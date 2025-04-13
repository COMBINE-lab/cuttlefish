
#include "Edge_Matrix.hpp"
#include "globals.hpp"
#include "parlay/parallel.h"

#include <cstddef>
#include <filesystem>
#include <algorithm>
#include <cassert>


namespace cuttlefish
{

template <uint16_t k>
Edge_Matrix<k>::Edge_Matrix(std::size_t part_count, const std::string& path):
      vertex_part_count_(part_count)
    , path(path)
    , row_to_read(part_count + 1, 0)
    , col_to_read(part_count + 1)
{
    // TODO: fix better policy.
    if((part_count & (part_count - 1)) != 0)
    {
        std::cerr << "Vertex-partition count needs to be a power of 2. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    edge_matrix.resize(part_count + 1);

    for(std::size_t i = 0; i <= part_count; ++i)
    {
        const auto row_dir = path + "/" + std::to_string(i);
        std::filesystem::create_directories(row_dir);
        for(std::size_t j = 0; j <= part_count; ++j)
            if(j < i || (i == 0 && j == 0))
                edge_matrix[i].emplace_back();
            else
                edge_matrix[i].emplace_back(row_dir + "/" + std::to_string(j));
    }

    for(std::size_t i = 1; i <= part_count; ++i)
        col_to_read[i] = i + 1;
}


template <uint16_t k>
Edge_Matrix<k>::Edge_Matrix(const cereal::BinaryInputArchive&):
      vertex_part_count_()
{}


template <uint16_t k>
void Edge_Matrix<k>::close()
{
    parlay::parallel_for(1, edge_matrix.size(),
    [&](const auto i)
    {
        parlay::parallel_for(i, edge_matrix.size(),
        [&](const auto j)
        {
            edge_matrix[i][j].close();
        }, 1);
    }, 1);
}


template <uint16_t k>
void Edge_Matrix<k>::read_diagonal_block(const std::size_t j, std::vector<Discontinuity_Edge<k>>& buf) const
{
    edge_matrix[j][j].load(buf);
}


template <uint16_t k>
std::size_t Edge_Matrix<k>::read_diagonal_block(const std::size_t j, Buffer<Discontinuity_Edge<k>>& buf) const
{
    buf.reserve_uninit(edge_matrix[j][j].size());
    return edge_matrix[j][j].load(buf.data());
}


template <uint16_t k>
std::size_t Edge_Matrix<k>::read_column_buffered(const std::size_t j, Buffer<Discontinuity_Edge<k>>& buf) const
{
    if(row_to_read[j] >= j)
        return 0;

    const auto to_read = edge_matrix[row_to_read[j]][j].size();
    buf.reserve_uninit(to_read);
    const auto read = edge_matrix[row_to_read[j]][j].load(buf.data());
    assert(read == to_read);
    row_to_read[j]++;

    return read;
}


template <uint16_t k>
std::size_t Edge_Matrix<k>::read_row_buffered(const std::size_t i, Buffer<Discontinuity_Edge<k>>& buf) const
{
    if(col_to_read[i] > vertex_part_count_)
        return 0;

    const auto to_read = edge_matrix[i][col_to_read[i]].size();
    buf.reserve_uninit(to_read);
    const auto read = edge_matrix[i][col_to_read[i]].load(buf.data());
    assert(read == to_read);
    col_to_read[i]++;

    return read;
}


template <uint16_t k>
std::size_t Edge_Matrix<k>::read_block(std::size_t i, std::size_t j, Buffer<Discontinuity_Edge<k>>& buf) const
{
    assert(i <= vertex_part_count_ && j <= vertex_part_count_);
    assert(i <= j); // Upper-diagonal matrix.
    buf.reserve_uninit(edge_matrix[i][j].size());
    return edge_matrix[i][j].load(buf.data());
}


template <uint16_t k>
std::size_t Edge_Matrix<k>::read_block_buffered(const std::size_t x, const std::size_t y, Buffer<Discontinuity_Edge<k>>& buf, const std::size_t n) const
{
    assert(x <= vertex_part_count_ && y <= vertex_part_count_);
    assert(x <= y); // Upper-diagonal matrix.
    return edge_matrix[x][y].read_buffered(buf, n);
}


template <uint16_t k>
void Edge_Matrix<k>::reset_read()
{
    for(std::size_t i = 0; i < edge_matrix.size(); ++i)
        for(std::size_t j = i; j < edge_matrix[i].size(); ++j)
            edge_matrix[i][j].reset_read();
}


template <uint16_t k>
std::size_t Edge_Matrix<k>::row_size(const std::size_t i) const
{
    assert(i <= vertex_part_count_);

    std::size_t sz = 0;
    for(std::size_t j = std::max(1lu, i); j <= vertex_part_count_; ++j)
        sz += edge_matrix[i][j].size();

    return sz;
}


template <uint16_t k>
std::size_t Edge_Matrix<k>::col_size(const std::size_t j) const
{
    assert(j <= vertex_part_count_);

    std::size_t sz = 0;
    for(std::size_t i = 0; i <= j; ++i)
        sz += edge_matrix[i][j].size();

    return sz;
}


template <uint16_t k>
std::size_t Edge_Matrix<k>::block_size(const std::size_t i, const std::size_t j) const
{
    return edge_matrix[i][j].size();
}


template <uint16_t k>
std::size_t Edge_Matrix<k>::size() const
{
    std::size_t sz = 0;
    for(std::size_t i = 0; i <= vertex_part_count_; ++i)
        sz += row_size(i);

    return sz;
}


template <uint16_t k>
std::size_t Edge_Matrix<k>::max_block_size() const
{
    std::size_t max_sz = 0;
    for(std::size_t i = 0; i <= vertex_part_count_; ++i)
        for(std::size_t j = std::max(1lu, i); j <= vertex_part_count_; ++j)
            max_sz = std::max(max_sz, edge_matrix[i][j].size());

    return max_sz;
}


template <uint16_t k>
void Edge_Matrix<k>::remove_block(const std::size_t i, const std::size_t j)
{
    edge_matrix[i][j].remove();
}


template <uint16_t k>
std::size_t Edge_Matrix<k>::RSS() const
{
    std::size_t bucket_bytes = 0;
    std::for_each(edge_matrix.cbegin(), edge_matrix.cend(), [&](const auto& row)
    {
        std::for_each(row.cbegin(), row.cend(), [&](const auto& blk){ bucket_bytes += blk.RSS(); });
    });

    return bucket_bytes;
}

}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Edge_Matrix)
