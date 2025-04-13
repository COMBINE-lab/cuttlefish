
#ifndef EDGE_MATRIX_HPP
#define EDGE_MATRIX_HPP



#include "Discontinuity_Edge.hpp"
#include "Ext_Mem_Bucket.hpp"
#include "utility.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/string.hpp"
#include "cereal/archives/binary.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>


namespace cuttlefish
{

// =============================================================================
// Blocked edge-matrix of a discontinuity-graph of `k`-mers.
template <uint16_t k>
class Edge_Matrix
{
private:

    const std::size_t vertex_part_count_;   // Number of vertex-partitions in the graph; it needs to be a power of 2.
    const std::string path; // File-path prefix to the external-memory blocks of the matrix.
    std::vector<std::vector<Ext_Mem_Bucket_Concurrent<Discontinuity_Edge<k>>>> edge_matrix; // Blocked edge matrix.
    // TODO: do the cells need padding, or do the pads within the concurrent buckets suffice?

    mutable std::vector<std::size_t> row_to_read;   // `j`'th entry contains the row of the next block to read from column `j`.
    mutable std::vector<std::size_t> col_to_read;   // `i`'th entry contains the column of the next block to read from row `i`.


public:

    // Constructs a blocked edge-matrix for `part_count` vertex-partitions. The
    // partition-count needs to be a power of 2.
    Edge_Matrix(std::size_t part_count, const std::string& path);

    // Dummy constructor required for `cereal` deserialization to work for
    // objects containing this matrix.
    Edge_Matrix(const cereal::BinaryInputArchive&);

    // Returns the number of vertex-partitions in the graph.
    std::size_t vertex_part_count() const { return vertex_part_count_; }

    // Returns the partition ID for the k-mer `kmer`.
    std::size_t partition(const Kmer<k>& kmer) const;

    // Adds the discontinuity-edge `({(u, s_u), (v, s_v)}, w, b)` to the matrix.
    // `b_idx` is the index of the corresponding unitig in its bucket.
    // `u_is_phi` and `v_is_phi` denote whether the `u` and the `v` endpoints
    // are ϕ, respectively.
    // TODO: why not using refs here for the k-mers? Check.
    void add(Kmer<k> u, side_t s_u, Kmer<k> v, side_t s_v, weight_t w, uint16_t b, std::size_t b_idx, bool u_is_phi, bool v_is_phi);

    // Serializes and closes the matrix. Edges should not be added anymore
    // after this has been invoked.
    void close();

    // Reads the edges from the `[j, j]`'th block into `buf`.
    void read_diagonal_block(std::size_t j, std::vector<Discontinuity_Edge<k>>& buf) const;

    // Reads the edges from the `[j, j]`'th block into `buf` and returns the
    // number of edges read.
    std::size_t read_diagonal_block(std::size_t j, Buffer<Discontinuity_Edge<k>>& buf) const;

    // Reads a chunk of edges from the column `j` into `buf`. Returns the count
    // of edges read. If `0` is returned, then the column has been depleted.
    // NB: this does not read the blocks in the diagonal.
    std::size_t read_column_buffered(std::size_t j, Buffer<Discontinuity_Edge<k>>& buf) const;

    // Reads a chunk of edges from the `(x, y)`'th block into `buf`, and returns
    // the count of edges read. If it is `0`, then the block has been depleted.
    std::size_t read_block_buffered(std::size_t x, std::size_t y, Buffer<Discontinuity_Edge<k>>& buf, std::size_t n) const;

    // Resets the read-status of each worker for the entire matrix.
    void reset_read();

    // Reads a chunk of edges from the row `i` into `buf`. Returns the count of
    // edges read. If `0` is returned, then the column has been depleted.
    // NB: this does not read the blocks in the diagonal.
    std::size_t read_row_buffered(std::size_t i, Buffer<Discontinuity_Edge<k>>& buf) const;

    // Reads the edges from the `[i, j]`'th block into `buf`. Returns the count
    // of edges read.
    std::size_t read_block(std::size_t i, std::size_t j, Buffer<Discontinuity_Edge<k>>& buf) const;

    // Returns the number of edges stored in row `i`.
    std::size_t row_size(std::size_t i) const;

    // Returns the number of edges stored in column `j`.
    std::size_t col_size(std::size_t j) const;

    // Returns the number of edges in the `[i, j]`'th block.
    std::size_t block_size(std::size_t i, std::size_t j) const;

    // Returns the number of edges stores in the matrix.
    std::size_t size() const;

    // Returns the maximum block-size of the matrix.
    std::size_t max_block_size() const;

    // Removes the `[i, j]`'th block.
    void remove_block(std::size_t i, std::size_t j);

    // Returns the resident set size of the space-dominant components of this
    // matrix.
    std::size_t RSS() const;

    // (De)serializes the matrix from / to the `cereal` archive `archive`.
    template <typename T_archive_> void serialize(T_archive_& archive);
};


template <uint16_t k>
inline std::size_t Edge_Matrix<k>::partition(const Kmer<k>& kmer) const
{
    return (kmer.to_u64() & (vertex_part_count_ - 1)) + 1;
}


template <uint16_t k>
inline void Edge_Matrix<k>::add(const Kmer<k> u, const side_t s_u, const Kmer<k> v, const side_t s_v, const weight_t w, const uint16_t b, const std::size_t b_idx, const bool u_is_phi, const bool v_is_phi)
{
    // TODO: add batched insertion per worker instead of locking a cell at every insertion.
    // TODO: almost surely the lack of batched insertion is the scalability bottleneck.

    auto p = u_is_phi ? 0 : partition(u);
    auto q = v_is_phi ? 0 : partition(v);

    if(p <= q)
        edge_matrix[p][q].emplace(u, s_u, v, s_v, w, b, b_idx, u_is_phi, v_is_phi, side_t::back);
    else    // p and q needs to be swapped along with the edge endpoints
        edge_matrix[q][p].emplace(v, s_v, u, s_u, w, b, b_idx, v_is_phi, u_is_phi, side_t::front);
}


template <uint16_t k>
template <typename T_archive_>
inline void Edge_Matrix<k>::serialize(T_archive_& archive)
{
    archive(type::mut_ref(vertex_part_count_), type::mut_ref(path), edge_matrix, row_to_read, col_to_read);
}

}



#endif
