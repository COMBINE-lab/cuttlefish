
#include "Super_Kmer_Bucket.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <cassert>


namespace cuttlefish
{

template <bool Colored_>
Super_Kmer_Bucket<Colored_>::Super_Kmer_Bucket(const uint16_t k, const uint16_t l, const std::string& path, const std::size_t chunk_cap):
      path_(path)
    , output(path_, std::ios::binary)
    , size_(0)
    , chunk_cap(chunk_cap)
    , chunk(k, l, chunk_cap)
    , bytes_(0)
    , compressed_bytes_(0)
{}


template <bool Colored_>
void Super_Kmer_Bucket<Colored_>::flush_chunk()
{
    if(!chunk.empty())
    {
        assert(output);
        // chunk.serialize(output);
        cmp_bytes.push_back(chunk.serialize_compressed(output));
        chunk_sz.push_back(chunk.size());

        bytes_ += chunk.bytes();
        compressed_bytes_ += (cmp_bytes.back().first + cmp_bytes.back().second);

        chunk.clear();
    }
}


template <bool Colored_>
void Super_Kmer_Bucket<Colored_>::close()
{
    flush_chunk();
    output.close();
}


template <bool Colored_>
void Super_Kmer_Bucket<Colored_>::remove()
{
    chunk.free();
    force_free(chunk_sz);
    force_free(cmp_bytes);

    if(output.is_open())
        output.close();

    if(!output || !remove_file(path_))
    {
        std::cerr << "Error removing file at " << path_ << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <bool Colored_>
Super_Kmer_Bucket<Colored_>::Iterator::Iterator(const Super_Kmer_Bucket& B):
      B(B)
    , input(B.path_, std::ios::in | std::ios::binary)
    , idx(0)
    , chunk_start_idx(0)
    , chunk_end_idx(0)
    , chunk_id(0)
{
    assert(input);
}


// TODO: inline.
template <bool Colored_>
std::size_t Super_Kmer_Bucket<Colored_>::Iterator::read_chunk()
{
    assert(chunk_end_idx < B.size());
    assert(chunk_id < B.chunk_sz.size());
    const auto super_kmers_to_read = B.chunk_sz[chunk_id];

    // B.chunk.deserialize(input, super_kmers_to_read);
    B.chunk.deserialize_decompressed(input, super_kmers_to_read, B.cmp_bytes[chunk_id]);
    chunk_id++;

    return super_kmers_to_read;
}


template <bool Colored_>
std::size_t Super_Kmer_Bucket<Colored_>::RSS() const
{
    return chunk.RSS() + memory::RSS(chunk_sz) + memory::RSS(cmp_bytes);
}

}



// Template-instantiations for the required instances.
template class cuttlefish::Super_Kmer_Bucket<false>;
template class cuttlefish::Super_Kmer_Bucket<true>;
