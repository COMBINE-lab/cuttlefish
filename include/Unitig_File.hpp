
#ifndef UNITIG_FILE_HPP
#define UNITIG_FILE_HPP



#include "Virtual_File.hpp"
#include "Maximal_Unitig_Scratch.hpp"
#include "globals.hpp"
#include "utility.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/string.hpp"
#include "cereal/archives/binary.hpp"

#include <cstdint>
#include <cstddef>
#include <limits>
#include <cstring>
#include <vector>
#include <string>
#include <utility>
#include <ios>
#include <fstream>
#include <cstdlib>
#include <cassert>


namespace cuttlefish
{

// =============================================================================
// Unitig-file writer manager.
class Unitig_File_Writer
{
private:

    static constexpr std::size_t in_memory_bytes = 16lu * 1024; // 16 KB.
    static constexpr auto in_memory_len = in_memory_bytes / sizeof(uni_len_t);

    const std::string file_path;    // Path to the file for the unitig content.

    std::vector<char> buf;  // In-memory buffer for the unitig content. TODO: replace with `Buffer`.

    std::size_t total_sz;   // Total size of the added unitig content.
    std::vector<uni_len_t> len; // Lengths of the unitigs in the file.  TODO: replace with `Buffer`.
    std::size_t unitig_c;   // Number of unitigs added.
    std::ofstream output; // The unitig file.
    std::ofstream output_len;   // The lengths file.


    // Returns path to the file containing the lengths of the unitigs.
    const std::string length_file_path() const { return file_path + std::string(".len"); }

    // Flushes the in-memory unitig content to external memory.
    void flush_unitigs();

    // Flushes the in-memory unitig lengths to external memory.
    void flush_lengths();


public:

    // Constructs a unitig-writer to the file at path `file_path`.
    Unitig_File_Writer(const std::string& file_path);

    // Dummy constructor required for `cereal` deserialization to work for
    // objects containing this writer.
    Unitig_File_Writer();

    // Returns the total size of the added unitig content.
    auto size() const { return total_sz; }

    // Returns the number of unitigs added.
    auto unitig_count() const { return unitig_c; }

    // Adds the unitig content in the range `[beg, end)` to the writer.
    template <typename T_it_> void add(T_it_ beg, T_it_ end);

    // Adds the unitig content split in the ranges `[beg_1, end_1)` and `[beg_2,
    // end_2)` to the writer.
    template <typename T_it_> void add(T_it_ beg_1, T_it_ end_1, T_it_ beg_2, T_it_ end_2);

    // Closes the stream.
    void close();

    // Returns the resident set size of the space-dominant components of this
    // writer.
    std::size_t RSS() const;

    // Serializes the writer to the `cereal` archive `archive`.
    template <typename T_archive_> void save(T_archive_& archive) const;

    // Deserializes the writer from the `cereal` archive `archive`.
    template <typename T_archive_> void load(T_archive_& archive);
};


// =============================================================================
// Unitig-file reader manager.
class Unitig_File_Reader
{
private:

    static constexpr std::streamsize in_memory_bytes = 16lu * 1024; // 16 KB.

    const std::string file_path;    // Path to the file with the unitig content.

    std::vector<char> buf;  // In-memory buffer for the unitig content. TODO: replace with `Buffer`.
    std::vector<uni_len_t> uni_len; // Sizes of the unitigs in the current buffer.  TODO: replace with `Buffer`.

    std::ifstream input;    // The unitigs-file.
    Virtual_File<uni_len_t> len;    // The lengths-file.

    std::size_t buf_idx;    // Index into the unitig-buffer for the next unitig to read-in.
    uni_idx_t uni_idx_in_file;  // Index of the next unitig to read from file.
    uni_idx_t uni_idx_in_mem;   // Index of the next unitig to parse from buffer.

    const uni_idx_t unitig_count_;  // Number of unitigs in the file.
    uni_idx_t unitig_parsed_;   // Number of unitigs parsed.
    std::size_t total_sz;   // Total size of the read unitig content.


    // Returns path to the file containing the lengths of the unitigs.
    const std::string length_file_path() const { return file_path + std::string(".len"); }


public:

    // Constructs a unitig-reader for the file at path `file_path`.
    Unitig_File_Reader(const std::string& file_path);

    // Returns the number of unitigs in the file.
    auto unitig_count() const { return unitig_count_; }

    // Reads the next unitig into `unitig` and returns its length iff there were
    // unitigs remaining to be read. Returns 0 otherwise.
    std::size_t read_next_unitig(Buffer<char>& unitig);

    // Removes the unitig-files.
    void remove_files();
};


// =============================================================================
// Distributor of unitig-write operations over multiple write-managers.
class Unitig_Write_Distributor
{
private:

    const std::size_t writer_count; // Number of write-managers.
    std::vector<Padded<Unitig_File_Writer>> writer; // Collection of the different write-managers.
    const std::size_t worker_count; // Number of workers doing the writings.
    const std::size_t writer_per_worker;    // Number of write-managers dedicated to a worker.
    std::vector<Padded<std::size_t>> next_writer;   // `next_writer[w]` contains the relative-ID of the next writer-manager to be used by worker `w`.

    std::vector<Padded<Unitig_File_Writer>> mtig_writer;    // Collection of write-managers for trivially maximal unitigs.

public:

    // Constructs a unitig-writer distributor to `writer_count` write-managers
    // for `worker_count` workers. The files are at path-prefix `path_pref`.
    // `trivial_mtigs` denotes whether trivially maximal unitigs are to be
    // written or not.
    Unitig_Write_Distributor(const std::string& path_pref, std::size_t writer_count, std::size_t worker_count, bool trivial_mtigs);

    // Dummy constructor required for `cereal` deserialization to work for
    // objects containing this distributor.
    Unitig_Write_Distributor(const cereal::BinaryInputArchive&);

    // Adds the unitig content in the scratch `maximal_unitig` to the writer
    // for the `w_id`'th worker. Returns the pair `(b, idx)` such that `b` is
    // the ID of the bucket where the unitig is put in at the index `idx`.
    template <uint16_t k> std::pair<std::size_t, std::size_t> add(std::size_t w_id, const Maximal_Unitig_Scratch<k>& maximal_unitig);

    // Adds the k-mer content in `kmer` to the writer for the `w_id`'th worker.
    // Returns the pair `(b, idx)` such that `b` is the ID of the bucket where
    // the unitig is put in at the index `idx`.
    template <uint16_t k> std::pair<std::size_t, std::size_t> add(std::size_t w_id, const Kmer<k>& kmer);

    // Adds the trivially maximal unitig in the scratch `maximal_unitig` to the
    // writer for the `w_id`'th worker. Returns the pair `(b, idx)` such that
    // `b` is the ID of the bucket where the unitig is put in at the index
    // `idx`.
    template <uint16_t k> std::pair<std::size_t, std::size_t> add_trivial_mtig(std::size_t w_id, const Maximal_Unitig_Scratch<k>& maximal_unitig);

    // Returns the total number of buckets used for writing.
    auto bucket_count() const { return  writer.size() + mtig_writer.size(); }

    // Closes the unitig-writer streams.
    void close();

    // Returns the resident set size of the space-dominant components of this
    // distributor.
    std::size_t RSS() const;

    // (De)serializes the distributor from / to the `cereal` archive `archive`.
    template <typename T_archive_> void serialize(T_archive_& archive);
};


template <typename T_it_>
inline void Unitig_File_Writer::add(const T_it_ beg, const T_it_ end)
{
    len.push_back(end - beg);
    buf.insert(buf.end(), beg, end);
    total_sz += (end - beg);
    unitig_c++;
    assert((end - beg) <= std::numeric_limits<uni_len_t>::max());

    if(buf.size() >= in_memory_bytes)
        flush_unitigs();

    if(len.size() >= in_memory_len)
        flush_lengths();
}


template <typename T_it_>
inline void Unitig_File_Writer::add(const T_it_ beg_1, const T_it_ end_1, const T_it_ beg_2, const T_it_ end_2)
{
    const auto l = (end_1 - beg_1) + (end_2 - beg_2);
    len.push_back(l);
    buf.insert(buf.end(), beg_1, end_1);
    buf.insert(buf.end(), beg_2, end_2);
    total_sz += l;
    unitig_c++;
    assert(l <= std::numeric_limits<uni_len_t>::max());

    if(buf.size() >= in_memory_bytes)
        flush_unitigs();

    if(len.size() >= in_memory_len)
        flush_lengths();
}


inline void Unitig_File_Writer::flush_unitigs()
{
    output.write(reinterpret_cast<const char*>(buf.data()), buf.size() * sizeof(decltype(buf)::value_type));
    if(!output)
    {
        std::cerr << "Error writing of unitig content to " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    buf.clear();
}


inline void Unitig_File_Writer::flush_lengths()
{
    output_len.write(reinterpret_cast<const char*>(len.data()), len.size() * sizeof(decltype(len)::value_type));
    if(!output_len)
    {
        std::cerr << "Error writing of unitig lengths to " << length_file_path() << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    len.clear();
}


template <typename T_archive_>
inline void Unitig_File_Writer::save(T_archive_& archive) const
{
    archive(file_path, buf, total_sz, len, unitig_c);
}


template <typename T_archive_>
inline void Unitig_File_Writer::load(T_archive_& archive)
{
    archive(type::mut_ref(file_path), buf, total_sz, len, unitig_c);

    if(!file_path.empty())
    {
        output.open(file_path, std::ios::binary | std::ios::app);
        output_len.open(length_file_path(), std::ios::binary | std::ios::app);

        if(!output || !output_len)
        {
            std::cerr << "Error opening unitig-files at " << file_path << " and " << length_file_path() << ". Aborting.\n";
            std::exit(EXIT_FAILURE);
        }
    }
}


inline std::size_t Unitig_File_Reader::read_next_unitig(Buffer<char>& unitig)
{
    if(buf_idx == buf.size())   // Buffer has been parsed completely; try a re-read.
    {
        assert(uni_idx_in_mem == uni_len.size());

        if(uni_idx_in_file == unitig_count_)    // All unitigs have been read-off.
            return 0;

        std::streamsize bytes_to_read = 0;
        uni_len.clear();

        while(uni_idx_in_file < unitig_count_ && bytes_to_read < in_memory_bytes)
        {
            uni_len.push_back(len[uni_idx_in_file]);
            bytes_to_read += uni_len.back();

            uni_idx_in_file++;
        }

        assert(bytes_to_read > 0);
        buf.resize(bytes_to_read);
        input.read(buf.data(), bytes_to_read * sizeof(char));
        assert(input.gcount() == bytes_to_read);
        if(!input)
        {
            std::cerr << "Error reading of unitig content from " << file_path << ". Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        buf_idx = 0;
        uni_idx_in_mem = 0;
    }


    const auto len = uni_len[uni_idx_in_mem];
    unitig.reserve_uninit(len);
    std::memcpy(unitig.data(), buf.data() + buf_idx, len * sizeof(decltype(buf)::value_type));

    buf_idx += len;
    uni_idx_in_mem++;
    unitig_parsed_++;
    total_sz += len;
    assert(buf_idx <= buf.size());

    return len;
}


template <uint16_t k>
inline std::pair<std::size_t, std::size_t> Unitig_Write_Distributor::add(const std::size_t w_id, const Maximal_Unitig_Scratch<k>& maximal_unitig)
{
    auto& next = next_writer[w_id].unwrap();
    const auto writer_id = w_id * writer_per_worker + next + 1; // +1 as edge-partition 0 conceptually contains edges without any associated lm-tig (i.e. with weight > 1)
    // TODO: instead of rolling `writer_id` in every add, do it in batches (say of size 64) to enhance cache-performance.
    const auto range_size = (w_id < worker_count - 1 ? writer_per_worker : writer_count - (w_id * writer_per_worker));
    assert(range_size > 0);
    assert(next < range_size);

    if(++next == range_size)
        next = 0;

    const auto& u_f = maximal_unitig.unitig_label(side_t::front);
    const auto& u_b = maximal_unitig.unitig_label(side_t::back);
    assert(writer_id < writer.size());
    auto& w = writer[writer_id].unwrap();
    const auto idx = w.unitig_count();
    w.add(u_f.cbegin(), u_f.cend(), u_b.cbegin() + k, u_b.cend());

    return {writer_id, idx};
}


template <uint16_t k>
inline std::pair<std::size_t, std::size_t> Unitig_Write_Distributor::add(std::size_t w_id, const Kmer<k>& kmer)
{
    auto& next = next_writer[w_id].unwrap();
    const auto writer_id = w_id * writer_per_worker + next + 1; // +1 as edge-partition 0 conceptually contains edges without any associated lm-tig (i.e. with weight > 1)
    const auto range_size = (w_id < worker_count - 1 ? writer_per_worker : writer_count - (w_id * writer_per_worker));
    assert(range_size > 0);
    assert(next < range_size);

    if(++next == range_size)
        next = 0;

    // TODO: more efficient? Reuse the string?
    std::string label;
    kmer.get_label(label);
    assert(writer_id < writer.size());
    auto& w = writer[writer_id].unwrap();
    const auto idx = w.unitig_count();
    w.add(label.cbegin(), label.cend());

    return {writer_id, idx};
}


template <uint16_t k>
inline std::pair<std::size_t, std::size_t> Unitig_Write_Distributor::add_trivial_mtig(const std::size_t w_id, const Maximal_Unitig_Scratch<k>& maximal_unitig)
{
    assert(w_id < mtig_writer.size());
    const auto& u_f = maximal_unitig.unitig_label(side_t::front);
    const auto& u_b = maximal_unitig.unitig_label(side_t::back);
    auto& w = mtig_writer[w_id].unwrap();
    const auto idx = w.unitig_count();
    w.add(u_f.cbegin(), u_f.cend(), u_b.cbegin() + k, u_b.cend());

    return {writer.size() + w_id, idx};
}


template <typename T_archive_>
inline void Unitig_Write_Distributor::serialize(T_archive_& archive)
{
    archive(type::mut_ref(writer_count), writer, type::mut_ref(worker_count), type::mut_ref(writer_per_worker), next_writer, mtig_writer);
}

}



#endif
