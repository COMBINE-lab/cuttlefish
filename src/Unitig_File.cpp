
#include "Unitig_File.hpp"
#include "utility.hpp"

#include <cstdlib>
#include <algorithm>
#include <filesystem>


namespace cuttlefish
{

Unitig_File_Writer::Unitig_File_Writer(const std::string& file_path):
      file_path(file_path)
    , total_sz(0)
    , unitig_c(0)
{
    if(!file_path.empty())
    {
        output.open(file_path, std::ios::binary);
        output_len.open(length_file_path(), std::ios::binary);

        if(!output || !output_len)
        {
            std::cerr << "Error opening unitig-files at " << file_path << " and " << length_file_path() << ". Aborting.\n";
            std::exit(EXIT_FAILURE);
        }
    }
}


Unitig_File_Writer::Unitig_File_Writer(): Unitig_File_Writer(std::string())
{}


void Unitig_File_Writer::close()
{
    if(!buf.empty())
        flush_unitigs();

    if(!len.empty())
        flush_lengths();

    output.close();
    output_len.close();
}


std::size_t Unitig_File_Writer::RSS() const
{
    return memory::RSS(buf) + memory::RSS(len);
}


Unitig_File_Reader::Unitig_File_Reader(const std::string& file_path):
      file_path(file_path)
    , input(file_path)
    , len(length_file_path().c_str())
    , buf_idx(0)
    , uni_idx_in_file(0)
    , uni_idx_in_mem(0)
    , unitig_count_(len.size())
    , unitig_parsed_(0)
    , total_sz(0)
{}


void Unitig_File_Reader::remove_files()
{
    if(!remove_file(file_path) || !remove_file(length_file_path()))
    {
        std::cerr << "Error removing file(s) at " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


Unitig_Write_Distributor::Unitig_Write_Distributor(const std::string& path_pref, const std::size_t writer_count, const std::size_t worker_count, const bool trivial_mtigs):
      writer_count(writer_count)
    , worker_count(worker_count)
    , writer_per_worker(writer_count / worker_count)
    , next_writer(worker_count, 0)
{
    assert(writer_count >= worker_count);

    std::filesystem::create_directories(path_pref);
    writer.reserve(1 + writer_count);
    writer.emplace_back(std::string()); // Edge-partition 0 is symbolic, for edges that do not have any associated lm-tig (i.e. has weight > 1).
    for(std::size_t b = 1; b <= writer_count; ++b)
        writer.emplace_back(path_pref + "/" + std::to_string(b));

    if(trivial_mtigs)
    {
        mtig_writer.reserve(worker_count);
        for(std::size_t w = 0; w < worker_count; ++w)
            mtig_writer.emplace_back(path_pref + "/" + std::to_string(writer.size() + w));
    }
}


Unitig_Write_Distributor::Unitig_Write_Distributor(const cereal::BinaryInputArchive&):
      writer_count()
    , worker_count()
    , writer_per_worker()
{}


void Unitig_Write_Distributor::close()
{
    for(auto& w : writer)
        w.unwrap().close();

    if(!mtig_writer.empty())
        for(auto& w : mtig_writer)
            w.unwrap().close();
}


std::size_t Unitig_Write_Distributor::RSS() const
{
    std::size_t writer_bytes = 0;
    std::for_each(writer.cbegin(), writer.cend(), [&](const auto& v){ writer_bytes += v.unwrap().RSS(); });
    std::for_each(mtig_writer.cbegin(), mtig_writer.cend(), [&](const auto& v){ writer_bytes += v.unwrap().RSS(); });

    return writer_bytes;
}

}
