
#ifndef VIRTUAL_PREFIX_FILE_HPP
#define VIRTUAL_PREFIX_FILE_HPP



#include <cstdint>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <iostream>


// A class to imitate the KMC3 prefix-file access as if it's in memory.
// Note: only linear indexing is supported.
class Virtual_Prefix_File
{
private:

	static constexpr std::size_t buffer_elem_count = (1 << 21);	// Number of prefixes to be kept in memory buffer at a time.
	std::size_t prefix_file_elem_count;	// Size of the KMC3 prefix-file (*.kmc_pre) in elements (i.e. 64-bit prefixes).
	std::vector<uint64_t> prefix_file_buf;	// The in-memory prefix-file buffer.

	std::size_t prefix_chunk_start_index;	// The index into the prefix-file where the prefix chunk currently loaded into memory starts.
	std::size_t prefix_chunk_end_index;	// The (non-inclusive) index into the prefix-file where the prefix chunk currently loaded into memory ends.

	uint64_t total_kmers;	// Total number of k-mers in the KMC3 database.
	std::FILE* fp;	// File handle to the KMC3 prefix-file.


	// Reads in as much data as possible from the prefix-file into the in-memory buffer,
	// and returns the number of elements read.
	std::size_t read_prefixes();


public:

	// Constructs an empty virtual file buffer.
	Virtual_Prefix_File();

	// Invalidate move and copy constructors, and copy-assignment operators.
	Virtual_Prefix_File(Virtual_Prefix_File&& rhs) = delete;
	Virtual_Prefix_File(const Virtual_Prefix_File& rhs) = delete;
	Virtual_Prefix_File& operator=(const Virtual_Prefix_File& rhs) = delete;
	Virtual_Prefix_File& operator=(Virtual_Prefix_File& rhs) = delete;

	// Destructs the virtual file.
	~Virtual_Prefix_File();

	// Initializes the file buffer with the file handle `fptr` that is supposed to contain
	// `prefix_count` number of prefixes as its prefix-content, and the associated KMC3
	// database has `kmer_count` number of k-mers.
	void init(std::FILE*& fptr, uint64_t prefix_count, uint64_t kmer_count);

	// Returns the data at index `idx` of the prefix-file.
	uint64_t operator[](std::size_t idx);

	// Returns the size of the buffer in bytes.
	static constexpr std::size_t memory();
};


inline std::size_t Virtual_Prefix_File::read_prefixes()
{
	const std::size_t elems_to_read = std::min(prefix_file_elem_count - prefix_chunk_end_index, buffer_elem_count);
	const std::size_t bytes_to_read = elems_to_read * sizeof(uint64_t);
	const std::size_t bytes_read = std::fread(prefix_file_buf.data(), 1, bytes_to_read, fp);
	
	if(bytes_read != bytes_to_read)
	{
		std::cerr << "Error reading the KMC database prefix file. Aborting.\n";
		std::exit(EXIT_FAILURE);
	}

	return elems_to_read;
}


inline uint64_t Virtual_Prefix_File::operator[](const std::size_t idx)
{
	if(idx >= prefix_file_elem_count)
		return total_kmers + 1;

	if(idx == prefix_file_elem_count - 1)
		return total_kmers;

	if(idx >= prefix_chunk_end_index)
	{
		prefix_chunk_start_index = idx;
		prefix_chunk_end_index = idx + read_prefixes();
	}

	return prefix_file_buf[idx - prefix_chunk_start_index];
}


inline constexpr std::size_t Virtual_Prefix_File::memory()
{
	return buffer_elem_count * sizeof(uint64_t);
}



#endif
