
#include "kmc_api/Virtual_Prefix_File.hpp"


Virtual_Prefix_File::Virtual_Prefix_File():
	prefix_file_elem_count(0),
	prefix_chunk_start_index(0),
	prefix_chunk_end_index(0),
	total_kmers(0),
	fp(nullptr)
{}


Virtual_Prefix_File::~Virtual_Prefix_File()
{
	if(fp)
	{ 
		std::fclose(fp); 
		fp = nullptr;
	}
}


void Virtual_Prefix_File::init(std::FILE*& fptr, const uint64_t prefix_count, const uint64_t kmer_count)
{
	// *Take ownership* of `fptr`.
	fp = fptr;
	fptr = NULL;
	
	prefix_file_elem_count = prefix_count;
	total_kmers = kmer_count;

	// Allocate the prefix-file buffer.
	prefix_file_buf.reserve(buffer_elem_count);

	// Read in some prefix-file data, and initialize the virtual indices into the prefix-file.
	prefix_chunk_start_index = 0;
	prefix_chunk_end_index = read_prefixes();
}
