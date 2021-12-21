
#include "kmc_api/Virtual_Prefix_File.hpp"
#include "kmc_api/kmer_defs.h"


Virtual_Prefix_File::Virtual_Prefix_File():
	prefix_file_elem_count(0),
	lut_area_size_in_bytes(0),
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


void Virtual_Prefix_File::init(std::FILE*& fptr, const uint64_t lut_area_bytes, const uint64_t kmer_count)
{
	// *Take ownership* of `fptr`.
	fp = fptr;
	fptr = NULL;
	
	// Skip the first 4 bytes of header to get to the start of the prefixes.
	my_fseek(fp, +4, SEEK_CUR);
	lut_area_size_in_bytes = lut_area_bytes;
	prefix_file_elem_count = (lut_area_size_in_bytes + 8) / sizeof(uint64_t);	// What's that extra 1 element for? KMC3 comment: reads without 4 bytes of a header_offset (and without markers)
	total_kmers = kmer_count;

	// Allocate the prefix-file buffer.
	prefix_file_buf.reserve(buffer_elem_count);

	// Read in some prefix-file data, and initialize the virtual indices into the prefix-file.
	prefix_chunk_start_index = 0;
	prefix_chunk_end_index = read_prefixes();
}


void Virtual_Prefix_File::init_kmc1(std::FILE*& fptr, const uint64_t elems, const uint64_t kmer_count)
{
	// *Take ownership* of `fptr`.
	fp = fptr;
	fptr = NULL;
	
	// Set metadata.
	lut_area_size_in_bytes = elems * sizeof(uint64_t);
	prefix_file_elem_count = elems;
	total_kmers = kmer_count;

	// Allocate the prefix-file buffer.
	prefix_file_buf.reserve(buffer_elem_count);

	// Read in some prefix-file data, and initialize the virtual indices into the prefix-file.
	prefix_chunk_start_index = 0;
	prefix_chunk_end_index = read_prefixes();
}
