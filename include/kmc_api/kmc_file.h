/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  Cuttlefish support: Jamshed Khan, Rob Patro

  Version: 3.1.1
  Date   : 2019-05-19
 */

#ifndef _KMC_FILE_H
#define _KMC_FILE_H

#include "kmer_defs.h"
#include "kmer_api.h"
#include <array>
#include <string>
#include <vector>
#include <unistd.h>
#include <utility>


// A class to imitate the KMC3 prefix-file access as if it's in memory.
// Note: only linear indexing is supported.
class Virtual_Prefix_File
{
private:

	static constexpr size_t buffer_elem_count = (1 << 16);	// Number of prefixes to be kept in memory buffer at a time.
	static constexpr size_t buffer_sz = buffer_elem_count * sizeof(uint64_t); // Size of buffer in bytes: 512KB. TODO: try small benchmarking for this size.
	size_t prefix_file_elem_count;	// Size of the KMC3 prefix-file (*.kmc_pre) in elements (i.e. 64-bit prefixes).
	std::array<uint64_t, buffer_elem_count> prefix_file_buf;	// The in-memory prefix-file buffer.

	uint64_t lut_area_size_in_bytes;	// From KMC3.
	size_t prefix_chunk_start_index;	// The index into the prefix-file where the prefix chunk currently loaded into memory starts.
	size_t prefix_chunk_end_index;	// The (non-inclusive) index into the prefix-file where the prefix chunk currently loaded into memory ends.

	uint64_t total_kmers;	// Total number of k-mers in the KMC3 database.
	FILE* fp;	// File handle to the KMC3 prefix-file.


	// Reads in as much data as possible from the prefix-file into the in-memory buffer,
	// and returns the number of elements read.
	size_t read_prefixes();


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
	// `lut_area_bytes` amount of bytes for its prefix-content, and the associated KMC3
	// database has `kmer_count` number of k-mers.
	void init(FILE*& fptr, uint64_t lut_area_bytes, uint64_t kmer_count);

	// Returns the data at index `idx` of the prefix-file.
	uint64_t operator[](size_t idx);
};


inline Virtual_Prefix_File::Virtual_Prefix_File():
	prefix_file_elem_count(0),
	lut_area_size_in_bytes(0),
	prefix_chunk_start_index(0),
	prefix_chunk_end_index(0),
	total_kmers(0),
	fp(nullptr)
{}


inline Virtual_Prefix_File::~Virtual_Prefix_File()
{
	if(fp)
	{ 
		std::fclose(fp); 
		fp = nullptr;
	}
}


inline void Virtual_Prefix_File::init(FILE*& fptr, const uint64_t lut_area_bytes, const uint64_t kmer_count)
{
	// *Take ownership* of `fptr`.
	fp = fptr;
	fptr = NULL;
	
	// Skip the first 4 bytes of header to get to the start of the prefixes.
	my_fseek(fp, +4, SEEK_CUR);
	lut_area_size_in_bytes = lut_area_bytes;
	prefix_file_elem_count = (lut_area_size_in_bytes + 8) / sizeof(uint64_t);	// What's that extra 1 element for? KMC3 comment: reads without 4 bytes of a header_offset (and without markers)
	total_kmers = kmer_count;

	// Read in some prefix-file data, and initialize the virtual indices into the prefix-file.
	prefix_chunk_start_index = 0;
	prefix_chunk_end_index = read_prefixes();
}


inline size_t Virtual_Prefix_File::read_prefixes()
{
	const size_t elems_to_read = std::min(prefix_file_elem_count - prefix_chunk_end_index, buffer_elem_count);
	const size_t bytes_to_read = elems_to_read * sizeof(uint64_t);
	const size_t bytes_read = std::fread(prefix_file_buf.data(), 1, bytes_to_read, fp);
	
	if(bytes_read != bytes_to_read)
	{
		std::cerr << "Error reading the KMC database prefix file. Aborting.\n";
		std::exit(EXIT_FAILURE);
	}

	return elems_to_read;
}


inline uint64_t Virtual_Prefix_File::operator[](const size_t idx)
{
	if(idx >= prefix_file_elem_count) 
		return total_kmers + 1;
	
	if(idx >= prefix_chunk_end_index)
	{
		prefix_chunk_start_index = idx;
		prefix_chunk_end_index = idx + read_prefixes();
	}

	return prefix_file_buf[idx - prefix_chunk_start_index];
}


struct CKMCFileInfo
{
	uint32 kmer_length;
	uint32 mode;
	uint32 counter_size;
	uint32 lut_prefix_length;
	uint32 signature_len;	
	uint32 min_count;
	uint64 max_count;
	bool both_strands;
	uint64 total_kmers;
};

// Forward declare Cuttlefish's k-mer class; required to parse KMC raw binary k-mers to Cuttlefish format.
template <uint16_t k> class Kmer;

class CKMCFile
{
	enum open_mode {closed, opened_for_RA, opened_for_listing};
	open_mode is_opened;
	uint64 suf_file_left_to_read = 0; // number of bytes that are yet to read in a listing mode
	uint64 suffix_file_total_to_read = 0; // number of bytes that constitutes records in kmc_suf file
	bool end_of_file;

	FILE *file_pre;
	FILE *file_suf;

	uint64* prefix_file_buf;
	Virtual_Prefix_File prefix_virt_buf;	// Virtual file to read over the prefix file in a buffered manner; for Cuttlefish.
	uint64 prefix_file_buf_size;
	uint64 prefix_index;			// The current prefix's index in an array "prefix_file_buf", readed from *.kmc_pre
	uint32 single_LUT_size;			// The size of a single LUT (in no. of elements)

	uint32* signature_map;
	uint32 signature_map_size;
	
	uchar* sufix_file_buf;
	uint64 sufix_number;			// The sufix's number to be listed
	uint64 index_in_partial_buf;	// The current byte's number in an array "sufix_file_buf", for listing mode

	uint32 kmer_length;
	uint32 mode;
	uint32 counter_size;
	uint32 lut_prefix_length;
	uint32 signature_len;
	uint32 min_count;
	uint64 max_count;
	uint64 total_kmers;
	bool both_strands;

	uint32 kmc_version;
	uint32 sufix_size;		// sufix's size in bytes 
	uint32 sufix_rec_size;  // sufix_size + counter_size

	uint32 original_min_count;
	uint64 original_max_count;

	static uint64 part_size; // the size of a block readed to sufix_file_buf, in listing mode 
	
	bool BinarySearch(int64 index_start, int64 index_stop, const CKmerAPI& kmer, uint64& counter, uint32 pattern_offset);

	// Open a file, recognize its size and check its marker. Auxiliary function.
	bool OpenASingleFile(const std::string &file_name, FILE *&file_handler, uint64 &size, char marker[]);	

	// Recognize current parameters. Auxiliary function.
	bool ReadParamsFrom_prefix_file_buf(uint64 &size, bool load_pref_file = true, bool init_pref_buf = true);	

	// Reload a contents of an array "sufix_file_buf" for listing mode. Auxiliary function. 
	void Reload_sufix_file_buf();

	// Implementation of GetCountersForRead for kmc1 database format for both strands
	bool GetCountersForRead_kmc1_both_strands(const std::string& read, std::vector<uint32>& counters);

	// Implementation of GetCountersForRead for kmc1 database format without choosing canonical k-mer
	bool GetCountersForRead_kmc1(const std::string& read, std::vector<uint32>& counters);		

	using super_kmers_t = std::vector<std::tuple<uint32, uint32, uint32>>;//start_pos, len, bin_no
	void GetSuperKmers(const std::string& transformed_read, super_kmers_t& super_kmers);

	// Implementation of GetCountersForRead for kmc2 database format for both strands
	bool GetCountersForRead_kmc2_both_strands(const std::string& read, std::vector<uint32>& counters);

	// Implementation of GetCountersForRead for kmc2 database format
	bool GetCountersForRead_kmc2(const std::string& read, std::vector<uint32>& counters);
public:
		
	CKMCFile();
	~CKMCFile();

	// Open files *.kmc_pre & *.kmc_suf, read them to RAM, close files. *.kmc_suf is opened for random access
	bool OpenForRA(const std::string &file_name);

	// Open files *kmc_pre & *.kmc_suf, read *.kmc_pre to RAM, *.kmc_suf is buffered
	bool OpenForListing(const std::string& file_name);

	// Open files `*kmc_pre` & `*.kmc_suf`, read `*.kmc_pre` to RAM; `*.kmc_suf` is not buffered internally.
	bool open_for_cuttlefish_listing(const std::string& file_name);

	// Open files `*kmc_pre` & `*.kmc_suf`, and read KMC DB parameters to RAM.
	bool read_parameters(const std::string& file_name);

	// Returns the size of a suffix-record in disk (in bytes); i.e. suffix-size plus counter-size.
	uint32_t suff_record_size() const;

	// Returns the current prefix (i.e. the next one to be potentially parsed).
	uint64_t curr_prefix() const;

	// Returns the current suffix's index (i.e. the next one to be parsed).
	uint64_t curr_suffix_idx() const;

	// Reads up-to `max_bytes_to_read` bytes worth of raw suffix records into the buffer `suff_buf`.
	// The prefixes corresponding to these suffixes are read into `pref_buf`, in the form
	// <prefix, #corresponding_suffix>. Returns the number of suffixes read. `0` is returned if
	// error(s) occurred during the read.
	uint64_t read_raw_suffixes(uint8_t* suff_buf, std::vector<std::pair<uint64_t, uint64_t>>& pref_buf, size_t max_bytes_to_read);

	// Parses a raw binary k-mer from the `buf_idx`'th byte onward of the buffer `suff_buf`, into
	// the Cuttlefish k-mer object `kmer`. `prefix_it` points to a pair of the form <prefix, abundance>
	// where "abundance" is the count of remaining k-mers to be parsed having this "prefix". The
	// iterator is adjusted accordingly for the next parse operation from the buffers.
	template <uint16_t k> void parse_kmer_buf(std::vector<std::pair<uint64_t, uint64_t>>::iterator& prefix_it, const uint8_t* suff_buf, size_t buf_idx, Kmer<k>& kmer) const;

	// Return next kmer in CKmerAPI &kmer. Return its counter in float &count. Return true if not EOF
	bool ReadNextKmer(CKmerAPI &kmer, float &count);

	inline bool ReadNextKmer(CKmerAPI &kmer, uint64 &count); //for small k-values when counter may be longer than 4bytes
	
	inline bool ReadNextKmer(CKmerAPI &kmer, uint32 &count);

	inline bool ReadNextKmer(CKmerAPI &kmer);
	// Release memory and close files in case they were opened 
	bool Close();

	// Set the minimal value for a counter. Kmers with counters below this theshold are ignored
	bool SetMinCount(uint32 x);

	// Return a value of min_count. Kmers with counters below this theshold are ignored 
	uint32 GetMinCount(void);

	// Set the maximal value for a counter. Kmers with counters above this theshold are ignored
	bool SetMaxCount(uint32 x);

	// Return a value of max_count. Kmers with counters above this theshold are ignored 
	uint64 GetMaxCount(void);
	
	//Return true if kmc was run without -b switch.
	bool GetBothStrands(void);

	// Return the total number of kmers between min_count and max_count
	uint64 KmerCount(void);

	// Return the length of kmers
	uint32 KmerLength(void);

	// Set initial values to enable listing kmers from the begining. Only in listing mode
	bool RestartListing(void);

	// Return true if all kmers are listed
	bool Eof(void);

	// Return true if kmer exists. In this case return kmer's counter in count
	bool CheckKmer(CKmerAPI &kmer, float &count);

	bool CheckKmer(CKmerAPI &kmer, uint32 &count);

	bool CheckKmer(CKmerAPI &kmer, uint64 &count);

	// Return true if kmer exists
	bool IsKmer(CKmerAPI &kmer);

	// Set original (readed from *.kmer_pre) values for min_count and max_count
	void ResetMinMaxCounts(void);

	// Get current parameters from kmer_database
	bool Info(uint32 &_kmer_length, uint32 &_mode, uint32 &_counter_size, uint32 &_lut_prefix_length, uint32 &_signature_len, uint32 &_min_count, uint64 &_max_count, uint64 &_total_kmers);
	
	// Get current parameters from kmer_database
	bool Info(CKMCFileInfo& info) const;

	// Get counters for all k-mers in read
	bool GetCountersForRead(const std::string& read, std::vector<uint32>& counters);
	bool GetCountersForRead(const std::string& read, std::vector<float>& counters);
	private:
		uint32 count_for_kmer_kmc1(CKmerAPI& kmer);
		uint32 count_for_kmer_kmc2(CKmerAPI& kmer, uint32 bin_start_pos);
};

//-----------------------------------------------------------------------------------------------
// Read next kmer
// OUT: kmer - next kmer
// OUT: count - kmer's counter
// RET: true - if not EOF
//-----------------------------------------------------------------------------------------------
inline bool CKMCFile::ReadNextKmer(CKmerAPI &kmer)
{
	uint64 prefix_mask = (1 << 2 * lut_prefix_length) - 1; //for kmc2 db

	if(is_opened != opened_for_listing) {
		return false;
	}

	//do 
	{
		if(end_of_file)
			return false;
		
		if(sufix_number == prefix_file_buf[prefix_index + 1]) {
			prefix_index++;
			while (prefix_file_buf[prefix_index] == prefix_file_buf[prefix_index + 1])
				prefix_index++;
		}
	
		uint32 off = (sizeof(prefix_index) * 8) - (lut_prefix_length * 2) - kmer.byte_alignment * 2;
			
		uint64 temp_prefix = (prefix_index & prefix_mask) << off;	// shift prefix towards MSD. "& prefix_mask" necessary for kmc2 db format

		memset(kmer.kmer_data, 0, kmer.no_of_rows * sizeof(kmer.kmer_data[0]));

		kmer.kmer_data[0] = temp_prefix;			// store prefix in an object CKmerAPI

		//for(uint32 i = 1; i < kmer.no_of_rows; i++) 
		//	kmer.kmer_data[i] = 0;
		

		//read sufix:
		uint32 row_index = 0;
 		uint64 suf = 0;
	
		off = off - 8;
				
		if(index_in_partial_buf+sufix_size < part_size) {
			// unroll?
			for(uint32_t a = 0; a < sufix_size; a ++) {
				suf = sufix_file_buf[index_in_partial_buf++];
				suf <<= off;
				kmer.kmer_data[row_index] |= suf;
				if (off == 0) { //the end of a word in kmer_data 
					off = 64 ;
					row_index++;
				} 
				off -=8; 
			}
		} else {
			for (uint32 a = 0; a < sufix_size; a++) {
				if (index_in_partial_buf == part_size)
					Reload_sufix_file_buf();

				suf = sufix_file_buf[index_in_partial_buf++];
				suf <<= off;
				kmer.kmer_data[row_index] |= suf;

				if (off == 0) { //the end of a word in kmer_data
					off = 64;
					row_index++;
				}
				off -= 8;
			}
		}

		//read counter:
		if(index_in_partial_buf == part_size) {
			Reload_sufix_file_buf();
		}

		// counter MUST be a single byte for this to work
		// auto count = sufix_file_buf[index_in_partial_buf++];
		index_in_partial_buf++;
		sufix_number++;
		if(sufix_number == total_kmers) {
			end_of_file = true;
		}
	}// while((count < min_count) || (count > max_count));

	return true;
}

//-----------------------------------------------------------------------------------------------
// Read next kmer
// OUT: kmer - next kmer
// OUT: count - kmer's counter
// RET: true - if not EOF
//-----------------------------------------------------------------------------------------------
bool CKMCFile::ReadNextKmer(CKmerAPI &kmer, uint32 &count)
{
	uint64 prefix_mask = (1 << 2 * lut_prefix_length) - 1; //for kmc2 db

	if(is_opened != opened_for_listing)
		return false;
	do
	{
		if(end_of_file)
			return false;
		
		if(sufix_number == prefix_file_buf[prefix_index + 1]) 
		{
			prefix_index++;
						
			while (prefix_file_buf[prefix_index] == prefix_file_buf[prefix_index + 1])
				prefix_index++;
		}
	
		uint32 off = (sizeof(prefix_index) * 8) - (lut_prefix_length * 2) - kmer.byte_alignment * 2;
			
		uint64 temp_prefix = (prefix_index & prefix_mask) << off;	// shift prefix towards MSD. "& prefix_mask" necessary for kmc2 db format
		
		kmer.kmer_data[0] = temp_prefix;			// store prefix in an object CKmerAPI

		for(uint32 i = 1; i < kmer.no_of_rows; i++)
			kmer.kmer_data[i] = 0;

		//read sufix:
		uint32 row_index = 0;
 		uint64 suf = 0;
	
		off = off - 8;
				
 		for(uint32 a = 0; a < sufix_size; a ++)
		{
			if(index_in_partial_buf == part_size)
				Reload_sufix_file_buf();
						
			suf = sufix_file_buf[index_in_partial_buf++];
			suf = suf << off;
			kmer.kmer_data[row_index] = kmer.kmer_data[row_index] | suf;

			if (off == 0)				//the end of a word in kmer_data
			{
					off = 56;
					row_index++;
			}
			else
					off -=8;
		}
	
		//read counter:
		if(index_in_partial_buf == part_size)
			Reload_sufix_file_buf();
		
		count = sufix_file_buf[index_in_partial_buf++];

		for(uint32 b = 1; b < counter_size; b++)
		{
			if(index_in_partial_buf == part_size)
				Reload_sufix_file_buf();
			
			uint32 aux = 0x000000ff & sufix_file_buf[index_in_partial_buf++];
			aux = aux << 8 * ( b);
			count = aux | count;
		}
			
		sufix_number++;
	
		if(sufix_number == total_kmers)
			end_of_file = true;

		if (mode != 0)
		{
			float float_counter;
			memcpy(&float_counter, &count, counter_size);
			if ((float_counter < min_count) || (float_counter > max_count))
				continue;
			else
				break;
		}

	}
	while((count < min_count) || (count > max_count));

	return true;
}

//-----------------------------------------------------------------------------------------------
// Read next kmer
// OUT: kmer - next kmer
// OUT: count - kmer's counter
// RET: true - if not EOF
//-----------------------------------------------------------------------------------------------
inline bool CKMCFile::ReadNextKmer(CKmerAPI &kmer, uint64 &count)
{
	uint64 prefix_mask = (1 << 2 * lut_prefix_length) - 1; //for kmc2 db

	if (is_opened != opened_for_listing)
		return false;
	do
	{
		if (end_of_file)
			return false;

		if (sufix_number == prefix_file_buf[prefix_index + 1])
		{
			prefix_index++;

			while (prefix_file_buf[prefix_index] == prefix_file_buf[prefix_index + 1])
				prefix_index++;
		}

		uint32 off = (sizeof(prefix_index)* 8) - (lut_prefix_length * 2) - kmer.byte_alignment * 2;

		uint64 temp_prefix = (prefix_index & prefix_mask) << off;	// shift prefix towards MSD. "& prefix_mask" necessary for kmc2 db format

		kmer.kmer_data[0] = temp_prefix;			// store prefix in an object CKmerAPI

		for (uint32 i = 1; i < kmer.no_of_rows; i++)
			kmer.kmer_data[i] = 0;

		//read sufix:
		uint32 row_index = 0;
		uint64 suf = 0;

		off = off - 8;

		for (uint32 a = 0; a < sufix_size; a++)
		{
			if (index_in_partial_buf == part_size)
				Reload_sufix_file_buf();

			suf = sufix_file_buf[index_in_partial_buf++];
			suf = suf << off;
			kmer.kmer_data[row_index] = kmer.kmer_data[row_index] | suf;

			if (off == 0)				//the end of a word in kmer_data
			{
				off = 56;
				row_index++;
			}
			else
				off -= 8;
		}

		//read counter:
		if (index_in_partial_buf == part_size)
			Reload_sufix_file_buf();

		count = sufix_file_buf[index_in_partial_buf++];

		for (uint32 b = 1; b < counter_size; b++)
		{
			if (index_in_partial_buf == part_size)
				Reload_sufix_file_buf();

			uint64 aux = 0x000000ff & sufix_file_buf[index_in_partial_buf++];
			aux = aux << 8 * (b);
			count = aux | count;
		}

		sufix_number++;

		if (sufix_number == total_kmers)
			end_of_file = true;

	} while ((count < min_count) || (count > max_count));

	return true;
}



inline uint64_t CKMCFile::read_raw_suffixes(uint8_t* const suff_buf, std::vector<std::pair<uint64_t, uint64_t>>& pref_buf, const size_t max_bytes_to_read)
{
	if(is_opened != opened_for_listing)
		return 0;

	const size_t max_suff_count = max_bytes_to_read / suff_record_size();
	uint64_t suff_read_count = 0;	// Count of suffixes to be read into the buffer `suff_buf`.
	pref_buf.clear();

	while(!end_of_file)
	{
		if(prefix_virt_buf[prefix_index] > total_kmers)
			break;

		// This conditional might be removable, by fixing the last entry of `prefix_file_buf` to `total_kmers` during its initialization.
		// TODO: Check if setting `prefix_file_buf[last_data_index]` to `total_kmers` instead of `total_kmers + 1` (current scheme) breaks stuffs.
		const uint64_t suff_id_next = (prefix_virt_buf[prefix_index + 1] > total_kmers ? total_kmers : prefix_virt_buf[prefix_index + 1]);
		// const uint64_t suff_id_next = std::min(prefix_file_buf[prefix_index + 1], total_kmers);

		// There are this many k-mers with the prefix `prefix_index`.
		const uint64_t suff_to_read = suff_id_next - sufix_number;
		if(suff_to_read > 0)
		{
			const uint64_t prev_sufix_number = sufix_number;
			
			if(suff_read_count + suff_to_read <= max_suff_count)
			{
				suff_read_count += suff_to_read;
				sufix_number += suff_to_read;
				pref_buf.emplace_back(prefix_index, sufix_number - prev_sufix_number);

				if(sufix_number == total_kmers)
					end_of_file = true;
			}
			else
			{
				sufix_number += (max_suff_count - suff_read_count);
				suff_read_count = max_suff_count;
				pref_buf.emplace_back(prefix_index, sufix_number - prev_sufix_number);

				break;
			}

		}

		prefix_index++;
	}

	const size_t bytes_to_read = suff_read_count * suff_record_size();
	const size_t bytes_read = std::fread(suff_buf, 1, bytes_to_read, file_suf);
	if(bytes_read != bytes_to_read)
		return 0;

	suf_file_left_to_read -= bytes_read;

	return suff_read_count;
}


inline uint32_t CKMCFile::suff_record_size() const
{
	return sufix_rec_size;
}


inline uint64_t CKMCFile::curr_prefix() const
{
	return prefix_index;
}


inline uint64_t CKMCFile::curr_suffix_idx() const
{
	return sufix_number;
}


template <uint16_t k>
inline void CKMCFile::parse_kmer_buf(std::vector<std::pair<uint64_t, uint64_t>>::iterator& prefix_it, const uint8_t* const suff_buf, size_t buf_idx, Kmer<k>& kmer) const
{
	static constexpr uint16_t NUM_INTS = (k + 31) / 32;
	uint64_t kmc_data[NUM_INTS]{};

	// Check if we have exhausted the currrent prefix.
	if(prefix_it->second == 0)
		++prefix_it;

	const uint64_t prefix = prefix_it->first;
	prefix_it->second--;

	
	// TODO: make some of these constant class-fields, to avoid repeated calculations.
	const uint64_t prefix_mask = (1 << 2 * lut_prefix_length) - 1; //for kmc2 db
	constexpr uint8_t byte_alignment = (k % 4 != 0 ? 4 - (k % 4) : 0);
	uint32_t off = (sizeof(prefix) * 8) - (lut_prefix_length * 2) - byte_alignment * 2;
	const uint64_t temp_prefix = (prefix & prefix_mask) << off;	// shift prefix towards MSD. "& prefix_mask" necessary for kmc2 db format

	// Store prefix in a KMC alignment (differs in endianness from Cuttlefish's).
	kmc_data[0] = temp_prefix;


	// Parse suffix.
	uint32_t row_idx{0};
	uint64_t suff{0};

	off -= 8;
	for(uint32 a = 0; a < sufix_size; a++)
	{			
		suff = suff_buf[buf_idx++];
		suff = suff << off;
		kmc_data[row_idx] = kmc_data[row_idx] | suff;

		if(off == 0)				//the end of a word in kmer_data
		{
			off = 56;
			row_idx++;
		}
		else
			off -= 8;
	}

	// Skip counter.
	// buf_idx += counter_size;

	// Parse KMC raw-binary k-mer data to Cuttlefish's k-mer format.
	kmer.from_KMC_data(kmc_data);
}

#endif

// ***** EOF
