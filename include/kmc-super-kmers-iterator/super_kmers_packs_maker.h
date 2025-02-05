#pragma once
#include "defs.h"
#include <list>

#include "bin_api.h"

class CExpanderPackQueue
{
	uint64_t _start = 0;
	std::list<simple_pack_data_t> l;

	std::mutex mtx;
public:
	CExpanderPackQueue(std::list<simple_pack_data_t>& expander_pack)
	{
		l = std::move(expander_pack);
	}
	bool Pop(uint64_t& n_super_kmers, uint64_t& start, uint64_t& end)
	{
		std::lock_guard<std::mutex> lck(mtx);
		if (l.empty())
			return false;
		start = _start;
		_start += l.front().end_pos;
		end = _start;

		n_super_kmers = l.front().n_super_kmers;

		l.pop_front();
		return true;
	}
};

class SuperKmersPacksMaker {
	bin_reader_super_kmers_packs_maker_queue_t& in;
	super_kmers_packs_maker_super_kmer_iterator_queue_t& out;
	CExpanderPackQueue pack_queue;
public:
	SuperKmersPacksMaker(bin_reader_super_kmers_packs_maker_queue_t& in, 
		super_kmers_packs_maker_super_kmer_iterator_queue_t& out,
		std::list<simple_pack_data_t>& expander_pack) :
		in(in),
		out(out),
		pack_queue(expander_pack) {

	}
	void Process() {
		uint64_t n_super_kmers;
		uint64_t start, end;

		reader_data_t in_pack;
		size_t in_pack_pos = 0;
		while (pack_queue.Pop(n_super_kmers, start, end)) {

			super_kmers_packer_data_t out_pack(n_super_kmers, end - start);
			
			//while outpack is not filled
			while (out_pack.data.size() < out_pack.data.capacity()) {

				//if input pack is finished, get next
				if (in_pack_pos == in_pack.data.size()) {

					in_pack_pos = 0;
					if (!in.pop(in_pack)) {
						std::cerr << "Error: Inconsistent data, there should be more data read from the input file\n";
						exit(1);
					}		
				}

				size_t yet_to_fill = out_pack.data.capacity() - out_pack.data.size();
				size_t left_in_current_pack = in_pack.data.size() - in_pack_pos;

				size_t to_copy = std::min(yet_to_fill, left_in_current_pack);
				std::copy(in_pack.data.data() + in_pack_pos, in_pack.data.data() + in_pack_pos + to_copy, out_pack.data.data() + out_pack.data.size());
				in_pack_pos += to_copy;
				out_pack.data.resize(out_pack.data.size() + to_copy);
				
			}
			out.push(std::move(out_pack));
		}
		out.mark_completed();
	}
};
