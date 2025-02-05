#ifndef _BITMEMORY_H
#define _BITMEMORY_H

#include <utility>
#include <memory>
#include <vector>
#include <cstdint>
#include <cassert>

namespace refresh
{
	// **************************************************************
	template<template<typename, typename...> typename STORAGE_T = std::vector>
	class basic_bit_memory
	{
	protected:
		STORAGE_T<uint64_t> stream;
		uint64_t buffer{};
		size_t in_buffer_pos{};

	public:
		basic_bit_memory() = default;

		basic_bit_memory(basic_bit_memory<STORAGE_T>&&) = default;
		basic_bit_memory& operator=(basic_bit_memory<STORAGE_T>&&) = default;

		basic_bit_memory(STORAGE_T<uint64_t>&& storage) :
			stream(std::move(storage))
		{

		}

		bool empty() {
			return stream.empty() && (in_buffer_pos == 0);
		}

		void clear() {
			stream.clear();
			in_buffer_pos = 0;
			buffer = 0;
		}

		void shrink_to_fit() {
			stream.shrink_to_fit();
		}
	};

	// **************************************************************
	template<template<typename, typename...> typename STORAGE_T = std::vector>
	class bit_memory_writer : public basic_bit_memory<STORAGE_T>
	{
		using basic_bit_memory<STORAGE_T>::buffer;
		using basic_bit_memory<STORAGE_T>::stream;
		using basic_bit_memory<STORAGE_T>::in_buffer_pos;

		void store_bits_imp(uint64_t x, size_t no_bits)
		{
			uint64_t b = x;
			b <<= in_buffer_pos;
			in_buffer_pos += no_bits;
			buffer += b;

			if (in_buffer_pos > 64)
			{
				stream.emplace_back(buffer);
				in_buffer_pos -= 64;

				buffer = x >> (no_bits - in_buffer_pos);
			}
			else if (in_buffer_pos == 64)
			{
				stream.emplace_back(buffer);
				in_buffer_pos = 0;
				buffer = 0;
			}
		}

	public:
		bit_memory_writer() = default;

		bit_memory_writer(bit_memory_writer<STORAGE_T>&&) = default;
		bit_memory_writer& operator=(bit_memory_writer<STORAGE_T>&&) = default;

		bit_memory_writer(STORAGE_T<uint64_t>&& storage) : basic_bit_memory<STORAGE_T>(std::move(storage)) {

		}

		void store_byte(uint8_t x) {
			store_bits_imp(x, 8);
		}

		void store_2bytes(uint16_t x) {
			store_bits_imp(x, 16);
		}

		void store_4bytes(uint32_t x) {
			store_bits_imp(x, 32);
		}

		void store_8bytes(uint64_t x) {
			store_bits_imp(x, 64);
		}

		// Assume that only no_bits lower bits in x are used
		void store_bits(uint64_t x, size_t no_bits) {
			store_bits_imp(x, no_bits);
		}

		// Mask higher (than no_bits) bits
		void store_bits_secure(uint64_t x, size_t no_bits) {
			x &= ~0ull >> (64 - no_bits);

			store_bits_imp(x, no_bits);
		}

		void flush_buffer() {
			if (in_buffer_pos == 0)
				return;

			stream.emplace_back(buffer);

			buffer = 0;
			in_buffer_pos = 0;
		}

		void flush_byte() {
			uint64_t t = in_buffer_pos & 7;

			if (t)
			{
				in_buffer_pos += 8 - t;

				if (in_buffer_pos == 64)
				{
					stream.emplace_back(buffer);

					buffer = 0;
					in_buffer_pos = 0;
				}
			}
		}

		size_t capacity() const {
			return 64 * stream.capacity() + 64;
		}

		size_t size() const {
			return 64 * stream.size() + in_buffer_pos;
		}

		//return number of bits that may be added without memory reallocation
		//such that memory reallocation is not needed even for flush_buffer
		//the assumption is that caller wants to avoid these reallocations
		//(for example when mem reallocs are not possible in the underlying storage)
		size_t free_bits() const {
			assert(stream.capacity() > stream.size() || in_buffer_pos == 0);

			return (stream.capacity() - stream.size()) * 64 - in_buffer_pos;
		}

		void serialize_fast(STORAGE_T<uint64_t>& vec) {
			flush_buffer();
			vec = std::move(stream);

			stream.clear();
		}

		void serialize(STORAGE_T<uint64_t>& vec) {
			flush_buffer();
			vec = stream;
		}

		void serialize(STORAGE_T<uint8_t>& vec) {
			vec.clear();
			vec.reserve(8 * stream.size() + 8);

			for (auto x : stream)
			{
				vec.emplace_back((uint8_t)(x & 0xffull));		x >>= 8;
				vec.emplace_back((uint8_t)(x & 0xffull));		x >>= 8;
				vec.emplace_back((uint8_t)(x & 0xffull));		x >>= 8;
				vec.emplace_back((uint8_t)(x & 0xffull));		x >>= 8;
				vec.emplace_back((uint8_t)(x & 0xffull));		x >>= 8;
				vec.emplace_back((uint8_t)(x & 0xffull));		x >>= 8;
				vec.emplace_back((uint8_t)(x & 0xffull));		x >>= 8;
				vec.emplace_back((uint8_t)(x & 0xffull));
			}

			uint64_t x = buffer;

			for (size_t i = 0; i < in_buffer_pos; i += 8)
			{
				vec.emplace_back((uint8_t)(x & 0xffull));		x >>= 8;
			}
		}

		//!!! may be unsafe, use with cautious
		//! invalidates the object, caller is responsible to call clear() after consuming returned memory
		void serialize_inplace(uint8_t*& data, size_t& size) {
			// consider checking hardware endianness and use some fast functions/intrinsic to swap bytes
			size_t full_elems = stream.size();
			if (in_buffer_pos)
				stream.emplace_back(); // make room for bytes from buffer

			data = reinterpret_cast<uint8_t*>(stream.data());

			size = 0;
			for (size_t i = 0; i < full_elems; ++i)
			{
				auto x = stream[i]; //important to make a copy because the content of the same memory is altered below!
				data[size++] = (uint8_t)(x & 0xffull);	x >>= 8;
				data[size++] = (uint8_t)(x & 0xffull);	x >>= 8;
				data[size++] = (uint8_t)(x & 0xffull);	x >>= 8;
				data[size++] = (uint8_t)(x & 0xffull);	x >>= 8;
				data[size++] = (uint8_t)(x & 0xffull);	x >>= 8;
				data[size++] = (uint8_t)(x & 0xffull);	x >>= 8;
				data[size++] = (uint8_t)(x & 0xffull);	x >>= 8;
				data[size++] = (uint8_t)(x & 0xffull);
			}

			uint64_t x = buffer;
			for (size_t i = 0; i < in_buffer_pos; i += 8)
			{
				data[size++] = (uint8_t)(x & 0xffull);	x >>= 8;
			}

			assert(full_elems * 8 + (in_buffer_pos + 7) / 8 == size);
			assert((uint8_t*)data == (uint8_t*)stream.data());
		}

		void serialize(uint8_t* vec) {
			for (auto x : stream)
			{
				*vec++ = (uint8_t)(x & 0xffull);		x >>= 8;
				*vec++ = (uint8_t)(x & 0xffull);		x >>= 8;
				*vec++ = (uint8_t)(x & 0xffull);		x >>= 8;
				*vec++ = (uint8_t)(x & 0xffull);		x >>= 8;
				*vec++ = (uint8_t)(x & 0xffull);		x >>= 8;
				*vec++ = (uint8_t)(x & 0xffull);		x >>= 8;
				*vec++ = (uint8_t)(x & 0xffull);		x >>= 8;
				*vec++ = (uint8_t)(x & 0xffull);
			}

			uint64_t x = buffer;

			for (size_t i = 0; i < in_buffer_pos; i += 8)
			{
				*vec++ = (uint8_t)(x & 0xffull);		x >>= 8;
			}
		}
	};

	// **************************************************************
	template<template<typename, typename...> typename STORAGE_T = std::vector>
	class bit_memory_reader : public basic_bit_memory<STORAGE_T> {
		using basic_bit_memory<STORAGE_T>::buffer;
		using basic_bit_memory<STORAGE_T>::stream;
		using basic_bit_memory<STORAGE_T>::in_buffer_pos;
		const uint64_t masks[65] =
		{ 0,
			0x1ull, 0x3ull, 0x7ull, 0xfull,
			0x1full, 0x3full, 0x7full, 0xffull,
			0x1ffull, 0x3ffull, 0x7ffull, 0xfffull,
			0x1fffull, 0x3fffull, 0x7fffull, 0xffffull,
			0x1ffffull, 0x3ffffull, 0x7ffffull, 0xfffffull,
			0x1fffffull, 0x3fffffull, 0x7fffffull, 0xffffffull,
			0x1ffffffull, 0x3ffffffull, 0x7ffffffull, 0xfffffffull,
			0x1fffffffull, 0x3fffffffull, 0x7fffffffull, 0xffffffffull,
			0x1ffffffffull, 0x3ffffffffull, 0x7ffffffffull, 0xfffffffffull,
			0x1fffffffffull, 0x3fffffffffull, 0x7fffffffffull, 0xffffffffffull,
			0x1ffffffffffull, 0x3ffffffffffull, 0x7ffffffffffull, 0xfffffffffffull,
			0x1fffffffffffull, 0x3fffffffffffull, 0x7fffffffffffull, 0xffffffffffffull,
			0x1ffffffffffffull, 0x3ffffffffffffull, 0x7ffffffffffffull, 0xfffffffffffffull,
			0x1fffffffffffffull, 0x3fffffffffffffull, 0x7fffffffffffffull, 0xffffffffffffffull,
			0x1ffffffffffffffull, 0x3ffffffffffffffull, 0x7ffffffffffffffull, 0xfffffffffffffffull,
			0x1fffffffffffffffull, 0x3fffffffffffffffull, 0x7fffffffffffffffull, 0xffffffffffffffffull
		};

		size_t in_vec_pos{};

		uint64_t read_bits_imp(size_t n)
		{
			uint64_t x = buffer >> in_buffer_pos;
			size_t no_buf_bits = 64 - in_buffer_pos;

			in_buffer_pos += n;

			if (in_buffer_pos >= 64 && in_vec_pos < stream.size())
			{
				buffer = stream[in_vec_pos++];
				in_buffer_pos -= 64;

				if (no_buf_bits < 64)
					x += buffer << no_buf_bits;
			}

			x &= masks[n];

			return x;
		}

	public:
		bit_memory_reader() = default;
		bit_memory_reader(bit_memory_reader<STORAGE_T>&&) = default;
		bit_memory_reader& operator=(bit_memory_reader<STORAGE_T>&&) = default;

		bit_memory_reader(STORAGE_T<uint64_t>&& storage) : basic_bit_memory<STORAGE_T>(std::move(storage)) {
			buffer = stream[in_vec_pos++];
		}

		//mkokot_TODO: probably to be removed
		void set_storage(STORAGE_T<uint64_t>&& storage)
		{
			stream = std::move(storage);
		}

		void assign(const STORAGE_T<uint64_t>& vec)
		{
			stream = vec;
			in_vec_pos = 0;

			buffer = stream[in_vec_pos++];
		}

		void assign(STORAGE_T<uint64_t>&& vec)
		{
			stream = std::move(vec);
			in_vec_pos = 0;

			buffer = stream[in_vec_pos++];
		}

		void assign(const STORAGE_T<uint8_t>& vec)
		{
			size_t i;
			uint64_t x;

			in_vec_pos = 0;
			stream.clear();
			stream.reserve((vec.size() + 7) / 8);

			for (i = 0; i + 7 < vec.size();)
			{
				x = ((uint64_t)vec[i++]);
				x += ((uint64_t)vec[i++]) << 8;
				x += ((uint64_t)vec[i++]) << 16;
				x += ((uint64_t)vec[i++]) << 24;
				x += ((uint64_t)vec[i++]) << 32;
				x += ((uint64_t)vec[i++]) << 40;
				x += ((uint64_t)vec[i++]) << 48;
				x += ((uint64_t)vec[i++]) << 56;

				stream.emplace_back(x);
			}

			if (i < vec.size())
			{
				x = 0;

				for (size_t j = 0; i < vec.size(); j += 8)
					x += ((uint64_t)vec[i++]) << j;

				stream.emplace_back(x);
			}
			buffer = stream[in_vec_pos++];
		}

		void assign(uint8_t* p, size_t n)
		{
			size_t i;
			uint64_t x;

			in_vec_pos = 0;
			stream.clear();
			stream.reserve((n + 7) / 8);

			for (i = 0; i + 7 < n; i += 8)
			{
				x = ((uint64_t)*p++);
				x += ((uint64_t)*p++) << 8;
				x += ((uint64_t)*p++) << 16;
				x += ((uint64_t)*p++) << 24;
				x += ((uint64_t)*p++) << 32;
				x += ((uint64_t)*p++) << 40;
				x += ((uint64_t)*p++) << 48;
				x += ((uint64_t)*p++) << 56;

				stream.emplace_back(x);
			}

			if (i < n)
			{
				x = 0;

				for (size_t j = 0; i < n; ++i, j += 8)
					x += ((uint64_t)*p++) << j;

				stream.emplace_back(x);
			}

			buffer = stream[in_vec_pos++];
		}

		uint8_t read_byte() {
			return (uint8_t)read_bits_imp(8);
		}

		uint16_t read_2bytes() {
			return (uint16_t)read_bits_imp(16);
		}

		uint32_t read_4bytes() {
			return (uint32_t)read_bits_imp(32);
		}

		uint64_t read_8bytes() {
			return (uint64_t)read_bits_imp(64);
		}

		uint64_t read_bits(size_t n) {
			return (uint64_t)read_bits_imp(n);
		}

		size_t size() const {
			return 64 * stream.size();
		}

		void setpos(size_t n)
		{
			in_vec_pos = n / 64;
			in_buffer_pos = n % 64;
		}
	};
}

#endif