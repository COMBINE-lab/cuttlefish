#ifndef _SERIALIZATION
#define _SERIALIZATION

#include <cassert>
#include <exception>
#include <type_traits>
#include <istream>
#include <ostream>

namespace refresh
{
	const uint32_t REFRESH_BUILD_SERIALIZATION = 1;
	
	namespace detail {

		
		template<size_t N_BYTES>
		struct UnsignedIntImpl {};

		template<>
		struct UnsignedIntImpl<1> {
			using type = uint8_t;
		};

		template<>
		struct UnsignedIntImpl<2> {
			using type = uint16_t;
		};

		template<>
		struct UnsignedIntImpl<4> {
			using type = uint32_t;
		};

		template<>
		struct UnsignedIntImpl<8> {
			using type = uint64_t;
		};

		template<size_t N_BYTES>
		using UnsignedInt = typename UnsignedIntImpl<N_BYTES>::type;


		//callback must be compatibile with void store_byte(uint8_t byte);
		template<typename T, typename STORE_BYTE_CALLBACK>
		void read_bytes_little_endian_impl(const T& value, const STORE_BYTE_CALLBACK& store_byte)
		{
			static_assert(std::is_unsigned<T>::value, "This only works for unsigned");
			for (uint32_t i = 0; i < sizeof(value); ++i)
			{
				uint8_t byte = value >> (8 * i);
				store_byte(byte);
			}
		}

		//callback must be compatibile with uint8_t read_byte();
		template<typename T, typename READ_BYTE_CALLBACK>
		void write_bytes_little_endian_impl(T& value, const READ_BYTE_CALLBACK& read_byte)
		{
			static_assert(std::is_unsigned<T>::value, "This only works for unsigned");
			value = T{};
			for (uint32_t i = 0; i < sizeof(value); ++i)
			{
				T byte = read_byte();
				value |= byte << (8 * i);
			}
		}

		template<typename UNSIGNED_T, typename BASE_T>
		UNSIGNED_T change_type(const BASE_T& val) {
			static_assert(sizeof(UNSIGNED_T) == sizeof(BASE_T));
			return *reinterpret_cast<const UNSIGNED_T*>(&val);
		}
	}

	//callback must be compatibile with void store_byte(uint8_t byte);
	template<typename T, typename STORE_BYTE_CALLBACK>
	void read_bytes_little_endian(const T& value, const STORE_BYTE_CALLBACK& store_byte)
	{
		if constexpr (std::is_unsigned_v<T>)
			detail::read_bytes_little_endian_impl(value, store_byte);
		else
		{
			using unsigned_type = detail::UnsignedInt<sizeof(T)>;
			auto converted_val = detail::change_type<unsigned_type>(value);
			detail::read_bytes_little_endian_impl(converted_val, store_byte);
		}
	}

	//callback must be compatibile with uint8_t read_byte();
	template<typename T, typename READ_BYTE_CALLBACK>
	void write_bytes_little_endian(T& value, const READ_BYTE_CALLBACK& read_byte)
	{
		if constexpr (std::is_unsigned_v<T>)
			detail::write_bytes_little_endian_impl(value, read_byte);
		else
		{
			using unsigned_type = detail::UnsignedInt<sizeof(T)>;
			unsigned_type res;
			detail::write_bytes_little_endian_impl(res, read_byte);
			value = detail::change_type<T>(res);
		}
	}
	
	template<typename T>
	void serialize_little_endian(const T& value, std::ostream& out)
	{
		read_bytes_little_endian(value, [&out](uint8_t byte) {
			out.write(reinterpret_cast<const char*>(&byte), 1);
		});
	}
	template<typename T>
	void load_little_endian(T& value, std::istream& in)
	{
		write_bytes_little_endian(value, [&in]() {
			uint8_t byte;
			in.read(reinterpret_cast<char*>(&byte), 1);
			return byte;
		});
	}

	template<typename T>
	void serialize_little_endian(const std::vector<T>& vec, std::ostream& out) {
		uint64_t size = vec.size();
		serialize_little_endian(size, out);
		for (auto& x: vec)
			serialize_little_endian(x, out);
	}
	template<typename T>
	void load_little_endian(std::vector<T>& vec, std::istream& in) {
		uint64_t size;
		load_little_endian(size, in);
		vec.resize(size);
		for (auto& x : vec)
			load_little_endian(x, in);
	}
} // namespace refresh

#endif // _SERIALIZATION
