#ifndef _MEMORY_CHUNK
#define _MEMORY_CHUNK

#include <cassert>
#include <exception>
#include <iterator>
#include <stdexcept>

namespace refresh
{
	const uint32_t REFRESH_BUILD_MEMORY_CHUNK = 1;
	template<typename T>
	class memory_chunk
	{
	public:
		using size_type = std::size_t;
		using value_type = T;
		using reference = value_type&;
		using const_reference = const value_type&;
		using difference_type = std::ptrdiff_t;
		using iterator = value_type*;
		using const_iterator = const value_type*;
		using reverse_iterator = std::reverse_iterator<iterator>;
		using const_reverse_iterator = std::reverse_iterator<const_iterator>;
		using pointer = value_type*;
		using const_pointer = const value_type*;

	private:
		T* memory{};
		size_type _size{};
		size_type _capacity{};
	public:
		memory_chunk() = default; //default memory view is view of nothing
		memory_chunk(const memory_chunk& rhs) = delete;
		memory_chunk<T>& operator=(const memory_chunk<T>& rhs) = delete;
		memory_chunk(memory_chunk&& rhs) noexcept :
			memory(rhs.memory),
			_size(rhs._size),
			_capacity(rhs._capacity)
		{
			rhs._size = rhs._capacity = 0;
			rhs.memory = nullptr;
		}

		//kinda unsafe, copy ctor and assignmend are disabled to avoid creating copy by incident, in some cases it is convinient to have a copy, but this is a caller responsibility to manage it appropriatelly
		memory_chunk<T> get_copy() const
		{
			auto copy = memory_chunk<T>(memory, _capacity);
			copy._size = _size;
			return copy;
		}
		memory_chunk<T>& operator=(memory_chunk<T>&& rhs) noexcept
		{
			if (this == &rhs)
				return *this;
			memory = rhs.memory;
			_size = rhs._size;
			_capacity = rhs._capacity;

			rhs._size = rhs._capacity = 0;
			rhs.memory = nullptr;

			return *this;
		}

		memory_chunk(T* memory, size_type capacity) :
			memory(memory),
			_capacity(capacity)
		{

		}
		void clear()
		{
			_size = 0;
		}

		void resize(size_type new_size)
		{
			assert(new_size <= _capacity);
			_size = new_size;
		}

		//does nothing
		void reserve(size_type new_capacity)
		{
			assert(new_capacity <= _capacity);
		}

		size_type size() const
		{
			return _size;
		}

		size_type capacity() const
		{
			return _capacity;
		}

		size_type max_size() const
		{
			return _capacity;
		}

		void change_capacity(size_type new_capacity)
		{
			assert(new_capacity <= _capacity);
			_capacity = new_capacity;
		}

		reference operator[](size_type pos)
		{
			assert(pos < _size);
			return memory[pos];
		}

		const_reference operator[](size_type pos) const
		{
			assert(pos < _size);
			return memory[pos];
		}

		reference at(size_type pos)
		{
			if (pos >= size)
			{
				throw std::out_of_range("invalid index");
			}
			return memory[pos];
		}

		const_reference at(size_type pos) const
		{
			if (pos >= size)
			{
				throw std::out_of_range("invalid index");
			}
			return memory[pos];
		}

		void push_back(const T& val)
		{
			assert(_size < _capacity);
			emplace_back(val);
		}

		void push_back(T&& val)
		{
			assert(_size < _capacity);
			emplace_back(std::move(val));
		}

		template< class... Args >
		reference emplace_back(Args&&... args)
		{
			assert(_size < _capacity);
			return *new (memory + _size++) T(std::forward<Args>(args)...);
		}
		void pop_back()
		{
			assert(_size);
			--size;
		}

		bool empty() const
		{
			return !_size;
		}

		reference front()
		{
			assert(_size);
			return memory[0];
		}

		reference back()
		{
			assert(_size);
			return memory[_size - 1];
		}

		const_reference front() const
		{
			assert(_size);
			return memory[0];
		}

		const_reference back() const
		{
			assert(_size);
			return memory[_size - 1];
		}

		iterator begin()
		{
			return memory;
		}

		iterator end()
		{
			return memory + _size;
		}

		const_iterator begin() const
		{
			return memory;
		}

		const_iterator end() const
		{
			return memory + _size;
		}

		const_iterator cbegin()
		{
			return memory;
		}

		const_iterator cend()
		{
			return memory + size;
		}

		const_iterator cbegin() const
		{
			return memory;
		}

		const_iterator cend() const
		{
			return memory + size;
		}

		reverse_iterator rbegin()
		{
			return reverse_iterator(memory + _size);
		}

		reverse_iterator rend()
		{
			return reverse_iterator(memory);
		}

		const_reverse_iterator rbegin() const
		{
			return const_reverse_iterator(memory + _size);
		}

		const_reverse_iterator rend() const
		{
			return const_reverse_iterator(memory);
		}

		const_reverse_iterator crbegin()
		{
			return const_reverse_iterator(memory + _size);
		}

		const_reverse_iterator crend()
		{
			return const_reverse_iterator(memory);
		}

		const_reverse_iterator crbegin() const
		{
			return const_reverse_iterator(memory + _size);
		}

		const_reverse_iterator crend() const
		{
			return const_reverse_iterator(memory);
		}

		pointer data()
		{
			return memory;
		}

		const_pointer data() const
		{
			return memory;
		}

		//returns end fragment of size rhs.capacity() - new_capacity
		//rhs becomes smaller: have new_capacity as capacity()
		//size() must <= new_capacity
		memory_chunk<T> split(size_t new_capacity)
		{
			assert(size() <= new_capacity);
			memory_chunk<T> res(data() + new_capacity, capacity() - new_capacity);
			change_capacity(new_capacity);
			return res;
		}
	};

	// Get new representation 
	// It is the caller responsibility to assure that casting TYPE_FROM* -> TYPE_TO* makes sense
	// typical use case: have object over some complex type and want to have it as a simpler type (for example as uint8_t) to keep simpler template version
	// sizeof(U) % sizeof(T) must = 0 OR sizeof(T) % sizeof(U) must = 0
	// Instead of this function there could be a template move ctor and move assignment, but I think it is to risky to allow this under the hood. Here caller must explicitly call this function
	template<typename TYPE_TO, typename TYPE_FROM>
	memory_chunk<TYPE_TO> convert_type(memory_chunk<TYPE_FROM>&& rhs)
	{
		decltype(rhs.size()) size;
		decltype(rhs.capacity()) capacity;

		if constexpr (sizeof(TYPE_TO) > sizeof(TYPE_FROM))
		{
			static_assert(sizeof(TYPE_TO) % sizeof(TYPE_FROM) == 0, "sizeof(TYPE_TO) % sizeof(TYPE_FROM) != 0");
			constexpr size_t X = sizeof(TYPE_TO) / sizeof(TYPE_FROM);

			size = rhs.size() / X;
			capacity = rhs.capacity() / X;
		}
		else
		{
			static_assert(sizeof(TYPE_FROM) % sizeof(TYPE_TO) == 0, "sizeof(TYPE_FROM) % sizeof(TYPE_TO) != 0");
			constexpr size_t X = sizeof(TYPE_FROM) / sizeof(TYPE_TO);
			size = rhs.size() * X;
			capacity = rhs.capacity() * X;
		}
		memory_chunk<TYPE_TO> res(reinterpret_cast<TYPE_TO*>(rhs.data()), capacity);
		res.resize(size);
		rhs = memory_chunk<TYPE_FROM>();
		return res;
	}

	template<typename TYPE_TO, typename TYPE_FROM>
	memory_chunk<TYPE_TO> convert_type(memory_chunk<TYPE_FROM>&& rhs, size_t new_capacity, size_t new_size = 0)
	{
		assert(sizeof(TYPE_TO) * new_capacity <= sizeof(TYPE_FROM) * rhs.capacity());
		assert(new_size <= new_capacity);

		memory_chunk<TYPE_TO> res(reinterpret_cast<TYPE_TO*>(rhs.data()), new_capacity);
		res.resize(new_size);
		rhs = memory_chunk<TYPE_FROM>();
		return res;
	}

} // namespace refresh

#endif // _MEMORY_CHUNK
