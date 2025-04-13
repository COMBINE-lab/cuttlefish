
#ifndef MIN_HEAP_HPP
#define MIN_HEAP_HPP



#include <vector>
#include <algorithm>


// =============================================================================


// A class implementing a Min Heap of `T_elem_` type of elements.
template <typename T_elem_>
class Min_Heap
{
private:

    std::vector<T_elem_> container; // The collection of elements in the heap.


public:

    // Constructs an empty min heap.
    Min_Heap();

    // Constructs a min heap with the provided elements at `elems`.
    Min_Heap(const std::vector<T_elem_>& elems);

    // Constructs a min heap with the provided elements at `elems`, emptying
    // `elems`.
    Min_Heap(std::vector<T_elem_>& elems);

    // Destructs the min heap.
    ~Min_Heap();

    // Initializes the heap with the elements from `elems`, emptying `elems`.
    void init_heap(std::vector<T_elem_>& elems);

    // Returns the minimum element of the heap, i.e. the top element.
    const T_elem_& top() const;

    // Returns `true` iff the heap is empty.
    bool empty() const;

    // Returns the size of the heap.
    std::size_t size() const;

    // Pushes the element `elem` to the heap.
    void push(const T_elem_& elem);

    // Pushes all the elements from `elems` to the heap.
    void push(const std::vector<T_elem_>& elems);

    // Removes the minimum element from the heap.
    void pop();
};


template <typename T_elem_>
inline Min_Heap<T_elem_>::Min_Heap()
{}


template <typename T_elem_>
inline Min_Heap<T_elem_>::Min_Heap(const std::vector<T_elem_>& elems):
    container(elems)
{
    std::make_heap(container.begin(), container.end(), std::greater<T_elem_>{});
}


template <typename T_elem_>
inline Min_Heap<T_elem_>::Min_Heap(std::vector<T_elem_>& elems)
{
    init_heap(elems);
}


template <typename T_elem_>
inline Min_Heap<T_elem_>::~Min_Heap()
{
    container.clear();
    container.shrink_to_fit();
}


template <typename T_elem_>
inline const T_elem_& Min_Heap<T_elem_>::top() const
{
    return container.front();
}


template <typename T_elem_>
inline bool Min_Heap<T_elem_>::empty() const
{
    return container.empty();
}


template <typename T_elem_>
inline std::size_t Min_Heap<T_elem_>::size() const
{
    return container.size();
}


template <typename T_elem_>
inline void Min_Heap<T_elem_>::init_heap(std::vector<T_elem_>& elems)
{
    container.swap(elems);
    std::make_heap(container.begin(), container.end(), std::greater<T_elem_>{});
}


template <typename T_elem_>
inline void Min_Heap<T_elem_>::push(const T_elem_& elem)
{
    container.emplace_back(elem);
    std::push_heap(container.begin(), container.end(), std::greater<T_elem_>{});
}


template <typename T_elem_>
inline void Min_Heap<T_elem_>::push(const std::vector<T_elem_>& elems)
{
    container.reserve(container.size() + elems.size());
    std::for_each(elems.begin(), elems.end(), [&](const T_elem_& elem) { push(elem); } );
}


template <typename T_elem_>
inline void Min_Heap<T_elem_>::pop()
{
    std::pop_heap(container.begin(), container.end(), std::greater<T_elem_>{});
    container.pop_back();
}



#endif
