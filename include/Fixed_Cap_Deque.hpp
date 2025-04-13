
#ifndef FIXED_CAP_DEQUE_HPP
#define FIXED_CAP_DEQUE_HPP



#include "utility.hpp"

#include <cstddef>


namespace cuttlefish
{

// Fixed-capacity deque of capacity at least `req_cap_`. UB when size exceeds
// capacity.
template <typename T_, std::size_t req_cap>
class Fixed_Cap_Deque
{
private:

    static constexpr std::size_t cap_ = ceil_pow_2(req_cap);    // Maximum capacity.
    static constexpr auto wrap_mask = cap_ - 1; // Mask to wrap indices around the deque.

    std::size_t front_; // Current index to the front.
    std::size_t back_;  // Current index to the back.
    T_ arr[cap_];   // Underlying memory pool.

    void grow_back() { back_ = (back_ + 1) & wrap_mask; }

    void grow_front() { front_ = (front_ - 1) & wrap_mask; }    // Under(over)-flow is defined for unsigned.

    void shrink_back() { back_ = (back_ - 1) & wrap_mask; }

    void shrink_front() { front_ = (front_ + 1) & wrap_mask; }

public:

    Fixed_Cap_Deque():
          front_(0)
        , back_(0)
    {}

    bool empty() const { return front_ == back_; }

    const T_& front() const { return arr[front_]; }

    const T_& back() const { return arr[(back_ - 1) & wrap_mask]; }

    void clear() { front_ = back_; }

    void push_back(const T_& val) { arr[back_] = val; grow_back(); }

    void push_front(const T_& val) { grow_front(); arr[front_] = val; }

    template <typename... Args>
    void emplace_back(Args&&... args) { new (arr + back_) T_(args...); grow_back(); }

    template <typename... Args>
    void emplace_front(Args&&... args) { grow_front(); new (arr + front_) T_(args...); }

    void pop_back() { shrink_back(); }

    void pop_front() { shrink_front(); }
};

}



#endif
