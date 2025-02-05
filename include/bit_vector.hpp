#pragma once

#include <cassert>
#include <vector>

#include "fulgor_util.hpp"

namespace fulgor {

struct bit_vector_builder {
    bit_vector_builder() { clear(); }

    void clear() {
        m_num_bits = 0;
        m_bits.clear();
        m_cur_word = nullptr;
    }

    uint64_t num_bits() const { return m_num_bits; }

    void resize(uint64_t num_bits) {
        m_num_bits = num_bits;
        m_bits.resize(util::num_64bit_words_for(num_bits), 0);
    }

    void reserve(uint64_t num_bits) { m_bits.reserve(util::num_64bit_words_for(num_bits)); }

    inline void set(uint64_t pos, bool b = true) {
        assert(pos < num_bits());
        uint64_t word = pos >> 6;
        uint64_t pos_in_word = pos & 63;
        m_bits[word] &= ~(uint64_t(1) << pos_in_word);
        m_bits[word] |= uint64_t(b) << pos_in_word;
    }

    void append_bits(uint64_t x, uint64_t len) {
        assert(len <= 64);
        assert(len == 64 or (x >> len) == 0);  // no other bits must be set
        if (len == 0) return;
        uint64_t pos_in_word = m_num_bits % 64;
        m_num_bits += len;
        if (pos_in_word == 0) {
            m_bits.push_back(x);
        } else {
            *m_cur_word |= x << pos_in_word;
            if (len > 64 - pos_in_word) { m_bits.push_back(x >> (64 - pos_in_word)); }
        }
        m_cur_word = &m_bits.back();
    }

    void append(bit_vector_builder const& bvb) {
        if (!bvb.num_bits()) return;
        uint64_t pos = m_bits.size();
        uint64_t shift = num_bits() % 64;
        m_num_bits = num_bits() + bvb.num_bits();
        m_bits.resize(util::num_64bit_words_for(m_num_bits));

        if (shift == 0) {  // word-aligned, easy case
            std::copy(bvb.m_bits.begin(), bvb.m_bits.end(), m_bits.begin() + ptrdiff_t(pos));
        } else {
            uint64_t* cur_word = &m_bits.front() + pos - 1;
            for (size_t i = 0; i < bvb.m_bits.size() - 1; ++i) {
                uint64_t w = bvb.m_bits[i];
                *cur_word |= w << shift;
                *++cur_word = w >> (64 - shift);
            }
            *cur_word |= bvb.m_bits.back() << shift;
            if (cur_word < &m_bits.back()) { *++cur_word = bvb.m_bits.back() >> (64 - shift); }
        }
        m_cur_word = &m_bits.back();
    }

    uint64_t const* data() const { return m_bits.data(); }
    std::vector<uint64_t>& bits() { return m_bits; }

private:
    uint64_t m_num_bits;
    std::vector<uint64_t> m_bits;
    uint64_t* m_cur_word;
};

struct bit_vector_iterator {
    bit_vector_iterator() : m_data(nullptr), m_num_64bit_words(0), m_pos(0), m_buf(0), m_avail(0) {}

    bit_vector_iterator(uint64_t const* data, uint64_t num_64bit_words, uint64_t pos = 0)
        : m_data(data), m_num_64bit_words(num_64bit_words) {
        at(pos);
    }

    void at(uint64_t pos) {
        m_pos = pos;
        m_buf = 0;
        m_avail = 0;
    }
    void at_and_clear_low_bits(uint64_t pos) {
        m_pos = pos;
        m_buf = m_data[pos / 64];
        m_buf &= uint64_t(-1) << (m_pos & 63);  // clear low bits
    }

    /* return the next l bits from the current position and advance by l bits */
    inline uint64_t take(uint64_t l) {
        assert(l <= 64);
        if (m_avail < l) fill_buf();
        uint64_t val;
        if (l != 64) {
            val = m_buf & ((uint64_t(1) << l) - 1);
            m_buf >>= l;
        } else {
            val = m_buf;
        }
        m_avail -= l;
        m_pos += l;
        return val;
    }

    uint64_t next() {
        unsigned long pos_in_word;
        uint64_t buf = m_buf;
        while (!util::lsb(buf, pos_in_word)) {
            m_pos += 64;
            buf = m_data[m_pos >> 6];
        }
        m_buf = buf & (buf - 1);  // clear LSB
        m_pos = (m_pos & ~uint64_t(63)) + pos_in_word;
        return m_pos;
    }

    /* skip all zeros from the current position and
    return the number of skipped zeros */
    inline uint64_t skip_zeros() {
        uint64_t zeros = 0;
        while (m_buf == 0) {
            m_pos += m_avail;
            zeros += m_avail;
            fill_buf();
        }
        uint64_t l = util::lsbll(m_buf);
        m_buf >>= l;
        m_buf >>= 1;
        m_avail -= l + 1;
        m_pos += l + 1;
        return zeros + l;
    }

    inline uint64_t position() const { return m_pos; }

    inline void fill_buf() {
        m_buf = get_next_word64();
        m_avail = 64;
    }

private:
    uint64_t get_next_word64() {
        uint64_t block = m_pos / 64;
        uint64_t shift = m_pos % 64;
        uint64_t word = m_data[block] >> shift;
        if (shift && block + 1 < m_num_64bit_words) { word |= m_data[block + 1] << (64 - shift); }
        return word;
    }

    uint64_t const* m_data;
    uint64_t m_num_64bit_words;
    uint64_t m_pos;
    uint64_t m_buf;
    uint64_t m_avail;
};

}  // namespace fulgor
