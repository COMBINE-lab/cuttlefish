#pragma once

#include "fulgor_util.hpp"
#include "bit_vector.hpp"
#include "integer_codes.hpp"

#include <cstdint>
#include <limits>

namespace fulgor {
    enum list_type { delta_gaps, bitmap, complement_delta_gaps, differential_list };
//struct hybrid {
    //static const bool meta_colored = false;
    //static const bool differential_colored = false;

    struct color_set_builder {
        color_set_builder() {}
        color_set_builder(uint64_t num_colors) { init(num_colors); }

        void init(uint64_t num_colors) {
            m_num_colors = num_colors;

            /* if list contains < sparse_set_threshold_size ints, code it with gaps+delta */
            m_sparse_set_threshold_size = 0.25 * m_num_colors;

            /* if list contains > very_dense_set_threshold_size ints, code it as a complementary set
               with gaps+delta */
            m_very_dense_set_threshold_size = 0.75 * m_num_colors;
            /* otherwise: code it as a bitmap of m_num_colors bits */

            /*
            std::cout << "m_num_colors: " << m_num_colors << std::endl;
            std::cout << "m_sparse_set_threshold_size " << m_sparse_set_threshold_size << std::endl;
            std::cout << "m_very_dense_set_threshold_size " << m_very_dense_set_threshold_size
                      << std::endl;

            m_bvb.reserve(8 * essentials::GB);
            m_offsets.push_back(0);

            m_num_lists = 0;
            m_num_total_integers = 0;
            */
        }

        void process(uint32_t const* colors, uint64_t list_size, bit_vector_builder& bvb_ints) {
            bvb_ints.clear();

            /* encode list_size */
            // constexpr uint64_t u64_limit = std::numeric_limits<uint64_t>::max() - 1;
            util::write_delta(bvb_ints, list_size);

            if (list_size < m_sparse_set_threshold_size) {
                uint32_t prev_val = colors[0];
                util::write_delta(bvb_ints, prev_val);
                for (uint64_t i = 1; i != list_size; ++i) {
                    uint32_t val = colors[i];
                    assert(val >= prev_val + 1);
                    util::write_delta(bvb_ints, val - (prev_val + 1));
                    prev_val = val;
                }
            } else if (list_size < m_very_dense_set_threshold_size) {
                // not sparse or dense --- encode as bitvector
                bit_vector_builder bvb_int_set;
                bvb_int_set.resize(m_num_colors);
                for (uint64_t i = 0; i != list_size; ++i) bvb_int_set.set(colors[i]);
                bvb_ints.append(bvb_int_set);
                //m_bvb.append(bvb_ints);
            } else {
                bool first = true;
                uint32_t val = 0;
                uint32_t prev_val = std::numeric_limits<uint32_t>::max();
                uint32_t written = 0;
                for (uint64_t i = 0; i != list_size; ++i) {
                    uint32_t x = colors[i];
                    while (val < x) {
                        if (first) {
                            util::write_delta(bvb_ints, val);
                            first = false;
                            ++written;
                        } else {
                            assert(val >= prev_val + 1);
                            util::write_delta(bvb_ints, val - (prev_val + 1));
                            ++written;
                        }
                        prev_val = val;
                        ++val;
                    }
                    assert(val == x);
                    val++;  // skip x
                }
                while (val < m_num_colors) {
                    assert(val >= prev_val + 1);
                    util::write_delta(bvb_ints, val - (prev_val + 1));
                    prev_val = val;
                    ++val;
                    ++written;
                }
                assert(val == m_num_colors);
                /* complementary_list_size = m_num_colors - list_size */
                assert(m_num_colors - list_size <= m_num_colors);
                assert(written == m_num_colors - list_size);
            }
            
            //m_offsets.push_back(m_bvb.num_bits());
            m_num_total_integers += list_size;
            //m_num_lists += 1;
            //if (m_num_lists % 500000 == 0) {
            //std::cout << "  processed " << m_num_lists << " lists" << std::endl;
            //}
        }

        /*
        void build(hybrid& h) {
            h.m_num_colors = m_num_colors;
            h.m_sparse_set_threshold_size = m_sparse_set_threshold_size;
            h.m_very_dense_set_threshold_size = m_very_dense_set_threshold_size;

            std::cout << "processed " << m_num_lists << " lists" << std::endl;
            std::cout << "m_num_total_integers " << m_num_total_integers << std::endl;
            assert(m_num_lists == m_offsets.size() - 1);

            h.m_offsets.encode(m_offsets.begin(), m_offsets.size(), m_offsets.back());
            h.m_colors.swap(m_bvb.bits());

            std::cout << "  total bits for ints = " << h.m_colors.size() * 64 << std::endl;
            std::cout << "  total bits per offsets = " << h.m_offsets.num_bits() << std::endl;
            std::cout << "  total bits = " << h.m_offsets.num_bits() + h.m_colors.size() * 64
                      << std::endl;
            std::cout << "  offsets: "
                      << static_cast<double>(h.m_offsets.num_bits()) / m_num_total_integers
                      << " bits/int" << std::endl;
            std::cout << "  lists: "
                      << static_cast<double>(h.m_colors.size() * 64) / m_num_total_integers
                      << " bits/int" << std::endl;
        }
        */

    private:
        uint32_t m_num_colors;
        uint32_t m_sparse_set_threshold_size;
        uint32_t m_very_dense_set_threshold_size;
        // uint64_t m_num_lists;
        uint64_t m_num_total_integers;

        //bit_vector_builder m_bvb;
        //std::vector<uint64_t> m_offsets;
    };
    
   
    //struct color_set_decoder {
    //    uint64_t m_num_colors;
    //    bit_vector_iterator m_it;

    //    color_set_decoder(bit_vector_iterator it, uint64_t num_colors) : m_it(it), m_num_colors(num_colors) {
    //    }

    //}

    struct forward_iterator {
        forward_iterator() {}

        forward_iterator(bit_vector_iterator it, uint64_t num_colors)
            : m_bitmap_begin(0)
            // , m_colors_begin(0)
            , m_num_colors(num_colors)
            , m_it(it) {
            m_sparse_set_threshold_size = 0.25 * m_num_colors;
            m_very_dense_set_threshold_size = 0.75 * m_num_colors;
            rewind();
        }

        void rewind() {
            m_pos_in_list = 0;
            m_pos_in_comp_list = 0;
            m_comp_list_size = 0;
            m_comp_val = -1;
            m_prev_val = -1;
            m_curr_val = 0;
            //m_it = bit_vector_iterator((m_ptr->m_colors).data(), (m_ptr->m_colors).size(),
            //                          m_colors_begin);
            m_size = util::read_delta(m_it);
            /* set m_type and read the first value */
            if (m_size < m_sparse_set_threshold_size) {
                //std::cerr << "m_size = " << m_size << " :: sparse set\n";
                m_type = list_type::delta_gaps;
                m_curr_val = util::read_delta(m_it);
            } else if (m_size < m_very_dense_set_threshold_size) {
                //std::cerr << "m_size = " << m_size << " :: middle set\n";
                m_type = list_type::bitmap;
                m_bitmap_begin = m_it.position();  // after m_size
                m_it.at_and_clear_low_bits(m_bitmap_begin);
                uint64_t pos = m_it.next();
                assert(pos >= m_bitmap_begin);
                m_curr_val = pos - m_bitmap_begin;
            } else {
                //std::cerr << "m_size = " << m_size << " :: dense set\n";
                m_type = list_type::complement_delta_gaps;
                m_comp_list_size = m_num_colors - m_size;
                if (m_comp_list_size > 0) m_comp_val = util::read_delta(m_it);
                next_comp_val();
            }
        }

        /* this is needed to annul the next_comp_val() done in the constructor
           if we want to iterate through the complemented set */
        void reinit_for_complemented_set_iteration() {
            assert(m_type == list_type::complement_delta_gaps);
            m_pos_in_comp_list = 0;
            m_prev_val = -1;
            m_curr_val = 0;
            m_it.at(0);//= bit_vector_iterator((m_ptr->m_colors).data(), (m_ptr->m_colors).size(),
                       //                m_colors_begin);
            util::read_delta(m_it); /* skip m_size */
            if (m_comp_list_size > 0) {
                m_comp_val = util::read_delta(m_it);
            } else {
                m_comp_val = m_num_colors;
            }
        }

        uint64_t value() const { return m_curr_val; }
        uint64_t comp_value() const { return m_comp_val; }
        uint64_t operator*() const { return value(); }

        void next() {
            if (m_type == list_type::complement_delta_gaps) {
                ++m_curr_val;
                if (m_curr_val >= m_num_colors) {  // saturate
                    m_curr_val = m_num_colors;
                    return;
                }
                next_comp_val();
            } else if (m_type == list_type::delta_gaps) {
                m_pos_in_list += 1;
                if (m_pos_in_list >= m_size) {  // saturate
                    m_curr_val = m_num_colors;
                    return;
                }
                m_prev_val = m_curr_val;
                m_curr_val = util::read_delta(m_it) + (m_prev_val + 1);
            } else {
                assert(m_type == list_type::bitmap);
                m_pos_in_list += 1;
                if (m_pos_in_list >= m_size) {  // saturate
                    m_curr_val = m_num_colors;
                    return;
                }
                uint64_t pos = m_it.next();
                assert(pos >= m_bitmap_begin);
                m_curr_val = pos - m_bitmap_begin;
            }
        }

        void next_comp() {
            ++m_pos_in_comp_list;
            if (m_pos_in_comp_list >= m_comp_list_size) {  // saturate
                m_comp_val = m_num_colors;
                return;
            }
            m_prev_val = m_comp_val;
            m_comp_val = util::read_delta(m_it) + (m_prev_val + 1);
        }

        void operator++() { next(); }

        /* update the state of the iterator to the element
           which is greater-than or equal-to lower_bound */
        void next_geq(const uint64_t lower_bound) {
            assert(lower_bound <= num_colors());
            if (m_type == list_type::complement_delta_gaps) {
                if (value() > lower_bound) return;
                next_geq_comp_val(lower_bound);
                m_curr_val = lower_bound + (m_comp_val == lower_bound);
            } else {
                while (value() < lower_bound) next();
            }
            assert(value() >= lower_bound);
        }

        uint32_t size() const { return m_size; }
        uint32_t num_colors() const { return m_num_colors; }
        int type() const { return m_type; }

    private:
        //hybrid const* m_ptr;
        uint64_t m_bitmap_begin;
        // uint64_t m_colors_begin;
        uint32_t m_num_colors;
        int m_type;
        uint32_t m_sparse_set_threshold_size;
        uint32_t m_very_dense_set_threshold_size;
 
        bit_vector_iterator m_it;
        uint32_t m_pos_in_list;
        uint32_t m_size;

        uint32_t m_pos_in_comp_list;
        uint32_t m_comp_list_size;

        uint32_t m_comp_val;
        uint32_t m_prev_val;
        uint32_t m_curr_val;

        void next_comp_val() {
            while (m_curr_val == m_comp_val) {
                ++m_curr_val;
                ++m_pos_in_comp_list;
                if (m_pos_in_comp_list >= m_comp_list_size) break;
                m_prev_val = m_comp_val;
                m_comp_val = util::read_delta(m_it) + (m_prev_val + 1);
            }
        }

        void next_geq_comp_val(const uint64_t lower_bound) {
            while (m_comp_val < lower_bound) {
                ++m_pos_in_comp_list;
                if (m_pos_in_comp_list >= m_comp_list_size) break;
                m_prev_val = m_comp_val;
                m_comp_val = util::read_delta(m_it) + (m_prev_val + 1);
            }
        }
    };
    

    //typedef forward_iterator iterator_type;

    //forward_iterator color_set(uint64_t color_set_id) const {
    //    assert(color_set_id < num_color_sets());
    //    uint64_t begin = m_offsets.access(color_set_id);
    //    return forward_iterator(this, begin);
    //}

    //uint32_t num_colors() const { return m_num_colors; }
    //uint64_t num_color_sets() const { return m_offsets.size() - 1; }

    //uint64_t num_bits() const {
    //    return (sizeof(m_num_colors) + sizeof(m_sparse_set_threshold_size) +
    //            sizeof(m_very_dense_set_threshold_size)) *
    //               8 +
    //           m_offsets.num_bits() + essentials::vec_bytes(m_colors) * 8;
    //}

    //void print_stats() const {
    //    uint64_t num_buckets = 10;
    //    assert(num_buckets > 0);
    //    uint64_t bucket_size = m_num_colors / num_buckets;
    //    std::vector<uint32_t> list_size_upperbounds;
    //    for (uint64_t i = 0, curr_list_size_upper_bound = bucket_size; i != num_buckets;
    //         ++i, curr_list_size_upper_bound += bucket_size) {
    //        if (i == num_buckets - 1) curr_list_size_upper_bound = m_num_colors;
    //        list_size_upperbounds.push_back(curr_list_size_upper_bound);
    //    }

    //    std::vector<uint64_t> num_bits_per_bucket;
    //    std::vector<uint64_t> num_lists_per_bucket;
    //    std::vector<uint64_t> num_ints_per_bucket;
    //    num_bits_per_bucket.resize(num_buckets, 0);
    //    num_lists_per_bucket.resize(num_buckets, 0);
    //    num_ints_per_bucket.resize(num_buckets, 0);

    //    const uint64_t num_lists = num_color_sets();
    //    uint64_t num_total_integers = 0;
    //    for (uint64_t color_set_id = 0; color_set_id != m_offsets.size() - 1; ++color_set_id) {
    //        uint64_t offset = m_offsets.access(color_set_id);
    //        bit_vector_iterator it(m_colors.data(), m_colors.size(), offset);
    //        uint32_t list_size = util::read_delta(it);
    //        uint64_t num_bits = m_offsets.access(color_set_id + 1) - offset;
    //        auto bucket_it = std::upper_bound(list_size_upperbounds.begin(),
    //                                          list_size_upperbounds.end(), list_size);
    //        if (bucket_it != list_size_upperbounds.begin() and *(bucket_it - 1) == list_size) {
    //            --bucket_it;
    //        }
    //        uint64_t bucket_index = std::distance(list_size_upperbounds.begin(), bucket_it);
    //        num_bits_per_bucket[bucket_index] += num_bits;
    //        num_lists_per_bucket[bucket_index] += 1;
    //        num_ints_per_bucket[bucket_index] += list_size;
    //        num_total_integers += list_size;
    //    }

    //    std::cout << "CCs SPACE BREAKDOWN:\n";
    //    uint64_t integers = 0;
    //    uint64_t bits = 0;
    //    const uint64_t total_bits = num_bits();
    //    for (uint64_t i = 0, curr_list_size_upper_bound = 0; i != num_buckets; ++i) {
    //        if (i == num_buckets - 1) {
    //            curr_list_size_upper_bound = m_num_colors;
    //        } else {
    //            curr_list_size_upper_bound += bucket_size;
    //        }
    //        if (num_lists_per_bucket[i] > 0) {
    //            uint64_t n = num_ints_per_bucket[i];
    //            integers += n;
    //            bits += num_bits_per_bucket[i];
    //            std::cout << "num. lists of size > " << (curr_list_size_upper_bound - bucket_size)
    //                      << " and <= " << curr_list_size_upper_bound << ": "
    //                      << num_lists_per_bucket[i] << " ("
    //                      << (num_lists_per_bucket[i] * 100.0) / num_lists
    //                      << "%) -- integers: " << n << " (" << (n * 100.0) / num_total_integers
    //                      << "%) -- bits/int: " << static_cast<double>(num_bits_per_bucket[i]) / n
    //                      << " -- "
    //                      << static_cast<double>(num_bits_per_bucket[i]) / total_bits * 100.0
    //                      << "\% of total space" << '\n';
    //        }
    //    }
    //    assert(integers == num_total_integers);
    //    assert(std::accumulate(num_lists_per_bucket.begin(), num_lists_per_bucket.end(),
    //                           uint64_t(0)) == num_lists);
    //    std::cout << "  colors: " << static_cast<double>(bits) / integers << " bits/int"
    //              << std::endl;
    //    std::cout << "  offsets: "
    //              << static_cast<double>((sizeof(m_num_colors) +
    //                                      sizeof(m_sparse_set_threshold_size) +
    //                                      sizeof(m_very_dense_set_threshold_size)) *
    //                                         8 +
    //                                     m_offsets.num_bits()) /
    //                     integers
    //              << " bits/int" << std::endl;
    //}

    //template <typename Visitor>
    //void visit(Visitor& visitor) {
    //    visit_impl(visitor, *this);
    //}

    //template <typename Visitor>
    //void visit(Visitor& visitor) const {
    //    visit_impl(visitor, *this);
    //}

/*private:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_num_colors);
        visitor.visit(t.m_sparse_set_threshold_size);
        visitor.visit(t.m_very_dense_set_threshold_size);
        visitor.visit(t.m_offsets);
        visitor.visit(t.m_colors);
    }

    uint32_t m_num_colors;
    uint32_t m_sparse_set_threshold_size;
    uint32_t m_very_dense_set_threshold_size;
    sshash::ef_sequence<false> m_offsets;
    std::vector<uint64_t> m_colors;
    */
//};

}  // namespace fulgor
