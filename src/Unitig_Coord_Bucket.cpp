
#include "Unitig_Coord_Bucket.hpp"
#include "globals.hpp"
#include "utility.hpp"

#include <cstdlib>
#include <algorithm>
#include <cassert>


namespace cuttlefish
{

template <uint16_t k>
Unitig_Coord_Bucket<k>::Unitig_Coord_Bucket(const std::string& path_pref):
      path_pref(path_pref)
    , coord_bucket(path_pref + ".coord", 8 * 1024)
    , label_bucket(path_pref + ".label", 8 * 1024)
    , size_(0)
    , label_len_(0)
{}


template <uint16_t k>
std::size_t Unitig_Coord_Bucket<k>::load_coords(Unitig_Coord<k, false>* const buf) const
{
    const auto b_sz = coord_bucket.load(buf);
    assert(b_sz == size_);

    return b_sz;
}


template <uint16_t k>
std::size_t Unitig_Coord_Bucket<k>::load_labels(char* const buf) const
{
    const auto len = label_bucket.load(buf);
    assert(len == label_len_);

    return len;
}


template <uint16_t k>
void Unitig_Coord_Bucket<k>::remove()
{
    coord_bucket.remove();
    label_bucket.remove();
}


template <uint16_t k, bool Colored_>
Unitig_Coord_Bucket_Concurrent<k, Colored_>::Unitig_Coord_Bucket_Concurrent(const std::string& path_pref):
      path_pref(path_pref)
    , flushed(0)
    , flushed_len(0)
    , flushed_color_c(0)
    , worker_buf(parlay::num_workers())
    , coord_os(coord_bucket_path(), std::ios::out | std::ios::binary)
    , label_os(label_bucket_path(), std::ios::out | std::ios::binary)
{
    if constexpr(Colored_)
        color_os.open(color_bucket_path(), std::ios::out | std::ios::binary);

    std::for_each(worker_buf.begin(), worker_buf.end(),
        [](auto& w_buf)
        {
            w_buf.unwrap().coord_buf.reserve(buf_sz_th / sizeof(Unitig_Coord<k, Colored_>));
            w_buf.unwrap().label_buf.reserve(buf_sz_th);
            if constexpr(Colored_)
                w_buf.unwrap().color_buf.reserve(buf_sz_th / sizeof(uint64_t));
        });
}


template <uint16_t k, bool Colored_>
Unitig_Coord_Bucket_Concurrent<k, Colored_>::Unitig_Coord_Bucket_Concurrent(Unitig_Coord_Bucket_Concurrent&& rhs):
      path_pref(std::move(rhs.path_pref))
    , flushed(std::move(rhs.flushed))
    , flushed_len(std::move(rhs.flushed_len))
    , flushed_color_c(std::move(rhs.flushed_color_c))
    , worker_buf(std::move(rhs.worker_buf))
    , coord_os(std::move(rhs.coord_os))
    , label_os(std::move(rhs.label_os))
    , color_os(std::move(rhs.color_os))
{}


template <uint16_t k, bool Colored_>
std::size_t Unitig_Coord_Bucket_Concurrent<k, Colored_>::size() const
{
    auto sz = flushed;
    std::for_each(worker_buf.cbegin(), worker_buf.cend(), [&](const auto& w_buf){ sz += w_buf.unwrap().coord_buf.size(); });

    return sz;
}


template <uint16_t k, bool Colored_>
std::size_t Unitig_Coord_Bucket_Concurrent<k, Colored_>::label_len() const
{
    auto len = flushed_len;
    std::for_each(worker_buf.cbegin(), worker_buf.cend(), [&](const auto& w_buf){ len += w_buf.unwrap().label_buf.size(); });

    return len;
}


template <uint16_t k, bool Colored_>
std::size_t Unitig_Coord_Bucket_Concurrent<k, Colored_>::color_count() const
{
    auto c = flushed_color_c;
    std::for_each(worker_buf.cbegin(), worker_buf.cend(), [&](const auto& w_buf){ c += w_buf.unwrap().color_buf.size(); });

    return c;
}


template <uint16_t k, bool Colored_>
std::size_t Unitig_Coord_Bucket_Concurrent<k, Colored_>::load_coords(Unitig_Coord<k, Colored_>* const buf) const
{
    const auto file_sz = load_file(coord_bucket_path(), reinterpret_cast<char*>(buf));
    assert(file_sz == flushed * sizeof(Unitig_Coord<k, Colored_>));
    (void)file_sz;

    auto coords_read = flushed; // Count of coordinates.
    std::size_t label_off_correct = flushed_len;    // Offset-correction factor of label-lengths for the in-memory coordinates.
    std::size_t color_off_correct = flushed_color_c;    // Offset-correction factor of colors for the in-memory coordinates.
    for(std::size_t i = 0; i < parlay::num_workers(); ++i)
    {
        const auto& w_buf = worker_buf[i].unwrap();
        const auto& coord_buf = w_buf.coord_buf;
        if(CF_LIKELY(!coord_buf.empty()))   // Conditional to avoid UB on `nullptr` being passed to `memcpy`.
        {
            std::memcpy(reinterpret_cast<char*>(buf + coords_read), reinterpret_cast<const char*>(coord_buf.data()), coord_buf.size() * sizeof(Unitig_Coord<k, Colored_>));

            // Offset-correct the un-flushed coordinates.
            const auto end = buf + coords_read;
            std::for_each(end, end + coord_buf.size(),
            [&](auto& v)
            {
                v.label_idx_ += label_off_correct;
                if constexpr(Colored_)
                    v.color_idx_ += color_off_correct;
            });

            label_off_correct += w_buf.label_buf.size();
            if constexpr(Colored_)
                color_off_correct += w_buf.color_buf.size();

            coords_read += coord_buf.size();
        }
    };

    return coords_read;
}


template <uint16_t k, bool Colored_>
std::size_t Unitig_Coord_Bucket_Concurrent<k, Colored_>::load_labels(char* const buf) const
{
    auto len = load_file(label_bucket_path(), buf); // Length of the label data (i.e. dump-string of the bucket).
    assert(len == flushed_len);

    std::for_each(worker_buf.cbegin(), worker_buf.cend(),
        [&](const auto& w_buf)
        {
            const auto& label_buf = w_buf.unwrap().label_buf;
            if(CF_LIKELY(!label_buf.empty()))
                std::memcpy(buf + len, label_buf.data(), label_buf.size());

            len += label_buf.size();
        });

    return len;
}


template <uint16_t k, bool Colored_>
std::size_t Unitig_Coord_Bucket_Concurrent<k, Colored_>::load_colors(Unitig_Color* const buf) const
{
    const auto file_sz = load_file(color_bucket_path(), reinterpret_cast<char*>(buf));
    assert(file_sz == flushed_color_c * sizeof(Unitig_Color));
    (void)file_sz;

    auto colors_read = flushed_color_c;
    std::for_each(worker_buf.cbegin(), worker_buf.cend(),
    [&](const auto& w_buf)
    {
        const auto& color_buf = w_buf.unwrap().color_buf;
        if(CF_LIKELY(!color_buf.empty()))   // Conditional to avoid UB on `nullptr` being passed to `memcpy`.
            std::memcpy(reinterpret_cast<char*>(buf + colors_read), reinterpret_cast<const char*>(color_buf.data()), color_buf.size() * sizeof(Unitig_Color));

        colors_read += color_buf.size();
    });

    return colors_read;
}


template <uint16_t k, bool Colored_>
void Unitig_Coord_Bucket_Concurrent<k, Colored_>::remove()
{
    coord_os.close();
    label_os.close();
    if constexpr(Colored_)
        color_os.close();

    bool err = !coord_os || !label_os || !remove_file(coord_bucket_path()) || !remove_file(label_bucket_path());
    if constexpr(Colored_)
        err = err || !color_os || !remove_file(color_bucket_path());

    if(err)
    {
        std::cerr << "Error removing files at path-prefix " << path_pref << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    force_free(worker_buf);
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Unitig_Coord_Bucket)
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_PER_BOOL, cuttlefish::Unitig_Coord_Bucket_Concurrent)
