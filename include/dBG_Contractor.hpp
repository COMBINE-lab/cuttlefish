
#ifndef DBG_CONTRACTOR_HPP
#define DBG_CONTRACTOR_HPP



#include "Discontinuity_Graph.hpp"
#include "Path_Info.hpp"
#include "Ext_Mem_Bucket.hpp"
#include "Async_Logger_Wrapper.hpp"
#include "Output_Sink.hpp"
#include "Character_Buffer.hpp"
#include "Build_Params.hpp"
#include "Data_Logistics.hpp"
#include "utility.hpp"
#include "globals.hpp"

#include <cstdint>
#include <vector>


namespace cuttlefish
{

// =============================================================================
// Compacted de Bruijn graph constructor.
template <uint16_t k>
class dBG_Contractor
{
public:

    typedef Obj_Path_Info_Pair<Kmer<k>, k> vertex_path_info_t;  // A vertex and its path-info.
    typedef Obj_Path_Info_Pair<uni_idx_t, k> unitig_path_info_t;    // A locally-maximal unitig's index in its bucket and its path-info.

    typedef Ext_Mem_Bucket_Concurrent<vertex_path_info_t> p_v_bucket_t;
    typedef Ext_Mem_Bucket_Concurrent<unitig_path_info_t> p_e_bucket_t;

    typedef std::vector<Padded<p_v_bucket_t>> P_v_t;
    typedef std::vector<Padded<p_e_bucket_t>> P_e_t;

    typedef Async_Logger_Wrapper sink_t;
    typedef Character_Buffer<sink_t> op_buf_t;
    typedef std::vector<Padded<op_buf_t>> op_buf_list_t;

private:

    const Build_Params params;  // Required parameters (wrapped inside).
    const Data_Logistics logistics; // Data logistics manager for the algorithm execution.

    P_v_t P_v;  // `P_v[j]` contains path-info for vertices in partition `j`.
    P_e_t P_e;  // `P_e[b]` contains path-info for edges induced by unitigs in bucket `b`.

    Output_Sink<sink_t> output_sink;    // Sink for the output maximal unitigs.

    // 128 KB (soft limit) worth of maximal unitig records (FASTA) can be retained in memory per worker, at most, before flushes.
    op_buf_list_t op_buf;   // Worker-specific output buffers.


    // Contracts the compacted de Bruijn graph from the parameters provided in
    // the constructor. `Colored_` determines whether to color the compacted
    // graph or not.
    template <bool Colored_> void construct();

    // Opens the containers for path-info of vertices.
    void open_p_v();

    // Releases the containers of path-info of vertices.
    void release_p_v();

    // Opens the containers for path-info of edges.
    void open_p_e();

    // Releases the containers of path-info of edges.
    void release_p_e();


public:

    dBG_Contractor(const dBG_Contractor&) = delete;
    dBG_Contractor& operator=(const dBG_Contractor&) = delete;

    // Constructs a compacted de Bruijn graph constructor with the parameters
    // required for the construction wrapped in `params`.
    dBG_Contractor(const Build_Params& params);

    // Contracts the compacted de Bruijn graph from the parameters provided in
    // the constructor.
    void construct();
};

}



#endif
