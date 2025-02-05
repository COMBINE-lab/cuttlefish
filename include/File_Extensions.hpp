
#ifndef FILE_EXTENSIONS_HPP
#define FILE_EXTENSIONS_HPP



namespace cuttlefish
{
    // File extensions for the data structures and files output by the algorithm.
    namespace file_ext
    {
        constexpr char edges_ext[] = ".cf_E";
        constexpr char vertices_ext[] = ".cf_V";
        constexpr char hash_ext[] = ".cf_hf";
        constexpr char buckets_ext[] = ".cf_hb";
        constexpr char unipaths_ext[] = ".fa";
        constexpr char json_ext[] = ".json";
        constexpr char temp[] = ".cf_op";

        // For reference dBGs only:

        constexpr char gfa1_ext[] = ".gfa1";
        constexpr char gfa2_ext[] = ".gfa2";
        constexpr char seg_ext[] = ".cf_seg";
        constexpr char seq_ext[] = ".cf_seq";

        // For cuttlefish 3.

        constexpr char atlas_ext[] = "_Atlas";
        constexpr char edge_matrix_ext[] = "_E";
        constexpr char edge_blk_ext[] = ".blk";
        constexpr char lmtig_bucket_ext[] = "_lmtig";
        constexpr char compressed_diagonal_ext[] = "_D";
        constexpr char vertex_p_inf_bucket_ext[] = "_P_v";
        constexpr char edge_p_inf_bucket_ext[] = "_P_e";
        constexpr char unitig_coord_bucket_ext[] = "_U";
        constexpr char color_rel_bucket_ext[] = "_C_rel";


        // For k-mer index.

        static constexpr char idx_file_ext[] = ".idx";
        static constexpr char min_inst_file_ext[] = ".mins";
    }
}



#endif
