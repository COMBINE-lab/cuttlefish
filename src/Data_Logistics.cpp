
#include "Data_Logistics.hpp"
#include "Build_Params.hpp"
#include "File_Extensions.hpp"
#include "utility.hpp"


Data_Logistics::Data_Logistics(const Build_Params& build_params):
    params(build_params)
{}


const std::vector<std::string> Data_Logistics::input_paths_collection() const
{
    return params.sequence_input().seqs();
}


const std::string Data_Logistics::working_dir_path() const
{
    return dirname(params.output_prefix());
}


const std::string Data_Logistics::edge_db_path() const
{
#ifdef CF_DEVELOP_MODE
    if(!params.edge_db_path().empty())
        return params.edge_db_path();
#endif

    return params.working_dir_path() + filename(params.output_prefix()) + cuttlefish::file_ext::edges_ext;
}


const std::string Data_Logistics::vertex_db_path() const
{
#ifdef CF_DEVELOP_MODE
    if(!params.vertex_db_path().empty())
        return params.vertex_db_path();
#endif

    return params.working_dir_path() + filename(params.output_prefix()) + cuttlefish::file_ext::vertices_ext;
}


const std::string Data_Logistics::output_file_path() const
{
    return params.output_file_path();
}


const std::string Data_Logistics::atlas_path() const
{
    return params.working_dir_path() + filename(params.output_prefix()) + cuttlefish::file_ext::atlas_ext;
}


const std::string Data_Logistics::edge_matrix_path() const
{
    return params.output_prefix() + cuttlefish::file_ext::edge_matrix_ext;
}


const std::string Data_Logistics::lmtig_buckets_path() const
{
    return params.output_prefix() + cuttlefish::file_ext::lmtig_bucket_ext;
}


const std::string Data_Logistics::compressed_diagonal_path() const
{
    return params.output_prefix() + cuttlefish::file_ext::compressed_diagonal_ext;
}


const std::string Data_Logistics::color_rel_bucket_path() const
{
    return params.working_dir_path() + filename(params.output_prefix()) + cuttlefish::file_ext::color_rel_bucket_ext;
}


const std::string Data_Logistics::vertex_path_info_buckets_path() const
{
    return params.working_dir_path() + filename(params.output_prefix()) + cuttlefish::file_ext::vertex_p_inf_bucket_ext;
}


const std::string Data_Logistics::edge_path_info_buckets_path() const
{
    return params.output_prefix() + cuttlefish::file_ext::edge_p_inf_bucket_ext;
}


const std::string Data_Logistics::unitig_coord_buckets_path() const
{
    return params.working_dir_path() + filename(params.output_prefix()) + cuttlefish::file_ext::unitig_coord_bucket_ext;
}
