
#include "Data_Logistics.hpp"
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
