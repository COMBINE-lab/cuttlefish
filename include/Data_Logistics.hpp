
#ifndef DATA_LOGISTICS_HPP
#define DATA_LOGISTICS_HPP



#include <string>
#include <vector>


class Build_Params;


// =============================================================================
// A class to govern the logistical policies regarding the various data used—
// either as input, output, or temporary—during the lifetime of Cuttlefish.
class Data_Logistics
{
private:

    const Build_Params& params;    // The construction parameters passed to Cuttlefish.


public:

    // Constructs a logistics manager object for the parameters in `params`.
    Data_Logistics(const Build_Params& build_params);

    // Returns the collection of file paths that are input to Cuttlefish.
    const std::vector<std::string> input_paths_collection() const;

    // Returns the path prefix for temporary files used by Cuttlefish.
    const std::string working_dir_path() const;

    // Returns the path prefix to the edge database being used by Cuttlefish.
    const std::string edge_db_path() const;

    // Returns the path prefix to the vertex database being used by Cuttlefish.
    const std::string vertex_db_path() const;

    // Returns the path to the final output file by Cuttlefish.
    const std::string output_file_path() const;

    // Returns the directory to the atlases.
    const std::string atlas_path() const;

    // Returns the directory to the edge-matrix produced by Cuttlefish.
    const std::string edge_matrix_path() const;

    // Returns the directory to the buckets for lm-tigs produced by Cuttlefish.
    const std::string lmtig_buckets_path() const;

    // Returns the directory to the edges introduced in diagonal blocks
    // contraction by Cuttlefish.
    const std::string compressed_diagonal_path() const;

    // Returns the directory to the color-relationship buckets.
    const std::string color_rel_bucket_path() const;

    // Returns the directory of the buckets storing path-information of
    // vertices.
    const std::string vertex_path_info_buckets_path() const;

    // Returns the directory of the buckets storing path-information of edges.
    const std::string edge_path_info_buckets_path() const;

    // Returns path prefix to the unitig-coordinate buckets produced in map-
    // reduce by Cuttlefish.
    const std::string unitig_coord_buckets_path() const;
};



#endif
