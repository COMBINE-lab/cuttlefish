
#ifndef DATA_LOGISTICS_HPP
#define DATA_LOGISTICS_HPP



#include "Build_Params.hpp"

#include <string>
#include <vector>


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
};



#endif
