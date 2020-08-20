
#ifndef REFERENCE_INPUT_HPP
#define REFERENCE_INPUT_HPP



#include <string>
#include <vector>


class Reference_Input
{
private:

    const std::vector<std::string> ref_paths_;  // Collection of paths to raw references.
    const std::vector<std::string> list_paths_; // Collection of paths to lists containing reference file paths.
    const std::vector<std::string> dir_paths_;  // Collection of paths to directories containing reference files.


public:

    // Constructs a collection of input references.
    Reference_Input(const std::vector<std::string>& refs,
                    const std::vector<std::string>& lists,
                    const std::vector<std::string>& dirs):
        ref_paths_(refs),
        list_paths_(lists),
        dir_paths_(dirs)
    {}


    // Returns the collection of paths to raw references.
    const std::vector<std::string>& ref_paths() const
    {
        return ref_paths_;
    }


    // Returns the collection of paths to lists containing reference file paths.
    const std::vector<std::string>& list_paths() const
    {
        return list_paths_;
    }


    // Returns the collection of paths to directories containing reference files.
    const std::vector<std::string>& dir_paths() const
    {
        return dir_paths_;
    }
};



#endif
