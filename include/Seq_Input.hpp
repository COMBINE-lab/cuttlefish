
#ifndef SEQUENCE_INPUT_HPP
#define SEQUENCE_INPUT_HPP



#include <string>
#include <vector>


// A class to pack the input sequences.
class Seq_Input
{
private:

    const std::vector<std::string> seq_paths_;  // Collection of paths to raw sequences.
    const std::vector<std::string> list_paths_; // Collection of paths to lists containing sequence file paths.
    const std::vector<std::string> dir_paths_;  // Collection of paths to directories containing sequence files.


public:

    // Constructs a collection of input sequences.
    Seq_Input(  const std::vector<std::string>& seqs,
                const std::vector<std::string>& lists,
                const std::vector<std::string>& dirs):
        seq_paths_(seqs),
        list_paths_(lists),
        dir_paths_(dirs)
    {}


    // Returns the collection of paths to raw sequences.
    const std::vector<std::string>& seq_paths() const
    {
        return seq_paths_;
    }


    // Returns the collection of paths to lists containing sequence file paths.
    const std::vector<std::string>& list_paths() const
    {
        return list_paths_;
    }


    // Returns the collection of paths to directories containing sequence files.
    const std::vector<std::string>& dir_paths() const
    {
        return dir_paths_;
    }


    // Returns whether the sequence collection is empty or not.
    bool empty() const
    {
        return seq_paths_.empty() && list_paths_.empty() && dir_paths_.empty();
    }
};



#endif
