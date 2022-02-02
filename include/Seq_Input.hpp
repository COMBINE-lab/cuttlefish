
#ifndef SEQ_INPUT_HPP
#define SEQ_INPUT_HPP



#include <string>
#include <vector>
#include <optional>


// A class to pack the input sequences.
class Seq_Input
{
private:

    const std::vector<std::string> seq_paths_;  // Collection of paths to raw sequences.
    const std::vector<std::string> list_paths_; // Collection of paths to lists containing sequence file paths.
    const std::vector<std::string> dir_paths_;  // Collection of paths to directories containing sequence files.

    static const std::vector<std::string> empty_collection; // A representative empty collection of sequences.


public:

    // Constructs a collection of input sequences.
    Seq_Input(const std::vector<std::string>& seqs, const std::vector<std::string>& lists, const std::vector<std::string>& dirs);

    // Constructs a collection of input sequences.
    Seq_Input(const std::optional<std::vector<std::string>>& seqs, const std::optional<std::vector<std::string>>& lists, const std::optional<std::vector<std::string>>& dirs);

    // Returns the collection of paths to raw sequences.
    const std::vector<std::string>& seq_paths() const;

    // Returns the collection of paths to lists containing sequence file paths.
    const std::vector<std::string>& list_paths() const;

    // Returns the collection of paths to directories containing sequence files.
    const std::vector<std::string>& dir_paths() const;

    // Returns the collection of all the input sequences.
    const std::vector<std::string> seqs() const;

    // Returns whether the sequence collection is empty or not.
    bool empty() const;
};



#endif
