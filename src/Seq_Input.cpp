
#include "Seq_Input.hpp"


Seq_Input::Seq_Input(   const std::vector<std::string>& seqs,
                        const std::vector<std::string>& lists,
                        const std::vector<std::string>& dirs):
    seq_paths_(seqs),
    list_paths_(lists),
    dir_paths_(dirs)
{}


const std::vector<std::string>& Seq_Input::seq_paths() const
{
    return seq_paths_;
}


const std::vector<std::string>& Seq_Input::list_paths() const
{
    return list_paths_;
}


const std::vector<std::string>& Seq_Input::dir_paths() const
{
    return dir_paths_;
}


bool Seq_Input::empty() const
{
    return seq_paths_.empty() && list_paths_.empty() && dir_paths_.empty();
}
