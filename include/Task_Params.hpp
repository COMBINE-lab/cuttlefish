
#ifndef TASK_PARAMS_HPP
#define TASK_PARAMS_HPP



#include "globals.hpp"

#include <cstddef>
#include <cstdint>


// Wrapper over the parameters for the classification task.
struct Classification_Task_Params
{
    const char* seq;
    size_t seq_len;
    size_t left_end;
    size_t right_end;

    
    Classification_Task_Params() {}

    Classification_Task_Params(const char* const seq, const size_t seq_len, const size_t left_end, const size_t right_end):
        seq(seq), seq_len(seq_len), left_end(left_end), right_end(right_end)
    {}
};


// Wrapper over the parameters for the output task.
struct Output_Task_Params
{
    uint16_t thread_id;
    const char* seq;
    size_t seq_len;
    size_t left_end;
    size_t right_end;


    Output_Task_Params() {}

    Output_Task_Params(const uint16_t thread_id, const char* const seq, const size_t seq_len, const size_t left_end, const size_t right_end):
        thread_id(thread_id), seq(seq), seq_len(seq_len), left_end(left_end), right_end(right_end)
    {}
};


// Wrapper over the parameters for the DFA-states computation and the maximal unitigs extraction tasks for read-dBGs.
struct Read_dBG_Compaction_Params
{
    uint16_t thread_id;


    Read_dBG_Compaction_Params() {}

    Read_dBG_Compaction_Params(const uint16_t thread_id):
        thread_id(thread_id)
    {}
};



#endif
