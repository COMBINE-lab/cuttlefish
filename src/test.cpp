
#include "Kmer.hpp"
#include "kseq/kseq.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <zlib.h>
#include <cstring>


// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(int, read);


void test_kseq(const char* fileName)
{
    // STEP 2: open the file handler
    FILE* fp = fopen(fileName, "r");

    // STEP 3: initialize seq
    kseq_t* seq = kseq_init(fileno(fp));


    // STEP 4: read sequence
    int c = 0;
    uint64_t maxLen = 0;
    uint64_t len = 0;
    while(true)
    {
        if(kseq_read(seq) < 0)
            break;

        c++;

        std::cout << "Name: " << seq->name.s << "\n";
        if(seq->comment.l)
            std::cout << "Comment: " << seq->comment.s << "\n";

        std::cout << "Seq: " << seq->seq.s << "\n";
        std::cout << "seq.m: " << seq->seq.m << "\n";
        std::cout << "seq.l: " << seq->seq.l << "\n";

        len += seq->seq.l;
        maxLen = seq->seq.l;
        
        
        if(seq->qual.l)
            std::cout << "Quality: " << seq->qual.s << "\n";
    }


    std::cout << "Line count: " << c << "\n";
    std::cout << "Max line length: " << maxLen << "\n";
    std::cout << "Total reference length: " << len << "\n";

    kseq_destroy(seq);
    fclose(fp);
}


int main(int argc, char** argv)
{
    const char* fileName = argv[1];

    test_kseq(fileName);

    return 0;
}
