
#include "CdBG.hpp"
#include "kseq/kseq.h"

#include <chrono>
#include "zlib.h"


// Declare the type of file handler and the read() function.
// Required for FASTA/FASTQ file reading using the kseq library.
KSEQ_INIT(int, read);


void CdBG::output_maximal_unitigs(const std::string& output_file)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    // Open the file handler for the FASTA / FASTQ file containing the reference.
    FILE* input = fopen(ref_file.c_str(), "r");
    if(input == NULL)
    {
        std::cerr << "Error opening input file " << ref_file << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    // Initialize the parser.
    kseq_t* parser = kseq_init(fileno(input));


    // Open the output file.
    std::ofstream output(output_file.c_str(), std::ofstream::out);
    if(!output)
    {
        std::cerr << "Error opening output file " << output_file << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    // Parse sequences one-by-one, and output each unique maximal unitig encountered through them.
    uint32_t seqCount = 0;
    while(kseq_read(parser) >= 0)
    {
        const char* seq = parser->seq.s;
        const size_t seq_len = parser->seq.l;

        std::cout << "Processing sequence " << ++seqCount << ", with length " << seq_len << ".\n";

        // Nothing to process for sequences with length shorter than `k`.
        if(seq_len < k)
            break;


        // Look for contiguous subsequences of k-mers without the placeholder nucleotide 'N'.
        size_t kmer_idx = 0;
        while(kmer_idx <= seq_len - k)
        {
            kmer_idx = search_valid_kmer(seq, kmer_idx, seq_len - k);

            // No valid k-mer remains in the sequence.
            if(kmer_idx > seq_len - k)
                break;

            // Process a maximal valid contiguous subsequence, and advance to the index following it.
            kmer_idx = output_maximal_unitigs(seq, seq_len, kmer_idx, output);
        }
    }

    
    // Close the output file.
    if(output.fail())
        std::cerr << "Errors had been encountered for the output stream.\n";

    output.close();


    // Close the parser and the input file.
    kseq_destroy(parser);
    fclose(input);


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done outputting the maximal unitigs. Time taken = " << elapsed_seconds << " seconds.\n";
}


size_t CdBG::output_maximal_unitigs(const char* seq, const size_t seq_len, const size_t start_idx, std::ofstream& output)
{
    size_t kmer_idx = start_idx;

    // assert(kmer_idx <= seq_len - k);

    Annotated_Kmer curr_annot_kmer(cuttlefish::kmer_t(seq, kmer_idx), kmer_idx, Vertices);

    // The k-mer is an isolated one, so is a maximal unitig by itself.
    if(kmer_idx + k == seq_len || seq[kmer_idx + k] == 'N')
        output_unitig(seq, curr_annot_kmer, curr_annot_kmer, output);
    else    // At least two adjacent k-mers are present from the index `kmer_idx`.
    {
        Annotated_Kmer prev_annot_kmer;
        Annotated_Kmer next_annot_kmer = curr_annot_kmer;
        next_annot_kmer.roll_to_next_kmer(seq[kmer_idx + k], Vertices);

        Annotated_Kmer unipath_start_kmer;


        // Process the first k-mer of this subsequence.

        // A maximal unitig starts at the beginning of a maximal valid subsequence.
        unipath_start_kmer = curr_annot_kmer;
        if(is_unipath_end(curr_annot_kmer.vertex_class, curr_annot_kmer.dir, next_annot_kmer.vertex_class, next_annot_kmer.dir))
            output_unitig(seq, unipath_start_kmer, curr_annot_kmer, output);


        // Process the internal k-mers of this subsequence.
        for(kmer_idx++; kmer_idx < seq_len - k && seq[kmer_idx + k] != 'N'; ++kmer_idx)
        {
            prev_annot_kmer = curr_annot_kmer;
            curr_annot_kmer = next_annot_kmer;
            next_annot_kmer.roll_to_next_kmer(seq[kmer_idx + k], Vertices);


            if(is_unipath_start(curr_annot_kmer.vertex_class, curr_annot_kmer.dir, prev_annot_kmer.vertex_class, prev_annot_kmer.dir))
                unipath_start_kmer = curr_annot_kmer;

            if(is_unipath_end(curr_annot_kmer.vertex_class, curr_annot_kmer.dir, next_annot_kmer.vertex_class, next_annot_kmer.dir))
                output_unitig(seq, unipath_start_kmer, curr_annot_kmer, output);
        }

        prev_annot_kmer = curr_annot_kmer;
        curr_annot_kmer = next_annot_kmer;


        // Process the last k-mer of this subsequence.

        if(is_unipath_start(curr_annot_kmer.vertex_class, curr_annot_kmer.dir, prev_annot_kmer.vertex_class, prev_annot_kmer.dir))
            unipath_start_kmer = curr_annot_kmer;

        // A maximal unitig ends at the ending of a maximal valid subsequence.
        output_unitig(seq, unipath_start_kmer, curr_annot_kmer, output);
    }


    // Return the non-inclusive ending index of the processed contiguous subsequence.
    return kmer_idx + k;
}


bool CdBG::is_unipath_start(const cuttlefish::Vertex_Class vertex_class, const cuttlefish::kmer_dir_t dir, const cuttlefish::Vertex_Class prev_kmer_class, const cuttlefish::kmer_dir_t prev_kmer_dir) const
{
    if(vertex_class == cuttlefish::Vertex_Class::multi_in_multi_out)
        return true;

    if(dir == cuttlefish::FWD)
    {
        if(vertex_class == cuttlefish::Vertex_Class::multi_in_single_out)
            return true;
    }
    else    // dir == cuttlefish::BWD
        if(vertex_class == cuttlefish::Vertex_Class::single_in_multi_out)
            return true;


    // assert(kmer_idx > 0);


    if(prev_kmer_class == cuttlefish::Vertex_Class::multi_in_multi_out)
        return true;

    if(prev_kmer_dir == cuttlefish::FWD)
    {
        if(prev_kmer_class == cuttlefish::Vertex_Class::single_in_multi_out)
            return true;
    }
    else    // prev_kmer_dir == cuttlefish::BWD
        if(prev_kmer_class == cuttlefish::Vertex_Class::multi_in_single_out)
            return true;

    
    return false;
}


bool CdBG::is_unipath_end(const cuttlefish::Vertex_Class vertex_class, const cuttlefish::kmer_dir_t dir, const cuttlefish::Vertex_Class next_kmer_class, const cuttlefish::kmer_dir_t next_kmer_dir) const
{
    if(vertex_class == cuttlefish::Vertex_Class::multi_in_multi_out)
        return true;

    if(dir == cuttlefish::FWD)
    {
        if(vertex_class == cuttlefish::Vertex_Class::single_in_multi_out)
            return true;
    }
    else    // dir == cuttlefish::BWD
        if(vertex_class == cuttlefish::Vertex_Class::multi_in_single_out)
            return true;


    // assert(kmer_idx < ref.length() - k);


    if(next_kmer_class == cuttlefish::Vertex_Class::multi_in_multi_out)
        return true;

    if(next_kmer_dir == cuttlefish::FWD)
    {
        if(next_kmer_class == cuttlefish::Vertex_Class::multi_in_single_out)
            return true;
    }
    else    // next_kmer_dir == cuttlefish::BWD
        if(next_kmer_class == cuttlefish::Vertex_Class::single_in_multi_out)
            return true;


    return false;
}


void CdBG::output_unitig(const char* seq, const Annotated_Kmer& start_kmer, const Annotated_Kmer& end_kmer, std::ofstream& output)
{
    // This is to avoid race conditions that may arise while multi-threading.
    // If two threads try to output the same unitig at the same time but
    // encounter it in the opposite orientations, then data races may arise.
    // For a particular unitig, always query the same well-defined canonical flanking
    // k-mer, irrespective of which direction the unitig may be traversed at.
    const cuttlefish::kmer_t min_flanking_kmer = std::min(start_kmer.canonical, end_kmer.canonical);
    Kmer_Hash_Entry_API hash_table_entry = Vertices[min_flanking_kmer];
    Vertex_Encoding& vertex_encoding = hash_table_entry.get_vertex_encoding();

    if(vertex_encoding.is_outputted())
        return;
    

    vertex_encoding = vertex_encoding.outputted();

    // If the hash table update is successful, only then this thread may output this unitig.
    if(Vertices.update(hash_table_entry))
        write_path(seq, start_kmer.idx, end_kmer.idx, start_kmer.kmer < end_kmer.rev_compl, output);
}


void CdBG::write_path(const char* seq, const uint32_t start_kmer_idx, const uint32_t end_kmer_idx, const bool in_forward, std::ofstream& output) const
{
    if(in_forward)
        for(uint32_t idx = start_kmer_idx; idx <= end_kmer_idx + k - 1; ++idx)
            output << seq[idx];
    else
    {
        // To avoid underflow of unsigned integers, the flanking indices are incremented by 1.
        uint32_t idx = end_kmer_idx + k;
        while(idx > start_kmer_idx)
        {
            idx--;
            output << complement(seq[idx]);
        }
    }

    output << "\n";
}
