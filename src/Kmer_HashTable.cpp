
#include "Kmer_HashTable.hpp"


void Kmer_HashTable::clear()
{
    hash.clear();
}


void Kmer_HashTable::print_hash_table() const
{
    for(auto key_val: hash)
        std::cout << key_val.first << " : " << key_val.second.decode() << "\n";
}
