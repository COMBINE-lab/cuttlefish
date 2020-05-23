
#include "Kmer.hpp"

#include <iostream>
#include <string>

int main()
{
    // freopen("input.txt", "r", stdin);

    std::string label;
    while(std::cin >> label)
    {
        Kmer::set_k(label.length());
        Kmer kmer(label);

        std::cout << "Kmer     : " << kmer.string_label() << " (" << kmer.int_label() << ")\n";
        std::cout << "RC       : " << kmer.reverse_complement() << "\n";
        std::cout << "Canonical: " << kmer.canonical() << "\n";
    }


    return 0;
}