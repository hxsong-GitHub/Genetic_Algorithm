//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XX  FindCombinations.cpp
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "FindCombinations.h"

#include <vector>
#include <stdexcept>
#include <algorithm>

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XX  FindPermuations()
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

void FindPermuations(std::vector<long> vec, std::vector<std::vector<long> > & perms) 
{
    perms.resize(0);

    do                          // find all permutations
    {
        perms.push_back(vec);
        
    } while ( std::next_permutation(vec.begin(), vec.end()) );
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XX  FindCombinations()
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

void FindCombinations(long m, long n, std::vector<std::vector<long> > & combs)  // combinations of length m from {1,..,n}
{
    combs.resize(0);
    
    if(m > n) throw std::runtime_error("FindCombinations: m is greater than n");
    
    std::vector<long> comb(m + 1, 0);             
    
    std::vector<bool> v(n);                       // boolean mask
    std::fill(v.begin() + n - m, v.end(), true);  // initially m end bits set  

    do 
    {
        long k = 0;
        for (long i = 0; i != n; ++i) 
        {
            if (v[i]) 
            {
                ++k;
                comb[k] = i + 1;                  // always has 0 in start position
            }                                     // (because seat 0 is always occupied)
        }
        
        combs.push_back(comb);
        
    } while (std::next_permutation(v.begin(), v.end())); // iterates through all combinations
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XX  end of file
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

















