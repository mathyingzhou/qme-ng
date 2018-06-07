/*
 * Copyright (c) 2012, Grégoire Dupont, Matthieu Pérotin
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "greenSizeHash.hpp"
GreenSizeHash::GreenSizeHash() {}
GreenSizeHash:: ~GreenSizeHash() {}

//Add (size,1) to the entries with key std::vector<int> &mutations and its initial subsequences and change the multiplicities map accordingly.
void GreenSizeHash::increment(std::vector<int> &mutations, uint64_t size)
{
    std::map<uint64_t,mpz_class> temp;
    temp[size] = 1;
    addSizes(mutations,temp);
}

//Add std::map<uint64_t,mpz_class>& sizes to the entries with key std::vector<int> &mutations and its initial subsequences and change the multiplicities map accordingly.
void GreenSizeHash::addSizes(std::vector<int> &mutations, std::map<uint64_t,mpz_class>& sizes)
{
    uint64_t i;
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    std::map<uint64_t,mpz_class> temp;
    std::map<uint64_t,mpz_class>::iterator map_it;
    std::string mutations_str;
    strhash::iterator strhash_it;
    mul_hash::iterator mul_hash_it;
    std::pair<mul_hash::iterator, bool> ins_it;
#ifdef DEBUG
    std::cout << "Before adding sizes" << std::endl;
    std::cout << "print Green Size:" << std::endl;
    printGreenSize();
#endif
    // 1. For each mutation subchain
    for(i=0;i<mutations.size();i++) {
        // Empty the ss object;
        ss.clear();ss.str("");
        ss << mutations_str << mutations[i]+1 << " ";
        mutations_str = ss.str();//Update the mutation string by adding one more mutation (with 1 added)
        // 2. Lookup the mutations_str 
        strhash_it = green_size.find(mutations_str);
        if(strhash_it != green_size.end())
        {
            // The subchain already has a multiplicities list
            temp = *(strhash_it->second);//a map
            mul_hash_it = multiplicities.find(temp);
            mul_hash_it->second -= 1;
            if(mul_hash_it->second == 0)
            {
                multiplicities.erase(temp);
            }
        }
        //Find the corresponding map of the mutation sequence in the existing green_size. If it is found then delete it once from the multiplicities map...aka replacement
        
        // If the subchain does not already have a multiplicities subchain then we do nothing
        
        // Now, add sizes to the old corresponding map
        for(map_it=sizes.begin();map_it!=sizes.end();map_it++)
        {
            temp[map_it->first]+=map_it->second;
        }
        //temp is basically defined as
        //3. Lookup the newly created map
        mul_hash_it = multiplicities.find(temp);
        if(mul_hash_it != multiplicities.end())
        {
            // The created map already exists !
            // a. increment refs
            mul_hash_it->second += 1;
            // b. update
            if(strhash_it != green_size.end())
            {
                strhash_it->second = &(mul_hash_it->first);
            }
            else
            {
                green_size[mutations_str] = &(mul_hash_it->first);
            }
        }
        else
        {
            // The created map does not yet exist !
            // a. insert
            ins_it = multiplicities.insert(std::make_pair(temp,1));
            // b. update
            green_size[mutations_str] = &(ins_it.first->first);
        } 
        temp.clear();
    }
#ifdef DEBUG
    std::cout << "After adding sizes" << std::endl;
    std::cout << "print Green Size:" << std::endl;
    printGreenSize();
#endif
}

//Given a modified mutation string produce the cardinality/size of the domain of the corresponding map.
uint64_t GreenSizeHash::GreenSizesGetSize(std::string s)
{
    strhash::iterator strhash_it;
    strhash_it = green_size.find(s);
    if(strhash_it == green_size.end())
    {
        return 0;
    }
    else
    {
        return strhash_it->second->size();
    }

}

//Given a modified mutation string and a uint64_t, produce the int associated to it.
//This is in essence an evaluation of the form green_size(s)(v).
mpz_class GreenSizeHash::GreenSize(std::string s,uint64_t v)
{
    strhash::iterator strhash_it;
    strhash_it = green_size.find(s);
    if(strhash_it == green_size.end())
    {
        return 0;
    }
    else
    {
        return strhash_it->second->at(v);
    }

}

//Given a map<uint64_t,mpz_class> print the mapped integers. This is also basically evaluation.
mpz_class GreenSizeHash::MultiplicitiesRefs(std::map<uint64_t,mpz_class> &rmap)
{
    mul_hash::iterator it;
    it = multiplicities.find(rmap);
    if(it == multiplicities.end())
    {
        return -1;
    }
    else
    {
        return it->second;
    }
}

//print green_size
void GreenSizeHash::printGreenSize()
{
    strhash::iterator strhash_it;
    std::map<uint64_t,mpz_class>::iterator map_it;
    std::map<uint64_t,mpz_class> tmp;
    
    for(strhash_it = green_size.begin(); strhash_it != green_size.end() ; strhash_it++)
    {
        std::cout << strhash_it->first << std::endl;
        std::cout << "\t";
        tmp = *(strhash_it->second);
        for(map_it=tmp.begin() ; map_it != tmp.end() ; map_it++)
        {
            std::cout << "(" << map_it->first << ";" << (map_it->second).get_str() << ")" ;
        }
        std::cout << std::endl;
    }

}

//print multiplicities
void GreenSizeHash::printMultiplicities()
{
    mul_hash::iterator mul_hash_it;
    std::map<uint64_t,mpz_class>::iterator map_it;
    std::map<uint64_t,mpz_class> tmp;
    
    for(mul_hash_it = multiplicities.begin(); mul_hash_it != multiplicities.end() ; mul_hash_it++)
    {
        tmp = mul_hash_it->first;
        for(map_it=tmp.begin() ; map_it != tmp.end() ; map_it++)
        {
            std::cout << "(" << map_it->first << ";" << (map_it->second).get_str() << ")" ;
        }
        std::cout <<std::endl;
        std::cout << "\t" << mul_hash_it->second << std::endl;
    }
}
