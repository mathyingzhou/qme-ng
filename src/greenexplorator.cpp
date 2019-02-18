/*
 * Copyright (c) 2007-2012, Grégoire Dupont, Matthieu Pérotin
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
#include "greenexplorator.hpp"

GreenExplorator::GreenExplorator()
{
    numGreen=0;
    minLength=std::numeric_limits<uint64_t>::max();
    maxLength=0;
    truncated = 0;
    isomorphTest = true;
    dumpTruncated = false;
    infCut = 0;
    depthCut = 0;
}
GreenExplorator::~GreenExplorator()
{
}
void GreenExplorator::printTree()
{
    std::list<IceQuiver>::iterator i;
    int j=0;
    for(i=c.begin();i!=c.end();i++)
    {
        std::cout << j++ << ":";
        i->printMutations(0);
    }
}

void GreenExplorator::clearC()
{
    c.clear();
}
        
int GreenExplorator::generateMutations(IceQuiver &pe)
{
    //int i;
    int ret;
    int vertex;
    //int create=0;
    unsigned int size;
    std::string mutations_str="";
    std::string filename = "";
    std::vector<int> mutations_v;
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    IceQuiver p = pe;
    strhash::iterator iit;
    std::map<uint64_t,mpz_class>::iterator mult_it;
    std::map<uint64_t,mpz_class> *mult;
    vecivect::iterator vi;
    #ifdef DEBUG
    std::cout << "genMut: Working on "; pe.printMutationsE(0);
    #endif
    vertex = pe.getNextFiniteGreenVertex();
    if(vertex == -1) {
        return 1;//No more finite green vertex
    }
    if(vertex==p.lastMutation())
    {
        // No need to mutate twice on the same vertex
        return 1;
    }
    if(p.getMutationsSize() >= this->max_depth) {
        return 4;//Cut the sequence
    }
    ret = p.mutate(vertex);
    if (ret == 0) {
        // if mutate returned 0, then non-BHIT violating infinity was detected
        return 0;
    }
    if (ret == 1) {
        //A maximal green sequence/tail was detected or infinity was detected!
        size = p.getMutationsSize();
        mult = p.getMultiplicityMap();
        printVectors(p.getAdmittedCVectors());
        if (!p.getAdmittedCVectors().empty()) {
            for (vi = p.admittedCVectors.begin(); vi != p.admittedCVectors.end(); vi++) {
                if (!belongsTo(*vi, admissibleCVectors))
                    admissibleCVectors.push_back(*vi);
            }
        }
        printVectors(p.getAdmittedGVectors());
        if (!p.getAdmittedGVectors().empty()) {
            for (vi = p.admittedGVectors.begin(); vi != p.admittedGVectors.end(); vi++) {
                if (!belongsTo(*vi, admissibleGVectors))
                    admissibleGVectors.push_back(*vi);
            }
        }
        printVectors(p.getAdmittedPVectors());
        if (!p.getAdmittedPVectors().empty()) {
            for (vi = p.admittedPVectors.begin(); vi != p.admittedPVectors.end(); vi++) {
                if (!belongsTo(*vi, admissiblePVectors))
                    admissiblePVectors.push_back(*vi);
            }
        }
        #ifdef DEBUG
        printVectors(admissibleCVectors);
        printVectors(admissibleGVectors);
        printVectors(admissiblePVectors);
        std::cout << "MGS found!" << std::endl;
        p.printMutationsE(0);
        #endif
        for(mult_it = mult->begin(); mult_it!=mult->end(); mult_it++) { 
            mgsInfo[mult_it->first] += mult_it->second;
#ifdef DEBUG
            std::cout << mult_it->first << "," << mult_it->second << std::endl;
#endif
        }
        
        if (isomorphTest) {
            mutations_v = p.getMutations();
            gsh.increment(mutations_v,size);
        }
        if(size > maxLength) { 
            maxLength = size; 
           // std::cerr << "M:";
            //p.printMutationsE(1);
            //std::cerr << "Q:";
            //p.printMutationsE(0);
        }
        if(size < minLength) { 
            minLength = size; 
            //std::cerr << "m:";
            //p.printMutationsE(1);
            //std::cerr << "Q:";
            //p.printMutationsE(0);
        }

        numGreen+=1;
        if((numGreen % 100000) == 0) {
//#ifdef DEBUG
            std::cerr << "S (" << numGreen << "):";
            p.printMutationsE(0);
//#endif
        }
        if(isomorphTest) {insertInCemetary(p,cemetary);}
        return 0;
    }
    if (isomorphTest) {
        return insertInList(p);
    }
    else {
        c.push_back(p);
        return 2;
    }
    /*#ifdef DEBUG
    std::cout << "Fin du travail avec "; p.printMutations(0);    
    #endif*/
}

int GreenExplorator::insertInList(IceQuiver &pe)
{
    std::list<IceQuiver>::iterator ri; 
    std::list<IceQuiver>::reverse_iterator rxi;
    std::map<uint64_t,mpz_class>::iterator it_mul;
    std::map<uint64_t,mpz_class> *mul_map;
    std::map<uint64_t,mpz_class> tmp;
    std::map<uint64_t,mpz_class> *green_sizes;
    std::map<uint64_t,mpz_class>::iterator it;
    strhash::iterator iit;
    std::string mutations_str="";
    std::vector<int> mutations_v;
    vecivect::iterator vi;
    //int i;
    uint64_t pe_size,n_size;
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    // ri is an  iterator, it browses the list from the beginning
    //1.If pe is in c then add the multiplicity of pe to the already existing copy of it (*ri).
    for(ri=c.begin();ri!=c.end();ri++)
    {
        if(this->myIsomorphismNauty(pe,*ri))
        {
            #ifdef DEBUG
                pe.printMutations(0); 
                std::cout << "Is isomorphic to  ";
                (*ri).printMutations(0); 
                std::cout << "\n";
            #endif
            // This principal extension is still to be considered
            // increase its multiplicity
            (*ri).addMultiplicity(pe); //(!!)
            break;
        }
    }
    if( ri != c.end()) { return 0;}//No addition in this case. Just change the multiplicity.
    for(rxi=cemetary.rbegin();rxi!=cemetary.rend();rxi++)
    {
        if(this->myIsomorphismNauty(pe,*rxi))
        {
            #ifdef DEBUG
                pe.printMutations(0); 
                std::cout << "Is isomorphic to (C) ";
                (*rxi).printMutations(0); 
                std::cout << "\n";
            #endif
            // This quiver is isomorph to a cemetary quiver

            mutations_str = (*rxi).getMutationsString();
            if(gsh.GreenSizesGetSize(mutations_str) != 0) {
                //So each mutation sequence that is an initial subsequence of an MGS is stored in gsh? Huh!
                // The quiver is isomorphic to a quiver leading to an maximal green sequence
                // Update the length list ! (Huh! What's this for!)

                // 2. For all the quiver sizes attainable with the quiver
                green_sizes=gsh.getGreenSizes(mutations_str);
                if (!pe.getAdmittedCVectors().empty()) {
                    for (vi = pe.admittedCVectors.begin(); vi != pe.admittedCVectors.end(); vi++) {
                        if (!belongsTo(*vi, admissibleCVectors))
                            admissibleCVectors.push_back(*vi);
                    }
                }
                std::cout << "admissibleCVecs:" << std::endl;
                printVectors(admissibleCVectors);
                if (!pe.getAdmittedGVectors().empty()) {
                    for (vi = pe.admittedGVectors.begin(); vi != pe.admittedGVectors.end(); vi++) {
                        if (!belongsTo(*vi, admissibleGVectors))
                            admissibleGVectors.push_back(*vi);
                    }
                }
                std::cout << "admissibleGVecs:" << std::endl;
                printVectors(admissibleGVectors);
                if (!pe.getAdmittedPVectors().empty()) {
                    for (vi = pe.admittedPVectors.begin(); vi != pe.admittedPVectors.end(); vi++) {
                        if (!belongsTo(*vi, admissiblePVectors))
                            admissiblePVectors.push_back(*vi);
                    }
                }
                std::cout << "admissiblePVecs:" << std::endl;
                printVectors(admissiblePVectors);
                for(it=green_sizes->begin();it!=green_sizes->end();it++) {
                    // For all the sizes multiplicities
                    mul_map = pe.getMultiplicityMap();
#ifdef DEBUG
                    //pe.printMultiplicityMap();
#endif
                    for(it_mul=mul_map->begin();it_mul != mul_map->end();it_mul++) {
                        pe_size=it_mul->first;
                        n_size=it->first-(*rxi).getMutationsSize()+pe_size;//MGS size
                        mgsInfo[n_size]+=pe.getMultiplicity(it_mul->first)* it->second;
                        //it_mul->first is the length of a mutation sequence from B0 to pe. Hence pe.getMultiplicity(it_mul->first) is literally just its corresponding multiplicity.
                        #ifdef DEBUG
                        std::cout << "Insert MGS of size: " <<  n_size <<
                                     " with multiplicity " << "("<<pe.getMultiplicity(it_mul->first) << "*" <<it->second<<") ="
                                             << pe.getMultiplicity(it_mul->first) *it->second 
                                             << std::endl;
                        #endif
                    }
                    n_size=it->first-(*rxi).getMutationsSize()+pe.getMutationsSize();
                    tmp[n_size] = it->second;
                }
                
                // 3. Update the quiver list
                mutations_v = pe.getMutations();
                gsh.addSizes(mutations_v,tmp);//new mutation sequence and (n_size, it->second) pair
            }
            break;
            
        }
    }
    // if ri went all the way through the end, then the quiver is not
    // ismomorph to any quiver in already in the list, so we add it
    //4.Add it if it is not in the cemetary.
    if(rxi==cemetary.rend())
    {
        c.push_back(pe);
        #ifdef DEBUG
        std::cout << "InsertInList: Adding mutation sequence ";
            pe.printMutations(0);
        #endif
        // Insertion done, return 2;
        return 2;
        
    }
    // No insertion done, return 0
    return 0;
}

//Insert a quiver into the cemetary (i.e. list of explored quivers) if it is not already there.
int GreenExplorator::insertInCemetary(IceQuiver &pe, std::list<IceQuiver> &cem)
{
    std::list<IceQuiver>::iterator ri; 
    std::list<IceQuiver>::reverse_iterator rxi;
#ifdef DEBUG
    fprintf(stderr,"insertInCemetary\n");
#endif

    // ri is an  iterator, it browses the list from the beginning
    for(ri=cem.begin();ri!=cem.end();ri++)
    {
        if(this->myIsomorphismNauty(pe,*ri))
        {
            #ifdef DEBUG
                pe.printMutations(0); 
                std::cout << "is isomorphic to  ";
                (*ri).printMutations(0); 
                std::cout << "\n";
            #endif
            break;
            
        }
    }
    // if ri went all the way through the end, then the quiver is not
    // ismomorphic to any quiver in already in the cemetary, so we add it
    if(ri==cem.end())
    {
        pe.semiDestroy();
        cem.push_back(pe);
        #ifdef DEBUG
        std::cout << "Adding to cemetary: ";
            pe.printMutations(0);
        #endif
        // Insertion done, return 2;
        return 2;
        
    }
    // No insertion done, return 0
    return 0;
}

int GreenExplorator::greenExploration(IceQuiver pe)
{
    int index = 0;
    int ret;
    int max;
    std::map<uint64_t,mpz_class>::iterator it;
    mpz_class total=0;
    std::list<IceQuiver>::iterator pei;
    std::list<IceQuiver>::iterator peitest;
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    std::string filename;
    vecivect::iterator vsi;
    std::vector<mpz_class>::iterator veci;
    uint64_t cutPending = 0;
    // Initial population of the list
    insertInList(pe);
    if (pe.getN() == 2) {//When n=1
        std::cout << 1 << "\t=>\t" << 1 << std::endl;
        std::cout << "Total: " << 1 << std::endl;
        return 1;
    }
    pei = c.begin();
//#ifdef DEBUG
 //   for(peitest = c.begin();peitest != c.end();peitest++) {
   //     fprintf(stderr, "Now printing all mutation sequences in the list c.\n");
     //   peitest->printMutationsE(0);
    //}
//#endif
    while(generateMutations(*pei) != 1){
        //The purpose of this loop is to completely exhaust all initial green sequences starting from a quiver
        //Then it can be removed from the list of to-be-determined quivers.
#ifdef DEBUG
        fprintf(stderr, "Now pei points to ");
        pei->printMutationsE(0);
#endif
        pei=c.begin();
#ifdef DEBUG
        fprintf(stderr, "Now printing all mutation sequences in the list c.\n");
        for(peitest = c.begin();peitest != c.end();peitest++) {
            peitest->printMutationsE(0);
        }
#endif
    };
//#ifdef DEBUG
  //  fprintf(stderr, "Now this weird process is finally over.\nNow pei points to ");
    //pei->printMutationsE(0);
    //fprintf(stderr, "Now printing all mutation sequences in the list c.\n");
    //for(peitest = c.begin();peitest != c.end();peitest++) {
      //  peitest->printMutationsE(0);
    //}
//#endif
    if (isomorphTest) {
        insertInCemetary(*pei,cemetary);
    }
    c.erase(pei);
    pei = c.begin();
    for(index=c.size();index>=1;index--) {
        while(generateMutations(*pei) != 1){
            pei=c.begin();
//#ifdef DEBUG
  //          fprintf(stderr, "Now pei points to ");
    //        pei->printMutationsE(0);
      //      fprintf(stderr, "Now printing all mutation sequences in the list c.\n");
        //    for(peitest = c.begin();peitest != c.end();peitest++) {
         //       peitest->printMutationsE(0);
          //  }
//#endif
        };
        if (isomorphTest) {
            insertInCemetary(*pei,cemetary);
        }
        c.erase(pei);
        pei = c.begin();
//#ifdef DEBUG
  //      fprintf(stderr, "Now pei points to ");
    //    pei->printMutationsE(0);
      //  fprintf(stderr, "Now printing all mutation sequences in the list c.\n");
        //for(peitest = c.begin();peitest != c.end();peitest++) {
          //  peitest->printMutationsE(0);
  //      }
//#endif
    }
    pei=c.end();
    pei--;
    // Main loop
    while(c.size()!=0) {
    #ifdef DEBUG
    std::cout << "C.size: " << c.size() << " Cemetary.size: " << cemetary.size() << "\t\t";
    std::cout << "Working on "; (*pei).printMutations(0);
    (*pei).print();
    #endif
        ret = generateMutations(*pei);
        switch(ret) {
            case 4:
                // Branch cut
            case 3:
                // Infinity detected on the branch
                // Cut the branch !
                if(dumpTruncated) {
                    ss.clear();
                    if(ret == 3) {
                        ss << "DInf_" << infCut << ".quiv";
                    }
                    if(ret == 4) {
                        ss << "DTrunc_" << depthCut << ".quiv";
                    }
                    ss >> filename;
                    (*pei).toFile(filename.c_str());
                }
                if(ret == 3) {
                    infCut++;
                }
                if(ret == 4) {
                    this->truncated = 1;
                    //std::cout << "Cut sequence: " << std::endl;
                    //pei->printMutationsE(0);
                    depthCut ++;
                }
                if(isomorphTest) {
                    if(((*pei).getMultiplicityMap())->size() > 0) {
                        cutPending++; 
                    }
                }
                
            case 1:
                if (isomorphTest) {
                    insertInCemetary(*pei,cemetary);
                }
                c.erase(pei);
            case 2:
                pei=c.end();pei--;
                break;
                // Nothing to be done for case 0
                // a green quiver was detected
                // and the list is not empty... some
                //mutations remains to be explored !
        }
    }
    // Print results
    for ( it=mgsInfo.begin() ; it != mgsInfo.end(); it++ ) {
        std::cout << (*it).first << "\t=>\t" << (*it).second << std::endl;
        total +=(*it).second;
    }
    std::cout << "Total: " << total << std::endl;
    if(this->truncated == 1) {
        std::cout << "Exploration Truncated at depth: " << max_depth << std::endl;
    }
    if(this->infCut > 0) {
        std::cout << "Num branches cut because of BHIT violations: " << infCut << std::endl;
    }
    if(this->depthCut > 0) {
        std::cout << "Num branches cut because of excessive depth: " << depthCut << std::endl;
    }
    if(cutPending > 0) {
        std::cout << "WARNING: "<< cutPending << " branches were cut with pending isomorphs." << std::endl;
        std::cout << "All max green sequences < "<< max_depth << " may not have been found." << std::endl;
    }
    //Print admissible c-vectors
    std::cout << "Final CVec list:" << std::endl;
    printVectors(admissibleCVectors);
    //vecProcs(pe, admissibleCVectors);
    std::cout << "Final GVec list:" << std::endl;
    printVectors(admissibleGVectors);
    std::cout << "Final PVec list:" << std::endl;
    printVectors(admissiblePVectors);
    std::cout << "CProc" << std::endl;
    pe.vecProcs(admissibleCVectors);
    std::cout << "PProc" << std::endl;
    pe.vecProcs(admissiblePVectors);
    if(total)
        return mgsInfo.rbegin()->first;
    else
        return 0;
}

bool GreenExplorator::myIsomorphismNauty(IceQuiver &a, IceQuiver &b)
{
    int i;
    //int n = a.getN();
    int nbNautyVert;
    std::map<mpz_class, mpz_class> *mul_a;
    std::map<mpz_class, mpz_class> *mul_b;
    std::map<mpz_class, mpz_class>::iterator ita,itb,mulend;

    graph *c1;
    graph *c2;


    // These two calls must be placed before hand
    // They are responsible for initializing all the other variables of
    // objects (nbVerticesNauty, multiplicities...)
    c2 = (graph *)a.oldGetNautyGraph();
    c1 = (graph *)b.oldGetNautyGraph();

    // This is very different from the number of vertices
    // This is the number of vertices in the nauty graph
    // which depends on the edge multiplicities
    nbNautyVert = a.getNbVertexsNauty();
    if(nbNautyVert!=b.getNbVertexsNauty()) { 
        return false;
    }

    // If here, the number of verticies is the same
    // We must ensure that the actual multiplicity map is
    // the same
  
    mul_a = a.getMultiplicitiesMap();
    mul_b = b.getMultiplicitiesMap();
    
    if(mul_a->size() != mul_b->size())
    { 
        return false;
    }
    for(ita=mul_a->begin(),itb=mul_b->begin(),mulend=mul_a->end();ita!=mulend;ita++,itb++)
    {
        if(ita->first  != itb->first)  {
            return false;
        }
        if(ita->second != itb->second) {
            return false;
        }
    }
    
    // Now that we are certain the multiplicities are the same,
    // compare nauty map 
    for(i=0;i<nbNautyVert;i++)
    {
        if(c1[i] != c2[i])
        {
            return false;
        }
    }

    // All possibilities are exhausted
    // The two graphs are Isomorphs !
    return true;
}
