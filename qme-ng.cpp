/*
 * Copyright (c) 2005-2012, Grégoire Dupont, Matthieu Pérotin
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

#include <boost/program_options.hpp>
#include "greenexplorator.hpp"
#include "greenfinder.hpp"
#include "quiver.hpp"
#include "mutexploratorSeq.hpp"
#include "mutexplorator.hpp"
#include "time.h"
#include <stdint.h>

int main(int argc, char **argv)
{
    //int **test;
    mpz_class p;
    int max_depth;
    int min_depth;
    int mutationClassSize;
    //bool iso=true;
    int random_tries;
    std::vector<int> size;
    std::string type;
    std::string orientation;
    Quiver *quiver;
    GreenExplorator ge;
    GreenFinder *gf;
    MutExploratorSeq *explorator;
    IceQuiver *pt;
    std::map<std::string,int> type_quiver;
    std::map<std::string,int>::const_iterator it_map;
    /* Initializing type_quiver map */
    type_quiver["A"] = A;
    type_quiver["D"] = D;
    type_quiver["E"] = E;
    type_quiver["ATILDE"] = ATILDE;
    type_quiver["X"] = X;
    type_quiver["DTILDE"] = DTILDE;
    type_quiver["ETILDE"] = ETILDE;
    type_quiver["SPORADIC"] = SPORADIC;
    type_quiver["UNAMED"] = UNAMED;
    type_quiver["E_ELIPTIC"] = E_ELIPTIC;

    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "Print this help")
        ("file,f", boost::program_options::value<std::string>(), "Exchange matrix file (n * n)")
        /*("eefile,e", boost::program_options::value<std::string>(), "Extended exchange matrix file (2n * n)")*/
        ("iqfile,i", boost::program_options::value<std::string>(), "Exchange matrix of ice quiver file (2n * 2n)")
        ("type,t", boost::program_options::value<std::string>(), "Quiver type (A, D, E, ATILDE, DTILDE, ETILDE, SPORADIC, UNAMED, E_ELIPTIC, X)")
        ("size,s", boost::program_options::value< std::vector<int> >(),"Quiver size (must be used with type)")
        ("orientation,o", boost::program_options::value<std::string>(), "Orientation (a sequence of l and r to denote orientation, if orientation is not given the default orientation will be used, can not be used for E_ELIPTIC or X)")
        ("green,g", "Green exploration")
        ("one", boost::program_options::value<int>(&random_tries)->default_value(0), "Find one green suite, give number of tries")
        ("p", boost::program_options::value<mpz_class>(&p)->default_value(0),"P param")
        ("max_depth", boost::program_options::value<int>(&max_depth)->default_value(INT_MAX),"Max exploration depth")
        ("min_depth", boost::program_options::value<int>(&min_depth)->default_value(0),"Min exploration depth")
        ("no-iso,n", "Isomorph discrimination")
        ("dump-class,c", boost::program_options::value<std::string>(), "Dump Mutation Class")
        ("dump-trunk,k", "Dump truncated quivers")
    ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);
    ge.setP(p);
    ge.setMaxDepth(max_depth);
    if(vm.count("no-iso"))
    {
        ge.setIsomorphTest(false);
    }
    if(vm.count("dump-trunk"))
    {
        ge.setDumpTruncated(true);
    }
    try
    {
        //Get the ice quiver
        if(vm.count("file") && vm.count("iqfile")) {throw Exception("Only one type of file can be given !");}
        if(vm.count("file"))
        {
                    quiver = new Quiver(vm["file"].as<std::string>().c_str());
        }
        //Generate quiver from type, size and orientation
        else if(vm.count("type"))
        {
            if(vm.count("size"))
            {
                size = vm["size"].as< std::vector<int> >();
                type = vm["type"].as<std::string>();
                if(vm.count("orientation")) {
                    orientation = vm["orientation"].as<std::string>();
                }
                else {
                    orientation = "";//If orientation is not given then we use the default orientation.
                }
                it_map = type_quiver.find(type);
                if(it_map == type_quiver.end())
                {
                    // Unknown type
                    std::cerr << desc << "\n";
                    return 1;
                }
                else
                {
                    //Type & size
                    if(size.size() == 1) {
                        if ((it_map->second == E_ELIPTIC || it_map->second == X) && vm.count("orientation")) {
                            throw Exception("E_ELIPTIC and X do not allow non-default orientations!");
                        }
                        quiver = new Quiver(it_map->second,size[0], orientation);
                    }
                    else {//absurd size length
                        std::cerr << desc << "\n";
                        return 1;
                    }
                }    
            }
            else
            {
                //No size
                std::cerr << desc << "\n";
                return 1;
            }
        }
        else
        {
            if(!vm.count("iqfile")) {
                std::cerr << desc << "\n";
                return 1;
            }
        }
    } catch (Exception e)
    {
        std::cout << e.m_Msg << "\n";
        return 1;
    }
    //Find maximal green sequences.
    if(vm.count("green"))
    {
        if(vm.count("iqfile")) { pt = new IceQuiver(vm["iqfile"].as<std::string>().c_str());}
        else  {
            pt = new IceQuiver(*quiver);
        }
        //Print the 2n * 2n exchange matrix of the ice quiver
        pt->print();
        pt->generateGreenVertices();
        try {
            if(random_tries == 0) {
                //Get the number of MGS
                ge.greenExploration(*pt);
            }
            else
            //--one, get a random MGS
            {
                gf = new GreenFinder(*pt, p, min_depth, max_depth);
                gf->find(random_tries);
                delete gf;
            }
        } catch (Exception e)
        {
            std::cout << e.m_Msg << "\n";
            delete quiver;
            return 1;
        }
    }
    else
    //Get mutation class size
    {
        explorator = new MutExploratorSeq();
        try {
        mutationClassSize = explorator->isomorphismExplorator(*quiver,5000);
        std::cout << "Set Size:" << mutationClassSize << "\n";
        explorator->getNbNeighboursMax();
        if(vm.count("dump-class")) {
            explorator->dumpFiles((vm["dump-class"].as<std::string>()).c_str());
        }
        if(explorator->isAcyclic()) {
            std::cout << "The mutation class is acyclic !\n";
        }
        else
        {
            std::cout << "The mutation class is not acyclic !\n";
        }
        } catch (Exception e)
        {

            std::cout << e.m_Msg << "\n";
            delete explorator;
            return 1;
        }
        delete explorator;
    }
    //Free the memory
    if(!vm.count("iqfile")) {
        delete quiver;
    }
    delete pt;
    return 0;
}
