/*
 * Copyright (c) 2011-2012, Grégoire Dupont, Matthieu Pérotin
 * 2018, Ying Zhou
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
#ifndef PRINCIPALEXT_H
#define PRINCIPALEXT_H

#include <iostream>
#include <utility>
#include <algorithm>
#include <vector>
#include <map>
#include <string>
#include <cstdlib>
#include <boost/tokenizer.hpp>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdint.h>
#include <gmpxx.h>
#include "Exception.h"
#include "quiver.hpp"
#include "nautinv.h"
#include "linalgext.hpp"
#include <set>
#include <cmath>


typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

class IceQuiver
{
	public:
		IceQuiver(Quiver c);
		IceQuiver(const char* file);
		~IceQuiver();
		IceQuiver(const IceQuiver &ca);
		int mutate(int k);
		void print();
		bool infinite(mpz_class p);
		void toFile(const char* filename);
		void printMutations(int s);
		void printMutationsE(int s);
		void genGraph();
        void printMultiplicityMap();

		/* Getters and Setters */
		void setM(int i, int j, mpz_class val);
		/*
		But: Getter pour la matrix d'incidence
		Entrée: 2 entiers i et j
		Sortie: un entier
		Précondition: i et j compris entre 0 et n-1
		PostCondition: return M[i][j] si i et j conformes, sinon jette une exception

		*/
		inline mpz_class getM(int i, int j)
		{
			if(i<n && j < n && i>=0 && j>=0)
				return M[i*n+j];
			else
			{
				throw Exception("DOMAIN_ERROR: getM");
			}
		}
		inline int getN()
		{
			return n;
		}
		inline int getNbVertexsNauty()
        {
            return nbVerticesNauty;
        }
		inline int lastMutation()
		{
			return (mutations.size()>0)?mutations[mutations.size()-1]:-1;
		}
		void generateGreenVertices();
        void generateFiniteGreenVertices();
		int getNextFiniteGreenVertex();
		int getRandomFiniteGreenVertex();
		void forceGreenVertex(int s);
		inline bool getGraphIsUpToDate()
		{
			return graphIsUpToDate;
		}
		inline void setMultiplicity(uint64_t size, mpz_class mul)
		{
			this->multiplicity[size]=mul;
		}

		inline mpz_class getMultiplicity(uint64_t size)
		{
			if(this->multiplicity.find(size) == this->multiplicity.end())
            {
                return 0;
            }
            else
            {
			    return this->multiplicity[size];
            }
		}

		inline mpz_class incMultiplicity(uint64_t size)
		{
            this->multiplicity[size] += 1;
			return this->multiplicity[size];
		}

		inline mpz_class addMultiplicity(uint64_t size,mpz_class value)
		{
            this->multiplicity[size] += value;
			return this->multiplicity[size];
		}
		void addMultiplicity(IceQuiver &p);
        inline std::map<uint64_t, mpz_class> *getMultiplicityMap()
        {
            return &multiplicity;
        }
        
        inline std::map<mpz_class, mpz_class> *getMultiplicitiesMap()
        {
            return &multiplicities;
        }

        inline std::vector<int> getMutations(void)
        {
            return mutations;
        }
        inline std::string getMutationsString(void)
        {
            return mutationString;
        }
        inline uint64_t getMutationsSize(void)
        {
            return mutationsSize;
        }
        inline vecivect getAdmittedCVectors(void)
        {
            return admittedCVectors;
        }
        inline vecivect getAdmittedGVectors(void)
        {
            return admittedGVectors;
        }
        inline vecivect getAdmittedPVectors(void)
        {
            return admittedPVectors;
        }
        inline imat getCartan(void) {
            return cartan;
        }
        inline imat getEuler(void) {
            return euler;
        }
        inline imat getPhi(void) {
            return phi;
        }
        inline imat getPhiInverse(void) {
            return phiin;
        }
        inline bool isCKnown(void) {
            return is_c_known;
        }
        inline int prepdeg(ivect v) {
            return mFinDeg(phi, v);
        }
        inline int preideg(ivect v) {
            return mFinDeg(phiin, v);
        }
        inline bool isprep(ivect v) {
            return (mFinDeg(phi, v) != -1);
        }
        inline bool isprei(ivect v) {
            return (mFinDeg(phiin, v) != -1);
        }
        inline bool isreg(ivect v) {
            return (mFinDeg(phi, v) == -1) && (mFinDeg(phiin, v) == -1);
        }
        inline bool isTauSincere(ivect v) {
            return isSincere(minvec(phi, phiin, v));
        }
        inline bool isEulerOne(ivect v) {
            return (calcEulerForm(v) == 1);
        }
        inline bool sameComponent(ivect v1, ivect v2) {
            return sameMOrbit(phi, getQuasiSimpleQuot(v1), getQuasiSimpleQuot(v2));
        }
        int calcEulerForm(ivect);
		Quiver *getQuiver(void);
		graph *oldGetNautyGraph();
        std::string mutationsToString();
        void shiftMultiplicity();
        void unshiftMultiplicity();
        void semiDestroy();
        int calcQuasiLength(ivect);
        int vecProc(ivect);
        void vecProcs(vecivect);
        ivect getQuasiSimpleQuot(ivect);
        inline ivect getDegZeroQuasiSimple(ivect v) {
            return minvec(phi, phiin, getQuasiSimpleQuot(v));
        }
        inline bool isComponentSincere(ivect v) {
            return isTauSincere(getQuasiSimpleQuot(v));
        }
    
        vecivect admittedCVectors;
        vecivect admittedGVectors;
        vecivect admittedPVectors;
        imat cartan;
        imat euler;
        imat phi;
        imat phiin;
        bool is_c_known;
	private:
		mpz_class *M;
		int n;
		int nbVerticesNauty;
    //We need a different term for "multiplicities"
		std::map<uint64_t,mpz_class> multiplicity;//Documents the mutation length and the multiplicity.
        std::map<mpz_class,mpz_class> multiplicities;
		std::vector<int> mutations;
		std::vector<int> greenVertices;
        std::vector<int> finiteGreenVertices;
        graph nautyG[MAXN*MAXM];
		graph nautyGC[MAXN*MAXM];
		set *gv;
		bool graphIsUpToDate;
        bool semiFreed;
        std::string mutationString;
        uint64_t mutationsSize;
};
#endif
