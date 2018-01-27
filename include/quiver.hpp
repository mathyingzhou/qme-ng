/*
 * Copyright (c) 2007-2012, Joris Calabrese, 
 *                          Grégoire Dupont, 
 *                          Matthieu Pérotin
 *                2018, Ying Zhou
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
#ifndef CARQUOIS_H
#define CARQUOIS_H

#include <iostream>
#include <utility>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include <boost/tokenizer.hpp>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <limits.h>
#include "Exception.h"
#define MAXN 100
#include "nauty.h"

#define A 0
#define D 1
#define E 2
#define ATILDE 3
#define DTILDE 4
#define ETILDE 5
#define SPORADIC 6
#define UNAMED 7
#define E_ELIPTIC 8
#define X 9
#define B 10
#define C 11
#define F 12
#define G 13
#define ATILDE14 14
#define BTILDE 15
#define CTILDE 16
#define BCTILDE 17
#define BDTILDE 18
#define CDTILDE 19
#define FTILDE41 20
#define FTILDE42 21
#define GTILDE21 22
#define GTILDE22 23

typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

class Quiver
{
	public:

		/* Constructeurs et Destructeurs */
		Quiver(int n);
		Quiver();
        Quiver(int type, int nbVertices, std::string orientation);
		Quiver(const Quiver &ca);
		Quiver(const Quiver &ca, int k);
		Quiver(int ** matrix, int n,int indice);
		Quiver(const char* file);
		~Quiver();

		/* Algorithmes */
		void mutate(int k);
		bool infinite();
		void genGraph();
		bool testInfiniEmpirique(int mutations);
		void toFile(const char* filename);
		void semiDestroy();
		int getNbNeighboursMax();
		int isConnected();
		bool aUneDouble(int vertex);
		bool isOrientedThreeCycle(int i, int j, int k);
		bool isCyclic();

		/* Printing */
		void print();
		void printMutations();

		/* Getters and Setters */
		void setM(int i, int j, int val);
		/*
		But: Getter pour la matrix d'incidence
		Entrée: 2 entiers i et j
		Sortie: un entier
		Précondition: i et j compris entre 0 et n-1
		PostCondition: return M[i][j] si i et j conformes, sinon jette une exception

		*/
        //get b_{i+1,j+1}
		inline int getM(int i, int j)
		{
			if(i<n && j < n && i>=0 && j>=0)
				return M[i][j];
			else
			{
				throw Exception("DOMAIN_ERROR: getM");
			}
		}
		inline int getN()
		{
			return n;
		}
		graph *getNautyGraph();
		inline int getScore()
		{
			return score;
		}
		inline int lastMutation()
		{
			return (mutations.size()>0)?mutations[mutations.size()-1]:-1;
		}
		inline int graphIsAJour()
		{
			return graphIsUpToDate;
		}
		std::string getMutations();

		//inline int getNextI()
		//{
			//return nextI;
		//}
		//inline int getNextJ()
		//{
			//return nextJ;
		//}
		//inline void setNextI(int i)
		//{
			//nextI = i;
		//}
		//inline void setNextJ(int j)
		//{
			//nextJ = j;
		//}
		
		/* Getter avancés pour raisonnements locaux */
		std::vector<int> getNeighbours(int vertex);
		std::vector<int> getNeighboursDoubles(int vertex);
		std::vector<int> getNeighboursSimples(int vertex);
		std::vector<int> getNeighboursSimplesPredecesseurs(int vertex);
		std::vector<int> getNeighboursSimplesSuccesseurs(int vertex);
		std::vector<int> getVertexsDoubleEdgeEntrante();
		int getVertexOrigineDoubleEdge(int i);
		std::vector<int> getVertexsPasDDoubleEdge();
		int getNbNeighboursSimplesPredecesseurs(int vertex);
		int getNbNeighboursSimplesSuccesseurs(int vertex);

	private:
		int **M;//The exchange matrix
		int n;//The number of vertices
		int absVal(int k);//Take the absolute value
		int graphIsUpToDate;
		std::vector<int> mutations;
		int score;
		void genScore();
		int semifree;
		int nbNeighboursMax;
		int connected;
		//int nextI;
		//int nextJ;
	    graph nautyG[MAXN*MAXM];
		graph nautyGC[MAXN*MAXM];
		set *gv;
       	bool exploreCycle(int *,int);

};
#endif
