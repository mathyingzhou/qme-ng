/*
 * Copyright (c) 2007-2012, Joris Calabrese, 
 *                          Grégoire Dupont, 
 *                          Matthieu Pérotin
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
/* Ying Zhou's notes:
* In essence this quiver (in French Carquois) file does not really depend on anything..except for Exception.h which does not
* count...
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
#define SPORADIQUE 6 /*Sporadic?*/
#define UNAMED 7
#define E_ELIPTIQUE 8 /*E_Elliptic?*/
#define ATILDEALT 9 

typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

class Carquois
{
	public:

		/* Constructeurs et Destructeurs/Constructors and Destructors */
		Carquois(int n);/*Construct a quiver given number of vertices?*/
		Carquois();/*Construct a quiver given nothing?*/
		Carquois(int type, int nbSommets, int orientation=0);/*Sommets = sources?*/
		Carquois(const Carquois &ca);/*Quiver copying?*/
		Carquois(const Carquois &ca, int k);/*Copy a part of a quiver?*/
		Carquois(int ** matrice, int n,int indice);/*Given a matrice = matrix, num of vertices and what..*/
		Carquois(const char* file);/*Read a quiver from .qmu?*/
		~Carquois();/*Destructor*/

		/* Algorithmes/Algorithms */
		void mutate(int k);/*Do mutation at vertex k?*/
		bool infinite();/*Is this green sequence infinite?*/
		void genGraph();/*Generate a graph (or quiver?)*/
		bool testInfiniEmpirique(int mutations);/*Test whether a green sequence is infinite empirically?*/
		void toFile(const char* filename);/*Put a quiver into a file?*/
		void semiDestroy();/*Semi-destroy. Not sure what this means.*/
		int getNbVoisinsMax();/*get Nb adjacent Max?*/
		int estConnexe();/*Is connected*/
		bool aUneDouble(int i);/*a double? What's this?*/
		bool troisCycleOriente(int i, int j, int k);/*Oriented 3-cycle?*/
		bool cyclique();/*Is cyclic?*/

		/* Affichage/Viewing?*/
		void affiche();/*Post (to the screen)?*/
		void printMutations();

		/* Getters et Setters/Get and set functions */
		void setM(int i, int j, int val);/*Set M[i][j] as val.*/
		/*
		But: Getter pour la matrice d'incidence
		Entrée: 2 entiers i et j
		Sortie: un entier
		Précondition: i et j compris entre 0 et n-1
		PostCondition: return M[i][j] si i et j conformes, sinon jette une exception

		*/
		inline int getM(int i, int j)
		{
			if(i<n && j < n && i>=0 && j>=0)
				return M[i][j];
			else
			{
				throw Exception("ERREUR_DOMAINE: getM");
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
		inline int graphEstAJour()
		{
			return graphAJour;
		}
		std::string getMutations();

		inline int getNextI()
		{
			return nextI;
		}
		inline int getNextJ()
		{
			return nextJ;
		}
		inline void setNextI(int i)
		{
			nextI = i;
		}
		inline void setNextJ(int j)
		{
			nextJ = j;
		}
		
		/* Getter avancés pour raisonnements locaux */
		std::vector<int> getVoisins(int sommet);
		std::vector<int> getVoisinsDoubles(int sommet);
		std::vector<int> getVoisinsSimples(int sommet);
		std::vector<int> getVoisinsSimplesPredecesseurs(int sommet);
		std::vector<int> getVoisinsSimplesSuccesseurs(int sommet);
		std::vector<int> getSommetsArreteDoubleEntrante();
		int getSommetOrigineArreteDouble(int i);
		std::vector<int> getSommetsPasDArreteDouble();
		int getNbVoisinsSimplesPredecesseurs(int sommet);
		int getNbVoisinsSimplesSuccesseurs(int sommet);

	private:
		int **M;
		int n;
		int valeurAbs(int k);
		int graphAJour;
		std::vector<int> mutations;
		int score;
		void genScore();
		int semifree;
		int nbVoisinsMax;
		int connexe;
		int nextI;
		int nextJ;
	    graph nautyG[MAXN*MAXM];
		graph nautyGC[MAXN*MAXM];
		set *gv;
       	bool exploreCycle(int *,int);

};
#endif
