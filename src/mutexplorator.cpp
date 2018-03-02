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
#include "mutexplorator.hpp"
MutExplorator::MutExplorator()
{
}
MutExplorator::~MutExplorator()
{
}
void MutExplorator::printTree()
{
	std::vector<Quiver>::iterator i;
	int j=0;
	for(i=c.begin();i!=c.end();i++)
	{
		std::cout << j++ << ":";
		i->printMutations();
	}
}
Quiver MutExplorator::getRepresentant()
{
	std::vector<Quiver>::iterator i,j;
	int meilleurScore = INT_MAX;
	for(i=c.begin();i!=c.end();i++)
	{
		if (i->getScore() < meilleurScore)
		{
			j=i;
			meilleurScore = i->getScore();
		}
		
	}
	return *j;
}
void MutExplorator::clearC()
{
	c.clear();
}
		
void MutExplorator::generateMutations(Quiver quiver)
{
	int i;
	Quiver c = quiver;
	#ifdef DEBUG
	std::cout << "Working on "; c.printMutations();
	#endif
	for(i=0;i<quiver.getN();i++)
	{
		if(i==c.lastMutation())
		{
			// No need to mutate twice on the same vertex
			continue;
		}
		if(i<c.lastMutation() && quiver.getM(i,c.lastMutation()) == 0)
		{
			/* A very difficult optimization here... I should have commented
			 * 3 years ago !
             * when two vertices are not connected then the mutation function is commutative
			 * Hence, if u and v are two vertices, and M[u][v] = 0
			 * µ_u.µ_v (Q) = µ_v. µ_u (Q)
			 *
			 * Let (S, y, x) be the sequence of mutations applied to Quiver so far
			 * We are right now wondering if we should add the Quiver (S,y,x,i)
			 * with (i < x)
			 * 
			 * It must be noted that before adding (S, y, x) to the Quiver list, 
			 * (S, y, i) was considered and analyzed (as i < x). For this quiver, two cases could
			 * have happened:
			 *   1. it was dismissed as isomorph to another quiver
			 *   2. it was added as a new quiver
			 * In any case, it is now present in the quiver list, as well as its mutations, in
			 * particular, (S, y, i, x)
			 * 
			 * i and x being unconnected, (S, y, i, x) and (S, y, x, i) are isomorphs, there
			 * is no need to consider this mutation ! 
			 */
			continue;
		}
		c.mutate(i);
		#ifdef DEBUG
		std::cout << "Analyzing ";
		c.printMutations();
		//c.print();
		#endif
		insertInList(&c);
		c.mutate(i);
	}
	#ifdef DEBUG
	std::cout << "We have finished the work on "; c.printMutations();	
	#endif
}
void MutExplorator::insertInList(Quiver *quiver)
{
	std::vector<Quiver>::iterator i; 
	std::vector<Quiver>::reverse_iterator ri;
	quiver->genGraph();

	if(modeComparaison == 1)
	{
		i=c.begin();
		if(myIsomorphismNauty(quiver,&(*i)))
		{
			std::cout << "The quiver belongs to the mutation class!";
			std::cout << "\t Mutations:";
			quiver->printMutations();
			throw Exception("Done.");
		
		}
	}
	// ri is a reverse iterator, it browse the list from the end
	for(ri=c.rbegin();ri!=c.rend();ri++)
	{
        if(myIsomorphismNauty(quiver,&(*ri)))
        {
            #ifdef DEBUG
                quiver->printMutations(); 
                std::cout << "Is isomorphic to  "; 
                (*ri).printMutations(); 
                std::cout << "\n";
            #endif
            break;
            
        }
	}
	// if ri went all the way through the end, then the quiver is not
	// ismomorph to any quiver in already in the list, so we add it
	if(ri==c.rend())
	{
		c.push_back(*quiver);
		#ifdef DEBUG
			std::cout << "Adding ";
			quiver->printMutations();
			std::cout << "\n";
		#endif
		
	}
}

/* Détection d'isomorphismes en utilisant nauty (algo polynomial de detection
   d'isomorphisme*/
bool MutExplorator::myIsomorphismNauty(Quiver *a, Quiver *b)
{
	int i;
	int n = a->getN();

	graph *c1;
	graph *c2;

	c1 = (graph *)b->oldGetNautyGraph();
	c2 = (graph *)a->oldGetNautyGraph();
	for(i=0;i<2*n;i++)
	{
		if(c1[i] !=c2[i])
		{
			return false;
		}
	}
	return true;
}
int MutExplorator::getModeCmp()
{
	return modeComparaison;
}

bool MutExplorator::isInC(Quiver *quiver)
{
	std::vector<Quiver>::reverse_iterator ri;
	quiver->genGraph();
	for(ri=c.rbegin();ri!=c.rend();ri++)
	{
		if(!ri->isGraphUpToDate())
			ri->genGraph();
		if(ri->getN()!=quiver->getN())
			break;
		else 
		{
			if(myIsomorphismNauty(quiver,&(*ri)))
			{
				#ifdef DEBUG
					quiver->printMutations(); 
					std::cout << "is isomorphic to  ";
					(*ri).printMutations(); 
					std::cout << "\n";
				#endif
				return true;
				
			}
		}
	}
	if(ri==c.rend())
	{
		return false;
		
	}
    return false;
}
