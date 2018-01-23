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

#include "mutexploratorSeq.hpp"
#include <sstream>

MutExploratorSeq::MutExploratorSeq()
{
	modeComparaison = 0;
	resume = 0;
}

MutExploratorSeq::~MutExploratorSeq()
{

}

int MutExploratorSeq::isomorphismExplorator(Quiver depart, int granularite)
{
	
	int stop=0;
	int infinite=0;
	Quiver ct = depart;
	//int i;
	if(resume == 0)
	{
		ct.testInfiniEmpirique(50*ct.getN());
		index = 0;
		depart.genGraph();
		c.push_back(depart);	
	}
	else
	{
		std::cout << "Resuming at pos: " << index << "\n";
		std::cout << "Nb Quivers in queue: " << c.size() << "\n";
		if(index == c.size())
		{
			return c.size();
		}
	}

	while(!stop)
	{
		while(!stop)
		{
			if(c[index].infinite())
			{
				stop=1; infinite=1;
				throw Exception("The mutation class is infinite" + c[index].getMutations());
			}
			
			#ifdef DEBUG
			std::cout << "================" << "\n";
			std::cout << "Position Index: " << index << "\n";
			std::cout << "================" << "\n";
			#endif
			generateMutations(c[index]);
			index++;
			if(index % granularite == 0)
				checkpoint();
			if(index == c.size())
				stop = 1;
		}
		stop = 1;
	}
	checkpoint();
	return c.size();
}

int MutExploratorSeq::isInTheMutationClass(Quiver test, Quiver depart)
{
	
	int stop=0;
	int infinite=0;
	if(resume == 0)
	{
		index = 1;
		modeComparaison = 1;
		if(test.getN() != depart.getN())
		{
			throw Exception("Les deux quiver n'ont pas le meme nombre de vertexs !");	
		}
		c.push_back(test);
		insertInList(&depart);
	}
	else
	{
		std::cout << "Resuming at pos: " << index << "\n";
		std::cout << "Nb Quivers in queue: " << c.size() << "\n";
		if(index == c.size())
		{
			return c.size()-1;
		}
	}
	while(!stop)
	{
		while(!stop)
		{
			if(c[index].infinite())
			{
				stop=1; infinite=1;
				throw Exception("The mutation class is infinite");
			}
			
			generateMutations(c[index]);
			index++;
			if(index % 5000 == 0)
				checkpoint();
			if(index == c.size())
				stop = 1;
		}
		stop = 1;
	}
	return c.size()-1;
	
}

int MutExploratorSeq::checkpoint()
{
	std::vector<Quiver>::iterator itC;
	std::ofstream ofs ("Out.txt");
	int i,j,n;

	std::cout << "Checkpointing at pos: " << index << "\n";
	std::cout << "Nb Quivers in queue: " << c.size() << "\n";

	ofs << index << "\n";
	itC=c.begin();
	n=(*itC).getN();
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			ofs << (*itC).getM(i,j) << "\t";
		}
		ofs << "\n";
	}
	
	if(modeComparaison == 1)
	{
		ofs << "--\n";
		itC++;
		for(i=0;i<n;i++)
		{
			for(j=0;j<n;j++)
			{
				ofs << (*itC).getM(i,j) << "\t";
			}
			ofs << "\n";
		}
	}
	ofs << "-\n";
	itC++;
	while(itC!=c.end())
	{
		ofs << (*itC).getMutations();
		itC++;
	}
	ofs << "I\n";
	ofs.close();
	return 1;
}

int MutExploratorSeq::dumpFiles(const char* prefix)
{
	std::vector<Quiver>::iterator itC;
	std::ofstream *ofs;
	std::stringstream *filename;
	int i,j,num,n;

	std::cout << "Dumping " << c.size() << " quivers !\n";

	itC=c.begin();
	n=(*itC).getN();
	num = 0;
	while(itC!=c.end())
	{
		filename = new std::stringstream;
		*filename << prefix << "_" << num++ << ".quiver"; 
		ofs = new std::ofstream((filename->str()).c_str()); 
		*ofs << "[";
		n=(*itC).getN();
		for(i=0;i<n;i++)
		{
			*ofs << "[";
			for(j=0;j<n;j++)
			{
				*ofs << (*itC).getM(i,j);
				if(j != (n-1)) { *ofs << " ";}
			}
			*ofs << "]";
			if(i != (n-1)) { *ofs << ",\n";}
		}
		*ofs << "]\n";
		itC++;
		ofs->close();
		delete ofs;
		delete filename;
	}
	return 1;
}

void MutExploratorSeq::reprise(const char* file)
{
	std::string contents,line;
	std::ifstream f(file);
	boost::char_separator<char> sep(",.[] \t;");
	std::vector<int> val;
	std::istringstream *iss;
	tokenizer::iterator tok_iter;
	tokenizer *tokens;
	int i,j;
	unsigned int n;
	Quiver *quiver, *c1;
	resume = 1;
	contents = "";
	std::cout << "Resuming from " << file << "...\n";
	if(!f)
		throw Exception("Impossible d'ouvrir le fichier !");
	
	// Retrieving index
	std::getline(f,line);
	iss = new std::istringstream(line);
	*iss >> index;	
	// Retrieving first Quiver
	std::getline (f, line);
	while(line != "-" && line != "--")
	{
		contents +=line;
		std::getline(f,line);
	}
	tokens = new tokenizer(contents, sep);
    	for (tokenizer::iterator tok_iter = tokens->begin();
         tok_iter != tokens->end(); ++tok_iter)
	{
		iss = new std::istringstream(*tok_iter);
		*iss >> i;
		delete iss;
		val.push_back(i);
	}
	delete tokens;
	n=(int)sqrt(val.size());
	if(n*n != val.size())
		throw Exception("Fichier mal formé !");

	quiver = new Quiver(n);
	for(i=0;i<(int)n;i++)
		for(j=0;j<(int)n;j++)
			quiver->setM(i,j,val[i*n+j]);
	quiver->print();
	quiver->genGraph();
	c.push_back(*quiver);
	// First Quiver retrieved
	if(line == "--")
	{
		delete quiver;
		std::cout << "Mode Comparaison...\n";
		modeComparaison = 1;
		std::getline(f,line);
		contents = "";
		val.clear();
		while(line != "-" && line != "--")
		{
			contents +=line;
			std::getline(f,line);
		}
		tokens = new tokenizer(contents, sep);
	    	for (tokenizer::iterator tok_iter = tokens->begin();
	         tok_iter != tokens->end(); ++tok_iter)
		{
			iss = new std::istringstream(*tok_iter);
			*iss >> i;
			delete iss;
			val.push_back(i);
		}
		delete tokens;
		n=(int)sqrt(val.size());
		if(n*n != val.size())
			throw Exception("Fichier mal formé !");

		quiver = new Quiver(n);
		for(i=0;i<(int)n;i++)
			for(j=0;j<(int)n;j++)
				quiver->setM(i,j,val[i*n+j]);
		quiver->print();
		quiver->genGraph();
		c.push_back(*quiver);
	}



	j=2;
	while(!f.eof())
	{
		std::getline(f,line);
		if(line == "I")
			break;
		j++;
		val.clear();
		tokens = new tokenizer(line, sep);
    		for (tok_iter = tokens->begin();tok_iter != tokens->end(); ++tok_iter)
		{
			iss = new std::istringstream(*tok_iter);
			*iss >> i;
			delete iss;
			val.push_back(i);
		}
		delete tokens;
		for(i=0;i<(int)val.size();i++)
		{
			quiver->mutate(val[i]);
		}
		c1 = new Quiver(*quiver);
		c1->genGraph();
		c.push_back(*c1);
		delete c1;
		for(i=val.size()-1;i>=0;i--)
		{
			quiver->mutate(val[i]);
		}
		
	}
	f.close ();
}

int MutExploratorSeq::getNbNeighboursMax()
{
	std::vector<Quiver>::iterator itC;
	int nbNeighboursMax=0;
	int valTemp;
	int i=0,index = 0;
	for(itC=c.begin();itC!=c.end();itC++)
	{
		valTemp = itC->getNbNeighboursMax();
		if( valTemp > nbNeighboursMax)
		{
			index = i;
			nbNeighboursMax = valTemp;
		}
		i++;
	}
	std::cout << "nbNeighbours Max: " << nbNeighboursMax << " found in " << c[index].getMutations() << "\n";
	return nbNeighboursMax;
	
}

bool MutExploratorSeq::isAcyclic()
{
	std::vector<Quiver>::iterator itC;
	for(itC=c.begin();itC!=c.end();itC++)
	{
		if(itC->isCyclic() == false)
		{
			std::cout << "No Cycles in " << itC->getMutations() << "\n";
			return true;
		}
	}
	return false;
}
