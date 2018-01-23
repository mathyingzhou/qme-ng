/*
 * Copyright (c) 2011-2012, Grégoire Dupont, Matthieu Pérotin
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

#include "iceQuiver.hpp"

IceQuiver::IceQuiver(const IceQuiver &ca)
{
    int i,j;
    n=ca.n;
    semiFreed = ca.semiFreed;
    mutationsSize = ca.mutationsSize;
    this->mutations = ca.mutations;
    if(!semiFreed) {
        this->M=new mpz_class[n*n];
        for(i=0;i<n;i++)
            for(j=0;j<n;j++)
                M[i*n+j]=mpz_class(ca.M[i*n+j]);
    }
    mutationString = ca.mutationString;
    this->greenVertices = ca.greenVertices;    
    this->graphAJour = ca.graphAJour;
    multiplicity=ca.multiplicity;
    if(ca.graphAJour)
    {
        this->nbVertexsNauty = ca.nbVertexsNauty;
        multiplicities = ca.multiplicities;
        for(i=0;i<this->nbVertexsNauty;i++)
            nautyGC[i]=ca.nautyGC[i];
    }
    
}

IceQuiver::IceQuiver(Quiver c)
{
    int i,j;
    this->n = 2 * c.getN();
    this->M=new mpz_class[n*n];
    for(i=0;i<n/2;i++)
    {
        for(j=0;j<n/2;j++)
        {
                M[i*n+j]=mpz_class(c.getM(i,j));
        }
        M[i*n+i+n/2] = mpz_class(1);
        M[(i+n/2)*n+i] = mpz_class(-1);
    }
    this->graphAJour = false;
    multiplicity[0] = 1;
    semiFreed = false;
    mutationsSize = 0;
}

IceQuiver::IceQuiver(const char *file)
{
    std::string contents,line;
    std::ifstream f(file);
    boost::char_separator<char> sep(",[] \t;");
    std::vector<mpz_class> val;
    std::istringstream *iss;
    int i,j;
    mpz_class m;
    unsigned int n;
    bool qmu = false;
    if(!f)
        throw Exception("ERROR: cannot open file !");
    std::getline (f, line);
    if(line == "//Number of points")
    {
        std::cout << "qmu file format detected !" << std::endl;
        qmu = true;
        while(line != "//Matrix")
            std::getline(f,line);
        std::getline(f,line);
        std::getline(f,line);
    }
    
    while(!f.eof())
    {
        
        contents+=line;
        std::getline(f,line);
        // Keller's new file format ends matrix definition by "Traffic lights"
        // the old file format ends it with "Points"
        // We keep both tests for retro compatibility
        if(qmu && (line == "//Traffic lights" || line == "//Points"))
        {
            break;
        }
        
    }
    tokenizer tokens(contents, sep);
    for (tokenizer::iterator tok_iter = tokens.begin();
         tok_iter != tokens.end(); ++tok_iter)
    {
        iss = new std::istringstream(*tok_iter);
        *iss >> m;
        delete iss;
        val.push_back(m);
    }
    f.close ();
    n=(unsigned int)sqrt(val.size());
    if(n*n != val.size())
    {
        throw Exception("Bad file format !");
    }
    
    this->M=new mpz_class[n*n];
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
                M[i*n+j]=val[i*n+j];
        }
    }
    this->n=n;
    this->graphAJour = false;
    multiplicity[0] = 1;
    semiFreed = false;
    mutationsSize = 0;
}

IceQuiver::~IceQuiver()
{
    if(!semiFreed) {
        delete[] this->M;
    }
}

void IceQuiver::semiDestroy()
{
    delete[] this->M;
    this->genGraph();
    multiplicity.clear();
    greenVertices.clear();
    mutationString = this->mutationsToString();
    mutations.clear();
    semiFreed = true;
}


/*
But: Appliquer la fonction de mutation sur un vertex du graphe
Entrée: en entier k correspondant à un des vertexs du graphe
Sortie: Néant
PréCondition: k est un vertex du graphe (=> k>=0 et k<n)
PostCondition: La fonction \mu_k est appliquée au quiver
*/

int IceQuiver::mutate(int k, mpz_class p)
{
    int i,j;
    int lastMutatedVertex;
    //int infinite = 0;
    // On ne fait rien si k ne correspond pas à un vertex du graphe
    if(k<0 || k>= this->n)
        return -1;
            
    // Application de la fonction mu        
    
    for(i=0;i<n;i++)
    {
        if (i==k) continue;
        for(j=0;j<n;j++)
        {
            if (j==k) continue;
            M[i*n+j] = M[i*n+j] + (abs(M[i*n+k])*M[k*n+j] + M[i*n+k]*abs(M[k*n+j]))/2;
        }
    }
    for(i=0;i<n;i++)
    {
        M[i*n+k]=-M[i*n+k];
        M[k*n+i]=-M[k*n+i];
    }

    if(p != 0 && this->infinite(p)) {
        //throw Exception("The mutation class is infinite ! " + getMutations());
        return 0;
    }

    /*
        On met à jour les mutations qui ont été appliquées sur le quiver
        Si la mutation appliquée est la même que la dernière qui avait été appliquée:
            alors on a appliqué deux fois la même, ce qui revient à ne pas l'appliquer, on l'efface de la liste
        Sinon
            on ajoute la mutation à la liste des mutations déjà appliquées
        
    */
    
    if(!mutations.empty())
    {
        lastMutatedVertex = mutations.back();
        if(lastMutatedVertex == k)
        {
            mutations.pop_back();
            this->unshiftMultiplicities();
            mutationsSize--;

        }
        else
        {
            mutations.push_back(k);
            this->shiftMultiplicities();
            mutationsSize++;
        }
    }
    else
    {
        mutations.push_back(k);
        this->shiftMultiplicities();
        mutationsSize++;
    }
    greenVertices.clear();
    this->generateGreenVertices();
    this->graphAJour = false;
    mutationString = "";
    if(greenVertices.size() == 0) { return 1;}
    else { return 2;}    
}

void IceQuiver::setM(int i, int j, mpz_class val)
{
    if(i<n && j < n && i>=0 && j>=0)
    {
        M[i*n+j]=val;
    }
    else
        throw Exception("DOMAIN_ERROR: setM");
    
    
}


bool IceQuiver::infinite(mpz_class p)
{
    int i,j;
    int N = this->getN();
    int n = N/2;
    if(this->getN()>2)
    {
        for(i=0;i<n;i++)
        {
            for(j=n;j<N;j++)
            {
                if(abs(M[i*N+j]) > p)
                {
                    return true;
                }
            }
        }
    } 
    return false;
}



void IceQuiver::printMutations(int s)
{
    std::vector<int>::iterator i;
    if(s == 0) {
        if(mutationString== "") {
            mutationString = this->mutationsToString();
        }
        std::cout << mutationString << std::endl;
    }
    else {
        std::cout << mutations.size() << std::endl;
    }
}

void IceQuiver::printMutationsE(int s)
{
    std::vector<int>::iterator i;
    if(s == 0) {
        if(mutations.empty()) std::cerr << "-" << std::endl;
        else
        {
            if(mutationString== "") {
                mutationString = this->mutationsToString();
            }
            std::cerr << mutationString << std::endl;
        }
    }
    else {
        std::cerr << mutations.size() << std::endl;
    }
}



void IceQuiver::generateGreenVertices()
{
    int i,j,c;
    greenVertices.clear();
    for(i=0;i<n/2;i++)
    {
        c=0;
        for(j=0;j<n/2;j++)
        {
            if(M[(n/2+j)*n+i] > 0) {
                c=1;
                break;
            }    
        }
        if(c==0) {
            greenVertices.push_back(i);
        }
    }
}

void IceQuiver::forceGreenVertex(int s)
{
    greenVertices.push_back(s);
}


int IceQuiver::getNextGreenVertex()
{
    int ret;
    if(greenVertices.size() == 0) { return -1;}
    #ifdef DEBUG
        std::vector<int>::iterator it;
        std::cout << "vertexsV: ";
        for(it=greenVertices.begin();it!=greenVertices.end();it++) { std::cout << (*it) << ",";}
        std::cout << std::endl;
    #endif
    ret = greenVertices.back();
    greenVertices.pop_back();
    return ret;
}

int IceQuiver::getRandomGreenVertex()
{
    int ret;
    int pos;
    if(greenVertices.size() == 0) { return -1;}
    #ifdef DEBUG
        std::vector<int>::iterator it;
        std::cout << "vertexsV: ";
        for(it=greenVertices.begin();it!=greenVertices.end();it++) { std::cout << (*it) << ",";}
        std::cout << std::endl;
    #endif
    if(greenVertices.size() == 0) { return -1;}
    pos = ((double)rand() / RAND_MAX)*(greenVertices.size());
    ret = greenVertices[pos];
    greenVertices.erase(greenVertices.begin()+pos);
    return ret;
}

Quiver *IceQuiver::getQuiver(void)
{
    return new Quiver(n);
}

/*
But: Afficher la matrix d'incidence sur la sortie standard
Entrée: Néant
Sortie: Néant
Précondition: Néant
PostCondition: La matrix d'incidence est affichée sur la sortie standard
*/
void IceQuiver::print()
{
    int i,j;
    for(i=0;i<this->n;i++)
    {
        for(j=0;j<this->n;j++)
        {
            std::cout << M[i*n+j] << "\t";
        }
        std::cout << "\n";
    }
    std::cout << "\n\n";
}

/* Cette fonction génère les structures qui vont bien pour les appels à Nauty */
void IceQuiver::genGraph()
{
    int i,j,m;
    int lab1[MAXN],ptn[MAXN],orbits[MAXN];
    std::map<mpz_class,mpz_class> multiplicities_index;
    std::map<mpz_class,mpz_class>::iterator mul_it;
    mpz_class nbSN_tmp = 0;

    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    setword workspace[5*MAXM];

    if(!this->graphAJour)
    {
    multiplicities.clear();
    // 1. Count multiplicities > 1 to get number of extra verticies
        for(i=0;i<this->getN();i++)
        {
            for(j=0;j<this->getN();j++)
            {
                if(this->getM(i,j)>1)
                {
                    multiplicities[this->getM(i,j)]+=1;
                }
            }
        }
    // 2. On with graph construction...

        nbVertexsNauty = 0;
        nbSN_tmp = 0;
        for(mul_it=multiplicities.begin();mul_it!=multiplicities.end();mul_it++)
        {
            multiplicities_index[mul_it->first]=nbVertexsNauty;
            nbSN_tmp += mul_it->second;
            if(nbSN_tmp.fits_sint_p()) {
                nbVertexsNauty = nbVertexsNauty + mul_it->second.get_si();
            }
            else
            {
                throw Exception("Wrap nbVertexsNauty!");
            }
        }

        nbSN_tmp += this->n;
        if(nbSN_tmp.fits_sint_p()) {
            nbVertexsNauty += this->n;
        }
        else
        {
            throw Exception("Wrap nbVertexsNauty!");
        }

        m=(nbVertexsNauty + WORDSIZE - 1)/WORDSIZE;

        /* Si on trouve une valeur strictement positive dans la matrix d'incidence, alors on ajoute une arrête dans notre graphe */
        for(i=0;i<nbVertexsNauty;i++)
        {
            gv=GRAPHROW(nautyG,i,m);
            EMPTYSET(gv,m);
        }
        for(i=0;i<nbVertexsNauty;i++)
        {
            lab1[i]=i;
            ptn[i]=1;
        }
        ptn[n/2-1]=0;
        ptn[n-1]=0;
        

        for(i=0;i<this->getN();i++)
        {
            /* On ajoute les fausses arrêtes entre le layer 0 et le layer 1 */
            for(j=0;j<this->getN();j++)
            {
                /* multiplicité de 1 */
                if(this->getM(i,j) <= 0) { continue;}
                if(this->getM(i,j)==1)
                {
                    gv=GRAPHROW(nautyG,i,m);
                    ADDELEMENT(gv,j);
                }
                else
                {
                    gv=GRAPHROW(nautyG,i,m);
                    nbSN_tmp = multiplicities_index[this->getM(i,j)] + this->n;
                    if(nbSN_tmp.fits_sint_p()) {
                        ADDELEMENT(gv,nbSN_tmp.get_si());
                    }
                    else
                    {
                        throw Exception("Wrap Vertex Nauty 1!");
                    }
                    gv=GRAPHROW(nautyG,nbSN_tmp.get_si(),m);
                    ADDELEMENT(gv,j);
                    multiplicities_index[this->getM(i,j)]++;
                }
            }
        }
        options.getcanon = TRUE;
        options.digraph = TRUE;
        options.defaultptn = FALSE;
        options.invarproc = adjacencies;
        options.mininvarlevel = 0;
        options.maxinvarlevel = 99;
        nauty_check(WORDSIZE,m,nbVertexsNauty,NAUTYVERSIONID);
        for(mul_it=multiplicities.begin();mul_it!=multiplicities.end();mul_it++)
        {
            nbSN_tmp = n-1+multiplicities_index[mul_it->first];
            if(nbSN_tmp.fits_sint_p()) {
                ptn[nbSN_tmp.get_si()] = 0;
            }
            else
            {
                throw Exception("Wrap Vertex Nauty 1!");
            }
        }
        
        

        nauty(nautyG,lab1,ptn,NULL,orbits,&options,&stats,
                                  workspace,5*MAXM,m,nbVertexsNauty,nautyGC);
        this->graphAJour=true;    
    }
}


graph *IceQuiver::getNautyGraph()
{
   // Quiver *quiver=NULL;
   // graph *g;
   // int i;
    if(!this->graphAJour)
    {
        this->genGraph();
        this->graphAJour=true;
    }
    return  (graph *)&nautyGC;
}

std::string IceQuiver::mutationsToString()
{
    std::string mutations_string;
    std::vector<int>::iterator i;
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    if(mutations.empty()) {ss << "-";}
    else
    {
        for(i=mutations.begin();i!=mutations.end();i++) {
            ss << *i+1 << " ";
        }
    }
    return ss.str(); 
}

void IceQuiver::addMultiplicity(IceQuiver &p)
{
    std::map<uint64_t,mpz_class> *mul = p.getMultiplicityMap();
    std::map<uint64_t,mpz_class>::iterator it;
    for(it=mul->begin();it!=mul->end();it++)
    {
        this->multiplicity[it->first] += it->second;
    }    
}

void IceQuiver::shiftMultiplicities()
{
    std::map<uint64_t,mpz_class>::reverse_iterator rit;
    std::map<uint64_t,mpz_class> nm;
    for(rit=multiplicity.rbegin();rit!=multiplicity.rend();rit++)
    {
       nm[rit->first + 1] = rit->second;
    }
    multiplicity = nm;
}

void IceQuiver::unshiftMultiplicities()
{
    std::map<uint64_t,mpz_class>::reverse_iterator rit;
    std::map<uint64_t,mpz_class> nm;
    for(rit=multiplicity.rbegin();rit!=multiplicity.rend();rit++)
    {
       nm[rit->first - 1] = rit->second;
    }
    multiplicity = nm;
}

void IceQuiver::toFile(const char* filename)
{
    int i,j;
    int n = this->getN();
    std::ofstream fichierSortie(filename);
    if(!fichierSortie)
        throw Exception("ERROR: cannot open output file !");
    fichierSortie << "[";
    for(i=0;i<n;i++)
    {
        fichierSortie << "[";
        for(j=0;j<n;j++)
        {
            fichierSortie << this->M[i*n+j];
            if(j!=n - 1)
                fichierSortie << ",";
        }
        fichierSortie << "]"  ;
        if(i!=n -1)
            fichierSortie << "," << std::endl;
    }
    fichierSortie << "]" << std::endl;
    fichierSortie.close();
}
