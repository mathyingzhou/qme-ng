/*
 * Copyright (c) 2011-2012, Grégoire Dupont, Matthieu Pérotin
 *               2018, Ying Zhou
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
    this->finiteGreenVertices = ca.finiteGreenVertices;
    this->graphIsUpToDate = ca.graphIsUpToDate;
    this->admittedCVectors = ca.admittedCVectors;
    this->admittedGVectors = ca.admittedGVectors;
    this->admittedPVectors = ca.admittedPVectors;
    multiplicity=ca.multiplicity;
    if(ca.graphIsUpToDate)
    {
        this->nbVerticesNauty = ca.nbVerticesNauty;
        multiplicities = ca.multiplicities;
        for(i=0;i<this->nbVerticesNauty;i++)
            nautyGC[i]=ca.nautyGC[i];
    }
    //std::cout << "Meow" << std::endl;
    is_c_known = ca.is_c_known;
    if (is_c_known) {
        cartan.setlength(n/2, n/2);
        euler.setlength(n/2, n/2);
        phi.setlength(n/2, n/2);
        phiin.setlength(n/2, n/2);
        cartan = ca.cartan;
        euler = ca.euler;
        phi = ca.phi;
        phiin = ca.phiin;
    }
}

IceQuiver::IceQuiver(Quiver c)
{
    int i,j;
    mat me, cmt, tcartan, tphi, tphiin;
    alglib::matinvreport rep;
    alglib::ae_int_t info;
    this->n = 2 * c.getN();
    this->M=new mpz_class[n*n];
    cartan.setlength(n/2, n/2);
    euler.setlength(n/2, n/2);
    phi.setlength(n/2, n/2);
    phiin.setlength(n/2, n/2);
    tphi.setlength(n/2, n/2);
    tphiin.setlength(n/2, n/2);
    tcartan.setlength(n/2, n/2);
    me.setlength(n/2, n/2);
    cmt.setlength(n/2, n/2);
    for(i=0;i<n/2;i++)
    {
        for(j=0;j<n/2;j++)
        {
            M[i*n+j]=mpz_class(c.getM(i,j));
            if (i == j) {
                cartan[i][j] = 1;
                me[i][j] = -1;
            }
            else {
                cartan[i][j] = c.getM(i, j) < 0 ? c.getM(i, j) : 0;
                me[i][j] = -cartan[i][j];
            }
        }
        M[i*n+i+n/2] = mpz_class(1);//Let's use the BDP convention on what "green" means.
        M[(i+n/2)*n+i] = mpz_class(-1);
    }
    this->graphIsUpToDate = false;
    multiplicity[0] = 1;
    semiFreed = false;
    mutationsSize = 0;
    euler = cartan;
    tcartan = rizem(cartan);
    alglib::rmatrixinverse(tcartan, n/2, info, rep);
    cartan = intizem(tcartan);
    alglib::rmatrixtranspose(n/2, n/2, tcartan, 0, 0, cmt, 0, 0);
    //std::cout << cmt.tostring(0) << me.tostring(0) << tphi.tostring(0) << std::endl;
    alglib::rmatrixgemm(n/2,n/2,n/2,1,cmt,0,0,0,me,0,0,0,0,tphi,0,0);
    phi = intizem(tphi);
    tphiin = tphi;
    alglib::rmatrixinverse(tphiin, n/2, info, rep);
    phiin = intizem(tphiin);
    is_c_known = true;
    std::cout << "Cartan matrix: " << cartan.tostring() << std::endl;
    std::cout << "Euler matrix: " << euler.tostring() << std::endl;
    std::cout << "Coxeter matrix: " << phi.tostring() << std::endl;
    std::cout << "Inverse of Coxeter matrix: " << phiin.tostring() << std::endl;
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
    this->graphIsUpToDate = false;
    multiplicity[0] = 1;
    semiFreed = false;
    mutationsSize = 0;
    is_c_known = false;
}

IceQuiver::~IceQuiver()
{
    if(!semiFreed) {
        delete[] this->M;
    }
}

//Generate the graph and the printing output. Destroy everything else
void IceQuiver::semiDestroy()
{
    delete[] this->M;
    this->genGraph();
    multiplicity.clear();
    greenVertices.clear();
    finiteGreenVertices.clear();
    admittedCVectors.clear();
    admittedGVectors.clear();
    mutationString = this->mutationsToString();
    mutations.clear();
    semiFreed = true;
    is_c_known = false;
}

/*
But: Appliquer la fonction de mutation sur un vertex du graphe
Entrée: en entier k correspondant à un des vertexs du graphe
Sortie: Néant
PréCondition: k est un vertex du graphe (=> k>=0 et k<n)
PostCondition: La fonction \mu_k est appliquée au quiver
*/

int IceQuiver::mutate(int k)
{
    int i,j;
    int lastMutatedVertex;
    ivect cVec, gVec, pVec;
    imat cMat, gMat, pMat;
    mat tcMat, tgMat, tpMat;
    alglib::real_2d_array temp;
    alglib::ae_int_t info;
    alglib::matinvreport rep;
    cVec.setlength(this->getN()/2);
    gVec.setlength(this->getN()/2);
    pVec.setlength(this->getN()/2);
    cMat.setlength(this->getN()/2, this->getN()/2);
    gMat.setlength(this->getN()/2, this->getN()/2);
    tgMat.setlength(this->getN()/2, this->getN()/2);
    pMat.setlength(this->getN()/2, this->getN()/2);
    tpMat.setlength(this->getN()/2, this->getN()/2);
    temp.setlength(this->getN()/2, this->getN()/2);
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
    
    for(i=0;i<n/2;i++)
    {
        for(j=0;j<n/2;j++) {
            cMat(i,j) = M[(i+n/2)*n+j].mpz_class::get_si();
            if (j == k) {
                cVec(i) = cMat(i,j);
            }
        }
    }
    //std::cout << "CMat" << std::endl;
    //std::cout << cMat.tostring() << std::endl;
    tcMat = rizem(cMat);
    alglib::rmatrixtranspose(this->getN()/2,this->getN()/2,tcMat,0,0,tgMat,0,0);
    alglib::rmatrixinverse(tgMat, info, rep);
    //alglib::rmatrixtranspose(this->getN()/2,this->getN()/2,gMat,0,0,temp,0,0);
    gMat = intizem(tgMat);
    //std::cout << "GMat" << std::endl;
    //std::cout << gMat.tostring() << std::endl;
    for (i = 0; i < n/2; i++) {
        gVec[i] = gMat(i, k);
    }
    //std::cout << "PMat" << std::endl;
    //std::cout << cartan.tostring(0) << std::endl;
    alglib::rmatrixgemm(n/2,n/2,n/2,1,rizem(cartan),0,0,0,tgMat,0,0,0,0,tpMat,0,0);
    pMat = intizem(tpMat);
    //std::cout << pMat.tostring() << std::endl;
    for (i = 0; i < n/2; i++) {
        pVec[i] = pMat(i, k) >= 0 ? pMat(i, k) : -pMat(i, k);
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
            this->unshiftMultiplicity();
            mutationsSize--;
            //TODO: Modify this part to accommodate red muts.
        }
        else
        {
            mutations.push_back(k);
            this->shiftMultiplicity();
            mutationsSize++;
            if (!belongsTo(cVec, admittedCVectors))
                admittedCVectors.push_back(cVec);
            if (!belongsTo(gVec, admittedGVectors))
                admittedGVectors.push_back(gVec);
            if (!belongsTo(pVec, admittedPVectors)) {
                //std::cout << "Meow" << std::endl;
                admittedPVectors.push_back(pVec);
            }
            
#ifdef DEBUG
            std::cout << "Printing c-vectors" << std::endl;
            printVector(cVec);
            std::cout << "Printing admitted c-vectors" << std::endl;
            printVectors(admittedCVectors);
            std::cout << "Printing g-vectors" << std::endl;
            printVector(gVec);
            std::cout << "Printing admitted g-vectors" << std::endl;
            printVectors(admittedGVectors);
            std::cout << "Printing p-vectors" << std::endl;
            printVector(pVec);
            std::cout << "Printing admitted p-vectors" << std::endl;
            printVectors(admittedPVectors);
#endif
        }
    }
    else
    {
        mutations.push_back(k);
        this->shiftMultiplicity();
        mutationsSize++;
        if (!belongsTo(cVec, admittedCVectors))
            admittedCVectors.push_back(cVec);
        if (!belongsTo(gVec, admittedGVectors))
            admittedGVectors.push_back(gVec);
        if (!belongsTo(pVec, admittedPVectors))
            admittedPVectors.push_back(pVec);
#ifdef DEBUG
        std::cout << "Printing c-vectors" << std::endl;
        printVector(cVec);
        std::cout << "Printing admitted c-vectors" << std::endl;
        printVectors(admittedCVectors);
        std::cout << "Printing g-vectors" << std::endl;
        printVector(gVec);
        std::cout << "Printing admitted g-vectors" << std::endl;
        printVectors(admittedGVectors);
        std::cout << "Printing p-vectors" << std::endl;
        printVector(pVec);
        std::cout << "Printing admitted p-vectors" << std::endl;
        printVectors(admittedPVectors);
#endif
    }
    this->generateGreenVertices();
    this->generateFiniteGreenVertices();
    this->graphIsUpToDate = false;
    mutationString = "";
    if(greenVertices.size() == 0) { return 1;}//1 = maximal green tail
    else if (finiteGreenVertices.size() == 0) {
#ifdef DEBUG
        std::cout << "Inf cut!" << std::endl;
        printMutationsE(0);
#endif
        return 0;
    }//0 = Banned vertices (BHIT & generalizations)
    else { return 2;} //No MGT yet
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

//This method is going to be deleted because it is completely absurd.
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

//If s is 0 print the mutation string. Otherwise print the length of the mutation sequence (no mutations = 0)
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

//If s is 0 print the mutation string. Otherwise print the length of the mutation sequence(no mutations = 0)
//There is only difference between this method and the method above, namely whether to print an empty mutation sequence at all.
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

//Actually used green sequences
void IceQuiver::generateFiniteGreenVertices()
{
    int i,j,c;
    finiteGreenVertices.clear();
    for(i=0;i<n/2;i++)
    {
        c=0;
        for(j=0;j<n/2;j++)
        {
            if(M[(n/2+j)*n+i] > 0 || (M[i*n+j] < 0 && M[i*n+j] * M[j*n+i] < -3)) {
                c=1;
                break;
            }    
        }
        if(c==0) {
            finiteGreenVertices.push_back(i);
        }
    }
#ifdef DEBUG
    std::vector<int>::iterator it;
    std::cout << "genFinGrV: Available finite green vertices: ";
    for(it=finiteGreenVertices.begin();it!=finiteGreenVertices.end();it++) { std::cout << (*it + 1) << ",";}
    std::cout << std::endl;
#endif
}

void IceQuiver::forceGreenVertex(int s)
{
    greenVertices.push_back(s);
}

//Get the last (WTF?) finite green vertex and remove it from the list of finite green vertices
int IceQuiver::getNextFiniteGreenVertex()
{
    int ret;
    if(finiteGreenVertices.size() == 0) { return -1;}
    #ifdef DEBUG
        std::vector<int>::iterator it;
        std::cout << "Available finite green vertices: ";
        for(it=finiteGreenVertices.begin();it!=finiteGreenVertices.end();it++) { std::cout << (*it + 1) << ",";}
        std::cout << std::endl;
    #endif
    ret = finiteGreenVertices.back();
    finiteGreenVertices.pop_back();
    greenVertices.pop_back();//greenVertices is not useful in all senses other than their size.
    return ret;
}

int IceQuiver::getRandomFiniteGreenVertex()
{
    int ret;
    int pos;
    if(finiteGreenVertices.size() == 0) { return -1;}
    #ifdef DEBUG
        std::vector<int>::iterator it;
        std::cout << "Available finite green vertices: ";
        for(it=finiteGreenVertices.begin();it!=finiteGreenVertices.end();it++) { std::cout << (*it) << ",";}
        std::cout << std::endl;
    #endif
    if(finiteGreenVertices.size() == 0) { return -1;}
    pos = ((double)rand() / RAND_MAX)*(finiteGreenVertices.size());
    ret = finiteGreenVertices[pos];
    finiteGreenVertices.erase(finiteGreenVertices.begin()+pos);
    greenVertices.erase(greenVertices.begin()+pos);
    return ret;
}

Quiver *IceQuiver::getQuiver(void)
{
    return new Quiver(n);
}

/*
Print the extended exchange matrix.
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

    if(!this->graphIsUpToDate)
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

        nbVerticesNauty = 0;
        nbSN_tmp = 0;
        for(mul_it=multiplicities.begin();mul_it!=multiplicities.end();mul_it++)
        {
            multiplicities_index[mul_it->first]=nbVerticesNauty;
            nbSN_tmp += mul_it->second;
            if(nbSN_tmp.fits_sint_p()) {
                nbVerticesNauty = nbVerticesNauty + mul_it->second.get_si();
            }
            else
            {
                throw Exception("Wrap nbVerticesNauty!");
            }
        }

        nbSN_tmp += this->n;
        if(nbSN_tmp.fits_sint_p()) {
            nbVerticesNauty += this->n;
        }
        else
        {
            throw Exception("Wrap nbVerticesNauty!");
        }

        m=(nbVerticesNauty + WORDSIZE - 1)/WORDSIZE;

        /* Si on trouve une valeur strictement positive dans la matrix d'incidence, alors on ajoute une arrête dans notre graphe */
        for(i=0;i<nbVerticesNauty;i++)
        {
            gv=GRAPHROW(nautyG,i,m);
            EMPTYSET(gv,m);
        }
        for(i=0;i<nbVerticesNauty;i++)
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
        nauty_check(WORDSIZE,m,nbVerticesNauty,NAUTYVERSIONID);
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
                                  workspace,5*MAXM,m,nbVerticesNauty,nautyGC);
        this->graphIsUpToDate=true;    
    }
}


graph *IceQuiver::oldGetNautyGraph()
{
   // Quiver *quiver=NULL;
   // graph *g;
   // int i;
    if(!this->graphIsUpToDate)
    {
        this->genGraph();
        this->graphIsUpToDate=true;
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
//Add multiplicity maps of *this and p considering the maps to be from uint64_t to Z with values at all the numbers that aren't mentioned 0.
void IceQuiver::addMultiplicity(IceQuiver &p)
{
    std::map<uint64_t,mpz_class> *mul = p.getMultiplicityMap();
    std::map<uint64_t,mpz_class>::iterator it;
    for(it=mul->begin();it!=mul->end();it++)
    {
        this->multiplicity[it->first] += it->second;
    }    
}
//I'm not sure what shift/unshiftMultiplicity are actually good for.
//Use g(x) = f(x-1)
void IceQuiver::shiftMultiplicity()
{
    std::map<uint64_t,mpz_class>::reverse_iterator rit;
    std::map<uint64_t,mpz_class> nm;
    for(rit=multiplicity.rbegin();rit!=multiplicity.rend();rit++)
    {
       nm[rit->first + 1] = rit->second;
    }
    multiplicity = nm;
}
//Use g(x) = f(x+1)
void IceQuiver::unshiftMultiplicity()
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
    std::ofstream outputFile(filename);
    if(!outputFile)
        throw Exception("ERROR: cannot open output file !");
    outputFile << "[";
    for(i=0;i<n;i++)
    {
        outputFile << "[";
        for(j=0;j<n;j++)
        {
            outputFile << this->M[i*n+j];
            if(j!=n - 1)
                outputFile << ",";
        }
        outputFile << "]"  ;
        if(i!=n -1)
            outputFile << "," << std::endl;
    }
    outputFile << "]" << std::endl;
    outputFile.close();
}

void IceQuiver::printMultiplicityMap() {
    std::map<uint64_t,mpz_class>::iterator it;
    std::cout << "Printing the Multiplicity Map of ";
    printMutationsE(0);
    for(it=multiplicity.begin();it!=multiplicity.end();it++) {
        std::cout << it->first << "," << it->second << std::endl;
    }
}
int IceQuiver::calcQuasiLength(ivect v) {
    int a = this->getN()/2;
    int max = a - 2;
    mat id = ident(a), temp;
    mat ph = rizem(phi);
    mat z;
    vect ll;
    int i = 2, j = 1;
    alglib::ae_int_t info;
    alglib::matinvreport rep;
    z.setlength(a, a);
    temp.setlength(a, a);
    ll.setlength(a);
    temp = id;
    if (!isreg(v)) {
        return 0;
    }
    //std::cout << "id " << id.tostring(3) << std::endl;
    //std::cout << "ph " << ph.tostring(3) << std::endl;
    for (;i <= max ;i++) {
        //std::cout << "temp " << temp.tostring(3) << std::endl;
        alglib::rmatrixgemm(a,a,a,1,temp,0,0,0,ph,0,0,0,1,id,0,0);
        temp = id;
        //std::cout << "i is " << i << std::endl;
        //std::cout << "temp " << temp.tostring(2) << std::endl;
        id = ident(a);
        z = rizem(intizem(temp));
        //std::cout << "z " << z.tostring(3) << std::endl;
        if (alglib::rmatrixdet(z) == 0) {
            //std::cout << "The mat is singular!" << std::endl;
            continue;
        }
            
        alglib::rmatrixinverse(z, a, info, rep);
        ll = matvecmult(z, rizev(v));
        //std::cout << "z " << z.tostring(3) << std::endl;
        //std::cout << "ll is" << ll.tostring(6) << std::endl;
        if (vecisint(ll) && ivecisnonneg(intizev(ll)) && isEulerOne(intizev(ll))) {
            j = i;
            //std::cout << "Get out at j = " << j << std::endl;
        }
    }
    return j;
}
int IceQuiver::calcEulerForm(ivect v) {
    imat temp = mult(vec2matr(v), euler);
    imat temp2 = mult(temp, vec2matc(v));
    return temp2[0][0];
}

int IceQuiver::vecProc(ivect v) {
    int p = prepdeg(v);
    int i = preideg(v);
    if (calcEulerForm(v) != 1) {
        std::cout << v.tostring() << " is not the dimension vector of any indecomposable rigid module!" << p << std::endl;
        return -1;
    }
    if(p >= 0) {
        std::cout << v.tostring() << " is preprojective of degree " << p << std::endl;
    }
    if (i >= 0) {
        std::cout << v.tostring() << " is preinjective of degree " << i << std::endl;
    }
    if (p < -1) {
        std::cout << v.tostring() << " is regular in standard stable tube of rank " << -p << " and has quasi-length " << calcQuasiLength(v) << std::endl;
        return 1;
    }
    if (p == -1 && i == -1) {
        std::cout << v.tostring() << " is regular of quasi-length " << calcQuasiLength(v);
        if (isComponentSincere(v))
            std::cout << " in a sincere component ";
        else
            std::cout << " in a nonsincere component ";
        std::cout << "with deg 0 quasi-simple " << getDegZeroQuasiSimple(v).tostring() << std::endl;
        return 2;
    }
    return 0;
}

void IceQuiver::vecProcs(vecivect vs) {
    vecivect::iterator vi, vj;
    vecivect reg_list;
    ivect qs;
    int i;
    for (vi = vs.begin(); vi != vs.end(); vi++) {
        i = vecProc(*vi);
        if (i == 2) {
            qs = getDegZeroQuasiSimple(*vi);
            for (vj = reg_list.begin(); vj != reg_list.end(); vj++) {
                if (sameMOrbit(phi, qs, *vj)) {
                    std::cout << qs.tostring() << " and " << vj->tostring() << " are in the same component with degree difference " << sameMOrbit(phi, *vj, qs) << " !" << std::endl;
                    break;
                }
            }
            if (vj == reg_list.end()) {
                reg_list.push_back(qs);
            }
        }
    }
    std::cout << "Component list" << std::endl;
    for (vj = reg_list.begin(); vj != reg_list.end(); vj++) {
        std::cout << vj->tostring() << std::endl;
    }
}

ivect IceQuiver::getQuasiSimpleQuot(ivect v) {
    int len = calcQuasiLength(v);
    alglib::ae_int_t info;
    alglib::matinvreport rep;
    imat mpom = mppo(phi,len - 1);
    mat tm = rizem(mpom);
    ivect qs;
    alglib::rmatrixinverse(tm, n/2, info, rep);
    qs = intizev(matvecmult(tm, rizev(v)));
    return qs;
}

/*int IceQuiver::calcTameM() {
    if (ab == )
}*/


