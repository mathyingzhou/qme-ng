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

 
#include "quiver.hpp"
#include "Exception.h"
/*
But: Constructeur
Entrée: un entier n, nombre de vertexs du graphe associé
Sortie: Néant
Précondition: Néant
PostCondition: les structures de données de l'objet sont allouées, la matrix d'incidence est initialisée à 0, pour tout i, pour tout j
*/
Quiver::Quiver(int n)
{
    int i,j;
    this->M=(int **)malloc(n*sizeof(int *));
    for(i=0;i<n;i++)
        (this->M)[i]=(int *)malloc(n*sizeof(int));
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
        {
            M[i][j]=0;
        }        
    this->n=n;
    this->graphAJour=0;
    semifree=0;
    nbNeighboursMax = -1;
    connected = -1;
    nextI=1;
    nextJ=0;
}

/*
But: Constructeur par recopie
Entrée: une référence sur un objet Quiver
Sortie: Néant
Précondition: Néant
PostCondition: les structures de données et les données de l'objet sont sont recopiées
*/

Quiver::Quiver(const Quiver &ca)
{
    int i,j;
    n=ca.n;
    if(ca.semifree == 0)
    {
        this->M=(int **)malloc(n*sizeof(int *));
        for(i=0;i<n;i++)
            (this->M)[i]=(int *)malloc(n*sizeof(int));
        for(i=0;i<n;i++)
            for(j=0;j<n;j++)
                M[i][j]=ca.M[i][j];
        this->semifree = 0;
    }
    else
    {
        this->semifree = 1;
    }
    this->graphAJour=ca.graphAJour;
    this->mutations = ca.mutations;
    this->score = ca.score;
    if(ca.graphAJour)
    {
        for(i=0;i<2*n;i++)
            nautyGC[i]=ca.nautyGC[i];
    }
    
    this->nbNeighboursMax = ca.nbNeighboursMax;
    this->connected = ca.connected;
    this->nextI = ca.nextI;
    this->nextJ = ca.nextJ;
    
}

/*
But: Constructeur par recopie et augmentation
Entrée: une référence sur un objet Quiver, un entier k
Sortie: Néant
Précondition: Néant
PostCondition: les structures de données et les données de l'objet sont recopiées dans un Quiver de taille n+k
    /!\     Le Quiver construit n'est pas connected !
        Il est nécessaire d'appeler connecter(i) pour obtenir un Quiver connected obtenu en reliant le vertex i à
            TOUS les autres.
*/

Quiver::Quiver(const Quiver &ca, int k)
{
    int i,j;
    n=ca.n+k;
    if(ca.semifree == 0)
    {
        this->M=(int **)malloc(n*sizeof(int *));
        for(i=0;i<n;i++)
            (this->M)[i]=(int *)malloc(n*sizeof(int));
        for(i=0;i<ca.n;i++)
            for(j=0;j<ca.n;j++)
                M[i][j]=ca.M[i][j];
        for(i=0;i<n;i++)
        {
            M[i][n-1]=0;
            M[n-1][i]=0;
        }
        this->semifree = 0;
    }
    else
    {
        this->semifree = 1;
    }

    this->graphAJour=0;
    semifree=0;
    nbNeighboursMax = -1;
    connected = -1;
    nextI=1;
    nextJ=0;
}

/*
But: Constructeur à partir d'une matrix
Entrée: une matrix
Sortie: Néant
Précondition: Néant
PostCondition: les structures de données et les données de l'objet sont initialisées
*/

Quiver::Quiver(int ** mat_quiver, int nbVertices, int indice)
{
    int i,j;
    n=nbVertices;
    this->M=(int **)malloc(n*sizeof(int *));
    
    for(i=0;i<n;i++)
        (this->M)[i]=(int *)malloc(n*sizeof(int));
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            M[i][j]=mat_quiver[i][j];
    //this->n=n; //Assign field to itself
    this->graphAJour=0;
    this->genScore();
    semifree=0;
    nbNeighboursMax = -1;
    connected = -1;
    nextI=1;
    nextJ=0;
}

/* But construire des quiver types
 *     Entrée: type
 *     A 0
    D 1
    E 2
    ATILDE 3
    DTILDE 4
    ETILDE 5
    SPORADIC 6
    UNAMED 7
    nbVertexs: le nombre de vertexs à construire
    Sortie: Néant
    Précondition: type d'un type défini cf ci-dessus et .h
    PostCondition: les structures de données et les données de l'objet sont désallouées
*/
Quiver::Quiver(int type, int nbVertices, int orientation)
{
    int i;
    
    // Dans tous les cas le graphe ne sera pas à jour
    this->graphAJour=0;
    this->semifree=0;
    switch(type)
    {
        case A:
            n=nbVertices;
            if(n<2)
                throw Exception("ERROR: not enough vertices");
            
            this->M=(int **)calloc(n,sizeof(int *));
            for(i=0;i<n;i++)
                (this->M)[i]=(int *)calloc(n,sizeof(int));
            for(i=0;i<n-1;i++)
            {
                M[i][i+1]=1;
                M[i+1][i]=-1;
            }
            
        break;
        case D:
            n=nbVertices;
            if(n<4)
                throw Exception("ERROR: not enough vertices");
            n=nbVertices;
            this->M=(int **)calloc(n,sizeof(int *));
            for(i=0;i<n;i++)
                (this->M)[i]=(int *)calloc(n,sizeof(int));
            for(i=0;i<n-2;i++)
            {
                M[i][i+1]=1;
                M[i+1][i]=-1;
            }
            M[n-3][n-1]=1;
            M[n-1][n-3]=-1;
            
        break;
        case E:
            n=nbVertices;
            switch(n)
            {
                case 6:
                case 7:
                case 8:
                    this->M=(int **)calloc(n,sizeof(int *));
                    for(i=0;i<n;i++)
                        (this->M)[i]=(int *)calloc(n,sizeof(int));
                    for(i=0;i<n-2;i++)
                    {
                        M[i][i+1]=1;
                        M[i+1][i]=-1;
                    }
                    M[n-4][n-1]=1;
                    M[n-1][n-4]=-1;
                break;
                default:
                    throw Exception("ERROR: bad vertex number asked");
            }    
        break;
        case ATILDE:
            n=nbVertices+orientation;
            if(n<2)
                throw Exception("ERROR: not enough vertices");
                
            this->M=(int **)calloc(n,sizeof(int *));
            for(i=0;i<n;i++)
                (this->M)[i]=(int *)calloc(n,sizeof(int));
            for(i=0;i<nbVertices;i++)
            {
                M[i][i+1]=1;
                M[i+1][i]=-1;
            }
            for(i=1;i<orientation;i++)
            {
                M[nbVertices+orientation-i][nbVertices+orientation-(i+1)]=1;
                M[nbVertices+orientation-(i+1)][nbVertices+orientation-i]=-1;
            }
            M[0][n-1] = 1;
            M[n-1][0]=-1;
            #ifdef DEBUG
            this->print();
            #endif
        break;
        case ATILDEALT:
            n=nbVertices+orientation;
            if(n<2)
                throw Exception("ERROR: not enough vertices");
                
            this->M=(int **)calloc(n,sizeof(int *));
            for(i=0;i<n;i++)
                (this->M)[i]=(int *)calloc(n,sizeof(int));
            for(i=0;i<n-1;i++)
            {
                M[i][i+1]=1-(2*(i%2));
                M[i+1][i]=-1+(2*(i%2));
            }
            M[n-1][0]=-1;
            M[0][n-1]=1;
            #ifdef DEBUG
            this->print();
            #endif
        break;
        case DTILDE:
            if(nbVertices<3)
                throw Exception("ERROR: not enough vertices");
            n=nbVertices+1;
            this->M=(int **)calloc(n,sizeof(int *));
            for(i=0;i<n;i++)
                (this->M)[i]=(int *)calloc(n,sizeof(int));
            for(i=1;i<n-2;i++)
            {
                M[i][i+1]=1;
                M[i+1][i]=-1;
            }
            M[n-3][n-1]=1;
            M[n-1][n-3]=-1;
            M[0][2]=1;
            M[2][0]=-1;
        break;
        case ETILDE:
            switch(nbVertices)
            {
                case 6:
                    n=nbVertices+1;
                    this->M=(int **)calloc(n,sizeof(int *));
                    for(i=0;i<n;i++)
                        (this->M)[i]=(int *)calloc(n,sizeof(int));

                    M[0][2] = -1;
                    M[1][3] = -1;
                    M[2][0] = 1;
                    M[2][5] = -1;
                    M[3][1] = 1;
                    M[3][5] = -1;
                    M[4][5] = 1;
                    M[4][6] = -1;
                    M[5][2] = 1;
                    M[5][3] = 1;
                    M[5][4] = -1;
                    M[6][4] = 1;                
                    break;
                case 7:
                    n=nbVertices+1;
                    this->M=(int **)calloc(n,sizeof(int *));
                    for(i=0;i<n;i++)
                        (this->M)[i]=(int *)calloc(n,sizeof(int));
                    
                    M[0][2] = -1;
                    M[1][3] = -1;
                    M[2][0] = 1;
                    M[2][4] = -1;
                    M[3][1] = 1;
                    M[3][5] = -1;
                    M[4][2] = 1;
                    M[4][7] = -1;
                    M[5][3] = 1;
                    M[5][7] = -1;
                    M[6][7] = 1;
                    M[7][4] = 1;
                    M[7][5] = 1;
                    M[7][6] = -1;
                    

                    break;
                case 8:
                    n=nbVertices+1;
                    this->M=(int **)calloc(n,sizeof(int *));
                    for(i=0;i<n;i++)
                        (this->M)[i]=(int *)calloc(n,sizeof(int));
                    M[0][2] = -1;
                    M[0][7] = 1;
                    M[1][3] = -1;
                    M[2][0] = 1;
                    M[2][4] = -1;
                    M[3][1] = 1;
                    M[3][6] = -1;
                    M[4][2] = 1;
                    M[4][6] = -1;
                    M[5][6] = 1;
                    M[6][3] = 1;
                    M[6][4] = 1;
                    M[6][5] = -1;
                    M[7][0] = -1;
                    M[7][8] = 1;
                    M[8][7] = -1;
                    
                break;
                default:
                    throw Exception("ERROR: bad vertex number asked");
            }
        break;
        case SPORADIC:
            n=nbVertices;
            if (n==3)
            {
                this->M=(int **)calloc(n,sizeof(int *));
                for(i=0;i<n;i++)
                    (this->M)[i]=(int *)calloc(n,sizeof(int));
                for(i=0;i<n-1;i++)
                {
                    M[i][i+1] = 2;
                    M[i+1][i] = -2;
                }
                M[n-1][0] = 2;
                M[0][n-1] = -2;
            }
            else if (n==4)
            {
                this->M=(int **)calloc(n,sizeof(int *));
                for(i=0;i<n;i++)
                    (this->M)[i]=(int *)calloc(n,sizeof(int));
                M[0][1] = -1;
                M[1][0] = 1;
                
                M[0][3] = -1;
                M[3][0] = 1;
                
                M[3][1] = 1;
                M[1][3]= -1;
                
                M[2][1] = 1;
                M[1][2] = -1;
                
                M[2][3] = 1;
                M[3][2] = -1;
                
                M[0][2] = 2;
                M[2][0] = -2;
                
            }
            else
            {
                throw Exception("ERROR: SPORADIC asked but n != 3 or 4");
            }
            
        break;
        case UNAMED:
            n=nbVertices;
            if (n>=5)
            {
                this->M=(int **)calloc(n,sizeof(int *));
                for(i=0;i<n;i++)
                    (this->M)[i]=(int *)calloc(n,sizeof(int));
                for(i=0;i<n-2;i++)
                {
                    M[i][i+1]=1;
                    M[i+1][i]=-1;
                }
                M[n-3][n-1]=1;
                M[n-1][n-3]=-1;
                M[0][1]=2;
                M[1][0]=-2;
            }
            else
                throw Exception("number of vertices too small");
        break;
        case  E_ELIPTIC:
            if(nbVertices>=6 && nbVertices<9)
            {
                n=nbVertices+2;
                this->M=(int **)calloc(n,sizeof(int *));
                for(i=0;i<n;i++)
                    (this->M)[i]=(int *)calloc(n,sizeof(int));
            }
            
            if(nbVertices==6)
            {
                M[0][1]=1;
                M[1][0]=-1;
                
                M[1][2]=1;
                M[2][1]=-1;
                            
                M[2][3]=2;
                M[3][2]=-2;
                
                M[1][3]=-1;
                M[3][1]=1;
                
                M[3][4]=1;
                M[4][3]=-1;
                
                M[5][3]=-1;
                M[3][5]=1;
                
                M[4][2]=1;
                M[2][4]=-1;
                
                M[5][6]=1;
                M[6][5]=-1;
                
                M[4][7]=1;
                M[7][4]=-1;
                
                M[5][2]=1;
                M[2][5]=-1;
    
                
            }
            else if(nbVertices==7)
            {
                
                M[0][1]=1;
                M[1][0]=-1;
                
                M[1][2]=1;
                M[2][1]=-1;
                            
                M[2][3]=2;
                M[3][2]=-2;
                
                M[1][3]=-1;
                M[3][1]=1;
                
                M[3][4]=1;
                M[4][3]=-1;
                
                M[5][3]=-1;
                M[3][5]=1;
                
                M[4][2]=1;
                M[2][4]=-1;
                
                M[5][6]=1;
                M[6][5]=-1;
                
                M[6][7]=1;
                M[7][6]=-1;
                
                M[5][2]=1;
                M[2][5]=-1;
                
                M[0][8]=1;
                M[8][0]=-1;
            }
            else if(nbVertices==8)
            {
                M[0][1]=1;
                M[1][0]=-1;
                M[1][2]=1;
                M[1][3]=-1;
                M[2][1]=-1;
                M[2][3]=2;
                M[2][4]=-1;
                M[2][5]=-1;
                M[3][1]=1;
                M[3][2]=-2;
                M[3][4]=1;
                M[3][5]=1;
                M[4][2]=1;
                M[4][3]=-1;
                M[5][2]=1;
                M[5][3]=-1;
                M[5][6]=1;
                M[6][5]=-1;
                M[6][7]=1;
                M[7][6]=-1;
                M[7][8]=1;
                M[8][7]=-1;
                M[8][9]=1;
                M[9][8]=-1;
            }
            else
            {
                throw Exception("Eliptique asked but wrong vertex number (must be 6,7 or 8)");
            }
            break;
        default:
            throw Exception("Bad type");
    } // end switch
    //Generate score
    this->genScore();
    nbNeighboursMax = -1;
    // This quiver is connected
    connected = 1;
    nextI=1;
    nextJ=0;

}
//Import quivers from files, both .qmu files and other files.
Quiver::Quiver(const char *file)
{
    std::string contents,line;
    std::ifstream f(file);
    boost::char_separator<char> sep(",[] \t;");
    std::vector<int> val;
    std::istringstream *iss;
    int i,j;
    unsigned int n;
    bool qmu = false;
    semifree=0;
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
        *iss >> i;
        delete iss;
        val.push_back(i);
    }
    f.close ();
    n=(unsigned int)sqrt(val.size());
    if(n*n != val.size())
    {
        throw Exception("Bad file format !");
    }
    
    this->M=(int **)malloc(n*sizeof(int *));
    for(i=0;i<n;i++)
    {
        (this->M)[i]=(int *)malloc(n*sizeof(int));
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            M[i][j]=val[i*n + j];
        }        
    }
    this->n=n;
    this->graphAJour=0;
    nbNeighboursMax = -1;
    connected = -1;
    nextI=1;
    nextJ=0;
}

/*
But: Destructeur par défaut
Entrée: Néant
Sortie: Néant
Précondition: Néant
PostCondition: les structures de données de l'objet sont désallouées
*/
//Destructor. Free the memory.
Quiver::~Quiver()
{
    int i;
    if(semifree==0)
    {
        for(i=0;i<this->n;i++)
            free((this->M)[i]);
        free(this->M);
    }
}

//Compute NbNeighboursMax, score and connectedness and then get rid of the matrix.
void Quiver::semiDestroy()
{
    int i;
    if(semifree==0)
    {
        getNbNeighboursMax(); // Penser à faire ça... sinon...
        genScore();
        isConnected();
        semifree=1;
        for(i=0;i<this->n;i++)
            free((this->M)[i]);
        free(this->M);
    }
}
//Default constructor. Doesn't do any initialization
Quiver::Quiver()
{

}

//Print the exchange matrix
void Quiver::print()
{
    int i,j;
    if(semifree==0)
    {
        for(i=0;i<this->n;i++)
        {
            for(j=0;j<this->n;j++)
            {
                std::cout << M[i][j] << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
    else
        std::cout << "The quiver has been semiFreed" << std::endl;
}

//Do mutation at k. Here counting starts at 0.
void Quiver::mutate(int k)
{
    int i,j;
    int lastMutatedVertex;
    // On ne fait rien si k ne correspond pas à un vertex du graphe
    if(k<0 || k>= this->n)
        return;
        
    //Création d'une matrix temporaire, copie de la matrix d'incidence
    /*int Mp[this->n][this->n];
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            Mp[i][j]=M[i][j];
    */        
    // Application de la fonction mu        
    
    for(i=0;i<n;i++)
    {
        if (i==k) continue;
        for(j=0;j<n;j++)
        {
            if (j==k) continue;
            M[i][j] = M[i][j] + (absVal(M[i][k])*M[k][j] + M[i][k]*absVal(M[k][j]))/2;
        }
    }
    for(i=0;i<n;i++)
    {
        M[i][k]=-M[i][k];
        M[k][i]=-M[k][i];
    }
    //If the exchange matrix is of infinite type then print an exception and getMutations().
    if(this->infinite())
        throw Exception("Mutation class is infinite ! " + getMutations());
    //Set graphAJour to 0.
    if(this->graphAJour)
    {
        this->graphAJour=0;
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
            //If the last mutation was done at k then doing another mutation cancels that mutation.
        }
        else
        {
            mutations.push_back(k);
        }
    }
    else
    {
        mutations.push_back(k);
    }
    //Compute score
    this->genScore();
}

//Take the absolute value of an integer

int Quiver::absVal(int k)
{
    if(k < 0)
        return -k;
    else
        return k;
}

//set M[i][j] to be val.
void Quiver::setM(int i, int j, int val)
{
    if(i<n && j < n && i>=0 && j>=0)
    {
        M[i][j]=val;
        //Restore graphAJour and connected to indeterminate after changing the exchange matrix
        if(this->graphAJour)
        {
            this->graphAJour = 0;
        }
        if(this->connected)
        {
            this->connected = -1;
        }
    }
    else
        throw Exception("DOMAIN_ERROR: setM");
    
    
}
/* Cette fonction génère les structures qui vont bien pour les appels à Nauty */
void Quiver::genGraph()
{
    int i,j,m,nbVertexsNauty;
    int lab1[MAXN],ptn[MAXN],orbits[MAXN];
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    setword workspace[5*MAXM];
    
    if(!this->graphAJour)
    {
        
        nbVertexsNauty = 2 * this->getN();
        m=(nbVertexsNauty + WORDSIZE - 1)/WORDSIZE;

        /* Si on trouve une valeur strictement positive dans la matrix d'incidence, alors on ajoute une arrête dans notre graphe */
        for(i=0;i<this->getN();i++)
        {
            gv=GRAPHROW(nautyG,i+this->getN(),m);
            EMPTYSET(gv,m);
            
            gv=GRAPHROW(nautyG,i,m);
            EMPTYSET(gv,m);
            /* On ajoute les fausses arrêtes entre le layer 0 et le layer 1 */
            ADDELEMENT(gv,i+this->getN());
            for(j=0;j<this->getN();j++)
            {
                /* multiplicité de 1 */
                if(this->getM(i,j)==1)
                {
                    gv=GRAPHROW(nautyG,i,m);
                    ADDELEMENT(gv,j);
                }
                else
                {
                    if(this->getM(i,j)==2)
                    {
                        gv=GRAPHROW(nautyG,i+this->getN(),m);
                        ADDELEMENT(gv,j+this->getN());
                    }
                }
            }
        }
        options.getcanon = TRUE;
        options.digraph = TRUE;
        options.defaultptn = FALSE;
        nauty_check(WORDSIZE,m,nbVertexsNauty,NAUTYVERSIONID);
        
        for(i=0;i<2*n;i++)
        {
            lab1[i]=i;
            ptn[i]=1;
        }
        ptn[n-1]=0;
        ptn[2*n-1]=0;
        
        
        nauty(nautyG,lab1,ptn,NULL,orbits,&options,&stats,
                                  workspace,5*MAXM,m,nbVertexsNauty,nautyGC);
        this->graphAJour=1;    
    }
}

/* Cette fonction teste deux propriétés du Quiver:
 *
 * - Si le quiver a une arrête au moins triple alors il n'est pas de mutation finie (Cas 0)
 * - Si un vertex x est connecté à un vertex v par une arrête double, et si v est
 *    aussi connecté par une arrête double à un autre vertex, alors le quiver n'est pas de
 *    mutation finie (Cas 1)
 *
 * Cette fonction renvoie vrai si l'un des deux cas présents ci-dessus est trouvé, et faux sinon.
 * Attention: ce n'est pas parceque cette fonction renvoie faux que le quiver est nécessairement de mutation finie !
 */
bool Quiver::infinite()
{
    int i,j,compteur;
    if(this->getN()>2)
    {
        for(i=0;i<this->getN();i++)
        {
            compteur = 0;
            for(j=0;j<this->getN();j++)
            {
                if(this->getM(i,j) >=2 || this->getM(i,j) <=-2)
                {
                    compteur++;
                }
                if(compteur == 2) /* Cas 1 le vertex i est connecté à 2 autres vertexs par une arrête double */
                {
                    return true;
                }
                if(this->getM(i,j) >=3) /* Cas 0 */
                {
                    return true;
                }
            }
        }
    } 
    return false;
}



graph *Quiver::getNautyGraph()
{
    if(!this->graphAJour)
    {
        this->genGraph();
    }
    return  (graph *)&nautyGC;
}

//Print the mutation sequence
//Ying Zhou's change: This time we want to start from 1.
void Quiver::printMutations()
{
    std::vector<int>::iterator i;
    if(mutations.empty()) std::cout << "-\n";
    else
    {
        for(i=mutations.begin();i!=mutations.end();i++)
        {
            std::cout << *i + 1;//In mathematics we want to start from 1.
            std::cout <<".";
        }
    std::cout << "\n";
    }
}

/* Calcule le "score" d'un quiver, utile pour sortir un "bon" représentant
 * du groupe de mutation
 */
void Quiver::genScore()
{
    int i,j;
    score=0;
    
    #ifdef SCORE1
    //Score 1 is literally the sum of all negative elements in the exchange matrix.
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            if (M[i][j]>0)
            {
                score-=M[i][j];
            }
    #else
    //Score Default (or Score 0) is literally 0 minus the count of all i j such that M[i][j]=2.
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            if (M[i][j]==2)
            {
                score-=1;
            }
    #endif
}

/* Effectue mutations mutations et regarde si le degré des arcs explose
 * 
 */
bool Quiver::testInfiniEmpirique(int mutations)
{
    int i;
    srand(time(NULL));
    // on travaille sur une copie
    Quiver t = *this;
    if(t.infinite())
        return true;
    for(i=0;i<mutations;i++)
    {
        t.mutate((int) (n * (rand() / (RAND_MAX + 1.0))));
        if(t.infinite())
            return true;
    }
    return false;
}
//Print the exchange matrix to a file.
void Quiver::toFile(const char* filename)
{
    int i,j;
    std::ofstream fichierSortie(filename);
    if(!fichierSortie)
        throw Exception("ERROR: cannot open output file !");
    fichierSortie << "[";
    for(i=0;i<this->getN();i++)
    {
        fichierSortie << "[";
        for(j=0;j<this->getN();j++)
        {
            fichierSortie << this->M[i][j];
            if(j!=this->getN() - 1)
                fichierSortie << ",";
        }
        fichierSortie << "]"  ;
        if(i!=this->getN()-1)
            fichierSortie << "," << std::endl;
    }
    fichierSortie << "]" << std::endl;
    fichierSortie.close();
}

//Get the mutation sequence as a string
//Note that here we start from 0.
std::string Quiver::getMutations()
{
    std::vector<int>::iterator i;
    std::stringstream out;

    if(mutations.empty()) out << "-\n";
    else
    {
        for(i=mutations.begin();i!=mutations.end();i++)
        {
            out << *i;
            out << ".";
        }
    out << "\n";
    }
    return out.str();
}

//Get the maximum number of vertices adjacent to any vertex.

int Quiver::getNbNeighboursMax()
{
    int i,j;
    int nbNeighboursTemp;
    if(nbNeighboursMax == -1)
    {
        for(i=0;i<n;i++)
        {
            nbNeighboursTemp=0;
            for(j=0;j<n;j++)
            {
                if(this->M[i][j] != 0)
                {
                    nbNeighboursTemp += 1;
                }    
            }
            if(nbNeighboursTemp > nbNeighboursMax)
            {
                nbNeighboursMax = nbNeighboursTemp;
            }
        }    

    }
    return nbNeighboursMax;
}
/* Warshall's algorithm. This calculates whether the quiver is connected.
 -1: connectedness is undetermined.
 0: the quiver is not connected.
 1: the quiver is connected.*/
int Quiver::isConnected()
{
    bool mat[n][n];
    int i,j,k;
    if(connected != -1)
        return connected;
    // On travaille sur une copie du graphe */
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            if(M[i][j] != 0)
                mat[i][j] = true;
            else
                mat[i][j] = false;
#ifdef DEBUG
    for(i=0;i<this->n;i++)
    {
        for(j=0;j<this->n;j++)
        {
            std::cout << mat[i][j] << "\t";
        }
        std::cout << "\n";
    }
    std::cout << "\n\n";
#endif

    for(i = 0; i < n; i++)
       for(j = 0; j < n; j++)
          for(k = 0; k < n; k++)
             mat[j][k] = mat[j][k] || (mat[j][i] && mat[i][k]);
#ifdef DEBUG
    for(i=0;i<this->n;i++)
    {
        for(j=0;j<this->n;j++)
        {
            std::cout << mat[i][j] << "\t";
        }
        std::cout << "\n";
    }
    std::cout << "\n\n";
#endif
/* Si la fermeture transitive est une clique, le graphe est connected */
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            if(mat[i][j]!=true)
            {
                connected = 0;
                return 0;
            }
    connected = 1;
    return 1;
}

/* But: Dire si un vertex i a une arrête double ou non
 * Entrée: i un vertex du quiver
 * Sortie: vraie si i est connecté par une arrête double, faux sinon
 * Précondition: i est compris entre 0 et n
 * PostCondition: néant
 */
/* Is vertex incident to some double edge?
 TODO: This method needs to be modified for valued quivers.*/
bool Quiver::aUneDouble(int vertex)
{
    int i;
    if(vertex < n && vertex >= 0)
    {
        for(i=0;i<n;i++)
        {
            if(M[i][vertex] == 2 || M[i][vertex] == -2)
                return true;
        }
        return false;
    }
    else
    {
        throw new Exception("aUneDouble: the argument is not a vertex");
    }
}

/*
But: Dire si trois vertexs du graphe forment un 3-cycle orienté ou non
Entrée: 3 entiers correspondant à 3 vertexs du quiver
Sortie: Vrai si les trois vertexs forment un 3-cycle orienté, Faux sinon.
Précondition: i, j et k sont compris entre 0 et n et sont tous les trois distincts
PostCondition: néant
*/
bool Quiver::isOrientedThreeCycle(int i, int j, int k)
{
    /* On vérifie que les entrées sont bien des vertexs du graphe */
    if(i<n && j < n && k < n && i>=0 && j>=0 && k>=0)
    {
        if(i != j && i != k && j != k)
        {
            if( M[i][j] != 0 && M[j][k] != 0 && M[k][i] != 0 )
            {
                if( (M[i][j] > 0 && M[j][k] > 0 && M[k][i] > 0) || (M[i][j] < 0 && M[j][k] < 0 && M[k][i] < 0) )
                    return true; // Toutes les arrêtes sont dans le même sens
                else
                    return false;
            }
            else
            {
                // Les trois vertexs ne forment pas un cycle !
                return false;
            }
        }
        else
        {
            // Les trois vertexs ne sont pas différents !
            return false;
        }
    }
    else
    {
        throw new Exception("ERROR, Quiver::isOrientedThreeCycle: One of the three arguments is not a vertex !");
    }
}

//Does the quiver contain oriented cycles?
bool Quiver::isCyclic()
{
    int visited[n];
    int i;
    for(i=0;i<n;i++)
    {
        visited[i] = 0;
    }
    for(i=0;i<n;i++)
    {
        if(visited[i] == 0) {
            if (exploreCycle(visited,i)) {
                return true;
            }
        }
    }
    return false;
}

bool Quiver::exploreCycle(int *visited,int i)
{
    int j;
    visited[i] = 1;
    for(j=0;j<n;j++)
    {
        if(M[i][j] > 0)
        {
            if(visited[j] == 1) {return true;}
            else 
            { 
                if (visited[j] == 0)
                {
                    if(exploreCycle(visited,j))
                    {
                        return true;
                    }
                }
            }
        }
    }
    visited[i] = 2;
    return false;
}

//get all neighbours of a given vertex
std::vector<int> Quiver::getNeighbours(int vertex)
{
    std::vector<int> neighbours;
    int i;
    for(i=0;i<n;i++)
    {
        if(M[i][vertex] != 0)
            neighbours.push_back(i);
    }
    return neighbours;
}

//get all neighbors with a (k,2) edge.
std::vector<int> Quiver::getNeighboursDoubles(int vertex)
{
    std::vector<int> neighbours;
    int i;
    for(i=0;i<n;i++)
    {
        if(M[i][vertex] == 2 || M[i][vertex] == -2)
            neighbours.push_back(i);
    }
    return neighbours;
}

//get all neighbors with a (k,1) edge.
std::vector<int> Quiver::getNeighboursSimples(int vertex)
{
    std::vector<int> neighbours;
    int i;
    for(i=0;i<n;i++)
    {
        if(M[i][vertex] == 1 || M[i][vertex] == -1)
            neighbours.push_back(i);
    }
    return neighbours;
}

//get all predecessors with a (k,1) edge
std::vector<int> Quiver::getNeighboursSimplesPredecesseurs(int vertex)
{
    std::vector<int> neighbours;
    int i;
    for(i=0;i<n;i++)
    {
        if(M[i][vertex] == -1)
            neighbours.push_back(i);
    }
    return neighbours;
}

//get all sucessors with a (k,1) edge
std::vector<int> Quiver::getNeighboursSimplesSuccesseurs(int vertex)
{
    std::vector<int> neighbours;
    int i;
    for(i=0;i<n;i++)
    {
        if(M[i][vertex] == 1)
            neighbours.push_back(i);
    }
    return neighbours;
}

//count the number of (k,1) predecessors
int Quiver::getNbNeighboursSimplesPredecesseurs(int vertex)
{
    int neighbours=0;
    int i;
    for(i=0;i<n;i++)
    {
        if(M[i][vertex] == -1)
            neighbours++;
    }
    return neighbours;
}

//count the number of (k,1) successors
int Quiver::getNbNeighboursSimplesSuccesseurs(int vertex)
{
    int neighbours=0;
    int i;
    for(i=0;i<n;i++)
    {
        if(M[i][vertex] == 1)
            neighbours++;
    }
    return neighbours;
}

//get all sinks of (k,2) arrows
std::vector<int> Quiver::getVertexsDoubleEdgeEntrante()
{
    std::vector<int> res;
    int i,j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            if(M[i][j]==2)
            {
                res.push_back(j);
                break;
            }    
        }
    }
    return res;
}
//get all vertices incident to (k,2) arrows
std::vector<int> Quiver::getVertexsPasDDoubleEdge()
{
    std::vector<int> res;
    int i,j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            if(M[i][j]==2 || M[i][j] == -2)
            {
                break;
            }    
        }
        res.push_back(i);
    }
    return res;
}

//Get the first vertex that is a source of a (k,2) arrow to i. If impossible return -1.
int Quiver::getVertexOrigineDoubleEdge(int i)
{
    int j;
    for(j=0;j<n;j++)
    {
        if(M[j][i] == 2)
            return j;
    }
    return -1;
}
