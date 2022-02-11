#include <bits/stdc++.h>

using namespace std;

class Graf
{
protected:
    int n;
    int m;
    ///structuri de retinere a muchiilor
    vector<vector<int>> muchii; //liste de adiacenta- array de vectori
    vector<tuple<int,int,int>> muchii_costuri_lista; //array de muchii cu costuri (e1,e2,cost)
    vector<vector<pair<int,int>>> muchii_costuri_adiacenta; //liste de adiacenta- array de vectori pentru muchii cu costuri (nod,cost)
    vector<vector<int>> matrice_ponderi;

protected:
    int Reprezentant(int x, vector<int>& rad);

public:
    Graf(int noduri, int muchii):n(noduri), m(muchii) {}
    void CitireMuchiiCosturiLista(ifstream& fin);

    static bool HavelHakimi(vector<int>& grade);

    //pentru paduri de multimi disjuncte
    void Reuneste(int x, int y, vector<int>& h, vector<int>& rad);
    bool Verifica(int x, int y, vector<int>& rad);

    int CompConexe();
    vector<int> BFS(int start); //returneaza distantele de la start la fiecare nod
    void DFS(int start, vector<vector<int>>& muchii_curente, int viz[]);
    vector<int> Dijkstra(int s);
};

class GrafOrientat: public Graf
{
private:
    void DFSTimpi(int nod,vector<vector<int>>& muchii_curente, int viz[], stack<int>& timpi_final); //memoreaza si timpii de final
    void DFSCtc (int nod, vector<vector<int>>& muchii_curente, int viz[], int nr_comp, vector<vector<int>>& comp_tare_con);//memoreaza si comp tare conexe
    vector<vector<int>> GrafTranspus();
    bool ExistaLantNesaturat(int s,int d, vector<int>& tata, vector<vector<int>>& muchii,
                             vector<vector<int>>& c, vector<vector<int>>& f, vector<int>& viz);

public:
    GrafOrientat(int noduri, int muchii):Graf(noduri,muchii) {}
    void CitireMuchiiAdiacenta(ifstream& fin);
    void CitireMuchiiCosturiAdiacenta(ifstream& fin);
    void CitireMatricePonderi(ifstream& fin);
    void CitireMuchiiCosturiAdiacentaInvers(ifstream& fin);

    pair<vector<vector<int>>, int> CTC(); //alg. lui Kosaraju
    vector<int> SortareTopologica();
    vector<int> BellmanFord(int s);
    vector<vector<int>> RoyFloyd();
    int FluxMaxim(int s, int d);
    int CicluHamiltonian();
};

class GrafNeorientat: public Graf
{
private:
    vector<pair<int,int>> capete_muchii;
    vector<vector<int>> muchii_adiacenta_indici; //fiecare nod va avea o lista cu indicii muchiilor de care apartine

    void DFSCritice(int nod, vector<vector<int>>& muchii, int viz[], vector<pair<int,int>>& muchii_critice, int nivel[], int nivel_min[] );
    void DFSBiconex(int nod, vector<vector<int>>& muchii, int viz[], int nivel[], int nivel_min[], stack<int>& noduri_traversate, vector<vector<int>>& componente);
    void DFSDiametru(int nod, vector<vector<int>>& muchii, int viz[], int distante[]);
    bool DFSCuplaj(int nod, vector<int>& viz, vector<int>& pereche_stanga, vector<int>& pereche_dreapta);

public:
    GrafNeorientat(int noduri, int muchii):Graf(noduri,muchii) {}
    void CitireMuchiiAdiacenta(ifstream& fin);
    void CitireMuchiiCosturiAdiacenta(ifstream& fin);
    void CitireMuchii(ifstream& fin); //fiecare muchie are un indice (se retine sursa si destinatia ei)
    void CitireMuchiiAdiacentaCuplaj(ifstream& fin, int nr1);

    vector<pair<int,int>> MuchiiCritice();
    vector<vector<int>> ComponenteBiconexe();
    vector<tuple<int,int,int>> APM(); //alg. lui Kruskal
    int DiametruArbore();
    vector<int> CicluEulerian();
    vector<pair<int,int>> CuplajMaxim(int nr1, int nr2);
};

bool GrafNeorientat::DFSCuplaj(int nod, vector<int>& viz, vector<int>& pereche_stanga, vector<int>& pereche_dreapta)
{
    if(viz[nod])
        return 0;

    viz[nod]=1;
    for(int i=0; i<(int)muchii[nod].size(); ++i) //iau vecinii nodului curent (noduri din dreapta)
    {
        int vecin=muchii[nod][i];
        if(!pereche_stanga[vecin] || DFSCuplaj(pereche_stanga[vecin], viz, pereche_stanga, pereche_dreapta))
            //daca nodul din dreapta nu are pereche in stanga sau are pereche in stanga care continua cu lant
        {
            pereche_stanga[vecin]=nod;
            pereche_dreapta[nod]=vecin;
            return 1; //s-a putut forma lant alternant
        }
    }
    return 0; //nu s-a gasit lant
}

vector<pair<int,int>> GrafNeorientat::CuplajMaxim(int nr1, int nr2)
{
    ///Algoritmul Hopcroft-Karp
    ///Complexitate: O(rad(n)*m), m=muchii, n=noduri
    //nr1= cardinal multime stanga, nr2= cardinal multime dreapta
    vector<pair<int,int>> sol;
    vector<int> pereche_dreapta(nr1+1, 0), pereche_stanga(nr2+1,0);  //retin perechile nodurilor din stanga (perechea dreapta),
    //respectiv din dreapta (perechea stanga); 0=fara pereche

    vector<int> viz(nr1+1); //retin daca am trecut prin acest nod in parcurgere (doar noduri din stanga)
    bool lant_gasit=1;

    //actualizez rezultatul pana nu mai pot adauga lanturi noi
    while(lant_gasit) //cat timp mai pot gasi drumuri("lanturi alternante") care sa plece din nod fara pereche din stanga,
        //sa ajunga tot in nod fara pereche din stanga
        //si sa treaca prin orice fel de nod
    {
        for(int i=0; i<=nr1; ++i)
            viz[i]=0;

        lant_gasit=0; //presupun ca nu mai exista lant de adaugat
        for(int i=1; i<=nr1; ++i)
            if(!pereche_dreapta[i] && DFSCuplaj(i, viz, pereche_stanga, pereche_dreapta))
                //iau toate nodurile din stanga care nu au pereche in dreapta si pot forma un lant alternant pornind din ele
                lant_gasit=1;
    }

    //formez solutia: iau nodurile din stanga care au pereche in dreapta
    for(int i=1; i<=nr1; ++i)
        if(pereche_dreapta[i])
            sol.push_back(make_pair(i, pereche_dreapta[i]));

    return sol;
}

int GrafOrientat::CicluHamiltonian()
{
    ///Programare dinamica, complexitate: O(m*2^n)
    int putere=1<<n; //2^n, toate nodurile=putere-1
    int matrice_cost[putere][n];
    //matrice_cost[i][j]= costul minim in lantul 0-j care contine exact nodurile marcate cu 1
    //in reprezentarea binara a lui i

    const int inf=100000000;
    int sol=inf;

    for(int i=0; i<putere; ++i)
        for(int j=0; j<n; ++j)
            matrice_cost[i][j]=inf;

    //consider ciclul ca pornind din 0 (un nod random)
    matrice_cost[1][0]=0; //lantul cu nodul 0 are costul 0

    //formez lanturi de cost minim de la 0 la 1,2,...,n-1
    for(int i=0; i<putere; ++i)//am lantul format din nodurile din reprez. lui i
        for(int j=0; j<n; ++j)//caut nodurile din lantul format
            if((1<<j) & i) //verific daca nodul j se afla printre cele din lantul format
            {
                for(int k=0; k<(int)muchii_costuri_adiacenta[j].size(); ++k) //iau vecinii care intra in nodul j
                {
                    int vecin=muchii_costuri_adiacenta[j][k].first;
                    int cost=muchii_costuri_adiacenta[j][k].second;

                    if((1<<vecin) & i) //verific daca nodul vecin se afla deja printre cele din lantul format
                    {
                        //actualizez costul minim curent
                        //prin a verifica daca folosirea nodului 'vecin' ca intermediar e profitabila

                        //compar costul precedent cu cel obtinut prin eliminarea nodului j din lantul care se termina in 'vecin'
                        //+costul muchiei de la 'vecin' la j
                        matrice_cost[i][j]= min(matrice_cost[i][j], matrice_cost[i & ~(1<<j)][vecin]+cost);
                    }
                }
            }

    for(int i=0; i<(int)muchii_costuri_adiacenta[0].size(); ++i)
        //iau nodurile care intra in 0 (care inchid ciclul) si aflu lungimea minima a lanturilor 0->nod + cost nod->0
    {
        int nod=muchii_costuri_adiacenta[0][i].first;
        int cost=muchii_costuri_adiacenta[0][i].second;
        sol=min(sol,matrice_cost[putere-1][nod]+cost);
    }

    if(sol!=inf)
        return sol;
    return 0;
}

vector<int> GrafNeorientat::CicluEulerian()
{
    ///Plec de la ideea ca orice graf Eulerian poate fi descompus ca reuniune de cicluri disjuncte.
    ///Construiesc un ciclu oarecare si intercalez recursiv ciclurile formate de muchiile ramase in graf.
    ///Teorema Leonhard Euler: un multigraf este Eulerian (admite un ciclu Eulerian) daca si numai daca este conex si toate nodurile sale au grad par

    vector<int> sol;
    vector<bool> eliminata(m,0); //pentru muchii
    stack<int> S;

    for(int i=1; i<=n; ++i)
        if(muchii_adiacenta_indici[i].size()%2==1
                || muchii_adiacenta_indici[i].size()==0) //nod izolat/cu grad impar=>nu poate fi graf eulerian
        {
            sol.push_back(-1);
            return sol;
        }

    //sigur este eulerian
    S.push(1);//plec din nodul 1

    while(!S.empty())
    {
        int nod=S.top();

        if(!muchii_adiacenta_indici[nod].empty()) //mai are vecini
        {
            //iau un vecin oarecare
            int muchie=muchii_adiacenta_indici[nod].back(); //iau indicele ultimei muchii ca sa elimin usor
            muchii_adiacenta_indici[nod].pop_back();

            if(!eliminata[muchie]) //marchez muchiile eliminate ca sa stiu atunci cand vin din celalalt capat al ei
            {
                eliminata[muchie]=1;
                //adaug celalalt capat in stiva
                if(nod==capete_muchii[muchie].first)
                    S.push(capete_muchii[muchie].second);
                else
                    S.push(capete_muchii[muchie].first);

            }
        }
        else
        {
            sol.push_back(nod);
            S.pop();
        }
    }

    return sol;
}

vector<vector<int>> GrafOrientat::RoyFloyd()
{
    ///complexitate: O(n^3)

    vector<vector<int>> drumuri_minime;
    const int valmax=100000000; //folosita la initializare in locurile unde nu exista muchie pentru a nu folosi costul 0
    drumuri_minime.resize(n);

    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            if(matrice_ponderi[i][j])
                drumuri_minime[i].push_back(matrice_ponderi[i][j]);
            else
                drumuri_minime[i].push_back(valmax);

    //pentru orice pereche de noduri (i,j) incerc sa optimizez drumul folosind un intermediar k
    for(int k=0; k<n; ++k)
        for(int i=0; i<n; ++i)
            for(int j=0; j<n; ++j)
                if(drumuri_minime[i][j]>drumuri_minime[i][k]+drumuri_minime[k][j])
                    drumuri_minime[i][j]=drumuri_minime[i][k]+drumuri_minime[k][j];

    return drumuri_minime;
}

void GrafNeorientat::DFSDiametru(int nod, vector<vector<int>>& muchii_curente, int viz[], int dist[])
{
    ///complexitate: O(n+m)
    viz[nod]=1;
    for(int i=0; i<(int)muchii_curente[nod].size(); ++i)
        if(!viz[muchii_curente[nod][i]])
        {
            dist[muchii_curente[nod][i]]=dist[nod]+1;
            DFSDiametru(muchii_curente[nod][i], muchii_curente, viz,dist);
        }
}

int GrafNeorientat::DiametruArbore()
{
    ///complexitate: 2*DFS cu m=n-1 => O(n)

    int dist[n+1];
    int viz[n+1];
    int maxim=-1, nod=1;

    for(int i=1; i<=n; ++i)
        viz[i]=0;

    //incep din nodul 1 (un nod random) sa fac dfs pentru a gasi cea mai indepartata frunza
    dist[1]=1;
    DFSDiametru(1, muchii, viz, dist);

    for(int i=1; i<=n; ++i)
        if(dist[i]>maxim)
            maxim=dist[i], nod=i;

    //continui cu un dfs din nodul cel mai indepartat gasit pentru a gasi celalalt capat al lantului maxim
    for(int i=1; i<=n; ++i)
        viz[i]=0;

    dist[nod]=1;
    DFSDiametru(nod, muchii, viz, dist);

    maxim=-1;
    for(int i=1; i<=n; ++i)
        if(dist[i]>maxim)
            maxim=dist[i];

    return maxim;

}
vector<int> GrafOrientat::BellmanFord(int s)
{
    ///complexitate: O(n*m), cu optimizare cu coada (Shortest Path Faster Algorithm)
    //folosind muchii_costuri_adiacenta

    bool in_coada[n+1];
    bool ok=1; //1-nu are cicluri negative, 0-are

    vector<int> d;
    int nr_incr[n+1]; //retine pentru fiecare nod numarul de incrementari pentru a detecta cicluri negative

    const int inf=250000005; //m*costmax muchie
    queue<int> noduri_modificate; //nodurile pentru care are sens sa incerc relaxarea vecinilor
    //(distantele acestora s-au modificat la pasul anterior)

    d.resize(n+1);

    for(int i=1; i<=n; ++i)
    {
        nr_incr[i]=0;
        in_coada[i]=0;
        if(i!=s)
            d[i]=inf;
        else d[i]=0;
    }

    noduri_modificate.push(s); //consider sursa ca fiind importanta la inceput
    in_coada[s]=1;

    //iau fiecare muchie care porneste dintr-un nod care s-a modificat anterior
    while(!noduri_modificate.empty() && ok)
    {
        int nod1=noduri_modificate.front();
        noduri_modificate.pop();
        in_coada[nod1]=0;

        for(int i=0; i<(int)muchii_costuri_adiacenta[nod1].size(); ++i)
        {
            int nod2=get<0>(muchii_costuri_adiacenta[nod1][i]);
            int cost=get<1>(muchii_costuri_adiacenta[nod1][i]);
            //am acum muchia nod1-nod2 cu costul ei

            if(d[nod2] > d[nod1]+cost) //pot relaxa muchia
            {
                d[nod2]=d[nod1]+cost;
                nr_incr[nod2]++;

                if(!in_coada[nod2]) //daca nu era in coada nodurilor importante (care s-au actualizat), il adaug
                {
                    noduri_modificate.push(nod2);
                    in_coada[nod2]=1;

                }
                if(nr_incr[nod2]>=n) //apare ciclu negativ pentru ca sunt prea multe incrementari (bucla infinita)
                {
                    ok=0;
                    d.clear();
                    break;
                }
            }
        }
    }
    return d;
}

vector<int> Graf::Dijkstra(int s)
{
    ///complexitate: O(m * logn)
    //folosind muchii_costuri_adiacenta

    bool viz[n+1]; //pentru a marca nodurile prin care trecem
    vector<int> d; //distanta de la s la celelalte noduri
    const int inf = 250005;

    //min heap cu sortare dupa distanta (primul element)
    priority_queue <pair<int,int>, vector<pair<int,int>>, greater<pair<int, int>>> dist_nod;

    //initializarea
    d.resize(n+1);
    for(int i=1; i<=n; ++i)
    {
        viz[i]=0;

        if(i!=s)
            d[i]=inf;
        else
        {
            d[i]=0;
            dist_nod.push(make_pair(d[i],i)); //adaug doar nodul sursa in heap initial
        }
    }

    while (!dist_nod.empty())
    {
        int u=get<1>(dist_nod.top()); //extrage nodul cu eticheta minima (distanta minima)
        dist_nod.pop();

        if(!viz[u])
        {
            viz[u]=1;

            for(int i=0; i<(int)muchii_costuri_adiacenta[u].size(); ++i) //caut nodurile v pentru care exista muchia uv
            {
                int v,cost;

                v=get<0>(muchii_costuri_adiacenta[u][i]);
                cost=get<1>(muchii_costuri_adiacenta[u][i]); //dintre u si v

                //verific daca pot relaxa muchia
                if(d[v] > d[u] + cost)
                {
                    d[v]=d[u]+cost; //actualizez distanta
                    dist_nod.push(make_pair(d[v],v));

                }

            }
        }
    }
    return d;
}

void GrafNeorientat::DFSBiconex(int nod, vector<vector<int>>& muchii, int viz[], int nivel[], int nivel_min[], stack<int>& noduri_traversate, vector<vector<int>>& componente)
{
    int fiu;

    viz[nod]=1;
    noduri_traversate.push(nod);

    nivel_min[nod]=nivel[nod];

    for(int j=0; j<(int)muchii[nod].size(); ++j)
    {
        fiu=muchii[nod][j];

        if(!viz[fiu])
        {
            nivel[fiu]=nivel[nod]+1;
            DFSBiconex(fiu, muchii, viz, nivel,nivel_min,noduri_traversate,componente);

            //formula A de la muchii critice
            //actualizez nivelul minim al nodului curent deoarece fiul poate urca la un nod de deasupra celui curent
            nivel_min[nod]=min(nivel_min[nod], nivel_min[fiu]);

            //verific daca parintele e nod critic
            if(nivel_min[fiu]>=nivel[nod])
            {
                vector<int> temp;
                //elimin din stiva pana la fiu, apoi fiul si nodul curent
                //nu se elimina pana la nodul curent direct deoarece se mai pot afla noduri
                //de la alta componenta biconexa intre cel curent si fiu

                while(noduri_traversate.top()!=fiu)
                {
                    temp.push_back(noduri_traversate.top());
                    noduri_traversate.pop();
                }
                temp.push_back(fiu);
                noduri_traversate.pop();

                temp.push_back(nod);
                //nu elimin nodul tata deoarece el poate face parte si din alta componenta biconexa

                componente.push_back(temp);

            }

        }
        else //muchie de intoarcere
            if(nivel[fiu]<nivel[nod]-1)
                //formula B de la muchii critice
                nivel_min[nod]=min(nivel_min[nod], nivel[fiu]);

    }
}

vector<vector<int>> GrafNeorientat::ComponenteBiconexe()
{
    ///asemanator cu muchiile critice
    ///complexitate: O(n+m)

    vector<vector<int>> componente;
    stack<int> noduri_traversate;
    int nivel[n+1];
    int nivel_min[n+1];
    int viz[n+1];

    for(int i=1; i<=n; ++i)
        viz[i]=0;

    nivel[1]=1;

    DFSBiconex(1,muchii,viz,nivel,nivel_min,noduri_traversate,componente);

    return componente;

}

void GrafNeorientat::DFSCritice(int nod, vector<vector<int>>& muchii,int viz[], vector<pair<int,int>>& muchii_critice, int nivel[], int nivel_min[])
{
    int fiu;

    viz[nod]=1;
    //setez nivelul minim ca fiind cel curent
    nivel_min[nod]=nivel[nod];

    for(int j=0; j<(int)muchii[nod].size(); ++j)
    {
        fiu=muchii[nod][j];

        if(!viz[fiu])
        {
            nivel[fiu]=nivel[nod]+1; //actualizez nivelul fiului
            DFSCritice(fiu, muchii,viz, muchii_critice,nivel,nivel_min);

            //la intoarcerea din dfs, actualizez dupa formula A (curs)
            nivel_min[nod]=min(nivel_min[nod], nivel_min[fiu]);

            //verific daca e muchie critica
            if(nivel_min[fiu]>nivel[nod])
                muchii_critice.push_back(make_pair(nod,fiu));

        }
        else //muchie de intoarcere
            if(nivel[fiu]<nivel[nod])
                //formula B (curs)
                //actualizez nivelul minim al nodului curent deoarece
                //fiul este deasupra celui curent
                nivel_min[nod]=min(nivel_min[nod], nivel[fiu]);
    }
}

vector<pair<int,int>> GrafNeorientat::MuchiiCritice()
{
    ///complexitate: O(n+m)
    //pentru fiecare nod calculez nivelul si nivelul minim la care poate ajunge
    //actualizez nivelul minim in functie de caz si verific daca gasesc muchie critica

    int nivel[n+1];
    int nivel_min[n+1];
    vector<pair<int,int>> muchii_critice;
    int viz[n+1];

    for(int i=1; i<=n; ++i)
        viz[i]=0;

    nivel[1]=1;

    DFSCritice(1, muchii, viz, muchii_critice, nivel, nivel_min);

    return muchii_critice;
}

vector<int> GrafOrientat::SortareTopologica()
{
    ///complexitate O(n+m)

    int grad_intern[n], vf;
    queue<int> sortare;
    vector<int> sortare_top;

    for(int i=1; i<=n; i++)
        grad_intern[i]=0;
    //calculez gradele interne pentru noduri

    for(int i=1; i<=n; ++i)
        for(int j=0; j<(int)muchii[i].size(); ++j)
            grad_intern[muchii[i][j]]+=1;

    //adaug toate nodurile cu grad intern 0 intr-o coada
    for(int i=1; i<=n; ++i)
        if(grad_intern[i]==0)
            sortare.push(i);


    while(!sortare.empty())
    {
        //extrag primul element din coada
        vf=sortare.front();
        sortare.pop();

        sortare_top.push_back(vf);

        //scad gradele vecinilor (elimin nodul curent in mod fictiv)
        //adaug in coada nodurile cu noul grad intern 0
        for(int j=0; j<(int)muchii[vf].size(); ++j)
        {
            grad_intern[muchii[vf][j]]--;
            if(grad_intern[muchii[vf][j]]==0)
                sortare.push(muchii[vf][j]);
        }
    }

    return sortare_top;
}

vector<vector<int>> GrafOrientat::GrafTranspus()
{
    vector<vector<int>> muchii_transp;
    muchii_transp.resize(n+1);

    for(int i=1; i<=n; ++i)
        for(int j=0; j<(int)muchii[i].size(); ++j)
            muchii_transp[muchii[i][j]].push_back(i);

    return muchii_transp;
}

void GrafOrientat::DFSTimpi(int nod, vector<vector<int>>& muchii_curente, int viz[], stack<int>& timpi_final)
{
    viz[nod]=1;
    for(int i=0; i<(int)muchii_curente[nod].size(); ++i)
        if(!viz[muchii_curente[nod][i]])
            DFSTimpi(muchii_curente[nod][i], muchii_curente, viz, timpi_final);

    timpi_final.push(nod);//am terminat cu un nod, il adaug in stiva
}

void GrafOrientat::DFSCtc (int nod, vector<vector<int>>& muchii_curente, int viz[], int nr_comp, vector<vector<int>>& comp_tare_con)
{
    viz[nod]=1;
    //adaug nodul la componenta tare conexa respectiva
    comp_tare_con[nr_comp].push_back(nod);

    for(int i=0; i<(int)muchii_curente[nod].size(); ++i)
        if(!viz[muchii_curente[nod][i]])
            DFSCtc(muchii_curente[nod][i], muchii_curente, viz, nr_comp, comp_tare_con);

}

pair<vector<vector<int>>, int> GrafOrientat::CTC()
{
    /// alg. lui Kosaraju
    ///complexitate: O(n+m)
    //returneaza componentele tare conexe si numarul lor

    int nr=0; //numarul de ctc
    stack<int> timpi_final;
    vector<vector<int>> muchii_transp; //liste de adiacenta- vector de vectori
    vector<vector<int>> comp_tare_con; //ctc, poz1- elem comp1, poz2-elem comp2 etc.

    comp_tare_con.resize(n+1);

    int viz[n+1];
    for(int i=1; i<=n; ++i)
        viz[i]=0;

    //pas1: creez graful transpus (muchii inversate)
    muchii_transp=GrafTranspus();

    //pas2: dfs in care retin timpii de final(stiva) pe graful initial
    for(int i=1; i<=n; ++i)
        if(!viz[i])
            DFSTimpi(i, muchii, viz, timpi_final);

    for(int i=1; i<=n; ++i)
        viz[i]=0;

    //pas3:dfs in functie de timpii de final pe graful transpus
    while(!timpi_final.empty())
    {
        int vf=timpi_final.top();
        timpi_final.pop();
        if(!viz[vf])
        {
            nr++;
            DFSCtc(vf, muchii_transp, viz, nr, comp_tare_con);
        }
    }

    return make_pair(comp_tare_con,nr);
}

int Graf::CompConexe()
{
    int viz[n+1];

    for(int i=1; i<=n; ++i)
        viz[i]=0;

    int nr_comp=0;
    for(int i=1; i<=n; ++i)
        if(!viz[i])
        {
            DFS(i, muchii, viz);
            nr_comp++;
        }
    return nr_comp;
}

void Graf::DFS(int nod, vector<vector<int>>& muchii_curente, int viz[])
{
    ///complexitate: O(n+m)

    viz[nod]=1;
    for(int i=0; i<(int)muchii_curente[nod].size(); ++i)
        if(!viz[muchii_curente[nod][i]])
            DFS(muchii_curente[nod][i], muchii_curente, viz);
}

vector<int> Graf::BFS(int start)
{
    ///complexitate: O(n+m)

    queue<int> coada;
    vector<int> viz;
    vector<int> dist;

    viz.resize(n+1);
    dist.resize(n+1);

    for(int i=1; i<=n; ++i)
        dist[i]=-1, viz[i]=0;

    viz[start]=1;
    dist[start]=0;
    coada.push(start);

    while(!coada.empty())
    {
        int vf=coada.front();
        coada.pop();
        for(int i=0; i<(int)muchii[vf].size(); ++i)
            if(!viz[muchii[vf][i]])
            {
                viz[muchii[vf][i]]=1;
                dist[muchii[vf][i]]=dist[vf]+1;
                coada.push(muchii[vf][i]);
            }
    }

    return dist;
}

void GrafOrientat::CitireMuchiiCosturiAdiacentaInvers(ifstream& fin)
{
    ///pastrez muchiile invers (vreau sa aflu cine intra intr-un nod nu cine iese din el)
    int n1,n2,c;
    muchii_costuri_adiacenta.resize(n+1);

    for(int i=0; i<m; ++i)
    {
        fin>>n1>>n2>>c;
        muchii_costuri_adiacenta[n2].push_back(make_pair(n1,c));
    }
}

void GrafOrientat::CitireMatricePonderi(ifstream& fin)
{
    int muchii=0,x;
    matrice_ponderi.resize(n);

    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
        {
            fin>>x;
            matrice_ponderi[i].push_back(x);
            if(matrice_ponderi[i][j])
                muchii++;
        }

    m=muchii;
}

void GrafOrientat::CitireMuchiiAdiacenta(ifstream& fin)
{
    int n1,n2;
    muchii.resize(n+1);

    for(int i=0; i<m; ++i)
    {
        fin>>n1>>n2;
        muchii[n1].push_back(n2);
    }
}

void GrafNeorientat::CitireMuchiiAdiacenta(ifstream& fin)
{
    int n1,n2;
    muchii.resize(n+1);

    for(int i=0; i<m; ++i)
    {
        fin>>n1>>n2;
        muchii[n1].push_back(n2);
        muchii[n2].push_back(n1);
    }
}

void GrafNeorientat::CitireMuchiiAdiacentaCuplaj(ifstream& fin, int nr1)
{
    //primeste ca parametru si numarul de noduri din multimea stanga
    //vor fi adaugate muchii orientate pentru simplificarea aflarii cuplajului

    int n1,n2;
    muchii.resize(nr1+1);

    for(int i=0; i<m; ++i)
    {
        fin>>n1>>n2;
        muchii[n1].push_back(n2);
    }
}

void GrafOrientat::CitireMuchiiCosturiAdiacenta(ifstream& fin)
{
    int n1,n2,c;
    muchii_costuri_adiacenta.resize(n+1);

    for(int i=0; i<m; ++i)
    {
        fin>>n1>>n2>>c;
        muchii_costuri_adiacenta[n1].push_back(make_pair(n2,c));
    }
}

void GrafNeorientat::CitireMuchiiCosturiAdiacenta(ifstream& fin)
{
    int n1,n2,c;
    muchii_costuri_adiacenta.resize(n+1);

    for(int i=0; i<m; ++i)
    {
        fin>>n1>>n2>>c;
        muchii_costuri_adiacenta[n1].push_back(make_pair(n2,c));
        muchii_costuri_adiacenta[n2].push_back(make_pair(n1,c));
    }
}

void Graf::CitireMuchiiCosturiLista(ifstream& fin)
{
    int n1,n2,c;

    for(int i=0; i<m; ++i)
    {
        fin>>n1>>n2>>c;
        muchii_costuri_lista.push_back(make_tuple(n1,n2,c));
    }
}

void GrafNeorientat::CitireMuchii(ifstream& fin)
{
    int a,b;

    muchii_adiacenta_indici.resize(n+1);

    for(int i=0; i<m; ++i)
    {
        //muchia cu numarul i
        fin>>a>>b;
        capete_muchii.push_back(make_pair(a,b));
        muchii_adiacenta_indici[a].push_back(i);
        muchii_adiacenta_indici[b].push_back(i);

    }
}

bool Graf::HavelHakimi(vector<int>& grade)
{
    ///complexitate: O(n^2 logn)

    int suma=0;
    int n=grade.size();

    for(int i=0; i<n; ++i)
        suma+=grade[i];

    //daca suma e impara sau unul din grade>n-1 nu e corect
    if(suma%2==1)
        return 0;

    for(int i=0; i<n; ++i)
        if(grade[i]>n-1)
            return 0;

    //sortez descrescator gradele
    sort(grade.begin(),grade.end(), greater<int> ());
    while(grade[0]!=0)
    {
        //decrementez cu 1 gradele de la poz 1 la grade[0]
        for(int i=1; i<=grade[0]; ++i)
        {
            grade[i]--;
            if(grade[i]<0)
                return 0;
        }
        grade[0]=0;
        sort(grade.begin(),grade.end(), greater<int> ());
    }
    return 1;
}

bool GrafOrientat::ExistaLantNesaturat(int s,int d, vector<int>& tata, vector<vector<int>>& muchii,
                                       vector<vector<int>>& c, vector<vector<int>>& f, vector<int>& viz)
{
    queue<int> C;

    //golesc viz
    viz.clear();
    viz.resize(n+1,0);

    C.push(s);
    viz[s]=1;
    tata[d]=0;

    //construiesc drum pana la d
    while(!C.empty() && tata[d]==0)
    {
        int v=C.front();
        C.pop();

        //iau vecinii
        for(int i=0; i<(int)muchii[v].size(); ++i)
        {
            int vecin=muchii[v][i];

            if(!viz[vecin] && c[v][vecin]>f[v][vecin])
            {
                C.push(vecin);
                viz[vecin]=1;
                tata[vecin]=v;
            }
        }
    }

    if(tata[d])
        return true; //exista drum pentru ca destinatia are un tata
    else
        return false;
}

int GrafOrientat::FluxMaxim(int s, int d)
{
    ///algoritmul Ford-Fulkerson

    vector<int> tata;
    vector<vector<int>> muchii; //consider graful neorientat
    vector<vector<int>> c; //matrice capacitate
    vector<vector<int>> f; //matrice flux
    vector<int> viz; //pentru a retine nodurile vizitate in bfs
    int flux_total=0;

    muchii.resize(n+1);
    tata.resize(n+1,0);
    c.resize(n+1,vector<int>(n+1,0));
    f.resize(n+1,vector<int>(n+1,0));

    //initializare (muchii, capacitate)
    for(int i=1; i<=n; ++i)
        for(int j=0; j<(int)muchii_costuri_adiacenta[i].size(); ++j)
        {
            int nod=get<0>(muchii_costuri_adiacenta[i][j]);
            int cost=get<1>(muchii_costuri_adiacenta[i][j]);
            muchii[i].push_back(nod);
            muchii[nod].push_back(i);
            c[i][nod]=cost;
        }

    while(ExistaLantNesaturat(s,d,tata,muchii,c,f,viz))
    {
        //iau toate drumurile care intra in d
        for(int i=0; i<(int)muchii[d].size(); ++i)
        {
            int v=muchii[d][i];
            if(f[v][d]!=c[v][d] && viz[v]) //muchie valida
            {
                tata[d]=v;
                int val_minima=110005;//un numar suficient de mare

                //aflu cu cat pot modifica fluxul pe drumul de la s la d (pornind din d) => capacitatile reziduale
                for (int i=d; i!=s; i=tata[i])
                    if(c[tata[i]][i]-f[tata[i]][i]<val_minima)
                        val_minima=c[tata[i]][i]-f[tata[i]][i];

                if(val_minima!=0)
                    //parcurg din nou drumul pentru a actualiza fluxurile
                {
                    for (int i=d; i!=s; i=tata[i])
                    {
                        f[tata[i]][i]+=val_minima;
                        f[i][tata[i]]-=val_minima;
                    }
                    flux_total+=val_minima; //actualizez fluxul final
                }
            }
        }
    }

    return flux_total;
}

int Graf::Reprezentant(int x, vector<int>& rad)
{
    /// complexitate: O(h arbore), O(1) amortizat

    int r=x; //radacina (reprezentantul)
    int aux;

    while(rad[r]!=r)
        r=rad[r];

    //compresia drumurilor (unesc nodurile direct de radacina)
    while(rad[x]!=x)
    {
        aux=rad[x];
        rad[x]=r;
        x=aux;
    }
    return r;
}

void Graf::Reuneste(int x, int y, vector<int>& h, vector<int>& rad)
{
    ///complexitate: O(1) (dupa aflarea anterioara a reprez. lui x si y)

    int rx=Reprezentant(x, rad);
    int ry=Reprezentant(y, rad);

    if(h[rx]>h[ry]) //unesc arborele mai mic la cel mai mare
    {
        rad[ry]=rx;
    }
    else
    {
        rad[rx]=ry;

        if(h[rx]==h[ry])
            h[ry]++;
    }
}

bool Graf::Verifica(int x, int y, vector<int>& rad)
{
    ///complexitate: O(1) (dupa aflarea anterioara a reprez. lui x si y)

    int rx=Reprezentant(x, rad);
    int ry=Reprezentant(y, rad);

    return (rx==ry);
}

vector<tuple<int,int,int>> GrafNeorientat::APM()
{
    ///Algoritmul lui Kruskal
    ///complexitate: O(Mlog*N + Mlog2M)
    //lucrez cu  muchii_costuri_lista

    int nr_sel=0; //numar muchii selectate
    int cost_total=0;

    vector<int> h; //h[i]=inaltimea maxima (de la compresie) a arborelui care il cuprinde pe i
    vector<int> rad;//tata[i]=radacina arborelui in care se afla i

    vector<tuple<int,int,int>> muchii_apm;

    //initializare
    h.resize(n+1);
    rad.resize(n+1);

    for(int i=1; i<=n; ++i)
        rad[i]=i, h[i]=1;

    //sortez muchiile crescator dupa cost
    sort(muchii_costuri_lista.begin(), muchii_costuri_lista.end(),
         [](const tuple<int, int, int>& c1, const tuple<int, int, int>& c2)
    {
        return get<2>(c1) < get<2>(c2);
    });

    for(int i=0; i<m && nr_sel<n-1; ++i)
    {
        int nod1=get<0>(muchii_costuri_lista[i]);
        int nod2=get<1>(muchii_costuri_lista[i]);
        int cost=get<2>(muchii_costuri_lista[i]);

        if(!Verifica(nod1,nod2,rad)) //nu fac parte din  aceeasi componenta conexa
        {
            nr_sel++;
            Reuneste(nod1, nod2, h, rad);
            cost_total+=cost;
            muchii_apm.push_back(make_tuple(nod1,nod2,cost));
        }
    }

    return muchii_apm;
}

void PbComponenteBiconexe()
{
    ifstream fin("biconex.in");
    ofstream fout("biconex.out");

    int n,m;
    vector<vector<int>> componente;

    fin>>n>>m;
    GrafNeorientat g(n,m);
    g.CitireMuchiiAdiacenta(fin);

    componente=g.ComponenteBiconexe();

    fout<<componente.size()<<"\n";
    for(int i=0; i<(int)componente.size(); ++i)
    {
        for(int j=0; j<(int)componente[i].size(); ++j)
            fout<<componente[i][j]<<" ";
        fout<<"\n";
    }

    fin.close();
    fout.close();
}

void PbBFS()
{
    ifstream fin("bfs.in");
    ofstream fout("bfs.out");

    int n,m,s;
    vector<int> dist;
    fin>>n>>m>>s;

    GrafOrientat g(n,m);
    g.CitireMuchiiAdiacenta(fin);

    dist=g.BFS(s);

    for(int i=1; i<=n; ++i)
        fout<<dist[i]<<" ";

    fout<<"\n";

    fin.close();
    fout.close();
}

void PbBellmanford()
{
    ifstream fin("bellmanford.in");
    ofstream fout("bellmanford.out");

    int n,m;
    vector<int> d;
    fin>>n>>m;

    GrafOrientat g(n,m);
    g.CitireMuchiiCosturiAdiacenta(fin);

    d=g.BellmanFord(1);

    if(d.size())
    {
        for(int i=2; i<=n; ++i)
            fout<<d[i]<<" ";
        fout<<"\n";
    }
    else
        fout<<"Ciclu negativ!\n";

    fin.close();
    fout.close();
}

void PbDisjoint()
{
    ifstream fin("disjoint.in");
    ofstream fout("disjoint.out");

    int n,m;
    int cod,x,y;
    vector<int> h; //h[i]=inaltimea maxima (de la compresie) a arborelui care il cuprinde pe i
    vector<int> rad;//tata[i]=radacina arborelui in care se afla

    fin>>n>>m;

    //initializare
    h.resize(n+1);
    rad.resize(n+1);
    for(int i=1; i<=n; ++i)
        rad[i]=i, h[i]=1;

    GrafOrientat g(n,0);

    for(int i=0; i<m; ++i)
    {
        fin>>cod>>x>>y;
        if(cod==1) //union
            g.Reuneste(x,y,h,rad);

        else //find
        {
            bool ok=g.Verifica(x,y,rad);
            if(ok)
                fout<<"DA\n";
            else
                fout<<"NU\n";

        }
    }

    fin.close();
    fout.close();
}

void PbDijkstra()
{
    ifstream fin("dijkstra.in");
    ofstream fout("dijkstra.out");

    int n,m,s=1;
    const int inf = 250005;
    vector<int> d;

    fin>>n>>m;

    GrafOrientat g(n,m);
    g.CitireMuchiiCosturiAdiacenta(fin);

    d=g.Dijkstra(s);  //distantele de la nodul start la celelalte noduri

    for(int i=1; i<=n; ++i)
        if(i!=s)
        {
            if(d[i]!=inf)
                fout<<d[i]<<" ";
            else fout<<0<<" ";
        }

    fout<<"\n";

    fin.close();
    fout.close();
}

void PbHavelHakimi()
{
    ifstream fin("havelhakimi.in");
    ofstream fout("havelhakimi.out");

    vector<int> grade;
    int grad;

    //citire grade
    while(fin>>grad)
        grade.push_back(grad);

    bool rez=Graf::HavelHakimi(grade);
    if(rez)
        fout<<"DA\n";
    else
        fout<<"NU\n";

    fin.close();
    fout.close();
}

void PbAPM()
{
    ifstream fin("apm.in");
    ofstream fout("apm.out");

    int n,m;
    vector<tuple<int,int,int>> muchii_apm;
    int cost_total=0;

    fin>>n>>m;

    GrafNeorientat g(n,m);
    g.CitireMuchiiCosturiLista(fin);

    muchii_apm=g.APM();

    for(int i=0; i<(int)muchii_apm.size(); ++i)
        cost_total+=get<2>(muchii_apm[i]);

    fout<<cost_total<<"\n";
    fout<<muchii_apm.size()<<"\n";

    for(int i=0; i<(int)muchii_apm.size(); ++i)
        fout<<get<0>(muchii_apm[i])<<" "<<get<1>(muchii_apm[i])<<"\n";

    fin.close();
    fout.close();
}

void PbMuchiiCritice()
{
    ifstream fin("critice.in");
    ofstream fout("critice.out");

    int n,m;
    vector<pair<int,int>> muchii_critice;
    fin>>n>>m;

    GrafNeorientat g(n,m);
    g.CitireMuchiiAdiacenta(fin);

    muchii_critice=g.MuchiiCritice();

    for(int i=0; i<(int)muchii_critice.size(); ++i)
    {
        fout<<muchii_critice[i].first<<" "<<muchii_critice[i].second<<"\n";
    }

    fin.close();
    fout.close();
}

void PbCTC()
{
    ifstream fin("ctc.in");
    ofstream fout("ctc.out");

    int n,m;
    pair<vector<vector<int>>, int> res;
    vector<vector<int>> comp_tare_con;
    int nr_comp;

    fin>>n>>m;
    GrafOrientat g(n,m);
    g.CitireMuchiiAdiacenta(fin);

    res=g.CTC();
    comp_tare_con=get<0>(res);
    nr_comp=get<1>(res);

    fout<<nr_comp<<"\n";
    for(int i=1; i<=nr_comp; ++i)
    {
        for(int j=0; j<(int)comp_tare_con[i].size(); ++j)
            fout<<comp_tare_con[i][j]<<" ";
        fout<<"\n";
    }

    fin.close();
    fout.close();
}

void PbCompConexe()
{
    ifstream fin("dfs.in");
    ofstream fout("dfs.out");

    int n,m;
    fin>>n>>m;

    GrafNeorientat g(n,m);
    g.CitireMuchiiAdiacenta(fin);

    fout<<g.CompConexe();

    fin.close();
    fout.close();
}

void PbSortareTopologica()
{
    ifstream fin("sortaret.in");
    ofstream fout("sortaret.out");

    int n,m;
    vector<int> sortare_top;
    fin>>n>>m;

    GrafOrientat g(n,m);
    g.CitireMuchiiAdiacenta(fin);

    sortare_top=g.SortareTopologica();

    if((int)sortare_top.size()!=n) //graful are cicluri
        fout<<"Nu exista sortare topologica!\n";
    else
        for(int i=0; i<(int)sortare_top.size(); ++i)
            fout<<sortare_top[i]<<" ";

    fout<<"\n";

    fin.close();
    fout.close();
}

void PbDiametruArbore()
{
    ifstream fin("darb.in");
    ofstream fout("darb.out");

    int n;
    fin>>n;
    GrafNeorientat g(n,n-1);
    g.CitireMuchiiAdiacenta(fin);

    fout<<g.DiametruArbore();

    fin.close();
    fout.close();
}

void PbRoyFloyd()
{
    ifstream fin("royfloyd.in");
    ofstream fout("royfloyd.out");

    int n;
    const int valmax=100000000;
    vector<vector<int>> drumuri_minime;
    fin>>n;

    GrafOrientat g(n,0);
    g.CitireMatricePonderi(fin);

    drumuri_minime=g.RoyFloyd();
    for(int i=0; i<n; ++i)
    {
        for(int j=0; j<n; ++j)
            if(drumuri_minime[i][j]!=valmax && i!=j)
                fout<<drumuri_minime[i][j]<<" ";
            else
                fout<<0<<" ";
        fout<<"\n";
    }

    fin.close();
    fout.close();
}

void PbFluxMaxim()
{
    ifstream fin("maxflow.in");
    ofstream fout("maxflow.out");

    int n,m;
    fin>>n>>m;
    GrafOrientat g(n,m);
    g.CitireMuchiiCosturiAdiacenta(fin);

    fout<<g.FluxMaxim(1,n);

    fin.close();
    fout.close();
}

void PbCicluEuler()
{
    ifstream fin("ciclueuler.in");
    ofstream fout("ciclueuler.out");

    int n,m;
    vector<int> sol;
    fin>>n>>m;

    GrafNeorientat g(n,m);
    g.CitireMuchii(fin);

    sol=g.CicluEulerian();

    if(sol.size()==1 && sol[0]==-1)
        fout<<-1; //nu exista ciclu
    else
        for(int i=0; i<(int)sol.size()-1; ++i) //nu afisam si nodul care inchide ciclul
            fout<<sol[i]<<" ";

    fout<<"\n";

    fin.close();
    fout.close();
}

void PbCuplajMaxim()
{
    ifstream fin("cuplaj.in");
    ofstream fout("cuplaj.out");

    int n,m,e;
    vector<pair<int,int>> sol;
    fin>>n>>m>>e;

    GrafNeorientat g(n+m,e);
    g.CitireMuchiiAdiacentaCuplaj(fin, n);

    sol=g.CuplajMaxim(n, m);

    fout<<sol.size()<<"\n";
    for(int i=0; i<(int)sol.size(); ++i)
        fout<<sol[i].first<<" "<<sol[i].second<<"\n";

    fin.close();
    fout.close();
}

void PbHamilton()
{
    ifstream fin("hamilton.in");
    ofstream fout("hamilton.out");

    int n,m;
    fin>>n>>m;

    GrafOrientat g(n,m);
    g.CitireMuchiiCosturiAdiacentaInvers(fin);
    int sol=g.CicluHamiltonian();

    if(sol==0)
        fout<<"Nu exista solutie\n";
    else
        fout<<sol<<"\n";

    fin.close();
    fout.close();
}

int main()
{
    PbCompConexe();

    return 0;
}

