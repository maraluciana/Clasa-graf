#include <iostream>
#include <fstream>
#include <stdio.h>
#include <queue>
#include <vector>
#include <utility>
#include <bits/stdc++.h>

using namespace std;

ifstream fin;
ofstream fout;

#define NMAX 100001

class Graf{
private:
    int nrNod, nrMuch;
    bool orientat;
    vector<vector<int>> listaAd;

public:
    Graf(int nrNoduri = 0, int nrMuchii = 0, bool eOrientat = false)
    {
        this->nrNod = nrNoduri;
        this->nrMuch = nrMuchii;
        this->orientat = eOrientat;
    }

    ~Graf()
    {
        this->nrNod = 0;
        this->nrMuch = 0;
        listaAd.clear();
    }

    friend class Probleme;

    void set_nrNod(int &);
    void set_nrMuch(int &);
    int get_nrNod();
    int get_nrMuch();

//tema 1 - BF-DF-si-aplicatii
    void bfs(int, vector<int>&);
    void dfs(int, vector<int>&, int, stack<int>&,vector<int>&, vector<vector<int>>);
    bool havel_hakimi(vector<int>&, int);
//tema 2 - Drumuri-minime-si-APM
    void init(vector<int>&,vector<int>&);
    int reprez(int,vector<int>&);
    void unite(int,int,vector<int>&,vector<int>&, vector<pair <int, int>>& muchii_apm);
    void apm_kruskall(int&, vector<pair<int, pair<int, int>>>&,vector<int>&,vector<int>&, vector<pair <int, int>>&);
    void dijkstra(int , vector<int>& , list<pair<int, int> > *muchii_dij);
    bool bellman_ford(int, vector<int>&, list<pair<int, int> > *muchii_dij);
//tema 3 - Maxflow-Royfloyd-Darb
    int darb(int,int&);
    void roy_floyd(vector<vector<int>>&);
    bool bfs_flow(int, int, vector<int>&, vector<vector<int>>&, vector<vector<int>>&, vector<int>&);
    int max_flow(vector<int>&, vector<vector<int>>&,vector<vector<int>>&,vector<int>&);
//tema 4 - Ciclu Eulerian
    bool ciclueurian(vector<int>&, int, list<pair<int, int> > *lista);
};

class Probleme{
public:
//tema 1 - BF-DF-si-aplicatii
    void bfs_infoarena();
    void dfs_infoarena();
    void havel_hakimi();
    void sortare_top_infoarena();
    void ctc_infoarena();
//tema 2 - Drumuri-minime-si-APM
    void apm_infoarena();
    void disjoint_infoarena();
    void dijkstra_infoarena();
    void bellman_ford_infoarena();
//tema 3 - Maxflow-Royfloyd-Darb
    void darb_infoarena();
    void roy_floyd_infoarena();
    void maxflow_infoarena();
//tema 4 - Ciclu Eulerian
    void ciclueulerian_infoarena();
};


void Graf::set_nrNod(int &n) {nrNod = n;}

void Graf::set_nrMuch(int &m) {nrMuch = m;}

int Graf::get_nrNod() {return nrNod;}

int Graf::get_nrMuch() {return nrMuch;}

bool Graf::ciclueurian(vector<int>& ciclu, int nod, list<pair<int, int>> *lista)
{
    for(int i = 1; i <= nrNod; ++i)
        if (lista[i].size() % 2 == 1)
            return false;

    stack<int> st;
    int viz_edge[nrMuch+1] = {0};

    st.push(nod);

    while (!st.empty())
    {
        int x = st.top();

        if(!lista[x].empty())
        {

            int e = lista[x].back().second;
            int vecin = lista[x].back().first;

            lista[x].pop_back();

            if (!viz_edge[e])
            {
                viz_edge[e] = 1;
                st.push(vecin);
            }
        }
        else
        {
            st.pop();
            ciclu.push_back(x);
        }
    }

    return true;
}

bool Graf::bellman_ford(int startNod, vector<int>& dist, list<pair<int, int> > *muchii_dij)
{
    int in_queue[nrNod+1]={0};
    int viz[nrNod+1]={0};

    for(int i = 1; i <= nrNod; i++){
        dist[i] = INT_MAX;}

	queue<int> q;

	q.push(startNod);
	in_queue[startNod] = 1;
    dist[startNod] = 0;


	while(!q.empty())
    {
        int x = q.front();
        q.pop();

		in_queue[x] = 0;
		viz[x]++;

		if(viz[x] > nrNod)
            return false;

		for(auto i : muchii_dij[x])
        {
            int y = i.first;
            int cost = i.second;

			if(dist[x] + cost < dist[y])
            {
				dist[y] = dist[x] + cost;

				if(!in_queue[x])
				{
					q.push(y);
                    in_queue[y] = 1;
                }
			}
		}
	}

    return true;
}


void Graf::dijkstra(int startNod, vector<int>& dist, list<pair<int, int> > *muchii_dij)
{
    priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> pq;

    int viz[nrNod+1]={0};

    for(int i = 1; i <= nrNod; i++)
        dist[i]=INT_MAX;

    pq.push(make_pair(0, startNod));
    dist[startNod] = 0;

    while (!pq.empty())
    {
        int x = pq.top().second;
        pq.pop();

        if(viz[x])
            continue;
        else
            viz[x] = 1;

        for (auto i : muchii_dij[x])
        {
            int y = i.first;
            int cost = i.second;

            if (dist[y] > dist[x] + cost)
            {
                dist[y] = dist[x] + cost;
                pq.push(make_pair(dist[y], y));
            }
        }
    }
}

void Graf::init(vector<int>& tata,vector<int>& dim)
{
	for(int i = 1; i <= nrNod; i++)
    {
        tata[i] = i;
        dim[i] = 1;
    }
}

int Graf::reprez(int x, vector<int>& tata)
{
	if(tata[x] == x)
        return x;
	return reprez(tata[x], tata);
}

void Graf::unite(int x,int y,vector<int>& tata,vector<int>& dim, vector<pair <int, int>>& muchii_apm)
{
	int repx = reprez(x, tata), repy = reprez(y, tata);

	if (repx == repy)
        return;
	if (dim[repx] >= dim[repy])
	{
		tata[repy] = repx;
		dim[repx] += dim[repy];
	}
	else
	{
		tata[repx] = repy;
		dim[repy] += dim[repx];
	}


    muchii_apm.push_back(make_pair(x,y));
}


void Graf::apm_kruskall(int& cost, vector<pair<int, pair<int, int>>>& muchii_cost, vector<int>& tata, vector<int>& dim, vector<pair <int, int>>& muchii_apm)
{
    init(tata, dim);

	sort(muchii_cost.begin(), muchii_cost.end());

	for(auto m : muchii_cost)
    {
		if(reprez(m.second.first, tata) != reprez(m.second.second, tata))
		{
            unite(m.second.first, m.second.second, tata, dim, muchii_apm);
            cost += m.first;
        }
    }
}

bool Graf::havel_hakimi(vector<int>& seq, int n)
{
    int sum = 0;

    for(int i = 0; i < n; i++)
        sum += seq[i];

    if(sum % 2 != 0)
    {
        return false;
    }

    for(auto i:seq)
        if(i >= n)
        {
            return false;
        }

    sort(seq.begin(), seq.end(), greater<int>());

    while(seq.front() > 0 && seq.back() >= 0)
    {
        n = seq.front();
        seq.erase(seq.begin());

        for(int i = 0; i < n; i++)
            seq[i] -= 1;

        sort(seq.begin(), seq.end(), greater<int>());
    }

    if(seq.back() < 0)
        return false;

    else if(seq.front() == 0)
        return true;

    return 0;
}

void Graf::dfs(int nod, vector<int>& viz, int c, stack<int>& st, vector<int>& v, vector<vector<int>> lista)
{
    viz[nod] = c;

    for(auto vecin:lista[nod])
        if(!viz[vecin])
        {
            viz[vecin] = c;
            v.push_back(vecin);
            dfs(vecin, viz, c, st, v, lista);
        }

    st.push(nod);
}

void Graf::bfs(int S, vector<int>& dist)
{
    queue<int> st;
    int nod;

    st.push(S);
    dist[S] = 0;

    while(!st.empty())
    {
        nod = st.front();
        st.pop();

        for(auto vecin:listaAd[nod])
        {
            if(dist[vecin] == -1)
            {
                st.push(vecin);
                dist[vecin] = dist[nod] + 1;
            }
        }
    }
}

int Graf::darb(int nod, int& diam)
{
    int viz[NMAX+1], d[NMAX+1], cap;
    queue<int> q;

    for(int i = 1; i <= NMAX; i++)
    {
        d[i] = 0;
        viz[i] = 0;
    }

    q.push(nod);
    d[nod] = 1;
    viz[nod] = 1;

    int x;

    while(!q.empty())
    {
        x = q.front();

        for(auto vecin:listaAd[x])
        {
            if(!viz[vecin])
            {
                q.push(vecin);
                viz[vecin] = 1;

                d[vecin] = d[x] + 1;

                diam = d[vecin];
                cap = vecin;
            }
        }

        q.pop();
    }

    return cap;
}

void Graf::roy_floyd(vector<vector<int>>& A)
{
    for(int k = 1; k <= nrNod; k++)
        for(int i = 1; i <= nrNod; i++)
            for(int j = 1; j <= nrNod; j++)
                if((i!=j) && A[i][k] && A[k][j] && (A[i][j] > A[i][k] + A[k][j] || !A[i][j]))
                    A[i][j] = A[i][k] + A[k][j];
}

bool Graf::bfs_flow(int s, int fin, vector<int>& t, vector<vector<int>>& c, vector<vector<int>>& f, vector<int>& viz)
{
    int x;

    for(int i = 1; i <= nrNod; i++)
        viz[i] = 0;

    queue<int> q;
    q.push(s);
    viz[s] = 1;

    while(!q.empty())
    {
        x = q.front();
        q.pop();

        if (x == fin)
            return true;

        for(auto vecin:listaAd[x])
        {
            if (!viz[vecin] && (c[x][vecin] != f[x][vecin]))
            {
                q.push(vecin);
                t[vecin] = x;
                viz[vecin] = 1;
            }
        }
    }

    return false;

}

int Graf::max_flow(vector<int>& t, vector<vector<int>>& c, vector<vector<int>>& f, vector<int>& viz)
{
    int rasp = 0, p, nod;
    int path_flow;

    while(bfs_flow(1, nrNod, t, c, f,viz))
    {
        for(auto vecin: listaAd[nrNod])
        {
            if((c[vecin][nrNod] != f[vecin][nrNod]) && viz[vecin])
            {
                t[nrNod] = vecin;

                path_flow = INT_MAX;

                for(nod = nrNod; nod != 1; nod = t[nod])
                {
                    p = t[nod];
                    path_flow = min(path_flow, c[p][nod] - f[p][nod]);
                }

                for(nod = nrNod; nod != 1; nod = t[nod])
                {
                    p = t[nod];
                    f[p][nod] += path_flow;
                    f[nod][p] -= path_flow;
                }

                rasp += path_flow;
            }
        }
    }

    return rasp;

}

void Probleme::bfs_infoarena()
{
    fin.open("bfs.in");
    fout.open("bfs.out");

    //N - nr noduri ; M - nr muchii; S - start
    int N, M, S, x, y;

    fin >> N >> M >> S;

    vector<int> dist;

    dist.resize(N+2);

    for(int i = 1; i <= N; i++)
        dist[i] = -1;


    Graf G(N,M,true);

    G.listaAd.resize(N+1);

    for(int i = 0; i < M; i++)
    {
        fin >> x >> y;
        G.listaAd[x].push_back(y);

    }

    G.bfs(S, dist);

    for(int i = 1; i <= N; i++)
        fout << dist[i] << ' ';

}

void Probleme::darb_infoarena()
{
    fin.open("darb.in");
    fout.open("darb.out");

    int n;

    fin >> n;

    Graf G(n,n-1,false);

    int cap, diam;

    int x, y;

    G.listaAd.resize(n+1);

    for(int i = 0; i < (n-1); i++)
    {
        fin >> x >> y;
        G.listaAd[x].push_back(y);
        G.listaAd[y].push_back(x);
    }

    cap = G.darb(1, diam);
    G.darb(cap, diam);

    fout << diam;
}

void Probleme::roy_floyd_infoarena()
{
    fin.open("royfloyd.in");
    fout.open("royfloyd.out");

    int n;
    vector<vector<int>> A;

    fin >> n;

    Graf G(n,0,true);

    int x;
    A.resize(n+1);

    for(int i = 1; i <= n; i++)
    {
        A[i].resize(n+1);
        for(int j = 1; j <= n; j++)
        {
            fin >> x;
            A[i][j] = x;
        }
    }


    G.roy_floyd(A);

    for(int i = 1; i <= n; i++)
    {
        for(int j = 1; j <= n; j++)
        {
            fout << A[i][j] << ' ';
        }

        fout << endl;
    }

}

void Probleme::maxflow_infoarena()
{
    fin.open("maxflow.in");
    fout.open("maxflow.out");

    int n, m, rasp = 0;

    vector<vector<int>> c,f;
    vector<int> viz,t;

    fin >> n >> m;

    t.resize(n+2);
    viz.resize(n+2);

    f.resize(n+2);

    for(int i = 1; i <= n; i++)
    {
        f[i].resize(n+2);
    }

    for(int i = 1; i <= n; i++)
        for(int j = 1; j <= n; j++)
            f[i][j] = 0;

    Graf G (n,m,true);

    int x,y,cap;

    G.listaAd.resize(n+2);
    c.resize(n+2);

    for(int i = 1; i <= n; i++)
    {
        c[i].resize(n+2);
    }

    for(int i = 1; i <= n; i++)
        for(int j = 1; j <= n; j++)
            c[i][j] = 0;

    for(int i = 0; i < m; i++)
    {
        fin >> x >> y >> cap;

        G.listaAd[x].push_back(y);
        G.listaAd[y].push_back(x);

        c[x][y] += cap;
        c[y][x] = 0;
    }


    rasp = G.max_flow(t,c,f,viz);

    fout << rasp;
}

void Probleme::dfs_infoarena()
{
    fin.open("dfs.in");
    fout.open("dfs.out");


    int N, M, c = 0, X, Y;

    fin >> N >> M;

    vector<int> viz;

    viz.resize(N+1);

    for(int i = 1; i <= N; i++)
        viz[i] = 0;

    Graf G (N,M,false);

    G.listaAd.resize(N+1);

    for(int i = 0; i < M; i++)
    {
        fin >> X >> Y;
        G.listaAd[X].push_back(Y);
        G.listaAd[Y].push_back(X);
    }

    stack<int> st;
    vector<vector<int>> ctc; vector<int> v;

    for(int i = 1; i <= N; i++)
        if(!viz[i])
        {
            c++;
            G.dfs(i,viz,c,st,v,G.listaAd);
        }

    fout << c;
}

void Probleme::havel_hakimi()
{
    fin.open("hh.in");
    fout.open("hh.out");

    int n, x;
    bool rasp;
    vector<int> seq;

    fin >> n;

    Graf G (0,0,false);

    for(int i = 0; i < n; i++)
    {
        fin >> x;
        seq.push_back(x);
    }

    rasp = G.havel_hakimi(seq, n);

    if(rasp)
        fout << "Da";
    else
        fout << "Nu";

}

void Probleme::sortare_top_infoarena()
{
    fin.open("sortaret.in");
    fout.open("sortaret.out");

    stack<int> st;
    int N, M;

    fin >> N >> M;

    int x,y;
    vector<int> viz;

    viz.resize(N+1);

    Graf G(N,M,true);

    for(int i = 1; i <= N; i++)
        viz[i] = 0;

    G.listaAd.resize(N+1);

    for(int i = 0; i < M; i++)
    {
        fin >> x >> y;
        G.listaAd[x].push_back(y);
    }

    vector<vector<int>> ctc; vector<int> v;

    for(int i = 1; i <= N; i++)
        if(!viz[i])
            G.dfs(i,viz,1,st,v,G.listaAd);

    while(!st.empty())
    {
        fout << st.top() << ' ';
        st.pop();
    }

}

void Probleme::ctc_infoarena()
{
    fin.open("ctc.in");
    fout.open("ctc.out");


    stack<int> st;
    vector<vector<int>> listaT;
    vector<vector<int>> ctc;
    vector<int> viz;
    vector<int> vizT;


    int N, M, c = 0;

    int x, y;

    fin >> N >> M;

    Graf G(N,M,true);

    G.listaAd.resize(N+1); listaT.resize(N+1); ctc.resize(N+1); viz.resize(N+1); vizT.resize(N+1);

    for(int i = 1; i <= N; i++)
    {
        viz[i] = 0;
        vizT[i] = 0;
    }


    for(int i = 0; i < M; i++)
    {
        fin >> x >> y;
        G.listaAd[x].push_back(y);
        listaT[y].push_back(x);
    }

    vector<int> vtemp;

    for(int i = 1; i <= N; i++)
        if(!viz[i])
            G.dfs(i,viz,1,st,vtemp,G.listaAd);

    stack<int> temp;
    vector<int> v;

    while(!st.empty())
    {
        x = st.top();
        st.pop();

        if(!vizT[x])
        {
            c++;
            v.push_back(x);
            G.dfs(x,vizT,c,temp,v,listaT);
            ctc[c].swap(v);
        }
    }

    fout << c << endl;

    for(int i = 1; i <= c; i++)
        if(!ctc[i].empty())
        {
            for(auto comp:ctc[i])
                fout << comp << ' ';
            fout << endl;
        }
}

void Probleme::apm_infoarena()
{
    fin.open("apm.in");
    fout.open("apm.out");

    int N, M;

    fin >> N >> M;

    Graf G(N, M, true);

    int x, y, c, cost = 0;
    pair<int, pair<int, int>> p;// cost(f) - x(sf) - y(ss)
    vector<pair<int, pair<int, int>>> muchii_cost;
    vector<pair <int, int>> muchii_apm;
    vector<int> tata;
    vector<int> dim;

    tata.resize(N+1); dim.resize(N+1);

    for(int i = 0; i < M; i++)
    {
        fin >> x >> y >> c;
        p.first = c;
        p.second.first = x;
        p.second.second = y;

        muchii_cost.push_back(p);
    }

    G.apm_kruskall(cost, muchii_cost, tata, dim, muchii_apm);

    fout << cost << endl;
    int n = muchii_apm.size();
    fout <<  n << endl;
    for(auto m : muchii_apm)
    {
        fout << m.first << ' ' << m.second << endl;
    }
}

void Probleme::disjoint_infoarena()
{
    fin.open("disjoint.in");
    fout.open("disjoint.out");

    int N, M, x, y, cod;
    vector<pair <int, int>> muchii_apm;
    vector<int> tata;
    vector<int> dim;

    tata.resize(N+1); dim.resize(N+1);

    fin >> N >> M;

    Graf G(N,M,false);

    G.init(tata,dim);

    for(int i = 0; i < M; i++)
    {
        fin >> cod >> x >> y;

        if(cod == 1)
            G.unite(x,y,tata,dim,muchii_apm);
        else
        {
            if(G.reprez(x,tata) == G.reprez(y,tata))
                fout << "DA" << endl;
            else
                fout << "NU" << endl;
        }
    }
}

void Probleme::dijkstra_infoarena()
{
    fin.open("dijkstra.in");
    fout.open("dijkstra.out");

    int N, M;

    fin >> N >> M;

    Graf G(N, M);

    int x, y, c;
    list<pair<int, int> > *muchii_dij;
    muchii_dij = new list<pair<int,int>> [N + 1];
    vector<int> dist; dist.resize(N+1);


    for(int i = 0; i < M; i++)
    {
        fin >> x >> y >> c;
        muchii_dij[x].push_back(make_pair(y,c));
    }

    G.dijkstra(1, dist, muchii_dij);

    for(int i = 2; i <= N; i++)
        if(dist[i] != INT_MAX)
            fout << dist[i] << ' ';
        else
            fout << 0 << ' ';
}

void Probleme::bellman_ford_infoarena()
{
    fin.open("bellmanford.in");
    fout.open("bellmanford.out");


    int N, M;
    bool rasp;

    fin >> N >> M;

    Graf G(N, M);

    int x, y, c;
    list<pair<int, int> > *muchii_dij;
    muchii_dij = new list<pair<int,int>> [N + 1];
    vector<int> dist; dist.resize(N+1);


    for(int i = 0; i < M; i++)
    {
        fin >> x >> y >> c;
        muchii_dij[x].push_back(make_pair(y,c));
    }

    rasp = G.bellman_ford(1,dist,muchii_dij);

    if(rasp)
    {
        for(int i = 2; i <= N; i++)
        {
            if(dist[i] != INT_MAX)
                fout << dist[i] << ' ';
            else
                fout << 0 << ' ';
        }
    }
    else
        fout << "Ciclu negativ!";
}

void Probleme::ciclueulerian_infoarena()
{
    fin.open("ciclueuler.in");
    fout.open("ciclueuler.out");

    vector<int> ciclu;
    bool rasp;
    int N, M, x, y;

    fin >> N >> M;

    Graf G(N, M);

    list<pair<int, int> > *lista;
    lista = new list<pair<int,int>> [N + 1];

    for(int i = 0; i < M; i++)
    {
        fin >> x >> y;

        lista[x].push_back(make_pair(y, i));
        lista[y].push_back(make_pair(x, i));
    }

    rasp = G.ciclueurian(ciclu, 1, lista);

    if(!rasp)
        fout << "-1";
    else
    {
        for (auto i = ciclu.begin(); i != ciclu.end(); i++)
            fout << *i << " ";
    }
}

int main()
{
    Probleme p;
    int cod;

    cout << "Tasteaza numarul care corespunde problemei dorite" << endl << endl;
    cout << "1.Diametrul unui arbore" << endl;
    cout << "2.Floyd-Warshall/Roy-Floyd" << endl;
    cout << "3.Flux maxim" << endl;
    cout << "4.BFS - Parcurgere in latime" << endl;
    cout << "5.Parcurgere DFS - componente conexe" << endl;
    cout << "6.Havel-Hakimi" << endl;
    cout << "7.Sortare topologica" << endl;
    cout << "8.Arbore partial de cost minim" << endl;
    cout << "9.Paduri de multimi disjuncte" << endl;
    cout << "10.Algoritmul lui Dijkstra" << endl;
    cout << "11.Algoritmul Bellman-Ford" << endl;
    cout << "12.Componente tare conexe" << endl;
    cout << "13.Ciclu Eulerian" << endl << endl;

    cin >> cod;

    switch(cod)
    {
        case 1:
        {
            p.darb_infoarena();
            break;
        }
        case 2:
        {
            p.roy_floyd_infoarena();
            break;
        }
        case 3:
        {
            p.maxflow_infoarena();
            break;
        }
        case 4:
        {
            p.bfs_infoarena();
            break;
        }
        case 5:
        {
            p.dfs_infoarena();
            break;
        }
        case 6:
        {
            p.havel_hakimi();
            break;
        }
        case 7:
        {
            p.sortare_top_infoarena();
            break;
        }
        case 8:
        {
            p.apm_infoarena();
            break;
        }
        case 9:
        {
            p.disjoint_infoarena();
            break;
        }
        case 10:
        {
            p.dijkstra_infoarena();
            break;
        }
        case 11:
        {
            p.bellman_ford_infoarena();
            break;
        }
        case 12:
        {
            p.ctc_infoarena();
            break;
        }
        case 13:
        {
            p.ciclueulerian_infoarena();
            break;
        }
        default:
            cout << "Trebuie un numar din lista :(" << endl;
    }

    p.ciclueulerian_infoarena();

    fin.close();
    fout.close();

    return 0;
}
