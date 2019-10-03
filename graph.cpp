#include "graph.h"
#include <algorithm>

graph::graph()
{
    //ctor
}

graph::~graph()
{
    //dtor
    adjlist.clear();
    adjlist_rev.clear();
}

void graph::read_graph(const char * file_name, bool directed)
{
    //freopen("out.txt","w",stdout);
    is_directed = directed;
    FILE * pFile;
    pFile = fopen(file_name, "rb");
    fread(&n_vertices, sizeof(int32_t), 1, pFile);
    fread(&n_edges, sizeof(int64_t), 1, pFile);
    deg = vector<int32_t>(n_vertices + 1);
    fread(&deg[1], sizeof(int32_t), n_vertices, pFile);
    adjlist = vector < vector <int32_t> >(n_vertices+1,vector <int32_t>(0));
    adjlist_rev = vector < vector <int32_t> >(n_vertices+1,vector <int32_t>(0));
    for (int32_t i = 1; i <= n_vertices; ++i)
    {
        adjlist[i] = vector <int32_t> (deg[i]);
        fread(&adjlist[i][0], sizeof(int32_t), deg[i], pFile);
    }
    if (is_directed)
    {
        for (int32_t i = 1; i <= n_vertices; ++i)
        {
            for (int32_t j = 0; j < adjlist[i].size(); ++j)
                adjlist_rev[adjlist[i][j]].push_back(i);
        }
    }
	else
	{
		for (int32_t i = 1; i <= n_vertices; ++i)
        {
            for (int32_t j = 0; j < adjlist[i].size(); ++j)
                adjlist_rev[adjlist[i][j]].push_back(i);
        }
		for (int32_t i = 1; i <= n_vertices; ++i)
			for (int32_t j = 0; j < adjlist_rev[i].size(); ++j)
                adjlist[i].push_back(adjlist_rev[i][j]);
	}
	for (int32_t i = 1; i <= n_vertices; ++i)
	{
		sort(adjlist[i].begin(),adjlist[i].end());
		adjlist[i].erase(unique(adjlist[i].begin(), adjlist[i].end()), adjlist[i].end());
	}
    for (int32_t i = 1; i <= n_vertices; ++i)
        deg[i] = adjlist[i].size();
	// cerr << "Loaded -- number of nodes: " << n_vertices << endl;
	// cerr << "       -- number of edges: " << n_edges << endl;
}

vector <int32_t> graph::get_adjlist(int32_t u)
{
    return adjlist[u];
}

vector <int32_t> graph::get_adjlist_rev(int32_t u)
{
    if (is_directed) return adjlist_rev[u];
    return adjlist[u];
}

int32_t graph::get_deg(int32_t u)
{
    return adjlist[u].size();
}

int32_t graph::get_deg_rev(int32_t u)
{
    if (is_directed) return adjlist_rev[u].size();
    return adjlist[u].size();
}

int32_t graph::get_n_vertices()
{
    return n_vertices;
}

int64_t graph::get_probability(int32_t u)
{
    return probability[u];
}

void graph::calculate_probability()
{
    probability = vector <int64_t> (n_vertices+1,0);
    for (int32_t i = 1; i <= n_vertices; ++i)
    {
        probability[i] = (int64_t)deg[i] * (deg[i]-1);
        probability[i] += probability[i-1];
    }
// pair<int,int> p;
	// for (int i = 1; i <= n_vertices; ++i)
		// for (int j = 0; j < adjlist[i].size(); ++j)
			// if (i < adjlist[i][j])
			// {
				// p = make_pair(i,adjlist[i][j]);
				// edges_set.insert(p);
			// }
}

bool graph::is_edge(int u, int v)
{
	if (u > v) swap (u,v);
	if (edges_set.count(make_pair(u,v)) > 0) return true;
	return false;
}

