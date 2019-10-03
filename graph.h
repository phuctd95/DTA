#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>

using namespace std;

class graph
{
    private:
        bool is_directed;
        int32_t n_vertices;
        int64_t n_edges;
        vector < vector <int32_t> > adjlist;
        vector < vector <int32_t> > adjlist_rev;
        vector <int32_t> deg;
        vector <int64_t> probability;
		set <pair<int,int> > edges_set;

    public:
        graph();
        virtual ~graph();
        //read graph from file
        void read_graph(const char * file_name, bool directed);
        vector <int32_t> get_adjlist(int32_t u);
        vector <int32_t> get_adjlist_rev(int32_t u);
        int32_t get_deg(int32_t u);
        int32_t get_deg_rev(int32_t u);
        int32_t get_n_vertices();
        int64_t get_probability(int32_t u);
        void calculate_probability();
		bool is_edge(int u, int v);
};

#endif // GRAPH_H
