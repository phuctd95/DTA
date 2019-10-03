#ifndef SAMPLER_H
#define SAMPLER_H
#include <queue>
#include "graph.h"
#include "xorshift.h"

class Sampler
{
    private:
        xorshift rand_gen;
        vector <int32_t> path;
        vector <int32_t> shortest_path_from_s;
        vector <double> n_paths_from_s;
        vector <int32_t> shortest_path_from_t;
        vector <double> n_paths_from_t;
        queue <int32_t> queue_s;
        queue <int32_t> queue_t;
        int64_t sum_deg_s, sum_deg_t;
        vector <int32_t> neigh;
        bool stop;
		int32_t k_path;
        vector <int32_t> layer;
		vector <bool> explored;
        graph &g;
        void get_random_pair(int32_t &s, int32_t &t);
        int32_t get_random_vertex();
        void init();
        void add_s(int32_t u);
        void add_t(int32_t u);
        int32_t get_front_s();
        int32_t get_front_t();
        void add_from_s();
        void add_from_t();
        void find_shortest_path(int32_t &s, int32_t &t);
        void get_random_shortest_path();
        void get_coverage_centrality_sample();
        bool get_random_triangle();

    public:
        Sampler(graph &gr);
        virtual ~Sampler();
        int32_t s,t;
		void init_k_path(int32_t init_k_path);
        void init_seed(int32_t p_init, vector <uint64_t> &s_init);
        // sample a random pair of vertices called s and t(choosen uniformly from set of vertices)
        // repeat this step until s and t are connected the shortest path from s to t is bigger than 1
        // return a path, which is choosen uniformly for set of shortest path from s to t
        // n_fails is the number of fail sampling
        vector <int32_t> bc_polling(int32_t &n_fails, int &zcover, vector<bool> & seedMark, vector<bool> & cSeedMark);		
        // sample s and t same as get_random_shortest_path function, but return the coverage centrality instead
        vector <int32_t> cc_polling(int32_t &n_fails, int &zcover, vector<bool> & seedMark, vector<bool> & cSeedMark);

        // sample a set of 3 vertices (use the method from "Triadic Measures on Graphs: The Power of Wedge Sampling")
        // repeat this step until they are a triangle
        // return a vector contain 3 elements is 3 vertices above
        // n_fails is the number of fail sampling
        // note: you need to run graph::calculate_probability before run this function
        vector <int32_t> triangle_polling(int32_t &n_fails, int &zcover, vector<bool> & seedMark, vector<bool> & cSeedMark);

        // sample a triangle , return 3 edges of this triangle
        vector < pair<int32_t, int32_t> > triangle_edge_polling(int32_t &n_fails, int &zcover, vector<bool> & seedMark, vector<bool> & cSeedMark);
        // sample a path have length less than or equal to k
        // use the method from "K-Path Centrality: A New Centrality Measure in Social Networks"
        // n_fails always equal to 0
        vector <int32_t> k_path_polling(int32_t &n_fails, int &zcover, vector<bool> & seedMark, vector<bool> & cSeedMark);
        vector <int32_t> ds2_polling(int32_t &n_fails, int &zcover, vector<bool> & seedMark, vector<bool> & cSeedMark);
};

#endif // SAMPLER_H
