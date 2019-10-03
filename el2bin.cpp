/*
* Functionality: convert from a graph file in weighted edge list to a binary file
* Syntax:
	./el2bin <graph file input> <binary graph output>

* The graph file input must follow the following format:
	<number of nodes> <number of edges>
	<first node of edge 1> <second node of edge 1> <weight of edge 1>
	...
	<first node of the last edge> <second node of the last edge> <weight of the last edge>

* The binary graph output will be used in our SSA/DSSA algorithms for fast reading

* Author: Hung T. Nguyen (hungnt@vcu.edu)
*/

#include <cstdio>
#include <fstream>
#include <random>
#include "GLib.hpp"
#include <cmath>
#include <cstring>
#include <cstring>
#include "xorshift.h"
#include <iostream>
using namespace std;

xorshift gen_rand;
void to_int(string s, int& u, int& v)
{
	int w;
	// stringstream sss(s);
	// sss >> u >> v;
	sscanf(s.c_str(), "%d %d %d", &u, &v, &w); 
}

int main(int argc, char ** argv)
{
	vector <uint64_t> s_init;
	vector <int> mmap(100000000,-1);
	for (int i = 0; i < 16; ++i)
		s_init.push_back(rand());
	gen_rand.init_seed(rand()%16,s_init);
	ifstream in(argv[1]);
	srand(time(NULL));
	int n = 0,u,v,x;
	long long m = 0;
	float w;

	vector<int> degree(100000000+1,0);
	vector<vector<int> > eList(100000000+1);
	vector<vector<float> > weight(100000000+1);
	vector<float> weightR(100000000+1,0);
	string s;
	while (getline(in,s))
	{
		if (s[0] != 'a') continue;
		s.erase(0,1);	
		to_int(s,u,v);	
		if (mmap[u] == -1)
		{
			mmap[u] = n++;			
		}
		if (mmap[v] == -1)
		{
			mmap[v] = n++;
		}
		u = mmap[u];
		v = mmap[v];
		x = gen_rand.get_max(3);		
		if (x == 0) w = 0.1;
		if (x == 1) w = 0.01;
		if (x == 2) w = 0.001;
		++m;
		// cout << u << " " << v << " " << w << endl;
		degree[v]++;
		eList[v].push_back(u);
		weight[v].push_back(w);
		weightR[u] += 1;
	}
	
	in.close();

		vector<size_t> idx(n);
	printf("%d %lld\n", n, m);
	FILE * pFile;
	pFile = fopen(argv[2],"wb");
	fwrite(&n, sizeof(int), 1, pFile);
	fwrite(&m, sizeof(long long), 1, pFile);

		for (int i = 0; i < n; ++i){
		idx[i] = i;
	}
	vector<int> inv_idx(n);
	for (int i = 0; i < n; ++i){
		inv_idx[idx[i]]	= i;
	}
	
	vector<int> iTmp(n);
	
	for (int i = 0; i < n; ++i){
		iTmp[i] = degree[idx[i]+1];
	}
	
	// Write node degrees
	fwrite(&iTmp[0], sizeof(int), n, pFile);
	
	for (int i = 1; i <= n; ++i){
		// Write neighbors
		for (unsigned int j = 0; j < eList[idx[i-1]+1].size(); ++j){
			iTmp[j] = inv_idx[eList[idx[i-1]+1][j]-1]+1;
		}
		fwrite(&iTmp[0], sizeof(int), eList[idx[i-1]+1].size(), pFile);
	}

	for (int i = 1; i <= n; ++i){
		// Write weights
				fwrite(&weight[idx[i-1] + 1][0], sizeof(float), weight[idx[i-1]+1].size(), pFile);
		}

	fclose(pFile);	
	return 1;
}
