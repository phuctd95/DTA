#ifndef _RWGRAPH_H
#define _RWGRAPH_H
#include <vector>
#include <random>
#include "sfmt/SFMT.h"
#include "mappedHeap.hpp"
#include "HeapData.hpp"
#include "StepwiseHeap.h"
#include <ctime>
#include <cmath>
#include <functional>

typedef uint32_t UI;
// typedef uint64_t ULL;

#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_set_num_threads(int t) { return 1;}
inline omp_int_t omp_get_thread_num() { return 0;}
#endif

// Graph class defines fundamental operators on graph
class Graph
{
	friend class IM_ICSampler;
	friend class IM_LTSampler;
	private:
		UI UI_MAX = 4294967295U;
//		ULL ULL_MAX = 18446744073709551615ULL;

		// number of nodes
		unsigned int numNodes;
		// number of edges
		unsigned int numEdges;
		// adjacency list
		std::vector<std::vector<int> > adjList;
		std::vector<int> node_deg;
		std::vector<std::vector<UI> > weights;
	
	public:
		Graph();
		// get a vector of neighbours of node u
		const std::vector<int> & operator [] (int u) const;
		// return weights of neighbours of node u
		const std::vector<UI> & getWeight(int u) const;

		// get a vector of neighbours of node u
		const std::vector<int> & operator [] (int u);
		// return weights of neighbours of node u
		const std::vector<UI> & getWeight(int u);

		// get degree of node u
		int getDegree(int u) const;
		// get size of the graph
		int getSize() const;
		// get number of edges
		int getEdge() const;
		// read graph from a file
		void readGraphLT(const char * filename);
		// read graph from a file
		void readGraphIC(const char * filename);
		// write the graph to file
		void writeToFile(const char * filename);
};

class HyperGraph
{
	private:
		// store the edges that a node is incident to
		std::vector<std::vector<int> > node_edge;
		// store hyperedges
		std::vector<std::vector<int> > edge_node;
		// complete copy of all generated hyperedge
		std::vector<std::vector<int> > all_edge_node;
		std::vector<std::vector<int> > all_node_edge;

		unsigned int numNodes;
		std::vector<bool> markRemovedEdges;
		int deletedSize, totalSize, fullSize, maxSize;
		bool do_check;

	public:
		Stepwise_heap sheap;
		HyperGraph();
		void reinitialize(unsigned int k, unsigned int b);
		HyperGraph(unsigned int n, bool c, unsigned int k, unsigned int b);
		void addHyperedge(std::vector<int> &e);
		void addHyperedge_all(std::vector<int> &e);
		void updateDeg();
		void updateEdge();
		void addEdge(std::vector<int> & edge);
		void addEdgeD(std::vector<int> & edge);
		const std::vector<int> & getEdge(int e) const;
		const std::vector<int> & getEdge(int e);
		const std::vector<int> & getNode(int n) const;
		const std::vector<int> & getNode(int n);
		void clearEdges();
		int getSketchSize();
		double increaseT(std::vector<int> & seeds, int k, std::vector<double> &ve, double alpha);
		double decreaseT(std::vector<int> & seeds, int k, std::vector<double> &ve, double alpha);
		double dualBound(std::vector<int> & seeds, int k);

		double ILPBound(std::vector<int> & seeds, int k);

		void selectNode(int v);
		int getMaxDegree();
		int getMaxDegreeNode();
		int getSumK();		
		void cleanRemovedEdges();
		void checkErrors(){
			int tsize=0,dsize=0;
			std::vector<int> ldegree(numNodes+1,0);
			for (int i = 1; i < numNodes+1; ++i){
				tsize += node_edge[i].size();
				ldegree[i] = node_edge[i].size();
                		for (unsigned int j = 0; j < node_edge[i].size(); ++j){
		                        if (markRemovedEdges[node_edge[i][j]]){
						dsize++;
						ldegree[i]--;
                		        }
		                }
        		}
			if (dsize != deletedSize){
				cout << "Error deleted size: " << dsize << " " << deletedSize << endl;
				exit(0);
			}
			if (tsize != totalSize){
                                cout << "Error total edge size: " << tsize << " " << totalSize << endl;
				exit(0);
                        }
			if (!sheap.checkDegree(ldegree)){
				cout << "Error matching degree " << endl;
				exit(0);
			}
		}
		

//                bool calculateInfluence(std::function< std::vector<int32_t>(int32_t &,int &, std::vector<bool> &) > sampler, std::vector<int> & seeds, int t, double deg, float epsilon, float delta, long long int numSuccessSamples, long long int numTotalSamples, int iter, double normalized);

//                bool MSz(std::function< std::vector<int32_t>(int32_t &,int &, std::vector<bool> &, std::vector<bool>&) > sampler, std::vector<int> & seeds, int t, double & deg, float epsilon, float delta, long long int & numSuccessSamples, long long int & numTotalSamples, double normalized, double c, double thresholdz, std::vector<bool> &cseedmark);

//                void DTA(std::function< std::vector<int32_t>(int32_t &,int &, std::vector<bool> &, std::vector<bool> &) > sampler, std::vector<int> & seeds, int t, float epsilon, float delta, double normalized);
//                void SRA(std::function< std::vector<int32_t>(int32_t &,int &, std::vector<bool> &, std::vector<bool> &) > sampler, std::vector<int> & seeds, int t, float epsilon, float delta, double normalized);

//                long long addHyperedge(std::function< std::vector<int32_t>(int32_t &,int &, std::vector<bool> &) > sampler, int t, long long num);
		void buildSeedSet(std::vector<int> & seeds, int k, std::vector<double> &degree);

//		void runMaxCover_SRA(std::function< std::vector<int32_t>(int32_t &,int &, std::vector<bool> &) > sampler, double epsilon, double delta, int k, int t, double normalized, int zscale);
//		bool runMaxCover_SSAz(std::function< std::vector<int32_t>(int32_t &,int &, std::vector<bool> &) > sampler, double epsilon, double delta, int k, int t, double normalized, int zscale);


//		void generateSamplesToBuffer(std::function< std::vector<int32_t>(int32_t &,int &, std::vector<bool> &) > sampler, int id);
};

#endif
