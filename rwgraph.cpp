#include "rwgraph.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <queue>
#include <cstdlib>
#include <unistd.h>
#include <sstream>
#include <ilcplex/ilocplex.h>

ILOSTLBEGIN

using namespace std;

Graph::Graph()
{
}

const vector<int> & Graph::operator [] (int u) const
{
	return adjList[u];
}


const vector<int> & Graph::operator [] (int u)
{
	return adjList[u];
}


const vector<UI> & Graph::getWeight (int u) const
{
        return weights[u];
}

const vector<UI> & Graph::getWeight (int u)
{
        return weights[u];
}

/*
* get degree of node u
*/
int Graph::getDegree(int u) const
{
	return adjList[u].size();
}

/*
* get the number of nodes
*/
int Graph::getSize() const
{
	return numNodes;
}

/*
* get the number of edges
*/
int Graph::getEdge() const
{
	return numEdges;
}

/*
* read binary graph input for LT model
* difference between LT and IC is for LT we accumulate the weights for fast choosing a random node
*/
void Graph::readGraphLT(const char* filename)
{
	FILE * pFile;
        pFile = fopen(filename, "rb");
        fread(&numNodes, sizeof(int), 1, pFile);
        fread(&numEdges, sizeof(long long), 1, pFile);
        node_deg=vector<int>(numNodes + 1);
        fread(&node_deg[1], sizeof(int), numNodes, pFile);

        vector<int> a;
        vector<UI> b;
        adjList.push_back(a);
        weights.push_back(b);

        for (unsigned int i = 1; i <= numNodes; ++i){
                vector<int> tmp(node_deg[i]);
                fread(&tmp[0], sizeof(int), node_deg[i], pFile);

                adjList.push_back(tmp);
        }

	for (unsigned int i = 1; i <= numNodes; ++i){
		vector<float> tmp(node_deg[i] + 1, 0);
                vector<UI> tmp1(node_deg[i] + 1, 0);
                fread(&tmp[1], sizeof(float), node_deg[i], pFile);

                for(int j = 1;j < node_deg[i] + 1; ++j){
                        tmp[j] += tmp[j-1];
                        if (tmp[j] >= 1){
                                tmp1[j] = UI_MAX;
                        } else {
                                tmp1[j] = tmp[j]*UI_MAX;
                        }
                }

                weights.push_back(tmp1);
                node_deg[i]++;
        }
}

/*
* read input graph for IC model
*/
void Graph::readGraphIC(const char* filename)
{
    	FILE * pFile;
    	pFile = fopen(filename, "rb");
    	fread(&numNodes, sizeof(int), 1, pFile);
    	fread(&numEdges, sizeof(long long), 1, pFile);
    	node_deg=vector<int>(numNodes + 1);
    	fread(&node_deg[1], sizeof(int), numNodes, pFile);

	vector<int> a;
	vector<UI> b;
    	adjList.push_back(a);
    	weights.push_back(b);

        for (unsigned int i = 1; i <= numNodes; ++i){
                vector<int> tmp(node_deg[i]);
                fread(&tmp[0], sizeof(int), node_deg[i], pFile);
                adjList.push_back(tmp);
        }

        for (unsigned int i = 1; i <= numNodes; ++i){
                vector<float> tmp(node_deg[i] + 1, 0);
		vector<UI> tmp1(node_deg[i] + 1, 0);
                fread(&tmp[1], sizeof(float), node_deg[i], pFile);

                for(int j = 1;j < node_deg[i] + 1; ++j){
                        tmp1[j] = tmp[j]*UI_MAX;
                }

		if (tmp1[node_deg[i]] <= 0)
                        tmp1[node_deg[i]] = UI_MAX;

                weights.push_back(tmp1);
        }
}

void Graph::writeToFile(const char * filename)
{/*
	ofstream output(filename);
	for (unsigned int i = 0; i < numNodes; ++i){
		for (unsigned int j = 0; j < adjList[i].size(); ++j){
			if (adjList[i][j] > i){
				output << adjList[i][j] << " " << i << " " << weights[i][j] << endl;
			}
		}
	}
	output.close();
*/
}

HyperGraph::HyperGraph()
{
	cout << "An empty hypergraph is created!" << endl;
	sheap = Stepwise_heap(1,1,1);
}

HyperGraph::HyperGraph(unsigned int n, bool c, unsigned int k, unsigned int b)
{
	node_edge = vector<vector<int> >(n+1);
	numNodes = n;

	do_check = c;
	if (do_check) all_node_edge = vector<vector<int> >(n+1);

	sheap = Stepwise_heap(n,k,b);
	deletedSize = 0;
	totalSize = 0;
	maxSize = 0;
}

void HyperGraph::reinitialize(unsigned int k, unsigned int b)
{
	node_edge = vector<vector<int> >(numNodes+1);
	edge_node = vector<vector<int> >();
	markRemovedEdges = vector<bool>();
	deletedSize = 0;
	totalSize = 0;
	//maxSize = 0;
	all_edge_node = vector<vector<int> >();
        if (do_check) all_node_edge = vector<vector<int> >(numNodes+1);

	sheap = Stepwise_heap(numNodes,k,b);
}

void HyperGraph::addHyperedge(vector<int> &e)
{
	unsigned int edgesize = edge_node.size();
	unsigned int num2 = e.size();
	totalSize += num2;
	maxSize = max(maxSize,totalSize);	

	for (unsigned int j = 0; j < num2; ++j){
		sheap.increase(e[j]);
		node_edge[e[j]].push_back(edgesize);
	}
        edge_node.push_back(e);
	markRemovedEdges.push_back(false);

        if (do_check){
                for (unsigned int j = 0; j < num2; ++j){
                        all_node_edge[e[j]].push_back(all_edge_node.size());
                }
                all_edge_node.push_back(e);
        }
	//checkErrors();
}

void HyperGraph::addHyperedge_all(vector<int> &e)
{
	unsigned int num2 = e.size();
	for (unsigned int j = 0; j < num2; ++j){
                all_node_edge[e[j]].push_back(all_edge_node.size());
        }
        all_edge_node.push_back(e);
}

void HyperGraph::selectNode(int v)
{
        for (int i:node_edge[v]){
		if (!markRemovedEdges[i]){
			deletedSize += edge_node[i].size();
			markRemovedEdges[i] = true;
			for (int j:edge_node[i]){
				sheap.decrease(j);
			}
		}
	}
	if (deletedSize > (totalSize+numNodes)/3+1){
		cout << "Clean removed edges" << endl;
		cleanRemovedEdges();
	}
	sheap.decrease(v);
	//checkErrors();
}

void HyperGraph::cleanRemovedEdges()
{
	//checkErrors();
        for (int i = 1; i < numNodes+1; ++i){
                for (unsigned int j = 0; j < node_edge[i].size(); ++j){
                        if (markRemovedEdges[node_edge[i][j]]){
                                node_edge[i][j] = node_edge[i].back();
                                node_edge[i].pop_back();
                                j--;
                        }
                }
        }
	totalSize -= deletedSize;
	deletedSize = 0;
	//checkErrors();
}

int HyperGraph::getSketchSize()
{
	return maxSize;
}


int HyperGraph::getMaxDegree()
{
	return sheap.getMaxDegree();
}

int HyperGraph::getMaxDegreeNode()
{
	return sheap.getTop();
}

int HyperGraph::getSumK()
{
	return sheap.getSumK();
}

// not completed or used
void HyperGraph::updateDeg(){
	/*
	unsigned int num=edge_node.size();
	unsigned int num2,num3,num4,etmp;
	for (unsigned int i = curEdge; i < num; ++i){
		num2 = edge_node[i].size();
		for (unsigned int j = 0; j < num2; ++j){
			if (j >= edge_node[i].size()){
				cout << j << " " << edge_node[i].size() << " " << edge_node[i][j] << endl << endl;
			}
			nodedegree[edge_node[i][j]]++;
			if (nodedegree[edge_node[i][j]] >= z){
				zcover += nodedegree[edge_node[i][j]];
				num3 = edge_node[i][j];
				for (unsigned int k = 0; k <= j; ++k){
					nodedegree[edge_node[i][k]]--;
					node_edge[edge_node[i][k]].pop_back();
				}
				for (int k = 0; k < nodedegree[num3]; ++k){
					etmp = node_edge[num3][k];
					num4 = edge_node[etmp].size();
					for (unsigned int l = 0; l < num4; ++l){
						nodedegree[edge_node[etmp][l]]--;
						node_edge[edge_node[etmp][l]].pop_back();
					}
					edge_node[etmp] = edge_node[num-1];
					num--;
					i--;
					edge_node.pop_back();
				}

				selectedSeed.push_back(num3);
				seedMark[num3] = true;
			} else {
				node_edge[edge_node[i][j]].push_back(i);
			}
		}
	}
	curEdge = edge_node.size();
	*/
}

void HyperGraph::updateEdge(){
	//curEdge = edge_node.size();
}

/*
* Add a hyperedge into the hypergraph
*/
void HyperGraph::addEdge(vector<int> & edge)
{
	edge_node.push_back(edge);
	unsigned int ind = edge_node.size() - 1;
	for (unsigned int i = 0; i < edge.size(); ++i)
		node_edge[edge[i]].push_back(ind);
}

/*
* Add a hyperedge into the hypergraph while keeping track of the node with max degree
*/
void HyperGraph::addEdgeD(vector<int> & edge)
{
	/*
        edge_node.push_back(edge);
        int ind = edge_node.size() - 1;
        for (unsigned int i = 0; i < edge.size(); ++i){
                node_edge[edge[i]].push_back(ind);
		if (node_edge[edge[i]].size() > maxDegree)
			maxDegree = node_edge[edge[i]].size();
	}
	*/
}

/*
* get an edge from the hypergraph
*/
const vector<int> & HyperGraph::getEdge(int e) const{
	return edge_node[e];
}

const vector<int> & HyperGraph::getEdge(int e){
	return edge_node[e];
}

/*
* get the list of hyperedges incident to node n
*/
const vector<int> & HyperGraph::getNode(int n) const{
	return node_edge[n];
}

const vector<int> & HyperGraph::getNode(int n){
	return node_edge[n];
}

/*
* remove all the hyperedges
*/
void HyperGraph::clearEdges()
{
	edge_node.clear();
	node_edge.clear();
	cout << "clear edges!" << endl;
	//maxDegree = 0;
}

/*
bool HyperGraph::calculateInfluence(std::function< std::vector<int32_t>(int32_t &,int &, std::vector<bool> &) > sampler, vector<int> & seeds, int t, double deg, float epsilon, float delta, long long int numSuccessSamples, long long int numTotalSamples, int iter, double normalized)
{
	double dualB = 0;
	if (do_check){
                  dualB = dualBound(numSeeds);
        }
	
	long long countSuccess = 0;
	long long int countTotal = 0;
	unsigned k = seeds.size();
	vector<unsigned int> link(numNodes+1, k);
	double f = (log(7/delta)+lgamma(numNodes+1)-lgamma(k+1)-lgamma(numNodes-k+1))*numNodes/(k*log(7*log2(numNodes)/delta));
	double lambda1 = 1+(1+2.0/9)*(2+1.0/3)*log(3*log2(f)/delta)*9;
	cout << "Seed size: " << k << " " << seeds.size() << " " << numSuccessSamples << endl;
	double degree=0;
	for (unsigned int i = 0; i < seeds.size();++i){
		link[seeds[i]] = i;
	}
	vector<bool> maxSeed(t, false);

	int32_t n_fail = 0;

	omp_set_num_threads(t);
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		vector<int> tmp;
		vector<bool> seedMark2 = seedMark;

                while(countSuccess < numSuccessSamples){
			maxSeed[id] = false;
                        tmp = sampler(n_fail,zcover,seedMark2);
			numTotalEdge += n_fail+1;
			numSuccessEdge += 1;
			countTotal += n_fail+1;
			for (int i:tmp){
				if (link[i] < k){
					maxSeed[id] = true;
					break;
				}
			}
                        #pragma omp critical
                        {
                                  if (tmp.size() > 0){
                                            addHyperedge(tmp);
                                  }
                                  countSuccess += 1;
                                  if (maxSeed[id]){
                                            degree++;
                                  }
                        }
               }
        }

//	double approx = pow((1-1.0/k),k);
	cout << "Degree: " << sheap.getMaxDegree() << " " << degree << " " << countSuccess << " " << countTotal << " " << numSuccessSamples << " " << numTotalSamples << endl;

      if (degree >= lambda1){
		double epsilon_1 = (deg*normalized/numTotalSamples)/(degree*normalized/countTotal) - 1;
		cout << "Epsilon_1 = " << epsilon_1 << endl;
		double epsilon_2 = epsilon*sqrt(normalized*(1+1.0/3)/(degree*normalized*pow(2,iter-1)/countTotal));
		cout << "Epsilon_2 = " << epsilon_2 << " " << epsilon*sqrt(normalized*(1+epsilon)/(degree*normalized*pow(2,iter-1)/countTotal)) << " " << pow(2,iter-1) << " " << pow(3,iter-1) << endl;
		double epsilon_3 = epsilon*sqrt(normalized*(1+1.0/3)*(1-1/exp(1)-epsilon)/((1+epsilon/3)*degree*normalized*pow(2,iter-1)/countTotal));
		cout << "Epsilon_3 = " << epsilon_3 << endl;
		cout << "Epsilon_t = " << (epsilon_1 + epsilon_2 + epsilon_1*epsilon_2)*(1-1/exp(1)-epsilon) + epsilon_3*(1-1/exp(1)) << endl;
		
		if (do_check){
                        dualB = deg/dualB;
			cout << "Optimality Guarantee: " << dualB - (epsilon_1 + epsilon_2 + epsilon_1*epsilon_2)*(dualB-epsilon) - epsilon_3*dualB << " " << dualB << " " << deg << " " << deg/dualB << endl;
		}

              if ((epsilon_1 + epsilon_2 + epsilon_1*epsilon_2)*(1-1/exp(1)-epsilon) + epsilon_3*(1-1/exp(1)) <= epsilon){
			cout << "Estimated value: " << (degree*normalized/countTotal) << endl;
                        return true;
              }
      }
      return false;
}
*/

/*
bool HyperGraph::MSz(std::function< std::vector<int32_t>(int32_t &,int &, std::vector<bool> &, std::vector<bool> &) > sampler, std::vector<int> & seeds, int t, double & deg, float epsilon, float delta, long long int & numSuccessSamples, long long int & numTotalSamples, double normalized, double c, double zthreshold, vector<bool> & cseedmark)
{
	z = zthreshold;
	cout << "z = " << z << endl;
	eager_thres = z/numSeeds;
	double checkpoint = 1;
	long long countSuccess = 0;
        long long int countTotal = 0;
        int k = numSeeds;
	double epsilon0 = epsilon;
	if (epsilon < 1.0/3){
		epsilon0 = 1.0/3;
	}

	double f = (log(7/delta)+lgamma(numNodes+1)-lgamma(k+1)-lgamma(numNodes-k+1))/(log(7*log2(numNodes)/delta));
        double lambda1 = 1+(1+epsilon0)*(2+2*epsilon0/3)*log(3*log2(f)/delta)/(epsilon0*epsilon0);

        bool maxSeed = false;

        int32_t n_fail = 0;
	int degree = 0;

	while (getNumSelected() < k){
		vector<int> tmp;

                while(countTotal < checkpoint && getNumSelected() < k){
                       	maxSeed = false;
			tmp = sampler(n_fail,zcover,seedMark,cseedmark);
			if (tmp.size() == 0){
				break;
			}
			for (int i:tmp){
				if (cseedmark[i]){
                	                maxSeed = true;
                        	        break;
                                }
	                }
                        if (tmp.size() > 0){
                               	addHyperedge(tmp);
			}
        	        numTotalEdge += n_fail+1;
	                numSuccessEdge += 1;
        	        countTotal += n_fail+1;
	                countSuccess += 1;
        	        if (maxSeed){
                	        degree++;
                        }
               	}

        	cout << "Print degree and total samples: " << zcover << " " << degree << " " << countSuccess << " " << countTotal << " " << deg << " " << numTotalSamples << endl;

      		if (degree >= lambda1){
	                double epsilon_1 = (deg*1.0/numTotalSamples)/(degree*1.0/countTotal) - 1;
        	        cout << "Epsilon_1 = " << epsilon_1 << " " << deg << " " << numTotalSamples << " " << degree << " " << countTotal << endl;
                	double epsilon_2 = sqrt((1+epsilon0)*log(3*log2(zmax/zmin)*log(Tmax)/log(1+c))/(degree));
	                cout << "Epsilon_2 = " << epsilon_2 << endl;
        	        double epsilon_3 = sqrt((1+epsilon0)*(1-1/exp(1)-epsilon)*countTotal*log(3*log2(zmax/zmin)*Tmax)/((1+epsilon/3)*degree*numTotalSamples));
                	cout << "Epsilon_3 = " << epsilon_3 << endl;
	                cout << "Epsilon_t = " << (epsilon_1 + epsilon_2 + epsilon_1*epsilon_2)*(1-1/exp(1)-epsilon) + epsilon_3*(1-1/exp(1)) << endl;

              		if ((epsilon_1 + epsilon_2 + epsilon_1*epsilon_2)*(1-1/exp(1)-epsilon) + epsilon_3*(1-1/exp(1)) <= epsilon){
				cout << "Estimated value: " << (degree*normalized/countTotal) << endl;
				numSuccessSamples = numSuccessEdge;
			        numTotalSamples = countTotal;
				deg = degree;
                        	return true;
        		}
      		}
		while (checkpoint <= countTotal)
			checkpoint = checkpoint*(1+c);
	}

	for (unsigned int i = 0; i < seeds.size();++i){
                link[seeds[i]] = k;
        }

	seeds = getSelectedSeeds();
	deg = 0;

	if (getNumSelected() < k){
		vector<double> degreek(k+1,0);
		buildSeedSet(seeds,k-getNumSelected(),degreek);
		deg += zcover + degreek[k-getNumSelected()];
	} else {
		for (unsigned int i = 0; i < seeds.size();++i){
                	link[seeds[i]] = i;
        	}

		if (zcover + sheap.getSumK() < z){
                       	vector<int> tmp;

        	        while(zcover + sheap.getSumK() < z){
                                maxSeed = false;
                       	        tmp = sampler(n_fail,zcover,seedMark,cseedmark);
                               	if (tmp.size() == 0){
                                       	break;
	                        }

                               	for (int i:tmp){
                                       	if (link[i] < k){
                                               	maxSeed = true;
						deg++;
	                                        break;
        	                        }
                                }
        	                countSuccess += 1;
					
        	              	numTotalEdge += n_fail+1;
	                        numSuccessEdge += 1;
        	       	        countTotal += n_fail+1;
                                if (maxSeed){
                                        degree++;
                      	        } else {
					for (int i:tmp){
						sheap.increase(i);
		                        }

				}
                       	}
                }
		deg += zcover;
	}
	numSuccessSamples = numSuccessEdge;
	numTotalSamples = countTotal;
        cout << "Print degree and total samples: " << zcover << " " << degree << " " << countSuccess << " " << countTotal << " " << deg << " " << numTotalSamples << endl;
	return false;
}

void HyperGraph::SRA(std::function< std::vector<int32_t>(int32_t &,int &, std::vector<bool> &, std::vector<bool> &) > sampler, std::vector<int> & seeds, int t, float epsilon, float delta, double normalized)
{
        double deg = 1;
        long long int numSuccessSamples = 1;
        long long int numTotalSamples = 1;
        long long int samples = 0;
        clock_t start = clock();
        vector<bool> cseedmark(numNodes,false);
	double c = 0.1;

        MSz(sampler,seeds,t,deg,epsilon,delta,numSuccessSamples,numTotalSamples,normalized,c,zmax,cseedmark);
        cout << "Seed Nodes: ";
        for (unsigned int i = 0; i < seeds.size(); ++i){
                cout << seeds[i] << " ";
        }
        cout << endl;
        printf("Influence: %0.2lf\n",deg*normalized/numTotalSamples);
        cout << "Time: " << (float)(clock()-start)/CLOCKS_PER_SEC << endl;
        cout << "Memory: " << getMemValue()/1024.0 << endl;
        cout << "Samples: " << numTotalSamples << endl;

        std::vector<std::vector<int> >().swap(node_edge);
        std::vector<std::vector<int> >().swap(edge_node);
        if (do_check){
                std::vector<std::vector<int> >().swap(all_node_edge);
                std::vector<std::vector<int> >().swap(all_edge_node);
        }
        std::vector<int>().swap(nodedegree);
        std::vector<int>().swap(selectedSeed);
        std::vector<bool>().swap(seedMark);
}

void HyperGraph::DTA(std::function< std::vector<int32_t>(int32_t &,int &, std::vector<bool> &, std::vector<bool> &) > sampler, std::vector<int> & seeds, int t, float epsilon, float delta, double normalized)
{
	double deg = 1;
	long long int numSuccessSamples = 1;
	long long int numTotalSamples = 1;
	long long int samples = 0;
	double c = 0.1;
	clock_t start = clock();
	cout << "Try z = " << zmin << endl;
	vector<bool> cseedmark(numNodes,false);
	vector<int> seed_t;

	MSz(sampler,seeds,t,deg,epsilon,delta,numSuccessSamples,numTotalSamples,normalized,c,zmin,cseedmark);
	int numIter = ceil(log2(zmax/zmin));
	for (int i = 1; i < numIter; ++i){
		samples += numTotalSamples;
		cout << "Try z = " << pow(2,i)*zmin << endl;
		reinitialize();
		for (int j = 0; j < seed_t.size(); ++j){
			cseedmark[seed_t[j]] = false;
		}
		for (int j = 0; j < seeds.size(); ++j){
			cseedmark[seeds[j]] = true;
		}
		seed_t = seeds;
		if (MSz(sampler,seeds,t,deg,epsilon,delta,numSuccessSamples,numTotalSamples,normalized,c,pow(2,i)*zmin,cseedmark)){
			break;
		}
	}
	
	cout << "Seed Nodes: ";
        for (unsigned int i = 0; i < seeds.size(); ++i){
                cout << seeds[i] << " ";
        }
        cout << endl;
        printf("Influence: %0.2lf\n",deg*normalized/numTotalSamples);
        cout << "Time: " << (float)(clock()-start)/CLOCKS_PER_SEC << endl;
        cout << "Memory: " << getMemValue()/1024.0 << endl;
	cout << "Samples: " << samples + numTotalSamples << endl;

        std::vector<std::vector<int> >().swap(node_edge);
        std::vector<std::vector<int> >().swap(edge_node);
	if (do_check){
		std::vector<std::vector<int> >().swap(all_node_edge);
	        std::vector<std::vector<int> >().swap(all_edge_node);
	}
        std::vector<int>().swap(nodedegree);
        std::vector<int>().swap(selectedSeed);
        std::vector<bool>().swap(seedMark);
}
*/

void HyperGraph::buildSeedSet(vector<int> & seeds, int k, vector<double> &degree)
{
        long long i;
        unsigned int j,l,maxInd;
        vector<int> e, nList;

        vector<int> nodeDegree(numNodes,0);
        vector<int> indx(numNodes,0);
        for (j = 0; j < numNodes; ++j){
                indx[j] = j;
                nodeDegree[j] = node_edge[j+1].size();
        }

        InfCost<int> hd(&nodeDegree[0]);
        MappedHeap<InfCost<int> > heap(indx,hd);
        long long numEdge = edge_node.size();

        vector<bool> edgeMark(numEdge, false);
        vector<bool> nodeMark(numNodes+1, true);

        double totalCost = 0;

        i=1;
        degree[0] = 0;
        cout << "Top degrees: ";
        while(totalCost < k && !heap.empty()){
                maxInd = heap.pop()+1;
                nodeMark[maxInd] = false;
                cout << maxInd << " " << nodeDegree[maxInd-1] << " ";

                totalCost++;

                e = node_edge[maxInd];

                degree[i] = degree[i-1]+nodeDegree[maxInd-1];

                seeds.push_back(maxInd);
                for (j = 0; j < e.size(); ++j){
                        if (edgeMark[e[j]]){
                                continue;
                        }

                        nList = edge_node[e[j]];
                        for (l = 0; l < nList.size(); ++l){
                                nodeDegree[nList[l]-1]--;
                                if (nodeMark[nList[l]]){
                                        heap.heapify(nList[l]-1);
                                }
                        }
                        edgeMark[e[j]] = true;
                }
                i++;
	}
        cout << endl << endl;

        vector<int>().swap(nodeDegree);
        vector<int>().swap(e);
        vector<int>().swap(nList);
        vector<int>().swap(indx);
        vector<bool>().swap(edgeMark);
}

double HyperGraph::increaseT(vector<int> & seeds, int k, vector<double> &ve, double alpha)
{
        vector<int> e, nList;
        unsigned int i,j;
        long long numEdge = all_edge_node.size();
        vector<int> deg(numNodes,0);
        vector<double> sumleft(numNodes,0);
        vector<double> ratio(numNodes,0);
        vector<bool> emark(numEdge,false);
        for (i = 0; i < numNodes; ++i){
                e = all_node_edge[i+1];
                for (int j:e){
                        sumleft[i] += ve[j];
                        if (ve[j] < 1-exp(-6)){
                                deg[i]++;
                                emark[j] = true;
                        }
                }
        }
        vector<int> indx(numNodes,0);
        for (j = 0; j < numNodes; ++j){
                indx[j] = j;
                if (deg[j] == 0){
                        ratio[j] = 0;
                } else {
                        ratio[j] = (alpha-sumleft[j])/deg[j];
                }
        }
        MinData<double> hd(&ratio[0]);
        MappedHeap<MinData<double> > heap(indx,hd);
        double maxincrease = 0;
        while (!heap.empty()){
                int cnode = heap.top();
                if (deg[cnode] == 0){
                        heap.pop();
                        continue;
                }

                e = all_node_edge[cnode+1];
                if (ratio[cnode] == 0){
                        heap.pop();
                        for (int i:e){
                                if (emark[i]){
                                        emark[i] = false;
                                        nList = all_edge_node[i];
                                        for (int j:nList){
                                                j--;
                                                deg[j]--;
                                                if (deg[j] == 0){
                                                        ratio[j] = 0;
                                                } else if (deg[j] > 0){
                                                        ratio[j] = (alpha-sumleft[j])/deg[j];
                                                }
                                                heap.heapify(j);
                                        }
                                }
                        }
                } else {
                        maxincrease = ratio[cnode];
                        for (int i:e){
                                if (emark[i]){
                                        if (maxincrease > (1-ve[i])){
                                                maxincrease = 1-ve[i];
                                        }
                                }
                        }
                        for (int i:e){
                                if (emark[i]){
                                        ve[i] += maxincrease;
                                        nList = all_edge_node[i];
                                        for (int j:nList){
                                                sumleft[j-1] += maxincrease;
						if (sumleft[j-1] > alpha){
	//						cout << "WARNING: exceed alpha! " << sumleft[j-1] << " " << alpha << " " <<  sumleft[j-1] - maxincrease << " " << ve[i] << " " << ratio[cnode] << endl;
						}
                                        }
                                }
                        }
                        if (maxincrease == ratio[cnode]){
                                heap.pop();
                                for (int i:e){
                                        if (emark[i]){
                                                emark[i] = false;
                                                nList = all_edge_node[i];
                                                for (int j:nList){
                                                        j--;
                                                        deg[j]--;
                                                        if (deg[j] == 0){
                                                                ratio[j] = 0;
                                                        } else if (deg[j] > 0){
                                                                ratio[j] = (alpha-sumleft[j])/deg[j];
                                                        }
                                                        heap.heapify(j);
                                                }
                                        }
                                }
                        } else {
                                for (int i:e){
                                        if (emark[i] && ve[i] >= 1-exp(-6)){
                                                emark[i] = false;
                                                nList = all_edge_node[i];
                                                for (int j:nList){
                                                        j--;
                                                        deg[j]--;
                                                        if (deg[j] == 0){
                                                                ratio[j] = 0;
                                                        } else if (deg[j] > 0){
                                                                ratio[j] = (alpha-sumleft[j])/deg[j];
                                                        }
                                                        heap.heapify(j);
                                                }
                                        } else if (emark[i]){
                                                nList = all_edge_node[i];
                                                for (int j:nList){
                                                        j--;
                                                        //deg[j]--;
                                                        if (deg[j] == 0){
                                                                ratio[j] = 0;
                                                        } else if (deg[j] > 0){
                                                                ratio[j] = (alpha-sumleft[j])/deg[j];
                                                        }
                                                        heap.heapify(j);
                                                }
                                        }
                                }
                        }
                }
        }

        double obj = 0;
        for (unsigned int i = 0; i < numEdge; ++i){
                obj += (1-ve[i]);
        }
        return obj;
}

double HyperGraph::decreaseT(vector<int> & seeds, int k, vector<double> &ve, double alpha)
{
        unsigned int vsize = ve.size();
        vector<int> vdegree(vsize,0);
        vector<double> cvalues(numNodes,0);
        vector<int> e,nList;

        double dtmp = 0;
        for (unsigned int i = 0; i < numNodes; ++i){
                cvalues[i] = 0;
                dtmp = 0;
                e = all_node_edge[i+1];
                for (int j:e){
                        dtmp += ve[j];
                }
                if (dtmp - alpha > 0){
                        cvalues[i] = dtmp - alpha;
                        for (int j:e){
                                if (ve[j] > 0){
                                        vdegree[j]++;
                                }
                        }
                }
        }

        vector<int> indx(vsize,0);
        for (unsigned int j = 0; j < vsize; ++j){
                indx[j] = j;
        }

        InfCost<int> hd(&vdegree[0]);
        MappedHeap<InfCost<int> > heap(indx,hd);

        int maxid = 0;
        double maxdecrease = 0;

        while (vdegree[heap.top()] > 0){
                maxid = heap.top();
                maxdecrease = ve[maxid];
                if (ve[maxid] <= 0){
                        heap.pop();
                        continue;
                }
                nList = all_edge_node[maxid];
                for (int j:nList){
                        if (cvalues[j-1] > 0 && cvalues[j-1] < maxdecrease){
                                maxdecrease = cvalues[j-1];
                        }
                }
                for (int j:nList){
                        if (cvalues[j-1] > 0){
                                cvalues[j-1] -= maxdecrease;
                                if (cvalues[j-1] <= exp(-10)){
                                        e = all_node_edge[j];
                                        for (int k:e){
                                                if (ve[k] > 0){
                                                        vdegree[k]--;
                                                        heap.heapify(k);
                                                }
                                        }
                                }
                        }
                }

                ve[maxid]-= maxdecrease;
                if (ve[maxid] <= exp(-10)){
                        vdegree[maxid] = 0;
                }
                heap.heapify(maxid);
        }

        double obj = 0;
        for (unsigned int i = 0; i < vsize; ++i){
                obj += (1-ve[i]);
        }
        return obj;
}

double HyperGraph::ILPBound(vector<int> & seeds, int k)
{
	double bound = 0;
	long long numEdge = all_edge_node.size();
	vector<int> e, nList;
	clock_t start = clock();

        IloEnv env;

        try{
                IloModel model(env);
                IloNumVarArray x(env);
                IloRangeArray cons(env);

                //int timeLimit = 60*60;
                //int threads = 8;

                IloCplex cplex(model);
                //cplex.setParam(IloCplex::TiLim,timeLimit);
                //cplex.setParam(IloCplex::Threads, threads);

                x.add(IloNumVar(env,0));
                for (int i = 1; i <= numEdge; ++i){
                        x.add(IloNumVar(env,0,1));
                }

                IloExpr ob(env);
                ob += k*x[0];
                for (int i = 1; i <= numEdge; ++i){
                        ob += 1-x[i];
                }
                model.add(IloMinimize(env,ob));

                int maxconssize = 0;
                int curconssize = 0;

                for (unsigned int i = 1; i <= numNodes; i++){
                        IloExpr constraint(env);
                        constraint += x[0];
                        e = all_node_edge[i];
                        curconssize = 0;
                        for (int j : e){
                                curconssize++;
                                constraint += (-1)*x[j+1];
                        }
                        if (curconssize > maxconssize)
                                maxconssize = curconssize;
                        cons.add(IloRange(env, 0, constraint));
                //        cons.add(IloRange(env, 0, constraint));
                }
                env.out() << "**** Max constraint size: " << maxconssize << endl << endl;
                model.add(cons);

                cplex.solve();
//                IloNumArray vals(env);
                env.out() << "Solution status = " << cplex.getStatus() << endl;
		bound = cplex.getObjValue();
                env.out() << "Solution value  = " << cplex.getObjValue() << endl;
//                env.out() << "**** Exact dual bound: " << degree[k]*1.0/cplex.getObjValue() << endl;
                env.out() << "Optimality gap = " << cplex.getMIPRelativeGap() << endl;
//                cplex.getValues(vals, x);
//                env.out() << "**** Variables: ";
//                for (unsigned int i = 0; i <= numEdge; ++i){
//                        env.out() << vals[i] << " ";
//                }
//                env.out() << endl << endl;
        } catch (IloException &e){
                cerr << "Concert exception caught: " << e << endl;
        } catch (...){
                cerr << "Unknown exception caught" << endl;
        }

        env.end();
	cout << "Time for computing ILP bound: " << (float)(clock()-start)/CLOCKS_PER_SEC << "s" << endl;

	return bound;
}

double HyperGraph::dualBound(vector<int> & seeds, int k)
{
	/*vector<int> seeds;
	vector<double> degree(k+1,0);
        long long i;
        unsigned int j,l,maxInd;
        vector<int> e, nList;

        vector<int> nodeDegree(numNodes,0);
        vector<int> indx(numNodes,0);
        for (j = 0; j < numNodes; ++j){
                indx[j] = j;
                nodeDegree[j] = all_node_edge[j+1].size();
        }

        InfCost<int> hd(&nodeDegree[0]);
        MappedHeap<InfCost<int> > heap(indx,hd);
        long long numEdge = all_edge_node.size();

        vector<bool> edgeMark(numEdge, false);
        vector<bool> nodeMark(numNodes+1, true);

        double totalCost = 0;

        i=1;
        degree[0] = 0;
        cout << "Top degrees: ";
        while(totalCost < k && !heap.empty()){
                maxInd = heap.pop()+1;
                nodeMark[maxInd] = false;
                cout << maxInd << " " << nodeDegree[maxInd-1] << " ";

                totalCost++;

                e = all_node_edge[maxInd];

                degree[i] = degree[i-1]+nodeDegree[maxInd-1];

                seeds.push_back(maxInd);
                for (j = 0; j < e.size(); ++j){
                        if (edgeMark[e[j]]){
                                continue;
                        }

                        nList = all_edge_node[e[j]];
                        for (l = 0; l < nList.size(); ++l){
                                nodeDegree[nList[l]-1]--;
                                if (nodeMark[nList[l]]){
                                        heap.heapify(nList[l]-1);
                                }
                        }
                        edgeMark[e[j]] = true;
                }
                i++;
        }
        cout << endl << endl;
	
	maxInd = heap.pop()+1;
	*/
        long long numEdge = all_edge_node.size();
	double nextDeg = getMaxDegree();
	vector<int> e, nList;

        vector<double> ve(numEdge,1);
        for (int i:seeds){
                e = all_node_edge[i];
                for (int j:e){
                        ve[j] = 0;
                }
        }

        double objbefore = 0;
        for (int i = 0; i < numEdge; ++i){
                objbefore += (1-ve[i]);
        }

//        cout << endl <<"**** Bound before: " << degree[k]*1.0/(objbefore+k*nextDeg) << " " << endl << endl;
        clock_t start = clock();

        double increase = 0;
        double stepsize = nextDeg;
        double bestvalue = objbefore+k*nextDeg;
        vector<double> ve_t;
        do{
                do{
                        ve_t = ve;
                        increase = bestvalue - increaseT(seeds,k,ve_t,nextDeg+stepsize)-k*(nextDeg+stepsize);
/*        		for (int i = 0; i < numNodes; ++i){
		                cursum = 0;
                		e = all_node_edge[i+1];
		                for (int j:e){
                		        cursum += ve_t[j];
		                }
                		if (cursum > nextDeg+stepsize + exp(-6)){
		                        cout << "$$$$$$$$$$$$$$ Alert: " << cursum << " " << nextDeg+stepsize << endl;
                		}
		        }*/
                        cout << "Increase " << increase << endl;
                        if (increase < 0){
				ve_t = ve;
                                cout << stepsize << " " << bestvalue << " " << nextDeg << endl;
                                increase = bestvalue - decreaseT(seeds,k,ve_t,nextDeg-stepsize)-k*(nextDeg-stepsize);
			/*        for (int i = 0; i < numNodes; ++i){
			                cursum = 0;
			                e = all_node_edge[i+1];
			                for (int j:e){
                        			cursum += ve_t[j];
			                }
			                if (cursum > nextDeg-stepsize+exp(-6)){
                        			cout << "$$$$$$$$$$$$$$ Alert: " << cursum << " " << nextDeg -stepsize << endl;
			                }
        			}*/
                                cout << "Decrease " << increase << " " << stepsize << endl;
                                if (increase > 0){
                                        nextDeg -= stepsize;
                                        ve = ve_t;
                                        bestvalue -= increase;
                                }
                        } else {
                                nextDeg += stepsize;
                                ve = ve_t;
				bestvalue -= increase;
                        }
                } while (increase > exp(-6));
                stepsize /= 2.0;
        } while (stepsize > 0.000001);

        vector<double>().swap(ve_t);

        cout << "Time for increase/decrease bound: " << (float)(clock()-start)/CLOCKS_PER_SEC << "s" << endl;
        //cout << endl << "**** Bound by increase/decrease T: " << degree[k]*1.0/bestvalue << " " << endl << endl;

        // vector<int>().swap(nodeDegree);
        vector<int>().swap(e);
        vector<int>().swap(nList);
        // vector<int>().swap(indx);
        // vector<bool>().swap(edgeMark);
	
	return bestvalue;
}
/*
long long HyperGraph::addHyperedge(std::function< std::vector<int32_t>(int32_t &,int &, std::vector<bool> &) > sampler, int t, long long num)
{
        omp_set_num_threads(t);
            
        long long iter = 0;
        int c = 100;

        #pragma omp parallel
        {
                vector<int> tmp;
		int32_t n_fail = 0;
                while (iter < num){
                        for (int i = 0; i < c; ++i){
                                tmp = sampler(n_fail,zcover,seedMark);
				numTotalEdge += n_fail+1;
				numSuccessEdge += 1;
                                if (tmp.size() > 0){
                                      	#pragma omp critical
					{
                                              	addHyperedge(tmp);
					}
                                }
                        }

                        #pragma omp atomic
                        iter += c;
                }
        }
        return getNumEdge();
}

void HyperGraph::runMaxCover_SRA(std::function< std::vector<int32_t>(int32_t &,int &, std::vector<bool> &) > sampler, double epsilon, double delta, int k, int t, double normalized, int zscale)
{
	vector<double> degree(k+1,0);

        vector<int> seeds;

        double f = (log(6/delta)+lgamma(numNodes+1)-lgamma(k+1)-lgamma(numNodes-k+1))*normalized/(k*log(6*log2(numNodes)/delta));
        double lambda1 = (2+2/3)*log(3*log2(f)/delta);
        long long int totalSamples = (long long int)lambda1;
        cout << lambda1 << " " << totalSamples << endl;

        double nmax = (2+2*epsilon/3)*(lgamma(numNodes+1)-lgamma(k+1)-lgamma(numNodes-k+1) + log(7/delta))*normalized/(epsilon*epsilon*k);

        clock_t start = clock();
        cout << totalSamples << " " << nmax << " " << lgamma(numNodes+1) << " " << lgamma(k+1) << " " << lgamma(numNodes-k+1) << endl;

	do{
                cout << "Adding " << totalSamples << " hyperedge" << endl;
                addHyperedge(sampler,t,totalSamples);
                seeds = getSelectedSeeds();
                buildSeedSet(seeds,k-getNumSelected(),degree);
                degree[k] = degree[k-getNumSelected()]+getZCover();
                totalSamples = numSuccessEdge;
        } while (degree[k] < lambda);

        cout << "Seed Nodes: ";
	for (unsigned int i = 0; i < seeds.size(); ++i){
                cout << seeds[i] << " ";
        }
        cout << endl;
        printf("Influence: %0.2lf\n",degree[k]*normalized/numTotalEdge);
        cout << "Time: " << (float)(clock()-start)/CLOCKS_PER_SEC << endl;
        cout << "Memory: " << getMemValue()/1024.0 << endl;

	std::vector<std::vector<int> >().swap(node_edge);
	std::vector<std::vector<int> >().swap(edge_node);
        std::vector<int>().swap(nodedegree);
        std::vector<int>().swap(selectedSeed);
        std::vector<bool>().swap(seedMark);
}

bool HyperGraph::runMaxCover_SSAz(std::function< std::vector<int32_t>(int32_t &,int &, std::vector<bool> &) > sampler, double epsilon, double delta, int k, int t, double normalized, int zscale)
{
	bool stop = true;

	z /= zscale;

	vector<double> degree(k+1,0);

        vector<int> seeds;

        double f = (log(6/delta)+lgamma(numNodes+1)-lgamma(k+1)-lgamma(numNodes-k+1))*normalized/(k*log(6*log2(numNodes)/delta));
        double lambda1 = (2+2/3)*log(3*log2(f)/delta)/(epsilon*epsilon);
        long long int totalSamples = (long long int)lambda1;
        cout << lambda1 << " " << totalSamples << endl;

        int iter = 1;

        addHyperedge(sampler,t,totalSamples);

        double nmax = (2+2*epsilon/3)*(lgamma(numNodes+1)-lgamma(k+1)-lgamma(numNodes-k+1) + log(7/delta))*normalized/(epsilon*epsilon*k);

        clock_t start = clock();
        cout << totalSamples << " " << nmax << " " << lgamma(numNodes+1) << " " << lgamma(k+1) << " " << lgamma(numNodes-k+1) << endl;
        while (degree[k] < lambda){
                totalSamples = numTotalEdge;
                seeds = getSelectedSeeds();
                buildSeedSet(seeds,k-getNumSelected(),degree);
                degree[k] = degree[k-getNumSelected()]+getZCover();
		if (degree[k] >= z*k/2 && getNumSelected() > 0){
                        stop = false;
                        break;
                }
                cout << "Total Samples: " << totalSamples << " " << degree[k] << " " << numSuccessEdge << " " << numTotalEdge << endl;
                if (calculateInfluence(sampler,seeds,t,degree[k],epsilon,delta,numSuccessEdge,numTotalEdge,iter,normalized)){
                        break;
                }
                iter++;
        }
	if (stop){
	        cout << "Seed Nodes: ";
		for (unsigned int i = 0; i < seeds.size(); ++i){
                	cout << seeds[i] << " ";
	        }
        	cout << endl;
	        printf("Influence: %0.2lf\n",degree[k]*normalized/totalSamples);
        	cout << "Time: " << (float)(clock()-start)/CLOCKS_PER_SEC << endl;
	        cout << "Memory: " << getMemValue()/1024.0 << endl;

		std::vector<std::vector<int> >().swap(node_edge);
		std::vector<std::vector<int> >().swap(edge_node);
        	std::vector<int>().swap(nodedegree);
	        std::vector<int>().swap(selectedSeed);
        	std::vector<bool>().swap(seedMark);
	}
	return stop;
}

double HyperGraph::getZ()
{
	return z;
}
*/

/*
* convert from an integer to a string
*/
string intToStr(int i) {
        stringstream ss;
        ss << i;
        return ss.str();
}

/*
* convert from a strong to an integer
*/
unsigned int strToInt(string s) {
        unsigned int i;
        istringstream myStream(s);

        if (myStream>>i) {
                return i;
        } else {
                cout << "String " << s << " is not a number." << endl;
                return atoi(s.c_str());
        }
        return i;
}
