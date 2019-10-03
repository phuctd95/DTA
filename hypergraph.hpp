#ifndef _HYPERGRAPH_H_
#define _HYPERGRAPH_H_

#include "rwgraph.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstring>
#include <random>
#include "mappedHeap.hpp"
#include "HeapData.hpp"

using namespace std;

/*
* find seed nodes procedure using greedy algorithm
*/
void buildSeedSet(HyperGraph & hg, vector<int> & seeds, unsigned int n, int k, vector<double> &degree)
{	
	double ob;
	long long i;
	unsigned int j,l,maxInd;
	vector<int> e, nList;

	vector<int> nodeDegree(n,0);
	vector<int> indx(n,0);
	for (j = 0; j < n; ++j){
		indx[j] = j;
		nodeDegree[j] = hg.getNode(j+1).size();
	}

	InfCost<int> hd(&nodeDegree[0]);
	MappedHeap<InfCost<int> > heap(indx,hd);
	long long numEdge = hg.getNumEdge();

	vector<bool> edgeMark(numEdge, false);
	vector<bool> nodeMark(n+1, true);
	
	double totalCost = 0;

	i=1;
	degree[0] = 0;
	cout << "Top degrees: ";
	while(totalCost < k && !heap.empty()){
		maxInd = heap.pop()+1;
		nodeMark[maxInd] = false;
		cout << maxInd << " " << nodeDegree[maxInd-1] << " ";

		totalCost++;

		e = hg.getNode(maxInd);
		
		degree[i] = degree[i-1]+nodeDegree[maxInd-1];
		
		seeds.push_back(maxInd);
		for (j = 0; j < e.size(); ++j){
			if (edgeMark[e[j]]){
				continue;
			}
	
			nList = hg.getEdge(e[j]);
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

	cout << "Degree of seeds: ";
	for (int i:seeds){
		cout << nodeDegree[i-1] << " ";
	}
	cout << endl;

	int leftover = 0;
        MappedHeap<InfCost<int> > heap2 = heap;
	cout << "Next degree: ";
        for (i = 0; i < k; ++i){
                if (heap2.empty()){
                        break;
                }
                maxInd = heap2.pop();
                leftover += nodeDegree[maxInd];
		cout << nodeDegree[maxInd] << " ";
        }
	cout << endl << endl;
        maxInd = heap.top()+1;
        int t = nodeDegree[maxInd-1];
        double dualbound = degree[k]+ k*nodeDegree[maxInd-1];
        cout << endl << endl << "**** Direct dual bound: " << degree[k]*1.0/dualbound << " " << degree[k] << " " << dualbound << " ****" << endl << endl;

        int maxedge = 0;
        int maxedgeid = 0;
        int tmpedge = 0;

	for (int i:seeds){
		e = hg.getNode(i);
		for (int nedge:e){
			maxedge = 0;
			nList = hg.getEdge(nedge);
			for (int enode:nList){
				if (maxedge < nodeDegree[enode-1]){
					maxedge = nodeDegree[enode-1];
				}
			}
			if (maxedge < t){
				dualbound -= 1;
				for (int enode:nList){
					nodeDegree[enode-1] += 1;
					if (nodeMark[enode]){
		                                heap.heapify(enode-1);
                        		}
				}
			}
		}
	}
	cout << endl << endl << "**** Improved dual bound: " << degree[k]*1.0/dualbound << " " << degree[k] << " " << dualbound << " ****" << endl << endl;

        while (true){
                maxedge = 0;
                maxInd = heap.top()+1;
                t = nodeDegree[maxInd-1];
                e = hg.getNode(maxInd);
                for (j = 0; j < e.size(); ++j){
                        if (edgeMark[e[j]]){
                                continue;
                        }
                        tmpedge = 0;
                        nList = hg.getEdge(e[j]);
                        for (l = 0; l < nList.size(); ++l){
                                if (nodeMark[nList[l]]){
                                        tmpedge++;
                                }
                        }
                        if (tmpedge > maxedge){
                                maxedge = tmpedge;
                                maxedgeid = e[j];
                        }
                }
                if (maxedge == 0){
                        break;
                }
                nList = hg.getEdge(maxedgeid);
                for (l = 0; l < nList.size(); ++l){
                        nodeDegree[nList[l]-1]--;
                        if (nodeMark[nList[l]]){
                                heap.heapify(nList[l]-1);
                        }
                }
                cout << "*** Heap: " << maxedge << " " << t << " " << nodeDegree[heap.top()] << " " << nodeDegree[maxInd-1] << " " << heap.top() << endl;
                edgeMark[maxedgeid] = true;
                if ((t-nodeDegree[heap.top()])*k < 1){
                        break;
                }
                cout << "Saving: " << (t-nodeDegree[heap.top()])*k << " " << maxedge << endl;
                dualbound = dualbound - (t-nodeDegree[heap.top()])*k + 1;
        }
        cout << endl << endl << "**** New dual bound: " << degree[k]*1.0/dualbound << " " << degree[k] << " " << dualbound << " ****" << endl << endl;
        vector<int>().swap(nodeDegree);
        vector<int>().swap(e);
        vector<int>().swap(nList);
        vector<int>().swap(indx);
        vector<bool>().swap(edgeMark);
        ob = degree[k]*1.0/(degree[k]+leftover);
        cout << endl << endl << "**** Online bound: " << ob << " ****" << endl << endl;
        double fixedbound = 1-pow((1-1.0/k),k);
        cout << "Fixed bound: " << fixedbound << endl;
        if (ob < fixedbound){
                ob = fixedbound;
        }
}

#endif
