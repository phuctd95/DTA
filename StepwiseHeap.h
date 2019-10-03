/* Stepwise Heap
 *
 *
 */

#ifndef STEPWISE_QUEUE_H_
#define STEPWISE_QUEUE_H_

// #define _DEBUG

#include<vector>
#include<iostream>

using namespace std;


class Stepwise_heap {
	private:
		int size, k;

		vector<int> order;
		vector<int> rorder;
		vector<int> buckethead;
	
		int sumK;

	public:
		vector<int> degree;
		// construct the stepwise heap with n elements to select top k nodes and threshold z of maximum value of each elememts		
		Stepwise_heap(){}
		Stepwise_heap(int n, int k, int z);
		// increase the value of element t by 1
		void increase(int t);
		// decrease the value of element t by 1
		void decrease(int t);
		// get the element with highest value
		int getTop();
		// get the highest value
		int getMaxDegree();
		// get sum of top k elements
		int getSumK();
		// remove the top node
		void removeMaxNode();
		// swap two elements t1 and t2
		void swap(int t1, int t2);
		// check matching degree
		bool checkDegree(vector<int> &d){
			for (int i = 1; i <= size; ++i){
				if (d[i] != degree[i]){
					cout << d[i] << " " << degree[i] << " " << i << endl;
					return false;
				}
			}
			return true;
		}

		void printOrder(){
			cout << "Order: ";
			for (int i = 0; i <= size; ++i){
				cout << order[i] << ":" << degree[order[i]] << " ";
			}
			cout << endl;
		}

		int checkOrder(){
                        for (int i = 0; i < size-1; ++i){
                                if (degree[order[i]] > degree[order[i+1]]){
					printOrder();
					cout << "Order error " << i << " " << order[i] << " " << degree[order[i]] << endl;
					return 1;
				}
                        }
			int topk = 0;
			for (int i = 0; i < k; ++i){
				topk += degree[order[size-i]];
			}
			if (topk != sumK){
				printOrder();
				cout << "Top k error " << topk << " " << sumK << endl;
				return 2;
			}	
                        return 0;
                }

		void printDegree(){
			cout << "Degree: ";
                        for (int i = 0; i <= size; ++i){
                                cout << degree[i] << " ";
                        }
                        cout << endl;
                }
};
#ifdef _DEBUG

class _CheckStepwiseHeap {
	public:
		_CheckStepwiseHeap(int n) {
			cout << "Created heap!" << endl;
			Stepwise_heap h(n, 2, 5);
			h.printOrder();
			h.printDegree();
			h.increase(1);
			h.printOrder();
			h.printDegree();
			h.increase(1);
			h.printOrder();
			h.printDegree();
			h.increase(6);
			h.printOrder();
			h.printDegree();
			h.decrease(1);
			h.printOrder();
			h.printDegree();
			h.decrease(1);
			h.printOrder();
			h.printDegree();
			
			cout << "Top degree: " << h.getMaxDegree() << " of node " << h.getTop() << endl;
		}
};
#endif

#endif
