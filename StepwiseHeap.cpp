#include "StepwiseHeap.h"

Stepwise_heap::Stepwise_heap(int n, int knode, int z)
{
	size = n;
	k = knode;

	degree = vector<int>(n+1,0);
	order = vector<int>(n+1,0);
	rorder = vector<int>(n+1,0);
	for (int i = 0; i < n+1; ++i){
		order[i] = i;
		rorder[i] = i;
	}

	buckethead = vector<int>(z+1,n+1);
	buckethead[0] = 0;
	sumK = 0;
}

void Stepwise_heap::swap(int t1, int t2)
{
	if (t1 != t2){
		order[rorder[t1]] += order[rorder[t2]];
		order[rorder[t2]] = order[rorder[t1]] - order[rorder[t2]];
		order[rorder[t1]] -= order[rorder[t2]];

		rorder[t1] += rorder[t2];
		rorder[t2] = rorder[t1] - rorder[t2];
		rorder[t1] -= rorder[t2];
	}
}

void Stepwise_heap::increase(int t)
{
	degree[t]++;
	swap(t,order[buckethead[degree[t]]-1]);
	buckethead[degree[t]]--;
	if (rorder[t] > size-k){
		sumK++;
	}
	#ifdef _DEBUG
        if (checkOrder()){
                cout << "Error increase()" << endl;
                exit(3);
        }
        #endif
}

void Stepwise_heap::decrease(int t)
{
	if (degree[t] > 0 && buckethead[degree[t]] > size-k){
		sumK--;
	}
	swap(t,order[buckethead[degree[t]]]);
	buckethead[degree[t]]++;
	degree[t]--;
	#ifdef _DEBUG
        if (checkOrder()){
                cout << "Error decrease()" << endl;
                exit(3);
        }
        #endif
}

int Stepwise_heap::getTop()
{
	/*int i = size;
	while (degree[order[i]] == -1){
		i--;
	}*/
	#ifdef _DEBUG
	if (checkOrder()){
		cout << "Error getTop" << endl;
		exit(3);
	}
	#endif
	return order[size];
}

int Stepwise_heap::getMaxDegree()
{
/*	int i = size;
        while (degree[order[i]] == -1){
                i--;
        }*/
	#ifdef _DEBUG
		if (checkOrder()){
                cout << "Error getMaxDegree" << endl;
                exit(3);
        }
        #endif
        return degree[order[size]];
}

int Stepwise_heap::getSumK()
{
	return sumK;
}

void Stepwise_heap::removeMaxNode()
{
	// to be implemented by calling decrease degree[top] times
}
