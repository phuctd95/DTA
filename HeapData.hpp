#ifndef _HEAP_DATA_H_
#define _HEAP_DATA_H_

template<class T1>
struct InfCost{
	T1 *v1;
public:
	InfCost(T1 *u1):v1(u1){};
	
	bool operator() (int &i, int &j) const{
		return v1[i] < v1[j];
	}
};

template<class T>
struct MinData{
        T *v;
public:
        MinData(T *u):v(u){};

        bool operator() (int &i, int &j) const{
                return v[i] > v[j];
        }
};

template<class T1>
struct MaxData{
        T1 *v1;
public:
        MaxData(T1 *u1):v1(u1){};

        bool operator() (int &i, int &j) const{
                return v1[i] < v1[j];
        }
};

#endif
