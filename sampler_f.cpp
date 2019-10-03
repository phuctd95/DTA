#include "sampler_f.h"
#include <fstream>
#include <sstream>
#include <iostream>
sampler_f::sampler_f()
{
	n = 0;
}
sampler_f::~sampler_f()
{

}

void sampler_f::readF(char* infile)
{
	ifstream fi(infile);
	string s;
	double x;
	while (getline(fi,s))
	{
		stringstream ss(s);
		f.push_back(vector<double>(0));
		while (ss >> x)
			f[n].push_back(x);
		m = f[n].size();
		++n;
	}
	cerr << n << " " << m << endl;
	fi.close();
}

void sampler_f::init_seed(int32_t p_init, vector <uint64_t> &s_init)
{
    rand_gen.init_seed(p_init,s_init);
}

vector <int> sampler_f::feature_polling(int32_t &n_fails, int &zcover, vector<bool> & seedMark, vector<bool> & cSeedMark)
{
	int s,t;
	s = rand_gen.get_max(n);
    while(1)
    {
        t = rand_gen.get_max(n);
        if (s != t) break;
    }
    vector <int> res;
    for (int i = 0; i <m; ++i)
        if (f[s][i] == 1 && f[t][i] == 1)
    	// if (f[s][i] == 1)
    		res.push_back(i+1);
    n_fails = 0;
	zcover = 0;
    bool t1 = false, t2 = false;
    for (int i : res){
        if (seedMark[i]){
            zcover = zcover | 2;
            t1 = true;
            if (t2)
                break;
        }
        if (cSeedMark[i]){
            zcover = zcover | 1;
            t2 = true;
            if (t1)
                break;
        }
    }
    return res;
}

vector <int> sampler_f::feature_probability_polling(int32_t &n_fails, int &zcover, vector<bool> & seedMark, vector<bool> & cSeedMark)
{
    int s,t;
    s = rand_gen.get_max(n);
    while(1)
    {
        t = rand_gen.get_max(n);
        if (s != t) break;
    }
    vector <int> res;
    double p;
    for (int i = 0; i <m; ++i)
    {
        p = f[s][i] * f[t][i];
        if (rand_gen.get_real() <= p)
            res.push_back(i+1);
    }
    n_fails = 0;
    zcover = 0;
    bool t1 = false, t2 = false;
    for (int i : res){
        if (seedMark[i]){
            zcover = zcover | 2;
            t1 = true;
            if (t2)
                break;
        }
        if (cSeedMark[i]){
            zcover = zcover | 1;
            t2 = true;
            if (t1)
                break;
        }
    }
    return res;
}