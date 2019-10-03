#ifndef SAMPLERF_H
#define SAMPLERF_H
#include <vector>
#include <string> 
#include "xorshift.h"

using namespace std;
class sampler_f
{
public:
	sampler_f();
	~sampler_f();
	void readF(char* infile);
	vector <int32_t> feature_polling(int32_t &n_fails, int &zcover, vector<bool> & seedMark, vector<bool> & cSeedMark);
	vector <int32_t> feature_probability_polling(int32_t &n_fails, int &zcover, vector<bool> & seedMark, vector<bool> & cSeedMark);
	void init_seed(int32_t p_init, vector <uint64_t> &s_init);
	int n,m;
private:
	xorshift rand_gen;
	vector <vector <double> > f;	
	
};





#endif
