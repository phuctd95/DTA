#ifndef _IMSAMPLER_H
#define _IMSAMPLER_H

#include "sfmt/SFMT.h"
#include <cmath>
#include <vector>
#include "rwgraph.h"

using namespace std;

class Graph;

class IM_ICSampler{
        private:
                std::vector<bool> visit;
                std::vector<int> visit_mark;
                int num_marked;
                sfmt_t sfmtSeed;
                Graph &g;
                int numNodes;
		bool do_check = false;

        public:
                IM_ICSampler(Graph &gr, bool c);
                std::vector<int32_t> polling(int32_t &n_fail, int &zcover, std::vector<bool> &seedMark, std::vector<bool> &cSeedMark);
};

class IM_LTSampler{
        private:
                std::vector<bool> visit;
                std::vector<int> visit_mark;
                int num_marked;
                sfmt_t sfmtSeed;
                Graph &g;
                int numNodes;
		bool do_check = false;
		
                inline int randIndex_bin(const std::vector<UI> &w, unsigned int si);
                inline int randIndex_lin(const std::vector<UI> &w, unsigned int si);

        public:
                IM_LTSampler(Graph &gr, bool c);
                std::vector<int32_t> polling(int32_t &n_fail, int &zcover, std::vector<bool> &seedMark, std::vector<bool> &cSeedMark);
};

#endif
