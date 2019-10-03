#include "IMSampler.h"

inline int IM_LTSampler::randIndex_lin(const vector<UI> &w, unsigned int si)
{
        UI ranNum = sfmt_genrand_uint32(&sfmtSeed);
        if (si <= 1 || ranNum > w[si - 1])
                return -1;

        for (unsigned int i = 1; i < si; ++i){
                if (ranNum <= w[i])
                        return i;
        }
        return -1;
}

inline int IM_LTSampler::randIndex_bin(const vector<UI> &w, unsigned int si)
{
        UI ran = sfmt_genrand_uint32(&sfmtSeed);
        if (si <= 1 || ran > w[si - 1])
                return -1;
        int left = 1;
        int right = si - 1;
        int prob;
        for (unsigned int i = 0; i < si; ++i){
                prob = (left + right)/2;
                if (w[prob - 1] > ran){
                        right = prob - 1;
                        continue;
                }
                if (w[prob] <= ran){
                        left = prob + 1;
                        continue;
                }
                break;
        }
        return prob;
}

IM_LTSampler::IM_LTSampler(Graph &gr, bool c):g(gr),do_check(c)
{
        numNodes = g.getSize();
        visit = vector<bool> (numNodes+1,false);
        visit_mark = vector<int> (numNodes+1,0);
        num_marked = 0;
        sfmt_init_gen_rand(&sfmtSeed, rand());
}

vector<int32_t> IM_LTSampler::polling(int32_t &n_fail, int & zcover, vector<bool> &seedMark, vector<bool> &cSeedMark)
{
	zcover = 0;
        n_fail = 0;
        int i;
        bool t1 = false;
	bool t2 = false;
        unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%numNodes+1;
        num_marked = 0;
        int ind;
        for (i = 0; i < numNodes; ++i){

                if (visit[cur] == true) break;
                visit[cur] = true;
                visit_mark[num_marked] = cur;
                num_marked++;

                if (seedMark[cur]){
                        t1=true;
			zcover = zcover | 2;
			if (!do_check && t2){
                        	break;
			}
                }
		if (cSeedMark[cur]){
			t2 = true;
			zcover = zcover | 1;
			if (!do_check && t1){
                                break;
                        }
		}

                if (g.weights[cur].size() >= 32)
                        ind = randIndex_bin(g.weights[cur],g.node_deg[cur]);
                else
                        ind = randIndex_lin(g.weights[cur],g.node_deg[cur]);

                if (ind == -1)
                        break;

                cur = g.adjList[cur][ind-1];
        }
	/*if (t1 && !t2){
	}
	if (!t1 && t2){
	}*/
        for (i = 0; i < num_marked; ++i){
                visit[visit_mark[i]]=false;
        }
        /*if (!do_check && t1 && t2){
                return vector<int32_t>();
        }
	*/
        return vector<int32_t>(visit_mark.begin(),visit_mark.begin()+num_marked);
}


IM_ICSampler::IM_ICSampler(Graph &gr, bool c):g(gr),do_check(c)
{
        numNodes = g.getSize();
        visit = vector<bool> (numNodes+1,false);
        visit_mark = vector<int> (numNodes+1,0);
        num_marked = 0;
        sfmt_init_gen_rand(&sfmtSeed, rand());
}

/*
 * * polling process under IC model
 * */
vector<int32_t> IM_ICSampler::polling(int32_t &n_fail, int &zcover, vector<bool> &seedMark, vector<bool> &cSeedMark)
{
	zcover = 0;
        n_fail = 0;
        int i;
        unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%(numNodes)+1;
        int curPos=0;
        num_marked=1;
        visit[cur] = true;
        visit_mark[0] = cur;
        bool t1 = false;
	bool t2 = false;
        while(curPos < num_marked){
                cur = visit_mark[curPos];
                curPos++;
                if (seedMark[cur]){
                        t1=true;
			zcover = zcover | 2;
			if (!do_check && t2){
	                        break;
			}
                }
		if (cSeedMark[cur]){
                        t2=true;
			zcover = zcover | 1;
                        if (!do_check && t1){
                                break;
                        }
                }

                const vector<UI> &w=g.weights[cur];
                const vector<int> &neigh = g.adjList[cur];
                for (i = 0; i < g.node_deg[cur]; ++i){
                        if (sfmt_genrand_uint32(&sfmtSeed) <  w[i+1]){
                                if (!visit[neigh[i]]){
                                        visit[neigh[i]] = true;
                                        visit_mark[num_marked]=neigh[i];
                                        num_marked++;
                                }
                        }
                }
        }
        for(i = 0; i < num_marked;++i){
                visit[visit_mark[i]]=false;
        }
        /*if (t && !do_check){
                num_marked = 0;
        }*/
        return vector<int32_t>(visit_mark.begin(),visit_mark.begin()+num_marked);
}
