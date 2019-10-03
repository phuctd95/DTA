#include "option.h"
#include "rwgraph.h"
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstring>
#include "sampler.h"
#include "IMSampler.h"

using namespace std;
using namespace std::placeholders;

void init_seeds(Sampler &sample)
{
    vector <uint64_t> s_init;
    for (int i = 0; i < 16; ++i)
        s_init.push_back(rand());
    sample.init_seed(rand()%16,s_init);
}

int main(int argc, char ** argv)
{
	srand(time(NULL));

	OptionParser op(argc, argv);
	if (!op.validCheck()){
		printf("Parameters error, please check the readme.txt file for correct format!\n");
		return -1;
	}
	char * inFile = op.getPara("-i");
	if (inFile == NULL){
		inFile = (char*)"network.bin";
	}

	char * outFile = op.getPara("-o");
	if (outFile == NULL){
		outFile = (char*)"network.seeds";
	}

	char * tmp = op.getPara("-epsilon");
	double epsilon = 0.1;
	if (tmp != NULL){
		epsilon = atof(tmp);
	}

	double delta = 0;
	tmp = op.getPara("-delta");
        if (tmp != NULL){
                delta = atof(tmp);
        }

	int k = 0;
	
	tmp = op.getPara("-k");
	if (tmp != NULL){
		k = atoi(tmp);
	}	

	int t = 1;
	tmp = op.getPara("-t");
	if (tmp != NULL){
		t = atoi(tmp);
	}
	
	char * prob = op.getPara("-prob");
	if (prob == NULL){
		cout << "Select an algorithm and run again!" << endl;
		return 1;
	}

	char * alg = op.getPara("-alg");
        if (alg == NULL){
                cout << "Select an algorithm and run again!" << endl;
                return 1;
        }
	
	int eager_mode = 1;
	char * emode = op.getPara("-emode");
	if (emode != NULL){
		if (strcmp(emode,"STRICT") == 0){
			eager_mode = 1;
		} else if (strcmp(emode,"EAGER") == 0){
			eager_mode = 2;
		} else if (strcmp(emode,"VEAGER") == 0){
			eager_mode = 3;
		} else {
			cout << "Invalid eager mode ====> set to default 'strict' mode" << endl;
		}
	}

	bool do_check = false;
	char * checked = op.getPara("-check");
	if (checked != NULL) do_check = true;

	if (strcmp(prob, "IM") == 0){
	        char * model = op.getPara("-m");
        	if (model == NULL)
                	model = (char *) "LT";

		if (strcmp(model, "LT") == 0){
			Graph g;
			g.readGraphLT(inFile);
			int n = g.getSize();
			HyperGraph hg(n,k,eager_mode,epsilon,delta,do_check);
			IM_LTSampler sampler(g,do_check);
			if (strcmp(alg, "SRA") == 0){
				hg.runMaxCover_SRA(std::bind(&IM_LTSampler::polling, &sampler, _1, _2, _3), epsilon,delta,k,t,n*1.0,1);
			} else if (strcmp(alg, "SSA_z") == 0){
				hg.runMaxCover_SSAz(std::bind(&IM_LTSampler::polling, &sampler, _1, _2, _3), epsilon,delta,k,t,n*1.0,1);
			} else if (strcmp(alg, "SSA_mul") == 0){
				double z = hg.getZ();
				for (int l = floor(log2(z))-2; l>= 0; l--){
					if (hg.runMaxCover_SSAz(std::bind(&IM_LTSampler::polling, &sampler, _1, _2, _3), epsilon,delta,k,t,n*1.0,l)){
						break;
					}
					hg = HyperGraph(n,k,eager_mode,epsilon,delta,do_check);
				}
			}
		} else if (strcmp(model, "IC") == 0){
			Graph g;
			g.readGraphIC(inFile);
			int n = g.getSize();
                        HyperGraph hg(n,k,eager_mode,epsilon,delta,do_check);
                        IM_ICSampler sampler(g,do_check);
			if (strcmp(alg, "SRA") == 0){
                        	hg.runMaxCover_SRA(std::bind(&IM_ICSampler::polling, &sampler, _1, _2, _3),epsilon,delta,k,t,n*1.0,1);
                        } else if (strcmp(alg, "SSA_z") == 0){
                                hg.runMaxCover_SSAz(std::bind(&IM_ICSampler::polling, &sampler, _1, _2, _3),epsilon,delta,k,t,n*1.0,1);
                        } else if (strcmp(alg, "SSA_mul") == 0){
                                double z = hg.getZ();
                                for (int l = floor(log2(z)); l>= 0; l--){
                                        if (hg.runMaxCover_SSAz(std::bind(&IM_ICSampler::polling, &sampler, _1, _2, _3), epsilon,delta,k,t,n*1.0,l)){
                                                break;
                                        }
					hg = HyperGraph(n,k,eager_mode,epsilon,delta,do_check);
                                }
                        }
		} else {
			printf("Incorrect model option!");
			return -1;
		}
	} else if (strcmp(prob, "BCM") == 0) {
		graph g;
		g.read_graph(inFile,true);
		Sampler sam(g);
		init_seeds(sam);
		int n = g.get_n_vertices();
		HyperGraph hg(n,k,eager_mode,epsilon,delta,do_check);
		
		if (strcmp(alg, "SRA") == 0){
                        hg.runMaxCover_SRA(std::bind(&Sampler::bc_polling, &sam, _1, _2, _3), epsilon,delta,k,t,n*1.0,1);
                } else if (strcmp(alg, "SSA_z") == 0){
                        hg.runMaxCover_SSAz(std::bind(&Sampler::bc_polling, &sam, _1, _2, _3), epsilon,delta,k,t,n*1.0,1);
                } else if (strcmp(alg, "SSA_mul") == 0){
                        double z = hg.getZ();
                        for (int l = floor(log2(z))-2; l>= 0; l--){
                                if (hg.runMaxCover_SSAz(std::bind(&Sampler::bc_polling, &sam, _1, _2, _3), epsilon,delta,k,t,n*1.0,l)){
                                        break;
                                }
				hg = HyperGraph(n,k,eager_mode,epsilon,delta,do_check);
                        }
                }
	} else if (strcmp(prob, "CCM") == 0) {
                graph g;
                g.read_graph(inFile,true);
                Sampler sam(g);
                init_seeds(sam);
                int n = g.get_n_vertices();
                HyperGraph hg(n,k,eager_mode,epsilon,delta,do_check);

		if (strcmp(alg, "SRA") == 0){
                        hg.runMaxCover_SRA(std::bind(&Sampler::cc_polling, &sam, _1, _2, _3), epsilon,delta,k,t,n*1.0,1);
                } else if (strcmp(alg, "SSA_z") == 0){
                        hg.runMaxCover_SSAz(std::bind(&Sampler::cc_polling, &sam, _1, _2, _3), epsilon,delta,k,t,n*1.0,1);
                } else if (strcmp(alg, "SSA_mul") == 0){
                        double z = hg.getZ();
                        for (int l = floor(log2(z))-2; l>= 0; l--){
                                if (hg.runMaxCover_SSAz(std::bind(&Sampler::cc_polling, &sam, _1, _2, _3), epsilon,delta,k,t,n*1.0,l)){
                                        break;
                                }
				hg = HyperGraph(n,k,eager_mode,epsilon,delta,do_check);
                        }
                }
	} else {
		cout << "Selected problem is invalid!" << endl;
	}
}
