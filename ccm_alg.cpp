#include "option.h"
#include "rwgraph.h"
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstring>
#include "getMem.hpp"
#include "sampler.h"
#include "IMSampler.h"

using namespace std;
using namespace std::placeholders;

HyperGraph hg;

std::vector<int> selectedSeed;
std::vector<bool> seedMark;
double z, zmax,zmin,Tmax;
int zcover;
long long zfull_sketch, zmax_sketch;

unsigned int numNodes;
unsigned int numSeeds;
double eager_thres;
double lambda;
long long full_sketch;
bool do_check = false;

void init_seeds(Sampler &sample)
{
    vector <uint64_t> s_init;
    for (int i = 0; i < 16; ++i)
        s_init.push_back(rand());
    sample.init_seed(rand()%16,s_init);
}

void initialize(unsigned int n, unsigned int numS, double epsilon, double delta, bool c, double normalized)
{
    seedMark = vector<bool>(n+1,false);
        numNodes = n;
        zcover = 0;
        numSeeds = numS;

        double lognk = lgamma(numNodes+1) - lgamma(numSeeds+1) - lgamma(numNodes-numSeeds-1);
        double epsilon2 = sqrt((1-1/exp(1))*log(lognk*7/delta))/((1-1/exp(1))*sqrt(log(7/delta))+sqrt((1-1/exp(1))*log(lognk*7/delta)))*epsilon;
        lambda = (2+2*epsilon2/3)*(lognk + log(7/delta))/(epsilon2*epsilon2);

        zmax = (1+epsilon)/(1-1/exp(1))*lambda;
        cout << "Lambda threshold: " << lambda << " " << zmax << endl;
        zmin = log(1/delta)/(epsilon*epsilon)/numSeeds;

        double f = (log(6/delta)+lgamma(numNodes+1)-lgamma(numSeeds+1)-lgamma(numNodes-numSeeds+1))*normalized/(numSeeds*log(6*log2(numNodes)/delta));

        Tmax = (2+2/3)*log(3*log2(f)/delta);
}

void reinitialize()
{
    selectedSeed = vector<int>();

        seedMark = vector<bool>(numNodes+1,false);
        zcover = 0;
        //full_sketch = 0;
    hg.reinitialize(numSeeds,ceil(zmax));
}

double lowerBound(double u, double i, double n, double a){
    double ci = u/i;
        double c = ci*n/i;
    double lbound = a-(a+1)*ci/(3+ci);
        double lb = (-(-2*a+2*ci/3+2*a*ci/3-2*c)-sqrt((-2*a+2*ci/3+2*a*ci/3-2*c)*(-2*a+2*ci/3+2*a*ci/3-2*c)-4*(a*a-2*a*ci/3)*(1-2*ci/3+2*c)))/(2*(1-2*ci/3+2*c));
        if (lbound > lb){
                lbound = lb;
        }
    return lbound;
}

double upperBound(double u, double i, double n, double a){
    double ci = u/i;
    double c = ci*n/i;
    double ubound = a + (a+1)*ci/(3-ci);
    cout << "Computing upper bound: " << ci << " " << c << " " << ubound << endl;
        double ub = (-(-2*a-2*ci/3-2*a*ci/3-2*c)+sqrt((-2*a-2*ci/3-2*a*ci/3-2*c)*(-2*a-2*ci/3-2*a*ci/3-2*c) - 4*(a*a+2*a*ci/3)*(1+2*ci/3+2*c)))/(2*(1+2*ci/3+2*c));
        if (ubound < ub){
                ubound = ub;
        }
    return ubound;
}

bool MSz(std::function< std::vector<int32_t>(int32_t &,int &, std::vector<bool> &, std::vector<bool>&) > sampler, std::vector<int> & seeds, int t, double & deg, float epsilon, float delta, long long int & numSuccessSamples, long long int & numTotalSamples, double normalized, double c, double thresholdz, std::vector<bool> &cseedmark, double ubound)
{
    z = thresholdz;
        cout << "z = " << z << endl;
        eager_thres = z/numSeeds;
        double checkpoint = 1;
        long long countSuccess = 0;
        long long max_sketch = 0;
        long long int countTotal = 0;
        int k = numSeeds;
        double epsilon0 = epsilon;
        cout << "init size" << hg.getSketchSize() << endl;
        if (epsilon < 1.0/3){
                epsilon0 = 1.0/3;
        }

        double f = (log(7/delta)+lgamma(numNodes+1)-lgamma(k+1)-lgamma(numNodes-k+1))/(log(7*log2(numNodes)/delta));
        double lambda1 = 1+(1+epsilon0)*(2+2*epsilon0/3)*log(3*log2(f)/delta)/(epsilon0*epsilon0);

        int32_t n_fail = 0;
        int degree = 0;
    int cc = 0;

        while (selectedSeed.size() < k){
                vector<int> tmp;

                while(countTotal < checkpoint && selectedSeed.size() < k){
                        tmp = sampler(n_fail,cc,seedMark,cseedmark);
                        full_sketch += tmp.size();
                        if (tmp.size() == 0){
                                break;
                        }
            //cout << cc << endl;
            
            if ((cc & 2) == 0){
                hg.addHyperedge(tmp);
            } else {
                zcover++;
            }
            if (cc & 1){
                degree++;
            }
                        
                        countTotal += n_fail+1;
                        countSuccess += 1;

            // Add hyperedge into the graph since it's not covered

            eager_thres = (z-zcover)*1.0/numSeeds;
            while (hg.getMaxDegree() >= eager_thres && selectedSeed.size() < numSeeds){
                        int newSeed = hg.getMaxDegreeNode();
                cout << "Select " << newSeed << endl;
                        zcover += hg.getMaxDegree();
                            eager_thres = (z-zcover)*1.0/numSeeds;

                        selectedSeed.push_back(newSeed);
                        seedMark[newSeed] = true;
                hg.selectNode(newSeed);
                }

            }

                cout << "Print degree and total samples: " << zcover << " " << degree << " " << countSuccess << " " << countTotal << " " << deg << " " << numTotalSamples << endl;
                max_sketch = hg.getSketchSize();
                cout << "Sketch reduce:" << max_sketch << " " << full_sketch << " " << double(max_sketch)/full_sketch << endl;
                if (double(zmax_sketch)/zfull_sketch < double(max_sketch)/full_sketch)
        {
            zmax_sketch = max_sketch;
            zfull_sketch = full_sketch;
        }
        double u = log(numNodes*(log(Tmax)/log(1+c)*log(zmax/zmin)));
        double lbound = lowerBound(u,countTotal,countTotal,degree*1.0/countTotal);
        cout << "Lower bound: " << lbound << endl;

        cout << "Assessment bound: " << lbound/ubound << " " << deg*2.0/z << endl;
        if (lbound/ubound >= 1-1/exp(1)-epsilon){
            cout << "Estimated value: " << (degree*normalized/countTotal) << endl;
                        numSuccessSamples = countSuccess;
                        numTotalSamples = countTotal;
                        deg = degree;
                        return true;
        }

/*                if (degree >= lambda1){
                        double epsilon_1 = (deg*1.0/numTotalSamples)/(degree*1.0/countTotal) - 1;
                        cout << "Epsilon_1 = " << epsilon_1 << " " << deg << " " << numTotalSamples << " " << degree << " " << countTotal << endl;
                        double epsilon_2 = sqrt((1+epsilon0)*log(3*log2(zmax/zmin)*log(Tmax)/log(1+c))/(degree));
                        cout << "Epsilon_2 = " << epsilon_2 << endl;
                        double epsilon_3 = sqrt((1+epsilon0)*(1-1/exp(1)-epsilon)*countTotal*log(3*log2(zmax/zmin)*Tmax)/((1+epsilon/3)*degree*numTotalSamples));
                        cout << "Epsilon_3 = " << epsilon_3 << endl;
                        cout << "Epsilon_t = " << (epsilon_1 + epsilon_2 + epsilon_1*epsilon_2)*(1-1/exp(1)-epsilon) + epsilon_3*(1-1/exp(1)) << endl;

                        if ((epsilon_1 + epsilon_2 + epsilon_1*epsilon_2)*(1-1/exp(1)-epsilon) + epsilon_3*(1-1/exp(1)) <= epsilon){
                                cout << "Estimated value: " << (degree*normalized/countTotal) << endl;
                                numSuccessSamples = countSuccess;
                                numTotalSamples = countTotal;
                                deg = degree;
                                return true;
                        }
                }*/
                while (checkpoint <= countTotal)
                        checkpoint = checkpoint*(1+c);
        }

        seeds = selectedSeed;
        deg = 0;

        if (selectedSeed.size() < k){
                vector<double> degreek(k+1,0);
                hg.buildSeedSet(seeds,k-selectedSeed.size(),degreek);
                deg += zcover + degreek[k-selectedSeed.size()];
        } else {
                if (zcover + hg.getSumK() < z){
                        vector<int> tmp;
                        while(zcover + hg.getSumK() < z){
                                tmp = sampler(n_fail,cc,seedMark,cseedmark);
                                full_sketch += tmp.size();
                                if (tmp.size() == 0){
                                        break;
                                }

                if (cc <= 1){
                                    hg.addHyperedge(tmp);
                            } else {
                                    zcover++;
                            }

                                countSuccess += 1;
                                countTotal += n_fail+1;
                        }
                }
                deg += zcover;
        }
        numSuccessSamples = countSuccess;
        numTotalSamples = countTotal;
        cout << "Print degree and total samples: " << zcover << " " << degree << " " << countSuccess << " " << countTotal << " " << deg << " " << numTotalSamples << endl;
        max_sketch = hg.getSketchSize();
        cout << "Sketch reduce:" << max_sketch << " " << full_sketch << " " << double(max_sketch)/full_sketch << endl;
        if (double(zmax_sketch)/zfull_sketch < double(max_sketch)/full_sketch)
        {
            zmax_sketch = max_sketch;
            zfull_sketch = full_sketch;
        }
        return false;
}

void DTA(std::function< std::vector<int32_t>(int32_t &,int &, std::vector<bool> &, std::vector<bool> &) > sampler, std::vector<int> & seeds, int t, float epsilon, float delta, double normalized)
{
    double deg = zmin;
        long long int numSuccessSamples = 1;
        long long int numTotalSamples = 1;
        long long int samples = 0;
        double c = 0.1;
        clock_t start = clock();
        cout << "Try z = " << zmin << endl;
        vector<bool> cseedmark(numNodes,false);
        vector<int> seed_t;
    double ubound = numNodes;
    double u = log(numNodes*(log(Tmax)/log(1+c)*log(zmax/zmin)));

        MSz(sampler,seeds,t,deg,epsilon,delta,numSuccessSamples,numTotalSamples,normalized,c,zmin,cseedmark, ubound);
        int numIter = ceil(log2(zmax/zmin));
    double zc = zmin;
        for (int i = 1; i < numIter; ++i){
        cout << "DEG and Z: " << deg << " " << z << endl;
                samples += numTotalSamples;
                reinitialize();
                for (int j = 0; j < seed_t.size(); ++j){
                        cseedmark[seed_t[j]] = false;
                }
                for (int j = 0; j < seeds.size(); ++j){
                        cseedmark[seeds[j]] = true;
                }
                seed_t = seeds;
        
        double N = ceil(pow(1+c,ceil(log(numTotalSamples)/log(1+c))));
        ubound = upperBound(u,numTotalSamples,N,zc/numTotalSamples);
                // cout << "Sketch reduce:" << zmax_sketch << " " << zfull_sketch << " " << double(zmax_sketch)/zfull_sketch << endl;

        cout << "Upper bound: " << ubound << endl;
        
        zc *= 2;
                cout << "Try z = " << zc << endl;

                if (MSz(sampler,seeds,t,deg,epsilon,delta,numSuccessSamples,numTotalSamples,normalized,c,zc,cseedmark,ubound)){
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
}

void SRA(std::function< std::vector<int32_t>(int32_t &,int &, std::vector<bool> &, std::vector<bool> &) > sampler, std::vector<int> & seeds, int t, float epsilon, float delta, double normalized)
{
    double deg = 1;
        long long int numSuccessSamples = 1;
        long long int numTotalSamples = 1;
        long long int samples = 0;
        clock_t start = clock();
        vector<bool> cseedmark(numNodes,false);
        double c = 0.1;
    double ubound = numNodes;

        MSz(sampler,seeds,t,deg,epsilon,delta,numSuccessSamples,numTotalSamples,normalized,c,zmax,cseedmark,ubound);
        cout << "Seed Nodes: ";
        for (unsigned int i = 0; i < seeds.size(); ++i){
                cout << seeds[i] << " ";
        }
        cout << endl;
        printf("Influence: %0.2lf\n",deg*normalized/numTotalSamples);
        cout << "Time: " << (float)(clock()-start)/CLOCKS_PER_SEC << endl;
        cout << "Memory: " << getMemValue()/1024.0 << endl;
        cout << "Samples: " << numTotalSamples << endl;
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

        char * alg = op.getPara("-alg");
        if (alg == NULL){
                cout << "Select an algorithm and run again!" << endl;
                return 1;
        }

        bool do_check = false;
        char * checked = op.getPara("-check");
        if (checked != NULL) do_check = true;

        vector<int> seeds(k,0);

        graph g;
        g.read_graph(inFile,false);
        Sampler sam(g);
        init_seeds(sam);
        int n = g.get_n_vertices();
        if (delta == 0) delta = 1.0/n;
        initialize(n,k,epsilon,delta,false,n*1.0);
        hg = HyperGraph(n,do_check,k,ceil(zmax));
        if (strcmp(alg, "SRA") == 0){
                SRA(std::bind(&Sampler::cc_polling, &sam, _1, _2, _3, _4), seeds, t, epsilon,delta,n*1.0);
        } else if (strcmp(alg, "DTA") == 0){
                DTA(std::bind(&Sampler::cc_polling, &sam, _1, _2, _3, _4),seeds,t,epsilon,delta,n*1.0);
        }
}
