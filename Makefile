PARA = -std=c++11 -Wall -O3 -g -pg
com: rwgraph.cpp im_alg.cpp bcm_alg.cpp ccm_alg.cpp el2bin.cpp option.cpp sampler.cpp xorshift.cpp graph.cpp IMSampler.cpp StepwiseHeap.cpp triangle_alg.cpp k_path_alg.cpp
	g++ -c mappedHeap.hpp -o mappedHeap.o $(PARA)
	g++ -c HeapData.hpp -o HeadData.o $(PARA)
	g++ -c StepwiseHeap.cpp -o StepwiseHeap.o $(PARA)
	g++ -c option.cpp -o option.o $(PARA)
	g++ -c rwgraph.cpp -o rwgraph.o -fopenmp $(PARA) -m64 -O -fPIC -fexceptions -DNDEBUG \
-DIL_STD -DILOSTRICTPOD -I/opt/ibm/ILOG/CPLEX_Studio1261/cplex/include \
-I/opt/ibm/ILOG/CPLEX_Studio1261/concert/include \
-L/opt/ibm/ILOG/CPLEX_Studio1261/cplex/lib/x86-64_linux/static_pic -lilocplex \
-lcplex -L/opt/ibm/ILOG/CPLEX_Studio1261/concert/lib/x86-64_linux/static_pic \
-lconcert -lm -lpthread
	g++ -c xorshift.cpp -o xorshift.o -fopenmp ${PARA}
	g++ -c graph.cpp -o graph.o -fopenmp ${PARA}
	g++ -c sampler_f.cpp -o sampler_f.o -fopenmp $(PARA)
	g++ -c IMSampler.cpp -o IMSampler.o -fopenmp ${PARA}
#	g++ ssaz.cpp rwgraph.o option.o sampler.o xorshift.o graph.o IMSampler.o -o ssaz -fopenmp $(PARA) sfmt/SFMT.c
	g++ feature.cpp rwgraph.o option.o sampler_f.o xorshift.o graph.o IMSampler.o StepwiseHeap.o -o feature -fopenmp $(PARA) sfmt/SFMT.c -m64 -O -fPIC -fexceptions -DNDEBUG \
-DIL_STD -DILOSTRICTPOD -I/opt/ibm/ILOG/CPLEX_Studio1261/cplex/include \
-I/opt/ibm/ILOG/CPLEX_Studio1261/concert/include \
-L/opt/ibm/ILOG/CPLEX_Studio1261/cplex/lib/x86-64_linux/static_pic -lilocplex \
-lcplex -L/opt/ibm/ILOG/CPLEX_Studio1261/concert/lib/x86-64_linux/static_pic \
-lconcert -lm -lpthread
	g++ feature_propability.cpp rwgraph.o option.o sampler_f.o xorshift.o graph.o IMSampler.o StepwiseHeap.o -o feature_propability -fopenmp $(PARA) sfmt/SFMT.c -m64 -O -fPIC -fexceptions -DNDEBUG \
-DIL_STD -DILOSTRICTPOD -I/opt/ibm/ILOG/CPLEX_Studio1261/cplex/include \
-I/opt/ibm/ILOG/CPLEX_Studio1261/concert/include \
-L/opt/ibm/ILOG/CPLEX_Studio1261/cplex/lib/x86-64_linux/static_pic -lilocplex \
-lcplex -L/opt/ibm/ILOG/CPLEX_Studio1261/concert/lib/x86-64_linux/static_pic \
-lconcert -lm -lpthread