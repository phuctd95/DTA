PARA = -std=c++11 -Wall -O3 -g -pg
cplex = -m64 -O -fPIC -fexceptions -DNDEBUG \
-DIL_STD -DILOSTRICTPOD -I/opt/ibm/ILOG/CPLEX_Studio1261/cplex/include \
-I/opt/ibm/ILOG/CPLEX_Studio1261/concert/include \
-L/opt/ibm/ILOG/CPLEX_Studio1261/cplex/lib/x86-64_linux/static_pic -lilocplex \
-lcplex -L/opt/ibm/ILOG/CPLEX_Studio1261/concert/lib/x86-64_linux/static_pic 
com: rwgraph.cpp im_alg.cpp ds_alg.cpp ccm_alg.cpp el2bin.cpp option.cpp sampler.cpp xorshift.cpp graph.cpp IMSampler.cpp StepwiseHeap.cpp
	g++ -c mappedHeap.hpp -o mappedHeap.o $(PARA)
	g++ -c HeapData.hpp -o HeadData.o $(PARA)
	g++ -c StepwiseHeap.cpp -o StepwiseHeap.o $(PARA)
	g++ -c option.cpp -o option.o $(PARA)
	g++ -c rwgraph.cpp -o rwgraph.o -fopenmp $(PARA) $(cplex) -lconcert -lm -lpthread
	g++ -c xorshift.cpp -o xorshift.o -fopenmp ${PARA}
	g++ -c graph.cpp -o graph.o -fopenmp ${PARA}
	g++ -c sampler.cpp -o sampler.o -fopenmp $(PARA)
	g++ -c IMSampler.cpp -o IMSampler.o -fopenmp ${PARA}
	g++ ccm_alg.cpp rwgraph.o option.o sampler_f.o xorshift.o graph.o IMSampler.o StepwiseHeap.o -o ccm_alg -fopenmp $(PARA) sfmt/SFMT.c $(cplex) -lconcert -lm -lpthread
	g++ ds_alg.cpp rwgraph.o option.o sampler_f.o xorshift.o graph.o IMSampler.o StepwiseHeap.o -o ds_alg -fopenmp $(PARA) sfmt/SFMT.c $(cplex) -lconcert -lm -lpthread
	g++ im_alg.cpp rwgraph.o option.o sampler_f.o xorshift.o graph.o IMSampler.o StepwiseHeap.o -o im_alg -fopenmp $(PARA) sfmt/SFMT.c $(cplex) -lconcert -lm -lpthread
