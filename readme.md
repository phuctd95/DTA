Complie using make file. You need to replace the cplex directory.

Run:
	<algorithm> -i <input-file> -k <size> -epsilon <epsilon> -delta <delta> -m <edge-modle> -alg DTA -u <function>
	algorithm: We implement 3 algorithms for 3 different problems:
		Dominating Set: ds_alg.cpp
		Influence Maximiation: im_alg.cpp
		Lanmark Selection: ccm_alg.cpp
	input-file: the input graph.
	size: the size of the seed set.
	epsilon: the error bound.
	delta: the probability bound.
	edge-modle: the edge model in influence maximiation problem; we support 2 models: IC and LT.
	function: the upper-bound function; we provide 3 upper-bound functions: '0' for requirement function, '1' for top-k and '2' for dual.


Example 
	./im_alg -i epinions.txt.infmax.e1.bin -k 10 -epsilon 0.01 -delta 0 -m IC -alg DTA -u 0
