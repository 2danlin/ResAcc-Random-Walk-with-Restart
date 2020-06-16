Readme (ResAcc)
=========================
This package contains the source codes for our proposed method ResAcc.
Please cite our paper below for usage:
------------------------
@inproceedings{lin2020index,
	title={Index-Free Approach with Theoretical Guarantee for Efficient  Random Walk with Restart Query},
	author={Lin, Dandan and Wong, Raymond Chi-Wing and Xie, Min and Wei, Victor Junqiu},
	booktitle={2020 IEEE 36th International Conference on Data Engineering (ICDE)},
	year={2020},
	organization={IEEE}
}
------------------------


Preparation before Usage (taking DBLP for example)
==========================
a. Dataset preparation
   1. Make a folder under "dataset": 
	mkdir dblp
   2. Put the file of dataset under the new folder with the name "dblp.txt" (Format of the file could be found in Appendix)

b. Generate the file of query nodes (the parameters of the command could be found in Appendix)
   	./ResAcc -d dblp -algo GEN_QUERY -n 1000

c. Generate the ground truth for each query node
   1. Make a folder under "real_rwr":
	mkdir dblp
   2. Generate the ground truth for this dataset (the parameters of the command could be found in Appendix)
	./ResAcc -d dblp -algo GEN_GROUND_TRUTH -n 20

d. Results preparation
   1. Make a folder under "estimated_rwr":
	mkdir dblp

Usage Step
==========
a. Compilation
	./compile.sh

b. Execution (the parameters of the command could be found in Appendix)
	./ResAcc -d dblp -algo RESACC -n 20 -h 2 -r 0.1 -rf 0.1

c. Output Files (Format of those files can be found in Appendix)
   1. "resacc.time" which captures the time statistics of the program
   2. "source.txt" which stores the estimated RWR score of each node w.r.t the "source" node
	where "source" is the node-id of a query node.


Appendix A. Parameters of this program
------------------------------------
This program contains 8 parameters.
The general command of this program is shown below:
	./ResAcc [parameters]

1.-d <dataset> : the name of dataset 
	 e.g., "-d dblp"

2.-algo <algorithm> : the name of algorithm (all the supported algorithms could be found in Appendix)
	 e.g., "-algo RESACC"

3.-n <node_count> : no. of source nodes
	 e.g., "-n 20" (the default value is 20)

4.-alpha <alpha> : the value of the restart probability
	 e.g., "-alpha 0.2" (the default value is 0.2)

5.-eps <epsilon> : the value of relative error to be guaranteed
	 e.g., "-epsilon 0.5" (the default value is 0.5)

6.-h <no. of hops> : the number of hops used in the h-HopFWD phase
	 e.g., "-h 2" (the default value is 2) 

7.-r <rmax_hop> : the residue threshold used in the h-HopFWD phase
	 e.g., "-r 0.1" (the default value is 0.1)
Note that the true residue threshold is rmax_hop multiplied by 10^{-8}

8.-rf <rmax_f> : the residue threshold used in the OMFWD phase
	 e.g., "-rf 0.1" (the default value is 0.1)
Note that the true residue threshold is rmax_f multiplied by 10^{-8}
 
Appendix B. Algorithms supported by this program
------------------------------------
1. "GEN_QUERY" : the algorithm for generating the query nodes

2. "GEN_GROUND_TRUTH" : the algorithm for generating the exact RWR scores for the query nodes

3. "RESACC" : the algorithm for estimating the RWR scores for the query nodes (proposed by us)


Appendix C. Format of dataset file (taking "dblp.txt" for example)
------------------------------------
1. The first line in the file is 
 	no. of nodes in this graph

2. The following lines after the first line is
 	the edges of this graph which are represented by the following format:
   <From-node> <To-node> (separated by space)


Appendix D. Format of the result file for SSRWR query 
------------------------------------------
Each line corresponds to a point and its RWR score w.r.t the source node.
The format of each line is 
   <node id> <rwr score> 
   ...

Appendix E. Format of "time.txt"
------------------------------------------
Each line corresponds to the query time cost
1. average time for SSRWR queries
2. average time for running the h-HopFWD phase
3. average time for running the OMFWD phase
4. average time for running the Remedy phase

