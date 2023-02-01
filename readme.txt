This is the repository for the Triangular Stability Maximization and Triangle Influence Maximization Project.

Sandwich:
	Triangle-O: OPIM-C-based algorithms
		OPIM: the original OPIM code (Baseline)
		HOPIM: Homologous TriangleIM-O
		GOPIM: General TriangleIM-O
		UOPIM: Component TriangleIM-O
	Triangle-S: D-SSA-based algorithms (√)
		DSSA: the original D-SSA code
		HDSSA: Homologous TriangleIM-S
		GDSSA: General TriangleIM-S
		DSSA: Component TriangleIM-O
	SandwichApproximation: Select the best solution among upper, lower and general.
JointBaking:
	Baking: Get the solution of upper/lower.
	Heuristics: Get the heuristic solution of General TriangleIM.
	SandwichApproximation: Select the best solution among upper, lower and general.
Validation: Evaluate the results of the algorithms
How to use?
0.  Use make in the algorithm folder.
/* Generate usable dataset files (the same as SSA/D-SSA) */

1.  Computing edge weights (probabilities) as described in the experiments:
	Enter Sandwich/TriangleIM-S/HDSSA:
	./format <input file> <output file> 1

	<input file>: the path to text file in edge list format with no weights: the first line contains the number of nodes n and number of edges m, each of the next m lines describes an edge following the format: <src> <dest>. Node index starts from 1. All self-linked edges should be removed in advance.
	<output file>: the path to text output file with edge probabilities
	The last parameter (1) means the input graph is considered as directed. (0) means undirected.
2.  Conversion from a text format to binary file
			Enter Sandwich/TriangleIM-S/HDSSA:
        	./el2bin <input file> <output file>

    	<input file>: the path to text file in weighted edge list format: the first line contains the number of nodes n and number of edges m, each of the next m lines describes an edge following the format: <src> <dest> <weight>. Node index starts from 1. 
    	<output file>: the path to binary output file
/* Compute seed set */
3.  Run the corresponding executable file.
        ./DSSA -i <dataset binary filepath> -o <output seeds> -k <the size of seed set> -triangle <the total number of triangle>
                -m <IC or LT> -delta <delta,1/n> -epsilon <epsilon,0.1>
		./OPIM -func=1 -gname=<dataset binary file> -alg=opim-c -mode=2 -seedsize=<the size of seed set> -eps=<epsilon,0.1> 
				-delta=<delta,1/n> -model=<IC or LT> -pdist=load -triangle=<the total number of triangle>  (For GOPIM, mode should be 0(vanilla) to avoid too small number of samples generated.)
		./HeuIC(LT) -i <dataset binary filename> -o <output seeds> -k <the size of seed set>
		./Sandwich -i <dataset binary filepath> -o <original seeds filepath> -u <upper seeds filepath> -l <lower seeds filepath> -f <output final seeds filepath> -k <the size of seed set> 
				-triangle <the total number of triangle> -m <IC or LT> -delta <delta,1/n> -epsilon <epsilon,0.1> -gamma=<gamma,0.1>
		For OPIM, dataset files should store at "/datasets" folder, otherwise argument.h should be changed.
/* Evaluate the results */
4. Run the corresponding executable file.
        ./verifyInf -i <dataset binary file> -o <seeds> -k <the size of seed set> -times <the number of samples to be generated>
        (Case study) ./verifyDead -i <dataset binary file> -o <seeds> -k <the size of seed set> -times <the number of samples to be generated> -nodescc <feature file of twitch-games>

The ablation experiment of HTriangleD-SSA can be referred to  Validation/HSSA_ForVerify/verifyInf.cpp
HeuIC can be used for directed triangle counting.
The codes only apply to the case of nodes with equal incoming edge probabilities, otherwise the codes need to be adjusted.
The executable files in the repository may not correspond exactly to their source code. *You should ``make'' first and then execute them.*
