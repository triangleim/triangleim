/**
 * @file OPIM.cpp
 * @brief This project implements the OPIM and OPIM-C for the following paper:
 * Jing Tang, Xueyan Tang, Xiaokui Xiao, Junsong Yuan, "Online Processing Algorithms for Influence Maximization," in Proc. ACM SIGMOD, 2018.
 *
 * @author Jing Tang (Nanyang Technological University)
 *
 * Copyright (C) 2018 Jing Tang and Nanyang Technological University. All rights reserved.
 *
 */

#include "stdafx.h"
#include "SFMT/dSFMT/dSFMT.c"
#include "alg.cpp"
#include <unistd.h>
/*
 * convert from an integer to a string
 */
std::string intToStr(int i)
{
	std::stringstream ss;
	ss << i;
	return ss.str();
}

/*
 * convert from a strong to an integer
 */
unsigned int strToInt(std::string s)
{
	unsigned int i;
	std::istringstream myStream(s);

	if (myStream >> i)
	{
		return i;
	}
	else
	{
		std::cout << "String " << s << " is not a number." << std::endl;
		return atoi(s.c_str());
	}
	return i;
}
float getCurrentMemoryUsage()
{

	std::string pid = intToStr(unsigned(getpid()));
	std::string outfile = "tmp_" + pid + ".txt";
	std::string command = "pmap " + pid + " | grep -i Total | awk '{print $2}' > " + outfile;
	system(command.c_str());

	std::string mem_str;
	std::ifstream ifs(outfile.c_str());
	std::getline(ifs, mem_str);
	ifs.close();

	mem_str = mem_str.substr(0, mem_str.size() - 1);
	float mem = (float)strToInt(mem_str);

	command = "rm " + outfile;
	system(command.c_str());

	return mem / 1024;

	return 0;
}
int main(int argc, char *argv[])
{
	// Randomize the seed for generating random numbers
	dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));
	const TArgument Arg(argc, argv);
	const std::string infilename = Arg._dir + "/" + Arg._graphname;
	if (Arg._func == 0 || Arg._func == 2)
	{
		// Format the graph
		GraphBase::format_graph(infilename, Arg._mode);
		if (Arg._func == 0)
			return 1;
	}

	std::cout << "---The Begin of " << Arg._outFileName << "---\n";
	Timer mainTimer("main");
	GraphInfo graphinfo;

	// Load the reverse graph
	GraphBase::load_graph(infilename, true, graphinfo, Arg._probDist, Arg._probEdge);
	graphinfo.numTriangles = Arg._triangle;

	if (Arg._model == LT)
	{
		// Normalize the propagation probabilities in accumulation format for LT cascade model for quickly generating RR sets
		to_normal_accum_prob(graphinfo.graph);
	}
	// Initialize a result object to record the results
	TResult tRes;
	TAlg tAlg(graphinfo, tRes);
	tAlg.set_cascade_model(Arg._model); // Set propagation model

	std::cout << "  ==>Graph loaded for RIS! total time used (sec): " << mainTimer.get_total_time() << '\n';
	int mode = 2; // Default is to use the minimum upper bound among all the rounds
	if (Arg._mode == "0" || Arg._mode == "vanilla")
	{
		mode = 0;
	}
	else if (Arg._mode == "1" || Arg._mode == "last")
	{
		mode = 1;
	}
	auto delta = Arg._delta;
	if (delta < 0)
		delta = 1.0 / graphinfo.numNodes;
	if (Arg._algName == "opim-c" || Arg._algName == "OPIM-C")
	{
		tAlg.opimc(Arg._seedsize, Arg._eps, delta, mode);
	}
	else if (Arg._algName == "opim" || Arg._algName == "OPIM")
	{
		tAlg.opim(Arg._seedsize, Arg._samplesize, delta, mode);
	}
	double mem = getCurrentMemoryUsage();
	std::cout << "Memory: " << mem << " MB" << std::endl;
	TIO::append_result("timesamplemem.tsv", tRes, Arg._resultFolder, mem);
	TIO::write_result(Arg._outFileName, tRes, Arg._resultFolder,graphinfo);
	TIO::write_order_seeds(Arg._outFileName, tRes, Arg._resultFolder, Arg._model == LT ? "LT" : "IC");
	std::cout << "---The End of " << Arg._outFileName << "---\n";
	return 0;
}
