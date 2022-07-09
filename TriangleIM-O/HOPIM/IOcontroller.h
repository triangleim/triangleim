#pragma once

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#include <stdio.h>
#endif

class IOcontroller
{
public:
	static void mkdir_absence(const char *outFolder)
	{
#if defined(_WIN32)
		CreateDirectoryA(outFolder, nullptr); // can be used on Windows
#else
		mkdir(outFolder, 0733); // can be used on non-Windows
#endif
	}

	/// Save a serialized file
	template <class T>
	static void save_file(const std::string filename, const T &output)
	{
		std::ofstream outfile(filename, std::ios::binary);
		if (!outfile.eof() && !outfile.fail())
		{
			StreamType res;
			serialize(output, res);
			outfile.write(reinterpret_cast<char *>(&res[0]), res.size());
			outfile.close();
			res.clear();
			std::cout << "Save file successfully: " << filename << '\n';
		}
		else
		{
			std::cout << "Save file failed: " + filename << '\n';
			exit(1);
		}
	}

	/// Load a serialized file
	template <class T>
	static void load_file(const std::string filename, T &input)
	{
		std::ifstream infile(filename, std::ios::binary);
		if (!infile.eof() && !infile.fail())
		{
			infile.seekg(0, std::ios_base::end);
			const std::streampos fileSize = infile.tellg();
			infile.seekg(0, std::ios_base::beg);
			std::vector<uint8_t> res(fileSize);
			infile.read(reinterpret_cast<char *>(&res[0]), fileSize);
			infile.close();
			input.clear();
			auto it = res.cbegin();
			input = deserialize<T>(it, res.cend());
			res.clear();
		}
		else
		{
			std::cout << "Cannot open file: " + filename << '\n';
			exit(1);
		}
	}

	static bool cmprforedge(Edge& e1,Edge& e2){
		return e1.first<e2.first;
	}

	static void load_file(const std::string filename, GraphInfo &input)
	{
		FILE *pFile;
		pFile = fopen(filename.c_str(), "rb");
		fread(&(input.numNodes), sizeof(int), 1, pFile);	   // Nodes
		fread(&(input.numEdges), sizeof(long long), 1, pFile); // Edges
		input.node_deg = Nodelist(input.numNodes + 1);
		fread(&(input.node_deg[1]), sizeof(int), input.numNodes, pFile); // in-degrees

		Nodelist d;
		input.spedgelist.push_back(d); // all edges

		input.graph = Graph(input.numNodes);
		input.revgraph = RevGraph(input.numNodes);
		// In-Neighbors
		for (unsigned int i = 1; i <= input.numNodes; ++i)
		{
			Nodelist tmp(input.node_deg[i]);
			fread(&tmp[0], sizeof(int), input.node_deg[i], pFile);
			for (auto nb : tmp)
			{
				input.graph[i - 1].push_back(Edge(nb - 1, 0.0));
			}
		}
		// Edge weights (to Longlong)
		for (unsigned int i = 1; i <= input.numNodes; ++i)
		{
			std::vector<float> tmp(input.node_deg[i] + 1, 0);
			fread(&tmp[1], sizeof(float), input.node_deg[i], pFile);

			for (int j = 1; j < input.node_deg[i] + 1; ++j)
			{
				input.graph[i - 1][j - 1].second = tmp[j];
			}
			std::sort(input.graph[i - 1].begin(),input.graph[i - 1].end(),cmprforedge);
		}
		



		input.node_deg_out = Nodelist(input.numNodes + 1);
		fread(&(input.node_deg_out[1]), sizeof(int), input.numNodes, pFile); // out-degrees
		// cout<<node_deg_out[1]<<" "<< node_deg_out[2]<<" "<< node_deg_out[3]<<" "<< node_deg_out[4]<<endl;
		//  Out-Neighbors
		for (unsigned int i = 1; i <= input.numNodes; ++i)
		{
			std::vector<int> tmp(input.node_deg_out[i]);
			fread(&tmp[0], sizeof(int), input.node_deg_out[i], pFile);
			for (auto nb : tmp)
			{
				input.revgraph[i - 1].push_back(nb-1);
			}
			std::sort(input.revgraph[i - 1].begin(),input.revgraph[i - 1].end());
			// for(auto t:tmp)
			//         cout<<t<<endl;
		}
		input.triweights = TriWeights(input.numEdges + 1);

		fread(&(input.triweights[0]), sizeof(ULL), input.numEdges + 1, pFile); // for weight sampling
		// for(auto t:weight_sampling)
		//         cout<<t<<endl;
		// all-edges
		for (unsigned int i = 1; i <= input.numEdges; ++i)
		{
			std::vector<uint32_t> tmp(3);
			fread(&tmp[0], sizeof(int), 3, pFile);
			// for(auto t:tmp)
			//         cout<<t<<" ";
			input.spedgelist.push_back(tmp);
			// cout<<endl;
		}
		input.numTriangles = 0;
	}
	/// Save graph structure to a file
	static void save_graph_struct(const std::string graphName, const Graph &vecGraph, const bool isReverse)
	{
		std::string postfix = ".vec.graph";
		if (isReverse)
			postfix = ".vec.rvs.graph";
		const std::string filename = graphName + postfix;
		save_file(filename, vecGraph);
	}

	/// Load graph structure from a file
	static void load_graph_struct(const std::string graphName, GraphInfo &graphinfo, const bool isReverse)
	{
		// std::string postfix = ".vec.graph";
		// if (isReverse) postfix = ".vec.rvs.graph";
		const std::string filename = graphName;
		load_file(filename, graphinfo);
	}

	/// Get out-file name
	static std::string get_out_file_name(const std::string graphName, const std::string algName, const int seedsize,
										 const std::string probDist, const float probEdge)
	{
		if (probDist == "UNI")
		{
			return graphName + "_" + algName + "_k" + std::to_string(seedsize) + "_" + probDist + std::to_string(probEdge);
		}
		return graphName + "_" + algName + "_k" + std::to_string(seedsize) + "_" + probDist;
	}
	static void append_result(const std::string &outFileName, const TResult &resultObj, const std::string &outFolder, double mem)
	{
		const auto runTime = resultObj.get_running_time();
		const auto RRsetsSize = resultObj.get_RRsets_size();

		mkdir_absence(outFolder.c_str());
		std::ofstream outFileNew(outFolder + "/" + outFileName,std::ios_base::app);
		if (outFileNew.is_open())
		{
			outFileNew << runTime <<"\t"<< RRsetsSize << "\t"<<mem<< std::endl;
			outFileNew.close();
		}
	}
	/// Print the results
	static void write_result(const std::string &outFileName, const TResult &resultObj, const std::string &outFolder, const GraphInfo& gi)
	{
		const auto approx = resultObj.get_approximation();
		const auto runTime = resultObj.get_running_time();
		const auto influence = resultObj.get_influence();
		const auto influenceOriginal = resultObj.get_influence_original();
		const auto seedSize = resultObj.get_seed_size();
		const auto RRsetsSize = resultObj.get_RRsets_size();

		std::cout << "   --------------------" << std::endl;
		std::cout << "  |Approx.: " << approx << std::endl;
		std::cout << "  |Time (sec): " << runTime << std::endl;
		std::cout << "  |Influence: " << influence <<", "<< influence/gi.numTriangles<< std::endl;
		std::cout << "  |Self-estimated influence: " << influenceOriginal << ", "<< influenceOriginal/gi.numTriangles << std::endl;
		std::cout << "  |#Seeds: " << seedSize << std::endl;
		std::cout << "  |#RR sets: " << RRsetsSize << std::endl;
		std::cout << "   --------------------" << std::endl;
		mkdir_absence(outFolder.c_str());
		std::ofstream outFileNew(outFolder + "/" + outFileName);
		if (outFileNew.is_open())
		{
			outFileNew << "Approx.: " << approx << std::endl;
			outFileNew << "Time (sec): " << runTime << std::endl;
			outFileNew << "Influence: " << influence << ", "<< influence/gi.numTriangles<<std::endl;
			outFileNew << "Self-estimated influence: " << influenceOriginal <<", "<< influenceOriginal/gi.numTriangles << std::endl;
			outFileNew << "#Seeds: " << seedSize << std::endl;
			outFileNew << "#RR sets: " << RRsetsSize << std::endl;
			outFileNew.close();
		}
	}

	/// Print the seeds
	static void write_order_seeds(const std::string &outFileName, const TResult &resultObj, const std::string &outFolder, const std::string &cascadeModel = "IC")
	{
		auto vecSeed = resultObj.get_seed_vec();
		mkdir_absence(outFolder.c_str());
		const auto outpath = outFolder + "/seed";
		mkdir_absence(outpath.c_str());
		std::ofstream outFile(outpath + "/seed_" + outFileName);
		for (auto i = 0; i < vecSeed.size(); i++)
		{
			outFile << vecSeed[i]+1 << '\n';
		}
		outFile.close();
	}
};

using TIO = IOcontroller;
using PIO = std::shared_ptr<IOcontroller>;
