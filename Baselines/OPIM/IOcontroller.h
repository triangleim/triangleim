#pragma once

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#endif

class IOcontroller
{
public:
	static void mkdir_absence(const char* outFolder)
	{
#if defined(_WIN32)
		CreateDirectoryA(outFolder, nullptr); // can be used on Windows
#else
		mkdir(outFolder, 0733); // can be used on non-Windows
#endif
	}

	/// Save a serialized file
	template <class T>
	static void save_file(const std::string filename, const T& output)
	{
		std::ofstream outfile(filename, std::ios::binary);
		if (!outfile.eof() && !outfile.fail())
		{
			StreamType res;
			serialize(output, res);
			outfile.write(reinterpret_cast<char*>(&res[0]), res.size());
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

	static bool cmprforedge(Edge& e1,Edge& e2){
		return e1.first<e2.first;
	}

	static void load_file(const std::string filename, Graph &input)
	{
		FILE *pFile;
		pFile = fopen(filename.c_str(), "rb");
		int numNodes;
		long long numEdges;
		fread(&(numNodes), sizeof(int), 1, pFile);	   // Nodes
		fread(&(numEdges), sizeof(long long), 1, pFile); // Edges
		Nodelist node_deg = Nodelist(numNodes + 1);
		fread(&(node_deg[1]), sizeof(int), numNodes, pFile); // in-degrees


		input = Graph(numNodes);
		// In-Neighbors
		for (unsigned int i = 1; i <= numNodes; ++i)
		{
			Nodelist tmp(node_deg[i]);
			fread(&tmp[0], sizeof(int), node_deg[i], pFile);
			for (auto nb : tmp)
			{
				input[i - 1].push_back(Edge(nb - 1, 0.0));
			}
		}
		// Edge weights (to Longlong)
		for (unsigned int i = 1; i <= numNodes; ++i)
		{
			std::vector<float> tmp(node_deg[i] + 1, 0);
			fread(&tmp[1], sizeof(float), node_deg[i], pFile);

			for (int j = 1; j < node_deg[i] + 1; ++j)
			{
				input[i - 1][j - 1].second = tmp[j];
			}
			std::sort(input[i - 1].begin(),input[i - 1].end(),cmprforedge);
		}
		
		
	}

	// /// Load a serialized file
	// template <class T>
	// static void load_file(const std::string filename, T& input)
	// {
	// 	std::ifstream infile(filename, std::ios::binary);
	// 	if (!infile.eof() && !infile.fail())
	// 	{
	// 		infile.seekg(0, std::ios_base::end);
	// 		const std::streampos fileSize = infile.tellg();
	// 		infile.seekg(0, std::ios_base::beg);
	// 		std::vector<uint8_t> res(fileSize);
	// 		infile.read(reinterpret_cast<char*>(&res[0]), fileSize);
	// 		infile.close();
	// 		input.clear();
	// 		auto it = res.cbegin();
	// 		input = deserialize<T>(it, res.cend());
	// 		res.clear();
	// 	}
	// 	else
	// 	{
	// 		std::cout << "Cannot open file: " + filename << '\n';
	// 		exit(1);
	// 	}
	// }

	/// Save graph structure to a file
	static void save_graph_struct(const std::string graphName, const Graph& vecGraph, const bool isReverse)
	{
		std::string postfix = ".vec.graph";
		if (isReverse) postfix = ".vec.rvs.graph";
		const std::string filename = graphName + postfix;
		save_file(filename, vecGraph);
	}

	/// Load graph structure from a file
	static void load_graph_struct(const std::string graphName, Graph& vecGraph, const bool isReverse)
	{
		// std::string postfix = ".vec.graph";
		// if (isReverse) postfix = ".vec.rvs.graph";
		const std::string filename = graphName;
		load_file(filename, vecGraph);
	}

	/// Get out-file name
	static std::string get_out_file_name(const std::string graphName, const std::string algName, const int seedsize,
		const std::string probDist, const float probEdge)
	{
		if (probDist == "UNI")
		{
			return graphName + "_" + algName + "_k" + std::to_string(seedsize) + "_" + probDist + std::
				to_string(probEdge);
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
	static void write_result(const std::string& outFileName, const TResult& resultObj, const std::string& outFolder, Graph& graph)
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
		std::cout << "  |Influence: " << influence <<", "<< influence/graph.size()<< std::endl;
		std::cout << "  |Self-estimated influence: " << influenceOriginal << ", "<< influenceOriginal/graph.size() << std::endl;
		std::cout << "  |#Seeds: " << seedSize << std::endl;
		std::cout << "  |#RR sets: " << RRsetsSize << std::endl;
		std::cout << "   --------------------" << std::endl;
		mkdir_absence(outFolder.c_str());
		std::ofstream outFileNew(outFolder + "/" + outFileName);
		if (outFileNew.is_open())
		{
			outFileNew << "Approx.: " << approx << std::endl;
			outFileNew << "Time (sec): " << runTime << std::endl;
			outFileNew << "Influence: " << influence << ", "<< influence/graph.size()<<std::endl;
			outFileNew << "Self-estimated influence: " << influenceOriginal <<", "<< influenceOriginal/graph.size() << std::endl;
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
