/*
 * @
 * @Date:
 * @LastEditors:  
 * @LastEditTime:
 * @Description: file content
 */
#ifndef _RWGRAPH_H
#define _RWGRAPH_H
#define WEIGHTEDCASCADE
#include <vector>
#include <random>
#include <unordered_map>
#include <set>
#include <algorithm>
#include "sfmt/SFMT.h"
#include "lib/segvcatch.h"
#include <string>
#include"hash.hpp"
typedef uint32_t UI;
typedef uint64_t ULL;
class BloomFilter
{
private:
	std::vector<uint> data;
public:
	BloomFilter();
	BloomFilter(std::vector<uint> d);
	BloomFilter operator|(const BloomFilter& b);
	BloomFilter& operator|=(const BloomFilter& b);
	BloomFilter operator&(const BloomFilter& b);
	BloomFilter& operator&=(const BloomFilter& b);
	bool operator==(const BloomFilter& b);
	bool isZero();
	void setData(std::vector<uint>& d);
	void genHash(std::vector<int>& d,std::vector<caluhash::HashFunction>& hs,size_t size,std::vector<int>& nodescc);
	void add(uint i);
	std::vector<uint>& getData();
	size_t getSize();
	const uint operator[](int u) const;
	const uint operator[](int u) ;
} ;
// Graph class defines fundamental operators on graph
class Graph
{
	friend class HyperGraph;

private:
	UI UI_MAX = 4294967295U;
	ULL ULL_MAX = 18446744073709551615ULL;

	// number of nodes
	unsigned int numNodes;
	// number of edges
	unsigned int numEdges;
	unsigned int numScc;
	unsigned int numFun;
	unsigned int indexWidth;
	// adjacency list
	std::vector<std::vector<int>> adjList; // in neighbors
	std::vector<std::vector<int>> edgeList;
	std::vector<int> node_deg;
	std::vector<uint> seeds;
	std::vector<int> node_scc;
	std::vector<std::vector<int>> adjOutList;
	std::vector<int> node_deg_out;
	std::vector<std::vector<UI>> weights;
	std::vector<ULL> weight_sampling;
	std::vector<BloomFilter> index;
	std::vector<caluhash::HashFunction> hs;
public:
		std::vector<uint> views;
		std::vector<uint> lifetimes;
		std::vector<uint> deads;
			std::vector<uint> parts;
	std::vector<uint> succs;
	std::vector<uint> initsuccs;
	std::vector<std::vector<int>> orders;
	std::vector<std::vector<int>> invalidorders;
	std::vector<double> rates, initsuccrates;
	unsigned long long numTriangles;
	// get a vector of neighbours of node u
	const std::vector<int> &operator[](int u) const;
	// return weights of neighbours of node u
	const std::vector<UI> &getWeight(int u) const;
	void readIndex(const char *filename);
	void readIndex2(const char *filename);
	// get a vector of neighbours of node u
	const std::vector<int> &operator[](int u);
	// return weights of neighbours of node u
	const std::vector<UI> &getWeight(int u);

	// get degree of node u
	int getDegree(int u) const;
	// get size of the graph
	int getSize() const;
	// get tri-size of the graph
	int getTriNumber() const;
	// get number of edges
	int getEdge() const;
	// get max weight
	ULL getMaxWeight() const;
	// read graph from a file
	void readGraphLT(const char *filename);
	// read graph from a file
	void readGraphIC(const char *filename);
	// write the graph to file
	void writeToFile(const char *filename);
};

class HyperGraph
{
private:
	// store the edges that a node is incident to
	std::vector<std::vector<int>> node_edge;
	// std::vector<std::vector<int>> node_edge_int;
	std::vector<std::vector<std::vector<int>>> node_edges;
	// store hyperedges
	std::vector<std::vector<std::vector<int>>> edge_node;
	unsigned int curEdge;
	unsigned int maxDegree;
	unsigned int numNodes;
	unsigned long long nemptyset;
	sfmt_t sfmtSeed;
	inline int randIndex_bin(const std::vector<UI> &w, unsigned int si);
	inline int randIndex_lin(const std::vector<UI> &w, unsigned int si);
	inline int randIndex_dart(const std::vector<UI> &w, unsigned int si);

public:
	HyperGraph(unsigned int n);
	std::vector<std::vector<int>> node_edge_int;
	ULL success;
	ULL allexpand;
	void updateDeg();
	void updateEdge();
	void addEdge(std::vector<std::vector<int>> &edge);
	void addEdgeD(std::vector<std::vector<int>> &edge);
	double emptyrate();
	int getMaxDegree();
	const std::vector<std::vector<int>> &getEdge(int e) const;
	const std::vector<std::vector<int>> &getEdge(int e);
	const std::vector<int> &getNode(int n) const;
	const std::vector<int> &getNodeIntersect(int n) const;
	const std::vector<int> &getNodeIntersect(int n) ;
	const std::vector<int> &getNode1(int n) const;
	const std::vector<int> &getNode1(int n) ;
	const std::vector<int> &getNode2(int n) const;
	const std::vector<int> &getNode2(int n) ;
	const std::vector<int> &getNode3(int n) const;
	const std::vector<int> &getNode3(int n) ;
	const std::vector<int> &getNode(int n);
	// get number of edges
	int getNumEdge() const;
	void clearEdges();
	void pollingLT1(Graph &g, std::vector<bool> &visit, std::vector<int> &mark_visit);
	bool pollingLT2(Graph &g, std::vector<unsigned int> &link, unsigned int k, std::vector<bool> &visit, std::vector<int> &mark_visit);
	bool pollingLT(Graph &g, std::vector<unsigned int> &link, unsigned int k, std::vector<bool> &visit, std::vector<int> &mark_visit);
	void pollingIC1(Graph &g, std::vector<bool> &visit, std::vector<int> &visit_mark);
	bool pollingIC2(Graph &g, std::vector<unsigned int> &link, unsigned int k, std::vector<bool> &visit, std::vector<int> &visit_mark);
	bool pollingIC(Graph &g, std::vector<unsigned int> &link, unsigned int k, std::unordered_map<int,bool>&visit, std::vector<int> &visit_mark);
	bool pollingIC3(Graph &g, std::vector<unsigned int> &link, unsigned int k, std::unordered_map<int,bool>&visit, std::vector<int> &visit_mark,int& node1,int& node2,int& node3);
	bool pollingLT3(Graph &g, std::vector<unsigned int> &link, unsigned int k, std::vector<bool> &visit, std::vector<int> &mark_visit,int& node1,int& node2,int& node3);
	
	
};

float getCurrentMemoryUsage();

#endif
