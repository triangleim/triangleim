#ifndef _HYPERGRAPH_H_
#define _HYPERGRAPH_H_

#include "rwgraph.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstring>
#include <random>
#include "mappedHeap.hpp"
#include "HeapData.hpp"

#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_set_num_threads(int t) { return 1;}
inline omp_int_t omp_get_thread_num() { return 0;}
#endif

using namespace std;

/*
* building the hypergraph procedure which generates hyperedges following LT model
*/
long long addHyperedge(Graph & g, HyperGraph & hg, int t, long long num, bool lt)
{
	int numNodes = g.getSize();

	omp_set_num_threads(t);
	
	long long iter = 0;
	int c = 100;

        #pragma omp parallel
	{
		vector<int> visit_mark(numNodes+1,0);
		vector<bool> visit(numNodes+1,false);
		vector<unsigned int> link;       		
		if (lt == 0){
			while (iter < num){
       	                	for (int i = 0; i < c; ++i){
					vector<int> he;
	                       	        hg.pollingLT1(g,visit,visit_mark);
        	               	}
				#pragma omp atomic
				iter += c;
			}
		} else {
			unordered_map<int,bool> visit;
			while (iter < num){
                                for (int i = 0; i < c; ++i){
                                        vector<int> he;
                                        hg.pollingIC1(g,visit,visit_mark);
										// cout<<i+iter<<" "<<flush;
                                }
								
                                #pragma omp atomic
                                iter += c;
                        }

		}
	}
	
	hg.updateDeg();
	return hg.getNumEdge();
}

/*
* find seed nodes procedure using greedy algorithm
*/
void buildSeedSet(HyperGraph & hg, vector<int> & seeds, unsigned int n, int k, vector<double> &degree)
{	
		long long i;
	unsigned int j, l, maxInd;
	vector<int> e, e1, e2, e3, nList1, nList2, nList3;

	vector<int> nodeDegree(n, 0);
	vector<vector<int>> processed = vector<vector<int>>(hg.node_edge_int);
	vector<int> indx(n, 0);
	for (j = 0; j < n; ++j)
	{
		indx[j] = j;
		nodeDegree[j] = hg.getNodeIntersect(j + 1).size(); // j+1出现在了多少RR集
	}

	InfCost<int> hd(&nodeDegree[0]);
	MappedHeap<InfCost<int>> heap(indx, hd);
	long long numEdge = hg.getNumEdge();

	// check if an edge is removed
	vector<bool> edgeMark(numEdge, false);
	vector<bool> edgeMark1(numEdge, false);
	vector<bool> edgeMark2(numEdge, false);
	vector<bool> edgeMark3(numEdge, false);
	vector<bool> nodeMark(n + 1, true);
	uint totalCover = 0;
	double totalCost = 0;
	
	i = 1;
	// building each seed at a time
	while (totalCost < k && !heap.empty())
	{
		maxInd = heap.pop() + 1;
		nodeMark[maxInd] = false;
		// cout<<nodeDegree[maxInd - 1];
		totalCost++;

		e = hg.getNode(maxInd);
		e1 = hg.getNode1(maxInd);
		e2 = hg.getNode2(maxInd);
		e3 = hg.getNode3(maxInd);
		degree[i] = degree[i - 1] + nodeDegree[maxInd - 1];
		set<int> influencedNodes;
		seeds.push_back(maxInd);

		for (j = 0; j < e.size(); ++j)
		{
			if (edgeMark1[e[j]] && edgeMark2[e[j]] && edgeMark3[e[j]])
			{
				continue;
			}
			nList1 = hg.getEdge(e[j])[0];
			for (auto node : nList1)
			{
				if (nodeMark[node])
					influencedNodes.insert(node);
			}
			nList2 = hg.getEdge(e[j])[1];
			for (auto node : nList2)
			{
				if (nodeMark[node])
					influencedNodes.insert(node);
			}
			nList3 = hg.getEdge(e[j])[2];
			for (auto node : nList3)
			{
				if (nodeMark[node])
					influencedNodes.insert(node);
			}
		}
		for (j = 0; j < e1.size(); ++j)
		{
			if (edgeMark1[e1[j]])
			{
				continue;
			}
			edgeMark1[e1[j]] = true;
		}
		for (j = 0; j < e2.size(); ++j)
		{
			if (edgeMark2[e2[j]])
			{
				continue;
			}
			edgeMark2[e2[j]] = true;
		}
		for (j = 0; j < e3.size(); ++j)
		{
			if (edgeMark3[e3[j]])
			{
				continue;
			}
			edgeMark3[e3[j]] = true;
		}
		// 
		// for (j = 0; j < e1.size(); ++j)
		// {
		// 	if (edgeMark1[e1[j]])
		// 	{
		// 		continue;
		// 	}
		// 	edgeMark1[e1[j]] = true;
		// 	nList2 = hg.getEdge(e1[j])[1];
		// 	for (auto node : nList2)
		// 	{
		// 		if (nodeMark[node])
		// 			influencedNodes.insert(node);
		// 	}
		// 	nList3 = hg.getEdge(e1[j])[2];
		// 	for (auto node : nList3)
		// 	{
		// 		if (nodeMark[node])
		// 			influencedNodes.insert(node);
		// 	}
		// }
		// for (j = 0; j < e2.size(); ++j)
		// {
		// 	if (edgeMark2[e2[j]])
		// 	{
		// 		continue;
		// 	}
		// 	edgeMark2[e2[j]] = true;
		// 	nList1 = hg.getEdge(e2[j])[0];
		// 	for (auto node : nList1)
		// 	{
		// 		if (nodeMark[node])
		// 			influencedNodes.insert(node);
		// 	}
		// 	nList3 = hg.getEdge(e2[j])[2];
		// 	for (auto node : nList3)
		// 	{
		// 		if (nodeMark[node])
		// 			influencedNodes.insert(node);
		// 	}
		// }
		// for (j = 0; j < e3.size(); ++j)
		// {
		// 	if (edgeMark3[e3[j]])
		// 	{
		// 		continue;
		// 	}
		// 	edgeMark3[e3[j]] = true;
		// 	nList1 = hg.getEdge(e3[j])[0];
		// 	for (auto node : nList1)
		// 	{
		// 		if (nodeMark[node])
		// 			influencedNodes.insert(node);
		// 	}
		// 	nList2 = hg.getEdge(e3[j])[1];
		// 	for (auto node : nList2)
		// 	{
		// 		if (nodeMark[node])
		// 			influencedNodes.insert(node);
		// 	}
		// }
		uint ni = 1;
		for (auto node : influencedNodes)
		{
			//cout << "\r" << i << "/" << k << "\t" << ni << "/" << influencedNodes.size();
			nodeDegree[node - 1] += degree[i - 1];
			for (auto edge : e)
			{
				if (!edgeMark[edge] && edgeMark1[edge] && edgeMark2[edge] && edgeMark3[edge] && (find(hg.getNode(node).begin(), hg.getNode(node).end(), edge) == hg.getNode(node).end()))
				{
					nodeDegree[node - 1]++;
					continue;
				}
				bool ni = !edgeMark[edge] &&
						  find(processed[node].begin(), processed[node].end(), edge) == processed[node].end() &&
						  (edgeMark1[edge] || (find(hg.getNode1(node).begin(), hg.getNode1(node).end(), edge) != hg.getNode1(node).end())) &&
						  (edgeMark2[edge] || (find(hg.getNode2(node).begin(), hg.getNode2(node).end(), edge) != hg.getNode2(node).end())) &&
						  (edgeMark3[edge] || (find(hg.getNode3(node).begin(), hg.getNode3(node).end(), edge) != hg.getNode3(node).end()));
				if (ni)
				{
					processed[node].push_back(edge);
					nodeDegree[node - 1]++;
				}
			}
			nodeDegree[node - 1] -= degree[i];
			heap.heapify(node - 1);
			ni++;
		}
		for (auto edge : e)
		{
			if (edgeMark1[edge] && edgeMark2[edge] && edgeMark3[edge] && !edgeMark[edge])
			{
				edgeMark[edge] = true;
				totalCover++;
			}
		}
		i++;
		// if (totalCost == k)
		// 	cout << "\n"
		// 		 << totalCover << endl;
	}

	vector<int>().swap(nodeDegree);
	vector<int>().swap(e);
	vector<int>().swap(indx);
	vector<bool>().swap(edgeMark1);
	vector<bool>().swap(edgeMark2);
	vector<bool>().swap(edgeMark3);
}

#endif
