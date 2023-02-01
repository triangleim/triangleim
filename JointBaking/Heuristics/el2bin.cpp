/*
* Functionality: convert from a graph file in weighted edge list to a binary file
* Syntax:
	./el2bin <graph file input> <binary graph output>

* The graph file input must follow the following format:
	<number of nodes> <number of edges>
	<first node of edge 1> <second node of edge 1> <weight of edge 1>
	...
	<first node of the last edge> <second node of the last edge> <weight of the last edge>

* The binary graph output will be used in our SSA/DSSA algorithms for fast reading

* 
*/
#include<iostream>
#include <cstdio>
#include <fstream>
#include <random>
#include "GLib.hpp"
#include <cmath>
#include <cstring>

using namespace std;

vector<int> intersection(vector<int> &v1,
                                      vector<int> &v2){
    vector<int> v3;

    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());

    set_intersection(v1.begin(),v1.end(),
                          v2.begin(),v2.end(),
                          back_inserter(v3));
    return v3;
}

int main(int argc, char ** argv)
{
	ifstream in(argv[1]);
	srand(time(NULL));
	int n,u,v;
	long long m;
	float w;
	in >> n >> m;
	printf("%d %lld\n", n, m);
	vector<int> degree(n+1,0);
	vector<vector<int> > eList(n+1);
	vector<vector<int> > edgeList(m+1); // [[ID,start,end,//triangles]]
	vector<vector<int> > eOutList(n+1);
	vector<vector<float> > weight(n+1);
	vector<float> weightR(n+1,0);
	vector<unsigned long long> weightSample;
	printf("Reading the graph!\n");

	for (long long i = 0; i < m; ++i){
		if(i%100000==0)
			cout<<"\r"<<i<<"/"<<m<<"\t"<<i/(float)m;
		in >> u >> v >> w;
		degree[v]++; // in degree
		eList[v].push_back(u); // in neighbors
		eOutList[u].push_back(v); // out neighbors
		weight[v].push_back(w); // weight of edge
		weightR[u] += 1; // out degree
		edgeList[i+1].push_back(i+1);
		edgeList[i+1].push_back(u);
		edgeList[i+1].push_back(v);
	}
	
	in.close();
	vector<long long> tmpWeight(m+1,0);
	for (long long i=1;i<m+1;++i){
		int u = edgeList[i][1];
		int v = edgeList[i][2];
		int tri = 0;
		// out-triangles
		tri += intersection(eList[u],eList[v]).size();
		tri += intersection(eOutList[u],eOutList[v]).size();
		tri += intersection(eOutList[u],eList[v]).size();
		tri += intersection(eList[u],eOutList[v]).size();
		//edgeList[i].push_back(tri);
		tmpWeight[i]=tri;
	}
    // for(auto t:tmpWeight)
    //             cout<<t<<endl;
	partial_sum(tmpWeight.begin(), tmpWeight.end(), back_inserter(weightSample)); // partial sum
    // for(auto t:weightSample)
    //             cout<<t<<endl;
    vector<size_t> idx(n);

	FILE * pFile;
	pFile = fopen(argv[2],"wb");
	// Number of vertexs
	fwrite(&n, sizeof(int), 1, pFile);
	// Number of edges
	fwrite(&m, sizeof(long long), 1, pFile);
	// Assume that id is 1 to i orderly
    for (int i = 0; i < n; ++i){
		idx[i] = i;
	}
	vector<int> inv_idx(n);
	for (int i = 0; i < n; ++i){
		inv_idx[idx[i]]	= i;
	}
	// Write in degree
	vector<int> iTmp(n);
	
	for (int i = 0; i < n; ++i){
		iTmp[i] = degree[idx[i]+1];
	}
	
	// Write node in-degrees
	fwrite(&iTmp[0], sizeof(int), n, pFile);
	
	for (int i = 1; i <= n; ++i){
		// Write in-neighbors
		for (unsigned int j = 0; j < eList[idx[i-1]+1].size(); ++j){
			iTmp[j] = inv_idx[eList[idx[i-1]+1][j]-1]+1; // in neighbors
		}
		fwrite(&iTmp[0], sizeof(int), eList[idx[i-1]+1].size(), pFile);
	}

	for (int i = 1; i <= n; ++i){
		// Write weights
                fwrite(&weight[idx[i-1] + 1][0], sizeof(float), weight[idx[i-1]+1].size(), pFile);
        }
	// Write out degree
	// vector<int> iTmp(n);
	
	for (int i = 0; i < n; ++i){
		iTmp[i] = weightR[idx[i]+1];
	}
	
	// Write node out-degrees
	fwrite(&iTmp[0], sizeof(int), n, pFile);
	
	for (int i = 1; i <= n; ++i){
		// Write out-neighbors
		for (unsigned int j = 0; j < eOutList[idx[i-1]+1].size(); ++j){
			iTmp[j] = inv_idx[eOutList[idx[i-1]+1][j]-1]+1; // in neighbors
		}
		fwrite(&iTmp[0], sizeof(int), eOutList[idx[i-1]+1].size(), pFile);
	}
	// Write sample weight
	

	fwrite(&weightSample[0], sizeof(unsigned long long), m+1, pFile);
	// write full adj list
	for (int i = 1; i <= m; ++i){
		// Write full adj list
        fwrite(&edgeList[i][0], sizeof(int), edgeList[i].size(), pFile);
    }
	fclose(pFile);
	printf("Done!\n");
	return 1;
}
