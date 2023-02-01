#include "rwgraph.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <queue>
#include <cstdlib>
#include <chrono>
#include <unistd.h>
#include <sstream>
#include <random>
#include <cmath>
using namespace std;
// BloomFilter::BloomFilter()
// {
// }
// BloomFilter::BloomFilter(vector<uint> d)
// {
//         data = vector<uint>(d);
// }
// BloomFilter BloomFilter::operator|(const BloomFilter &b)
// {
//         BloomFilter result;
//         size_t s = getSize();
//         for (uint i = 0; i < s; i++)
//         {
//                 result.add(data[i] | b[i]);
//         }
//         return result;
// }
// BloomFilter &BloomFilter::operator|=(const BloomFilter &b)
// {
//         size_t s = getSize();
//         for (uint i = 0; i < s; i++)
//         {
//                 data[i] |= b[i];
//         }
//         return *this;
// }
// BloomFilter BloomFilter::operator&(const BloomFilter &b)
// {
//         BloomFilter result;
//         size_t s = getSize();
//         for (uint i = 0; i < s; i++)
//         {
//                 result.add(data[i] & b[i]);
//         }
//         return result;
// }
// BloomFilter &BloomFilter::operator&=(const BloomFilter &b)
// {
//         size_t s = getSize();
//         for (uint i = 0; i < s; i++)
//         {
//                 data[i] &= b[i];
//         }
//         return *this;
// }
// bool BloomFilter::operator==(const BloomFilter &b)
// {
//         size_t s = getSize();
//         for (uint i = 0; i < s; i++)
//         {
//                 if (data[i] != b[i])
//                         return false;
//         }
//         return true;
// }
// bool BloomFilter::isZero()
// {
//         size_t s = getSize();
//         for (uint i = 0; i < s; i++)
//         {
//                 if (data[i] != 0)
//                         return false;
//         }
//         return true;
// }
// void BloomFilter::setData(std::vector<uint> &d)
// {
//         data = vector<uint>(d);
// }
// void BloomFilter::genHash(std::vector<int> &d, std::vector<caluhash::HashFunction> &hs, size_t size, std::vector<int> &nodescc)
// {
//         data = vector<uint>(size, 0);
//         set<int> sccs;
//         for (auto i : d)
//                 sccs.insert(nodescc[i]);
//         for (int i : sccs)
//         {
//                 for (auto h : hs)
//                 {
//                         uint hashv = h(i) % (size * 32);
//                         uint k = hashv / 32;
//                         uint r = hashv % 32;
//                         uint one = ((uint)1) << r;
//                         data[k] |= one;
//                 }
//         }
// }
// void BloomFilter::add(uint i)
// {
//         data.push_back(i);
// }
// vector<uint> &BloomFilter::getData()
// {
//         return data;
// }
// size_t BloomFilter::getSize()
// {
//         return data.size();
// }
// const uint BloomFilter::operator[](int u) const
// {
//         return data[u];
// }

// const uint BloomFilter::operator[](int u)
// {
//         return data[u];
// }
const vector<int> &Graph::operator[](int u) const
{
        return adjList[u];
}

const vector<int> &Graph::operator[](int u)
{
        return adjList[u];
}

const vector<UI> &Graph::getWeight(int u) const
{
        return weights[u];
}

const vector<UI> &Graph::getWeight(int u)
{
        return weights[u];
}

/*
* get degree of node u
*/
int Graph::getDegree(int u) const
{
        return adjList[u].size();
}

/*
* get the number of nodes
*/
int Graph::getSize() const
{
        return numNodes;
}
/*
* get the number of triangles, uncompleted
*/
int Graph::getTriNumber() const
{
        if (numTriangles > 0)
                return numTriangles;
        return numNodes;
}
/*
* get the number of edges
*/
int Graph::getEdge() const
{
        return numEdges;
}
/*
* get max weight
*/
ULL Graph::getMaxWeight() const
{
        return weight_sampling.back();
}
/*
* read binary graph input for LT model
* difference between LT and IC is for LT we accumulate the weights for fast choosing a random node
*/
void Graph::readGraphLT(const char *filename)
{
        FILE *pFile;
        pFile = fopen(filename, "rb");
        fread(&numNodes, sizeof(int), 1, pFile);       // Nodes
        fread(&numEdges, sizeof(long long), 1, pFile); // Edges
        node_deg = vector<int>(numNodes + 1);
        fread(&node_deg[1], sizeof(int), numNodes, pFile); // in-degrees

        vector<int> a, c, d;
        vector<UI> b;
        adjList.push_back(a);    // in-neighbors
        adjOutList.push_back(c); // out-neighbors
        edgeList.push_back(d);   //all edges
        weights.push_back(b);    // edge weights
                                 // In-Neighbors
        for (unsigned int i = 1; i <= numNodes; ++i)
        {
                vector<int> tmp(node_deg[i]);
                fread(&tmp[0], sizeof(int), node_deg[i], pFile);
                std::sort(tmp.begin(),tmp.end());
                adjList.push_back(tmp);
        }
        // Edge weights (to Longlong)
        for (unsigned int i = 1; i <= numNodes; ++i)
        {
                vector<float> tmp(node_deg[i] + 1, 0);
                vector<UI> tmp1(node_deg[i] + 1, 0);
                fread(&tmp[1], sizeof(float), node_deg[i], pFile);

                for (int j = 1; j < node_deg[i] + 1; ++j)
                {
                        tmp[j] += tmp[j - 1];
                        if (tmp[j] >= 1)
                        {
                                tmp1[j] = UI_MAX;
                        }
                        else
                        {
                                tmp1[j] = tmp[j] * UI_MAX;
                        }
                }

                weights.push_back(tmp1);
                node_deg[i]++; // what?
        }
        node_deg_out = vector<int>(numNodes + 1);
        fread(&node_deg_out[1], sizeof(int), numNodes, pFile); // out-degrees
        //cout<<node_deg_out[1]<<" "<< node_deg_out[2]<<" "<< node_deg_out[3]<<" "<< node_deg_out[4]<<endl;
        // Out-Neighbors
        for (unsigned int i = 1; i <= numNodes; ++i)
        {
                vector<int> tmp(node_deg_out[i]);
                fread(&tmp[0], sizeof(int), node_deg_out[i], pFile);
                adjOutList.push_back(tmp);
                // for(auto t:tmp)
                //         cout<<t<<endl;
        }
        weight_sampling = vector<ULL>(numEdges + 1);

        fread(&weight_sampling[0], sizeof(ULL), numEdges + 1, pFile); // for weight sampling
        // for(auto t:weight_sampling)
        //         cout<<t<<endl;
        // all-edges
        for (unsigned int i = 1; i <= numEdges; ++i)
        {
                vector<int> tmp(3);
                fread(&tmp[0], sizeof(int), 3, pFile);
                // for(auto t:tmp)
                //         cout<<t<<" ";
                edgeList.push_back(tmp);
                // cout<<endl;
        }
        numTriangles = 0;
}
void split(const string &s, vector<string> &tokens, const string &delimiters = " ")
{
        string::size_type lastPos = s.find_first_not_of(delimiters, 0);
        string::size_type pos = s.find_first_of(delimiters, lastPos);
        while (string::npos != pos || string::npos != lastPos)
        {
                tokens.push_back(s.substr(lastPos, pos - lastPos)); //use emplace_back after C++11
                lastPos = s.find_first_not_of(delimiters, pos);
                pos = s.find_first_of(delimiters, lastPos);
        }
}

// void Graph::readIndex(const char *filename, const char *filenc)
// {
//         ifstream in(filename);
//         string line;
//         vector<string> array;
//         getline(in, line);
//         split(line, array, " ");
//         numScc = atoi(array[0].c_str());
//         numFun = 0;
//         for (int i = 1; i < array.size(); i++)
//         {
//                 seeds.push_back(unsigned(atoi(array[i].c_str())));
//                 numFun += 1;
//         }

//         index = vector<BloomFilter>(numScc + 1);

//         while (!in.eof())
//         {
//                 string line;
//                 vector<string> array;
//                 getline(in, line);
//                 split(line, array, " ");
//                 for (int i = 1; i < array.size(); i++)
//                         index[atoi(array[0].c_str())].add(unsigned(atoi(array[i].c_str())));
//         }
//         indexWidth = index[1].getSize();
//         for (int i = 0; i < numFun; i++)
//         {
//                 unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//                 seeds.push_back(seed);
//                 std::mt19937_64 tmpg = std::mt19937_64(seed);
//                 hs.push_back(caluhash::HashFunction(tmpg, int(ceil(log2(indexWidth * 32)))));
//         }
//         cout << index[1].getData()[0] << index[1].getData()[1] << index[1].getData()[2] << index[1].getData()[3] << endl;
//         in.close();
//         ifstream in2(filenc);
//         int node, scc;
//         node_scc = vector<int>(numNodes + 1);
//         for (uint i = 1; i < numNodes + 1; i++)
//         {
//                 in2 >> scc >> node;
//                 node_scc[node] = scc;
//         }
//         in2.close();
// }
void Graph::readIndex(const char *filename)
{
        ifstream in(filename);
        string line;
        vector<string> array;
        getline(in, line);
        views=vector<uint>(numNodes+1);
        deads=vector<uint>(numNodes+1);
        lifetimes=vector<uint>(numNodes+1);
        // split(line, array, ",");
        // numScc = atoi(array[0].c_str());
        // numFun = 0;
        // for (int i = 1; i < array.size(); i++)
        // {
        //         seeds.push_back(unsigned(atoi(array[i].c_str())));
        //         numFun += 1;
        // }

        // index = vector<BloomFilter>(numScc + 1);

        while (!in.eof())
        {
                string line;
                vector<string> array;
                getline(in, line);
                if(line=="")
                        break;
                split(line, array, ",");
                // cout<<"\r"<<atoi(array[5].c_str())+1<<flush;
                views[atoi(array[5].c_str())+1]=(unsigned(atoi(array[0].c_str())));
                lifetimes[atoi(array[5].c_str())+1]=(unsigned(atoi(array[2].c_str())));
                deads[atoi(array[5].c_str())+1]=(unsigned(atoi(array[6].c_str())));
        }
        in.close();
        cout<<"done!"<<endl<<flush;
}
/*
* read input graph for IC model
*/
void Graph::readGraphIC(const char *filename)
{
        FILE *pFile;
        pFile = fopen(filename, "rb");
        fread(&numNodes, sizeof(int), 1, pFile);       // Nodes
        fread(&numEdges, sizeof(long long), 1, pFile); //Edges
        node_deg = vector<int>(numNodes + 1);
        fread(&node_deg[1], sizeof(int), numNodes, pFile); // in-degrees

        vector<int> a, c, d;
        vector<UI> b;
        adjList.push_back(a);    // in-neighbors
        adjOutList.push_back(c); // out-neighbors
        edgeList.push_back(d);   //all edges
        weights.push_back(b);    // edge weights

        // In-Neighbors
        for (unsigned int i = 1; i <= numNodes; ++i)
        {
                vector<int> tmp(node_deg[i]);
                fread(&tmp[0], sizeof(int), node_deg[i], pFile);
                std::sort(tmp.begin(),tmp.end());
                adjList.push_back(tmp);
        }
        // Edge weights (to Longlong)
        for (unsigned int i = 1; i <= numNodes; ++i)
        {
                vector<float> tmp(node_deg[i] + 1, 0);
                vector<UI> tmp1(node_deg[i] + 1, 0);
                fread(&tmp[1], sizeof(float), node_deg[i], pFile);

                for (int j = 1; j < node_deg[i] + 1; ++j)
                {
                        tmp1[j] = tmp[j] * UI_MAX;
                }

                if (tmp1[node_deg[i]] <= 0)
                        tmp1[node_deg[i]] = UI_MAX;

                weights.push_back(tmp1);
        }

        node_deg_out = vector<int>(numNodes + 1);
        fread(&node_deg_out[1], sizeof(int), numNodes, pFile); // out-degrees
        //cout<<node_deg_out[1]<<" "<< node_deg_out[2]<<" "<< node_deg_out[3]<<" "<< node_deg_out[4]<<endl;
        // Out-Neighbors
        for (unsigned int i = 1; i <= numNodes; ++i)
        {
                vector<int> tmp(node_deg_out[i]);
                fread(&tmp[0], sizeof(int), node_deg_out[i], pFile);
                adjOutList.push_back(tmp);
                // for(auto t:tmp)
                //         cout<<t<<endl;
        }
        weight_sampling = vector<ULL>(numEdges + 1);

        fread(&weight_sampling[0], sizeof(ULL), numEdges + 1, pFile); // for weight sampling
        // for(auto t:weight_sampling)
        //         cout<<t<<endl;
        // all-edges
        for (unsigned int i = 1; i <= numEdges; ++i)
        {
                vector<int> tmp(3);
                fread(&tmp[0], sizeof(int), 3, pFile);
                // for(auto t:tmp)
                //         cout<<t<<" ";
                edgeList.push_back(tmp);
                // cout<<endl;
        }
        numTriangles = 0;
}
void Graph::readIndex2(const char *filename)
{
        ifstream in(filename);
        string line;
        vector<string> array;
        getline(in, line);
                parts=vector<uint>(numNodes+1);
        succs=vector<uint>(numNodes+1);
        rates=vector<double>(numNodes+1);
        initsuccs=vector<uint>(numNodes+1);
        initsuccrates=vector<double>(numNodes+1);
        // split(line, array, ",");
        // numScc = atoi(array[0].c_str());
        // numFun = 0;
        // for (int i = 1; i < array.size(); i++)
        // {
        //         seeds.push_back(unsigned(atoi(array[i].c_str())));
        //         numFun += 1;
        // }

        // index = vector<BloomFilter>(numScc + 1);

        while (!in.eof())
        {
                string line;
                vector<string> array;
                getline(in, line);
                if(line=="")
                        break;
                split(line, array, "\t");
                // cout<<"\r"<<atoi(array[5].c_str())+1<<flush;
                parts[atoi(array[0].c_str())]=(unsigned(atoi(array[1].c_str())));
                succs[atoi(array[0].c_str())]=(unsigned(atoi(array[2].c_str())));
                initsuccs[atoi(array[0].c_str())]=(unsigned(atoi(array[4].c_str())));
                initsuccrates[atoi(array[0].c_str())]=atof(array[7].c_str());
                rates[atoi(array[0].c_str())]=atof(array[6].c_str());
        }
        in.close();
        // cout<<"done!"<<endl<<flush;
}
void Graph::writeToFile(const char *filename)
{ /*
	ofstream output(filename);
	for (unsigned int i = 0; i < numNodes; ++i){
		for (unsigned int j = 0; j < adjList[i].size(); ++j){
			if (adjList[i][j] > i){
				output << adjList[i][j] << " " << i << " " << weights[i][j] << endl;
			}
		}
	}
	output.close();
*/
}

// choose a random edge in LT model based on linear search
inline int HyperGraph::randIndex_lin(const vector<UI> &w, unsigned int si)
{
        UI ranNum = sfmt_genrand_uint32(&sfmtSeed);
        if (si <= 1 || ranNum > w[si - 1])
                return -1;

        for (unsigned int i = 1; i < si; ++i)
        {
                if (ranNum <= w[i])
                        return i;
        }
        return -1;
}

// choose a random live edge in LT model based on binary search
inline int HyperGraph::randIndex_bin(const vector<UI> &w, unsigned int si)
{
        UI ran = sfmt_genrand_uint32(&sfmtSeed);
        if (si <= 1 || ran > w[si - 1])
                return -1;
        int left = 1;
        int right = si - 1;
        int prob;
        for (unsigned int i = 0; i < si; ++i)
        {
                prob = (left + right) / 2;
                if (w[prob - 1] > ran)
                {
                        right = prob - 1;
                        continue;
                }
                if (w[prob] <= ran)
                {
                        left = prob + 1;
                        continue;
                }
                break;
        }
        return prob;
}

inline int HyperGraph::randIndex_dart(const vector<UI> &w, unsigned int si)
{
        int prob = 0;
        while (prob == 0)
        {
                UI ran = sfmt_genrand_uint32(&sfmtSeed) % si;
                UI ran2 = sfmt_genrand_uint32(&sfmtSeed);
                if (w[ran] >= ran2)
                        prob = ran;
        }
        return prob;
}

HyperGraph::HyperGraph(unsigned int n)
{
        sfmt_init_gen_rand(&sfmtSeed, rand());
        node_edge = vector<vector<int>>(n + 1);
        node_edge_int = vector<vector<int>>(n + 1);
        node_edges = vector<vector<vector<int>>>();
        node_edges.push_back(vector<vector<int>>(n + 1));
        node_edges.push_back(vector<vector<int>>(n + 1));
        node_edges.push_back(vector<vector<int>>(n + 1));
        maxDegree = 0;
        numNodes = n;
        curEdge = 0;
        nemptyset = 0;
        success = 0;
}

void HyperGraph::updateDeg()
{
        unsigned int num = edge_node.size(); //edge_node是RRset的集合
        for (unsigned int i = curEdge; i < num; ++i)
        {
                for (unsigned int k = 0; k < 3; ++k)
                {
                        unsigned int num2 = edge_node[i][k].size();
                        for (unsigned int j = 0; j < num2; ++j)
                        {
                                //cout<<"\r";
                                if (node_edge[edge_node[i][k][j]].empty() || node_edge[edge_node[i][k][j]][node_edge[edge_node[i][k][j]].size() - 1] != i)
                                        node_edge[edge_node[i][k][j]].push_back(i); //这个node出现在了哪些的rrset0中
                                node_edges[k][edge_node[i][k][j]].push_back(i);     //This node appears in which rr sets
                        }
                }
                //  cout<<"Here"<<flush;
                unsigned int num3 = edge_node[i][2].size();
                for (unsigned int j = 0; j < num3; ++j)
                {
                        if (find(edge_node[i][0].begin(), edge_node[i][0].end(), edge_node[i][2][j]) != edge_node[i][0].end() &&
                            find(edge_node[i][1].begin(), edge_node[i][1].end(), edge_node[i][2][j]) != edge_node[i][1].end())
                                node_edge_int[edge_node[i][2][j]].push_back(i); //This node appears in which rri set
                }
                // cout<<"Here"<<flush;
        }
        curEdge = edge_node.size();
        // cout<<"Here2"<<flush;
}

void HyperGraph::updateEdge()
{
        curEdge = edge_node.size();
}

/*
* Add a hyperedge into the hypergraph
*/
void HyperGraph::addEdge(vector<vector<int>> &edge)
{
        edge_node.push_back(edge);
        unsigned int ind = edge_node.size() - 1;
        for (unsigned int j = 0; j < 3; ++j)
                for (unsigned int i = 0; i < edge[j].size(); ++i)
                        node_edge[edge[j][i]].push_back(ind);
}

/*
* Add a hyperedge into the hypergraph while keeping track of the node with max degree
*/
void HyperGraph::addEdgeD(vector<vector<int>> &edge)
{
        edge_node.push_back(edge);
        int ind = edge_node.size() - 1;
        for (unsigned int j = 0; j < 3; ++j)
                for (unsigned int i = 0; i < edge[j].size(); ++i)
                {
                        node_edge[edge[j][i]].push_back(ind);
                        if (node_edge[edge[j][i]].size() > maxDegree)
                                maxDegree = node_edge[edge[j][i]].size();
                }
}

/*
* get an edge from the hypergraph
*/
const vector<vector<int>> &HyperGraph::getEdge(int e) const
{
        return edge_node[e];
}

const vector<vector<int>> &HyperGraph::getEdge(int e)
{
        return edge_node[e];
}

/*
* get the list of hyperedges incident to node n
*/
const vector<int> &HyperGraph::getNode(int n) const
{
        return node_edge[n]; //Appears in which samples 
}
const vector<int> &HyperGraph::getNodeIntersect(int n) const
{
        return node_edge_int[n]; //Appears in which samples (RRI) 
}
const vector<int> &HyperGraph::getNode(int n)
{
        return node_edge[n];
}
const vector<int> &HyperGraph::getNodeIntersect(int n)
{
        return node_edge_int[n]; //Appears in which samples (RRI) 
}
const vector<int> &HyperGraph::getNode1(int n) const
{
        return node_edges[0][n]; //Appears in which samples 
}
const vector<int> &HyperGraph::getNode1(int n)
{
        return node_edges[0][n]; //Appears in which samples 
}
const vector<int> &HyperGraph::getNode2(int n) const
{
        return node_edges[1][n]; //Appears in which samples 
}
const vector<int> &HyperGraph::getNode2(int n)
{
        return node_edges[1][n]; //Appears in which samples 
}
const vector<int> &HyperGraph::getNode3(int n) const
{
        return node_edges[2][n]; //Appears in which samples 
}
const vector<int> &HyperGraph::getNode3(int n)
{
        return node_edges[2][n]; //Appears in which samples 
}
/*
* get the number of hyperedges
*/
int HyperGraph::getNumEdge() const
{
        return edge_node.size();
}

/*
* get the maximum degree
*/
int HyperGraph::getMaxDegree()
{
        return maxDegree;
}

/*
* remove all the hyperedges
*/
void HyperGraph::clearEdges()
{
        edge_node.clear();
        node_edge.clear();
        cout << "clear edges!" << endl;
        maxDegree = 0;
}
vector<int> intersection(vector<int> &v1,
                         vector<int> &v2)
{
        vector<int> v3;
#ifndef WEIGHTEDCASCADE
        sort(v1.begin(), v1.end());
        sort(v2.begin(), v2.end());
#endif

        set_intersection(v1.begin(), v1.end(),
                         v2.begin(), v2.end(),
                         back_inserter(v3));
        return v3;
}
double HyperGraph::emptyrate()
{
        return (double)nemptyset / (double)edge_node.size();
}
struct joinULL
{
        uint a;
        uint b;
};
/*
* polling process under LT model
*/
bool HyperGraph::pollingLT2(Graph &g, vector<unsigned int> &link, unsigned int k, vector<bool> &visit, vector<int> &visit_mark)
{
        return false;
}

bool HyperGraph::pollingLT(Graph &g, vector<unsigned int> &link, unsigned int k, vector<bool> &visit, vector<int> &visit_mark)
{
   bool t = false, t1 = false, t2 = false, t3 = false;
    unsigned int i;
    unsigned int gSize = g.getSize();
    //unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%gSize+1;
    unsigned int num_marked = 0;
    unsigned int num_marked2 = 0;
    joinULL rands;
    rands.a = sfmt_genrand_uint32(&sfmtSeed);
    rands.b = sfmt_genrand_uint32(&sfmtSeed);
    ULL joinrands;
    joinrands = *(unsigned long long *)&rands;
    ULL random = joinrands % (g.getMaxWeight());
    //cout<<"Come there"<<endl;
    ULL cur = upper_bound(g.weight_sampling.begin(), g.weight_sampling.end(), random) - g.weight_sampling.begin(); //sample vertex
    //if (cur==4)
    //        cout<<cur<<endl;
    // Here the edge IDs need to be generated randomly with probability
    uint u, v, w, cur1, cur2, cur3;
    //cout<<g.edgeList.size()<<" "<<g.getMaxWeight()<<" "<<g.weight_sampling.size()<<endl;
    u = cur1 = g.edgeList[cur][1];
    v = cur2 = g.edgeList[cur][2];
    if (u == v)
        return t;
    vector<int> triangles;
    vector<int> tmp = intersection(g.adjList[cur1], g.adjList[cur2]);
    // Then traverse the triangles in which this edge is involved (here the edge may or may not be directed)
    triangles.insert(triangles.end(), tmp.begin(), tmp.end());
    tmp = intersection(g.adjOutList[cur1], g.adjOutList[cur2]);
    triangles.insert(triangles.end(), tmp.begin(), tmp.end());
    tmp = intersection(g.adjOutList[cur1], g.adjList[cur2]);
    triangles.insert(triangles.end(), tmp.begin(), tmp.end());
    tmp = intersection(g.adjList[cur1], g.adjOutList[cur2]);
    triangles.insert(triangles.end(), tmp.begin(), tmp.end());
    unordered_map<int, int> m;
    for (int j = 0; j < triangles.size(); ++j)
    {
        m[triangles[j]] += 1;
    }
    std::vector<int> keys;
    keys.reserve(m.size());
    std::vector<int> vals;
    vals.reserve(m.size());

    for (auto kv : m)
    {
        keys.push_back(kv.first);
        vals.push_back(kv.second);
    }
    vector<int> weightsample_node;
    partial_sum(vals.begin(), vals.end(), back_inserter(weightsample_node));
    w = cur3 = keys[upper_bound(weightsample_node.begin(), weightsample_node.end(), sfmt_genrand_uint32(&sfmtSeed) % (weightsample_node.back())) - weightsample_node.begin()];
    if (u == w || v == w)
        return t;

    cur = cur1;
    vector<int> rrset1, rrset2, rrset3;
    for (i = 0; i < gSize; ++i)
    {
        if (visit[cur] == true)
            break;
        visit[cur] = true;
        visit_mark[num_marked] = cur;
        rrset1.push_back(cur);
        num_marked++;
        num_marked2++;
        const vector<int> &neigh = g[cur];
        int ind;
        if (g.weights[cur].size() >= 32)
            ind = randIndex_bin(g.weights[cur], g.node_deg[cur]);
        else
            ind = randIndex_lin(g.weights[cur], g.node_deg[cur]);

        if (ind == -1)
        {
            break;
        }
        cur = neigh[ind - 1];
    }

    vector<bool> visit1 = vector<bool>(visit);
    int intersect_flag = -1;
    cur = cur2;
    for (i = 0; i < gSize; ++i)
    {
        if (visit1[cur] == true)
        {
            intersect_flag = find(visit_mark.begin(), visit_mark.begin() + num_marked, cur) - visit_mark.begin();
            break;
        }
        if (visit[cur] == true)
        {
            break;
        }
        rrset2.push_back(cur);
        visit[cur] = true;
        visit_mark[num_marked2] = cur;
        num_marked2++;
        const vector<int> &neigh = g[cur];
        int ind;
        if (g.weights[cur].size() >= 32)
            ind = randIndex_bin(g.weights[cur], g.node_deg[cur]);
        else
            ind = randIndex_lin(g.weights[cur], g.node_deg[cur]);

        if (ind == -1)
        {
            break;
        }
        cur = neigh[ind - 1];
    }
    if (intersect_flag >= 0)
        rrset2.insert(rrset2.end(), rrset1.begin() + intersect_flag, rrset1.end());
    cur = cur3;
    int intersect_flag2 = -1;
    int num_marked_2 = num_marked2;
    vector<bool> visit12 = vector<bool>(visit);
    for (i = 0; i < gSize; ++i)
    {
        if (visit12[cur] == true)
        {
            intersect_flag2 = find(visit_mark.begin(), visit_mark.begin() + num_marked_2, cur) - visit_mark.begin();
            break;
        }
        if (visit[cur] == true)
        {
            break;
        }
        rrset3.push_back(cur);
        visit[cur] = true;
        visit_mark[num_marked2] = cur;
        num_marked2++;
        const vector<int> &neigh = g[cur];
        int ind;
        if (g.weights[cur].size() >= 32)
            ind = randIndex_bin(g.weights[cur], g.node_deg[cur]);
        else
            ind = randIndex_lin(g.weights[cur], g.node_deg[cur]);

        if (ind == -1)
        {
            break;
        }
        cur = neigh[ind - 1];
    }
    if (intersect_flag2 >= 0)
    {
        if (intersect_flag2 >= num_marked)
        {
            rrset3.insert(rrset3.end(), rrset2.begin() + intersect_flag2 - num_marked, rrset2.end());
        }
        else
        {
            rrset3.insert(rrset3.end(), rrset1.begin() + intersect_flag2, rrset1.end());
        }
    }
    for (int node : rrset1)
    {
        if (link[node] < k)
        {
            t1 = true;
            break;
        }
    }
    for (int node : rrset2)
    {
        if (link[node] < k)
        {
            t2 = true;
            break;
        }
    }
    for (int node : rrset3)
    {
        if (link[node] < k)
        {
            t3 = true;
            break;
        }
    }
    for (i = 0; i < num_marked2; ++i)
    {
        visit[visit_mark[i]] = false;
    }
    return t1 && t2 && t3;
}
bool HyperGraph::pollingLT3(Graph &g, vector<unsigned int> &link, unsigned int k, vector<bool> &visit, vector<int> &visit_mark,int& node1,int& node2,int& node3)
{
   bool t = false, t1 = false, t2 = false, t3 = false;
    unsigned int i;
    unsigned int gSize = g.getSize();
    //unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%gSize+1;
    unsigned int num_marked = 0;
    unsigned int num_marked2 = 0;
    joinULL rands;
    rands.a = sfmt_genrand_uint32(&sfmtSeed);
    rands.b = sfmt_genrand_uint32(&sfmtSeed);
    ULL joinrands;
    joinrands = *(unsigned long long *)&rands;
    ULL random = joinrands % (g.getMaxWeight());
    //cout<<"Come there"<<endl;
    ULL cur = upper_bound(g.weight_sampling.begin(), g.weight_sampling.end(), random) - g.weight_sampling.begin(); //sample vertex
    //if (cur==4)
    //        cout<<cur<<endl;
    // Here the edge IDs need to be generated randomly with probability
    uint u, v, w, cur1, cur2, cur3;
    //cout<<g.edgeList.size()<<" "<<g.getMaxWeight()<<" "<<g.weight_sampling.size()<<endl;
    u = cur1 = g.edgeList[cur][1];
    v = cur2 = g.edgeList[cur][2];
    if (u == v)
        return t;
    vector<int> triangles;
    vector<int> tmp = intersection(g.adjList[cur1], g.adjList[cur2]);
    // Then traverse the triangles in which this edge is involved (here the edge may or may not be directed)
    triangles.insert(triangles.end(), tmp.begin(), tmp.end());
    tmp = intersection(g.adjOutList[cur1], g.adjOutList[cur2]);
    triangles.insert(triangles.end(), tmp.begin(), tmp.end());
    tmp = intersection(g.adjOutList[cur1], g.adjList[cur2]);
    triangles.insert(triangles.end(), tmp.begin(), tmp.end());
    tmp = intersection(g.adjList[cur1], g.adjOutList[cur2]);
    triangles.insert(triangles.end(), tmp.begin(), tmp.end());
    unordered_map<int, int> m;
    for (int j = 0; j < triangles.size(); ++j)
    {
        m[triangles[j]] += 1;
    }
    std::vector<int> keys;
    keys.reserve(m.size());
    std::vector<int> vals;
    vals.reserve(m.size());

    for (auto kv : m)
    {
        keys.push_back(kv.first);
        vals.push_back(kv.second);
    }
    vector<int> weightsample_node;
    partial_sum(vals.begin(), vals.end(), back_inserter(weightsample_node));
    w = cur3 = keys[upper_bound(weightsample_node.begin(), weightsample_node.end(), sfmt_genrand_uint32(&sfmtSeed) % (weightsample_node.back())) - weightsample_node.begin()];
    if (u == w || v == w)
        return t;
    node1=u;
    node2=v;
    node3=w;
    cur = cur1;
    vector<int> rrset1, rrset2, rrset3;
    for (i = 0; i < gSize; ++i)
    {
        if (visit[cur] == true)
            break;
        visit[cur] = true;
        visit_mark[num_marked] = cur;
        rrset1.push_back(cur);
        num_marked++;
        num_marked2++;
        const vector<int> &neigh = g[cur];
        int ind;
        if (g.weights[cur].size() >= 32)
            ind = randIndex_bin(g.weights[cur], g.node_deg[cur]);
        else
            ind = randIndex_lin(g.weights[cur], g.node_deg[cur]);

        if (ind == -1)
        {
            break;
        }
        cur = neigh[ind - 1];
    }

    vector<bool> visit1 = vector<bool>(visit);
    int intersect_flag = -1;
    cur = cur2;
    for (i = 0; i < gSize; ++i)
    {
        if (visit1[cur] == true)
        {
            intersect_flag = find(visit_mark.begin(), visit_mark.begin() + num_marked, cur) - visit_mark.begin();
            break;
        }
        if (visit[cur] == true)
        {
            break;
        }
        rrset2.push_back(cur);
        visit[cur] = true;
        visit_mark[num_marked2] = cur;
        num_marked2++;
        const vector<int> &neigh = g[cur];
        int ind;
        if (g.weights[cur].size() >= 32)
            ind = randIndex_bin(g.weights[cur], g.node_deg[cur]);
        else
            ind = randIndex_lin(g.weights[cur], g.node_deg[cur]);

        if (ind == -1)
        {
            break;
        }
        cur = neigh[ind - 1];
    }
    if (intersect_flag >= 0)
        rrset2.insert(rrset2.end(), rrset1.begin() + intersect_flag, rrset1.end());
    cur = cur3;
    int intersect_flag2 = -1;
    int num_marked_2 = num_marked2;
    vector<bool> visit12 = vector<bool>(visit);
    for (i = 0; i < gSize; ++i)
    {
        if (visit12[cur] == true)
        {
            intersect_flag2 = find(visit_mark.begin(), visit_mark.begin() + num_marked_2, cur) - visit_mark.begin();
            break;
        }
        if (visit[cur] == true)
        {
            break;
        }
        rrset3.push_back(cur);
        visit[cur] = true;
        visit_mark[num_marked2] = cur;
        num_marked2++;
        const vector<int> &neigh = g[cur];
        int ind;
        if (g.weights[cur].size() >= 32)
            ind = randIndex_bin(g.weights[cur], g.node_deg[cur]);
        else
            ind = randIndex_lin(g.weights[cur], g.node_deg[cur]);

        if (ind == -1)
        {
            break;
        }
        cur = neigh[ind - 1];
    }
    if (intersect_flag2 >= 0)
    {
        if (intersect_flag2 >= num_marked)
        {
            rrset3.insert(rrset3.end(), rrset2.begin() + intersect_flag2 - num_marked, rrset2.end());
        }
        else
        {
            rrset3.insert(rrset3.end(), rrset1.begin() + intersect_flag2, rrset1.end());
        }
    }
    for (int node : rrset1)
    {
        if (link[node] < k)
        {
            t1 = true;
            break;
        }
    }
    for (int node : rrset2)
    {
        if (link[node] < k)
        {
            t2 = true;
            break;
        }
    }
    for (int node : rrset3)
    {
        if (link[node] < k)
        {
            t3 = true;
            break;
        }
    }
    for (i = 0; i < num_marked2; ++i)
    {
        visit[visit_mark[i]] = false;
    }
    return t1 && t2 && t3;
}
void HyperGraph::pollingLT1(Graph &g, vector<bool> &visit, vector<int> &visit_mark)
{
        unsigned int i;
        unsigned int gSize = g.getSize();
        //unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%gSize+1;
        unsigned int num_marked = 0;
        unsigned int num_marked2 = 0;
        joinULL rands;
        rands.a = sfmt_genrand_uint32(&sfmtSeed);
        rands.b = sfmt_genrand_uint32(&sfmtSeed);
        ULL joinrands;
        joinrands = *(unsigned long long *)&rands;
        ULL random = joinrands % (g.getMaxWeight());
        //cout<<"Come there"<<endl;
        ULL cur = upper_bound(g.weight_sampling.begin(), g.weight_sampling.end(), random) - g.weight_sampling.begin(); //sample vertex
        //if (cur==4)
        //        cout<<cur<<endl;
        // Here the edge IDs need to be generated randomly with probability
        uint u, v, w, cur1, cur2, cur3;
        //cout<<g.edgeList.size()<<" "<<g.getMaxWeight()<<" "<<g.weight_sampling.size()<<endl;
        u = cur1 = g.edgeList[cur][1];
        v = cur2 = g.edgeList[cur][2];
        if (u == v)
                return;
        vector<int> triangles;
        vector<int> tmp = intersection(g.adjList[cur1], g.adjList[cur2]);
        // Then traverse the triangles in which this edge is involved (here the edge may or may not be directed)
        triangles.insert(triangles.end(), tmp.begin(), tmp.end());
        tmp = intersection(g.adjOutList[cur1], g.adjOutList[cur2]);
        triangles.insert(triangles.end(), tmp.begin(), tmp.end());
        tmp = intersection(g.adjOutList[cur1], g.adjList[cur2]);
        triangles.insert(triangles.end(), tmp.begin(), tmp.end());
        tmp = intersection(g.adjList[cur1], g.adjOutList[cur2]);
        triangles.insert(triangles.end(), tmp.begin(), tmp.end());
        unordered_map<int, int> m;
        for (int j = 0; j < triangles.size(); ++j)
        {
                m[triangles[j]] += 1;
        }
        std::vector<int> keys;
        keys.reserve(m.size());
        std::vector<int> vals;
        vals.reserve(m.size());

        for (auto kv : m)
        {
                keys.push_back(kv.first);
                vals.push_back(kv.second);
        }
        vector<int> weightsample_node;
        partial_sum(vals.begin(), vals.end(), back_inserter(weightsample_node));
        w = cur3 = keys[upper_bound(weightsample_node.begin(), weightsample_node.end(), sfmt_genrand_uint32(&sfmtSeed) % (weightsample_node.back())) - weightsample_node.begin()];
        if (u == w || v == w)
                return;

        cur = cur1;
        vector<int> rrset1, rrset2, rrset3;
        for (i = 0; i < gSize; ++i)
        {
                if (visit[cur] == true)
                        break;
                rrset1.push_back(cur);
                visit[cur] = true;
                visit_mark[num_marked] = cur;
                num_marked++;
                num_marked2++;
                const vector<int> &neigh = g[cur];
                int ind;
                if (g.weights[cur].size() >= 32)
                        ind = randIndex_bin(g.weights[cur], g.node_deg[cur]);
                else
                        ind = randIndex_lin(g.weights[cur], g.node_deg[cur]);

                if (ind == -1)
                {
                        break;
                }
                cur = neigh[ind - 1];
        }

        vector<bool> visit1 = vector<bool>(visit);
        int intersect_flag = -1;
        cur = cur2;
        for (i = 0; i < gSize; ++i)
        {
                if (visit1[cur] == true)
                {
                        intersect_flag = find(visit_mark.begin(), visit_mark.begin() + num_marked, cur) - visit_mark.begin();
                        break;
                }
                if (visit[cur] == true)
                {
                        break;
                }
                rrset2.push_back(cur);
                visit[cur] = true;
                visit_mark[num_marked2] = cur;
                num_marked2++;
                const vector<int> &neigh = g[cur];
                int ind;
                if (g.weights[cur].size() >= 32)
                        ind = randIndex_bin(g.weights[cur], g.node_deg[cur]);
                else
                        ind = randIndex_lin(g.weights[cur], g.node_deg[cur]);

                if (ind == -1)
                {
                        break;
                }
                cur = neigh[ind - 1];
        }

        cur = cur3;
        int intersect_flag21, intersect_flag22 = -1;
        vector<bool> visit12 = vector<bool>(visit);
        for (i = 0; i < gSize; ++i)
        {
                if (visit12[cur] == true)
                {
                        intersect_flag21 = find(rrset1.begin(), rrset1.end(), cur) - rrset1.begin();
                        if (intersect_flag21 == rrset1.size())
                        {
                                intersect_flag21 = -1;
                                intersect_flag22 = find(rrset2.begin(), rrset2.end(), cur) - rrset2.begin();
                        }
                        break;
                }
                if (visit[cur] == true)
                {
                        break;
                }
                rrset3.push_back(cur);
                visit[cur] = true;
                visit_mark[num_marked2] = cur;
                num_marked2++;
                const vector<int> &neigh = g[cur];
                int ind;
                if (g.weights[cur].size() >= 32)
                        ind = randIndex_bin(g.weights[cur], g.node_deg[cur]);
                else
                        ind = randIndex_lin(g.weights[cur], g.node_deg[cur]);

                if (ind == -1)
                {
                        break;
                }
                cur = neigh[ind - 1];
        }
        if (intersect_flag != -1)
        {
                rrset2.insert(rrset2.end(), rrset1.begin() + intersect_flag, rrset1.end());
        }
        if (intersect_flag21 != -1)
        {
                rrset3.insert(rrset3.end(), rrset1.begin() + intersect_flag21, rrset1.end());
        }
        if (intersect_flag22 != -1)
        {
                rrset3.insert(rrset3.end(), rrset2.begin() + intersect_flag22, rrset2.end());
        }
        vector<vector<int>> rrset;
        rrset.push_back(rrset1);
        rrset.push_back(rrset2);
        rrset.push_back(rrset3);
        //cout<<num_marked-intersect_flag<<" "<<num_marked<<" "<<intersect_flag<<endl;
        edge_node.push_back(rrset);
        for (i = 0; i < num_marked2; ++i)
        {
                visit[visit_mark[i]] = false;
        }
}

void HyperGraph::pollingIC1(Graph &g, vector<bool> &visit, vector<int> &visit_mark)
{
        int i;
        //cout<<"Come here"<<endl;
        joinULL rands;
        rands.a = sfmt_genrand_uint32(&sfmtSeed);
        rands.b = sfmt_genrand_uint32(&sfmtSeed);
        ULL joinrands;
        joinrands = *(unsigned long long *)&rands;
        ULL random = joinrands % (g.getMaxWeight());
        //cout<<"Come there"<<endl;
        ULL cur = upper_bound(g.weight_sampling.begin(), g.weight_sampling.end(), random) - g.weight_sampling.begin(); //sample vertex
        //if (cur==4)
        //        cout<<cur<<endl;
        // Here the edge IDs need to be generated randomly with probability
        uint u, v, w, cur1, cur2, cur3;
        //cout<<g.edgeList.size()<<" "<<g.getMaxWeight()<<" "<<g.weight_sampling.size()<<endl;
        u = cur1 = g.edgeList[cur][1];
        v = cur2 = g.edgeList[cur][2];
        if (u == v)
                return;
        vector<int> triangles;
        vector<int> tmp = intersection(g.adjList[cur1], g.adjList[cur2]);
        // Then traverse the triangles in which this edge is involved (here the edge may or may not be directed)
        triangles.insert(triangles.end(), tmp.begin(), tmp.end());
        tmp = intersection(g.adjOutList[cur1], g.adjOutList[cur2]);
        triangles.insert(triangles.end(), tmp.begin(), tmp.end());
        tmp = intersection(g.adjOutList[cur1], g.adjList[cur2]);
        triangles.insert(triangles.end(), tmp.begin(), tmp.end());
        tmp = intersection(g.adjList[cur1], g.adjOutList[cur2]);
        triangles.insert(triangles.end(), tmp.begin(), tmp.end());
        unordered_map<int, int> m;
        for (int j = 0; j < triangles.size(); ++j)
        {
                m[triangles[j]] += 1;
        }
        std::vector<int> keys;
        keys.reserve(m.size());
        std::vector<int> vals;
        vals.reserve(m.size());

        for (auto kv : m)
        {
                keys.push_back(kv.first);
                vals.push_back(kv.second);
        }
        vector<int> weightsample_node;
        partial_sum(vals.begin(), vals.end(), back_inserter(weightsample_node));
        w = cur3 = keys[upper_bound(weightsample_node.begin(), weightsample_node.end(), sfmt_genrand_uint32(&sfmtSeed) % (weightsample_node.back())) - weightsample_node.begin()];
        if (u == w || v == w)
                return;
        // Accelerated determination of empty sets
        // Degree-orientated
        // Build DFS graph
        unordered_map<int, vector<int>> ancestors;
        unordered_map<int, vector<int>> children;
        // Sampling 3 nodes
        int num_marked = 1;
        int curPos = 0;
        visit[cur1] = true;
        visit_mark[0] = cur1;
        while (curPos < num_marked)
        {
                cur1 = visit_mark[curPos];
                children[cur1].push_back(cur1);
                curPos++;
                const vector<UI> &w = g.getWeight(cur1);
                const vector<int> &neigh = g[cur1];
                uint step = 1;
                double u = sfmt_genrand_uint32(&sfmtSeed) / (double)4294967295U;
                if (w.size() > 1 && 4294967295U > w[1] && u > 0)
                        step = ceil(log(u) / log(1 - (double)w[1] / 4294967295U));
                // cout<<"\r"<<step<<" "<<u<<flush;
                for (i = step - 1; i < g.node_deg[cur1]; i += step)
                {
                        u = sfmt_genrand_uint32(&sfmtSeed) / (double)4294967295U;
                        if (4294967295U > w[1] && u > 0)
                                step = ceil(log(u) / log(1 - (double)w[1] / 4294967295U));
                        if (true)
                        { // Mark whether the edge is preserved or not by probability
                                allexpand++;
                                if (!visit[neigh[i]])
                                {
                                        visit[neigh[i]] = true;
                                        visit_mark[num_marked] = neigh[i];
                                        num_marked++; // BFS结构里节点访问完毕则停止
                                }
                                ancestors[neigh[i]].push_back(cur1);
                                if (ancestors.find(cur1) != ancestors.end())
                                        ancestors[neigh[i]].insert(ancestors[neigh[i]].end(), ancestors[cur1].begin(), ancestors[cur1].end());
                                sort(ancestors[neigh[i]].begin(), ancestors[neigh[i]].end());
                                ancestors[neigh[i]].erase(unique(ancestors[neigh[i]].begin(), ancestors[neigh[i]].end()), ancestors[neigh[i]].end());
                                for (auto anc : ancestors[cur1])
                                {
                                        if (find(children[anc].begin(), children[anc].end(), neigh[i]) == children[anc].end())
                                                children[anc].push_back(neigh[i]);
                                }
                                children[cur1].push_back(neigh[i]);
                        }
                }
        } // 
        // 
        // int numrr_1 = num_marked;
        // BloomFilter bf;
        // bf.genHash(rrset, g.hs, g.indexWidth, g.node_scc);
        if (!visit[cur2])
        {
                visit[cur2] = true;
                children[cur2].push_back(cur2);
                visit_mark[num_marked] = cur2;
                num_marked += 1;
                while (curPos < num_marked)
                {
                        cur2 = visit_mark[curPos];
                        curPos++;
                        const vector<UI> &w = g.getWeight(cur2);
                        const vector<int> &neigh = g[cur2];
                        uint step = 1;
                        double u = sfmt_genrand_uint32(&sfmtSeed) / (double)4294967295;
                        if (w.size() > 1 && 4294967295U > w[1] && u > 0)
                                step = ceil(log(u) / log(1 - (double)w[1] / 4294967295U));
                        for (i = step - 1; i < g.node_deg[cur2]; i += step)
                        {
                                u = sfmt_genrand_uint32(&sfmtSeed) / (double)4294967295U;
                                if (4294967295U > w[1] && u > 0)
                                        step = ceil(log(u) / log(1 - (double)w[1] / 4294967295U));

                                if (true)
                                { // Mark whether the edge is preserved or not by probability

                                        // if (abandoned[g.node_scc[neigh[i]]] || (g.index[g.node_scc[neigh[i]]] & bf).isZero())
                                        // {
                                        //         abandoned[g.node_scc[neigh[i]]] = true;
                                        //         success++;
                                        //         continue;
                                        // }
                                        allexpand++;
                                        if (!visit[neigh[i]])
                                        {
                                                visit[neigh[i]] = true;
                                                visit_mark[num_marked] = neigh[i];
                                                num_marked++; // BFS结构里节点访问完毕则停止
                                        }
                                        ancestors[neigh[i]].push_back(cur2);
                                        if (ancestors.find(cur2) != ancestors.end())
                                                ancestors[neigh[i]].insert(ancestors[neigh[i]].end(), ancestors[cur2].begin(), ancestors[cur2].end());
                                        sort(ancestors[neigh[i]].begin(), ancestors[neigh[i]].end());
                                        ancestors[neigh[i]].erase(unique(ancestors[neigh[i]].begin(), ancestors[neigh[i]].end()), ancestors[neigh[i]].end());
                                        for (auto anc : ancestors[cur2])
                                        {
                                                if (find(children[anc].begin(), children[anc].end(), neigh[i]) == children[anc].end())
                                                        children[anc].push_back(neigh[i]);
                                        }
                                        children[cur2].push_back(neigh[i]);
                                }
                        }
                }
        }
        // bf.genHash(rrset, g.hs, g.indexWidth, g.node_scc);
        if (!visit[cur3])
        {
                visit[cur3] = true;
                visit_mark[num_marked] = cur3;
                num_marked += 1;
                while (curPos < num_marked)
                {
                        cur3 = visit_mark[curPos];
                        children[cur3].push_back(cur3);
                        curPos++;
                        const vector<UI> &w = g.getWeight(cur3);
                        const vector<int> &neigh = g[cur3];
                        uint step = 1;
                        double u = sfmt_genrand_uint32(&sfmtSeed) / (double)4294967295;
                        if (w.size() > 1 && 4294967295U > w[1] && u > 0)
                                step = ceil(log(u) / log(1 - (double)w[1] / 4294967295U));
                        for (i = step - 1; i < g.node_deg[cur3]; i += step)
                        {
                                u = sfmt_genrand_uint32(&sfmtSeed) / (double)4294967295U;
                                if (4294967295U > w[1] && u > 0)
                                        step = ceil(log(u) / log(1 - (double)w[1] / 4294967295U));
                                if (true)
                                { // Mark whether the edge is preserved or not by probability

                                        // if (abandoned[g.node_scc[neigh[i]]] || (g.index[g.node_scc[neigh[i]]] & bf).isZero())
                                        // {
                                        //         abandoned[g.node_scc[neigh[i]]] = true;
                                        //         success++;
                                        //         continue;
                                        // }
                                        allexpand++;
                                        if (!visit[neigh[i]])
                                        {
                                                visit[neigh[i]] = true;
                                                visit_mark[num_marked] = neigh[i];
                                                num_marked++; // BFS结构里节点访问完毕则停止
                                        }
                                        ancestors[neigh[i]].push_back(cur3);
                                        if (ancestors.find(cur3) != ancestors.end())
                                                ancestors[neigh[i]].insert(ancestors[neigh[i]].end(), ancestors[cur3].begin(), ancestors[cur3].end());
                                        sort(ancestors[neigh[i]].begin(), ancestors[neigh[i]].end());
                                        ancestors[neigh[i]].erase(unique(ancestors[neigh[i]].begin(), ancestors[neigh[i]].end()), ancestors[neigh[i]].end());
                                        for (auto anc : ancestors[cur3])
                                        {
                                                if (find(children[anc].begin(), children[anc].end(), neigh[i]) == children[anc].end())
                                                        children[anc].push_back(neigh[i]);
                                        }
                                        children[cur3].push_back(neigh[i]);
                                }
                        }
                }
        }
        // 
        // TrIM
        //edge_node.push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+num_marked)); //

        vector<vector<int>> rrset;
        rrset.push_back(children[u]);
        rrset.push_back(children[v]);
        rrset.push_back(children[w]);
        edge_node.push_back(rrset);
        for (i = 0; i < num_marked; ++i)
        {
                visit[visit_mark[i]] = false;
        }
}

/*
* polling process under IC model
*/
bool HyperGraph::pollingIC2(Graph &g, vector<unsigned int> &link, unsigned int k, vector<bool> &visit, vector<int> &visit_mark)
{
        return false;
}

bool HyperGraph::pollingIC(Graph &g, vector<unsigned int> &link, unsigned int k, unordered_map<int,bool> &visit, vector<int> &visit_mark)
{
     int i;
    bool t = false, t1 = false, t2 = false, t3 = false;
    // unordered_map<int,int> visit_mark;
    // bool t = false;
    joinULL rands;
    rands.a = sfmt_genrand_uint32(&sfmtSeed);
    rands.b = sfmt_genrand_uint32(&sfmtSeed);
    ULL joinrands;
    joinrands = *(unsigned long long *)&rands;
    ULL random = joinrands % (g.getMaxWeight());
    ULL cur = upper_bound(g.weight_sampling.begin(), g.weight_sampling.end(), random) - g.weight_sampling.begin(); //sample vertex
    // Here the edge IDs need to be generated randomly with probability
    uint u, v, w, cur1, cur2, cur3;
    //cout<<g.edgeList.size()<<" "<<g.getMaxWeight()<<" "<<g.weight_sampling.size()<<endl;
    u = cur1 = g.edgeList[cur][1];
    v = cur2 = g.edgeList[cur][2];
    if (u == v)
        return t;
    vector<int> triangles;
    vector<int> tmp = intersection(g.adjList[cur1], g.adjList[cur2]);
    // Then traverse the triangles in which this edge is involved (here the edge may or may not be directed)
    triangles.insert(triangles.end(), tmp.begin(), tmp.end());
    tmp = intersection(g.adjOutList[cur1], g.adjOutList[cur2]);
    triangles.insert(triangles.end(), tmp.begin(), tmp.end());
    tmp = intersection(g.adjOutList[cur1], g.adjList[cur2]);
    triangles.insert(triangles.end(), tmp.begin(), tmp.end());
    tmp = intersection(g.adjList[cur1], g.adjOutList[cur2]);
    triangles.insert(triangles.end(), tmp.begin(), tmp.end());
    unordered_map<int, int> m;
    for (int j = 0; j < triangles.size(); ++j)
    {
        m[triangles[j]] += 1;
    }
    std::vector<int> keys;
    keys.reserve(m.size());
    std::vector<int> vals;
    vals.reserve(m.size());

    for (auto kv : m)
    {
        keys.push_back(kv.first);
        vals.push_back(kv.second);
    }
    vector<int> weightsample_node;
    partial_sum(vals.begin(), vals.end(), back_inserter(weightsample_node));
    w = cur3 = keys[upper_bound(weightsample_node.begin(), weightsample_node.end(), sfmt_genrand_uint32(&sfmtSeed) % (weightsample_node.back())) - weightsample_node.begin()];
    if (u == w || v == w)
        return t;
    // Accelerated determination of empty sets
    // Degree-orientated
    // Build DFS graph
    unordered_map<int, vector<int>> child;
    // Sampling 3 nodes
    // int num_marked = 1;
    // int curPos = 0;
    int stackTop = 0;
    // vector<bool> processed(visit.size(),false);
    unordered_map<int, bool> processed;
    visit[cur1] = true;
    visit_mark[0] = cur1;
    // vector<int> cache1;
    vector<int> rrset1;
    rrset1.push_back(cur1);
    while (stackTop >= 0)
    {
        cur1 = visit_mark[stackTop];
        if (processed[cur1])
        {
            stackTop--;
            continue;
        }
        processed[cur1] = true;
        // curPos++;
        const vector<UI> &w = g.getWeight(cur1);
        const vector<int> &neigh = g[cur1];
        uint step = 1;
        double u = sfmt_genrand_uint32(&sfmtSeed) / (double)4294967295;
        if (w.size() > 1 && 4294967295U > w[1] && u > 0)
            step = ceil(log(u) / log(1 - (double)w[1] / 4294967295U));
        int code = 0;
        for (i = step - 1; i < g.node_deg[cur1]; i += step)
        {
            u = sfmt_genrand_uint32(&sfmtSeed) / (double)4294967295U;
            if (4294967295U > w[1] && u > 0)
                step = ceil(log(u) / log(1 - (double)w[1] / 4294967295U));
            if (true)
            { // Mark whether the edge is preserved or not by probability
                allexpand++;
                child[cur1].push_back(neigh[i]);
                // prefix1[neigh[i]].push_back(prefix1[cur1][0] + to_string(code) + ".");
                if (!visit[neigh[i]])
                {
                    visit[neigh[i]] = true;
                    stackTop++;
                    visit_mark[stackTop] = neigh[i];
                    rrset1.push_back(neigh[i]);
                }
            }
        }

    } // 
    // 
    //int numrr_1 = num_marked;
    for (int node : rrset1)
    {
        if (link[node] < k)
        {
            t1 = true;
            break;
        }
    }
    vector<int> rrset2;
    unordered_map<int, bool> visit1(visit);
    set<int> inter_mark;
    if (!visit[cur2])
    {
        visit[cur2] = true;
        visit_mark[0] = cur2;
        rrset2.push_back(cur2);
        stackTop += 1;
        while (stackTop >= 0)
        {
            cur2 = visit_mark[stackTop];
            if (processed[cur2])
            {
                stackTop--;
                continue;
            }
            processed[cur2] = true;
            const vector<UI> &w = g.getWeight(cur2);
            const vector<int> &neigh = g[cur2];
            uint step = 1;
            double u = sfmt_genrand_uint32(&sfmtSeed) / (double)4294967295;
            if (w.size() > 1 && 4294967295U > w[1] && u > 0)
                step = ceil(log(u) / log(1 - (double)w[1] / 4294967295U));
            for (i = step - 1; i < g.node_deg[cur2]; i += step)
            {
                u = sfmt_genrand_uint32(&sfmtSeed) / (double)4294967295U;
                if (4294967295U > w[1] && u > 0)
                    step = ceil(log(u) / log(1 - (double)w[1] / 4294967295U));
                if (true)
                { // Mark whether the edge is preserved or not by probability
                    child[cur2].push_back(neigh[i]);
                    if (visit1[neigh[i]])
                    {
                        inter_mark.insert(neigh[i]);
                    }
                    if (!visit[neigh[i]])
                    {
                        visit[neigh[i]] = true;
                        stackTop++;
                        visit_mark[stackTop] = neigh[i];
                        rrset2.push_back(neigh[i]);
                    }
                }
            }
        }
    }
    else
    {
        inter_mark.insert(v);
    }
    int subsize = rrset1.size() + rrset2.size();
    {
        {
            unordered_map<int, bool> visit_;
            vector<UI> visit_mark_(subsize);
            int num_marked_ = 0;
            int curPos_ = 0;
            // UI cur__ =v;
            for (UI cur__ : inter_mark)
            {
                UI cur_ = cur__;
                if (!visit_[cur_])
                {
                    rrset2.push_back(cur_);
                    num_marked_ += 1;
                    visit_[cur_] = true;
                    visit_mark_[curPos_] = cur_;
                    while (curPos_ < num_marked_)
                    {
                        cur_ = visit_mark_[curPos_];
                        curPos_++;
                        const vector<int> &neigh = child[cur_];
                        // cout<<"\r"<<step<<" "<<u<<flush;
                        for (i = 0; i < child[cur_].size(); i += 1)
                        {
                            if (true)
                            { // Mark whether the edge is preserved or not by probability
                                if (!visit_[neigh[i]])
                                {
                                    visit_[neigh[i]] = true;
                                    visit_mark_[num_marked_] = neigh[i];
                                    rrset2.push_back(neigh[i]);
                                    num_marked_++; // BFS结构里节点访问完毕则停止
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for (int node : rrset2)
    {
        if (link[node] < k)
        {
            t2 = true;
            break;
        }
    }
    // cout<<tmpnm<< " "<<rrset.size()<<endl;
    vector<int> rrset3;
    unordered_map<int, bool> visit12(visit);
    set<int> inter_mark2;
    if (!visit[cur3])
    {
        rrset3.push_back(cur3);
        visit[cur3] = true;
        visit_mark[0] = cur3;
        stackTop += 1;
        while (stackTop >= 0)
        {
            cur3 = visit_mark[stackTop];
            if (processed[cur3])
            {
                stackTop--;
                continue;
            }
            processed[cur3] = true;
            const vector<UI> &w = g.getWeight(cur3);
            const vector<int> &neigh = g[cur3];
            uint step = 1;
            double u = sfmt_genrand_uint32(&sfmtSeed) / (double)4294967295;
            if (w.size() > 1 && 4294967295U > w[1] && u > 0)
                step = ceil(log(u) / log(1 - (double)w[1] / 4294967295U));
            for (i = step - 1; i < g.node_deg[cur3]; i += step)
            {
                u = sfmt_genrand_uint32(&sfmtSeed) / (double)4294967295U;
                if (4294967295U > w[1] && u > 0)
                    step = ceil(log(u) / log(1 - (double)w[1] / 4294967295U));
                if (true)
                { // Mark whether the edge is preserved or not by probability
                    child[cur3].push_back(neigh[i]);
                    if (visit12[neigh[i]])
                    {
                        inter_mark2.insert(neigh[i]);
                    }
                    if (!visit[neigh[i]])
                    {
                        visit[neigh[i]] = true;
                        stackTop++;
                        visit_mark[stackTop] = neigh[i];
                        rrset3.push_back(neigh[i]);
                    }
                }
            }
        }
    }
    else
    {
        inter_mark2.insert(w);
    }

    // 
    // HTrIM

    // set<int> rrset;
    // set<int> rrset2, rrset3;
    //repeatNode.insert(0);
    // repeatNode = true;
    if (true)
    {
        subsize += rrset3.size();

        {
            unordered_map<int, bool> visit_;
            vector<UI> visit_mark_(subsize);
            int num_marked_ = 1;
            int curPos_ = 0;
            // UI cur__ = w;
            for (UI cur__ : inter_mark2)
            {
                // cout<< curPos_<<num_marked_<<endl;
                UI cur_ = cur__;
                if (!visit_[cur_])
                {
                    // if (rrset.find(cur_) != rrset.end())
                    rrset3.push_back(cur_);
                    if (num_marked_ == curPos_)
                        num_marked_ += 1;
                    visit_[cur_] = true;
                    visit_mark_[curPos_] = cur_;
                    while (curPos_ < num_marked_)
                    {
                        cur_ = visit_mark_[curPos_];
                        curPos_++;
                        const vector<int> &neigh = child[cur_];
                        // cout<<"\r"<<step<<" "<<u<<flush;
                        for (i = 0; i < child[cur_].size(); i += 1)
                        {
                            if (true)
                            { // Mark whether the edge is preserved or not by probability
                                if (!visit_[neigh[i]])
                                {
                                    visit_[neigh[i]] = true;
                                    visit_mark_[num_marked_] = neigh[i];
                                    // if (rrset.find(cur_) != rrset.end())
                                    rrset3.push_back(neigh[i]);
                                    num_marked_++; // BFS结构里节点访问完毕则停止
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for (int node : rrset3)
    {
        if (link[node] < k)
        {
            t3 = true;
            break;
        }
    }

    // cout<<repeatNode.size()<<" "<<rrset.size()<<" "<<rrset2.size()<<endl;
    visit.clear();
    return t1 && t2 && t3;
}
bool HyperGraph::pollingIC3(Graph &g, vector<unsigned int> &link, unsigned int k, unordered_map<int,bool> &visit, vector<int> &visit_mark,int& node1,int& node2,int& node3)
{
     int i;
    bool t = false, t1 = false, t2 = false, t3 = false;
    // unordered_map<int,int> visit_mark;
    // bool t = false;
    joinULL rands;
    rands.a = sfmt_genrand_uint32(&sfmtSeed);
    rands.b = sfmt_genrand_uint32(&sfmtSeed);
    ULL joinrands;
    joinrands = *(unsigned long long *)&rands;
    ULL random = joinrands % (g.getMaxWeight());
    ULL cur = upper_bound(g.weight_sampling.begin(), g.weight_sampling.end(), random) - g.weight_sampling.begin(); //sample vertex
    // Here the edge IDs need to be generated randomly with probability
    uint u, v, w, cur1, cur2, cur3;
    //cout<<g.edgeList.size()<<" "<<g.getMaxWeight()<<" "<<g.weight_sampling.size()<<endl;
    u = cur1 = g.edgeList[cur][1];
    v = cur2 = g.edgeList[cur][2];
    if (u == v)
        return t;
    vector<int> triangles;
    vector<int> tmp = intersection(g.adjList[cur1], g.adjList[cur2]);
    // Then traverse the triangles in which this edge is involved (here the edge may or may not be directed)
    triangles.insert(triangles.end(), tmp.begin(), tmp.end());
    tmp = intersection(g.adjOutList[cur1], g.adjOutList[cur2]);
    triangles.insert(triangles.end(), tmp.begin(), tmp.end());
    tmp = intersection(g.adjOutList[cur1], g.adjList[cur2]);
    triangles.insert(triangles.end(), tmp.begin(), tmp.end());
    tmp = intersection(g.adjList[cur1], g.adjOutList[cur2]);
    triangles.insert(triangles.end(), tmp.begin(), tmp.end());
    unordered_map<int, int> m;
    for (int j = 0; j < triangles.size(); ++j)
    {
        m[triangles[j]] += 1;
    }
    std::vector<int> keys;
    keys.reserve(m.size());
    std::vector<int> vals;
    vals.reserve(m.size());

    for (auto kv : m)
    {
        keys.push_back(kv.first);
        vals.push_back(kv.second);
    }
    vector<int> weightsample_node;
    partial_sum(vals.begin(), vals.end(), back_inserter(weightsample_node));
    w = cur3 = keys[upper_bound(weightsample_node.begin(), weightsample_node.end(), sfmt_genrand_uint32(&sfmtSeed) % (weightsample_node.back())) - weightsample_node.begin()];
    if (u == w || v == w)
        return t;
    node1=u;
    node2=v;
    node3=w;    
    // Accelerated determination of empty sets
    // Degree-orientated
    // Build DFS graph
    unordered_map<int, vector<int>> child;
    // Sampling 3 nodes
    // int num_marked = 1;
    // int curPos = 0;
    int stackTop = 0;
    // vector<bool> processed(visit.size(),false);
    unordered_map<int, bool> processed;
    visit[cur1] = true;
    visit_mark[0] = cur1;
    // vector<int> cache1;
    vector<int> rrset1;
    rrset1.push_back(cur1);
    while (stackTop >= 0)
    {
        cur1 = visit_mark[stackTop];
        if (processed[cur1])
        {
            stackTop--;
            continue;
        }
        processed[cur1] = true;
        // curPos++;
        const vector<UI> &w = g.getWeight(cur1);
        const vector<int> &neigh = g[cur1];
        uint step = 1;
        double u = sfmt_genrand_uint32(&sfmtSeed) / (double)4294967295;
        if (w.size() > 1 && 4294967295U > w[1] && u > 0)
            step = ceil(log(u) / log(1 - (double)w[1] / 4294967295U));
        int code = 0;
        for (i = step - 1; i < g.node_deg[cur1]; i += step)
        {
            u = sfmt_genrand_uint32(&sfmtSeed) / (double)4294967295U;
            if (4294967295U > w[1] && u > 0)
                step = ceil(log(u) / log(1 - (double)w[1] / 4294967295U));
            if (true)
            { // Mark whether the edge is preserved or not by probability
                allexpand++;
                child[cur1].push_back(neigh[i]);
                // prefix1[neigh[i]].push_back(prefix1[cur1][0] + to_string(code) + ".");
                if (!visit[neigh[i]])
                {
                    visit[neigh[i]] = true;
                    stackTop++;
                    visit_mark[stackTop] = neigh[i];
                    rrset1.push_back(neigh[i]);
                }
            }
        }

    } // 
    // 
    //int numrr_1 = num_marked;
    for (int node : rrset1)
    {
        if (link[node] < k)
        {
            t1 = true;
            break;
        }
    }
    vector<int> rrset2;
    unordered_map<int, bool> visit1(visit);
    set<int> inter_mark;
    if (!visit[cur2])
    {
        visit[cur2] = true;
        visit_mark[0] = cur2;
        rrset2.push_back(cur2);
        stackTop += 1;
        while (stackTop >= 0)
        {
            cur2 = visit_mark[stackTop];
            if (processed[cur2])
            {
                stackTop--;
                continue;
            }
            processed[cur2] = true;
            const vector<UI> &w = g.getWeight(cur2);
            const vector<int> &neigh = g[cur2];
            uint step = 1;
            double u = sfmt_genrand_uint32(&sfmtSeed) / (double)4294967295;
            if (w.size() > 1 && 4294967295U > w[1] && u > 0)
                step = ceil(log(u) / log(1 - (double)w[1] / 4294967295U));
            for (i = step - 1; i < g.node_deg[cur2]; i += step)
            {
                u = sfmt_genrand_uint32(&sfmtSeed) / (double)4294967295U;
                if (4294967295U > w[1] && u > 0)
                    step = ceil(log(u) / log(1 - (double)w[1] / 4294967295U));
                if (true)
                { // Mark whether the edge is preserved or not by probability
                    child[cur2].push_back(neigh[i]);
                    if (visit1[neigh[i]])
                    {
                        inter_mark.insert(neigh[i]);
                    }
                    if (!visit[neigh[i]])
                    {
                        visit[neigh[i]] = true;
                        stackTop++;
                        visit_mark[stackTop] = neigh[i];
                        rrset2.push_back(neigh[i]);
                    }
                }
            }
        }
    }
    else
    {
        inter_mark.insert(v);
    }
    int subsize = rrset1.size() + rrset2.size();
    {
        {
            unordered_map<int, bool> visit_;
            vector<UI> visit_mark_(subsize);
            int num_marked_ = 0;
            int curPos_ = 0;
            // UI cur__ =v;
            for (UI cur__ : inter_mark)
            {
                UI cur_ = cur__;
                if (!visit_[cur_])
                {
                    rrset2.push_back(cur_);
                    num_marked_ += 1;
                    visit_[cur_] = true;
                    visit_mark_[curPos_] = cur_;
                    while (curPos_ < num_marked_)
                    {
                        cur_ = visit_mark_[curPos_];
                        curPos_++;
                        const vector<int> &neigh = child[cur_];
                        // cout<<"\r"<<step<<" "<<u<<flush;
                        for (i = 0; i < child[cur_].size(); i += 1)
                        {
                            if (true)
                            { // Mark whether the edge is preserved or not by probability
                                if (!visit_[neigh[i]])
                                {
                                    visit_[neigh[i]] = true;
                                    visit_mark_[num_marked_] = neigh[i];
                                    rrset2.push_back(neigh[i]);
                                    num_marked_++; // BFS结构里节点访问完毕则停止
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for (int node : rrset2)
    {
        if (link[node] < k)
        {
            t2 = true;
            break;
        }
    }
    // cout<<tmpnm<< " "<<rrset.size()<<endl;
    vector<int> rrset3;
    unordered_map<int, bool> visit12(visit);
    set<int> inter_mark2;
    if (!visit[cur3])
    {
        rrset3.push_back(cur3);
        visit[cur3] = true;
        visit_mark[0] = cur3;
        stackTop += 1;
        while (stackTop >= 0)
        {
            cur3 = visit_mark[stackTop];
            if (processed[cur3])
            {
                stackTop--;
                continue;
            }
            processed[cur3] = true;
            const vector<UI> &w = g.getWeight(cur3);
            const vector<int> &neigh = g[cur3];
            uint step = 1;
            double u = sfmt_genrand_uint32(&sfmtSeed) / (double)4294967295;
            if (w.size() > 1 && 4294967295U > w[1] && u > 0)
                step = ceil(log(u) / log(1 - (double)w[1] / 4294967295U));
            for (i = step - 1; i < g.node_deg[cur3]; i += step)
            {
                u = sfmt_genrand_uint32(&sfmtSeed) / (double)4294967295U;
                if (4294967295U > w[1] && u > 0)
                    step = ceil(log(u) / log(1 - (double)w[1] / 4294967295U));
                if (true)
                { // Mark whether the edge is preserved or not by probability
                    child[cur3].push_back(neigh[i]);
                    if (visit12[neigh[i]])
                    {
                        inter_mark2.insert(neigh[i]);
                    }
                    if (!visit[neigh[i]])
                    {
                        visit[neigh[i]] = true;
                        stackTop++;
                        visit_mark[stackTop] = neigh[i];
                        rrset3.push_back(neigh[i]);
                    }
                }
            }
        }
    }
    else
    {
        inter_mark2.insert(w);
    }

    // 
    // HTrIM

    // set<int> rrset;
    // set<int> rrset2, rrset3;
    //repeatNode.insert(0);
    // repeatNode = true;
    if (true)
    {
        subsize += rrset3.size();

        {
            unordered_map<int, bool> visit_;
            vector<UI> visit_mark_(subsize);
            int num_marked_ = 1;
            int curPos_ = 0;
            // UI cur__ = w;
            for (UI cur__ : inter_mark2)
            {
                // cout<< curPos_<<num_marked_<<endl;
                UI cur_ = cur__;
                if (!visit_[cur_])
                {
                    // if (rrset.find(cur_) != rrset.end())
                    rrset3.push_back(cur_);
                    if (num_marked_ == curPos_)
                        num_marked_ += 1;
                    visit_[cur_] = true;
                    visit_mark_[curPos_] = cur_;
                    while (curPos_ < num_marked_)
                    {
                        cur_ = visit_mark_[curPos_];
                        curPos_++;
                        const vector<int> &neigh = child[cur_];
                        // cout<<"\r"<<step<<" "<<u<<flush;
                        for (i = 0; i < child[cur_].size(); i += 1)
                        {
                            if (true)
                            { // Mark whether the edge is preserved or not by probability
                                if (!visit_[neigh[i]])
                                {
                                    visit_[neigh[i]] = true;
                                    visit_mark_[num_marked_] = neigh[i];
                                    // if (rrset.find(cur_) != rrset.end())
                                    rrset3.push_back(neigh[i]);
                                    num_marked_++; // BFS结构里节点访问完毕则停止
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for (int node : rrset3)
    {
        if (link[node] < k)
        {
            t3 = true;
            break;
        }
    }

    // cout<<repeatNode.size()<<" "<<rrset.size()<<" "<<rrset2.size()<<endl;
    visit.clear();
    return t1 && t2 && t3;
}
/*
* convert from an integer to a string
*/
string intToStr(int i)
{
        stringstream ss;
        ss << i;
        return ss.str();
}

/*
* convert from a strong to an integer
*/
unsigned int strToInt(string s)
{
        unsigned int i;
        istringstream myStream(s);

        if (myStream >> i)
        {
                return i;
        }
        else
        {
                cout << "String " << s << " is not a number." << endl;
                return atoi(s.c_str());
        }
        return i;
}

/*
* measure the consumed memory
*/
float getCurrentMemoryUsage()
{

        string pid = intToStr(unsigned(getpid()));
        string outfile = "tmp_" + pid + ".txt";
        string command = "pmap " + pid + " | grep -i Total | awk '{print $2}' > " + outfile;
        system(command.c_str());

        string mem_str;
        ifstream ifs(outfile.c_str());
        std::getline(ifs, mem_str);
        ifs.close();

        mem_str = mem_str.substr(0, mem_str.size() - 1);
        float mem = (float)strToInt(mem_str);

        command = "rm " + outfile;
        system(command.c_str());

        return mem / 1024;

        return 0;
}
