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
        unsigned int num2 = edge_node[i].size();
        for (unsigned int j = 0; j < num2; ++j)
        {
            node_edge[edge_node[i][j]].push_back(i); //这个node出现在了哪些的rrset中
        }
    }
    curEdge = edge_node.size();
}

void HyperGraph::updateEdge()
{
    curEdge = edge_node.size();
}

/*
* Add a hyperedge into the hypergraph
*/
void HyperGraph::addEdge(vector<int> &edge)
{
    edge_node.push_back(edge);
    unsigned int ind = edge_node.size() - 1;
    for (unsigned int i = 0; i < edge.size(); ++i)
        node_edge[edge[i]].push_back(ind);
}

/*
* Add a hyperedge into the hypergraph while keeping track of the node with max degree
*/
void HyperGraph::addEdgeD(vector<int> &edge)
{
    edge_node.push_back(edge);
    int ind = edge_node.size() - 1;
    for (unsigned int i = 0; i < edge.size(); ++i)
    {
        node_edge[edge[i]].push_back(ind);
        if (node_edge[edge[i]].size() > maxDegree)
            maxDegree = node_edge[edge[i]].size();
    }
}

/*
* get an edge from the hypergraph
*/
const vector<int> &HyperGraph::getEdge(int e) const
{
    return edge_node[e];
}

const vector<int> &HyperGraph::getEdge(int e)
{
    return edge_node[e];
}

/*
* get the list of hyperedges incident to node n
*/
const vector<int> &HyperGraph::getNode(int n) const
{
    return node_edge[n]; //出现在了哪些rrset
}

const vector<int> &HyperGraph::getNode(int n)
{
    return node_edge[n];
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

    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());

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
    unsigned int i;
    unsigned int gSize = g.getSize();
    bool t = false;
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
    // 这里需要按概率随机生成边的ID
    uint u, v, w, cur1, cur2, cur3;
    //cout<<g.edgeList.size()<<" "<<g.getMaxWeight()<<" "<<g.weight_sampling.size()<<endl;
    u = cur1 = g.edgeList[cur][1];
    v = cur2 = g.edgeList[cur][2];
    if (u == v)
        return t;
    vector<int> triangles;
    vector<int> tmp = intersection(g.adjList[cur1], g.adjList[cur2]);
    // 再遍历出这条边参与的三角形 （这里的边可以考虑方向，也可以不考虑）
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
    for (i = 0; i < gSize; ++i)
    {
        if (visit[cur] == true)
            break;
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
            //edge_node.push_back(vector<int>());
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
            nemptyset += 1;
            edge_node.push_back(vector<int>());
            for (i = 0; i < num_marked2; ++i)
            {
                visit[visit_mark[i]] = false;
            }
            return t;
        }
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
            nemptyset += 1;
            edge_node.push_back(vector<int>());
            for (i = 0; i < num_marked2; ++i)
            {
                visit[visit_mark[i]] = false;
            }
            return t;
        }
        cur = neigh[ind - 1];
    }

    cur = cur3;
    int intersect_flag2 = -1;
    vector<bool> visit12 = vector<bool>(visit);
    for (i = 0; i < gSize; ++i)
    {
        if (visit12[cur] == true)
        {
            intersect_flag2 = find(visit_mark.begin(), visit_mark.begin() + num_marked, cur) - visit_mark.begin();
            if (intersect_flag2 != num_marked && intersect_flag2 > intersect_flag)
                intersect_flag = intersect_flag2;
            break;
        }
        if (visit[cur] == true)
        {
            nemptyset += 1;
            edge_node.push_back(vector<int>());
            for (i = 0; i < num_marked2; ++i)
            {
                visit[visit_mark[i]] = false;
            }
            return t;
        }
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
            nemptyset += 1;
            edge_node.push_back(vector<int>());
            for (i = 0; i < num_marked2; ++i)
            {
                visit[visit_mark[i]] = false;
            }
            return t;
        }
        cur = neigh[ind - 1];
    }
    //cout<<num_marked-intersect_flag<<" "<<num_marked<<" "<<intersect_flag<<endl;
    edge_node.push_back(vector<int>(visit_mark.begin() + intersect_flag, visit_mark.begin() + num_marked));
    for (int i = intersect_flag; i < num_marked; ++i)
    {
        if (link[visit_mark[i]] < k)
        {
            t = true;
            break;
        }
    }
    for (i = 0; i < num_marked2; ++i)
    {
        visit[visit_mark[i]] = false;
    }
    return t;
}

bool HyperGraph::pollingLT(Graph &g, vector<unsigned int> &link, unsigned int k, vector<bool> &visit, vector<int> &visit_mark)
{
    unsigned int i;
    unsigned int gSize = g.getSize();
    bool t = false;
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
    // 这里需要按概率随机生成边的ID
    uint u, v, w, cur1, cur2, cur3;
    //cout<<g.edgeList.size()<<" "<<g.getMaxWeight()<<" "<<g.weight_sampling.size()<<endl;
    u = cur1 = g.edgeList[cur][1];
    v = cur2 = g.edgeList[cur][2];
    if (u == v)
        return t;
    vector<int> triangles;
    vector<int> tmp = intersection(g.adjList[cur1], g.adjList[cur2]);
    // 再遍历出这条边参与的三角形 （这里的边可以考虑方向，也可以不考虑）
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
    for (i = 0; i < gSize; ++i)
    {
        if (visit[cur] == true)
            break;
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
            //edge_node.push_back(vector<int>());
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
            //edge_node.push_back(vector<int>());
            for (i = 0; i < num_marked2; ++i)
            {
                visit[visit_mark[i]] = false;
            }
            return t;
        }
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
            //edge_node.push_back(vector<int>());
            for (i = 0; i < num_marked2; ++i)
            {
                visit[visit_mark[i]] = false;
            }
            return t;
        }
        cur = neigh[ind - 1];
    }

    cur = cur3;
    int intersect_flag2 = -1;
    vector<bool> visit12 = vector<bool>(visit);
    for (i = 0; i < gSize; ++i)
    {
        if (visit12[cur] == true)
        {
            intersect_flag2 = find(visit_mark.begin(), visit_mark.begin() + num_marked, cur) - visit_mark.begin();
            if (intersect_flag2 != num_marked && intersect_flag2 > intersect_flag)
                intersect_flag = intersect_flag2;
            break;
        }
        if (visit[cur] == true)
        {
            //edge_node.push_back(vector<int>());
            for (i = 0; i < num_marked2; ++i)
            {
                visit[visit_mark[i]] = false;
            }
            return t;
        }
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
            //edge_node.push_back(vector<int>());
            for (i = 0; i < num_marked2; ++i)
            {
                visit[visit_mark[i]] = false;
            }
            return t;
        }
        cur = neigh[ind - 1];
    }
    //cout<<num_marked-intersect_flag<<" "<<num_marked<<" "<<intersect_flag<<endl;
    //edge_node.push_back(vector<int>(visit_mark.begin()+intersect_flag,visit_mark.begin()+num_marked));
    for (int i = intersect_flag; i < num_marked; ++i)
    {
        if (link[visit_mark[i]] < k)
        {
            t = true;
            break;
        }
    }
    for (i = 0; i < num_marked2; ++i)
    {
        visit[visit_mark[i]] = false;
    }
    return t;
}

void HyperGraph::pollingLT1(Graph &g, vector<bool> &visit, vector<int> &visit_mark)
{
    // cout<<"Here!"<<flush;
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
    // 这里需要按概率随机生成边的ID
    uint u, v, w, cur1, cur2, cur3;
    //cout<<g.edgeList.size()<<" "<<g.getMaxWeight()<<" "<<g.weight_sampling.size()<<endl;
    u = cur1 = g.edgeList[cur][1];
    v = cur2 = g.edgeList[cur][2];
    if (u == v)
        return;
    vector<int> triangles;
    vector<int> tmp = intersection(g.adjList[cur1], g.adjList[cur2]);
    // 再遍历出这条边参与的三角形 （这里的边可以考虑方向，也可以不考虑）
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
    for (i = 0; i < gSize; ++i)
    {
        if (visit[cur] == true)
            break;
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
            edge_node.push_back(vector<int>());
            nemptyset += 1;
            for (i = 0; i < num_marked2; ++i)
            {
                visit[visit_mark[i]] = false;
            }
            return;
        }
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
            edge_node.push_back(vector<int>());
            nemptyset += 1;
            for (i = 0; i < num_marked2; ++i)
            {
                visit[visit_mark[i]] = false;
            }
            return;
        }
        cur = neigh[ind - 1];
    }

    cur = cur3;
    int intersect_flag2 = -1;
    vector<bool> visit12 = vector<bool>(visit);
    for (i = 0; i < gSize; ++i)
    {
        if (visit12[cur] == true)
        {
            intersect_flag2 = find(visit_mark.begin(), visit_mark.begin() + num_marked, cur) - visit_mark.begin();
            if (intersect_flag2 != num_marked && intersect_flag2 > intersect_flag)
                intersect_flag = intersect_flag2;
            break;
        }
        if (visit[cur] == true)
        {
            edge_node.push_back(vector<int>());
            nemptyset += 1;
            for (i = 0; i < num_marked2; ++i)
            {
                visit[visit_mark[i]] = false;
            }
            return;
        }
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
            edge_node.push_back(vector<int>());
            nemptyset += 1;
            for (i = 0; i < num_marked2; ++i)
            {
                visit[visit_mark[i]] = false;
            }
            return;
        }
        cur = neigh[ind - 1];
    }
    //cout<<num_marked-intersect_flag<<" "<<num_marked<<" "<<intersect_flag<<endl;
    edge_node.push_back(vector<int>(visit_mark.begin() + intersect_flag, visit_mark.begin() + num_marked));
    for (i = 0; i < num_marked2; ++i)
    {
        visit[visit_mark[i]] = false;
    }
}

void HyperGraph::pollingIC1(Graph &g, unordered_map<int, bool> &visit, vector<int> &visit_mark)
{
    // unordered_map<int,int> visit_mark;
    int i;
    bool t = false;
    joinULL rands;
    rands.a = sfmt_genrand_uint32(&sfmtSeed);
    rands.b = sfmt_genrand_uint32(&sfmtSeed);
    ULL joinrands;
    joinrands = *(unsigned long long *)&rands;
    ULL random = joinrands % (g.getMaxWeight());
    ULL cur = upper_bound(g.weight_sampling.begin(), g.weight_sampling.end(), random) - g.weight_sampling.begin(); //sample vertex
    // 这里需要按概率随机生成边的ID
    uint u, v, w, cur1, cur2, cur3;
    //cout<<g.edgeList.size()<<" "<<g.getMaxWeight()<<" "<<g.weight_sampling.size()<<endl;
    u = cur1 = g.edgeList[cur][1];
    v = cur2 = g.edgeList[cur][2];
    if (u == v)
        return;
    vector<int> triangles;
    vector<int> tmp = intersection(g.adjList[cur1], g.adjList[cur2]);
    // 再遍历出这条边参与的三角形 （这里的边可以考虑方向，也可以不考虑）
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
    // 加速判断空集
    // Degree启发
    unsigned int tmp1;
    if (g.adjList[u] < g.adjList[v])
    {
        tmp1 = u;
        u = v;
        v = tmp1;
    }
    if (g.adjList[v] < g.adjList[w])
    {
        tmp1 = v;
        v = w;
        w = tmp1;
    }
    if (g.adjList[u] < g.adjList[v])
    {
        tmp1 = u;
        u = v;
        v = tmp1;
    }
    cur1 = u;
    cur2 = v;
    cur3 = w;
    // 构建BFS图
    unordered_map<int, vector<int>> child;
    // 假设sample出三个点cur1,cur2,cur3
    int num_marked = 3;
    // int curPos = 0;
    int stackTop = 0;
    // vector<bool> processed(visit.size(),false);
    unordered_map<int, bool> processed;
    visit[cur1] = true;
    visit_mark[0] = cur1;
    unordered_map<int, int> inTime1;
    unordered_map<int, int> outTime1;
    // vector<int> inTime1(g.getSize()+1,0);
    // vector<int> outTime1(g.getSize()+1,0);
    // vector<vector<string>> prefix1(g.getSize() + 1);
    // prefix1[cur1].push_back(string("0."));
    // unordered_set<int> repeatNode;
    bool repeatNode = false;
    bool neighborCheck = false;
    // vector<int> cache1;
    int time = 1;
    while (stackTop >= 0)
    {
        cur1 = visit_mark[stackTop];
        if (processed[cur1])
        {
            stackTop--;
            outTime1[cur1] = time;
            time++;
            continue;
        }
        processed[cur1] = true;
        inTime1[cur1] = time;
        if (neighborCheck && cur1 == u)
        {
            repeatNode = true;
            for (auto c : child[cur1])
            {
                if (!visit[c])
                {
                    visit[c] = true;
                    stackTop++;
                    visit_mark[stackTop] = c;
                }
            }
            continue;
        }
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
            { // 按概率标注边是否保留
                allexpand++;
                child[cur1].push_back(neigh[i]);
                // prefix1[neigh[i]].push_back(prefix1[cur1][0] + to_string(code) + ".");
                if (!visit[neigh[i]])
                {
                    num_marked++;
                    visit[neigh[i]] = true;
                    stackTop++;
                    visit_mark[stackTop] = neigh[i];
                }
                else
                {
                    // repeatNode.insert(neigh[i]);
                    repeatNode = true;
                }
            }
        }
        time++;
        if (!neighborCheck)
        {
            neighborCheck = true;
            bool uAv = false, uAw = false;
            if (find(visit_mark.begin(), visit_mark.begin() + stackTop + 1, v) != visit_mark.begin() + stackTop + 1)
            {
                uAv = true;
            }
            if (find(visit_mark.begin(), visit_mark.begin() + stackTop + 1, cur3) != visit_mark.begin() + stackTop + 1)
            {
                uAw = true;
            }
            if (uAv || uAw)
            {
                visit.clear();
                processed[u] = false;
                // visit_mark.clear();
                if (uAv && uAw)
                {
                    visit[v] = true;
                    visit[cur3] = true;
                    visit_mark[0] = v;
                    visit_mark[1] = cur3;
                    stackTop = 1;
                }
                else
                {
                    visit[v] = uAv;
                    visit[cur3] = uAw;
                    if (uAv)
                        visit_mark[0] = v;
                    if (uAw)
                        visit_mark[0] = cur3;
                    stackTop = 0;
                }
            }
        }

    } // 无论有没有标注visited，被访问过的边不会再被访问（因为起点一定被visited了）。
    // 记录一下这里的num_mark->numrr_1
    // int numrr_1 = num_marked;
    // vector<int> rrset1(visit_mark.begin(), visit_mark.begin() + numrr_1);
    unordered_map<int, bool> visit1(visit);
    set<int> inter_mark;
    unordered_map<int, int> inTime2;
    unordered_map<int, int> outTime2;
    // vector<int> inTime2(g.getSize()+1,0);
    // vector<int> outTime2(g.getSize()+1,0);
    neighborCheck = false;
    if (!visit[cur2])
    {
        visit[cur2] = true;
        visit_mark[0] = cur2;
        stackTop += 1;
        while (stackTop >= 0)
        {
            cur2 = visit_mark[stackTop];
            if (processed[cur2])
            {
                stackTop--;
                outTime2[cur2] = time;
                if (visit1[cur2])
                {
                    inTime2[cur2] = time;
                }
                time++;
                continue;
            }
            processed[cur2] = true;
            inTime2[cur2] = time;
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
                { // 按概率标注边是否保留
                    child[cur2].push_back(neigh[i]);
                    if (visit1[neigh[i]])
                    {
                        inter_mark.insert(neigh[i]);
                        stackTop++;
                        visit_mark[stackTop] = neigh[i];
                    }
                    if (!visit[neigh[i]])
                    {

                        num_marked++;
                        visit[neigh[i]] = true;
                        stackTop++;
                        visit_mark[stackTop] = neigh[i];
                    }
                    else
                    {
                        if (!visit1[neigh[i]])
                            // repeatNode.insert(neigh[i]);
                            repeatNode = true;
                    }
                }
            }
            time++;
            if (!neighborCheck)
            {
                neighborCheck = true;
                bool vAw = false;
                if (find(visit_mark.begin(), visit_mark.begin() + stackTop + 1, cur3) != visit_mark.begin() + stackTop + 1)
                {
                    vAw = true;
                }
                if (vAw && visit1[cur3])
                {
                    inter_mark.clear();
                    inter_mark.insert(cur3);
                    repeatNode = false;
                    break;
                }
            }
        }
    }
    else
    {
        inter_mark.insert(v);
        inTime2[v] = outTime2[v] = ++time;
    }
    if (inter_mark.empty())
    {
        success += 1;
        visit.clear();
        // for (i = 0; i < num_marked; ++i)
        // {
        //         visit[visit_mark[i]] = false;
        // }
        nemptyset += 1;
        edge_node.push_back(vector<int>());
        return;
    }

    // cout<<tmpnm<< " "<<rrset.size()<<endl;
    unordered_map<int, bool> visit12(visit);
    set<int> inter_mark2;
    if (!visit[cur3])
    {
        visit[cur3] = true;
        visit_mark[0] = cur3;
        stackTop += 1;
        while (stackTop >= 0)
        {
            cur3 = visit_mark[stackTop];
            if (processed[cur3])
            {
                stackTop--;
                time++;
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
                { // 按概率标注边是否保留
                    child[cur3].push_back(neigh[i]);
                    if (visit12[neigh[i]])
                    {
                        inter_mark2.insert(neigh[i]);
                    }
                    if (!visit[neigh[i]])
                    {
                        num_marked++;
                        visit[neigh[i]] = true;
                        stackTop++;
                        visit_mark[stackTop] = neigh[i];
                    }
                }
            }
            time++;
        }
    }
    else
    {
        inter_mark2.insert(w);
    }
    if (inter_mark2.empty())
    {
        //success += 1;
        visit.clear();
        edge_node.push_back(vector<int>());
        nemptyset += 1;
        return;
    }
    // 重复上面的操作，得到三个rrset的vector。
    // HTrIM

    set<int> rrset;
    set<int> rrset2, rrset3;
    //repeatNode.insert(0);
    // repeatNode = true;
    if (repeatNode)
    {
        {
            unordered_map<int, bool> visit_;
            vector<UI> visit_mark_(num_marked);
            int num_marked_ = 0;
            int curPos_ = 0;
            // UI cur__ =v;
            for (UI cur__ : inter_mark)
            {
                UI cur_ = cur__;
                if (!visit_[cur_])
                {
                    rrset.insert(cur_);
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
                            { // 按概率标注边是否保留
                                if (!visit_[neigh[i]])
                                {
                                    visit_[neigh[i]] = true;
                                    visit_mark_[num_marked_] = neigh[i];
                                    rrset.insert(neigh[i]);
                                    num_marked_++; // BFS结构里节点访问完毕则停止
                                }
                            }
                        }
                    }
                }
            }
        }

        {
            unordered_map<int, bool> visit_;
            vector<UI> visit_mark_(num_marked);
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
                    rrset2.insert(cur_);
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
                            { // 按概率标注边是否保留
                                if (!visit_[neigh[i]])
                                {
                                    visit_[neigh[i]] = true;
                                    visit_mark_[num_marked_] = neigh[i];
                                    // if (rrset.find(cur_) != rrset.end())
                                    rrset2.insert(neigh[i]);
                                    num_marked_++; // BFS结构里节点访问完毕则停止
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        set<int> cand;
        for (auto im2 : inter_mark2)
            for (auto im : inter_mark)
            {
                if ((inTime1[im] >= inTime1[im2] && outTime1[im] <= outTime1[im2]) || (inTime2[im] >= inTime2[im2] && outTime2[im] <= outTime2[im2]))
                {
                    // cout<<inTime1[im] <<" "<< inTime1[im2] <<" "<< outTime1[im] <<" " <<outTime1[im2]<<";";
                    // cout<<inTime2[im] <<" "<< inTime2[im2] <<" "<< outTime2[im] <<" " <<outTime2[im2]<<endl;
                    cand.insert(im); //im2是im的祖先
                }
                if ((inTime1[im2] * outTime1[im2] != 0) && (inTime1[im] <= inTime1[im2] && outTime1[im] >= outTime1[im2])) //or (inTime2[im] <= inTime2[im2] && outTime2[im] >= outTime2[im2]))
                    cand.insert(im2);                                                                                      //im是im2祖先
            }
        unordered_map<int, bool> visit_;
        vector<UI> visit_mark_(num_marked);
        int num_marked_ = 1;
        int curPos_ = 0;
        // UI cur__ = w;
        for (UI cur__ : cand)
        {
            UI cur_ = cur__;
            if (!visit_[cur_])
            {
                // if (rrset.find(cur_) != rrset.end())
                rrset3.insert(cur_);
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
                        { // 按概率标注边是否保留
                            if (!visit_[neigh[i]])
                            {
                                visit_[neigh[i]] = true;
                                visit_mark_[num_marked_] = neigh[i];
                                // if (rrset.find(cur_) != rrset.end())
                                rrset3.insert(neigh[i]);
                                num_marked_++; // BFS结构里节点访问完毕则停止
                            }
                        }
                    }
                }
            }
        }
    }
    // cout<<repeatNode.size()<<" "<<rrset.size()<<" "<<rrset2.size()<<endl;
    vector<int> tmps1, tmps2, tmps3, result;
    // int tmpp1, tmpp2;
    if (!repeatNode)
    {
        // success += 1;
        tmps3 = vector<int>(rrset3.begin(), rrset3.end());
        result = tmps3;
        // tmpp1 = result.size();
    }
    else
    {
        tmps1 = vector<int>(rrset.begin(), rrset.end());
        tmps2 = vector<int>(rrset2.begin(), rrset2.end());
        result = intersection(tmps1, tmps2);
        // tmpp2 = result.size();
    }
    // if (tmpp1 > tmpp2)
    // {
    //         cout << tmpp1 << " " << tmpp2 << endl;
    //         for (auto i : tmps1)
    //                 cout << i << " " << inTime1[i] << " " << outTime1[i] << " " << inTime2[i] << " " << outTime2[i] << ";";
    //         cout << "\n";
    //         for (auto i : tmps2)
    //                 cout << i << " " << inTime1[i] << " " << outTime1[i] << " " << inTime2[i] << " " << outTime2[i] << ";";
    //         cout << "\n";
    //         for (auto i : tmps3)
    //                 cout << i << " " << inTime1[i] << " " << outTime1[i] << " " << inTime2[i] << " " << outTime2[i] << ";";
    //         cout << "\n";
    // }
    // for (i = 0; i < num_marked; ++i)
    // {
    //         visit[visit_mark[i]] = false;
    // }
    // visit = vector<bool>(visit.size(), false);
    if (result.empty())
        nemptyset += 1;
    visit.clear();
    edge_node.push_back(result);
    return;
}

/*
* polling process under IC model
*/
bool HyperGraph::pollingIC2(Graph &g, vector<unsigned int> &link, unsigned int k, unordered_map<int, bool> &visit, vector<int> &visit_mark)
{
    // unordered_map<int,int> visit_mark;
    int i;
    bool t = false;
    joinULL rands;
    rands.a = sfmt_genrand_uint32(&sfmtSeed);
    rands.b = sfmt_genrand_uint32(&sfmtSeed);
    ULL joinrands;
    joinrands = *(unsigned long long *)&rands;
    ULL random = joinrands % (g.getMaxWeight());
    ULL cur = upper_bound(g.weight_sampling.begin(), g.weight_sampling.end(), random) - g.weight_sampling.begin(); //sample vertex
    // 这里需要按概率随机生成边的ID
    uint u, v, w, cur1, cur2, cur3;
    //cout<<g.edgeList.size()<<" "<<g.getMaxWeight()<<" "<<g.weight_sampling.size()<<endl;
    u = cur1 = g.edgeList[cur][1];
    v = cur2 = g.edgeList[cur][2];
    if (u == v)
        return false;
    vector<int> triangles;
    vector<int> tmp = intersection(g.adjList[cur1], g.adjList[cur2]);
    // 再遍历出这条边参与的三角形 （这里的边可以考虑方向，也可以不考虑）
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
        return false;
    // 加速判断空集
    // Degree启发
    unsigned int tmp1;
    if (g.adjList[u] < g.adjList[v])
    {
        tmp1 = u;
        u = v;
        v = tmp1;
    }
    if (g.adjList[v] < g.adjList[w])
    {
        tmp1 = v;
        v = w;
        w = tmp1;
    }
    if (g.adjList[u] < g.adjList[v])
    {
        tmp1 = u;
        u = v;
        v = tmp1;
    }
    cur1 = u;
    cur2 = v;
    cur3 = w;
    // 构建BFS图
    unordered_map<int, vector<int>> child;
    // 假设sample出三个点cur1,cur2,cur3
    int num_marked = 3;
    // int curPos = 0;
    int stackTop = 0;
    // vector<bool> processed(visit.size(),false);
    unordered_map<int, bool> processed;
    visit[cur1] = true;
    visit_mark[0] = cur1;
    unordered_map<int, int> inTime1;
    unordered_map<int, int> outTime1;
    // vector<int> inTime1(g.getSize()+1,0);
    // vector<int> outTime1(g.getSize()+1,0);
    // vector<vector<string>> prefix1(g.getSize() + 1);
    // prefix1[cur1].push_back(string("0."));
    // unordered_set<int> repeatNode;
    bool repeatNode = false;
    bool neighborCheck = false;
    // vector<int> cache1;
    int time = 1;
    while (stackTop >= 0)
    {
        cur1 = visit_mark[stackTop];
        if (processed[cur1])
        {
            stackTop--;
            outTime1[cur1] = time;
            time++;
            continue;
        }
        processed[cur1] = true;
        inTime1[cur1] = time;
        if (neighborCheck && cur1 == u)
        {
            repeatNode = true;
            for (auto c : child[cur1])
            {
                if (!visit[c])
                {
                    visit[c] = true;
                    stackTop++;
                    visit_mark[stackTop] = c;
                }
            }
            continue;
        }
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
            { // 按概率标注边是否保留
                allexpand++;
                child[cur1].push_back(neigh[i]);
                // prefix1[neigh[i]].push_back(prefix1[cur1][0] + to_string(code) + ".");
                if (!visit[neigh[i]])
                {
                    num_marked++;
                    visit[neigh[i]] = true;
                    stackTop++;
                    visit_mark[stackTop] = neigh[i];
                }
                else
                {
                    // repeatNode.insert(neigh[i]);
                    repeatNode = true;
                }
            }
        }
        time++;
        if (!neighborCheck)
        {
            neighborCheck = true;
            bool uAv = false, uAw = false;
            if (find(visit_mark.begin(), visit_mark.begin() + stackTop + 1, v) != visit_mark.begin() + stackTop + 1)
            {
                uAv = true;
            }
            if (find(visit_mark.begin(), visit_mark.begin() + stackTop + 1, cur3) != visit_mark.begin() + stackTop + 1)
            {
                uAw = true;
            }
            if (uAv || uAw)
            {
                visit.clear();
                processed[u] = false;
                // visit_mark.clear();
                if (uAv && uAw)
                {
                    visit[v] = true;
                    visit[cur3] = true;
                    visit_mark[0] = v;
                    visit_mark[1] = cur3;
                    stackTop = 1;
                }
                else
                {
                    visit[v] = uAv;
                    visit[cur3] = uAw;
                    if (uAv)
                        visit_mark[0] = v;
                    if (uAw)
                        visit_mark[0] = cur3;
                    stackTop = 0;
                }
            }
        }

    } // 无论有没有标注visited，被访问过的边不会再被访问（因为起点一定被visited了）。
    // 记录一下这里的num_mark->numrr_1
    // int numrr_1 = num_marked;
    // vector<int> rrset1(visit_mark.begin(), visit_mark.begin() + numrr_1);
    unordered_map<int, bool> visit1(visit);
    set<int> inter_mark;
    unordered_map<int, int> inTime2;
    unordered_map<int, int> outTime2;
    neighborCheck = false;
    // vector<int> inTime2(g.getSize()+1,0);
    // vector<int> outTime2(g.getSize()+1,0);
    if (!visit[cur2])
    {
        visit[cur2] = true;
        visit_mark[0] = cur2;
        stackTop += 1;
        while (stackTop >= 0)
        {
            cur2 = visit_mark[stackTop];
            if (processed[cur2])
            {
                stackTop--;
                outTime2[cur2] = time;
                if (visit1[cur2])
                {
                    inTime2[cur2] = time;
                }
                time++;
                continue;
            }
            processed[cur2] = true;
            inTime2[cur2] = time;
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
                { // 按概率标注边是否保留
                    child[cur2].push_back(neigh[i]);
                    if (visit1[neigh[i]])
                    {
                        inter_mark.insert(neigh[i]);
                        stackTop++;
                        visit_mark[stackTop] = neigh[i];
                    }
                    if (!visit[neigh[i]])
                    {
                        num_marked++;
                        visit[neigh[i]] = true;
                        stackTop++;
                        visit_mark[stackTop] = neigh[i];
                    }
                    else
                    {
                        if (!visit1[neigh[i]])
                            // repeatNode.insert(neigh[i]);
                            repeatNode = true;
                    }
                }
            }
            time++;
            if (!neighborCheck)
            {
                neighborCheck = true;
                bool vAw = false;
                if (find(visit_mark.begin(), visit_mark.begin() + stackTop + 1, cur3) != visit_mark.begin() + stackTop + 1)
                {
                    vAw = true;
                }
                if (vAw && visit1[cur3])
                {
                    inter_mark.clear();
                    inter_mark.insert(cur3);
                    repeatNode = false;
                    break;
                }
            }
        }
    }
    else
    {
        inter_mark.insert(v);
        inTime2[v] = outTime2[v] = ++time;
    }
    if (inter_mark.empty())
    {
        success += 1;
        visit.clear();
        // for (i = 0; i < num_marked; ++i)
        // {
        //         visit[visit_mark[i]] = false;
        // }
        nemptyset += 1;
        edge_node.push_back(vector<int>());
        return false;
    }

    // cout<<tmpnm<< " "<<rrset.size()<<endl;
    unordered_map<int, bool> visit12(visit);
    set<int> inter_mark2;
    if (!visit[cur3])
    {
        visit[cur3] = true;
        visit_mark[0] = cur3;
        stackTop += 1;
        while (stackTop >= 0)
        {
            cur3 = visit_mark[stackTop];
            if (processed[cur3])
            {
                stackTop--;
                time++;
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
                { // 按概率标注边是否保留
                    child[cur3].push_back(neigh[i]);
                    if (visit12[neigh[i]])
                    {
                        inter_mark2.insert(neigh[i]);
                    }
                    if (!visit[neigh[i]])
                    {
                        num_marked++;
                        visit[neigh[i]] = true;
                        stackTop++;
                        visit_mark[stackTop] = neigh[i];
                    }
                }
            }
            time++;
        }
    }
    else
    {
        inter_mark2.insert(w);
    }
    if (inter_mark2.empty())
    {
        //success += 1;
        visit.clear();
        edge_node.push_back(vector<int>());
        nemptyset += 1;
        return false;
    }
    // 重复上面的操作，得到三个rrset的vector。
    // HTrIM

    set<int> rrset;
    set<int> rrset2, rrset3;
    //repeatNode.insert(0);
    // repeatNode = true;
    if (repeatNode)
    {
        {
            unordered_map<int, bool> visit_;
            vector<UI> visit_mark_(num_marked);
            int num_marked_ = 0;
            int curPos_ = 0;
            // UI cur__ =v;
            for (UI cur__ : inter_mark)
            {
                UI cur_ = cur__;
                if (!visit_[cur_])
                {
                    rrset.insert(cur_);
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
                            { // 按概率标注边是否保留
                                if (!visit_[neigh[i]])
                                {
                                    visit_[neigh[i]] = true;
                                    visit_mark_[num_marked_] = neigh[i];
                                    rrset.insert(neigh[i]);
                                    num_marked_++; // BFS结构里节点访问完毕则停止
                                }
                            }
                        }
                    }
                }
            }
        }

        {
            unordered_map<int, bool> visit_;
            vector<UI> visit_mark_(num_marked);
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
                    rrset2.insert(cur_);
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
                            { // 按概率标注边是否保留
                                if (!visit_[neigh[i]])
                                {
                                    visit_[neigh[i]] = true;
                                    visit_mark_[num_marked_] = neigh[i];
                                    // if (rrset.find(cur_) != rrset.end())
                                    rrset2.insert(neigh[i]);
                                    num_marked_++; // BFS结构里节点访问完毕则停止
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        set<int> cand;
        for (auto im2 : inter_mark2)
            for (auto im : inter_mark)
            {
                if ((inTime1[im] >= inTime1[im2] && outTime1[im] <= outTime1[im2]) || (inTime2[im] >= inTime2[im2] && outTime2[im] <= outTime2[im2]))
                {
                    // cout<<inTime1[im] <<" "<< inTime1[im2] <<" "<< outTime1[im] <<" " <<outTime1[im2]<<";";
                    // cout<<inTime2[im] <<" "<< inTime2[im2] <<" "<< outTime2[im] <<" " <<outTime2[im2]<<endl;
                    cand.insert(im); //im2是im的祖先
                }
                if ((inTime1[im2] * outTime1[im2] != 0) && (inTime1[im] <= inTime1[im2] && outTime1[im] >= outTime1[im2])) //or (inTime2[im] <= inTime2[im2] && outTime2[im] >= outTime2[im2]))
                    cand.insert(im2);                                                                                      //im是im2祖先
            }
        unordered_map<int, bool> visit_;
        vector<UI> visit_mark_(num_marked);
        int num_marked_ = 1;
        int curPos_ = 0;
        // UI cur__ = w;
        for (UI cur__ : cand)
        {
            UI cur_ = cur__;
            if (!visit_[cur_])
            {
                // if (rrset.find(cur_) != rrset.end())
                rrset3.insert(cur_);
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
                        { // 按概率标注边是否保留
                            if (!visit_[neigh[i]])
                            {
                                visit_[neigh[i]] = true;
                                visit_mark_[num_marked_] = neigh[i];
                                // if (rrset.find(cur_) != rrset.end())
                                rrset3.insert(neigh[i]);
                                num_marked_++; // BFS结构里节点访问完毕则停止
                            }
                        }
                    }
                }
            }
        }
    }
    // cout<<repeatNode.size()<<" "<<rrset.size()<<" "<<rrset2.size()<<endl;
    vector<int> tmps1, tmps2, tmps3, result;
    // int tmpp1, tmpp2;
    if (!repeatNode)
    {
        // success += 1;
        tmps3 = vector<int>(rrset3.begin(), rrset3.end());
        result = tmps3;
        // tmpp1 = result.size();
    }
    else
    {
        tmps1 = vector<int>(rrset.begin(), rrset.end());
        tmps2 = vector<int>(rrset2.begin(), rrset2.end());
        result = intersection(tmps1, tmps2);
        // tmpp2 = result.size();
    }
    // if (tmpp1 > tmpp2)
    // {
    //         cout << tmpp1 << " " << tmpp2 << endl;
    //         for (auto i : tmps1)
    //                 cout << i << " " << inTime1[i] << " " << outTime1[i] << " " << inTime2[i] << " " << outTime2[i] << ";";
    //         cout << "\n";
    //         for (auto i : tmps2)
    //                 cout << i << " " << inTime1[i] << " " << outTime1[i] << " " << inTime2[i] << " " << outTime2[i] << ";";
    //         cout << "\n";
    //         for (auto i : tmps3)
    //                 cout << i << " " << inTime1[i] << " " << outTime1[i] << " " << inTime2[i] << " " << outTime2[i] << ";";
    //         cout << "\n";
    // }
    //edge_node.push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+num_marked)); //求交集or把他们连起来。
    for (int node : result)
    {
        if (link[node] < k)
        {
            t = true;
            break;
        }
    }
    // for (i = 0; i < num_marked; ++i)
    // {
    //         visit[visit_mark[i]] = false;
    // }
    // visit = vector<bool>(visit.size(), false);
    if (result.empty())
        nemptyset += 1;
    visit.clear();
    edge_node.push_back(result);
    return t;
}

bool HyperGraph::pollingIC(Graph &g, vector<unsigned int> &link, unsigned int k, unordered_map<int, bool> &visit, vector<int> &visit_mark)
{
    // unordered_map<int,int> visit_mark;
    int i;
    bool t = false;
    joinULL rands;
    rands.a = sfmt_genrand_uint32(&sfmtSeed);
    rands.b = sfmt_genrand_uint32(&sfmtSeed);
    ULL joinrands;
    joinrands = *(unsigned long long *)&rands;
    ULL random = joinrands % (g.getMaxWeight());
    ULL cur = upper_bound(g.weight_sampling.begin(), g.weight_sampling.end(), random) - g.weight_sampling.begin(); //sample vertex
    // 这里需要按概率随机生成边的ID
    uint u, v, w, cur1, cur2, cur3;
    //cout<<g.edgeList.size()<<" "<<g.getMaxWeight()<<" "<<g.weight_sampling.size()<<endl;
    u = cur1 = g.edgeList[cur][1];
    v = cur2 = g.edgeList[cur][2];
    if (u == v)
        return false;
    vector<int> triangles;
    vector<int> tmp = intersection(g.adjList[cur1], g.adjList[cur2]);
    // 再遍历出这条边参与的三角形 （这里的边可以考虑方向，也可以不考虑）
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
        return false;
    // 加速判断空集
    // Degree启发
    unsigned int tmp1;
    if (g.adjList[u] < g.adjList[v])
    {
        tmp1 = u;
        u = v;
        v = tmp1;
    }
    if (g.adjList[v] < g.adjList[w])
    {
        tmp1 = v;
        v = w;
        w = tmp1;
    }
    if (g.adjList[u] < g.adjList[v])
    {
        tmp1 = u;
        u = v;
        v = tmp1;
    }
    cur1 = u;
    cur2 = v;
    cur3 = w;
    // 构建BFS图
    unordered_map<int, vector<int>> child;
    // 假设sample出三个点cur1,cur2,cur3
    // int num_marked = 1;
    // int curPos = 0;
    int stackTop = 0;
    // vector<bool> processed(visit.size(),false);
    unordered_map<int, bool> processed;
    visit[cur1] = true;
    visit_mark[0] = cur1;
    unordered_map<int, int> inTime1;
    unordered_map<int, int> outTime1;
    // vector<int> inTime1(g.getSize()+1,0);
    // vector<int> outTime1(g.getSize()+1,0);
    // vector<vector<string>> prefix1(g.getSize() + 1);
    // prefix1[cur1].push_back(string("0."));
    // unordered_set<int> repeatNode;
    bool repeatNode = false;
    int time = 0;
    while (stackTop >= 0)
    {
        cur1 = visit_mark[stackTop];
        if (processed[cur1])
        {
            stackTop--;
            outTime1[cur1] = time;
            time++;
            continue;
        }
        processed[cur1] = true;
        inTime1[cur1] = time;
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
            { // 按概率标注边是否保留
                allexpand++;
                child[cur1].push_back(neigh[i]);
                // prefix1[neigh[i]].push_back(prefix1[cur1][0] + to_string(code) + ".");
                if (!visit[neigh[i]])
                {
                    visit[neigh[i]] = true;
                    stackTop++;
                    visit_mark[stackTop] = neigh[i];
                }
                else
                {
                    // repeatNode.insert(neigh[i]);
                    repeatNode = true;
                }
            }
        }
        time++;

    } // 无论有没有标注visited，被访问过的边不会再被访问（因为起点一定被visited了）。
    // 记录一下这里的num_mark->numrr_1
    // int numrr_1 = num_marked;
    // vector<int> rrset1(visit_mark.begin(), visit_mark.begin() + numrr_1);
    unordered_map<int, bool> visit1(visit);
    set<int> inter_mark;
    unordered_map<int, int> inTime2;
    unordered_map<int, int> outTime2;
    // vector<int> inTime2(g.getSize()+1,0);
    // vector<int> outTime2(g.getSize()+1,0);
    if (!visit[cur2])
    {
        visit[cur2] = true;
        visit_mark[0] = cur2;
        stackTop += 1;
        while (stackTop >= 0)
        {
            cur2 = visit_mark[stackTop];
            if (processed[cur2])
            {
                stackTop--;
                outTime2[cur2] = time;
                if (visit1[cur2])
                {
                    inTime2[cur2] = time;
                }
                time++;
                continue;
            }
            processed[cur2] = true;
            inTime2[cur2] = time;
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
                { // 按概率标注边是否保留
                    child[cur2].push_back(neigh[i]);
                    if (visit1[neigh[i]])
                    {
                        inter_mark.insert(neigh[i]);
                        stackTop++;
                        visit_mark[stackTop] = neigh[i];
                    }
                    if (!visit[neigh[i]])
                    {
                        visit[neigh[i]] = true;
                        stackTop++;
                        visit_mark[stackTop] = neigh[i];
                    }
                    else
                    {
                        if (!visit1[neigh[i]])
                            // repeatNode.insert(neigh[i]);
                            repeatNode = true;
                    }
                }
            }
            time++;
        }
    }
    else
    {
        inter_mark.insert(v);
        inTime2[v] = outTime2[v] = ++time;
    }
    if (inter_mark.empty())
    {
        success += 1;
        visit.clear();
        // for (i = 0; i < num_marked; ++i)
        // {
        //         visit[visit_mark[i]] = false;
        // }
        return false;
    }

    // cout<<tmpnm<< " "<<rrset.size()<<endl;
    unordered_map<int, bool> visit12(visit);
    set<int> inter_mark2;
    if (!visit[cur3])
    {
        visit[cur3] = true;
        visit_mark[0] = cur3;
        stackTop += 1;
        while (stackTop >= 0)
        {
            cur3 = visit_mark[stackTop];
            if (processed[cur3])
            {
                stackTop--;
                time++;
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
                { // 按概率标注边是否保留
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
                    }
                }
            }
            time++;
        }
    }
    else
    {
        inter_mark2.insert(w);
    }
    if (inter_mark2.empty())
    {
        //success += 1;
        visit.clear();
        return false;
    }
    // 重复上面的操作，得到三个rrset的vector。
    // HTrIM

    set<int> rrset;
    set<int> rrset2, rrset3;
    //repeatNode.insert(0);
    // repeatNode = true;
    if (repeatNode)
    {
        {
            vector<bool> visit_(g.getSize() + 1, false);
            vector<UI> visit_mark_(g.getSize());
            int num_marked_ = 0;
            int curPos_ = 0;
            // UI cur__ =v;
            for (UI cur__ : inter_mark)
            {
                UI cur_ = cur__;
                if (!visit_[cur_])
                {
                    rrset.insert(cur_);
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
                            { // 按概率标注边是否保留
                                if (!visit_[neigh[i]])
                                {
                                    visit_[neigh[i]] = true;
                                    visit_mark_[num_marked_] = neigh[i];
                                    rrset.insert(neigh[i]);
                                    num_marked_++; // BFS结构里节点访问完毕则停止
                                }
                            }
                        }
                    }
                }
            }
        }

        {
            vector<bool> visit_(g.getSize() + 1, false);
            vector<UI> visit_mark_(g.getSize());
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
                    rrset2.insert(cur_);
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
                            { // 按概率标注边是否保留
                                if (!visit_[neigh[i]])
                                {
                                    visit_[neigh[i]] = true;
                                    visit_mark_[num_marked_] = neigh[i];
                                    // if (rrset.find(cur_) != rrset.end())
                                    rrset2.insert(neigh[i]);
                                    num_marked_++; // BFS结构里节点访问完毕则停止
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        set<int> cand;
        for (auto im2 : inter_mark2)
            for (auto im : inter_mark)
            {
                if ((inTime1[im] >= inTime1[im2] && outTime1[im] <= outTime1[im2]) || (inTime2[im] >= inTime2[im2] && outTime2[im] <= outTime2[im2]))
                {
                    // cout<<inTime1[im] <<" "<< inTime1[im2] <<" "<< outTime1[im] <<" " <<outTime1[im2]<<";";
                    // cout<<inTime2[im] <<" "<< inTime2[im2] <<" "<< outTime2[im] <<" " <<outTime2[im2]<<endl;
                    cand.insert(im); //im2是im的祖先
                }
                if ((inTime1[im2] * outTime1[im2] != 0) && (inTime1[im] <= inTime1[im2] && outTime1[im] >= outTime1[im2])) //or (inTime2[im] <= inTime2[im2] && outTime2[im] >= outTime2[im2]))
                    cand.insert(im2);                                                                                      //im是im2祖先
            }
        vector<bool> visit_(g.getSize() + 1, false);
        vector<UI> visit_mark_(g.getSize());
        int num_marked_ = 1;
        int curPos_ = 0;
        // UI cur__ = w;
        for (UI cur__ : cand)
        {
            UI cur_ = cur__;
            if (!visit_[cur_])
            {
                // if (rrset.find(cur_) != rrset.end())
                rrset3.insert(cur_);
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
                        { // 按概率标注边是否保留
                            if (!visit_[neigh[i]])
                            {
                                visit_[neigh[i]] = true;
                                visit_mark_[num_marked_] = neigh[i];
                                // if (rrset.find(cur_) != rrset.end())
                                rrset3.insert(neigh[i]);
                                num_marked_++; // BFS结构里节点访问完毕则停止
                            }
                        }
                    }
                }
            }
        }
    }
    // cout<<repeatNode.size()<<" "<<rrset.size()<<" "<<rrset2.size()<<endl;
    vector<int> tmps1, tmps2, tmps3, result;
    // int tmpp1, tmpp2;
    if (!repeatNode)
    {
        // success += 1;
        tmps3 = vector<int>(rrset3.begin(), rrset3.end());
        result = tmps3;
        // tmpp1 = result.size();
    }
    else
    {
        tmps1 = vector<int>(rrset.begin(), rrset.end());
        tmps2 = vector<int>(rrset2.begin(), rrset2.end());
        result = intersection(tmps1, tmps2);
        // tmpp2 = result.size();
    }
    // if (tmpp1 > tmpp2)
    // {
    //         cout << tmpp1 << " " << tmpp2 << endl;
    //         for (auto i : tmps1)
    //                 cout << i << " " << inTime1[i] << " " << outTime1[i] << " " << inTime2[i] << " " << outTime2[i] << ";";
    //         cout << "\n";
    //         for (auto i : tmps2)
    //                 cout << i << " " << inTime1[i] << " " << outTime1[i] << " " << inTime2[i] << " " << outTime2[i] << ";";
    //         cout << "\n";
    //         for (auto i : tmps3)
    //                 cout << i << " " << inTime1[i] << " " << outTime1[i] << " " << inTime2[i] << " " << outTime2[i] << ";";
    //         cout << "\n";
    // }
    //edge_node.push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+num_marked)); //求交集or把他们连起来。
    for (int node : result)
    {
        if (link[node] < k)
        {
            t = true;
            break;
        }
    }
    // for (i = 0; i < num_marked; ++i)
    // {
    //         visit[visit_mark[i]] = false;
    // }
    // visit = vector<bool>(visit.size(), false);
    visit.clear();
    return t;
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
