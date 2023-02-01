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
    edgeList.push_back(d);   // all edges
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
    // cout<<node_deg_out[1]<<" "<< node_deg_out[2]<<" "<< node_deg_out[3]<<" "<< node_deg_out[4]<<endl;
    //  Out-Neighbors
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
        tokens.push_back(s.substr(lastPos, pos - lastPos)); // use emplace_back after C++11
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
    fread(&numEdges, sizeof(long long), 1, pFile); // Edges
    node_deg = vector<int>(numNodes + 1);
    fread(&node_deg[1], sizeof(int), numNodes, pFile); // in-degrees

    vector<int> a, c, d;
    vector<UI> b;
    adjList.push_back(a);    // in-neighbors
    adjOutList.push_back(c); // out-neighbors
    edgeList.push_back(d);   // all edges
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
    // cout<<node_deg_out[1]<<" "<< node_deg_out[2]<<" "<< node_deg_out[3]<<" "<< node_deg_out[4]<<endl;
    //  Out-Neighbors
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
    unsigned int num = edge_node.size(); // edge_node是RRset的集合
    for (unsigned int i = curEdge; i < num; ++i)
    {
        unsigned int num2 = edge_node[i].size();
        for (unsigned int j = 0; j < num2; ++j)
        {
            node_edge[edge_node[i][j]].push_back(i); // 这个node出现在了哪些的rrset中
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
    return node_edge[n]; // 出现在了哪些rrset
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
    bool t = false;
    unsigned int gSize = g.getSize();
    joinULL rands;
    rands.a = sfmt_genrand_uint32(&sfmtSeed);
    rands.b = sfmt_genrand_uint32(&sfmtSeed);
    ULL joinrands;
    joinrands = *(unsigned long long *)&rands;
    ULL random = joinrands % (g.getMaxWeight());
    // cout<<"Come there"<<endl;
    ULL cur = upper_bound(g.weight_sampling.begin(), g.weight_sampling.end(), random) - g.weight_sampling.begin(); // sample vertex
    // if (cur==4)
    //         cout<<cur<<endl;
    //  这里需要按概率随机生成边的ID
    uint u, v, w, cur1, cur2, cur3;
    // cout<<g.edgeList.size()<<" "<<g.getMaxWeight()<<" "<<g.weight_sampling.size()<<endl;
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
    auto selected = rand() % 3;
    switch (selected)
    {
    case 0:
        cur = cur1;
        break;
    case 1:
        cur = cur2;
        break;
    case 2:
        cur = cur3;
        break;
    default:
        break;
    }
    unsigned int num_marked = 0;
    for (i = 0; i < gSize; ++i)
    {
        if (visit[cur] == true)
            break;
        visit[cur] = true;
        visit_mark[num_marked] = cur;
        num_marked++;
        if (link[cur] < k)
            t = true;
        int ind;
        if (g.weights[cur].size() >= 32)
            ind = randIndex_bin(g.weights[cur], g.node_deg[cur]);
        else
            ind = randIndex_lin(g.weights[cur], g.node_deg[cur]);

        if (ind == -1)
            break;

        cur = g.adjList[cur][ind - 1];
    }
    edge_node.push_back(vector<int>(visit_mark.begin(), visit_mark.begin() + num_marked));
    for (i = 0; i < num_marked; ++i)
    {
        visit[visit_mark[i]] = false;
    }
    return t;
}

bool HyperGraph::pollingLT(Graph &g, vector<unsigned int> &link, unsigned int k, vector<bool> &visit, vector<int> &visit_mark)
{
    unsigned int i;
    bool t = false;
    unsigned int gSize = g.getSize();
    joinULL rands;
    rands.a = sfmt_genrand_uint32(&sfmtSeed);
    rands.b = sfmt_genrand_uint32(&sfmtSeed);
    ULL joinrands;
    joinrands = *(unsigned long long *)&rands;
    ULL random = joinrands % (g.getMaxWeight());
    // cout<<"Come there"<<endl;
    ULL cur = upper_bound(g.weight_sampling.begin(), g.weight_sampling.end(), random) - g.weight_sampling.begin(); // sample vertex
    // if (cur==4)
    //         cout<<cur<<endl;
    //  这里需要按概率随机生成边的ID
    uint u, v, w, cur1, cur2, cur3;
    // cout<<g.edgeList.size()<<" "<<g.getMaxWeight()<<" "<<g.weight_sampling.size()<<endl;
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
    auto selected = rand() % 3;
    switch (selected)
    {
    case 0:
        cur = cur1;
        break;
    case 1:
        cur = cur2;
        break;
    case 2:
        cur = cur3;
        break;
    default:
        break;
    }

    unsigned int num_marked = 0;
    for (i = 0; i < gSize; ++i)
    {
        if (link[cur] < k)
        {
            t = true;
            break;
        }
        if (visit[cur] == true)
            break;
        visit[cur] = true;
        visit_mark[num_marked] = cur;
        num_marked++;
        int ind;
        if (g.weights[cur].size() >= 32)
            ind = randIndex_bin(g.weights[cur], g.node_deg[cur]);
        else
            ind = randIndex_lin(g.weights[cur], g.node_deg[cur]);

        if (ind == -1)
            break;
        cur = g.adjList[cur][ind - 1];
    }
    for (i = 0; i < num_marked; ++i)
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
    // unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%gSize+1;
    unsigned int num_marked = 0;
    joinULL rands;
    rands.a = sfmt_genrand_uint32(&sfmtSeed);
    rands.b = sfmt_genrand_uint32(&sfmtSeed);
    ULL joinrands;
    joinrands = *(unsigned long long *)&rands;
    ULL random = joinrands % (g.getMaxWeight());
    // cout<<"Come there"<<endl;
    ULL cur = upper_bound(g.weight_sampling.begin(), g.weight_sampling.end(), random) - g.weight_sampling.begin(); // sample vertex
    // if (cur==4)
    //         cout<<cur<<endl;
    //  这里需要按概率随机生成边的ID
    uint u, v, w, cur1, cur2, cur3;
    // cout<<g.edgeList.size()<<" "<<g.getMaxWeight()<<" "<<g.weight_sampling.size()<<endl;
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
    auto selected = rand() % 3;
    switch (selected)
    {
    case 0:
        cur = cur1;
        break;
    case 1:
        cur = cur2;
        break;
    case 2:
        cur = cur3;
        break;
    default:
        break;
    }
    for (i = 0; i < gSize; ++i)
    {
        if (visit[cur] == true)
            break;
        visit[cur] = true;
        visit_mark[num_marked] = cur;
        num_marked++;
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

    // cout<<num_marked-intersect_flag<<" "<<num_marked<<" "<<intersect_flag<<endl;
    edge_node.push_back(vector<int>(visit_mark.begin(), visit_mark.begin() + num_marked));
    for (i = 0; i < num_marked; ++i)
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
    ULL cur = upper_bound(g.weight_sampling.begin(), g.weight_sampling.end(), random) - g.weight_sampling.begin(); // sample vertex
    // 这里需要按概率随机生成边的ID
    uint u, v, w, cur1, cur2, cur3;
    // cout<<g.edgeList.size()<<" "<<g.getMaxWeight()<<" "<<g.weight_sampling.size()<<endl;
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
    auto selected = rand() % 3;
    switch (selected)
    {
    case 0:
        cur = cur1;
        break;
    case 1:
        cur = cur2;
        break;
    case 2:
        cur = cur3;
        break;
    default:
        break;
    }
    int num_marked = 1;
    int curPos = 0;
    visit[cur] = true;
    visit_mark[0] = cur;
    while (curPos < num_marked)
    {
        cur = visit_mark[curPos];
        curPos++;
        const vector<UI> &w = g.getWeight(cur);
        const vector<int> &neigh = g[cur];
        for (i = 0; i < g.node_deg[cur]; ++i)
        {
            if (sfmt_genrand_uint32(&sfmtSeed) < w[i + 1])
            {
                if (!visit[neigh[i]])
                {
                    visit[neigh[i]] = true;
                    visit_mark[num_marked] = neigh[i];
                    num_marked++;
                }
            }
        }
    }
    edge_node.push_back(vector<int>(visit_mark.begin(), visit_mark.begin() + num_marked));
    for (i = 0; i < num_marked; ++i)
    {
        visit[visit_mark[i]] = false;
    }
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
    ULL cur = upper_bound(g.weight_sampling.begin(), g.weight_sampling.end(), random) - g.weight_sampling.begin(); // sample vertex
    // 这里需要按概率随机生成边的ID
    uint u, v, w, cur1, cur2, cur3;
    // cout<<g.edgeList.size()<<" "<<g.getMaxWeight()<<" "<<g.weight_sampling.size()<<endl;
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
    auto selected = rand() % 3;
    switch (selected)
    {
    case 0:
        cur = cur1;
        break;
    case 1:
        cur = cur2;
        break;
    case 2:
        cur = cur3;
        break;
    default:
        break;
    }
    int curPos = 0;
    int num_marked = 1;
    visit[cur] = true;
    visit_mark[0] = cur;
    // bool t = false;
    while (curPos < num_marked)
    {
        cur = visit_mark[curPos];
        curPos++;
        if (link[cur] < k)
            t = true;
        const vector<UI> &w = g.getWeight(cur);
        const vector<int> &neigh = g[cur];
        for (i = 0; i < g.node_deg[cur]; ++i)
        {
            if (sfmt_genrand_uint32(&sfmtSeed) < w[i + 1])
            {
                if (!visit[neigh[i]])
                {
                    visit[neigh[i]] = true;
                    visit_mark[num_marked] = neigh[i];
                    num_marked++;
                }
            }
        }
    }
    edge_node.push_back(vector<int>(visit_mark.begin(), visit_mark.begin() + num_marked));

    for (i = 0; i < num_marked; ++i)
    {
        visit[visit_mark[i]] = false;
    }
    return t;
}

bool HyperGraph::pollingIC(Graph &g, vector<unsigned int> &link, unsigned int k, unordered_map<int, bool> &visit, vector<int> &visit_mark)
{
    // unordered_map<int,int> visit_mark;
    // unordered_map<int,int> visit_mark;
    int i;
    bool t = false;
    joinULL rands;
    rands.a = sfmt_genrand_uint32(&sfmtSeed);
    rands.b = sfmt_genrand_uint32(&sfmtSeed);
    ULL joinrands;
    joinrands = *(unsigned long long *)&rands;
    ULL random = joinrands % (g.getMaxWeight());
    ULL cur = upper_bound(g.weight_sampling.begin(), g.weight_sampling.end(), random) - g.weight_sampling.begin(); // sample vertex
    // 这里需要按概率随机生成边的ID
    uint u, v, w, cur1, cur2, cur3;
    // cout<<g.edgeList.size()<<" "<<g.getMaxWeight()<<" "<<g.weight_sampling.size()<<endl;
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
    auto selected = rand() % 3;
    switch (selected)
    {
    case 0:
        cur = cur1;
        break;
    case 1:
        cur = cur2;
        break;
    case 2:
        cur = cur3;
        break;
    default:
        break;
    }
    int curPos = 0;
    int num_marked = 1;
    visit[cur] = true;
    visit_mark[0] = cur;
    // bool t = false;
    while (curPos < num_marked)
    {
        cur = visit_mark[curPos];
        curPos++;
        if (link[cur] < k)
        {
            t = true;
            return t;
        }
        const vector<UI> &w = g.getWeight(cur);
        const vector<int> &neigh = g[cur];
        for (i = 0; i < g.node_deg[cur]; ++i)
        {
            if (sfmt_genrand_uint32(&sfmtSeed) < w[i + 1])
            {
                if (!visit[neigh[i]])
                {
                    visit[neigh[i]] = true;
                    visit_mark[num_marked] = neigh[i];
                    num_marked++;
                }
            }
        }
    }
    // edge_node.push_back(vector<int>(visit_mark.begin(), visit_mark.begin() + num_marked));

    for (i = 0; i < num_marked; ++i)
    {
        visit[visit_mark[i]] = false;
    }
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
