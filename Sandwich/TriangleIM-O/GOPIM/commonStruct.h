#pragma once
typedef uint32_t UI;
typedef uint64_t ULL;
bool app_stopped = false;
struct joinULL
{
    uint a;
    uint b;
};
/// Node list
typedef std::vector<uint32_t> Nodelist;
/// Edge structure, neighbor id and the edge weight
typedef std::pair<uint32_t, float> Edge;
/// Edgelist structure from one source/target node
typedef std::vector<Edge> Edgelist;
/// Graph structure
typedef std::vector<Edgelist> Graph;
/// Reverse Graph structure
typedef std::vector<std::vector<uint32_t>> RevGraph;
/// One forward reachable set
typedef std::vector<size_t> FRset;
/// A set of forward reachable sets
typedef std::vector<FRset> FRsets;
/// A set of forward reachable sets
typedef std::vector<FRsets> FRseqs;
/// One reverse reachable set
typedef std::vector<uint32_t> RRset;
/// A set of reverse reachable sets
typedef std::vector<RRset> RRseq;
/// A set of reverse reachable sequences
typedef std::vector<RRseq> RRcol;
/// One weight sampling vector for triple
typedef std::vector<uint64_t> TriWeights;
/// EdgeList for sample
typedef std::vector<std::vector<uint32_t>> SpEdgeList;
/// Cascade models: IC, LT
enum CascadeModel { IC, LT };

/// Node element with id and a property value
typedef struct NodeElement
{
	int id;
	double value;
} NodeEleType;

typedef struct GraphInformation
{
	int numNodes;
	long long numEdges;
	long long numTriangles;
	Graph graph;
	RevGraph revgraph;
	TriWeights triweights;
	SpEdgeList spedgelist;
	Nodelist node_deg;
	Nodelist node_deg_out;
} GraphInfo;
/// Smaller operation for node element
struct smaller
{
	bool operator()(const NodeEleType& Ele1, const NodeEleType& Ele2) const
	{
		return (Ele1.value < Ele2.value);
	}
};

/// Greater operation for node element
struct greater
{
	bool operator()(const NodeEleType& Ele1, const NodeEleType& Ele2) const
	{
		return (Ele1.value > Ele2.value);
	}
};
