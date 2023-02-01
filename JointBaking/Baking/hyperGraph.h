#pragma once
#include <algorithm>
#include <unordered_map>
#include <set>
class HyperGraph
{
private:
	/// __numV: number of nodes in the graph.
	uint32_t __numV;
	/// __numE: number of edges in the graph.
	size_t __numE;
	/// __numT: number of directed triangles in the graph.
	unsigned long long __numT;
	/// __numRRsets: number of RR sets.
	size_t __numRRsets = 0;
	// std::vector<bool> __vecVisitBool;
	std::unordered_map<uint, bool> __vecVisitBool;
	Nodelist __vecVisitNode;

	/// Initialization
	void init_hypergraph()
	{
		__numV = (uint32_t)_graphinfo.numNodes;
		// _graph=_graphinfo.graph;
		// for (auto& nbrs : _graph) __numE += nbrs.size();
		__numE = _graphinfo.numEdges;
		__numT = _graphinfo.numTriangles;
		_FRsets = FRsets(__numV);
		_FRsetsU = FRsets(__numV);
		// __vecVisitBool = std::vector<bool>(__numV);
		__vecVisitNode = Nodelist(__numV);
	}

public:
	GraphInfo &_graphinfo;
	/// _graph: reverse graph
	Graph &_graph;
	RevGraph &_revgraph;
	/// _FRsets: forward cover sets, _FRsets[i] is the node sets that node i can reach
	FRsets _FRsets;
	FRsets _FRsetsU;
	/// _RRsets: reverse cover sets, _RRsets[i] is the node set that can reach node i
	RRsets _RRsets;
	RRsets _RRsetsU;
	bool stopLower = false;
	bool stopUpper = false;
	/// _cascadeModel: the cascade model, default is IC
	CascadeModel _cascadeModel = IC;

	explicit HyperGraph(GraphInfo &graphinfo) : _graphinfo(graphinfo), _graph(graphinfo.graph), _revgraph(graphinfo.revgraph)
	{
		init_hypergraph();
	}

	/// Set cascade model
	void set_cascade_model(const CascadeModel model)
	{
		_cascadeModel = model;
	}

	/// Returns the number of nodes in the graph.
	uint32_t get_nodes() const
	{
		return __numV;
	}

	/// Returns the number of edges in the graph.
	size_t get_edges() const
	{
		return __numE;
	}
	/// Returns the number of edges in the graph.
	size_t get_triangles() const
	{
		return __numT;
	}
	/// Returns the number of RR sets in the graph.
	size_t get_RR_sets_size() const
	{
		return __numRRsets;
	}

	/// Get out degree
	std::vector<size_t> get_out_degree() const
	{
		std::vector<size_t> outDeg(__numV);
		for (auto &nbrs : _graph)
		{
			for (auto &nbr : nbrs)
			{
				outDeg[nbr.first]++;
			}
		}
		return outDeg;
	}
	std::vector<uint> intersection(std::vector<uint> &v1,
								   std::vector<uint> &v2)
	{
		std::vector<uint> v3;
		std::set_intersection(v1.begin(), v1.end(),
							  v2.begin(), v2.end(),
							  std::back_inserter(v3));
		return v3;
	}
	std::vector<uint> intersection(std::vector<Edge> &es1,
								   std::vector<uint> &v2)
	{
		std::vector<uint> v1(es1.size()), v3;
		int count = 0;
		for (auto &e1 : es1)
		{
			v1[count++] = e1.first;
		}
		std::set_intersection(v1.begin(), v1.end(),
							  v2.begin(), v2.end(),
							  std::back_inserter(v3));
		return v3;
	}
	std::vector<uint> intersection(std::vector<uint> &v2, std::vector<Edge> &es1)
	{
		std::vector<uint> v1(es1.size()), v3;
		int count = 0;
		for (auto &e1 : es1)
		{
			v1[count++] = e1.first;
		}
		std::set_intersection(v1.begin(), v1.end(),
							  v2.begin(), v2.end(),
							  std::back_inserter(v3));
		return v3;
	}
	std::vector<uint> intersection(std::vector<Edge> &es1,
								   std::vector<Edge> &es2)
	{
		std::vector<uint> v1(es1.size()), v2(es2.size()), v3;
		int count = 0;
		for (auto &e1 : es1)
		{
			v1[count++] = e1.first;
		}
		count = 0;
		for (auto &e2 : es2)
		{
			v2[count++] = e2.first;
		}
		std::set_intersection(v1.begin(), v1.end(),
							  v2.begin(), v2.end(),
							  std::back_inserter(v3));
		return v3;
	}
	/// Generate a set of n RR sets
	void build_n_RRsets(const size_t numSamples)
	{
		if (numSamples > SIZE_MAX)
		{
			std::cout << "Error:R too large" << std::endl;
			exit(1);
		}
		const auto prevSize = __numRRsets;
		__numRRsets = __numRRsets > numSamples ? __numRRsets : numSamples;
		for (auto i = prevSize; i < numSamples; i++)
		{

			build_one_RRset(i);
		}
	}

	/// Generate one RR set
	void build_one_RRset(const size_t hyperIdx)
	{
		// //std::cout << "\rHere!" << 0 << std::flush;
		int i;
		bool t = false;
		joinULL rands;
		rands.a = dsfmt_gv_genrand_uint32();
		rands.b = dsfmt_gv_genrand_uint32();
		ULL joinrands;
		joinrands = *(unsigned long long *)&rands;
		ULL random = joinrands % (_graphinfo.triweights.back());
		ULL cur = std::upper_bound(_graphinfo.triweights.begin(), _graphinfo.triweights.end(), random) - _graphinfo.triweights.begin(); // sample vertex
		// Here the edge IDs need to be generated randomly with probability
		uint u, v, w, cur1, cur2, cur3;
		// cout<<_graphinfo.spedgelist.size()<<" "<<g.getMaxWeight()<<" "<<_graphinfo.triweights.size()<<endl;
		u = cur1 = _graphinfo.spedgelist[cur][1] - 1;
		v = cur2 = _graphinfo.spedgelist[cur][2] - 1;
		if (u == v)
			return;
		std::vector<int> triangles;
		std::vector<uint> tmp = intersection(_graph[cur1], _graph[cur2]);
		// Then traverse the triangles in which this edge is involved (here the edge may or may not be directed)
		triangles.insert(triangles.end(), tmp.begin(), tmp.end());
		tmp = intersection(_revgraph[cur1], _revgraph[cur2]);
		triangles.insert(triangles.end(), tmp.begin(), tmp.end());
		tmp = intersection(_revgraph[cur1], _graph[cur2]);
		triangles.insert(triangles.end(), tmp.begin(), tmp.end());
		tmp = intersection(_graph[cur1], _revgraph[cur2]);
		triangles.insert(triangles.end(), tmp.begin(), tmp.end());
		std::unordered_map<int, int> m;
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
		std::vector<int> weightsample_node;
		partial_sum(vals.begin(), vals.end(), back_inserter(weightsample_node));
		w = cur3 = keys[std::upper_bound(weightsample_node.begin(), weightsample_node.end(), dsfmt_gv_genrand_uint32() % (weightsample_node.back())) - weightsample_node.begin()];
		if (u == w || v == w)
			return;
		// Accelerated determination of empty sets
		// Degree-orientated
		unsigned int tmp1;
		if (_graph[u].size() < _graph[v].size())
		{
			tmp1 = u;
			u = v;
			v = tmp1;
		}
		if (_graph[v].size() < _graph[w].size())
		{
			tmp1 = v;
			v = w;
			w = tmp1;
		}
		if (_graph[u].size() < _graph[v].size())
		{
			tmp1 = u;
			u = v;
			v = tmp1;
		}
		cur1 = u;
		cur2 = v;
		cur3 = w;
		auto selected = rand() % 3;
		auto forU = 0;
		switch (selected)
		{
		case 0:
			forU = u;
			break;
		case 1:
			forU = v;
			break;
		case 2:
			forU = w;
			break;
		default:
			break;
		}
		if (_cascadeModel == IC)
		{
			// Build DFS graph
			std::unordered_map<uint, std::vector<uint>> child;
			// Sampling 3 nodes
			int num_marked = 3;
			// int curPos = 0;
			int stackTop = 0;
			// std::vector<bool> processed(visit.size(),false);
			std::unordered_map<uint, bool> processed;
			__vecVisitBool[cur1] = true;
			__vecVisitNode[0] = cur1;
			std::unordered_map<uint, uint> inTime1;
			std::unordered_map<uint, uint> outTime1;
			// std::vector<int> inTime1(g.getSize()+1,0);
			// std::vector<int> outTime1(g.getSize()+1,0);
			// std::vector<std::vector<string>> prefix1(g.getSize() + 1);
			// prefix1[cur1].push_back(string("0."));
			// unordered_std::set<int> repeatNode;
			bool repeatNode = false;
			bool neighborCheck = false;
			// std::vector<int> cache1;
			int time = 1;
			while (stackTop >= 0)
			{
				cur1 = __vecVisitNode[stackTop];
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
						if (!__vecVisitBool[c])
						{
							__vecVisitBool[c] = true;
							stackTop++;
							__vecVisitNode[stackTop] = c;
						}
					}
					continue;
				}
				// curPos++;
				// const std::vector<UI> &w = g.getWeight(cur1);

				const std::vector<Edge> &neigh = _graph[cur1];

				uint step = 1;
				double u = dsfmt_gv_genrand_open_open();
				// std::cout<<neigh[0].second<< std::flush;
				if (neigh.size() > 0) // && 1 > neigh[0].second && neigh[0].second > 0 && u > 0)
				{
					step = ceil(log(u) / log(1 - neigh[0].second));
				}
				else
				{
					continue;
				}
				int code = 0;

				for (i = step - 1; i < _graphinfo.node_deg[cur1 + 1]; i += step)
				{
					u = dsfmt_gv_genrand_open_open();
					if (1 > neigh[0].second && u > 0)
						step = ceil(log(u) / log(1 - (double)neigh[0].second));
					if (true)
					{ // Mark whether the edge is preserved or not by probability
						child[cur1].push_back(neigh[i].first);
						// prefix1[neigh[i]].push_back(prefix1[cur1][0] + to_string(code) + ".");
						if (!__vecVisitBool[neigh[i].first])
						{
							num_marked++;
							__vecVisitBool[neigh[i].first] = true;
							stackTop++;
							__vecVisitNode[stackTop] = neigh[i].first;
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
					if (std::find(__vecVisitNode.begin(), __vecVisitNode.begin() + stackTop + 1, v) != __vecVisitNode.begin() + stackTop + 1)
					{
						uAv = true;
					}
					if (std::find(__vecVisitNode.begin(), __vecVisitNode.begin() + stackTop + 1, cur3) != __vecVisitNode.begin() + stackTop + 1)
					{
						uAw = true;
					}
					if (uAv || uAw)
					{
						__vecVisitBool.clear();
						processed[u] = false;
						// __vecVisitNode.clear();
						if (uAv && uAw)
						{
							__vecVisitBool[v] = true;
							__vecVisitBool[cur3] = true;
							__vecVisitNode[0] = v;
							__vecVisitNode[1] = cur3;
							stackTop = 1;
						}
						else
						{
							__vecVisitBool[v] = uAv;
							__vecVisitBool[cur3] = uAw;
							if (uAv)
								__vecVisitNode[0] = v;
							if (uAw)
								__vecVisitNode[0] = cur3;
							stackTop = 0;
						}
					}
				}

			} //
			if (!stopUpper && forU == u)
			{
				std::set<int> rrsetU;
				{
					UI cur__ = forU;
					std::unordered_map<int, bool> visit_;
					std::vector<UI> __vecVisitNode_(num_marked);
					int num_marked_ = 0;
					int curPos_ = 0;
					// UI cur__ =v;

					UI cur_ = cur__;
					if (!visit_[cur_])
					{
						rrsetU.insert(cur_);
						num_marked_ += 1;
						visit_[cur_] = true;
						__vecVisitNode_[curPos_] = cur_;
						while (curPos_ < num_marked_)
						{
							cur_ = __vecVisitNode_[curPos_];
							curPos_++;
							const std::vector<uint> &neigh = child[cur_];
							// cout<<"\r"<<step<<" "<<u<<flush;
							for (i = 0; i < child[cur_].size(); i += 1)
							{
								if (true)
								{ // Mark whether the edge is preserved or not by probability
									if (!visit_[neigh[i]])
									{
										visit_[neigh[i]] = true;
										__vecVisitNode_[num_marked_] = neigh[i];
										rrsetU.insert(neigh[i]);
										num_marked_++; // BFS结构里节点访问完毕则停止
									}
								}
							}
						}
					}
				}
				auto result = std::vector<uint>(rrsetU.begin(), rrsetU.end());
				_RRsetsU.push_back(result);
				for (auto res : result)
				{
					_FRsetsU[res].push_back(hyperIdx);
				}
			}
			// //std::cout << "\rHere!" << 2 << std::flush;
			// int numrr_1 = num_marked;
			// std::vector<int> rrset1(__vecVisitNode.begin(), __vecVisitNode.begin() + numrr_1);
			std::unordered_map<uint, bool> visit1(__vecVisitBool);
			std::set<int> inter_mark;
			std::unordered_map<uint, uint> inTime2;
			std::unordered_map<uint, uint> outTime2;
			neighborCheck = false;
			// std::vector<int> inTime2(g.getSize()+1,0);
			// std::vector<int> outTime2(g.getSize()+1,0);
			if (!__vecVisitBool[cur2])
			{

				__vecVisitBool[cur2] = true;
				__vecVisitNode[0] = cur2;
				stackTop += 1;
				while (stackTop >= 0)
				{
					cur2 = __vecVisitNode[stackTop];
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
					// const std::vector<UI> &w = g.getWeight(cur2);
					const std::vector<Edge> &neigh = _graph[cur2];
					uint step = 1;
					double u = dsfmt_gv_genrand_open_open();
					if (neigh.size() > 0) // && 1 > neigh[0].second && neigh[0].second > 0 && u > 0)
						step = ceil(log(u) / log(1 - neigh[0].second));
					else
						continue;
					int code = 0;
					for (i = step - 1; i < _graphinfo.node_deg[cur2 + 1]; i += step)
					{
						u = dsfmt_gv_genrand_open_open();
						if (1 > neigh[0].second && u > 0)
							step = ceil(log(u) / log(1 - (double)neigh[0].second));
						if (true)
						{ // Mark whether the edge is preserved or not by probability
							child[cur2].push_back(neigh[i].first);
							if (visit1[neigh[i].first])
							{
								inter_mark.insert(neigh[i].first);
								stackTop++;
								__vecVisitNode[stackTop] = neigh[i].first;
							}
							if (!__vecVisitBool[neigh[i].first])
							{
								num_marked++;
								__vecVisitBool[neigh[i].first] = true;
								stackTop++;
								__vecVisitNode[stackTop] = neigh[i].first;
							}
							else
							{
								if (!visit1[neigh[i].first])
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
						if (std::find(__vecVisitNode.begin(), __vecVisitNode.begin() + stackTop + 1, cur3) != __vecVisitNode.begin() + stackTop + 1)
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
			if (!stopUpper && forU == v)
			{
				std::set<int> rrsetU;
				{
					UI cur__ = forU;
					std::unordered_map<int, bool> visit_;
					std::vector<UI> __vecVisitNode_(num_marked);
					int num_marked_ = 0;
					int curPos_ = 0;
					// UI cur__ =v;

					UI cur_ = cur__;
					if (!visit_[cur_])
					{
						rrsetU.insert(cur_);
						num_marked_ += 1;
						visit_[cur_] = true;
						__vecVisitNode_[curPos_] = cur_;
						while (curPos_ < num_marked_)
						{
							cur_ = __vecVisitNode_[curPos_];
							curPos_++;
							const std::vector<uint> &neigh = child[cur_];
							// cout<<"\r"<<step<<" "<<u<<flush;
							for (i = 0; i < child[cur_].size(); i += 1)
							{
								if (true)
								{ // Mark whether the edge is preserved or not by probability
									if (!visit_[neigh[i]])
									{
										visit_[neigh[i]] = true;
										__vecVisitNode_[num_marked_] = neigh[i];
										rrsetU.insert(neigh[i]);
										num_marked_++; // BFS结构里节点访问完毕则停止
									}
								}
							}
						}
					}
				}
				auto result = std::vector<uint>(rrsetU.begin(), rrsetU.end());
				_RRsetsU.push_back(result);
				for (auto res : result)
				{
					_FRsetsU[res].push_back(hyperIdx);
				}
			}
			if (inter_mark.empty())
			{
				// success += 1;
				__vecVisitBool.clear();
				// for (i = 0; i < num_marked; ++i)
				// {
				//         __vecVisitBool[__vecVisitNode[i]] = false;
				// }
				// nemptyset += 1;
				_RRsets.push_back(std::vector<uint>());
				if (!stopUpper && forU == w)
				{
					std::set<int> rrset;
					int cur1 = w;
					std::unordered_map<uint, bool> __vecVisitBool;
					std::unordered_map<uint, uint> __vecVisitNode;
					{
						std::unordered_map<uint, std::vector<uint>> child;
						// Sampling 3 nodes
						int num_marked = 3;
						int stackTop = 0;
						std::unordered_map<uint, bool> processed;
						__vecVisitBool[cur1] = true;
						__vecVisitNode[0] = cur1;
						RRset rrset;
						// std::vector<int> cache1;
						// int time = 1;
						rrset.push_back(cur1);
						while (stackTop >= 0)
						{
							cur1 = __vecVisitNode[stackTop];
							if (processed[cur1])
							{
								stackTop--;

								continue;
							}
							processed[cur1] = true;
							// curPos++;
							// const std::vector<UI> &w = g.getWeight(cur1);

							const std::vector<Edge> &neigh = _graph[cur1];

							uint step = 1;
							double u = dsfmt_gv_genrand_open_open();
							// std::cout<<neigh[0].second<< std::flush;
							if (neigh.size() > 0) // && 1 > neigh[0].second && neigh[0].second > 0 && u > 0)
							{
								step = ceil(log(u) / log(1 - neigh[0].second));
							}
							else
							{
								continue;
							}
							int code = 0;

							for (i = step - 1; i < _graphinfo.node_deg[cur1 + 1]; i += step)
							{
								u = dsfmt_gv_genrand_open_open();
								if (1 > neigh[0].second && u > 0)
									step = ceil(log(u) / log(1 - (double)neigh[0].second));
								if (true)
								{ // Mark whether the edge is preserved or not by probability
									child[cur1].push_back(neigh[i].first);
									// prefix1[neigh[i]].push_back(prefix1[cur1][0] + to_string(code) + ".");
									if (!__vecVisitBool[neigh[i].first])
									{
										// num_marked++;
										__vecVisitBool[neigh[i].first] = true;
										stackTop++;
										__vecVisitNode[stackTop] = neigh[i].first;
										rrset.push_back(neigh[i].first);
									}
								}
							}

						} //
					}
					auto result = std::vector<uint>(rrset.begin(), rrset.end());
					_RRsetsU.push_back(result);
					for (auto res : result)
					{
						_FRsetsU[res].push_back(hyperIdx);
					}
				}
				return;
			}

			// cout<<tmpnm<< " "<<rrset.size()<<endl;
			std::unordered_map<uint, bool> visit12(__vecVisitBool);
			std::set<int> inter_mark2;
			if (!__vecVisitBool[cur3])
			{
				__vecVisitBool[cur3] = true;
				__vecVisitNode[0] = cur3;
				stackTop += 1;
				while (stackTop >= 0)
				{
					cur3 = __vecVisitNode[stackTop];
					if (processed[cur3])
					{
						stackTop--;
						time++;
						continue;
					}
					processed[cur3] = true;
					// const std::vector<UI> &w = g.getWeight(cur3);
					const std::vector<Edge> &neigh = _graph[cur3];
					uint step = 1;
					double u = dsfmt_gv_genrand_open_open();
					if (neigh.size() > 0) // && 1 > neigh[0].second && neigh[0].second > 0 && u > 0)
						step = ceil(log(u) / log(1 - neigh[0].second));
					else
						continue;
					int code = 0;
					for (i = step - 1; i < _graphinfo.node_deg[cur3 + 1]; i += step)
					{
						u = dsfmt_gv_genrand_open_open();
						if (1 > neigh[0].second && u > 0)
							step = ceil(log(u) / log(1 - (double)neigh[0].second));
						if (true)
						{ // Mark whether the edge is preserved or not by probability
							child[cur3].push_back(neigh[i].first);
							if (visit12[neigh[i].first])
							{
								inter_mark2.insert(neigh[i].first);
							}
							if (!__vecVisitBool[neigh[i].first])
							{
								num_marked++;
								__vecVisitBool[neigh[i].first] = true;
								stackTop++;
								__vecVisitNode[stackTop] = neigh[i].first;
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
			if (!stopUpper && forU == w)
			{
				std::set<int> rrsetU;
				{
					UI cur__ = forU;
					std::unordered_map<int, bool> visit_;
					std::vector<UI> __vecVisitNode_(num_marked);
					int num_marked_ = 0;
					int curPos_ = 0;
					// UI cur__ =v;

					UI cur_ = cur__;
					if (!visit_[cur_])
					{
						rrsetU.insert(cur_);
						num_marked_ += 1;
						visit_[cur_] = true;
						__vecVisitNode_[curPos_] = cur_;
						while (curPos_ < num_marked_)
						{
							cur_ = __vecVisitNode_[curPos_];
							curPos_++;
							const std::vector<uint> &neigh = child[cur_];
							// cout<<"\r"<<step<<" "<<u<<flush;
							for (i = 0; i < child[cur_].size(); i += 1)
							{
								if (true)
								{ // Mark whether the edge is preserved or not by probability
									if (!visit_[neigh[i]])
									{
										visit_[neigh[i]] = true;
										__vecVisitNode_[num_marked_] = neigh[i];
										rrsetU.insert(neigh[i]);
										num_marked_++; // BFS结构里节点访问完毕则停止
									}
								}
							}
						}
					}
				}
				auto result = std::vector<uint>(rrsetU.begin(), rrsetU.end());
				_RRsetsU.push_back(result);
				for (auto res : result)
				{
					_FRsetsU[res].push_back(hyperIdx);
				}
			}
			if (inter_mark2.empty())
			{
				// success += 1;
				__vecVisitBool.clear();
				_RRsets.push_back(std::vector<uint>());
				// nemptyset += 1;
				return;
			}
			//
			// HTrIM

			std::set<int> rrset;
			std::set<int> rrset2, rrset3;
			// repeatNode.insert(0);
			//  repeatNode = true;

			if (repeatNode)
			{
				{
					std::unordered_map<int, bool> visit_;
					std::vector<UI> __vecVisitNode_(num_marked);
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
							__vecVisitNode_[curPos_] = cur_;
							while (curPos_ < num_marked_)
							{
								cur_ = __vecVisitNode_[curPos_];
								curPos_++;
								const std::vector<uint> &neigh = child[cur_];
								// cout<<"\r"<<step<<" "<<u<<flush;
								for (i = 0; i < child[cur_].size(); i += 1)
								{
									if (true)
									{ // Mark whether the edge is preserved or not by probability
										if (!visit_[neigh[i]])
										{
											visit_[neigh[i]] = true;
											__vecVisitNode_[num_marked_] = neigh[i];
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
					std::unordered_map<int, bool> visit_;
					std::vector<UI> __vecVisitNode_(num_marked);
					int num_marked_ = 1;
					int curPos_ = 0;
					// UI cur__ = w;
					for (UI cur__ : inter_mark2)
					{
						// cout<< curPos_<<num_marked_<<endl;
						UI cur_ = cur__;
						if (!visit_[cur_])
						{
							// if (rrset.std::find(cur_) != rrset.end())
							rrset2.insert(cur_);
							if (num_marked_ == curPos_)
								num_marked_ += 1;
							visit_[cur_] = true;
							__vecVisitNode_[curPos_] = cur_;
							while (curPos_ < num_marked_)
							{
								cur_ = __vecVisitNode_[curPos_];
								curPos_++;
								const std::vector<uint> &neigh = child[cur_];
								// cout<<"\r"<<step<<" "<<u<<flush;
								for (i = 0; i < child[cur_].size(); i += 1)
								{
									if (true)
									{ // Mark whether the edge is preserved or not by probability
										if (!visit_[neigh[i]])
										{
											visit_[neigh[i]] = true;
											__vecVisitNode_[num_marked_] = neigh[i];
											// if (rrset.std::find(cur_) != rrset.end())
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
				std::set<int> cand;
				for (auto im2 : inter_mark2)
					for (auto im : inter_mark)
					{
						if ((inTime1[im] >= inTime1[im2] && outTime1[im] <= outTime1[im2]) || (inTime2[im] >= inTime2[im2] && outTime2[im] <= outTime2[im2]))
						{
							// cout<<inTime1[im] <<" "<< inTime1[im2] <<" "<< outTime1[im] <<" " <<outTime1[im2]<<";";
							// cout<<inTime2[im] <<" "<< inTime2[im2] <<" "<< outTime2[im] <<" " <<outTime2[im2]<<endl;
							cand.insert(im); // im2是im的祖先
						}
						if ((inTime1[im2] * outTime1[im2] != 0) && (inTime1[im] <= inTime1[im2] && outTime1[im] >= outTime1[im2])) // or (inTime2[im] <= inTime2[im2] && outTime2[im] >= outTime2[im2]))
							cand.insert(im2);																					   // im是im2祖先
					}
				std::unordered_map<uint, bool> visit_;
				std::vector<UI> __vecVisitNode_(num_marked);
				int num_marked_ = 1;
				int curPos_ = 0;
				// UI cur__ = w;
				for (UI cur__ : cand)
				{
					UI cur_ = cur__;
					if (!visit_[cur_])
					{
						// if (rrset.std::find(cur_) != rrset.end())
						rrset3.insert(cur_);
						if (num_marked_ == curPos_)
							num_marked_ += 1;
						visit_[cur_] = true;
						__vecVisitNode_[curPos_] = cur_;
						while (curPos_ < num_marked_)
						{
							cur_ = __vecVisitNode_[curPos_];
							curPos_++;
							const std::vector<uint> &neigh = child[cur_];
							// cout<<"\r"<<step<<" "<<u<<flush;
							for (i = 0; i < child[cur_].size(); i += 1)
							{
								if (true)
								{ // Mark whether the edge is preserved or not by probability
									if (!visit_[neigh[i]])
									{
										visit_[neigh[i]] = true;
										__vecVisitNode_[num_marked_] = neigh[i];
										// if (rrset.std::find(cur_) != rrset.end())
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
			std::vector<uint> tmps1, tmps2, tmps3, result;
			// int tmpp1, tmpp2;
			if (!repeatNode)
			{
				// success += 1;
				tmps3 = std::vector<uint>(rrset3.begin(), rrset3.end());
				result = tmps3;
				// tmpp1 = result.size();
			}
			else
			{
				tmps1 = std::vector<uint>(rrset.begin(), rrset.end());
				tmps2 = std::vector<uint>(rrset2.begin(), rrset2.end());
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
			// _RRsets.push_back(std::vector<int>(__vecVisitNode.begin(),__vecVisitNode.begin()+num_marked)); //
			// for (int node : result)
			// {
			// 	if (link[node] < k)
			// 	{
			// 		t = true;
			// 		break;
			// 	}
			// }
			// for (i = 0; i < num_marked; ++i)
			// {
			//         __vecVisitBool[__vecVisitNode[i]] = false;
			// }
			// visit = std::vector<bool>(visit.size(), false);
			// if (result.empty())
			// 	nemptyset += 1;
			// //std::cout << "\rHere!" << 4 << std::flush;
			__vecVisitBool.clear();
			for (auto res : result)
			{
				_FRsets[res].push_back(hyperIdx);
			}
			_RRsets.push_back(result);
		}
		else if (_cascadeModel == LT)
		{
			// std::cout << "\rHere!" << 1 << std::flush;
			unsigned int num_marked = 0;
			unsigned int num_marked2 = 0;
			cur = cur1;
			for (i = 0; i < __numV; ++i)
			{
				if (__vecVisitBool[cur] == true)
					break;
				__vecVisitBool[cur] = true;
				__vecVisitNode[num_marked] = cur;
				num_marked++;
				num_marked2++;

				const std::vector<Edge> &neigh = _graph[cur];

				size_t ind;
				ind = gen_random_node_by_weight_LT(neigh);

				if (ind == neigh.size() + 1)
				{
					// _RRsets.push_back(std::vector<int>());
					break;
				}
				cur = neigh[ind].first;
			}
			if (!stopUpper && forU == cur1)
			{
				auto result = std::vector<uint>(__vecVisitNode.begin(), __vecVisitNode.begin() + num_marked);
				for (auto res : result)
				{
					_FRsetsU[res].push_back(hyperIdx);
				}
				_RRsetsU.push_back(result);
			}
			// std::cout << "\rHere!" << 2 << std::flush;
			std::unordered_map<uint, bool> visit1 = std::unordered_map<uint, bool>(__vecVisitBool);
			int intersect_flag = -1;
			cur = cur2;
			auto num_marked2_backup = num_marked2;
			for (i = 0; i < __numV; ++i)
			{
				if (visit1[cur] == true)
				{
					intersect_flag = find(__vecVisitNode.begin(), __vecVisitNode.begin() + num_marked, cur) - __vecVisitNode.begin();
					break;
				}
				if (__vecVisitBool[cur] == true)
				{
					_RRsets.push_back(std::vector<uint>());
					for (i = 0; i < num_marked2; ++i)
					{
						__vecVisitBool[__vecVisitNode[i]] = false;
					}
					if (!stopUpper && forU == cur2)
					{
						auto result = std::vector<uint>(__vecVisitNode.begin() + num_marked2_backup, __vecVisitNode.begin() + num_marked2);
						for (auto res : result)
						{
							_FRsetsU[res].push_back(hyperIdx);
						}
						_RRsetsU.push_back(result);
					}
					if (!stopUpper && forU == cur3)
					{
						RRset rrset;
						unsigned int num_marked = 0;
						unsigned int num_marked2 = 0;
						std::unordered_map<uint, bool> __vecVisitBool;
						std::unordered_map<uint, uint> __vecVisitNode;
						int cur = cur3;
						int ind;
						for (i = 0; i < __numV; ++i)
						{
							if (__vecVisitBool[cur] == true)
								break;
							__vecVisitBool[cur] = true;
							__vecVisitNode[num_marked] = cur;
							rrset.push_back(cur);
							num_marked++;
							num_marked2++;

							const std::vector<Edge> &neigh = _graph[cur];

							size_t ind;
							ind = gen_random_node_by_weight_LT(neigh);

							if (ind == neigh.size() + 1)
							{
								// _RRsets.push_back(std::vector<int>());
								break;
							}
							cur = neigh[ind].first;
						}
						// std::cout << "\rHere!" << 2 << std::flush;
						_RRsetsU.push_back(rrset);
						for (auto res : rrset)
						{
							_FRsetsU[res].push_back(hyperIdx);
						}
					}
					return;
				}
				__vecVisitBool[cur] = true;
				__vecVisitNode[num_marked2] = cur;
				num_marked2++;
				const std::vector<Edge> &neigh = _graph[cur];
				size_t ind;
				ind = gen_random_node_by_weight_LT(neigh);

				if (ind == neigh.size() + 1)
				{
					_RRsets.push_back(std::vector<uint>());
					for (i = 0; i < num_marked2; ++i)
					{
						__vecVisitBool[__vecVisitNode[i]] = false;
					}
					if (!stopUpper &&forU == cur2)
					{
						auto result = std::vector<uint>(__vecVisitNode.begin() + num_marked2_backup, __vecVisitNode.begin() + num_marked2 + 1);
						for (auto res : result)
						{
							_FRsetsU[res].push_back(hyperIdx);
						}
						_RRsetsU.push_back(result);
					}
					if (!stopUpper &&forU == cur3)
					{
						RRset rrset;
						unsigned int num_marked = 0;
						unsigned int num_marked2 = 0;
						std::unordered_map<uint, bool> __vecVisitBool;
						std::unordered_map<uint, uint> __vecVisitNode;
						int cur = cur3;
						int ind;
						for (i = 0; i < __numV; ++i)
						{
							if (__vecVisitBool[cur] == true)
								break;
							__vecVisitBool[cur] = true;
							__vecVisitNode[num_marked] = cur;
							rrset.push_back(cur);
							num_marked++;
							num_marked2++;

							const std::vector<Edge> &neigh = _graph[cur];

							size_t ind;
							ind = gen_random_node_by_weight_LT(neigh);

							if (ind == neigh.size() + 1)
							{
								// _RRsets.push_back(std::vector<int>());
								break;
							}
							cur = neigh[ind].first;
						}
						// std::cout << "\rHere!" << 2 << std::flush;
						_RRsetsU.push_back(rrset);
						for (auto res : rrset)
						{
							_FRsetsU[res].push_back(hyperIdx);
						}
					}
					return;
				}
				cur = neigh[ind].first;
			}
			if (!stopUpper &&forU == cur2)
			{
				auto tmp = std::vector<uint>(__vecVisitNode.begin() + num_marked2_backup, __vecVisitNode.begin() + num_marked2);
				auto tmp2 = std::vector<uint>(__vecVisitNode.begin() + intersect_flag, __vecVisitNode.begin() + num_marked);
				tmp.insert(tmp.end(), tmp2.begin(), tmp2.end());
				std::vector<uint> result = tmp;
				for (auto res : result)
				{
					_FRsetsU[res].push_back(hyperIdx);
				}
				_RRsetsU.push_back(result);
			}
			// std::cout << "\rHere!" << 3 << std::flush;
			cur = cur3;
			int intersect_flag2 = -1;
			num_marked2_backup = num_marked2;
			std::unordered_map<uint, bool> visit12 = std::unordered_map<uint, bool>(__vecVisitBool);
			for (i = 0; i < __numV; ++i)
			{
				if (visit12[cur] == true)
				{
					intersect_flag2 = find(__vecVisitNode.begin(), __vecVisitNode.begin() + num_marked, cur) - __vecVisitNode.begin();
					if (intersect_flag2 != num_marked && intersect_flag2 > intersect_flag)
						intersect_flag = intersect_flag2;
					break;
				}
				if (__vecVisitBool[cur] == true)
				{
					_RRsets.push_back(std::vector<uint>());
					for (i = 0; i < num_marked2; ++i)
					{
						__vecVisitBool[__vecVisitNode[i]] = false;
					}
					if (forU == cur3)
					{
						auto result = std::vector<uint>(__vecVisitNode.begin() + num_marked2_backup, __vecVisitNode.begin() + num_marked2);
						for (auto res : result)
						{
							_FRsetsU[res].push_back(hyperIdx);
						}
						_RRsetsU.push_back(result);
					}
					return;
				}
				__vecVisitBool[cur] = true;
				__vecVisitNode[num_marked2] = cur;
				num_marked2++;
				const std::vector<Edge> &neigh = _graph[cur];
				size_t ind;
				ind = gen_random_node_by_weight_LT(neigh);

				if (ind == neigh.size() + 1)
				{
					_RRsets.push_back(std::vector<uint>());
					for (i = 0; i < num_marked2; ++i)
					{
						__vecVisitBool[__vecVisitNode[i]] = false;
					}
					if (forU == cur3)
					{
						auto result = std::vector<uint>(__vecVisitNode.begin() + num_marked2_backup, __vecVisitNode.begin() + num_marked2 + 1);
						for (auto res : result)
						{
							_FRsetsU[res].push_back(hyperIdx);
						}
						_RRsetsU.push_back(result);
					}
					return;
				}
				cur = neigh[ind].first;
			}
			if (!stopUpper &&forU == cur3)
			{
				auto tmp = std::vector<uint>(__vecVisitNode.begin() + num_marked2_backup, __vecVisitNode.begin() + num_marked2);
				auto tmp2 = std::vector<uint>(__vecVisitNode.begin() + intersect_flag2, __vecVisitNode.begin() + num_marked);
				tmp.insert(tmp.end(), tmp2.begin(), tmp2.end());
				std::vector<uint> result = tmp;
				for (auto res : result)
				{
					_FRsetsU[res].push_back(hyperIdx);
				}
				_RRsetsU.push_back(result);
			}
			// std::cout << "\rHere!" << 4 << std::flush;
			//  cout<<num_marked-intersect_flag<<" "<<num_marked<<" "<<intersect_flag<<endl;
			auto result = std::vector<uint>(__vecVisitNode.begin() + intersect_flag, __vecVisitNode.begin() + num_marked);
			_RRsets.push_back(result);
			__vecVisitBool.clear();
			for (auto res : result)
			{
				_FRsets[res].push_back(hyperIdx);
			}
		}

		//
		//
		// size_t numVisitNode = 0, currIdx = 0;
		// _FRsets[uStart].push_back(hyperIdx);
		// __vecVisitNode[numVisitNode++] = uStart;
		// __vecVisitBool[uStart] = true;
		// while (currIdx < numVisitNode)
		// {
		// 	const auto expand = __vecVisitNode[currIdx++];
		// 	if (_cascadeModel == IC)
		// 	{
		// 		for (auto &nbr : _graph[expand])
		// 		{
		// 			const auto nbrId = nbr.first;
		// 			if (__vecVisitBool[nbrId])
		// 				continue;
		// 			const auto randDouble = dsfmt_gv_genrand_open_close();
		// 			if (randDouble > nbr.second)
		// 				continue;
		// 			__vecVisitNode[numVisitNode++] = nbrId;
		// 			__vecVisitBool[nbrId] = true;
		// 			_FRsets[nbrId].push_back(hyperIdx);
		// 		}
		// 	}
		// 	else if (_cascadeModel == LT)
		// 	{
		// 		if (_graph[expand].size() == 0)
		// 			continue;
		// 		const auto nextNbrIdx = gen_random_node_by_weight_LT(_graph[expand]);
		// 		if (nextNbrIdx >= _graph[expand].size())
		// 			break; // No element activated
		// 		const auto nbrId = _graph[expand][nextNbrIdx].first;
		// 		if (__vecVisitBool[nbrId])
		// 			break; // Stop, no further node activated
		// 		__vecVisitNode[numVisitNode++] = nbrId;
		// 		__vecVisitBool[nbrId] = true;
		// 		_FRsets[nbrId].push_back(hyperIdx);
		// 	}
		// }
		// for (int i = 0; i < numVisitNode; i++)
		// 	__vecVisitBool[__vecVisitNode[i]] = false;
		// _RRsets.push_back(RRset(__vecVisitNode.begin(), __vecVisitNode.begin() + numVisitNode));
	}

	/// Evaluate the influence spread of a seed set on current generated RR sets
	double self_inf_cal(const Nodelist &vecSeed)
	{
		std::vector<bool> vecBoolVst = std::vector<bool>(__numRRsets);
		std::vector<bool> vecBoolSeed(__numV);
		for (auto seed : vecSeed)
			vecBoolSeed[seed] = true;
		for (auto seed : vecSeed)
		{
			for (auto node : _FRsets[seed])
			{
				vecBoolVst[node] = true;
			}
		}
		return 1.0 * std::count(vecBoolVst.begin(), vecBoolVst.end(), true) * __numT / __numRRsets;
	}
	double self_inf_calU(const Nodelist &vecSeed)
	{
		std::vector<bool> vecBoolVst = std::vector<bool>(__numRRsets);
		std::vector<bool> vecBoolSeed(__numV);
		for (auto seed : vecSeed)
			vecBoolSeed[seed] = true;
		for (auto seed : vecSeed)
		{
			for (auto node : _FRsetsU[seed])
			{
				vecBoolVst[node] = true;
			}
		}
		return 1.0 * std::count(vecBoolVst.begin(), vecBoolVst.end(), true) * __numT / __numRRsets;
	}
	/// Efficiently estimate the influence spread with sampling error epsilon within probability 1-delta
	double effic_inf_valid_algo(const Nodelist &vecSeed, const double delta = 1e-3, const double eps = 0.01)
	{
		const double c = 2.0 * (exp(1.0) - 2.0);
		const double LambdaL = 1.0 + 2.0 * c * (1.0 + eps) * log(2.0 / delta) / (eps * eps);
		size_t numHyperEdge = 0;
		size_t numCoverd = 0;
		std::vector<bool> vecBoolSeed(__numV);
		for (auto seed : vecSeed)
			vecBoolSeed[seed] = true;

		while (numCoverd < LambdaL)
		{
			numHyperEdge++;
			size_t numVisitNode = 0, currIdx = 0;
			const auto uStart = dsfmt_gv_genrand_uint32_range(__numV);
			if (vecBoolSeed[uStart])
			{
				// Stop, this sample is covered
				numCoverd++;
				continue;
			}
			__vecVisitNode[numVisitNode++] = uStart;
			__vecVisitBool[uStart] = true;
			while (currIdx < numVisitNode)
			{
				const auto expand = __vecVisitNode[currIdx++];
				if (_cascadeModel == IC)
				{
					for (auto &nbr : _graph[expand])
					{
						const auto nbrId = nbr.first;
						if (__vecVisitBool[nbrId])
							continue;
						const auto randDouble = dsfmt_gv_genrand_open_close();
						if (randDouble > nbr.second)
							continue;
						if (vecBoolSeed[nbrId])
						{
							// Stop, this sample is covered
							numCoverd++;
							goto postProcess;
						}
						__vecVisitNode[numVisitNode++] = nbrId;
						__vecVisitBool[nbrId] = true;
					}
				}
				else if (_cascadeModel == LT)
				{
					if (_graph[expand].size() == 0)
						continue;
					const auto nextNbrIdx = gen_random_node_by_weight_LT(_graph[expand]);
					if (nextNbrIdx >= _graph[expand].size())
						break; // No element activated
					const auto nbrId = _graph[expand][nextNbrIdx].first;
					if (__vecVisitBool[nbrId])
						break; // Stop, no further node activated
					if (vecBoolSeed[nbrId])
					{
						// Stop, this sample is covered
						numCoverd++;
						goto postProcess;
					}
					__vecVisitNode[numVisitNode++] = nbrId;
					__vecVisitBool[nbrId] = true;
				}
			}
		postProcess:
			for (auto i = 0; i < numVisitNode; i++)
				__vecVisitBool[__vecVisitNode[i]] = false;
		}
		return 1.0 * numCoverd * __numV / numHyperEdge;
	}

	/// Efficiently evaluate the influence spread of a seed set with a given number of RR sets to test
	double effic_inf_valid_algo_with_samplesize(const std::vector<uint32_t> &vecSeed, const size_t numSamples)
	{
		size_t numHyperEdge = 0;
		size_t numCoverd = 0;
		std::vector<bool> vecBoolSeed(__numV);
		for (auto seed : vecSeed)
			vecBoolSeed[seed] = true;
		while (++numHyperEdge < numSamples)
		{
			size_t numVisitNode = 0, currIdx = 0;
			const auto uStart = dsfmt_gv_genrand_uint32_range(__numV);
			if (vecBoolSeed[uStart])
			{
				// Stop, this sample is covered
				numCoverd++;
				continue;
			}
			__vecVisitNode[numVisitNode++] = uStart;
			__vecVisitBool[uStart] = true;
			while (currIdx < numVisitNode)
			{
				const auto expand = __vecVisitNode[currIdx++];
				if (_cascadeModel == IC)
				{
					for (auto &nbr : _graph[expand])
					{
						const auto nbrId = nbr.first;
						if (__vecVisitBool[nbrId])
							continue;
						const auto randDouble = dsfmt_gv_genrand_open_close();
						if (randDouble > nbr.second)
							continue;
						if (vecBoolSeed[nbrId])
						{
							// Stop, this sample is covered
							numCoverd++;
							goto postProcess;
						}
						__vecVisitNode[numVisitNode++] = nbrId;
						__vecVisitBool[nbrId] = true;
					}
				}
				else if (_cascadeModel == LT)
				{
					if (_graph[expand].size() == 0)
						continue;
					const auto nextNbrIdx = gen_random_node_by_weight_LT(_graph[expand]);
					if (nextNbrIdx >= _graph[expand].size())
						break; // No element activated
					const auto nbrId = _graph[expand][nextNbrIdx].first;
					if (__vecVisitBool[nbrId])
						break; // Stop, no further node activated
					if (vecBoolSeed[nbrId])
					{
						// Stop, this sample is covered
						numCoverd++;
						goto postProcess;
					}
					__vecVisitNode[numVisitNode++] = nbrId;
					__vecVisitBool[nbrId] = true;
				}
			}
		postProcess:
			for (auto i = 0; i < numVisitNode; i++)
				__vecVisitBool[__vecVisitNode[i]] = false;
		}
		return 1.0 * numCoverd * __numV / numHyperEdge;
	}

	/// Refresh the hypergraph
	void refresh_hypergraph()
	{
		for (auto i = __numRRsets; i--;)
		{
			RRset().swap(_RRsets[i]);
		}
		RRsets().swap(_RRsets);
		for (auto i = __numV; i--;)
		{
			FRset().swap(_FRsets[i]);
		}
		__numRRsets = 0;
	}

	/// Release memory
	void release_memory()
	{
		refresh_hypergraph();
		// std::vector<bool>().swap(__vecVisitBool);
		std::unordered_map<uint, bool>().swap(__vecVisitBool);
		Nodelist().swap(__vecVisitNode);
		FRsets().swap(_FRsets);
	}
};

using THyperGraph = HyperGraph;
using PHyperGraph = std::shared_ptr<THyperGraph>;
