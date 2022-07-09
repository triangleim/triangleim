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
		_FRsetsI = FRsets(__numV);
		_FRseqs.push_back(std::vector<std::vector<uint64_t>>(__numV));
		_FRseqs.push_back(std::vector<std::vector<uint64_t>>(__numV));
		_FRseqs.push_back(std::vector<std::vector<uint64_t>>(__numV));
		// __vecVisitBool = std::vector<bool>(__numV);
		__vecVisitNode = Nodelist(__numV);
	}

public:
	GraphInfo &_graphinfo;
	/// _graph: reverse graph
	Graph &_graph;
	RevGraph &_revgraph;
	/// _FRsets: forward cover sets, _FRsets[i] is the node sets that node i can reach
	FRsets _FRsets;	 // Forward to RR set
	FRsets _FRsetsI; // Forward to RRI set
	FRseqs _FRseqs;
	/// _RRsets: reverse cover sets, _RRsets[i] is the node set that can reach node i
	RRcol _RRcol;
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
			// std::cout << "\rHere!" << i << "/" << numSamples << std::flush;
			build_one_RRset(i);
		}
	}

	/// Generate one RR set
	void build_one_RRset(const size_t hyperIdx)
	{
		size_t i;
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

		if (_cascadeModel == IC)
		{
			// std::cout << "\rHereIC!" << 0 << std::flush;
			// Build DFS graph
			std::unordered_map<uint, std::vector<uint>> child;
			// Sampling 3 nodes
			int num_marked = 3;
			int stackTop = 0;
			std::unordered_map<uint, bool> processed;
			__vecVisitBool[cur1] = true;
			__vecVisitNode[0] = cur1;
			RRset rrset1, rrset2, rrset3;
			// std::vector<int> cache1;
			// int time = 1;
			rrset1.push_back(cur1);
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
							rrset1.push_back(neigh[i].first);
						}
					}
				}

			} //
			//
			// std::cout << "\rHere!" << 2 << std::flush;
			// int numrr_1 = num_marked;
			// std::vector<int> rrset1(__vecVisitNode.begin(), __vecVisitNode.begin() + numrr_1);
			std::unordered_map<uint, bool> visit1(__vecVisitBool);
			std::set<int> inter_mark;
			// std::vector<int> inTime2(g.getSize()+1,0);
			// std::vector<int> outTime2(g.getSize()+1,0);
			if (!__vecVisitBool[cur2])
			{

				__vecVisitBool[cur2] = true;
				__vecVisitNode[0] = cur2;
				rrset2.push_back(cur2);
				stackTop += 1;
				while (stackTop >= 0)
				{
					cur2 = __vecVisitNode[stackTop];
					if (processed[cur2])
					{
						stackTop--;
						continue;
					}
					processed[cur2] = true;
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
							}
							if (!__vecVisitBool[neigh[i].first])
							{
								// num_marked++;
								__vecVisitBool[neigh[i].first] = true;
								stackTop++;
								__vecVisitNode[stackTop] = neigh[i].first;
								rrset2.push_back(neigh[i].first);
							}
						}
					}
				}
			}
			else
			{
				inter_mark.insert(v);
			}
			size_t subsize = rrset1.size() + rrset2.size();

			{
				std::unordered_map<uint, bool> visit_;
				std::vector<UI> visit_mark_(subsize);
				int num_marked_ = 0;
				int curPos_ = 0;
				// UI cur__ =v;
				for (UI cur_ : inter_mark)
				{
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
							const std::vector<uint> &neigh = child[cur_];
							// cout<<"\r"<<step<<" "<<u<<flush;
							for (i = 0; i < child[cur_].size(); i += 1)
							{
								if (true)
								{ // 按概率标注边是否保留
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
				// std::cout<<"\r"<<subsize<<"\t"<<curPos_<<std::flush;
			}

			// cout<<tmpnm<< " "<<rrset.size()<<endl;
			// std::cout << "\rHere!" << 3 << std::flush;
			std::unordered_map<uint, bool> visit12(__vecVisitBool);
			std::set<int> inter_mark2;
			if (!__vecVisitBool[cur3])
			{
				rrset3.push_back(cur3);
				__vecVisitBool[cur3] = true;
				__vecVisitNode[0] = cur3;
				stackTop += 1;
				while (stackTop >= 0)
				{
					cur3 = __vecVisitNode[stackTop];
					if (processed[cur3])
					{
						stackTop--;
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
								// num_marked++;
								__vecVisitBool[neigh[i].first] = true;
								stackTop++;
								__vecVisitNode[stackTop] = neigh[i].first;
								rrset3.push_back(neigh[i].first);
							}
						}
					}
				}
			}
			else
			{
				inter_mark2.insert(w);
			}
			subsize += rrset3.size();

			{
				std::unordered_map<int, bool> visit_;
				std::vector<UI> visit_mark_(subsize);
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
							const std::vector<uint> &neigh = child[cur_];
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
										rrset3.push_back(neigh[i]);
										num_marked_++; // BFS结构里节点访问完毕则停止
									}
								}
							}
						}
					}
				}
			}

			// HTrIM

			RRseq result;
			result.push_back(rrset1);
			result.push_back(rrset2);
			result.push_back(rrset3);
			// //std::cout << "\rHere!" << 4 << std::flush;
			__vecVisitBool.clear();
			// for (auto res : result)
			// {
			// 	_FRsets[res].push_back(hyperIdx);
			// }
			_RRcol.push_back(result);
		}
		else if (_cascadeModel == LT)
		{
			// std::cout << "\rHere!" << 1 << std::flush;
			RRset rrset1, rrset2, rrset3;
			unsigned int num_marked = 0;
			unsigned int num_marked2 = 0;
			cur = cur1;
			for (i = 0; i < __numV; ++i)
			{
				if (__vecVisitBool[cur] == true)
					break;
				__vecVisitBool[cur] = true;
				__vecVisitNode[num_marked] = cur;
				rrset1.push_back(cur);
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
			std::unordered_map<uint, bool> visit1 = std::unordered_map<uint, bool>(__vecVisitBool);
			int intersect_flag = -1;
			cur = cur2;
			for (i = 0; i < __numV; ++i)
			{
				if (visit1[cur] == true)
				{
					intersect_flag = find(__vecVisitNode.begin(), __vecVisitNode.begin() + num_marked, cur) - __vecVisitNode.begin();
					break;
				}
				if (__vecVisitBool[cur] == true)
				{
					break;
				}
				rrset2.push_back(cur);
				__vecVisitBool[cur] = true;
				__vecVisitNode[num_marked2] = cur;
				num_marked2++;
				const std::vector<Edge> &neigh = _graph[cur];
				size_t ind;
				ind = gen_random_node_by_weight_LT(neigh);

				if (ind == neigh.size() + 1)
				{
					break;
				}
				cur = neigh[ind].first;
			}
			// std::cout << "\rHere!" << 3 << std::flush;
			if (intersect_flag >= 0)
				rrset2.insert(rrset2.end(), rrset1.begin() + intersect_flag, rrset1.end());
			cur = cur3;
			int intersect_flag2 = -1;
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
					break;
				}
				rrset3.push_back(cur);
				__vecVisitBool[cur] = true;
				__vecVisitNode[num_marked2] = cur;
				num_marked2++;
				const std::vector<Edge> &neigh = _graph[cur];
				size_t ind;
				ind = gen_random_node_by_weight_LT(neigh);

				if (ind == neigh.size() + 1)
				{
					break;
				}
				cur = neigh[ind].first;
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
			// std::cout << "\rHere!" << 4 << std::flush;
			//  cout<<num_marked-intersect_flag<<" "<<num_marked<<" "<<intersect_flag<<endl;
			// auto result = std::vector<uint>(__vecVisitNode.begin() + intersect_flag, __vecVisitNode.begin() + num_marked);
			RRseq result;
			result.push_back(rrset1);
			result.push_back(rrset2);
			result.push_back(rrset3);
			_RRcol.push_back(result);
			__vecVisitBool.clear();
			// for (auto res : result)
			// {
			// 	_FRsets[res].push_back(hyperIdx);
			// }
		}
		// std::cout << "\rJustWantoUpdate" << std::flush;
		i = hyperIdx;
		for (unsigned int k = 0; k < 3; ++k)
		{
			unsigned int num2 = _RRcol[i][k].size();

			for (unsigned int j = 0; j < num2; ++j)
			{
				// cout<<"\r";
				if (_FRsets[_RRcol[i][k][j]].empty() || _FRsets[_RRcol[i][k][j]][_FRsets[_RRcol[i][k][j]].size() - 1] != i)
					_FRsets[_RRcol[i][k][j]].push_back(i); //这个node出现在了哪些的rrset0中
				// std::cout<<_FRseqs[k].size()<<"\t"<<_RRcol[i][k][j]<<std::endl;
				_FRseqs[k][_RRcol[i][k][j]].push_back(i); //这个node出现在了哪些的rrset的哪一个中
			}
		}

		unsigned int num3 = _RRcol[i][2].size();
		for (unsigned int j = 0; j < num3; ++j)
		{
			if (find(_RRcol[i][0].begin(), _RRcol[i][0].end(), _RRcol[i][2][j]) != _RRcol[i][0].end() &&
				find(_RRcol[i][1].begin(), _RRcol[i][1].end(), _RRcol[i][2][j]) != _RRcol[i][1].end())

				_FRsetsI[_RRcol[i][2][j]].push_back(i); //这个node出现在了哪些的rrset的交集中
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
		std::vector<bool> vecBoolVst = std::vector<bool>(__numRRsets, false);
		std::vector<bool> vecBoolVst0 = std::vector<bool>(__numRRsets, false);
		std::vector<bool> vecBoolVst1 = std::vector<bool>(__numRRsets, false);
		std::vector<bool> vecBoolVst2 = std::vector<bool>(__numRRsets, false);
		// std::vector<bool> vecBoolSeed(__numV);
		// for (auto seed : vecSeed)
		// 	vecBoolSeed[seed] = true;
		for (auto seed : vecSeed)
		{
			for (auto node : _FRseqs[0][seed])
			{
				vecBoolVst0[node] = true;
			}
			for (auto node : _FRseqs[1][seed])
			{
				vecBoolVst1[node] = true;
			}
			for (auto node : _FRseqs[2][seed])
			{
				vecBoolVst2[node] = true;
			}
		}
		for( size_t i =0;i<__numRRsets;i++){
			vecBoolVst[i]=vecBoolVst0[i]&&vecBoolVst1[i]&&vecBoolVst2[i];
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
			RRseq().swap(_RRcol[i]);
		}
		RRcol().swap(_RRcol);
		for (auto i = __numV; i--;)
		{
			FRset().swap(_FRsets[i]);
			FRset().swap(_FRsetsI[i]);
		}
		_FRseqs.clear();
		_FRseqs.push_back(std::vector<std::vector<uint64_t>>(__numV));
		_FRseqs.push_back(std::vector<std::vector<uint64_t>>(__numV));
		_FRseqs.push_back(std::vector<std::vector<uint64_t>>(__numV));
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
