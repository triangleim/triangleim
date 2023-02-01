#include "option.h"
#include "hypergraph.hpp"
#include "sfmt/SFMT.h"
#include <numeric>
#include <iostream>
#include <ctime>
#include <cmath>
#include <signal.h>
using namespace std;


bool app_stopped = false;

void sigint_handler(int sig)
{
	if (sig == SIGINT)
	{
		// ctrl+c
		std::cout << "ctrl+c pressed!" << std::endl;
		app_stopped = true;
	}
}
vector<long long> heuristics;
bool compare(int a, int b)
{
	return heuristics[a] > heuristics[b];
}
vector<int> intersection2(vector<int> &v1,
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
int main(int argc, char **argv)
{
	// signal(SIGINT, sigint_handler);
	srand(time(NULL));

	OptionParser op(argc, argv);
	if (!op.validCheck())
	{
		printf("Parameters error, please check the readme.txt file for correct format!\n");
		return -1;
	}
	char *inFile = op.getPara("-i");
	if (inFile == NULL)
	{
		inFile = (char *)"network.bin";
	}

	char *outFile = op.getPara("-o");
	if (outFile == NULL)
	{
		outFile = (char *)"network.seeds";
	}

	char *model = op.getPara("-m");
	if (model == NULL)
		model = (char *)"IC";

	Graph g;
	if (strcmp(model, "LT") == 0)
	{
		g.readGraphLT(inFile);
	}
	else if (strcmp(model, "IC") == 0)
	{
		g.readGraphIC(inFile);
	}
	else
	{
		printf("Incorrect model option!");
		return -1;
	}

	int n = g.getSize();

	char *tmp = op.getPara("-epsilon");
	float epsilon = 0.1;
	if (tmp != NULL)
	{
		epsilon = atof(tmp);
	}
	tmp = op.getPara("-triangle");
	long long numTriangle = 0;
	if (tmp != NULL)
	{
		numTriangle = atoll(tmp);
	}
	float delta = 1.0 / n;
	tmp = op.getPara("-delta");
	if (tmp != NULL)
	{
		delta = atof(tmp);
	}

	float k = 1;

	tmp = op.getPara("-k");
	if (tmp != NULL)
	{
		k = atof(tmp);
	}

	int t = 1;
	tmp = op.getPara("-t");
	if (tmp != NULL)
	{
		t = atoi(tmp);
	}
	if (numTriangle > 0)
		g.numTriangles = numTriangle;
	long long nt = g.getTriNumber();
	HyperGraph hg(n);
	vector<double> degree(k + 1, 0);

	vector<int> seeds;
	//.?



	int mo = 0;
	if (strcmp(model, "IC") == 0)
		mo = 1;

	int iter = 1;


	sfmt_t sfmtSeed;
	sfmt_init_gen_rand(&sfmtSeed, rand());
	//cache convert
	unordered_map<int,unsigned long long> nt_node;
	long long tmpp;
	int a,b;
	unordered_map<int, unordered_map<int,unsigned long long>> nt_edge;
	for (int i = 1; i < g.edgeList.size() + 1; i++)
	{
		// pair<int, int> tmp;
		if(g.edgeList[i].size()<2)
			continue;
		a = g.edgeList[i][0];
		b = g.edgeList[i][1];
		tmpp =g.weight_sampling[i]-g.weight_sampling[i-1];
		nt_edge[a][b]=tmpp;
		nt_node[a] += tmpp;
		nt_node[b] += tmpp;
	}
	cout<<g.weight_sampling.size()<<" "<<g.edgeList.size()<<endl;
	cout<<"total triangle:"<<g.weight_sampling.back()/3<<endl;

	clock_t start = clock();
	heuristics = vector<long long>(n + 1, 0);
	for (int i = 1; i < n + 1; i++)
	{
		cout << "\r" << i << "/" << n <<flush;
		long long formtri = 0;
		vector<int> &neigh = g.adjOutList[i];
		vector<UI> weights(neigh.size());
		for (int d = 0; d < neigh.size(); d++)
		{
			weights[d] = (double)1 / g.getDegree(neigh[d]);
		}
		for (int c = 0; c < neigh.size(); c++)
		{
			if (sfmt_genrand_uint32(&sfmtSeed) / (double)4294967295U < weights[c])
			{
				formtri+=nt_edge[i][neigh[c]];
				formtri+=nt_edge[neigh[c]][i];
			}
		}
		formtri += nt_node[i] / 2;
		heuristics[i] = formtri;
	}
	vector<int> tmps(n + 1, 0);
	for (int i = 1; i <= n; i++)
	{
		tmps[i] = i;
	}
	sort(tmps.begin(), tmps.end(), compare);
	for (unsigned int i = 0; i < k; ++i)
	{
		seeds.push_back(tmps[i]);
	}
	cout << "Seed Nodes: ";
	ofstream out(outFile);
	for (unsigned int i = 0; i < seeds.size(); ++i)
	{
		cout << seeds[i] << " ";
		out << seeds[i] << endl;
	}
	out.close();
	cout << endl;
	cout << "Time: " << (float)(clock() - start) / CLOCKS_PER_SEC << "s" << endl;
	cout << "Memory: " << getCurrentMemoryUsage() << " MB" << endl;
}
