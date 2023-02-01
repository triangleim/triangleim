// This version fixes some of the problem in the SSA.cpp either for improvements or correctness
//
#include "option.h"
#include "hypergraph.hpp"
#include <iostream>
#include <ctime>
#include <cmath>

using namespace std;

bool calculateInfluence(Graph &g, HyperGraph &hg, vector<int> &seeds, int t, double expected, double epsilon_1, float epsilon_2, float delta, int m, long long int maxSamples, long long int &checkSam)
{
	long long counter = 0;
	int n = g.getSize();
	long long nt = g.getTriNumber();
	unsigned int k = seeds.size();
	double f = log2(nt * (log(2 / delta) + lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) / ((k / 3) * log(2 / delta))); // The IMM algorithm yields a bound here.
	//cout << f << endl;
	double upperBound = 1 + (2 + 2 * epsilon_2 / 3) * (1 + epsilon_2) * (log(3 / delta) + log(f)) / (epsilon_2 * epsilon_2); // Lambda_2
	vector<unsigned int> link(n + 1, k);
	int degree = 0;
	vector<float> benefits(n, 0);
	//cout << "Seeds: ";
	for (unsigned int i = 0; i < k; ++i)
	{
		link[seeds[i]] = i; //目前种子集的index映射
							//	cout << seeds[i] << " ";
	}
	//cout << endl;
	vector<bool> maxSeed(t, false);
	omp_set_num_threads(t);
#pragma omp parallel
	{
		vector<bool> visit(n + 1, false);
		vector<int> visit_mark(n, 0);
		int id = omp_get_thread_num(); // All algorithms in this work are designed for single-threaded use.
		if (m == 0)
		{
			cout << endl;
			while (counter <= maxSamples && degree < upperBound)
			{
				if (counter % 10000 == 0)
					cout << "\r" << counter << "/" << maxSamples << " " << degree << "/" << upperBound;
				if ((1 + epsilon_1) * (upperBound / counter) < expected)
					break;
				maxSeed[id] = hg.pollingLT(g, link, k, visit, visit_mark);
#pragma omp critical
				{
					counter++;
					if (maxSeed[id])
					{
						degree++;
					}
				}
			}
		}
		else
		{
			cout << endl;
			while (counter <= maxSamples && degree < upperBound)
			{
				if (counter % 10000 == 0)
					cout << "\r" << counter << "/" << maxSamples << " " << degree << "/" << upperBound;
				if ((1 + epsilon_1) * (upperBound / counter) < expected)
					break;
				maxSeed[id] = hg.pollingIC(g, link, k, visit, visit_mark);
#pragma omp critical
				{
					counter++;
					if (maxSeed[id])
					{ // aka. |Rj intersect S| or 1
						degree++;
					}
				}
			}
		}
	}
	checkSam += counter;
	//	cout << "Degree: " << degree << " " << counter << endl;
	cout << "\n"
		 << expected << " " << (1 + epsilon_1) * ((double)degree / (double)counter) << endl;
	if (expected <= (1 + epsilon_1) * ((double)degree / (double)counter))
	{
		return true;
	}
	return false;
}

int main(int argc, char **argv)
{
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
	char* indexFile = op.getPara("-index");
	// if (tmp != NULL)
	// {
	// 	t = atoi(tmp);
	// }
	char* nodesccFile= op.getPara("-nodescc");
	// if (tmp != NULL)
	// {
	// 	t = atoi(tmp);
	// }
	long long nt = g.getTriNumber();
	if(indexFile&&nodesccFile)
		g.readIndex(indexFile,nodesccFile);
	HyperGraph hg(n);
	vector<double> degree(k + 1, 0);
	// Improved version SSA
	double epsilon_1 = epsilon / (2 * (1 - 1 / exp(1) - epsilon));
	double epsilon_2 = epsilon / (3 * (1 - 1 / exp(1) - epsilon));
	double epsilon_3 = (epsilon - (epsilon_1 - epsilon_2 - epsilon_1 * epsilon_2) * (1 - 1 / exp(1) - epsilon)) / (1 - 1 / exp(1));
	cout << "epsilons: " << epsilon_1 << " " << epsilon_2 << " " << epsilon_3 << endl;

	double f = log2(nt * (log(2 / delta) + lgammal(n + 1) - lgammal(k + 1) - lgammal(n - k + 1)) / ((k / 3) * log(2 / delta))); //lgamma(n+1):log(n!), Bound from IMM
	double maxsamples = 8 * nt * (1 - 1 / exp(1)) * (log(6 / delta) + lgammal(n + 1) - lgammal(k + 1) - lgammal(n - k + 1)) / ((k / 3) * epsilon * epsilon);
	// log[n(log(2/delta)+log(n!)-log(k!)-log((n-k+1)!) / (k*log(2/delta)))]
	double degreeRequired = (2 + 2 * epsilon_3 / 3) * (1 + epsilon_1) * (1 + epsilon_2) * log(3 * f / delta) / (epsilon_3 * epsilon_3); // Lambda_1
																																		//	cout << f << " " << degreeRequired << " " << k << " " << n << endl;
	vector<int> seeds;

	long long int curSamples = (long long)degreeRequired;
	// cout << "Adding " << curSamples << " samples " << endl;
	long long int totalSamples = 0;

	int mo = 0;
	if (strcmp(model, "IC") == 0)
	{
		mo = 1;
	}
	clock_t start = clock();
	long long int checkSam = 0;
	cout << "Max samples:"
		 << " " << maxsamples << endl;
	while (true)
	{
		// try{
		seeds.clear();
		cout << "Adding " << curSamples << " samples " << endl<< flush;
		addHyperedge(g, hg, t, curSamples, mo);
		curSamples *= 2;
		totalSamples = hg.getNumEdge();
		cout << hg.emptyrate() << endl << flush;
		buildSeedSet(hg, seeds, n, k, degree);
		cout << "Current seeds:";
		for (unsigned int i = 0; i < seeds.size(); ++i)
		{
			cout << seeds[i] << " ";
		}

		cout << "\n"
			 << degree[k] << " " << degreeRequired << " " << (double)degree[k] / (double)hg.getNumEdge() <<" "<<hg.success/(double)hg.allexpand<< endl;
		//continue;
		if (degree[k] < degreeRequired)
		{ //degree[k]是覆盖的RR set个数
			continue;
		}
		if (totalSamples > maxsamples)
			break;
		//	cout << "Total Samples: " << totalSamples << " " << degreeRequired << endl;

		if (calculateInfluence(g, hg, seeds, t, (double)degree[k] / (double)hg.getNumEdge(), epsilon_1, epsilon_2, delta, mo, totalSamples * 2 * (1 + epsilon_2) / (1 - epsilon_2) * (epsilon_3 * epsilon_3) / (epsilon_2 * epsilon_2), checkSam))
		{
			break;
		}
		// }catch (std::exception& e)
		// {
		// 	std::cerr << "Exception catched : " << e.what() << std::endl;
		// }
	}
	cout << endl;
	cout << "\nTotal Samples: " << (totalSamples + checkSam) << endl;
	cout << "Seed Nodes: ";
	ofstream out(outFile);
	for (unsigned int i = 0; i < seeds.size(); ++i)
	{
		cout << seeds[i] << " ";
		out << seeds[i] << endl;
	}
	out.close();
	cout << endl
		 << endl;
	printf("Influence rate: %0.5lf\n", (double)degree[k] / (double)hg.getNumEdge());
	cout << "Time: " << (float)(clock() - start) / CLOCKS_PER_SEC << "s" << endl;
	cout << "Memory: " << getCurrentMemoryUsage() << " MB" << endl;
	cout << "\nTotal Samples: " << (totalSamples + checkSam) << endl;
}
