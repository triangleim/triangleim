#include "option.h"
#include "hypergraph.hpp"
#include "sfmt/SFMT.h"
#include <iostream>
#include <ctime>
#include <cmath>
#include <signal.h>
using namespace std;

bool calculateInfluence(HyperGraph &hg, Graph &g, vector<int> &seeds, int t, double &deg, float epsilon, float delta, int m, long long int maxSamples, int iter)
{
	long long counter = 0;
	long long nt = g.getTriNumber();
	int n = g.getSize();
	unsigned k = seeds.size();
	vector<unsigned int> link(n + 1, seeds.size());
	double f = (log(6 / delta) + lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) * nt / ((k / 3) * log(6 * log2(nt) / delta));
	double lambda1 = 1 + (1 + epsilon) * (2 + 2 * epsilon / 3) * log(3 * log2(f) / delta) / (epsilon * epsilon);
	double degree = 0;
	for (unsigned int i = 0; i < k; ++i)
	{
		link[seeds[i]] = i;
	}
	vector<bool> maxSeed(t, false);

	omp_set_num_threads(t);
#pragma omp parallel
	{
		vector<bool> visit(n + 1, false);
		vector<int> visit_mark(n, 0);
		int id = omp_get_thread_num();

		if (m == 0)
		{
			while (counter < maxSamples)
			{
				if (counter % 10000 == 0)
					cout << "\r" << counter << "/" << maxSamples << " " << degree << "/" << lambda1 << " ";
				maxSeed[id] = hg.pollingLT2(g, link, k, visit, visit_mark);
#pragma omp critical
				{
					counter += 1;
					if (maxSeed[id])
					{
						degree++;
					}
				}
			}
		}
		else
		{
			unordered_map<int,bool> visit;
			while (counter < maxSamples)
			{
				if (counter % 10000 == 0)
					cout << "\r" << counter << "/" << maxSamples << " " << degree << "/" << lambda1;
				maxSeed[id] = hg.pollingIC2(g, link, k, visit, visit_mark);
#pragma omp critical
				{
					counter += 1;
					if (maxSeed[id])
					{
						degree++;
					}
				}
			}
		}
	}
	//cout << "Degree: " << degree << "/" << lambda1<< endl;

	if (degree >= lambda1)
	{
		double epsilon_1 = (deg * nt / maxSamples) / (degree * nt / counter) - 1;
		cout << "Epsilon_1 = " << epsilon_1 << "; ";
		double epsilon_2 = epsilon * sqrt(nt * (1 + epsilon) / (degree * nt * pow(2, iter - 1) / counter));
		cout << "Epsilon_2 = " << epsilon_2 << " " << epsilon * sqrt(nt * (1 + epsilon) / (degree * nt * pow(2, iter - 1) / counter)) << " " << pow(2, iter - 1) << " " << pow(3, iter - 1) << "; ";
		double epsilon_3 = epsilon * sqrt(nt * (1 + epsilon) * (1 - 1 / exp(1) - epsilon) / ((1 + epsilon / 3) * degree * nt * pow(2, iter - 1) / counter));
		cout << "Epsilon_3 = " << epsilon_3 << "; ";
		double epsilon_t =(epsilon_1 + epsilon_2 + epsilon_1 * epsilon_2) * (1 - 1 / exp(1) - epsilon) + epsilon_3 * (1 - 1 / exp(1));
		cout << "Epsilon_t = " << epsilon_t << "; ";
		cout << "Epsilon: " << epsilon << endl;
		if (epsilon_t <= epsilon)
		{
			return true;
		}
	}

	hg.updateDeg();
	return false;
}
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
int main(int argc, char **argv)
{
	signal(SIGINT, sigint_handler);
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
	double f = (log(6 / delta) + lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) * nt / ((k / 3) * log(6 * log2(nt) / delta));
	double maxsamples = 8 * nt * (1 - 1 / exp(1)) * (log(6 / delta) + lgammal(n + 1) - lgammal(k + 1) - lgammal(n - k + 1)) / ((k / 3) * epsilon * epsilon);
	//double f = (nt*(log(6/delta)+lgammal(n+1)-lgammal(k+1)-lgammal(n-k+1))/((k/3)*log(6/delta))); //lgamma(n+1):log(n!), Bound from IMM
	double lambda = (2 + 2 * epsilon / 3) * log(3 * log2(f) / delta) / (epsilon * epsilon);

	long long int totalSamples = (long long int)lambda;
	cout << lambda << " " << totalSamples << endl;

	int mo = 0;
	if (strcmp(model, "IC") == 0)
		mo = 1;

	int iter = 1;

	addHyperedge(g, hg, t, totalSamples, mo);
	double nmax = (2 + 2 * epsilon / 3) * (lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1) + log(6 / delta)) * nt / (epsilon * epsilon * k);
	cout << "Max samples:"
		 << " " << maxsamples << " " << nmax << endl;
	clock_t start = clock();
	cout << totalSamples << " " << nmax << " " << lgamma(n + 1) << " " << lgamma(k + 1) << " " << lgamma(n - k + 1) << endl;

	while (totalSamples < nmax)
	{
		seeds.clear();
		totalSamples = hg.getNumEdge();
		cout << "Total Samples: " << totalSamples << endl;
		buildSeedSet(hg, seeds, n, k, degree);
		if (calculateInfluence(hg, g, seeds, t, degree[k], epsilon, delta, mo, totalSamples, iter))
		{
			break;
		}
		cout << "\nEmptyrate:"
			 << " " << hg.emptyrate() << endl;
		iter++;
		if (app_stopped)
		{
			break;
		}
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
	printf("Influence rate: %0.5lf\n", (double)degree[k] / totalSamples);
	cout << "Time: " << (float)(clock() - start) / CLOCKS_PER_SEC << "s" << endl;
	cout << "Memory: " << getCurrentMemoryUsage() << " MB" << endl;
	cout << "Total Samples: " << totalSamples << endl;
}
