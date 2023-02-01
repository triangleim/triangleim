// This version fixes some of the problem in the SSA.cpp either for improvements or correctness
//
#include "option.h"
#include "hypergraph.hpp"
#include <iostream>
#include <ctime>
#include <cmath>

using namespace std;

double calculateInfluence(Graph &g, HyperGraph &hg, vector<int> &seeds, int t, int times, int m, float gamma, float delta)
{
    long long counter = 0;
    int n = g.getSize();
    long long nt = g.getTriNumber();
    unsigned int k = seeds.size();
    int maxSamples = times;

    double upperBound = 1+ (1+gamma)*4*(exp(1)-2)*log(2/delta)/(gamma*gamma);
    // If enable improved SSA bounds:
    // double gamma_2 = gamma / (3 * (1 - 1 / exp(1) - gamma));
    // double f = log2(nt * (log(2 / delta) + lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) / ((k / 3) * log(2 / delta))); // The IMM algorithm yields a bound here.
    // double upperBound = 1 + (2 + 2 * gamma_2 / 3) * (1 + gamma_2) * (log(3 / delta) + log(f)) / (gamma_2 * gamma_2); // Lambda_2
    vector<unsigned int> link(n + 1, k);
    int degree = 0;
    vector<float> benefits(n, 0);
    // cout << "Seeds: ";
    for (unsigned int i = 0; i < k; ++i)
    {
        link[seeds[i]] = i; // 目前种子集的index映射
                            //	cout << seeds[i] << " ";
    }
    // cout << endl;
    vector<bool> maxSeed(t, false);
    omp_set_num_threads(t);
#pragma omp parallel
    {
        vector<bool> visit(n + 1, false);
        vector<int> visit_mark(n, 0);
        int id = omp_get_thread_num(); // All algorithms in this work are designed for single-threaded use.
        if (m == 0)
        {
            // cout << endl;
            while (counter <= maxSamples && degree < upperBound)
            {

                // cout << "\r" << counter << "/" << maxSamples << " " << degree << "/";
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
            // cout << endl;
            while (counter <= maxSamples && degree < upperBound)
            {
                unordered_map<int, bool> visit;
                // cout << "\r" << counter << "/" << maxSamples << " " << degree << "/";

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
    
    cout << counter<<","<<upperBound<<","<<degree <<endl;
    return degree / (double)counter;
}
void readSeeds(const char *filename, vector<int> &seeds, int k)
{
    ifstream in(filename);
    int tmp;
    for (int i = 0; i < k; i++)
    {
        in >> tmp;
        seeds.push_back(tmp);
    }
    in.close();
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

    char *upperFile = op.getPara("-u");
    if (upperFile == NULL)
    {
        upperFile = (char *)"network.seeds";
    }

    char *lowerFile = op.getPara("-l");
    if (lowerFile == NULL)
    {
        lowerFile = (char *)"network.seeds";
    }
    char *finalFile = op.getPara("-f");
    if (finalFile == NULL)
    {
        finalFile = (char *)"network.seeds";
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
    float gamma = 0.1;
    tmp = op.getPara("-gamma");
    if (tmp != NULL)
    {
        gamma = atof(tmp);
    }
    int times = 0;
    tmp = op.getPara("-times");
    if (tmp != NULL)
    {
        times = atoi(tmp);
    }
    if (numTriangle > 0)
        g.numTriangles = numTriangle;
    // char* indexFile = op.getPara("-index");
    // if (tmp != NULL)
    // {
    // 	t = atoi(tmp);
    // }
    // char* nodesccFile= op.getPara("-nodescc");
    // if (tmp != NULL)
    // {
    // 	t = atoi(tmp);
    // }
    long long nt = g.getTriNumber();
    // g.readIndex(indexFile,nodesccFile);
    clock_t start = clock();
    double sigma, mu, nu;
    {
        HyperGraph hg(n);
        vector<double> degree(k + 1, 0);
        // Improved version SSA
        // Improved version SSA Bound

        vector<int> seeds;
        readSeeds(outFile, seeds, k);
        // cout << "Adding " << curSamples << " samples " << endl;
        long long int totalSamples = 0;

        int mo = 0;
        if (strcmp(model, "IC") == 0)
        {
            mo = 1;
        }

        sigma = calculateInfluence(g, hg, seeds, t, times, mo, gamma, delta);
    }
    {
        HyperGraph hg(n);
        vector<double> degree(k + 1, 0);
        // Improved version SSA
        // Improved version SSA Bound

        vector<int> seeds;
        readSeeds(lowerFile, seeds, k);
        // cout << "Adding " << curSamples << " samples " << endl;
        long long int totalSamples = 0;

        int mo = 0;
        if (strcmp(model, "IC") == 0)
        {
            mo = 1;
        }

        mu = calculateInfluence(g, hg, seeds, t, times, mo, gamma, delta);
    }
    {
        HyperGraph hg(n);
        vector<double> degree(k + 1, 0);
        // Improved version SSA
        // Improved version SSA Bound

        vector<int> seeds;
        readSeeds(upperFile, seeds, k);
        // cout << "Adding " << curSamples << " samples " << endl;
        long long int totalSamples = 0;

        int mo = 0;
        if (strcmp(model, "IC") == 0)
        {
            mo = 1;
        }

        nu = calculateInfluence(g, hg, seeds, t, times, mo, gamma, delta);
    }
    // cout << "\nTotal Samples: " << (totalSamples + checkSam) << endl;
    // cout << "Seed Nodes: ";
    cout << "Selection Time: " << (float)(clock() - start) / CLOCKS_PER_SEC << "s" << endl;
    // cout << "Memory: " << getCurrentMemoryUsage() << " MB" << endl;
    if (mu >= nu && mu >= sigma)
    {
        vector<int> seeds;
        readSeeds(lowerFile, seeds, k);
        cout << "Seed Nodes: ";
        ofstream out(finalFile);
        for (unsigned int i = 0; i < seeds.size(); ++i)
        {
            cout << seeds[i] << " ";
            out << seeds[i] << endl;
        }
        printf("Influence rate: %0.5lf\n", mu);
    }
    else if (nu >= mu && nu >= sigma)
    {
        vector<int> seeds;
        readSeeds(upperFile, seeds, k);
        cout << "Seed Nodes: ";
        ofstream out(finalFile);
        for (unsigned int i = 0; i < seeds.size(); ++i)
        {
            cout << seeds[i] << " ";
            out << seeds[i] << endl;
        }
        printf("Influence rate: %0.5lf\n", nu);
    }
    else
    {
        vector<int> seeds;
        readSeeds(outFile, seeds, k);
        cout << "Seed Nodes: ";
        ofstream out(finalFile);
        for (unsigned int i = 0; i < seeds.size(); ++i)
        {
            cout << seeds[i] << " ";
            out << seeds[i] << endl;
        }
        printf("Influence rate: %0.5lf\n", sigma);
    }
}
