#include <iostream>
#include <iomanip>
#include <chrono>

#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>

#include "mantissaSegmentation.hpp"

using namespace std;

constexpr double d = 0.85;
constexpr double tol = 1e-7;
constexpr int maxIter = 100;

class SparseMatrixCOO
{
public:
    int numVertices;
    int numEdges;

    SparseMatrixCOO(string file)
    {
        ifstream f;
        f.open(file, std::ios_base::in);
        if(!f.is_open())
        {
            cerr << "error opening file in COO" << endl;
            exit(1);
        }
        string line;
        getline(f, line);

        f >> numVertices;
        f >> numEdges;

        source = new int[numEdges];
        destination = new int[numEdges];

        for(int i = 0; i < numEdges; ++i)
            f >> source[i] >> destination[i];
    }

    void calculateOutDegree(int outdeg[])
    {
        for(int i = 0; i < numEdges; ++i)
            ++outdeg[source[i]];
    }

    template<class ManSegBase>
    void iterate(double d, ManSegBase* prevPr, ManSegBase* newPr, int outdeg[], ManSegBase* contr)
    {
        int src, dest;
        for(int i = 0; i < numEdges; ++i)
        {
            src = source[i];
            dest = destination[i];
            newPr[dest] += d*(prevPr[src]/outdeg[src]);
        }
    }


private:
    int* source;
    int* destination;
};

class SparseMatrixCSR
{
public:
    int numVertices;
    int numEdges;

    SparseMatrixCSR(string file)
    {
        ifstream f;
        f.open(file, std::ios_base::in);
        if(!f.is_open())
        {
            cerr << "error opening file in CSR" << endl;
            exit(1);
        }

        string line;
        getline(f, line);

        f >> numVertices;
        f >> numEdges;

        index = new int[numVertices + 1];
        dest = new int[numEdges];

        index[0] = 0;
        int pos = 0;

        getline(f, line);
        for(int i = 0; i < numVertices; ++i)
        {
            getline(f, line);
            istringstream iss(line);
            vector<string> results(istream_iterator<string>{iss}, istream_iterator<string>());

            pos += (results.size() - 1);
            index[i + 1] = pos;

            for(int j = 1; j < results.size(); ++j)
                dest[index[i] + (j - 1)] = stoi(results[j]);
        }
    }

    void calculateOutDegree(int outdeg[])
    {
        for(int i = 0; i < numVertices; ++i)
            outdeg[i] = (index[i + 1] - index[i]);
    }

    template<class ManSegBase>
    void iterate(double d, ManSegBase* prevPr, ManSegBase* newPr, int outdeg[], ManSegBase* contr)
    {
        for(int i = 0; i < numVertices; ++i)
            contr[i] = d*(prevPr[i]/outdeg[i]);

        int curr, next;
        for(int i = 0; i < numVertices; ++i)
        {
            curr = index[i];
            next = index[i + 1];
            for(int j = curr; j < next; ++j)
                newPr[dest[j]] += contr[i];
        }
    }

private:
    int* index;
    int* dest;
};

class SparseMatrixCSC
{
public:
    int numVertices;
    int numEdges;

    SparseMatrixCSC(string file)
    {
        ifstream f;
        f.open(file, std::ios_base::in);
        if(!f.is_open())
        {
            cerr << "error opening file in CSC" << endl;
            exit(1);
        }

        string line;
        getline(f, line);

        f >> numVertices;
        f >> numEdges;

        index = new int[numVertices + 1];
        source = new int[numEdges];

        index[0] = 0;
        int pos = 0;

        getline(f, line);
        for(int i = 0; i < numVertices; ++i)
        {
            getline(f, line);
            istringstream iss(line);
            vector<string> results(istream_iterator<string>{iss}, istream_iterator<string>());

            pos += (results.size() - 1);
            index[i + 1] = pos;

            for(int j = 1; j < results.size(); ++j)
                source[index[i] + (j - 1)] = stoi(results[j]);
        }
    }

    void calculateOutDegree(int outdeg[])
    {
        for(int i = 0; i < numEdges; ++i)
            ++outdeg[source[i]];
    }

    template<class ManSegBase>
    void iterate(double d, ManSegBase* prevPr, ManSegBase* newPr, int outdeg[], ManSegBase* contr)
    {
        for(int i = 0; i < numVertices; ++i)
            contr[i] = d*(prevPr[i]/outdeg[i]);

        int curr, next;
        for(int i = 0; i < numVertices; ++i)
        {
            curr = index[i];
            next = index[i + 1];
            for(int j = curr; j < next; ++j)
                newPr[i] += contr[source[j]];
        }
    }

private:
    int* index;
    int* source;
};

template<class ManSegBase>
double sum(ManSegBase* a, int& n)
{
    double d = 0.0;
    double err = 0.0;
    for(int i = 0; i < n; ++i)
    {
        // does d += a[i] with high accuracy
        double temp = d;
        double y = a[i] + err;
        d = temp + y;
        err = temp - d;
        err += y;
    }
    return d;
}

template<class ManSegBase>
double normDiff(ManSegBase* a, ManSegBase* b, int& n)
{
    double d = 0.0;
    double err = 0.0;
    for(int i = 0; i < n; ++i)
    {
        // does d += abs(b[i] - a[i]) with high accuracy
        double temp = d;
        double y = abs(b[i] - a[i]) + err;
        d = temp + y;
        err = temp - d;
        err += y;
    }
    return d;
}

template<typename SparseMatrix>
int pageRank(SparseMatrix* matrix, std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::nanoseconds>& totalSt)
{
    auto tmStart = chrono::high_resolution_clock::now();

    int n = matrix->numVertices;
    ManSegBase<false> x_f(n); // pagerank
    ManSegBase<true> x_t = *x_f.updoot();
    
    ManSegBase<false> v_f(n);
    ManSegBase<true> v_t = *v_f.updoot();
    
    ManSegBase<false> y_f(n); // new pagerank
    ManSegBase<true> y_t = *y_f.updoot();
    int* outdeg = new int[n];
    ManSegBase<false> contr_f(n); // contribution for each vertex
    ManSegBase<true> contr_t = *contr_f.updoot();

    double delta = 2.0;
    int iter = 0;

    for(int i = 0; i < n; ++i)
    {
        x_f.setPair(i, 1.0 / (double)(n));
        v_f.setPair(i, 1.0 / (double)(n));
        y_f.setPair(i, 0.0);
        outdeg[i] = 0;
    }

    matrix->calculateOutDegree(outdeg);

    auto tmInit = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
    cout << "Initialisation: " << tmInit << " seconds" << endl;
    auto tmTotal = 0.0;
    tmStart = chrono::high_resolution_clock::now();

    auto iterateS = chrono::high_resolution_clock::now();
    while(iter < maxIter && delta > tol)
    {
        matrix->iterate(d, x_t, y_t, outdeg, contr_t);

        double w = 1.0 - sum(y_t, n);
        for(int i = 0; i < n; ++i)
            y_t[i] += w * v_t[i];

        delta = normDiff(&x_t, &y_t, n);
        ++iter;

        for(int i = 0; i < n; ++i)
        {
            x_t[i] = static_cast<double>(y_t[i]);
            y_t[i] = 0.0;
        }

        auto tmStep = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
        cout << "iteration " << iter << ": delta=" << delta << " xnorm=" << sum(&x_t, n) 
            << " time=" << tmStep << " seconds" << endl;

        tmStart = chrono::high_resolution_clock::now();
    }
    auto endT = chrono::high_resolution_clock::now();

    auto iterateT = chrono::duration_cast<chrono::nanoseconds>(endT - iterateS).count()*1e-9;
    auto totalT = chrono::duration_cast<chrono::nanoseconds>(endT - totalSt).count()*1e-9;

    cout << "\nTotal time:" << totalT << " seconds\n"
        << "Total iterate time:" << iterateT << " seconds\n"
        << "Total seq time:" << (totalT - iterateT) << " seconds" << endl;

    if(delta > tol)
        cerr << "error: solution has not converged" << endl;

    for(int i = 0; i < n; ++i)
        cerr << i << " " << x[i] << endl;

    return 0;
}

int main(int argc, char** argv)
{
    cout << setprecision(16);
    cerr << setprecision(16);

    if(argc < 3)
    {
        cerr << "usage: ./pagerank <format> <input_file>" << endl;
        return 1;
    }

    string format = argv[1];
    string inputFile = argv[2];

    cout << "Format: " << format 
    << "\nInput file: " << inputFile
    << endl;

    auto tmStart = chrono::high_resolution_clock::now();
    auto totalSt = tmStart;

    if(format.compare("CSR") == 0)
    {
        SparseMatrixCSR* matrix = new SparseMatrixCSR(inputFile);
        auto tmInput = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
        cout << "Reading input: " << tmInput << " seconds" << endl;

        return pageRank(matrix, totalSt);
    }
    else if(format.compare("CSC") == 0)
    {
        SparseMatrixCSC* matrix = new SparseMatrixCSC(inputFile);
        auto tmInput = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
        cout << "Reading input: " << tmInput << " seconds" << endl;

        return pageRank(matrix, totalSt);
    }
    else if(format.compare("COO") == 0)
    {
        SparseMatrixCOO* matrix = new SparseMatrixCOO(inputFile);
        auto tmInput = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
        cout << "Reading input: " << tmInput << " seconds" << endl;

        return pageRank(matrix, totalSt);
    }
    else
    {
        cerr << "Unknown format: " << format << endl;
        return 1;
    }

    return 0;
}
