#include <iostream>
#include <iomanip>
#include <random>

#include "mantissaSegmentation.hpp"

int main()
{
    int size = 10000;
    ManSegArray* arr = new ManSegArray(size);
    double* darr = new double[size];

    srand(1);
    for(int i = 0; i < size; i++)
    {
        double val = (static_cast<double>(rand()) / RAND_MAX);
        arr->setPair(i, val);
        darr[i] = val;
    }

    std::cout << std::setprecision(16);

    double sum = 0.0;
    for(int i = 0; i < size; i++)
        sum += arr->readHead(i);

    std::cout << "heads only =\t" << sum << std::endl;

    sum = 0.0;
    for(int i = 0; i < size; i++)
        sum += arr->readPair(i);

    std::cout << "heads + tails =\t" << sum << std::endl;

    sum = 0.0;
    for(int i = 0; i < size; i++)
        sum += darr[i];

    std::cout << "std double =\t" << sum << std::endl;

    return 0;
}