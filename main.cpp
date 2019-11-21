#include <iostream>
#include <iomanip>
#include <random>

#include "mantissaSegmentation.hpp"

void printBinary(double d) 
{ 
    unsigned long l = reinterpret_cast<unsigned long&>(d);
    for(int i = 63; i >= 0; --i)
    {
        std::cout << ((l >> i) & 1);
        if(i == 63 || i == 52)
            std::cout << " ";
        else if(i == 32)
            std::cout << "|";
    } 
    std::cout << std::endl; 
}

int main()
{
    srand(1);

    const size_t size = 1e8;
    ManSegArray* arr = new ManSegArray(size);
    double* darr = new double[size];

    for(size_t i = 0; i < size; i++)
    {
        double val = (static_cast<double>(rand()) / RAND_MAX);
        arr->set(i, val);
        darr[i] = val;
    }

    //     head sum     pair sum   std double sum
    double hsum = 0.0, psum = 0.0, dsum = 0.0;

    for(size_t i = 0; i < size; i++)
        hsum += (arr->read<ManSegHead>(i) * 0.25) / 3.0;

    for(size_t i = 0; i < size; i++)
        psum += (arr->read(i) * 0.25) / 3.0;

    for(size_t i = 0; i < size; i++)
        dsum += (darr[i] * 0.25) / 3.0;

    std::cout << std::fixed << std::setprecision(16);
    std:: cout << "heads only =\t" << hsum << "        ";
    printBinary(hsum);
    std::cout << "heads+tails =\t" << psum << "        ";
    printBinary(psum);
    std::cout << "std double =\t" << dsum << "        ";
    printBinary(dsum);

    delete arr;
    delete[] darr;

    return 0;
}