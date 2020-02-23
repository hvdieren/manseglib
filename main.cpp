#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>

#include <chrono>
#include <math.h>

#include "mantissaSegmentation.hpp"
// #include "mantissaSegmentation_s.hpp"
#include "quicksort.h"

using namespace ManSeg;

// print binary of double
void printBinary(double d)
{
    unsigned long l = *reinterpret_cast<unsigned long*>(&d);
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

// print binary of float
void printBinary(float f)
{
    unsigned int l = *reinterpret_cast<unsigned int*>(&f);
    for(int i = 31; i >= 0; --i)
    {
        std::cout << ((l >> i) & 1);
        if(i == 31 || i == 23)
            std::cout << " ";
    }
    std::cout << std::endl;
}

// print binary of int
void printBinary(int l)
{
    for(int i = 31; i >= 0; --i)
    {
        std::cout << ((l >> i) & 1);
        if(i % 8 == 0 && i > 0)
            std::cout << " ";
    }
}

int main(int argc, char** argv)
{
    std::cout << std::setprecision(16);

    using fArray = TwoSegArray<false>;
    using tArray = TwoSegArray<true>;

    const int size = 5;

    fArray f(size);
    tArray f1(size);

    // f.setPair(0, 0.1);
    // f.setPair(1, 1.4f);
    // f.setPair(2, 22);
    // f.setPair(3, 45L);
    // f.setPair(4, true);

    std::cout << "sizeof(std::uint_fast32) = " << sizeof(std::uint32_t)*8 << std::endl;

    return 0;
   
    f[0] = 0.38572985629756297361978;
    f[1] = 1.4f;
    f[2] = 22;
    f[3] = 45L;
    f[4] = true;
 
    f1[0] = f[0]; 
    f1[1] = f[1]; 
    f1[2] = f[2];
    f1[3] = f[3];
    f1[4] = f[4];

    for(int i = 0; i < size; ++i)
        std::cout << i << ":f -> " << f[i] << "\t\tf1 -> " << f1[i] << "\n";
    std::cout << "\n";

    return 0;
}
