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
    const size_t size = 10;
    srand(1);

    ManSegArray* arr = new ManSegBase<false>(size);
    ManSegArray* arr2 = new ManSegBase<false>(size);
    double* darr = new double[size];
    double* darr2 = new double[size];
    bool arbitraryPrecisionIncrease = false;

    for(size_t i = 0; i < size; ++i)
    {
        double val = (static_cast<double>(rand()) / RAND_MAX);

        darr[i] = val;
        arr->setPair(i, val);

        val += 2.0;

        arr2->setPair(i, val);
        darr2[i] = val;
    }

    //     head sum     pair sum   std double sum
    double hsum = 0.0, psum = 0.0, dsum = 0.0;

    std::cout << "starting sums" << std::endl;
    for(size_t i = 0; i < size; ++i)
        hsum += ((arr->read(i) + arr2->read(i)));

    arbitraryPrecisionIncrease = true;

    // TODO: swap precision -> true
    // not sure how you promote from ManSegArray<false> -> ManSegArray<true>
    if(arbitraryPrecisionIncrease)
    {
        arr = arr->updoot();
        arr2 = arr2->updoot();
    }

    for(size_t i = 0; i < size; ++i)
        psum += ((arr->read(i) + arr2->read(i)));

    for(size_t i = 0; i < size; ++i)
        dsum += (darr[i] + darr2[i]);

    std::cout << std::fixed << std::setprecision(16);
    std:: cout << "heads only =\t" << hsum << "        ";
    printBinary(hsum);
    std::cout << "heads+tails =\t" << psum << "        ";
    printBinary(psum);
    std::cout << "std double =\t" << dsum << "        ";
    printBinary(dsum);

    // delete arr, arr2;
    delete[] darr, darr2;

    return 0;
}