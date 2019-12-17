#include <iostream>
#include <iomanip>
#include <random>

#include "mantissaSegmentation.hpp"

void printBinary(double& d)
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

template<class ManSegBase>
void f(ManSegBase arr, ManSegBase arr2, double* darr, double* darr2, size_t size)
{
    //     mseg sum    std double sum
    double msum = 0.0, dsum = 0.0;

    // std::cout << "starting sums" << std::endl;
    for(size_t i = 0; i < size; ++i)
        msum += (arr[i] + arr2[i]);

    std:: cout << "manseg =\t" << msum << "        ";
    printBinary(msum);

    for(size_t i = 0; i < size; ++i)
        dsum += (darr[i] + darr2[i]);

    std::cout << "std double =\t" << dsum << "        ";
    printBinary(dsum);
}

void test1()
{
    // make calculations into templated function taking ManSegBase as type parameter
    std::cout << std::fixed << std::setprecision(16);
    const size_t size = 10000000UL;

    ManSegBase<false> arr(size);
    ManSegBase<false> arr2(size);
    // setup whatever precisions you need/want
    ManSegBase<true> tarr = arr.updoot();
    ManSegBase<true> tarr2 = arr2.updoot();

    double* darr = new double[size];
    double* darr2 = new double[size];

    srand(1);
    for(size_t i = 0; i < size; ++i)
    {
        double val = (static_cast<double>(rand()) / RAND_MAX);
        darr[i] = val;
        arr.setPair(i, val);

        val += 2.0;

        arr2.setPair(i, val);
        darr2[i] = val;
    }

    std::cout << "heads only" << std::endl;
    f(arr, arr2, darr, darr2, size);

    std::cout << "\npairs" << std::endl;
    f(tarr, tarr2, darr, darr2, size);

    delete[] darr, darr2;
}

void test2()
{
    std::cout << std::fixed << std::setprecision(16);
    // const size_t size = 10000000UL;
    const size_t size = 10UL;

    ManSegBase<true> arr(size);
    ManSegBase<true> arr2(size);
    ManSegBase<true> arr3(size);
    double d[10];
    double d2[10];
    double d3[10];

    srand(1);
    for(size_t i = 0; i < size; ++i)
    {
        double val = (static_cast<double>(rand()) / RAND_MAX);
        arr.setPair(i, val);
        d[i] = val;

        val += 2.0;

        arr2.setPair(i, val);
        d2[i] = val;
    }

    for(size_t i = 0; i < size; ++i)
    {
        double val = (static_cast<double>(rand()) / RAND_MAX);
        arr3[i] = arr[i];
        d3[i] = d[i];
    }

    for(size_t i = 0; i < size; ++i)
    {
        std::cout << "arr[" << i << "]=" << arr[i] << std::endl;
        std::cout << "arr2[" << i << "]=" << arr2[i] << std::endl;
        std::cout << "arr3[" << i << "]=" << arr3[i] << std::endl;
    }
    std::cout << "\n\n";
    for(size_t i = 0; i < size; ++i)
    {
        std::cout << "d[" << i << "]=" << d[i] << std::endl;
        std::cout << "d2[" << i << "]=" << d2[i] << std::endl;
        std::cout << "d3[" << i << "]=" << d3[i] << std::endl;
    }
}

int main()
{
    // test1(); // comparing manseg sum to std double sum & manseg conversion

    test2(); // checking array assignments

    return 0;
}

// TODO:
// need to fix assignment operators for ManSegHead/ManSegPair
// when accessing values via [] they do not allow assignment properly
// i.e. arr[i] = arr2[i] without casting
