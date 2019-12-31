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
    // the destructor is being called as arr is being passed by value
    // this deletes the references to the array
    f(arr, arr2, darr, darr2, size);

    std::cout << "heads only" << std::endl;
    f(arr, arr2, darr, darr2, size);

    // setup whatever precisions you need/want
    ManSegBase<true> tarr = arr.updoot();
    ManSegBase<true> tarr2 = arr2.updoot();

    std::cout << "\npairs" << std::endl;
    f(tarr, tarr2, darr, darr2, size);

    std::cout << "pairs" << std::endl;
    f(tarr, tarr2, darr, darr2, size);

    delete[] darr, darr2;

    arr.freeMemory();
    arr2.freeMemory();
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

    // copy values
    for(size_t i = 0; i < size; ++i)
    {
        double val = (static_cast<double>(rand()) / RAND_MAX);
        arr3[i] = arr[i];
        d3[i] = d[i];
    }

    for(size_t i = 0; i < size; ++i)
    {
        std::cout << "arr[" << i << "]=" << arr[i] << " == arr3[" << i << "]=" << arr3[i] << std::endl;
        std::cout << "arr2[" << i << "]=" << arr2[i] << std::endl;
        // std::cout << "arr3[" << i << "]=" << arr3[i] << std::endl;
    }
    std::cout << "\n\n";
    for(size_t i = 0; i < size; ++i)
    {
        std::cout << "d[" << i << "]=" << d[i] << " == d3[" << i << "]=" << d3[i] << std::endl;
        std::cout << "d2[" << i << "]=" << d2[i] << std::endl;
        // std::cout << "d3[" << i << "]=" << d3[i] << std::endl;
    }

    arr.freeMemory();
    arr2.freeMemory();
    arr3.freeMemory();
}

void test3()
{
    ManSegBase<false> x(5);
    ManSegBase<true>y = x.updoot();

    y[0] = 1.26;
    y.setPair(1, 1.26);
    y[2] = 2.5;
    x[3] = 1.22;
    y.set(4, 1.22);

    for(int i = 0; i < 5; ++i)
        std::cout << y[i] << std::endl;
}

void test4()
{
    ManSegBase<false> x(10);
    ManSegBase<true> y = x.updoot();

    double d[10];

    for(int i = 0; i < 10; ++i)
    {
        y[i] = (i+1);
        d[i] = (i+1);
        std::cout<<"y["<<i<<"]="<<y[i]<<std::endl;
    }
    std::cout << "\n\n";

    for(int i = 0; i < 10; ++i)
    {
        // +=
        d[i] += 0.85*(d[i]+d[0]/0.5);
        x[i] += 0.85*(x[i]+x[0]/0.5);
        // y[i] += 0.85*(y[i]+y[0]/0.5);
        
        // -=
        d[i] -= 0.85;
        x[i] -= 0.85;
        // y[i] -= 0.85;
        
        // *=
        d[i] *= 1.5;
        x[i] *= 1.5;
        // y[i] *= 1.5;
        
        // /=
        d[i] /= 2;
        x[i] /= 2;
        // y[i] /= 2;

        std::cout<<"d["<<i<<"]="<<d[i]<<std::endl;
        std::cout<<"x["<<i<<"]="<<x[i]<<std::endl<<std::endl;
        // std::cout<<"y["<<i<<"]="<<y[i]<<std::endl<<std::endl;

    }
}

int main()
{
    // test1(); // comparing manseg sum to std double sum & manseg conversion (static assigned)
    // test2(); // checking array to array assignments
    // test3(); // initialisation 
    test4();    // checking to other assignment works correctly (i.e. x[i] += 2 and such)

    return 0;
}

// TODO:
// Can probably rename ManSegBase to ManSegArray since there is no encapsuling thing anymore
// unless you want the wrapper thing you'd talked about which holds refs to low and high precision or whatever
// but that might not be so good if multiple precisions in future or whatever

// I think this still applies, but with pointer version of ManSegBase

// need to fix assignment operators for ManSegHead/ManSegPair
// when accessing values via [] they do not allow assignment properly
// i.e. arr[i] = arr2[i] without casting
