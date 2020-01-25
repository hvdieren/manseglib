#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>

#include <chrono>
#include <math.h>

#include "mantissaSegmentation.hpp"

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

// helper function for test1
template<class TwoSegmentArray>
void f(TwoSegmentArray arr, TwoSegmentArray arr2, double* darr, double* darr2, size_t size)
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

// test basic functionality
// test that heads and pairs operate correctly
void test1()
{
    // make calculations into templated function taking TwoSegmentArray as type parameter
    std::cout << std::fixed << std::setprecision(16);
    const size_t size = 10000000UL;

    TwoSegmentArray<false> arr(size);
    TwoSegmentArray<false> arr2(size);

    double* darr = new double[size];
    double* darr2 = new double[size];

    srand(1);
    for(size_t i = 0; i < size; ++i)
    {
        double val = (static_cast<double>(rand()) / RAND_MAX);
        darr[i] = val;
        arr.set(i, val);

        val += 2.0;

        arr2.set(i, val);
        darr2[i] = val;
    }

    std::cout << "heads only" << std::endl;
    // the destructor is being called as arr is being passed by value
    // this deletes the references to the array
    f(arr, arr2, darr, darr2, size);

    std::cout << "heads only" << std::endl;
    f(arr, arr2, darr, darr2, size);

    // setup whatever precisions you need/want
    TwoSegmentArray<true> tarr = arr.increasePrecision();
    TwoSegmentArray<true> tarr2 = arr2.increasePrecision();

    for(size_t i = 0; i < size; ++i)
    {
        tarr.set(i, darr[i]);
        tarr2.set(i, darr2[i]);
    }

    std::cout << "\npairs" << std::endl;
    f(tarr, tarr2, darr, darr2, size);

    std::cout << "pairs" << std::endl;
    f(tarr, tarr2, darr, darr2, size);

    delete[] darr, darr2;

    arr.del();
    arr2.del();
}

// test that assignment works as expected
// array to array assignment specifically
void test2()
{
    std::cout << std::fixed << std::setprecision(16);
    // const size_t size = 10000000UL;
    const size_t size = 10UL;

    TwoSegmentArray<true> arr(size);
    TwoSegmentArray<true> arr2(size);
    TwoSegmentArray<true> arr3(size);
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
        if((double)arr[i] != (double)arr3[i])
            std::cout << "arr[" << i << "]=" << arr[i] << " != arr3[" << i << "]=" << arr3[i] << std::endl;
    }
    for(size_t i = 0; i < size; ++i)
    {
        if(d2[i] != (double)arr2[i])
            std::cout << "d2[" << i << "]=" << d2[i] << " != arr2[" << i << "]=" << arr2[i] << std::endl;
    }

    arr.del();
    arr2.del();
    arr3.del();
}

// tests that other assignments work as expected
// i.e the same as set/setPair does
void test3()
{
    std::cout << std::setprecision(16) << std::endl;

    TwoSegmentArray<false> x(4);
    TwoSegmentArray<true>y = x.increasePrecision();

    y[0] = 1.26;
    y.setPair(1, 1.26);

    y[2] = 1.22;
    y.set(3, 1.22);

    if(y[0] != y[1])
        std::cout << "y[0] != y[1], values: y[0]=" << y[0] << " y[1]=" << y[1] << "\n";

    if(y[2] != y[3])
        std::cout << "y[2] != y[3], values: y[2]=" << y[2] << " y[3]=" << y[3] << "\n";


    x[0] = 1.26;
    x.set(1, 1.26);

    x[2] = 1.22;
    x.set(3, 1.22);

    if(x[0] != x[1])
        std::cout << "x[0] != x[1], values: x[0]=" << x[0] << " x[1]=" << x[1] << "\n";

    if(x[2] != x[3])
        std::cout << "x[2] != x[3], values: x[2]=" << x[2] << " x[3]=" << x[3] << "\n";

    x.del();
    y.del();
}

// helper function for test4
void fill(double* d, TwoSegmentArray<true>& y, int* v, const int& size)
{
    for(int i = 0; i < size; ++i)
    {
        y[i] = (i+1)/0.1;
        d[i] = (i+1)/0.1;
        v[i] = (i+1);
    }
}
constexpr int size = 10000;
constexpr double val = 0.85;
// tried to replicate what pagerank was doing
// in a small tester function for micro benchmarking
// not exactly the same, but somewhat close.
void test4(bool usePairs)
{
    TwoSegmentArray<false> x(size);                  // x = heads
    TwoSegmentArray<true> y = x.increasePrecision();            // y = pairs

    int* v = new int[size];
    double* d = new double[size];

    fill(d, y, v, size);

    // standard doubles
    auto dst = std::chrono::_V2::high_resolution_clock::now();

    for(int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j)
            d[i] += val*(d[size-1-j]/v[j]);

    auto ded = std::chrono::_V2::high_resolution_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::nanoseconds>(ded - dst).count()*1e-9;


    // using manseg
    if(!usePairs)
    {
        auto mxst = std::chrono::_V2::high_resolution_clock::now();

        for(int i = 0; i < size; ++i)
            for(int j = 0; j < size; ++j)
                x[i] += val*(x[size-1-j]/v[j]);

        auto mxed = std::chrono::_V2::high_resolution_clock::now();
        auto mt = std::chrono::duration_cast<std::chrono::nanoseconds>(mxed - mxst).count()*1e-9;

        std::cout << "std double time = " << dt << "s\n"
        << "heads time = " << mt << "s\n";
    }
    else
    {
        auto mxst = std::chrono::_V2::high_resolution_clock::now();

        for(int i = 0; i < size; ++i)
            for(int j = 0; j < size; ++j)
                y[i] += val*(y[size-1-j]/v[j]);

        auto mxed = std::chrono::_V2::high_resolution_clock::now();
        auto mt = std::chrono::duration_cast<std::chrono::nanoseconds>(mxed - mxst).count()*1e-9;

        std::cout << "std double time = " << dt << "s\n"
        << "pairs time = " << mt << "s\n";
    }

    delete[] v, d;
    y.del();
}

// reading SNAP file format test
void test5()
{
    using std::ifstream;
    using std::cout;
    using std::endl;

    ifstream file;
    file.open("../../graphs/Yahoo.txt");

    if(file.is_open())
    {
        for(int i = 0; i < 50; ++i)
        {
            std::string line;
            getline(file, line);
            cout << line << "\n";
        }

        file.close();
    }

    /*
        Outside of the initial couple of lines, looks similar to COO format.
        ****************************
        File format:
        #SNAP Yahoo
        #Nodes:123539879 Edges:123456789
        #FromNodeId\t   ToNodeId
        0 1
        0 12
        0 234
        1 25232
        *****************************
    */
}

// helper for test6 & test7
template<class TwoSegmentArray>
void t6Helper(TwoSegmentArray* ms, int size, int repeats)
{
    srand((unsigned int)(time(NULL)));
    for(int i = 0; i < repeats; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            double val = (static_cast<double>(rand()) / RAND_MAX);
            (*ms)[j] = val;
        }
    }
}
template<class TwoSegmentArray>
void t6Helper_v2(TwoSegmentArray* ms, int size, int repeats)
{
    srand((unsigned int)(time(NULL)));
    for(int i = 0; i < repeats; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            double val = (static_cast<double>(rand()) / RAND_MAX);
            ms->set(j, val);
        }
    }
}
void t6Helper2(double* d, int size, int repeats)
{
    srand((unsigned int)(time(NULL)));
    for(int i = 0; i < repeats; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            double val = (static_cast<double>(rand()) / RAND_MAX);
            d[j] = val;
        }
    }
}
// benchmark writing
void test6()
{
    const int size = 1000000;
    const int repetitions = 1000;
    TwoSegmentArray<false>* f = new TwoSegmentArray<false>(size);
    TwoSegmentArray<true>* t = new TwoSegmentArray<true>(size);
    double* d = new double[size];

    auto d_start = std::chrono::_V2::high_resolution_clock::now();
    t6Helper2(d, size, repetitions);
    auto d_end = std::chrono::_V2::high_resolution_clock::now();

    auto f_start = std::chrono::_V2::high_resolution_clock::now();
    t6Helper(f, size, repetitions);
    auto f_end = std::chrono::_V2::high_resolution_clock::now();

    auto f_startv2 = std::chrono::_V2::high_resolution_clock::now();
    t6Helper_v2(f, size, repetitions);
    auto f_endv2 = std::chrono::_V2::high_resolution_clock::now();

    auto p_start = std::chrono::_V2::high_resolution_clock::now();
    t6Helper(t, size, repetitions);
    auto p_end = std::chrono::_V2::high_resolution_clock::now();

    auto p_startv2 = std::chrono::_V2::high_resolution_clock::now();
    t6Helper_v2(t, size, repetitions);
    auto p_endv2 = std::chrono::_V2::high_resolution_clock::now();



    auto d_time = std::chrono::duration_cast<std::chrono::nanoseconds>(d_end - d_start).count()*1e-9;
    auto f_time = std::chrono::duration_cast<std::chrono::nanoseconds>(f_end - f_start).count()*1e-9;
    auto p_time = std::chrono::duration_cast<std::chrono::nanoseconds>(p_end - p_start).count()*1e-9;
    auto f_timev2 = std::chrono::duration_cast<std::chrono::nanoseconds>(f_endv2 - f_startv2).count()*1e-9;
    auto p_timev2 = std::chrono::duration_cast<std::chrono::nanoseconds>(p_endv2 - p_startv2).count()*1e-9;

    std::cout << "std double: " << d_time << "s\nheads only: " << f_time << "s\npairs:      " << p_time
                << "\nsheads func:      " << f_timev2 << "s\npairs func:    " << p_timev2 << "s" << std::endl;

    delete[] d, f, t;
}
// benchmark reading
template<class TwoSegmentArray>
double t7Helper(TwoSegmentArray* ms, int size, int repetitions)
{
    double t = 0.0;
    for(int i = 0; i < repetitions; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            t += (*ms)[j];
        }
    }
    return t;
}
template<class TwoSegmentArray>
double t7Helper_v2(TwoSegmentArray* ms, int size, int repetitions)
{
    double t = 0.0;
    for(int i = 0; i < repetitions; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            t += ms->read(j);
        }
    }
    return t;
}
double t7Helper2(double* d, int size, int repetitions)
{
    double t = 0.0;
    for(int i = 0; i < repetitions; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            t += d[j];
        }
    }
    return t;
}
void test7()
{
    const int size = 1000000;
    const int repetitions = 1000;
    TwoSegmentArray<false>* f = new TwoSegmentArray<false>(size);
    TwoSegmentArray<true>* t = new TwoSegmentArray<true>(size);
    double* d = new double[size];


    t6Helper2(d, size, 1);
    t6Helper(f, size, 1);
    t6Helper(t, size, 1);

    auto f_start = std::chrono::_V2::high_resolution_clock::now();
    double ftotal = t7Helper(f, size, repetitions);
    // double ftotal = t7Helper_v2(f, size, repetitions);
    auto f_end = std::chrono::_V2::high_resolution_clock::now();

    auto p_start = std::chrono::_V2::high_resolution_clock::now();
    double ptotal = t7Helper(t, size, repetitions);
    // double ptotal = t7Helper_v2(t, size, repetitions);
    auto p_end = std::chrono::_V2::high_resolution_clock::now();

    auto d_start = std::chrono::_V2::high_resolution_clock::now();
    double dtotal = t7Helper2(d, size, repetitions);
    auto d_end = std::chrono::_V2::high_resolution_clock::now();




    auto d_time = std::chrono::duration_cast<std::chrono::nanoseconds>(d_end - d_start).count()*1e-9;
    auto f_time = std::chrono::duration_cast<std::chrono::nanoseconds>(f_end - f_start).count()*1e-9;
    auto p_time = std::chrono::duration_cast<std::chrono::nanoseconds>(p_end - p_start).count()*1e-9;
    // auto f_timev2 = std::chrono::duration_cast<std::chrono::nanoseconds>(f_endv2 - f_startv2).count()*1e-9;
    // auto p_timev2 = std::chrono::duration_cast<std::chrono::nanoseconds>(p_endv2 - p_startv2).count()*1e-9;


    std::cout << "std double: " << d_time << "s\nheads only: " << f_time << "s\npairs:       " << p_time << "s\ntotals" << std::endl;
    std::cout << "std double: " << dtotal << "\nheads only: " << ftotal << "\npairs:       " << ptotal << std::endl;
    // std::cout << "\nheadsv2 only: " << f_timev2 << "s\npairs:       " << p_timev2 << "s" << std::endl;

    delete[] d, f, t;
}

// todo: replace with memcpy and whatever
// see if that provides some form of performance improvement

// spoiler: it does not.

int main(int argc, char** argv)
{
    if(argc < 2)
    {
        std::cerr << "need 1 argument" << std::endl;
        std::cerr << "usage: ./tests <test_case>" << std::endl;
        exit(1);
    }

    int t = 0;
    try
    {
        t = std::atoi(argv[1]);
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << "\nselecting default case...\n";
        t = 0;
    }

    switch(t)
    {
    case 1:
        std::cout << "testing manseg sum to standard double sum" << std::endl;
        test1(); // comparing manseg sum to std double sum & manseg conversion (static assigned)
        break;
    case 2:
        std::cout << "checking array type syntax and assignments" << std::endl;
        test2(); // checking array to array assignments
        break;
    case 3:
        std::cout << "checking initialisation and deconstruction" << std::endl;
        test3(); // initialisation
        break;
    case 4:
        std::cout << "nano benchmark that emulates pagerank, sort of, kind of (heads only)" << std::endl;
        test4(false); // checking to other assignment works correctly (i.e. x[i] += 2 and such)
        break;
    case 5:
        std::cout << "nano benchmark that emulates pagerank, sort of, kind of (using pairs)" << std::endl;
        test4(true);
        break;
    case 6:
        std::cout << "pico benchmark for sequential writes only" << std::endl;
        test6(); // write tests
        break;
    case 7:
        std::cout << "pico benchmark for sequential reads only" << std::endl;
        test7(); // read tests
        break;
    default:
        test1();
        std::cout << "\n";
        test2();
        std::cout << "\n";
        test3();
        std::cout << "\n";
        test4(false); // use heads only
        std::cout << "\n";
        test4(true); // use pairs
        std::cout << "\n";
        test6();
        std::cout << "\n";
        test7();
        std::cout << "\n";
    };

    return 0;
}

// TODO:
// Can probably rename TwoSegmentArray to TwoSegmentArray since there is no encapsuling thing anymore
// unless you want the wrapper thing you'd talked about which holds refs to low and high precision or whatever
// but that might not be so good if multiple precisions in future or whatever

// I think this still applies, but with pointer version of TwoSegmentArray
// dynamic version of TwoSegmentArray still doesn't like [] operator (need to wrap in *() to use)