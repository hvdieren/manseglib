#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>

#include <chrono>
#include <math.h>

// #include "../mantissaSegmentation.hpp"
#include "../mantissaSegmentation_dev.hpp"

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
template<class TwoSegArray>
void f(TwoSegArray arr, TwoSegArray arr2, double* darr, double* darr2, size_t size)
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
    // make calculations into templated function taking TwoSegArray as type parameter
    std::cout << std::fixed << std::setprecision(16);
    const size_t size = 10000000UL;

    TwoSegArray<false> arr(size);
    TwoSegArray<false> arr2(size);

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
    TwoSegArray<true> tarr = arr.createFullPrecision();
    TwoSegArray<true> tarr2 = arr2.createFullPrecision();

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

    TwoSegArray<true> arr(size);
    TwoSegArray<true> arr2(size);
    TwoSegArray<true> arr3(size);
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

    TwoSegArray<false> x(4);
    TwoSegArray<true>y = x.createFullPrecision();

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
void fill(double* d, double* dy, TwoSegArray<false>& x, TwoSegArray<false>& xy, 
    TwoSegArray<true>& y, TwoSegArray<true>& yy, int* v, const int& size)
{
    for(int i = 0; i < size; ++i)
    {
        v[i] = (i+1);

        x[i] = (double)(i+1)/0.1;
        xy[i] = 0.0;

        y[i] = (double)(i+1)/0.1;
        yy[i] = 0.0;
        
        d[i] = (double)(i+1)/0.1;
        dy[i] = 0.0;        
    }
}

constexpr int size = 64000;
constexpr double val = 0.85;
// tried to replicate what pagerank was doing
// in a small tester function for micro benchmarking
// not exactly the same, but somewhat close.
void test4_double_ver(double* d, double* dy, int* v)
{
    auto dst = std::chrono::_V2::high_resolution_clock::now();

    for(int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j)
            dy[i] += val*(d[size-1-j]/v[j]);

    auto ded = std::chrono::_V2::high_resolution_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::nanoseconds>(ded - dst).count()*1e-9;

    std::cout << "std double = " << dt << "s\n";
}
void test4_heads_ver(TwoSegArray<false>& x, TwoSegArray<false>& xy, int* v)
{
    auto mxst = std::chrono::_V2::high_resolution_clock::now();

    for(int i = 0; i < size; ++i)
    {
        for(int j = 0; j < size; ++j) xy[i] += val*(x[size-1-j]/v[j]);
    }

    auto mxed = std::chrono::_V2::high_resolution_clock::now();
    auto mt = std::chrono::duration_cast<std::chrono::nanoseconds>(mxed - mxst).count()*1e-9;

    std::cout  << "heads = " << mt << "s\n";
}
void test4_pairs_ver(TwoSegArray<true>& y, TwoSegArray<true>& yy, int* v)
{
    auto mxst = std::chrono::_V2::high_resolution_clock::now();

    for(int i = 0; i < size; ++i)
    {
        for(int j = 0; j < size; ++j) yy[i] += val*(y[size-1-j]/v[j]);
    }

    auto mxed = std::chrono::_V2::high_resolution_clock::now();
    auto mt = std::chrono::duration_cast<std::chrono::nanoseconds>(mxed - mxst).count()*1e-9;

    std::cout  << "pairs = " << mt << "s\n";
}
void test4_head_acc_ver(TwoSegArray<false>& x, TwoSegArray<false>& xy, int* v)
{
    auto mxst = std::chrono::_V2::high_resolution_clock::now();

    double* acc = new double[size];
    for(int i = 0; i < size; ++i)
    {
        for(int j = 0; j < size; ++j) acc[i] += val*(x[size-1-j]/v[j]);

        xy[i] = acc[i];
    }
    delete[] acc;

    auto mxed = std::chrono::_V2::high_resolution_clock::now();
    auto mt = std::chrono::duration_cast<std::chrono::nanoseconds>(mxed - mxst).count()*1e-9;

    std::cout  << "heads w/acc = " << mt << "s\n";
}
void test4_pairs_acc_ver(TwoSegArray<true>& y, TwoSegArray<true>& yy, int* v)
{
    auto mxst = std::chrono::_V2::high_resolution_clock::now();

    double* acc = new double[size];
    for(int i = 0; i < size; ++i)
    {
        for(int j = 0; j < size; ++j) acc[i] += val*(y[size-1-j]/v[j]);

        yy[i] = acc[i];
    }
    delete[] acc;

    auto mxed = std::chrono::_V2::high_resolution_clock::now();
    auto mt = std::chrono::duration_cast<std::chrono::nanoseconds>(mxed - mxst).count()*1e-9;

    std::cout  << "pairs w/acc = " << mt << "s\n";
}
void test4()
{
    TwoSegArray<false> x(size), xy(size);       // x = heads
    TwoSegArray<true> y(size), yy(size);        // y = pairs

    int* v = new int[size];
    double* d = new double[size];
    double* dy = new double[size];

    // fill arrays initially
    std::cout << "round 1 - read and write to manseg \n";
    fill(d, dy, x, xy, y, yy, v, size);

    // doubles
    test4_double_ver(d, dy, v);
    // heads
    test4_heads_ver(x, xy, v);
    // pairs
    test4_pairs_ver(y, yy, v);

    std::cout << "round 2 - read from manseg, write to double acc \n";
    // reset
    fill(d, dy, x, xy, y, yy, v, size);

    // doubles
    test4_double_ver(d, dy, v);
    // heads w acc
    test4_head_acc_ver(x, xy, v);
    // pairs w acc
    test4_pairs_acc_ver(y, yy, v);

    delete[] d, dy, v;
    x.del(); xy.del();
    y.del(); yy.del();
}

// reading SNAP file format test
void testSNAP()
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
template<class TwoSegArray>
void t5Helper(TwoSegArray* ms, int size, int repeats)
{
    srand((unsigned int)(time(NULL)));
    for(int i = 0; i < repeats; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            double val = (static_cast<double>(rand()) / RAND_MAX);
            (*ms)[(int)(rand()/size)] = val;
        }
    }
}
template<class TwoSegArray>
void t5Helper_v2(TwoSegArray* ms, int size, int repeats)
{
    srand((unsigned int)(time(NULL)));
    for(int i = 0; i < repeats; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            double val = (static_cast<double>(rand()) / RAND_MAX);
            ms->set((int)(rand()/size), val);
        }
    }
}
void t5Helper2(double* d, int size, int repeats)
{
    srand((unsigned int)(time(NULL)));
    for(int i = 0; i < repeats; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            double val = (static_cast<double>(rand()) / RAND_MAX);
            d[(int)(rand()/size)] = val;
        }
    }
}
// benchmark writing
void test5()
{
    const int size = 1000000;
    const int repetitions = 1000;
    TwoSegArray<false>* f = new TwoSegArray<false>(size);
    TwoSegArray<true>* t = new TwoSegArray<true>(size);
    double* d = new double[size];

    auto d_start = std::chrono::_V2::high_resolution_clock::now();
    t5Helper2(d, size, repetitions);
    auto d_end = std::chrono::_V2::high_resolution_clock::now();

    auto f_start = std::chrono::_V2::high_resolution_clock::now();
    t5Helper(f, size, repetitions);
    auto f_end = std::chrono::_V2::high_resolution_clock::now();

    auto f_startv2 = std::chrono::_V2::high_resolution_clock::now();
    t5Helper_v2(f, size, repetitions);
    auto f_endv2 = std::chrono::_V2::high_resolution_clock::now();

    auto p_start = std::chrono::_V2::high_resolution_clock::now();
    t5Helper(t, size, repetitions);
    auto p_end = std::chrono::_V2::high_resolution_clock::now();

    auto p_startv2 = std::chrono::_V2::high_resolution_clock::now();
    t5Helper_v2(t, size, repetitions);
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
template<class TwoSegArray>
double t6Helper(TwoSegArray* ms, int size, int repetitions)
{
    srand((unsigned int)(time(NULL)));
    double t = 0.0;
    for(int i = 0; i < repetitions; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            t += (*ms)[(int)(rand()/size)];
        }
    }
    return t;
}
template<class TwoSegArray>
double t6Helper_v2(TwoSegArray* ms, int size, int repetitions)
{
    srand((unsigned int)(time(NULL)));
    double t = 0.0;
    for(int i = 0; i < repetitions; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            t += ms->read((int)(rand()/size));
        }
    }
    return t;
}
double t6Helper2(double* d, int size, int repetitions)
{
    srand((unsigned int)(time(NULL)));
    double t = 0.0;
    for(int i = 0; i < repetitions; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            t += d[(int)(rand()/size)];
        }
    }
    return t;
}
void test6()
{
    const int size = 1000000;
    const int repetitions = 2000;
    TwoSegArray<false>* f = new TwoSegArray<false>(size);
    TwoSegArray<true>* t = new TwoSegArray<true>(size);
    double* d = new double[size];


    t5Helper2(d, size, 1);
    t5Helper(f, size, 1);
    t5Helper(t, size, 1);

    auto f_start = std::chrono::_V2::high_resolution_clock::now();
    double ftotal = t6Helper(f, size, repetitions);
    // double ftotal = t7Helper_v2(f, size, repetitions);
    auto f_end = std::chrono::_V2::high_resolution_clock::now();

    auto p_start = std::chrono::_V2::high_resolution_clock::now();
    double ptotal = t6Helper(t, size, repetitions);
    // double ptotal = t7Helper_v2(t, size, repetitions);
    auto p_end = std::chrono::_V2::high_resolution_clock::now();

    auto d_start = std::chrono::_V2::high_resolution_clock::now();
    double dtotal = t6Helper2(d, size, repetitions);
    auto d_end = std::chrono::_V2::high_resolution_clock::now();


    auto d_time = std::chrono::duration_cast<std::chrono::nanoseconds>(d_end - d_start).count()*1e-9;
    auto f_time = std::chrono::duration_cast<std::chrono::nanoseconds>(f_end - f_start).count()*1e-9;
    auto p_time = std::chrono::duration_cast<std::chrono::nanoseconds>(p_end - p_start).count()*1e-9;

    std::cout << "std double: " << d_time << "s\nheads only: " << f_time << "s\npairs:       " << p_time << "s\ntotals" << std::endl;
    std::cout << "std double: " << dtotal << "\nheads only: " << ftotal << "\npairs:       " << ptotal << std::endl;

    delete[] d, f, t;
}

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
        std::cout << "checking array type overloads" << std::endl;
        test2(); // checking array to array assignments
        break;
    case 3:
        std::cout << "ensure .set/setPair functions the same as overloaded operators" << std::endl;
        test3(); // .set/setPair functions / construction/deconstruction
        break;
    case 4:
        std::cout << "nano benchmark that emulates pagerank, sort of, kind of" << std::endl;
        test4(); // test basic math potential in a few ways
        break;
    case 5:
        std::cout << "pico benchmark for pseudo-random writes" << std::endl;
        test5();
        break;
    case 6:
        std::cout << "pico benchmark for pseudo-random reads" << std::endl;
        test6(); // read tests
        break;
    default:
        test1();
        std::cout << "\n";
        test2();
        std::cout << "\n";
        test3();
        std::cout << "\n";
        test4(); // use heads only
        std::cout << "\n";
        test5();
        std::cout << "\n";
        test6();
        std::cout << "\n";
    };

    return 0;
}
