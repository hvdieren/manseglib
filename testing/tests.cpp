#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>

#include <chrono>
#include <math.h>

// #include "../mantissaSegmentation.hpp"
// #include "../mantissaSegmentation_f.hpp"
#include "../mantissaSegmentation_dev.hpp"

using namespace ManSeg;

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

    HeadsArray arr(size);
    HeadsArray arr2(size);

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
    PairsArray tarr = arr.createFullPrecision();
    PairsArray tarr2 = arr2.createFullPrecision();

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

    PairsArray arr(size);
    PairsArray arr2(size);
    PairsArray arr3(size);
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

    HeadsArray x(4);
    PairsArray y = x.createFullPrecision();

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
void fill(double* d, double* dy, HeadsArray* x, HeadsArray* xy, 
    PairsArray* y, PairsArray* yy, int* v, const int& size)
{
    HeadsArray _x = *x, _xy = *xy;
    PairsArray _y = *y, _yy = *yy;
    for(int i = 0; i < size; ++i)
    {
        v[i] = (i+1);

        _x[i] = (double)(i+1)/0.1;
        _xy[i] = 0.0;

        _y[i] = (double)(i+1)/0.1;
        _yy[i] = 0.0;
        
        d[i] = (double)(i+1)/0.1;
        dy[i] = 0.0;        
    }
}

constexpr int size = 64000;
constexpr double val = 0.85;
// tried to replicate what pagerank was doing
// in a small tester function for micro benchmarking
// not exactly the same, but..
void test4_double_ver(double* d, double* dy, int* v, int* loc)
{
    auto dst = std::chrono::_V2::high_resolution_clock::now();

    for(int i = 0; i < size; ++i)
    {
        for(int j = 0; j < size; ++j) dy[i] += val*(d[loc[j]]/v[j]);
    }

    auto ded = std::chrono::_V2::high_resolution_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::nanoseconds>(ded - dst).count()*1e-9;

    std::cout << "std double = " << dt << "s\n";
}
void test4_heads_ver(HeadsArray* _x, HeadsArray* _xy, int* v, int* loc)
{
    HeadsArray x = *_x, xy = *_xy;
    auto mxst = std::chrono::_V2::high_resolution_clock::now();

    for(int i = 0; i < size; ++i)
    {
        for(int j = 0; j < size; ++j) xy[i] += val*(x[loc[j]]/v[j]);
    }

    auto mxed = std::chrono::_V2::high_resolution_clock::now();
    auto mt = std::chrono::duration_cast<std::chrono::nanoseconds>(mxed - mxst).count()*1e-9;

    std::cout  << "heads = " << mt << "s\n";
}
void test4_pairs_ver(PairsArray* _y, PairsArray* _yy, int* v, int* loc)
{
    PairsArray y = *_y, yy = *_yy;
    auto mxst = std::chrono::_V2::high_resolution_clock::now();

    for(int i = 0; i < size; ++i)
    {
        for(int j = 0; j < size; ++j) yy[i] += val*(y[loc[j]]/v[j]);
    }

    auto mxed = std::chrono::_V2::high_resolution_clock::now();
    auto mt = std::chrono::duration_cast<std::chrono::nanoseconds>(mxed - mxst).count()*1e-9;

    std::cout  << "pairs = " << mt << "s\n";
}
void test4_head_acc_ver(HeadsArray* _x, HeadsArray* _xy, int* v, int* loc)
{
    HeadsArray x = *_x, xy = *_xy;
    double* acc = new double[size];

    auto mxst = std::chrono::_V2::high_resolution_clock::now();
    for(int i = 0; i < size; ++i)
    {
        for(int j = 0; j < size; ++j) acc[i] += val*(x[loc[j]]/v[j]);

        xy[i] = acc[i];
    }
    auto mxed = std::chrono::_V2::high_resolution_clock::now();
    auto mt = std::chrono::duration_cast<std::chrono::nanoseconds>(mxed - mxst).count()*1e-9;

    delete[] acc;
    std::cout  << "heads w/acc = " << mt << "s\n";
}
void test4_pairs_acc_ver(PairsArray* _y, PairsArray* _yy, int* v, int* loc)
{
    PairsArray y = *_y, yy = *_yy;
    double* acc = new double[size];

    auto mxst = std::chrono::_V2::high_resolution_clock::now();
    for(int i = 0; i < size; ++i)
    {
        for(int j = 0; j < size; ++j) acc[i] += val*(y[loc[j]]/v[j]);

        yy[i] = acc[i];
    }
    auto mxed = std::chrono::_V2::high_resolution_clock::now();
    auto mt = std::chrono::duration_cast<std::chrono::nanoseconds>(mxed - mxst).count()*1e-9;

    delete[] acc;
    std::cout  << "pairs w/acc = " << mt << "s\n";
}
void test4()
{
    HeadsArray* x = new HeadsArray(size);
    HeadsArray* xy = new HeadsArray(size);       // x, xy = heads
    PairsArray* y = new PairsArray(size);
    PairsArray* yy = new PairsArray(size);       // y, yy = pairs

    int* v = new int[size];
    double* d = new double[size];
    double* dy = new double[size];              // d, dy = doubles

    srand(1992);
    // randomise data accesses
    int* loc = new int[size];
    for(int i = 0; i < size; ++i)
        loc[i] = (rand() % size);

    // fill arrays initially
    fill(d, dy, x, xy, y, yy, v, size);

    std::cout << "round 1 - read and write direct \n";
    // doubles
    test4_double_ver(d, dy, v, loc);
    // heads
    test4_heads_ver(x, xy, v, loc);
    // pairs
    test4_pairs_ver(y, yy, v, loc);

    std::cout << "round 2 - read from manseg, write to double accumulator \n";
    // reset
    fill(d, dy, x, xy, y, yy, v, size);

    // doubles
    test4_double_ver(d, dy, v, loc);
    // heads w acc
    test4_head_acc_ver(x, xy, v, loc);
    // pairs w acc
    test4_pairs_acc_ver(y, yy, v, loc);

    delete[] d, dy, v, loc;
    x->del(); xy->del();
    y->del(); yy->del();
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

// helper for test5 & test6
template<class TwoSegArray>
void t5Helper(TwoSegArray* ms, int size, int repeats, int* loc)
{
    for(int i = 0; i < repeats; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            double val = (static_cast<double>(rand()) / RAND_MAX);
            (*ms)[loc[j]] = val;
        }
    }
}
template<class TwoSegArray>
void t5Helper_v2(TwoSegArray* ms, int size, int repeats, int* loc)
{
    for(int i = 0; i < repeats; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            double val = (static_cast<double>(rand()) / RAND_MAX);
            ms->set(loc[j], val);
        }
    }
}
void t5Helper2(double* d, int size, int repeats, int* loc)
{
    for(int i = 0; i < repeats; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            double val = (static_cast<double>(rand()) / RAND_MAX);
            d[loc[j]] = val;
        }
    }
}
// benchmark writing
void test5()
{
    const int size = 1000000;
    const int repetitions = 1000;
    HeadsArray* f = new HeadsArray(size);
    PairsArray* t = new PairsArray(size);
    double* d = new double[size];
    
    srand(1992);
    int* loc = new int[size];
    for(int i = 0; i < size; ++i) loc[i] = rand() % size;

    auto d_start = std::chrono::_V2::high_resolution_clock::now();
    t5Helper2(d, size, repetitions, loc);
    auto d_end = std::chrono::_V2::high_resolution_clock::now();

    auto f_start = std::chrono::_V2::high_resolution_clock::now();
    t5Helper(f, size, repetitions, loc);
    auto f_end = std::chrono::_V2::high_resolution_clock::now();

    auto f_startv2 = std::chrono::_V2::high_resolution_clock::now();
    t5Helper_v2(f, size, repetitions, loc);
    auto f_endv2 = std::chrono::_V2::high_resolution_clock::now();

    auto p_start = std::chrono::_V2::high_resolution_clock::now();
    t5Helper(t, size, repetitions, loc);
    auto p_end = std::chrono::_V2::high_resolution_clock::now();

    auto p_startv2 = std::chrono::_V2::high_resolution_clock::now();
    t5Helper_v2(t, size, repetitions, loc);
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
double t6Helper(TwoSegArray* ms, int size, int repetitions, int* loc)
{
    double t = 0.0;
    for(int i = 0; i < repetitions; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            t += (*ms)[loc[j]];
        }
    }
    return t;
}
template<class TwoSegArray>
double t6Helper_v2(TwoSegArray* ms, int size, int repetitions, int* loc)
{
    double t = 0.0;
    for(int i = 0; i < repetitions; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            t += ms->read(loc[j]);
        }
    }
    return t;
}
double t6Helper2(double* d, int size, int repetitions, int* loc)
{
    srand((unsigned int)(time(NULL)));
    double t = 0.0;
    for(int i = 0; i < repetitions; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            t += d[loc[j]];
        }
    }
    return t;
}
void test6()
{
    const int size = 1000000;
    const int repetitions = 2000;
    HeadsArray* f = new HeadsArray(size);
    PairsArray* t = new PairsArray(size);
    double* d = new double[size];

    srand(1992);
    int* loc = new int[size];
    for(int i = 0; i < size; ++i) loc[i] = rand() % size;

    t5Helper2(d, size, 1, loc);
    t5Helper(f, size, 1, loc);
    t5Helper(t, size, 1, loc);

    auto f_start = std::chrono::_V2::high_resolution_clock::now();
    double ftotal = t6Helper(f, size, repetitions, loc);
    // double ftotal = t7Helper_v2(f, size, repetitions);
    auto f_end = std::chrono::_V2::high_resolution_clock::now();

    auto p_start = std::chrono::_V2::high_resolution_clock::now();
    double ptotal = t6Helper(t, size, repetitions, loc);
    // double ptotal = t7Helper_v2(t, size, repetitions);
    auto p_end = std::chrono::_V2::high_resolution_clock::now();

    auto d_start = std::chrono::_V2::high_resolution_clock::now();
    double dtotal = t6Helper2(d, size, repetitions, loc);
    auto d_end = std::chrono::_V2::high_resolution_clock::now();


    auto d_time = std::chrono::duration_cast<std::chrono::nanoseconds>(d_end - d_start).count()*1e-9;
    auto f_time = std::chrono::duration_cast<std::chrono::nanoseconds>(f_end - f_start).count()*1e-9;
    auto p_time = std::chrono::duration_cast<std::chrono::nanoseconds>(p_end - p_start).count()*1e-9;

    std::cout << "std double: " << d_time << "s\nheads only: " << f_time << "s\npairs:       " << p_time << "s\ntotals" << std::endl;
    std::cout << "std double: " << dtotal << "\nheads only: " << ftotal << "\npairs:       " << ptotal << std::endl;

    f->del(); t->del();
    delete[] d, f, t, loc;
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
        test4();
        std::cout << "\n";
        test5();
        std::cout << "\n";
        test6();
        std::cout << "\n";
    };

    return 0;
}
