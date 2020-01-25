#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>

#include <chrono>
#include <math.h>

#include "mantissaSegmentation.hpp"
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


// todo: replace with memcpy and whatever
// see if that provides some form of performance improvement

// spoiler: it does not.

// spoiler: changing other things to assignment/creation things also does not make an improvement - it is worse.
// not entirely sure if *reinterpret_cast<double>(&l) is the best way to do things..
// but also other things don't seem to work /shrug

using namespace std;

// continuous memory implementation of unsigned char* slow as balls bc memory locality sucks
// try splitting into heads and tails
// for increased head locality.
// pairs might suffer dramatically though.

// splitting things into heads and tails for locality made things even worse

int main()
{
    const int size = 1e6;
    int* a = new int[size];
    int* b = new int[size];

    srand(time(NULL));
    for(int i = 0; i < size; ++i)
    {
        a[i] = rand() % 50 + 10;
        b[i] = rand() % 50 + 10;
    }

    // std::cout << "a = ";
    // for(int i = 0; i < size; ++i)
    //     std::cout << a[i] << " ";
    // std::cout << "\n";

    // std::cout << "b = ";
    // for(int i = 0; i < size; ++i)
    //     std::cout << b[i] << " ";
    // std::cout << "\n\n do sort \n\n";

    // SORT
    quicksort_pair(a, b, 0, size-1);

    // std::cout << "a = ";
    // for(int i = 0; i < size; ++i)
    //     std::cout << a[i] << " ";
    // std::cout << "\n";

    // std::cout << "b = ";
    // for(int i = 0; i < size; ++i)
    //     std::cout << b[i] << " ";
    // std::cout << "\n";

    return 0;
}

int main1(int argc, char** argv)
{
    // unsigned char* c = new unsigned char[8];
    // cout << setprecision(16);
    // double d = 0.123456789987654321;
    // // unsigned char* c = reinterpret_cast<unsigned char*>(&d);

    // // for(int i = 0; i < 4; ++i) c[i] = 0;

    // unsigned char* cpy = reinterpret_cast<unsigned char*>(&d);

    // // memcpy(c+4, cpy+4, 4);
    // for(int i = 7; i >=4; --i)
    //     c[i] = cpy[i];

    // double e = *reinterpret_cast<double*>(c);

    // cout << "d = " << d << endl;
    // printBinary(d);
    // cout << endl;
    // cout << "e = " << e << endl;
    // printBinary(e);
    // cout << endl;

    // return 0;

    int size = 5;
    TwoSegmentArray<false> f(size);
    TwoSegmentArray<false> f1(size);
    f.set(0, 0.1);
    f.set(1, 1.4f);
    f.set(2, 22);
    f.set(3, 45L);
    f.set(4, true);

    // f.setPair(0, 0.1);
    // f.setPair(1, 1.4f);
    // f.setPair(2, 22);
    // f.setPair(3, 45L);
    // f.setPair(4, true);

    // f[0] = 0.1;
    // f[1] = 1.4f;
    // f[2] = 22;
    // f[3] = 45L;
    // f[4] = true;

    f1[0] = f[0]; // for some reason, now, instead of using an operator =
    f1[1] = f[1]; // this is using a copy constructor
    f1[2] = f[2]; // something which has not existed until this point - side effect of using unsigned char* it seems
    f1[3] = f[3];
    f1[4] = f[4];

    for(int i = 0; i < size; ++i)
        std::cout << i << ":f -> " << f[i] << "\t\tf1 -> " << f1[i] << "\n";
    std::cout << "\n";

    return 0;
}

// performance is legitimately awful, and i don't really know what it is that's wrong
// reads and writes are slow but how else do you combine unsigned ints to a double?
// unsigned int -> unsigned long | unsigned int -> double


// TODO:
// Can probably rename TwoSegmentArray to TwoSegmentArray since there is no encapsuling thing anymore
// unless you want the wrapper thing you'd talked about which holds refs to low and high precision or whatever
// but that might not be so good if multiple precisions in future or whatever

// dynamically allocated version of TwoSegmentArray requires the use of (*a)[] to get access to operator[] without
// using a->operator[]()
// ..not sure if there is a way around this.
