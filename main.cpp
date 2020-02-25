#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>

#include <chrono>
#include <math.h>

// #include "mantissaSegmentation.hpp"
#include "mantissaSegmentation_f.hpp"

using namespace ManSeg;

// print binary of double
// void printBinary(double d)
// {
//     unsigned long l = *reinterpret_cast<unsigned long*>(&d);
//     for(int i = 63; i >= 0; --i)
//     {
//         std::cout << ((l >> i) & 1);
//         if(i == 63 || i == 52)
//             std::cout << " ";
//         else if(i == 32)
//             std::cout << "|";
//     }
//     std::cout << std::endl;
// }

// // print binary of float
// void printBinary(float f)
// {
//     unsigned int l = *reinterpret_cast<unsigned int*>(&f);
//     for(int i = 31; i >= 0; --i)
//     {
//         std::cout << ((l >> i) & 1);
//         if(i == 31 || i == 23)
//             std::cout << " ";
//     }
//     std::cout << std::endl;
// }

// // print binary of int
// void printBinary(int l)
// {
//     for(int i = 31; i >= 0; --i)
//     {
//         std::cout << ((l >> i) & 1);
//         if(i % 8 == 0 && i > 0)
//             std::cout << " ";
//     }
// }

int main(int argc, char** argv)
{
    std::cout << std::setprecision(15);
    std::cout.setf(std::ios::fixed, std::ios::floatfield);

    using fArray = TwoSegArray<false>;
    using tArray = TwoSegArray<true>;

    const int size = 6;

    fArray f(size);
    tArray t(size);

    // f.setPair(0, 0.1);
    // f.setPair(1, 1.4f);
    // f.setPair(2, 22);
    // f.setPair(3, 45L);
    // f.setPair(4, true);

    f[0] = 0.385729789733887;
    f[1] = 1.4f;
    f[2] = 22;
    f[3] = 45L;
    f[4] = true;
    f[5] = 1.25e-308;
 
    t[0] = f[0]; 
    t[1] = f[1]; 
    t[2] = f[2];
    t[3] = f[3];
    t[4] = f[4];
    t[5] = f[5];

    for(int i = 0; i < size; ++i)
        std::cout << i << ":f -> " << f[i] << "\t\tt -> " << t[i] << "\n";
    std::cout << "\n";

    return 0;
}



// int main(int argc, char** argv)
// {
//     // what i want to do:
//     // take double, split into two parts (floats)
//     // take parts, reassemble

//     // todo: try converting mantissaSegmentation_f to use ints? again as underlying type but __m64 potentially (vector of ints)
//     // should then be able to convert directly from __m64 -> double, or from __m128d -> __m64 ->> at no cost <<-

//     double d = 1.2;
//     std::cout << "begin with = " << d << "\n";
    
//     __m128d init = {d, 0};
//     __m128 fl = _mm_castpd_ps(init);

//     for(int i = 0; i < 2; ++i)
//         std::cout << init[i] << "\n";
    
//     std::cout << "\n";

//     // fl[0] = tail, fl[1] = head
//     for(int i = 0; i < 4; ++i)
//         std::cout << fl[i] << "\n";

//     __m128d d_out = _mm_castps_pd(fl);
//     double res;
//     res = d_out[0];

//     std::cout << "res = " << res << "\n";    

//     return 0;
// }
