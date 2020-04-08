#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>

#include <chrono>
#include <math.h>

#include "util.h"
#include "../mantissaSegmentation_dev.hpp"

using namespace ManSeg;
using namespace std;

int main()
{
	cout << fixed << setprecision(16);
	const size_t size = 100UL;

	PairsArray arr(size);
    PairsArray arr2(size);
    PairsArray arr3(size);
    double d[size];
    double d2[size];
    double d3[size];

    srand(2517);
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

	int return_code = 0;
    for(size_t i = 0; i < size; ++i)
    {
        if((double)arr[i] != (double)arr3[i])
		{
            std::cout << "arr[" << i << "]=" << arr[i] << " != arr3[" << i << "]=" << arr3[i] << std::endl;
			return_code = 1;
		}
    }
    for(size_t i = 0; i < size; ++i)
    {
        if(d2[i] != (double)arr2[i])
		{
            std::cout << "d2[" << i << "]=" << d2[i] << " != arr2[" << i << "]=" << arr2[i] << std::endl;
			return_code = 1;
		}
    }

    arr.del();
    arr2.del();
    arr3.del();

	if(return_code == 0)
		cout << "test passed!" << endl;

	return return_code;
}
