#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>

#include <chrono>
#include <math.h>

#include "util.h"
#include "../manseglib.hpp"

using namespace ManSeg;
using namespace std;

int main()
{
	cout << fixed << setprecision(16);

	int return_code = 0;

	HeadsArray x(4);
    PairsArray y = x.createFullPrecision();

    y[0] = 1.26;
    y.setPair(1, 1.26);

    y[2] = 1.22;
    y.set(3, 1.22);

    if(y[0] != y[1])
	{
        std::cout << "y[0] != y[1], values: y[0]=" << y[0] << " y[1]=" << y[1] << "\n";
		return_code = 1;
	}
    if(y[2] != y[3])
	{
        std::cout << "y[2] != y[3], values: y[2]=" << y[2] << " y[3]=" << y[3] << "\n";
		return_code = 1;
	}

    x[0] = 1.26;
    x.set(1, 1.26);

    x[2] = 1.22;
    x.set(3, 1.22);

    if(x[0] != x[1])
	{
        std::cout << "x[0] != x[1], values: x[0]=" << x[0] << " x[1]=" << x[1] << "\n";
		return_code = 1;
	}

    if(x[2] != x[3])
	{
        std::cout << "x[2] != x[3], values: x[2]=" << x[2] << " x[3]=" << x[3] << "\n";
		return_code = 1;
	}

    x.del();
    y.del();

	if(return_code == 0)
		cout << "test passed!" << endl;

	return return_code;
}