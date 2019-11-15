#include <iostream>
#include "mantissaSegmentation.hpp"

int main()
{
    ManSegArray* arr = new ManSegArray(100);

    double d = 1.25;
    arr->setPair(0, d);

    arr->setPair(50, d);

    arr->print();

    return 0;
}