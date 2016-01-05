#include <iostream>
#include <fprotos.hxx>

int main (int argc, char * argv[])
{

    std::cout << "Area of circle with unit radius is:" << circleArea(1.0) << std::endl;
    std::cout << "Area of square with side length " << 2 << " is: " << squareArea(2.0) << std::endl;

    return 1.0;
}

