#include <scalar.H>
#include <label.H>
#include <Switch.H>
#include <word.H>
#include <cmath>
#include <iostream>

using namespace Foam;

int main (int argc, char * argv[])
{

    label i = 1;
    int   j = 2;
    label k = i;

    scalar pi1 = 3.14;
    double d   = M_PI - pi1;
    
    Switch eq = (d <= 0.0);
    
    word goodStr = "abc";
    word badStr  = "abc abc";
    
    return 0;
}

