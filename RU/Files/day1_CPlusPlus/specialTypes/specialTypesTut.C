#include <scalar.H>
#include <label.H>
#include <Switch.H>
#include <word.H>
#include <autoPtr.H>
#include <tmp.H>
#include <Xfer.H>
#include <cmath>
#include <iostream>

using namespace Foam;

scalar mult (scalar a, scalar b)
{
    return a*b;
}

int main (int argc, char * argv[])
{
    //pointers
    label i = 1;
    label& li = i;
    label *pi = &i;
    *pi++;
    autoPtr<label> ap;
    ap.reset(pi);
    label *pi2 = ap.ptr();
    
    //arrays
    scalar da[20];
    scalar * db;
    db = new scalar [20];
    for (label q=0; q<20; q++)
    {
	da[q] = q;
	*(db+q) = mult(da[q], q);
    }
    delete db;

    return 0;
}

