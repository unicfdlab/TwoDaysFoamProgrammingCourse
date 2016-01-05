#include <scalar.H>
#include <label.H>
#include <Switch.H>
#include <word.H>
#include <List.H>
#include <cmath>
#include <iostream>

using namespace Foam;

int main (int argc, char * argv[])
{
    //for example
    scalar sum = 0.0;
    label N = 10;
    for (label j=0; j<N; j++)
    {
	sum += j;
    }
    
    //do-while example
    label k = N;
    do
    {
	k--;
	sum -= k;
    }
    while ( k > 0);
    
    //if example
    if (sum == 0)
    {
	k = 0;
    }
    else if (sum > 0)
    {
	k = 1;
    }
    else
    {
	k = -1;
    }
    
    //switch example
    switch (k)
    {
	case 0:
	{
	    sum = 0.0;
	    break;
	}
	case 1:
	{
	    sum = 1.0;
	    break;
	}
	case -1:
	{
	    sum = -1.0;
	    break;
	}
	default:
	{
	    sum = VSMALL;
	}
    }
    
    //forAll example
    List<label> l(N);
    forAll (l, i)
    {
	l[i] = 2*i;
    }
    return 0;
}

