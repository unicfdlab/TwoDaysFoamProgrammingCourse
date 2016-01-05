/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    simplexTest

Description
    Tests the simplex algorithm using a stand-alone function paraboloidSin

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "SimplexNelderMead.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


class paraboloidSin
{
    scalar p0_, p1_, p2_, p3_, p4_;

public:

    paraboloidSin(scalar p0, scalar p1, scalar p2, scalar p3, scalar p4)
    :
        p0_(p0),
        p1_(p1),
        p2_(p2),
        p3_(p3),
        p4_(p4)
    {}

    label nArgs() const
    {
        return 3;
    }

    Tuple2<scalar, bool> operator()(const scalarField& xv) const
    {
        scalar x = xv[0];
        scalar y = xv[1];
        scalar z = xv[2];

        Info << "Calling paraboloidSin x = "
            << x << " y = " << y << " z = " << z << endl;

        //if ( z < 0.0 )
        //{
        //    return Tuple2<scalar,bool>(0.0, false);
        //}

        return Tuple2<scalar, bool>
        (
            p2_*sqr(x - p0_)
          + p3_*sqr(y - p1_)
          + 10*sqr(z)*(1 + Foam::sin(x))
          + p4_,
            true
        );
    }
};


int main(int argc, char *argv[])
{
    paraboloidSin f(1.0, 2.0, 10.0, 20.0, 30.0);

    // Starting point
    scalarField startPoint(f.nArgs());
    startPoint[0] = 5.0;
    startPoint[1] = 7.0;
    startPoint[2] = 0.0;

    scalarField lambda(f.nArgs());
    lambda[0] = 1.0;
    lambda[1] = 1.0;
    lambda[2] = 1.0;

    SimplexNelderMead<paraboloidSin> simplex(f, startPoint, lambda);

    Info<< "before iterate" << endl;

    label iter = 0;
    while ( simplex.size() > 0.001 )
    {
        simplex.iterate();

        Info << "iter = " << iter
            << " minPos = " << simplex.minCoord()
            << " v = " << simplex.min()
            << " size = " << simplex.size() << endl;

        iter++;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
