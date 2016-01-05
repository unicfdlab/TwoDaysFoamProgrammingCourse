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
    optimiserFoam

Description
    Generic optimisation solver using the simplex algorithm

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "SimplexNelderMead.H"
#include "objectiveFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createSimplex.H"

    label simplexIter = 0;

    while (simplex.size() > tolerance)
    {
        simplexIter++;
        simplex.iterate();

        Info << "simplex optimisation iteration = " << simplexIter
            << " minPos = " << simplex.minCoord()
            << " v = " << simplex.min()
            << " size = " << simplex.size() << endl;

        if (simplexIter > maxIter)
        {
            Info<< "Convergence not achieved" << endl;
            break;
        }
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
