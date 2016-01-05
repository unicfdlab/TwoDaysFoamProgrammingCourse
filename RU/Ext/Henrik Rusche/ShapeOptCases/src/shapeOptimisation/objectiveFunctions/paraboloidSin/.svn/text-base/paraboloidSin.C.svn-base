/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "paraboloidSin.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(paraboloidSin, 0);
    addToRunTimeSelectionTable(objectiveFunction, paraboloidSin, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::paraboloidSin::paraboloidSin
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    objectiveFunction(mesh, dict),
    // Note: read from private dictionary.  HJ, 4/Apr/2010
    p0_(readScalar(functionProperties().lookup("p0"))),
    p1_(readScalar(functionProperties().lookup("p1"))),
    p2_(readScalar(functionProperties().lookup("p2"))),
    p3_(readScalar(functionProperties().lookup("p3"))),
    p4_(readScalar(functionProperties().lookup("p4")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Tuple2<Foam::scalar, bool>
Foam::paraboloidSin::operator()
(
    const scalarField& xv
)
{
    scalar x = xv[0];
    scalar y = xv[1];
    scalar z = xv[2];

    Info<< "Calling paraboloidSin x = " << x
        << " y = " << y << " z = " << z << endl;

//     if ( z < 0.0 )
//     {
//         return Tuple2<scalar,bool>(0.0, false);
//     }

    return Tuple2<scalar,bool>
    (
        p2_*sqr(x - p0_)
      + p3_*sqr(y - p1_)
      + 10*sqr(z)*(1 + Foam::sin(x))
      + p4_,
        true
    );
}


// ************************************************************************* //
