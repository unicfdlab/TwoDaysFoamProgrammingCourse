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

Class
    paraboloidSin

Description
    Simple objective function for testing purposes

Author
    Henrik Rusche and Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    paraboloidSin.C

\*---------------------------------------------------------------------------*/

#ifndef paraboloidSin_H
#define paraboloidSin_H

#include "objectiveFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class paraboloidSin Declaration
\*---------------------------------------------------------------------------*/

class paraboloidSin
:
    public objectiveFunction
{
    // Private data

        //- Parameters
        scalar p0_, p1_, p2_, p3_, p4_;

public:

    //- Runtime type information
    TypeName("paraboloidSin");


    // Constructors

        //- Construct from dictionary
        paraboloidSin
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~paraboloidSin()
        {}


    // Member Functions

        //- Return number of arguments
        label nArgs() const
        {
            return 3;
        }

        //- Evaluate and return objective.
        //  Objective is optimised by minimisation of the scalar
        //  bool indicates successful evaluation
        virtual Tuple2<scalar, bool> operator()
        (
            const scalarField& xv
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
