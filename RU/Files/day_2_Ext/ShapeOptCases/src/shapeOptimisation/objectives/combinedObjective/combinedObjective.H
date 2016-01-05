/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010 Hrvoje Jasak
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
    combinedObjective

Description
    Combined objective as a weighted average of objectives

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    combinedObjective.C

\*---------------------------------------------------------------------------*/

#ifndef combinedObjective_H
#define combinedObjective_H

#include "objective.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class combinedObjective Declaration
\*---------------------------------------------------------------------------*/

class combinedObjective
:
    public objective
{
    // Private data

        //- List of objectives
        PtrList<objective> objectives_;

        //- Objective weights
        scalarField weights_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        combinedObjective(const combinedObjective&);

        //- Disallow default bitwise assignment
        void operator=(const combinedObjective&);


public:

    //- Runtime type information
    TypeName("combinedObjective");


    // Constructors

        //- Construct from dictionary
        combinedObjective
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~combinedObjective()
        {}


    // Member Functions

        //- Evaluate the objective
        virtual scalar evaluate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
