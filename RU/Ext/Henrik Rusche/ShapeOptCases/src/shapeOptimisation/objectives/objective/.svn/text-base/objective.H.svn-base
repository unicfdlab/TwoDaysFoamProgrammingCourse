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

Class
    objective

Description
    Virtual base class for data extraction objects from a simulation:
    This encodes the objective of optimisation, eg. minimum drag

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    objective.C
    newObjective.C

\*---------------------------------------------------------------------------*/

#ifndef objective_H
#define objective_H

#include "fvMesh.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class objective Declaration
\*---------------------------------------------------------------------------*/

class objective
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Objective properties dictionary
        dictionary objProperties_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        objective(const objective&);

        //- Disallow default bitwise assignment
        void operator=(const objective&);


protected:

    // Protected member functions

        //- Return mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return objective properties dictionary
        const dictionary& objProperties() const
        {
            return objProperties_;
        }


public:

    //- Runtime type information
    TypeName("objective");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            objective,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Constructors

        //- Construct from components
        objective
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Selectors

        //- Select constructed from mesh
        static autoPtr<objective> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~objective();


    // Member Functions

        //- Evaluate the objective
        virtual scalar evaluate() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
