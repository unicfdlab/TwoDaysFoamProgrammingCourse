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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    objectiveFunction

Description
    Objective function used in optimisation.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    objectiveFunction.C
    newObjectiveFunction.C

\*---------------------------------------------------------------------------*/

#ifndef objectiveFunction_H
#define objectiveFunction_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class objectiveFunction Declaration
\*---------------------------------------------------------------------------*/

class objectiveFunction
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Objective function properties dictionary
        dictionary functionProperties_;


    // Private Member Functions

        //- Disallow copy construct
        objectiveFunction(const objectiveFunction&);

        //- Disallow default bitwise assignment
        void operator=(const objectiveFunction&);


protected:

    // Protected member functions

        //- Return reference to mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return objective function properties dictionary
        const dictionary& functionProperties() const
        {
            return functionProperties_;
        }


public:

    //- Runtime type information
    TypeName("objectiveFunction");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            objectiveFunction,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Selectors

        //- Return a reference to the selected thermal model
        static autoPtr<objectiveFunction> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        objectiveFunction
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~objectiveFunction()
        {}


    // Member Functions

        //- Return number of arguments
        virtual label nArgs() const = 0;

        //- Evaluate and return objective.
        //  Objective is optimised by minimisation of the scalar
        //  bool indicates successful evaluation
        virtual Tuple2<scalar, bool> operator()
        (
            const scalarField& xv
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
