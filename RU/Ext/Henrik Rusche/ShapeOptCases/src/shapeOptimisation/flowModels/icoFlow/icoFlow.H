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
    icoFlow

Description
    icoFoam flow model

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    icoFlow.C

\*---------------------------------------------------------------------------*/

#ifndef icoFlow_H
#define icoFlow_H

#include "flowModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace flowModels
{

/*---------------------------------------------------------------------------*\
                        Class icoFlow Declaration
\*---------------------------------------------------------------------------*/

class icoFlow
:
    public flowModel
{
    // Private data

        //- Velocity field
        volVectorField U_;

        //- Pressure field
        volScalarField p_;

        //- Flux field
        surfaceScalarField phi_;

        //- Kinematic viscosity
        dimensionedScalar nu_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        icoFlow(const icoFlow&);

        //- Disallow default bitwise assignment
        void operator=(const icoFlow&);


public:

    //- Runtime type information
    TypeName("icoFlow");


    // Constructors

        //- Construct from components
        icoFlow
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~icoFlow()
        {}


    // Member Functions

        // Edit

            //- Evolve the flow model
            virtual void evolve();

            //- Has the flow model converged?
            virtual bool converged() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace flowModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
