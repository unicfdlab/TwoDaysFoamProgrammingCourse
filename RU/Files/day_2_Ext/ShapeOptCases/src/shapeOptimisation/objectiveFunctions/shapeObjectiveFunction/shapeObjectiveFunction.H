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
    shapeObjectiveFunction

Description
    Objective function which optimises a shape in a flow using RBF morphing
    and target property to minimise

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    shapeObjectiveFunction.C

\*---------------------------------------------------------------------------*/

#ifndef shapeObjectiveFunction_H
#define shapeObjectiveFunction_H

#include "objectiveFunction.H"
#include "RBFMeshMorph.H"
#include "flowModel.H"
#include "objective.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class shapeObjectiveFunction Declaration
\*---------------------------------------------------------------------------*/

class shapeObjectiveFunction
:
    public objectiveFunction
{
    // Private data

        //- Mesh morphing object
        RBFMeshMorph morph_;

        //- Point parametrisation of the mesh morphing independently
        //  in (x y z) directions.  Note that each point carries 3
        //  free parameters
        //  Indices of moving points are grouped under parameter sets
        // that are being changed in unison
        labelListList pointParametrisation_;

        //- Line parametrisation of the mesh morphing along lines
        //  Indices of moving points are grouped under parameter sets
        // that are being changed in unison
        labelListList lineParametrisation_;

        //- Flow model
        autoPtr<flowModel> flowPtr_;

        //- Minimisation function: this is the object of shape optimisation
        autoPtr<objective> objectivePtr_;


        // Simulation control data

            //- Max number of iterations to flow solver convergence
            label maxIter_;

            //- Convergence tolerance for the objective
            scalar objectiveTol_;

            //- Iteration span for objective tolerance check
            //  If the objective does not change by more than objectiveTol_
            //  over the last objectiveSpan_ iterations, the objective is
            //  considered to be converged
            label objectiveSpan_;

            //- Configuration offset.  Where to store current configuration
            const label configOffset_;


        // State data

            //- Configuration index.  Incremented by evaluation call from
            //  the optimiser
            label configIndex_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        shapeObjectiveFunction(const shapeObjectiveFunction&);

        //- Disallow default bitwise assignment
        void operator=(const shapeObjectiveFunction&);


public:

    //- Runtime type information
    TypeName("shapeObjective");

    // Static data members


    // Constructors

        //- Construct from dictionary
        shapeObjectiveFunction
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~shapeObjectiveFunction()
        {}


    // Member Functions

        //- Return number of arguments
        virtual label nArgs() const;

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
