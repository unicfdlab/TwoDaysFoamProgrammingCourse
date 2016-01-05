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
    RBFMeshMorph

Description
    Radial basis function mesh morphing tool

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    RBFMeshMorph.C

\*---------------------------------------------------------------------------*/

#ifndef RBFMeshMorph_H
#define RBFMeshMorph_H

#include "IOdictionary.H"
#include "polyMesh.H"
#include "Tuple2.H"
#include "RBFInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class RBFMeshMorph Declaration
\*---------------------------------------------------------------------------*/

class RBFMeshMorph
:
    public IOdictionary
{
    // Private data

        //- Mesh reference
        const polyMesh& mesh_;

        //- Moving patch names
        wordList movingPatches_;

        //- Control points governing mesh morphing
        vectorField controlPoints_;

        //- Motion bounds
        List<Tuple2<vector, vector> > motionBounds_;

        //- Motion mask
        scalarField motionMask_;

        //- Reference points
        vectorField referencePoints_;

        //- RBF interpolation
        RBFInterpolation interpolation_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        RBFMeshMorph(const RBFMeshMorph&);

        //- Disallow default bitwise assignment
        void operator=(const RBFMeshMorph&);

        //- Make control masks
        void makeControlMasks();


public:

    //- Runtime type information
    TypeName("RBFMeshMorph");


    // Constructors

        //- Construct from polyMesh
        RBFMeshMorph(const polyMesh& mesh);


    // Destructor

        virtual ~RBFMeshMorph();


    // Member Functions

        //- Return control points governing mesh morphing
        const vectorField& controlPoints() const
        {
            return controlPoints_;
        }

        //- Return motion of points given motion parameters (0 <= cpm <= 1)
        tmp<pointField> motion(const scalarField& cpm) const;

        //- Return motion of points given motion parameters (0 <= cpm <= 1)
        //  Used for independent parametrisation of vector directions
        tmp<pointField> motion(const vectorField& cpm) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
