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
    SimplexNelderMead

Description
    Simplex Algorithmn for minimisation based on GSL library.
    Function is provided as a template parameter function object, evaluated
    using operator()(const scalarField x) returning a Tuple2<scalar, bool>

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved.

SourceFiles
    SimplexNelderMead.C

\*---------------------------------------------------------------------------*/

#ifndef SimplexNelderMead_H
#define SimplexNelderMead_H

#include "PtrList.H"
#include "scalarField.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class SimplexNelderMead Declaration
\*---------------------------------------------------------------------------*/


template<class Func>
class SimplexNelderMead
{
    // Private data types

        class simplexCorner
        {
            //- Coordinates for evaluation
            scalarField coord_;

            //- Value at coordinates
            scalar value_;


        public:

            //- Construct from coordinates and value
            simplexCorner(const scalarField& coord, scalar value)
            :
                coord_(coord),
                value_(value)
            {}

            //- Reset coordinates and value
            void reset(const scalarField& coord, scalar value)
            {
                coord_ = coord;
                value_ = value;
            }

            //- Return coordinates
            const scalarField& coord() const
            {
                return coord_;
            }

            //- Return value
            scalar value() const
            {
                return value_;
            }
        };


    // Private data

        //- Reference to a function
        Func& f_;

        //- Points of the simplex;
        PtrList<simplexCorner> points_;

        //- Creates a new point by moving a simplex corner scaled by coeff
        //  (negative value represents mirroring by the middle point of
        //  the "other" corner points)
        autoPtr<simplexCorner> newFromWorst
        (
            const scalar coeff
        );

        //- Calculate the center of the simplex
        scalarField center() const;

        //- Contract the simplex in respect to best valued corner
        void contractByBest();

        //- Swap two points of the simplex
        void swapPoints(label i, label j);

        //- Sort the points partially by value, so that the smallest value
        //  comes first; seconds largest at nArgs()-1; largest at nArgs()
        void updateOrdering();


public:

    // Constructors

        //- Construct given a function, a reference point and range of lambda
        SimplexNelderMead
        (
            Func& f,
            const scalarField& p0,
            const scalarField& lambda
        );


    // Destructor - default


    // Member Functions

        //- Return the current minimum
        scalar min() const
        {
            return points_[0].value();
        };

        //- Return the current coordinates of the minimum
        const scalarField& minCoord() const
        {
            return points_[0].coord();
        }

        //- Make a simplex step
        void iterate();

        //- return size of the simplex
        scalar size();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SimplexNelderMead.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
