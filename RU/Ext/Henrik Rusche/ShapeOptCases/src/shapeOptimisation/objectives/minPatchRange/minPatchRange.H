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
    minPatchRange

Description
    Minimise range if a variable on a patch: optimisation objective

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    minPatchRange.C

\*---------------------------------------------------------------------------*/

#ifndef minPatchRange_H
#define minPatchRange_H

#include "objective.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class minPatchRange Declaration
\*---------------------------------------------------------------------------*/

class minPatchRange
:
    public objective
{
    // Private data

        //- Name of field
        word fieldName_;

        //- List of patch names to integrate
        wordList patchNames_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        minPatchRange(const minPatchRange&);

        //- Disallow default bitwise assignment
        void operator=(const minPatchRange&);


        //- Check patch names
        void checkPatchNames() const;

public:

    //- Runtime type information
    TypeName("minPatchRange");

    // Static data members


    // Constructors

        //- Construct from dictionary
        minPatchRange
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~minPatchRange()
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
