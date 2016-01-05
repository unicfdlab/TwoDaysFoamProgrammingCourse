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

\*---------------------------------------------------------------------------*/

#include "minPatchRange.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(minPatchRange, 0);
    addToRunTimeSelectionTable(objective, minPatchRange, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::minPatchRange::checkPatchNames() const
{
    forAll (patchNames_, pnI)
    {
        if (mesh().boundaryMesh().findPatchID(patchNames_[pnI]) < 0)
        {
            FatalErrorIn("void minPatchRange::checkPatchNames() const")
                << "Patch names " << patchNames_[pnI] << " not found." << nl
                << "Available patches are: " << mesh().boundaryMesh().names()
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::minPatchRange::minPatchRange
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    objective(mesh, dict),
    // Note: read from local dictionary
    fieldName_(objProperties().lookup("field")),
    patchNames_(objProperties().lookup("patches"))
{
    checkPatchNames();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::minPatchRange::evaluate() const
{
    scalar maxField = -GREAT;
    scalar minField = GREAT;

    // Lookup the variable as flux
    if (mesh().foundObject<surfaceScalarField>(fieldName_))
    {
        const surfaceScalarField& field =
            mesh().lookupObject<surfaceScalarField>(fieldName_);

        forAll (patchNames_, pnI)
        {
            const label patchIndex =
                mesh().boundaryMesh().findPatchID(patchNames_[pnI]);

            minField = Foam::min
            (
                minField,
                gMin
                (
                    field.boundaryField()[patchIndex]/
                    field.boundaryField()[patchIndex].patch().magSf()
                )
            );

            maxField = Foam::max
            (
                maxField,
                gMax
                (
                    field.boundaryField()[patchIndex]/
                    field.boundaryField()[patchIndex].patch().magSf()
                )
            );
        }
    }
    // Lookup the variable as a scalar field
    else if (mesh().foundObject<volScalarField>(fieldName_))
    {
        const volScalarField& field =
            mesh().lookupObject<volScalarField>(fieldName_);

        forAll (patchNames_, pnI)
        {
            const label patchIndex =
                mesh().boundaryMesh().findPatchID(patchNames_[pnI]);

            if (field.boundaryField()[patchIndex].size() > 0)
            {
                minField = Foam::min
                (
                    minField,
                    min(field.boundaryField()[patchIndex])
                );

                maxField = Foam::max
                (
                    maxField,
                    max(field.boundaryField()[patchIndex])
                );
            }
        }

        reduce(minField, minOp<scalar>());
        reduce(maxField, maxOp<scalar>());
    }
    else
    {
        FatalErrorIn("scalar minPatchRange::evaluate() const")
            << "Cannot find field " << fieldName_
            << " as flux or scalar field"
            << abort(FatalError);
    }

    Info<< "field range: (" << minField << ", " << maxField << ")" << endl;

    return maxField - minField;
}


// ************************************************************************* //
