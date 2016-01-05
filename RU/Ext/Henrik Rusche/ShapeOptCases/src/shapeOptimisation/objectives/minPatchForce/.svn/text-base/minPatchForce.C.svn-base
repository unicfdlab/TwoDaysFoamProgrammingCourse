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

#include "minPatchForce.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

#include "RASModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(minPatchForce, 0);
    addToRunTimeSelectionTable(objective, minPatchForce, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::minPatchForce::checkPatchNames() const
{
    forAll (patchNames_, pnI)
    {
        if (mesh().boundaryMesh().findPatchID(patchNames_[pnI]) < 0)
        {
            FatalErrorIn("void minPatchForce::checkPatchNames() const")
                << "Patch names " << patchNames_[pnI] << " not found." << nl
                << "Available patches are: " << mesh().boundaryMesh().names()
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::minPatchForce::minPatchForce
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    objective(mesh, dict),
    // Note: read from local dictionary
    patchNames_(objProperties().lookup("patches")),
    pName_(objProperties().lookupOrDefault<word>("pName", "p")),
    UName_(objProperties().lookupOrDefault<word>("UName", "U")),
    rhoRef_(readScalar(objProperties().lookup("rhoInf"))),
    useDirection_(objProperties().lookup("useDirection")),
    direction_(1, 0, 0)
{
    if (useDirection_)
    {
        direction_ = vector(objProperties().lookup("direction"));

        if (mag(direction_) > SMALL)
        {
            direction_ /= mag(direction_);
        }
        else
        {
            FatalErrorIn
            (
                "minPatchForce::minPatchForce\n"
                "(\n"
                "    const fvMesh& mesh,\n"
                "    const dictionary& dict\n"
                ")\n"
            )   << "Incorrect direction: zero magnitude.  "
                << "direction = " << direction_
                << abort(FatalError);
        }
    }

    checkPatchNames();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::minPatchForce::evaluate() const
{
    // Calculate force
    const volVectorField& U = mesh().lookupObject<volVectorField>(UName_);
    const volScalarField& p = mesh().lookupObject<volScalarField>(pName_);

    const incompressible::RASModel& ras =
        mesh().lookupObject<incompressible::RASModel>("RASProperties");

    volSymmTensorField devRhoReff = rhoRef_*ras.devReff();

    const volSymmTensorField::GeometricBoundaryField& devRhoReffb
        = devRhoReff.boundaryField();

    const surfaceVectorField::GeometricBoundaryField& Sfb =
        mesh().Sf().boundaryField();

    vector totalForce = vector::zero;

    forAll (patchNames_, pnI)
    {
        const label patchIndex =
            mesh().boundaryMesh().findPatchID(patchNames_[pnI]);

        totalForce += rhoRef_*
            gSum
            (
                Sfb[patchIndex]*p.boundaryField()[patchIndex]
              + (Sfb[patchIndex] & devRhoReffb[patchIndex])
            );
    }

    if (useDirection_)
    {
        return mag(direction_ & totalForce);
    }
    else
    {
        return mag(direction_);
    }
}


// ************************************************************************* //
