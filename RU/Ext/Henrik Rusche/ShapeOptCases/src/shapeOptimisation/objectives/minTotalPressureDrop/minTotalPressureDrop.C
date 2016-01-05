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

#include "minTotalPressureDrop.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(minTotalPressureDrop, 0);
    addToRunTimeSelectionTable(objective, minTotalPressureDrop, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::minTotalPressureDrop::checkPatchNames() const
{
    if
    (
        mesh().boundaryMesh().findPatchID(inletPatchName_) < 0
     || mesh().boundaryMesh().findPatchID(outletPatchName_) < 0
    )
    {
        FatalErrorIn("void minTotalPressureDrop::checkPatchNames() const")
            << "Patch names (" << inletPatchName_ << ", " << outletPatchName_
            << ") not found." << nl
            << "Available patches are: " << mesh().boundaryMesh().names()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::minTotalPressureDrop::minTotalPressureDrop
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    objective(mesh, dict),
    // Note: read from local dictionary
    pName_(objProperties().lookup("p")),
    UName_(objProperties().lookup("U")),
    inletPatchName_(objProperties().lookup("inlet")),
    outletPatchName_(objProperties().lookup("outlet")),
    rhoRef_(objProperties().lookup("rhoRef"))
{
    checkPatchNames();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::minTotalPressureDrop::evaluate() const
{
    const volScalarField& p = mesh().lookupObject<volScalarField>(pName_);
    const volVectorField& U = mesh().lookupObject<volVectorField>(UName_);


    const label inletIndex =
        mesh().boundaryMesh().findPatchID(inletPatchName_);

    const label outletIndex =
        mesh().boundaryMesh().findPatchID(outletPatchName_);

    scalar pTotAvgIn = rhoRef_.value()*
        gAverage
        (
            p.boundaryField()[inletIndex]
          + 0.5*magSqr(U.boundaryField()[inletIndex])
        );

    scalar pTotAvgOut = rhoRef_.value()*
        gAverage
        (
            p.boundaryField()[outletIndex]
          + 0.5*magSqr(U.boundaryField()[outletIndex])
        );

    Info<< "Total pressure drop: (" << pTotAvgIn << ", " << pTotAvgOut << ")"
        << endl;

    return pTotAvgIn - pTotAvgOut;
}


// ************************************************************************* //
