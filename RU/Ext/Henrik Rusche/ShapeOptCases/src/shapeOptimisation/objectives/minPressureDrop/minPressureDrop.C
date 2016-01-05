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

#include "minPressureDrop.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(minPressureDrop, 0);
    addToRunTimeSelectionTable(objective, minPressureDrop, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::minPressureDrop::checkPatchNames() const
{
    if
    (
        mesh().boundaryMesh().findPatchID(inletPatchName_) < 0
     || mesh().boundaryMesh().findPatchID(outletPatchName_) < 0
    )
    {
        FatalErrorIn("void minPressureDrop::checkPatchNames() const")
            << "Patch names (" << inletPatchName_ << ", " << outletPatchName_
            << ") not found." << nl
            << "Available patches are: " << mesh().boundaryMesh().names()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::minPressureDrop::minPressureDrop
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    objective(mesh, dict),
    // Note: read from local dictionary
    pName_(objProperties().lookup("p")),
    inletPatchName_(objProperties().lookup("inlet")),
    outletPatchName_(objProperties().lookup("outlet"))
{
    checkPatchNames();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::minPressureDrop::evaluate() const
{
    const volScalarField& p = mesh().lookupObject<volScalarField>(pName_);


    const label inletIndex =
        mesh().boundaryMesh().findPatchID(inletPatchName_);

    const label outletIndex =
        mesh().boundaryMesh().findPatchID(outletPatchName_);

    scalar pAvgIn = gAverage(p.boundaryField()[inletIndex]);
    scalar pAvgOut = gAverage(p.boundaryField()[outletIndex]);

    Info<< "Pressure drop: (" << pAvgIn << ", " << pAvgOut << ")" << endl;

    return pAvgOut - pAvgIn;
}


// ************************************************************************* //
