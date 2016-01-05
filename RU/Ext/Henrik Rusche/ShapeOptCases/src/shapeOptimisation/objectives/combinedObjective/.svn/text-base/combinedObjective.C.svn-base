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

#include "combinedObjective.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(combinedObjective, 0);
    addToRunTimeSelectionTable(objective, combinedObjective, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::combinedObjective::combinedObjective
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    objective(mesh, dict),
    objectives_(),
    weights_(dict.lookup("weights"))
{
    // Read combined objective
    PtrList<dictionary> objectiveDicts(objProperties().lookup("objectives"));

    if (weights_.size() != objectiveDicts.size())
    {
        FatalErrorIn
        (
            "combinedObjective::combinedObjective\n"
            "(\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "List of objectives and weights not of same size.  "
            << "objectives = " << objectiveDicts.size()
            << " weights = " << weights_.size()
            << abort(FatalError);
    }

    // Read objectives
    objectives_.setSize(objectiveDicts.size());

    forAll (objectiveDicts, objI)
    {
        objectives_.set
        (
            objI,
            objective::New(mesh, objectiveDicts[objI])
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::combinedObjective::evaluate() const
{
    scalarField co(objectives_.size(), 0);

    forAll (objectives_, objI)
    {
        co[objI] = objectives_[objI].evaluate();
    }

    scalar weightedObjective = sum(co*weights_);

    Info<< "objectives: " << co
        << " combined objective = " << weightedObjective << endl;

    return weightedObjective;
}


// ************************************************************************* //
