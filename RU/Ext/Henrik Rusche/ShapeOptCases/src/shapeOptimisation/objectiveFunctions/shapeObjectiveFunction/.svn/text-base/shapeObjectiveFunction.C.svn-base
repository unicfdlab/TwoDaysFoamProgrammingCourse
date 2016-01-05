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

#include "shapeObjectiveFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(shapeObjectiveFunction, 0);
    addToRunTimeSelectionTable
    (
        objectiveFunction,
        shapeObjectiveFunction,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::shapeObjectiveFunction::shapeObjectiveFunction
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    objectiveFunction(mesh, dict),
    morph_(mesh),
    // Note: read subdictionary
    pointParametrisation_(functionProperties().lookup("pointParametrisation")),
    lineParametrisation_(functionProperties().lookup("lineParametrisation")),
    flowPtr_(flowModel::New(mesh, functionProperties().subDict("flowModel"))),
    objectivePtr_
    (
        objective::New(mesh, functionProperties().subDict("objective"))
    ),
     // Simulation control data
    maxIter_(readLabel(functionProperties().lookup("maxIter"))),
    objectiveTol_(readScalar(functionProperties().lookup("objectiveTol"))),
    objectiveSpan_(readScalar(functionProperties().lookup("objectiveSpan"))),
    configOffset_(readLabel(functionProperties().lookup("configOffset"))),
    configIndex_(0)
{
    // Check tolerance
    if (objectiveTol_ < SMALL)
    {
        FatalErrorIn
        (
            "shapeObjectiveFunction::shapeObjectiveFunction\n"
            "(\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "Invalid objective tolerance = " << objectiveTol_
            << abort(FatalError);
    }

    // Check iteration counters
    if (objectiveSpan_ > maxIter_)
    {
        InfoIn
        (
            "shapeObjectiveFunction::shapeObjectiveFunction\n"
            "(\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "objectiveSpan is smaller than maxIter.  "
            << "Objective convergence shall not be checked."
            << endl;
    }


    // Check parametrisation

    if (pointParametrisation_.size() + lineParametrisation_.size() < 1)
    {
        FatalErrorIn
        (
            "shapeObjectiveFunction::shapeObjectiveFunction\n"
            "(\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "Empty parametrisation.  This is not allowed"
            << abort(FatalError);

    }

    // Check parametrisation for double entries
    boolList check(morph_.controlPoints().size(), false);

    forAll (pointParametrisation_, parI)
    {
        const labelList& curSet = pointParametrisation_[parI];

        if (min(curSet) < 0 || max(curSet) >= morph_.controlPoints().size())
        {
            FatalErrorIn
            (
                "shapeObjectiveFunction::shapeObjectiveFunction\n"
                "(\n"
                "    const fvMesh& mesh,\n"
                "    const dictionary& dict\n"
                ")\n"
            )   << "PointParametrisation out or range.  All entries should be "
                << "between zero and number of control points "
                << "without duplication" << nl
                << "Number of control points = "
                << morph_.controlPoints().size() << nl
                << "Offending set: " << curSet
                << abort(FatalError);
        }

        forAll (curSet, i)
        {
            if (check[curSet[i]] == true)
            {
                FatalErrorIn
                (
                    "shapeObjectiveFunction::shapeObjectiveFunction\n"
                    "(\n"
                    "    const fvMesh& mesh,\n"
                    "    const dictionary& dict\n"
                    ")\n"
                )   << "Duplicate point in pointParametrisation: " << curSet[i]
                    << abort(FatalError);
            }
            else
            {
                check[curSet[i]] = true;
            }
        }
    }

    forAll (lineParametrisation_, parI)
    {
        const labelList& curSet = lineParametrisation_[parI];

        if (min(curSet) < 0 || max(curSet) >= morph_.controlPoints().size())
        {
            FatalErrorIn
            (
                "shapeObjectiveFunction::shapeObjectiveFunction\n"
                "(\n"
                "    const fvMesh& mesh,\n"
                "    const dictionary& dict\n"
                ")\n"
            )   << "LineParametrisation out or range.  All entries should be "
                << "between zero and number of control points "
                << "without duplication" << nl
                << "Number of control points = "
                << morph_.controlPoints().size() << nl
                << "Offending set: " << curSet
                << abort(FatalError);
        }

        forAll (curSet, i)
        {
            if (check[curSet[i]] == true)
            {
                FatalErrorIn
                (
                    "shapeObjectiveFunction::shapeObjectiveFunction\n"
                    "(\n"
                    "    const fvMesh& mesh,\n"
                    "    const dictionary& dict\n"
                    ")\n"
                )   << "Duplicate point in lineParametrisation: " << curSet[i]
                    << abort(FatalError);
            }
            else
            {
                check[curSet[i]] = true;
            }
        }
    }

    Info<< "Paramatrisation check OK.  Size = " << nArgs() << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::shapeObjectiveFunction::nArgs() const
{
    return 3*pointParametrisation_.size() + lineParametrisation_.size();
}


Foam::Tuple2<Foam::scalar, bool>
Foam::shapeObjectiveFunction::operator()
(
    const scalarField& xv
)
{
    // Receive control motion and test
    Info<< "Evaluating shapeObjectiveFunction with controls = " << xv;

    if (xv.size() != nArgs())
    {
        FatalErrorIn
        (
            "shapeObjectiveFunction::operator()(const scalarField& xv)"
        )   << "Wrong size of controls: " << xv.size()
            << ".  Should be " << nArgs()
            << abort(FatalError);
    }

    // Reject parametrisation out of range
    if (min(xv) < 0 || max(xv) > 1)
    {
        // Parametrisation out of range
        Info<< "  Rejected: out of range." << endl;

        return Tuple2<scalar, bool>(0, false);
    }
    else
    {
        Info<< endl;
    }

    // Get non-constant access to mesh and runTime
    fvMesh& m = const_cast<fvMesh&>(mesh());
    Time& runTime = const_cast<Time&>(mesh().time());


    // Save state and solution
    if (configIndex_ > 0)
    {
        // Set time to offset and save data
        Info<< "Resetting time for configuration " << configIndex_
            << " to " << configOffset_ + configIndex_ << " for data dump";
        runTime.setTime(configOffset_ + configIndex_, configIndex_);
        runTime.writeNow();
        Info<< "... done" << endl;
    }

    configIndex_++;

    vectorField controlMotion(morph_.controlPoints().size(), vector::zero);

    label xvI = 0;

    // In point parametrisation, each vector component moves independently
    forAll (pointParametrisation_, parI)
    {
        const labelList& curSet = pointParametrisation_[parI];

        vector value = vector(xv[xvI], xv[xvI + 1], xv[xvI + 2]);
        xvI += 3;

        forAll (curSet, i)
        {
            controlMotion[curSet[i]] = value;
        }
    }

    forAll (lineParametrisation_, parI)
    {
        const labelList& curSet = lineParametrisation_[parI];

        vector value = xv[xvI]*vector::one;
        xvI++;

        forAll (curSet, i)
        {
            controlMotion[curSet[i]] = value;
        }
    }

    // Deform the mesh and reset motion, mesh and database
    {
        runTime.setTime(0, 0);
        m.movePoints(morph_.motion(controlMotion)());
        m.resetMotion();
        m.moving(false);
        m.checkMesh(true);

        scalarField objList(maxIter_, 0);

        Info<< "Iterating flow model" << endl;
        for (label iter = 0; iter < maxIter_; iter++)
        {
            // Evolve the flow model
            flowPtr_->evolve();

            // Record the objective
            objList[iter] = objectivePtr_->evaluate();

            // If there is sufficient span, check objective for convergence
            if (iter > objectiveSpan_)
            {
                scalar minInRange = GREAT;
                scalar maxInRange = -GREAT;

                for (label i = iter - objectiveSpan_; i < iter; i++)
                {
                    minInRange = min(minInRange, objList[i]);
                    maxInRange = max(maxInRange, objList[i]);
                }

                Info<< "Convergence of objective: ("
                    << minInRange << ", " << maxInRange
                    << ") = " << maxInRange - minInRange << ".";

                if (maxInRange - minInRange < objectiveTol_)
                {
                    Info<< "  Converged!" << endl;
                    break;
                }
                else
                {
                    Info<< endl;
                }
            }
        }
    }

    // Evaluate the objective
    scalar value = objectivePtr_->evaluate();
    Info << "objective value = " << value << endl;

    return Tuple2<scalar, bool>(value, true);
}


// ************************************************************************* //
