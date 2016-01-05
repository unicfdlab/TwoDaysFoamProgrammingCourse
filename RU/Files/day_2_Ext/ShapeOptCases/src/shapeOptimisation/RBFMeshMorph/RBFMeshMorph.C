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

\*---------------------------------------------------------------------------*/

#include "RBFMeshMorph.H"
#include "Time.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RBFMeshMorph, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::RBFMeshMorph::makeControlMasks()
{
    // Check definition of moving patches
    forAll (movingPatches_, patchI)
    {
        // Find the patch in boundary
        label patchIndex =
            mesh_.boundaryMesh().findPatchID(movingPatches_[patchI]);

        if (patchIndex < 0)
        {
            FatalErrorIn("void RBFMeshMorph::makeControlMasks()")
                << "Patch " << movingPatches_[patchI] << " not found.  "
                << "valid patch names: " << mesh_.boundaryMesh().names()
                << abort(FatalError);
        }
    }

    // Find static patches and fix them by setting motion mask to zero
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label nMarkedPoints = 0;

    forAll (patches, patchI)
    {
        const word& patchName = patches[patchI].name();

        // Find out if the patch is moving
        bool moving = false;

        forAll (movingPatches_, mpI)
        {
            if (movingPatches_[mpI] == patchName)
            {
                moving = true;
                break;
            }
        }

        // Check for empty or coupled patch type: it is always moving
        if (isA<emptyPolyPatch>(patches[patchI]) || patches[patchI].coupled())
        {
            moving = true;
        }

        if (!moving)
        {
            const labelList& mp = patches[patchI].meshPoints();

            forAll (mp, i)
            {
                motionMask_[mp[i]] = 0;
                nMarkedPoints++;
            }
        }
    }

    Info << "Total points on static boundaries: " << nMarkedPoints << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBFMeshMorph::RBFMeshMorph
(
    const polyMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "rbfMotionDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    movingPatches_(lookup("movingPatches")),
    controlPoints_(lookup("controlPoints")),
    motionBounds_(lookup("motionBounds")),
    motionMask_(mesh.allPoints().size(), 1.0),
    referencePoints_(mesh.allPoints()),
    interpolation_
    (
        subDict("interpolation"),
        controlPoints_,
        referencePoints_
    )
{
    if (controlPoints_.size() != motionBounds_.size())
    {
        FatalErrorIn("RBFMeshMorph::RBFMeshMorph(const polyMesh& mesh)")
            << "Inconsistent size of controlPoints and motionBounds" << nl
            << "controlPoints = " << controlPoints_.size()
            << " motionBounds = " << motionBounds_.size()
            << abort(FatalError);
    }

    makeControlMasks();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBFMeshMorph::~RBFMeshMorph()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::RBFMeshMorph::motion(const scalarField& cpm) const
{
    if (cpm.size() != motionBounds_.size())
    {
        FatalErrorIn("RBFMeshMorph::motion(const scalarField& cpm) const")
            << "Incorrect size of motion control parameters.  "
            << "cpm: " << cpm.size() << "; should be " << motionBounds_.size()
            << abort(FatalError);
    }

    vectorField controlMotion(motionBounds_.size(), vector::zero);

    forAll (controlMotion, pointI)
    {
        controlMotion[pointI] = motionBounds_[pointI].first()
          + cpm[pointI]*
            (
                motionBounds_[pointI].second()
              - motionBounds_[pointI].first()
            );
    }

    return referencePoints_
        + motionMask_*interpolation_.interpolate(controlMotion);
}


Foam::tmp<Foam::pointField>
Foam::RBFMeshMorph::motion(const vectorField& cpm) const
{
    if (cpm.size() != motionBounds_.size())
    {
        FatalErrorIn("RBFMeshMorph::motion(const vectorField& cpm) const")
            << "Incorrect size of motion control parameters.  "
            << "cpm: " << cpm.size() << "; should be " << motionBounds_.size()
            << abort(FatalError);
    }

    vectorField controlMotion(motionBounds_.size(), vector::zero);

    forAll (controlMotion, pointI)
    {
        controlMotion[pointI] = motionBounds_[pointI].first()
          + cmptMultiply
            (
                cpm[pointI],
                (
                    motionBounds_[pointI].second()
                  - motionBounds_[pointI].first()
                )
            );
    }

    return referencePoints_
        + motionMask_*interpolation_.interpolate(controlMotion);
}


// ************************************************************************* //
