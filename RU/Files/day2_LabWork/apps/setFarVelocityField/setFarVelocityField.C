/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    setFarVelocityField

Description
    Sets far velocity field to fixed value except thin region around selected object

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar minPatchDist (vector C, const fvMesh& mesh, label patchId)
{
    vectorField Cf = mesh.Cf().boundaryField()[patchId];
    scalar dist = 1 / VSMALL;
    forAll (Cf, i)
    {
	scalar cDist = mag(C - Cf[i]);
	if (cDist < dist)
	{
	    dist = cDist;
	}
    }

    return dist;
}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    
    // set startTime and endTime depending on -time and -latestTime options
    instantList Times = runTime.times();
    #include "checkTimeOptions.H"

    // Set time to start time
    runTime.setTime(Times[startTime], startTime);
    
    for (label i=startTime; i<endTime; i++)
    {
        //set time to current time dir
        runTime.setTime(Times[i], i);
        //update mesh
        polyMesh::readUpdateState state = mesh.readUpdate();
        
	//set volume fields
        forAll (mesh.C(), i)
        {
	    scalar dist = minPatchDist (mesh.C()[i], mesh, wallPatchId);
	    if (dist > minDist)
	    {
		U[i] = Ufar;
	    }
	    else
	    {
		U[i] = Unear;
	    }
        }
        
        Info << "Writing velocity field for Time = " << runTime.timeName() << endl;

        //write results
        U.write();
        //runTime.write();
    }

    Info<< "end" << endl;

    return 0;
}


// ************************************************************************* //
