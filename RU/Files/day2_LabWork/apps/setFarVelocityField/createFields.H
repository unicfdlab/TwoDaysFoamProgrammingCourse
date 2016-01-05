    Info<< "Reading field U\n" << endl;

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    //Open dictionary
    IOdictionary setFarVelocityFieldDict
    (
        IOobject
        (
            "setFarVelocityFieldDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    //read dictionary entries
    word wallPatchName;
    label wallPatchId = -1; 
    
    vector Ufar (0, 0, 0);
    vector Unear (0, 0, 0);
    scalar minDist = 0;
    
    setFarVelocityFieldDict.lookup("wallPatchName") >> wallPatchName;
    wallPatchId = mesh.boundary().findPatchID(wallPatchName);
    if (wallPatchId < 0)
    {
	FatalError
	    << "Can\'t find patch " << wallPatchName << nl
	    << exit(FatalError);
    }
    
    if (!isA<wallFvPatch>(mesh.boundary()[wallPatchId]))
    {
	FatalError
	    << "Patch " << wallPatchName << "is not a wall, exiting" << nl
	    << exit(FatalError);
    }
    
    setFarVelocityFieldDict.lookup("Ufar") >> Ufar;
    setFarVelocityFieldDict.lookup("Unear") >> Unear;
    setFarVelocityFieldDict.lookup("minDist") >> minDist;

//END_OF_FILE

