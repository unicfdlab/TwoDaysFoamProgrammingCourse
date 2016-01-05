    Info<< "Reading rbfMotionDict\n" << endl;

    RBFMeshMorph rbfMotion(mesh);

    label nSteps = readLabel(rbfMotion.lookup("nSteps"));
