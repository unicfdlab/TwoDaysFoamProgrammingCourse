
    // Read optimiserDict
    IOdictionary optimiserDict
    (
        IOobject
        (
            "optimiserDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Create objective function
    autoPtr<objectiveFunction> objective = objectiveFunction::New
    (
        mesh,
        optimiserDict.subDict("objectiveFunction")
    );

    // Read start point and lambda range
    dictionary simplexDict = optimiserDict.subDict("simplex");

    scalarField startPoint(simplexDict.lookup("startPoint"));
    scalarField lambda(simplexDict.lookup("lambda"));

    label maxIter(readLabel(simplexDict.lookup("maxIter")));
    scalar tolerance(readScalar(simplexDict.lookup("tolerance")));

    // Create simplex optimiser
    SimplexNelderMead<objectiveFunction> simplex
    (
        objective(),
        startPoint,
        lambda
    );
