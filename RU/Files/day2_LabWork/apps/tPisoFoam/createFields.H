    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

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

#   include "createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    dimensionedScalar rho0   ("rho0", dimensionSet(1, -3, 0, 0, 0), 1000);
    dimensionedScalar beta   ("beta", dimensionSet(1, -3, 0, -1, 0), -1.0);
    dimensionedScalar T0     ("T0", dimensionSet(0, 0, 0, 1, 0), 293.15);
    dimensionedScalar lambda ("lambda", dimensionSet(1, 1, -3, -1, 0), 0.6);
    dimensionedScalar Cp     ("Cp", dimensionSet(0, 2, -2, -1, 0), 4.5e+3);
    dimensionedScalar Prt    ("Prt",dimensionSet(0, 0, 0, 0, 0), 1.0);

    dictionary transportProperties = mesh.thisDb().lookupObject<dictionary>("transportProperties");

    transportProperties.lookup("rho0") >> rho0;
    transportProperties.lookup("beta") >> beta;
    transportProperties.lookup("T0") >> T0;
    transportProperties.lookup("lambda") >> lambda;
    transportProperties.lookup("Cp") >> Cp;
    transportProperties.lookup("Prt") >> Prt;

    volScalarField rho
    (
	IOobject
	(
	    "rho",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
	rho0 - beta*rho0*(T - T0)
    );

    surfaceScalarField rhoPhi = fvc::interpolate(rho)*phi;
    rhoPhi.rename("rho*phi");

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

