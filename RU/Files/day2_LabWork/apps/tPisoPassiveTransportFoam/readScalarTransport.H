IOdictionary ptDict
(
    IOobject
    (
	"passiveScalarTransportProperties",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    )
);

label nScalars = 0;

ptDict.lookup("nScalars") >> nScalars;

List<scalar> DiffCoeffs(nScalars);
ptDict.lookup("DC") >> DiffCoeffs;
forAll(DiffCoeffs, iCoeff)
{
    Info << "DC[" << iCoeff << "]=" << DiffCoeffs[iCoeff] << endl;
}

List<autoPtr<volScalarField> > CiPtr (nScalars);


forAll(CiPtr, iPtr)
{
    word nameCi ("C" + name(iPtr));

    CiPtr[iPtr].reset
    (
	new volScalarField
        (
	    IOobject
	    (
		nameCi,
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh
	)
    );
}
