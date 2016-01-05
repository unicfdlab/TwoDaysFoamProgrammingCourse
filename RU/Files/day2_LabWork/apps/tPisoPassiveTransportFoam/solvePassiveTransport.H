{
    forAll(CiPtr, iPtr)
    {
	volScalarField& Ci = CiPtr[iPtr]();
	
	dimensionedScalar ci ("Ci", dimensionSet(0, 2, -1, 0, 0, 0, 0), DiffCoeffs[iPtr]);
	volScalarField dc = rho*ci;
	
	solve
	(
	    fvm::ddt(rho, Ci)
	    +
	    fvm::div(rhoPhi,Ci)
	    - 
	    fvm::laplacian(dc,Ci)
        );
    }
}

