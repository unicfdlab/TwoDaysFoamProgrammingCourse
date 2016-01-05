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

Application
    tPisoFoam

Description
    Transient solver for incompressible non-isothermal flow.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readPISOControls.H"
#       include "CourantNo.H"

        // Pressure-velocity PISO corrector
        {
            // Momentum predictor
	    
	    surfaceScalarField muEff
	    (
		"muEff",
		fvc::interpolate
		(
		    rho*
		    (
			laminarTransport.nu()
			+ turbulence->nut()
		    )
		)
	    );

            fvVectorMatrix UEqn
            (
                fvm::ddt(rho, U)
              + fvm::div(rhoPhi, U)
	      - fvm::laplacian(muEff, U)
	      - (fvc::grad(U) & fvc::grad(muEff))
            );

            UEqn.relax();

            if (momentumPredictor)
            {
                solve(UEqn == -fvc::grad(p));
            }

            // --- PISO loop

            for (int corr=0; corr<nCorr; corr++)
            {
                volScalarField rUA = 1.0/UEqn.A();
                surfaceScalarField rhorAUf = fvc::interpolate(rho*rUA);

                U = rUA*UEqn.H();
                rhoPhi = fvc::interpolate(rho)*
			(
			       (fvc::interpolate(U) & mesh.Sf())
			     + fvc::ddtPhiCorr(rUA, rho, U, phi)
			);

                adjustPhi(rhoPhi, U, p);

                // Non-orthogonal pressure corrector loop
                for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                {
                    // Pressure corrector
                    fvScalarMatrix pEqn
                    (
                	  fvc::ddt(rho)
                	  +
                	  fvc::div(rhoPhi)
                	  -
                	  fvm::laplacian(rho*rUA, p)
                    );

                    pEqn.setReference(pRefCell, pRefValue);

                    if
                    (
                        corr == nCorr-1
                     && nonOrth == nNonOrthCorr
                    )
                    {
                        pEqn.solve(mesh.solver("pFinal"));
                    }
                    else
                    {
                        pEqn.solve();
                    }

                    if (nonOrth == nNonOrthCorr)
                    {
                        rhoPhi += pEqn.flux();
                        phi = rhoPhi / fvc::interpolate(rho);
                    }
                }

#               include "continuityErrs.H"
		
		U = rUA*fvc::reconstruct(rhoPhi/rhorAUf);
                //U -= rUA*fvc::grad(p);
                U.correctBoundaryConditions();
            }
        }
        
        #include "solvePassiveTransport.H"

	volScalarField kappaEff = lambda / Cp + rho*turbulence->nut();
	kappaEff.rename("kappaEff");
	fvScalarMatrix TEqn
	(
	      fvm::ddt(rho,T)
	    + fvm::div(rhoPhi,T)
	    - fvm::laplacian(kappaEff,T)
	);
	
	TEqn.solve();
	
	rho = rho0 - (beta * rho) *(T-T0);

        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
