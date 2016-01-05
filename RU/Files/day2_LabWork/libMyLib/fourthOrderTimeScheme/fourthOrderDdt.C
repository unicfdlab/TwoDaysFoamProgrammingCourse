#include "fourthOrderDdt.H"
#include "volFields.H"

Foam::fourthOrderDdt::fourthOrderDdt (const volScalarField& psi)
:
    psi_(psi),
    mesh_(psi.mesh()),
    cTime_(-1),
    dt_(4),
    oldPsi_(4),
    oldV_(4),
    delta_(4),
    k_(4)
{

    cTime_ = mesh_.time().value();

    forAll (dt_, i)
    {
	dt_[i] = mesh_.time().deltaT().value();
    }
    
    forAll (oldPsi_, i)
    {
	oldPsi_[i].resize(psi_.size());
	oldV_[i].resize(psi_.size());
	scalarField& coPsi = oldPsi_[i];
	scalarField& coV   = oldV_[i];
	forAll (coPsi, j)
	{
	    coPsi[j] = psi_[j];
	    coV[j]   = mesh_.V()[j];
	}
    }
    
    advanceInTime();
    
    makeCoeffs();
}

void Foam::fourthOrderDdt::advanceInTime()
{

    if (mesh_.time().value() <= cTime_)
    {
	return;
    }

    cTime_ = mesh_.time().value();
    
    //store old time deltas
    for (label ti=2; ti>=0; ti--)
    {
	dt_[ti] = dt_[ti+1];
    }
    
    dt_[3] = mesh_.time().deltaT().value();
    
    //store old field values and mesh volumes
    for (label ti=0; ti<=2; ti++)
    {
	oldPsi_[ti] = oldPsi_[ti+1];
	oldV_[ti] = oldV_[ti+1];
    }
    
    forAll (psi_, j)
    {
	oldPsi_[3][j] = psi_[j];
	oldV_[3][j] = mesh_.V()[j];
    }
}

void Foam::fourthOrderDdt::makeCoeffs()
{
    //build time delta sums
    forAll (delta_, i)
    {
	delta_[i] = 0.0;
	for (label j=3; j>=i; j--)
	{
	    delta_[i] += dt_[j];
	}
    }

    k_[0] = delta_[1]*delta_[2]*delta_[3] / 
    (
	(-delta_[3] + delta_[0])
	*
	delta_[0]
	*
	(-delta_[1] + delta_[0])
	*
	(-delta_[2] + delta_[0])
    );
    
    k_[1] = -
    delta_[0]*delta_[2]*delta_[3] /
    (
	(delta_[1] - delta_[3])
	*
	(-delta_[1] + delta_[0])
	*
	delta_[1]
	*
	(delta_[1] - delta_[2])
    );
    
    k_[2] = delta_[0]*delta_[3]*delta_[1] /
    (
	(delta_[2] - delta_[3])*
	(
	    delta_[1]*delta_[0]
	    -
	    delta_[2]*delta_[0]
	    +
	    delta_[2]*delta_[2]
	    -
	    delta_[2]*delta_[1]
	)
	*
	delta_[2]
    );
    
    k_[3] = -
    delta_[0]*delta_[2]*delta_[1] /
    (
	(
	    delta_[0]*delta_[2]*delta_[1]
	    -
	    delta_[0]*delta_[3]*delta_[1]
	    +
	    delta_[3]*delta_[3]*delta_[0]
	    -
	    delta_[0]*delta_[2]*delta_[3]
	    +
	    delta_[3]*delta_[3]*delta_[1]
	    -
	    delta_[1]*delta_[2]*delta_[3]
	    -
	    delta_[3]*delta_[3]*delta_[3]
	    +
	    delta_[3]*delta_[3]*delta_[2]
	)
	*
	delta_[3]
    );
}

const Foam::fvScalarMatrix& Foam::fourthOrderDdt::ddt()
{
    advanceInTime();
    
    makeCoeffs();
    
    mtrx_ = tmp<fvScalarMatrix> 
    (
        new fvScalarMatrix
        (
            psi_,
            psi_.dimensions()*dimVol / dimTime
        )
    );
    
    forAll (mtrx_().diag(), i)
    {
	mtrx_().diag()[i] = - (k_[0] + k_[1] + k_[2] + k_[3]) * mesh_.V()[i];
	for (label j=0; j<=3; j++)
	{
	    mtrx_().source()[i] -= k_[j]*oldV_[j][i]*oldPsi_[j][i];
	}
    }
    
    return mtrx_();
}

//END_OF_FILE

