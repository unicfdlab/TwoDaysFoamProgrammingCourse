#ifndef fourthOrderDdt_H
#define fourthOrderDdt_H

#include "fvMesh.H"
#include "fvMatrices.H"
#include "volFieldsFwd.H"

namespace Foam
{

class fourthOrderDdt
{

private:

    //-
    const volScalarField& psi_;
    
    //-
    const fvMesh& mesh_;
    
    //-
    scalar cTime_;
    
    //-
    List<scalar> dt_;
    
    //-
    List<scalarField> oldPsi_;
    
    //-
    List<scalarField> oldV_;
    
    //-
    List<scalar> delta_;
    
    //-
    List<scalar> k_;
    
    //-
    tmp<fvScalarMatrix> mtrx_;
    
private:

    //-
    fourthOrderDdt();
    
    //-
    fourthOrderDdt (const fourthOrderDdt&);
    
protected:

    //-
    void makeCoeffs();

    //-
    void advanceInTime();

public:

    //-
    fourthOrderDdt (const volScalarField& psi);

    //-
    const fvScalarMatrix& ddt();
};

}

#endif

//END_OF_FILE

