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

\*---------------------------------------------------------------------------*/

#include "fourierFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fourierFvPatchScalarField::fourierFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    N_(1),
    A_(List<scalar> (1, 1.0)),
    theta_(List<scalar> (1, 0.0)),
    freq_(1)
{}


Foam::fourierFvPatchScalarField::fourierFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    N_(dict.lookupOrDefault<label>("N", 1)),
    A_(dict.lookupOrDefault<List<scalar> >("A", List<scalar>(1,1.0))),
    theta_(dict.lookupOrDefault<List<scalar> >("theta",  List<scalar>(1,0.0))),
    freq_(dict.lookupOrDefault<scalar>("freq",1.0 ))
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(0.0);
    }
}


Foam::fourierFvPatchScalarField::fourierFvPatchScalarField
(
    const fourierFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    N_(ptf.N_),
    A_(ptf.A_),
    theta_(ptf.theta_),
    freq_(ptf.freq_)
{}


Foam::fourierFvPatchScalarField::fourierFvPatchScalarField
(
    const fourierFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    N_(tppsf.N_),
    A_(tppsf.A_),
    theta_(tppsf.theta_),
    freq_(tppsf.freq_)
{}


Foam::fourierFvPatchScalarField::fourierFvPatchScalarField
(
    const fourierFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    N_(tppsf.N_),
    A_(tppsf.A_),
    theta_(tppsf.theta_),
    freq_(tppsf.freq_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fourierFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
}


void Foam::fourierFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);
}


void Foam::fourierFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalar val = 0.0;
    scalar t = this->patch().boundaryMesh().mesh().time().value();
    forAll (A_, k)
    {
	val += A_[k] * cos (2*M_PI*scalar(k+1)*freq_*t + theta_[k]);
    }
    
    operator == (val);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::fourierFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("N") << N_ << token::END_STATEMENT << nl;
    os.writeKeyword("A") << A_ << token::END_STATEMENT << nl;
    os.writeKeyword("theta") << theta_ << token::END_STATEMENT << nl;
    os.writeKeyword("freq") << freq_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fourierFvPatchScalarField
    );
}

// ************************************************************************* //
