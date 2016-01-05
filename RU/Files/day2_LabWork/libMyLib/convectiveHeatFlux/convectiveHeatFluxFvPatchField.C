/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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

#include "convectiveHeatFluxFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "EulerDdtScheme.H"
#include "CrankNicholsonDdtScheme.H"
#include "backwardDdtScheme.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::convectiveHeatFluxFvPatchField<Type>::convectiveHeatFluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    lambda_(1e-3),
    alpha_(1000),
    Tref_(273.0)
{
    this->refValue() = pTraits<Type>::zero;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::convectiveHeatFluxFvPatchField<Type>::convectiveHeatFluxFvPatchField
(
    const convectiveHeatFluxFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    lambda_(ptf.lambda_),
    alpha_(ptf.alpha_),
    Tref_(ptf.Tref_)
{}


template<class Type>
Foam::convectiveHeatFluxFvPatchField<Type>::convectiveHeatFluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    lambda_(dict.lookupOrDefault<scalar>("lambda", 1e-3)),
    alpha_(dict.lookupOrDefault<scalar>("alpha", 1000)),
    Tref_(dict.lookupOrDefault<scalar>("Tref", 273))
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<Type>::operator=(this->patchInternalField());
    }

    this->refValue() = *this;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::convectiveHeatFluxFvPatchField<Type>::convectiveHeatFluxFvPatchField
(
    const convectiveHeatFluxFvPatchField& ptpsf
)
:
    mixedFvPatchField<Type>(ptpsf),
    lambda_(ptpsf.lambda_),
    alpha_(ptpsf.alpha_),
    Tref_(ptpsf.Tref_)
{}


template<class Type>
Foam::convectiveHeatFluxFvPatchField<Type>::convectiveHeatFluxFvPatchField
(
    const convectiveHeatFluxFvPatchField& ptpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptpsf, iF),
    lambda_(ptpsf.lambda_),
    alpha_(ptpsf.alpha_),
    Tref_(ptpsf.Tref_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::convectiveHeatFluxFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    scalarField C = lambda_ / alpha_ * this->patch().deltaCoeffs();
    
    this->valueFraction() = 1 + 1./C;
    
    this->refValue() = Type(pTraits<Type>::one)*Tref_;
    
    this->refGrad() = Type(pTraits<Type>::zero);

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::convectiveHeatFluxFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    os.writeKeyword("lambda") << lambda_ << token::END_STATEMENT << nl;
    
    os.writeKeyword("alpha") << alpha_ << token::END_STATEMENT << nl;
    
    os.writeKeyword("Tref") << Tref_ << token::END_STATEMENT << nl;

    this->writeEntry("value", os);
}


// ************************************************************************* //
