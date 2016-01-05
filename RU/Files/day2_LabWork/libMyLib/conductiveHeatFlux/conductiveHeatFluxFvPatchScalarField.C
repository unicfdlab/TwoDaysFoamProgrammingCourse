/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "conductiveHeatFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

conductiveHeatFluxFvPatchScalarField::
conductiveHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    q_(p.size(), 0.0),
    lambdaWall_(0.0)
{
}


conductiveHeatFluxFvPatchScalarField::
conductiveHeatFluxFvPatchScalarField
(
    const conductiveHeatFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    q_(ptf.q_, mapper),
    lambdaWall_(ptf.lambdaWall_)
{}


conductiveHeatFluxFvPatchScalarField::
conductiveHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    q_("q", dict, p.size()),
    lambdaWall_(dict.lookupOrDefault<scalar>("lambdaWall", 1e-6))
{
    if (dict.found("gradient"))
    {
	gradient() = Field<scalar> ("gradient", dict, p.size());
    }
    else
    {
	gradient() = 0;
    }
    
    if (dict.found("value"))
    {
	fvPatchField<scalar>::operator= (Field<scalar>("value", dict, p.size()));
    }
    else
    {
	fvPatchField<scalar>::operator=(patchInternalField());
    }
}


conductiveHeatFluxFvPatchScalarField::
conductiveHeatFluxFvPatchScalarField
(
    const conductiveHeatFluxFvPatchScalarField& thftpsf
)
:
    fixedGradientFvPatchScalarField(thftpsf),
    q_(thftpsf.q_),
    lambdaWall_(thftpsf.lambdaWall_)
{}


conductiveHeatFluxFvPatchScalarField::
conductiveHeatFluxFvPatchScalarField
(
    const conductiveHeatFluxFvPatchScalarField& thftpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(thftpsf, iF),
    q_(thftpsf.q_),
    lambdaWall_(thftpsf.lambdaWall_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void conductiveHeatFluxFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
    q_.autoMap(m);
}


void conductiveHeatFluxFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const conductiveHeatFluxFvPatchScalarField& thftptf =
        refCast<const conductiveHeatFluxFvPatchScalarField>
        (
            ptf
        );

    q_.rmap(thftptf.q_, addr);
}


void conductiveHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    gradient() = q_/ lambdaWall_;

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void conductiveHeatFluxFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    q_.writeEntry("q", os);
    os.writeKeyword("lambdaWall") << lambdaWall_ << token::END_STATEMENT << nl;
    gradient().writeEntry("gradient", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    conductiveHeatFluxFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //

