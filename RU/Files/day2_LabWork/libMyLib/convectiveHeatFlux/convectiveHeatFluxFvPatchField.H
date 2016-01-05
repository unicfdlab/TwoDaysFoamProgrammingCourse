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

Class
    Foam::convectiveHeatFluxFvPatchField

Description
    Robin B.C. for solving equation lambda*dT/dn + alpha(T-Tref)=0
    
    \verbatim
        wall1
        {
            type            convectiveHeatFlux;
            lambda	    1.0e-3;
            alpha	    1000;
            Tref	    293.0;
        }
    \endverbatim

SourceFiles
    convectiveHeatFluxFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef convectiveHeatFluxFvPatchField_H
#define convectiveHeatFluxFvPatchField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class advectiveFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class convectiveHeatFluxFvPatchField
:
    public mixedFvPatchField<Type>
{
protected:

    // Private data

        //- Temperature diffusivity
        scalar lambda_;

        //- Heat exchange coefficient
	scalar alpha_;

        //- Reference temperature
        scalar Tref_;

public:

    //- Runtime type information
    TypeName("convectiveHeatFlux");


    // Constructors

        //- Construct from patch and internal field
        convectiveHeatFluxFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        convectiveHeatFluxFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given convectiveHeatFluxFvPatchField
        //  onto a new patch
        convectiveHeatFluxFvPatchField
        (
            const convectiveHeatFluxFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        convectiveHeatFluxFvPatchField
        (
            const convectiveHeatFluxFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new convectiveHeatFluxFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        convectiveHeatFluxFvPatchField
        (
            const convectiveHeatFluxFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new convectiveHeatFluxFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "convectiveHeatFluxFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
