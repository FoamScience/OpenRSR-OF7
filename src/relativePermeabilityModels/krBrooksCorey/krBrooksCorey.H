/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::twoPhaseRelPermModels::krBrooksCorey

Description
    Two-Phase relative permeability model implementing Modified Brooks-Corey
    approach. 
    \f[
        k_{rw} = k_{rw, max}S_{w, norm}^n
    \f]
    \f[
        k_{rn} = k_{rn, max}(1-S_{w, norm})^n
    \f]
    Where 
    \f[
        S_{w, norm} = \frac{S_w, S_{w, irr}}{1-S_{w, irr}-S_{n, irr}}
    \f]
    \f[
        S_{w, norm} = 0\ if\ S_w<S_{w, irr}
    \f]
    \f[
        S_{w, norm} = 1\ if\ S_w>1-S_{n, irr}
    \f]
    \vartable
        k_r             | Relative permeability
        w,n             | Wetting and NonWetting indices respectively
        S_{w, norm}     | Normalized wetting-phase saturation
        S_{w, irr}      | Irreducible wetting phase saturation
        S_{n, irr}      | Irreducible wetting phase saturation
    \endtable

    Use:
    - Phase exponent (m) = 4: Fall back to 1954 Corey expressions for oil/gas 
Note
    Mostly suitable for water-oil systems.
    And the canonical phase is usually the wetting phase (water).

SourceFiles
    krBrooksCorey.C

\*---------------------------------------------------------------------------*/

#ifndef krBrooksCorey_H
#define krBrooksCorey_H

#include "VolatileDimensionedField.H"
#include "dimensionedScalarFwd.H"
#include "relPermModel.H"
#include "singleCellFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace twoPhaseRelPermModels
{

/*---------------------------------------------------------------------------*\
                        Class krBrooksCorey Declaration
\*---------------------------------------------------------------------------*/

template<class RockType>
class krBrooksCorey
:
    public relPermModel<RockType, 2>
{
    //- Name of the other phase
    word otherPhase_;

    //- Const-access to canonical saturation
    const volScalarField& alpha_;

    // A single-cell mesh
    autoPtr<singleCellFvMesh> oneCellMesh_;

    //- Residual saturation for canonical and the other phase
    dimensionedScalar Scr0_, Sor0_;
    VolatileDimensionedField<scalar, volMesh> Scr_, Sor_;

    //- Model exponent for canonical and the other phase
    dimensionedScalar mc0_, mo0_;
    VolatileDimensionedField<scalar, volMesh> mc_, mo_;

    //- End-point relative permeabilities
    dimensionedScalar krcMax0_, kroMax0_;
    VolatileDimensionedField<scalar, volMesh> krcMax_, kroMax_;

public:

    //- Runtime type information
    ClassName("BrooksCorey");

    // Constructors

        // Construct from components
        krBrooksCorey
        (
            const word& name,
            const dictionary& transportProperties,
            const RockType& rock
        );

        // Construct from copy
        krBrooksCorey(const krBrooksCorey&) = delete;

	// Destructor
	virtual ~krBrooksCorey();

    // Member Functions

		//- Correct relative-permeability fields
		virtual void correct();

    // Member operators

        //- Disallow default bitwise assignment
        krBrooksCorey& operator=(const krBrooksCorey&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoPhaseRelPermModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "krBrooksCorey.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
