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
    Foam::twoPhaseCapPressModels::pcBrooksCorey

Description
    Two-Phase capillary pressure model implementing Standard Brooks-Corey
    approach. 
    \f[
        P_c = P_{c0} S_{pc, norm}^{-\alpha}
    \f]

    The normalized saturation is given by:
    \f[
        S_{w, norm} = \frac{S_w - S_{pc, min}}{S_{pc, max}-S_{pc, min}}
    \f]
    Where \f$S_{pc, min} < S_{w, irr}\f$
    See Foam.twoPhaseCapPressModels.krBrooksCoery relative permeability model.
    \vartable
        P_c                   | Capillary pressure
        S_w                   | Canonical phase saturation
        S_{pc, norm}          | Normalized capillary saturation
        S_{pc, min}           | Min. capillary saturation
        S_{pc, max}           | Max. cappilary saturation
        \frac{1}{n}           | Pore-size distribution index        
    \endtable
Note
    The derivate \f$\frac{dP_c}{dS_w}\f$ is calculated theoretically.

SourceFiles
    pcBrooksCorey.C

\*---------------------------------------------------------------------------*/

#ifndef pcBrooksCorey_H
#define pcBrooksCorey_H

#include "VolatileDimensionedField.H"
#include "dimensionedScalarFwd.H"
#include "capPressModel.H"
#include "singleCellFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace twoPhaseCapPressModels
{

/*---------------------------------------------------------------------------*\
                        Class pcBrooksCorey Declaration
\*---------------------------------------------------------------------------*/

template<class RockType>
class pcBrooksCorey
:
    public capPressModel<RockType, 2>
{
    //- Const-access to canonical saturation
    const volScalarField& alpha_;

    // A single-cell mesh
    autoPtr<singleCellFvMesh> oneCellMesh_;

    //- Model coeffs
    dimensionedScalar pcSmin0_, pcSmax0_, pc00_;
    VolatileDimensionedField<scalar, volMesh> pcSmin_, pcSmax_, pc0_;

    //- Model exponent (1/n is the pore size distribution)
    dimensionedScalar n0_;
    VolatileDimensionedField<scalar, volMesh> n_;

public:

    //- Runtime type information
    ClassName("BrooksCorey");

    // Constructors

        // Construct from components
        pcBrooksCorey
        (
            const word& name,
            const dictionary& transportProperties,
            const RockType& rock
        );

        // Construct from copy
        pcBrooksCorey(const pcBrooksCorey&) = delete;

	// Destructor
	virtual ~pcBrooksCorey();

    // Member Functions

		//- Correct relative-permeability fields
		virtual void correct();

    // Member operators

        //- Disallow default bitwise assignment
        pcBrooksCorey& operator=(const pcBrooksCorey&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoPhaseCapPressModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pcBrooksCorey.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
