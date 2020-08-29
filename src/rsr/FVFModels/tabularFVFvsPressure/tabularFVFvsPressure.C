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

\*---------------------------------------------------------------------------*/

#include "tabularFVFvsPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace FVFModels {
        defineTypeNameAndDebug (tabularFVFvsPressure, 0);
        addToRunTimeSelectionTable
        (
            compressibleFVFModel, tabularFVFvsPressure, dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FVFModels::tabularFVFvsPressure::tabularFVFvsPressure
(
    const word& name,
    const dictionary& phaseDict,
    const fvMesh& mesh
)
    :
    FVFModel(name, phaseDict, mesh),
    pName_
    (
        phaseDict.subDict("FVFData").lookupOrDefault<word>("pName", "p")
    ),
    p_(mesh.lookupObject<volScalarField>(pName_)),
    rFVFseries_
    (
        basicInterpolationTable<scalar>::New
        (
            phaseDict.subDict("FVFData")
        )
    )
{
}


Foam::FVFModels::tabularFVFvsPressure::tabularFVFvsPressure
(
    const tabularFVFvsPressure& fvfModel
)
    :
    FVFModel(fvfModel),
    pName_(fvfModel.pName_),
    p_(fvfModel.p_),
    rFVFseries_
    (
        basicInterpolationTable<scalar>::New
        (
            fvfModel.phaseDict_.subDict("FVFData")
        )
    )
{
}


// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

void Foam::FVFModels::tabularFVFvsPressure::correct()
{
    // Read In FVF data
    forAll(p_.internalField(), celli)
    {
        scalar pInCell = p_.internalField()[celli];
        auto interpolatedValues = rFVFseries_->interpolate(pInCell);
        rFVF_[celli] = interpolatedValues[0];
        drFVFdP_[celli] = interpolatedValues[1];
    }

    // Correct Boundary Conditions
    rFVF_.correctBoundaryConditions();
    drFVFdP_.correctBoundaryConditions();
}

// ************************************************************************* //
