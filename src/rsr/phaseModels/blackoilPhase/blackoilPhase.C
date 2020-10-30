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

#include "blackoilPhase.H"
#include "fvcReconstruct.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace phases {
        defineTypeNameAndDebug (blackoilPhase, 0);
        addToRunTimeSelectionTable
        (
            phase, blackoilPhase, dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phases::blackoilPhase::blackoilPhase
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const mixtureType& mT
)
:
    phase(name, mesh, transportProperties, mT)
{
}


Foam::phases::blackoilPhase::blackoilPhase
(
    const blackoilPhase& ph
)
:
    phase(ph)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phases::blackoilPhase::~blackoilPhase()
{}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

void Foam::phases::blackoilPhase::correct()
{
    if (BModel_->isIncompressible()) return;
    // Correct FVFModel
    BModel_->correct();

    // Update current density
    forAll(rho_.internalField(), ci)
    {
        rho_[ci] = rhoSc_.value()*BModel_->rFVF()[ci];
    }

    // Reconstruct velocity to reflect flux
    U_ = fvc::reconstruct(phi());
    U_.correctBoundaryConditions();
}

// ************************************************************************* //
