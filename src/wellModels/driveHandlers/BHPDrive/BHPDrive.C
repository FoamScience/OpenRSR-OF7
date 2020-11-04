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

#include "BHPDrive.H"

namespace Foam 
{
namespace driveHandlers
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
BHPDrive<RockType, nPhases>::BHPDrive
(
    const word& name,
    const dictionary& driveDict,
    HashTable<autoPtr<wellSource<RockType, nPhases>>>& sources,
    sourceProperties& srcProps,
    HashPtrTable<fvScalarMatrix>& matrices
)
:
    driveHandler<RockType,nPhases>(name,driveDict,sources,srcProps,matrices),
    phases_(sources.toc())
{
    if (sources.empty())
    {
        FatalErrorInFunction
            << "No Well Source describers were passed to drive handler"
            << name << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
BHPDrive<RockType, nPhases>::~BHPDrive() {}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class RockType, int nPhases>
void BHPDrive<RockType, nPhases>::correct()
{
    // Get refs to mesh, time, pressure and phase matrix
    const fvMesh& mesh = this->wellSources_[phases_[0]]->rock().mesh();
    const scalar& timeValue = mesh.time().timeOutputValue();
    const label& opSign = this->srcProps_.operationSign();

    // Get interpolated value for imposed BHP
    scalar BHP = this->driveSeries_->interpolate(timeValue)[0];

    // Calculate sources for each phase
    forAll(phases_, pi)
    {
        // Do nothing if this is an injector working on the wrong phase
        if
        (
            this->srcProps_.operation() 
            == sourceProperties::operationHandling::injection
            and
            this->srcProps_.injPhase() != phases_[pi]
        )
        {
            continue;
        }

        fvScalarMatrix& phEqn = *(this->matrices_[phases_[pi]]);

        // Get well equation coefficients from well source describer
        this->wellSources_[phases_[pi]]->calculateCoeff0
        (
            this->coeffs_[0], this->srcProps_, this->cells_
        );
        this->wellSources_[phases_[pi]]->calculateCoeff1
        (
            this->coeffs_[1], this->srcProps_, this->cells_
        );
        this->wellSources_[phases_[pi]]->calculateCoeff2
        (
            this->coeffs_[2], this->srcProps_, this->cells_
        );

        // Loop through cells and add diagonal coeffs and matrix source
        forAll(this->cells_, ci)
        {
            const label cellID = this->cells_[ci];
            phEqn.diag()[cellID] += opSign * this->coeffs_[0][ci];
            phEqn.source()[cellID] += opSign
                * (this->coeffs_[1][ci]*BHP + this->coeffs_[2][ci]);
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace driveHandlers

} // End namespace Foam

// ************************************************************************* //
