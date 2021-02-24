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

#include "flowRateDrive.H"

namespace Foam 
{
namespace driveHandlers
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
flowRateDrive<RockType, nPhases>::flowRateDrive
(
    const word& name,
    const dictionary& driveDict,
    HashTable<autoPtr<wellSource<RockType, nPhases>>>& sources,
    sourceProperties& srcProps,
    HashPtrTable<fvScalarMatrix>& matrices
)
:
    driveHandler<RockType,nPhases>(name,driveDict,sources,srcProps,matrices),
    phase_ (driveDict.lookupOrDefault<word>("phase", "water"))
{
    if(phase_ != sources[phase_]->phaseName())
    {
        FatalErrorInFunction
            << "Drive handler " << name << " works on "
            << phase_ << " phase but got well source describer for phase"
            << sources[phase_]->phaseName() << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
flowRateDrive<RockType, nPhases>::~flowRateDrive() {}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class RockType, int nPhases>
void flowRateDrive<RockType, nPhases>::correct()
{
    // Do nothing if this is an injector working on the wrong phase
    if
    (
        this->srcProps_.operation() 
        == sourceProperties::operationHandling::injection
        and
        this->srcProps_.injPhase() != phase_
    )
    {
        return;
    }

    // Get refs to mesh, time, pressure and phase matrix
    const fvMesh& mesh = this->wellSources_[phase_]->rock().mesh();
    const scalar& timeValue = mesh.time().timeOutputValue();
    const volScalarField& p = mesh.lookupObject<volScalarField>("p");
    fvScalarMatrix& phEqn = *(this->matrices_[phase_]);
    const label& opSign = this->srcProps_.operationSign();

    // Get interpolated value for imposed phase flowrate
    scalar qt = opSign * this->driveSeries_->interpolate(timeValue)[0];

    // If well has one cell
    if (this->cells_.size() == 1)
    {
        phEqn.source()[this->cells_[0]] += qt;
        return;
    }

    // Get well equation coefficients from well source describer
    this->wellSources_[phase_]->calculateCoeff0
    (
        this->coeffs_[0], this->srcProps_, this->cells_
    );
    this->wellSources_[phase_]->calculateCoeff1
    (
        this->coeffs_[1], this->srcProps_, this->cells_
    );
    this->wellSources_[phase_]->calculateCoeff2
    (
        this->coeffs_[2], this->srcProps_, this->cells_
    );
    // Get well internal faces from sourceProperties
    const auto& iFaces = this->srcProps_.faces();

    // qi = ai * pi + bi * BHP + ci
    // BHP substituted from total flowrate eqn

    // Get sum of bi for the whole well on all processes:
    scalar bSum = sum(this->coeffs_[1]);
    reduce(bSum, sumOp<scalar>());
    Pstream::scatter(bSum);

    // Divide coeff1 by its global sum 
    auto& a = this->coeffs_[0];
    auto& c = this->coeffs_[2];
    auto b = this->coeffs_[1];
    forAll(b, i) { b[i] /= bSum; }

    // Get sum of ci for the whole well on all processes:
    scalar cSum = sum(c);
    reduce(cSum, sumOp<scalar>());
    Pstream::scatter(cSum);

    // Loop through cells and add diagonal coeffs and matrix source
    forAll(this->cells_, ci)
    {
        const label cellID = this->cells_[ci];
        phEqn.diag()[cellID] += a[ci]*(1-b[ci]);
        phEqn.source()[cellID] += c[ci] + b[ci]*(qt-cSum);
        forAll(this->cells_, cj)
        {
            const label cellj = this->cells_[cj];
            bool notNeiCell = findIndex(mesh.cellCells()[cellID], cellj) == -1
                           ? true : false;
            if (cellj != cellID and notNeiCell)
            {
                phEqn.source()[cellID] += 
                    -b[ci]*a[cj]*p[this->cells_[cj]];
            }
        }
    }

    // Loop through all internal faces and add off-diagonal coefficients
    forAll(iFaces, fi)
    {
        const label faceID = iFaces[fi];
        label in = findIndex(this->cells_, mesh.lduAddr().lowerAddr()[faceID]);
        label jn = findIndex(this->cells_, mesh.lduAddr().upperAddr()[faceID]);
        phEqn.upper()[faceID] += -b[jn]*a[in];
        if (phEqn.asymmetric()) phEqn.lower()[faceID] += -b[in]*a[jn];
    }
    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace driveHandlers

} // End namespace Foam

// ************************************************************************* //
