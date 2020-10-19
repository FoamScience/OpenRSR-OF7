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
    wellSource<RockType, nPhases>& source,
    sourceProperties& srcProps,
    HashTable<fvScalarMatrix>& matrices
)
:
    driveHandler<RockType,nPhases>(name, driveDict, source, srcProps, matrices),
    phase_ (driveDict.lookup("phase"))
{
    if(phase_ != source.phaseName())
    {
        FatalErrorInFunction
            << "Drive handler " << name << " works on "
            << phase_ << " phase but got well source describer for phase"
            << source.phaseName() << "." << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
flowRateDrive<RockType, nPhases>::~flowRateDrive() {}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class RockType, int nPhases>
void flowRateDrive<RockType, nPhases>::correct()
{
    // Get refs to mesh, time, pressure and phase matrix
    const fvMesh& mesh = this->wellSource_.rock().mesh();
    const scalar& timeValue = mesh.time().timeOutputValue();
    const volScalarField& p = mesh.lookupObject<volScalarField>("p");
    fvScalarMatrix& phEqn = this->matrices_[phase_];

    // Get well equation coefficients from well source describer
    this->wellSource_.calculateCoeff0
    (
        this->coeffs_[0], this->srcProps_, this->cells_
    );
    this->wellSource_.calculateCoeff1
    (
        this->coeffs_[1], this->srcProps_, this->cells_
    );
    this->wellSource_.calculateCoeff2
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

    // Get interpolated value for imposed phase flowrate
    scalar qt = this->driveSeries_->interpolate(timeValue)[0];

    // Loop through cells and add diagonal coeffs and matrix source
    forAll(this->cells_, ci)
    {
        const label cellID = this->cells_[ci];
        phEqn.diag()[cellID] += a[ci]*(1-b[ci]);
        phEqn.source()[cellID] += c[ci] + b[ci]*qt - b[ci]*cSum;
        forAll(this->cells_, cj)
        {
            const label cellj = this->cells_[cj];
            bool notNeiCell = findIndex(mesh.cellCells()[cellID], cellj) == -1
                           ? true : false;
            if (cellj != cellID and notNeiCell)
            {
                phEqn.source()[cellID] += -b[ci]*a[cj]*p[this->cells_[cj]];
            }
        }
    }

    // Loop through all internal faces and add off-diagonal coefficients
    forAll(iFaces, fi)
    {
        const label faceID = iFaces[fi];
        label in = findIndex(this->cells_, mesh.lduAddr().lowerAddr()[faceID]);
        label jn = findIndex(this->cells_, mesh.lduAddr().upperAddr()[faceID]);
        phEqn.lower()[faceID] += -b[in] * a[jn];
        if (phEqn.hasUpper()) phEqn.upper()[faceID] += -b[jn] * a[in];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace driveHandlers

} // End namespace Foam

// ************************************************************************* //
