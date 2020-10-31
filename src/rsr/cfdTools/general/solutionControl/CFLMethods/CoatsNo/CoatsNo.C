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

#include "CoatsNo.H"
#include "fvcSurfaceIntegrate.H"
#include "surfaceFieldsFwd.H"

namespace Foam 
{
namespace twoPhaseCFLMethods
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class RockType>
CoatsNo<RockType>::CoatsNo
(
    const word& name,
    const dictionary& algorithmProperties,
    const surfaceScalarField& phi,
    const wordList& phaseNames,
    const RockType& rock
)
:
    CFLMethod<RockType, 2>(name, algorithmProperties, phi, phaseNames, rock),
    cPhase_(rock.mesh().template lookupObject<phase>(phaseNames[0])),
    nPhase_(rock.mesh().template lookupObject<phase>(phaseNames[1])),
    fieldNames_
    (
        algorithmProperties.found("fieldNames")
        ? algorithmProperties.subDict("fieldNames")
        : dictionary()
    ),
    muRatio_(name+".muRatio", cPhase_.mu()/nPhase_.mu()),
    dPhi_ (rock.mesh().nCells(), 0.0)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class RockType>
CoatsNo<RockType>::~CoatsNo() {}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class RockType>
void CoatsNo<RockType>::correct()
{
    // Get refs to phase kr and pc fields
    const volScalarField& krc =
        this->rock_.mesh().template lookupObject<volScalarField>
        (
            fieldNames_.lookupOrAddDefault<word>
            (cPhase_.name()+".kr",cPhase_.name()+".kr")
        );
    const volScalarField& krn =
        this->rock_.mesh().template lookupObject<volScalarField>
        (
            fieldNames_.lookupOrAddDefault<word>
            (nPhase_.name()+".kr",nPhase_.name()+".kr")
        );
    const volScalarField& dkrcdS =
        this->rock_.mesh().template lookupObject<volScalarField>
        (
            fieldNames_.lookupOrAddDefault<word>
            (
                cPhase_.name()+".dkrdS",
                cPhase_.name()+".dkrdS("+cPhase_.name()+")"
            )
        );
    const volScalarField& dkrndS =
        this->rock_.mesh().template lookupObject<volScalarField>
        (
            fieldNames_.lookupOrAddDefault<word>
            (
                nPhase_.name()+".dkrdS",
                nPhase_.name()+".dkrdS("+cPhase_.name()+")"
            )
        );

    // Gravitational acceleration
    const UniformDimensionedField<vector>& g =
        this->rock_.mesh().template lookupObject
        <UniformDimensionedField<vector>>
        (
            fieldNames_.lookupOrAddDefault<word>("g", "g")
        );

    // face-interpolated absolute permeability
    using SurfacePermType = GeometricField
        <typename RockType::KcmptType, fvsPatchField, surfaceMesh>;
    const SurfacePermType& Kf = 
        this->rock_.mesh().template lookupObject<SurfacePermType>
        (
            fieldNames_.lookupOrAddDefault<word>("Kf", "Kf")
        );

    // Update viscosity ratio
    muRatio_ = cPhase_.mu()/nPhase_.mu();

    // Preparations
    dimensionedScalar rSmall("epsRate",dimVolume/dimTime, vSmall);
    surfaceScalarField magPhi = mag(this->phi_);
    scalarField symmPhaseKr = 
            muRatio_*Foam::pow(krn,2) + 2 * krc * krn + 1/muRatio_*Foam::pow(krc,2);

    // Inertia's contribution to fractional flux
    dPhi_ = (dkrcdS*krn - dkrndS*krc)/symmPhaseKr;

    // Gravity's contribution to fractional flux
    dPhi_-= this->rock_.K()
            * (nPhase_.rho()-cPhase_.rho())
            * fvc::surfaceSum(mag(this->rock_.mesh().Sf() & g))
            / fvc::surfaceSum(magPhi+rSmall)
            * (Foam::pow(krn,2)*dkrcdS/nPhase_.mu() + Foam::pow(krc,2)*dkrndS/cPhase_.mu())
            / symmPhaseKr;

    // Update CFL number
    this->CFLNo_ = this->rock_.mesh().time().deltaT()/this->rock_.porosity()
            * dPhi_ * fvc::surfaceSum(magPhi);

    // Capillarity's contribution to CFL Number
    word dpcName = fieldNames_.lookupOrAddDefault<word>
            ("dpc", cPhase_.name()+".dpcdS("+cPhase_.name()+")");
    if
    (
        this->rock_.mesh().template lookupClass<volScalarField>().found(dpcName)
    )
    {
        const volScalarField& dpc = 
            this->rock_.mesh().template lookupObject<volScalarField>
            (
                fieldNames_.lookupOrAddDefault<word>
                ("dpc", cPhase_.name()+".dpcdS("+cPhase_.name()+")")
            );

        this->CFLNo_+= this->rock_.mesh().time().deltaT()/this->rock_.porosity()
                * 2 * mag(dpc) * fvc::surfaceSum
                ( Kf * this->rock_.mesh().magSf() / mag(this->rock_.mesh().delta()) )
                * (krn*krc/(cPhase_.mu()*krn+nPhase_.mu()*krc));
    }

    // Divide by cell volume
    this->CFLNo_ /= this->rock_.mesh().V();

    // Report Findings
    Info << "CFL " << typeName_()
        << " mean: " << gAverage(this->CFLNo_)
        << " max: " << gMax(this->CFLNo_) << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoPhaseCFLMethods

} // End namespace Foam

// ************************************************************************* //
