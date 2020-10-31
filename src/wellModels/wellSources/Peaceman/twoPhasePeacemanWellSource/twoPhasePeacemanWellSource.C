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

#include "error.H"
#include "mathematicalConstants.H"
#include "twoPhasePeacemanWellSource.H"
#include "sourceProperties.H"

namespace Foam 
{
namespace wellSources
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class RockType>
twoPhasePeacemanWellSource<RockType>::twoPhasePeacemanWellSource
(
    const word& name,
    const phase& attachedPhase,
    const dictionary& wellSourceDict,
    const RockType& rock
)
:
    wellSource<RockType, 2>(name, attachedPhase, wellSourceDict, rock),
    peacemanWellSourceCore<RockType>(name, attachedPhase, wellSourceDict, rock),
    krModel_
    (
        relPermModel<RockType, 2>::getKrModel
        (
            attachedPhase.name(), rock.mesh()
        )
    ),
    pcModel_
    (
        capPressModel<RockType, 2>::getPcModel
        (
            attachedPhase.name(), rock.mesh()
        )
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class RockType>
twoPhasePeacemanWellSource<RockType>::~twoPhasePeacemanWellSource() {}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class RockType>
void twoPhasePeacemanWellSource<RockType>::calculateCoeff0
(
    scalarList &coeff0,
    sourceProperties& srcProps,
    const labelList& cellIDs
)
{
    word krName = relPermModel<RockType,2>::krName(this->phase_.name());
    const auto& kr = krModel_[krName];
    const auto& mu = this->phase_.mu();

    this->calculateWellIndex(cellIDs, srcProps);
    coeff0.resize(srcProps.wellIndex().size());
    forAll(coeff0, ci)
    {
        const label cellID = cellIDs[ci];
        coeff0[ci] = - srcProps.wellIndex()[ci] * kr[cellID] / mu[cellID];
    }
}

template<class RockType>
void twoPhasePeacemanWellSource<RockType>::calculateCoeff1
(
    scalarList &coeff1,
    sourceProperties& srcProps,
    const labelList& cellIDs
)
{
    word krName = relPermModel<RockType,2>::krName(this->phase_.name());
    const auto& kr = krModel_[krName];
    const auto& mu = this->phase_.mu();

    this->calculateWellIndex(cellIDs, srcProps);
    coeff1.resize(srcProps.wellIndex().size());
    forAll(coeff1, ci)
    {
        const label cellID = cellIDs[ci];
        coeff1[ci] = srcProps.wellIndex()[ci] * kr[cellID] / mu[cellID];
    }
}

template<class RockType>
void twoPhasePeacemanWellSource<RockType>::calculateCoeff2
(
    scalarList &coeff2,
    sourceProperties& srcProps,
    const labelList& cellIDs
)
{
    const auto& rho = this->phase_.rho();
    const auto& g = srcProps.g();
    scalar gg = (g && g).value()/(mag(g).value() + vSmall);
    scalar ZBH =
        (gg == 0) ? 0 : ((srcProps.gLowerCell().first() && g)/gg).value();

    this->calculateWellIndex(cellIDs, srcProps);
    coeff2.resize(srcProps.wellIndex().size());

    forAll(coeff2, ci)
    {
        const label cellID = cellIDs[ci];
        scalar cellZ = 
            (gg == 0) ? 0 : ((rho.mesh().C()[cellID] && g)/gg).value();
        // TODO: Consider adding capillary pressure support
        coeff2[ci] = - srcProps.wellIndex()[ci] * rho[cellID] * gg
            * (ZBH - cellZ);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wellSources

} // End namespace Foam

// ************************************************************************* //
