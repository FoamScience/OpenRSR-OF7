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

#include "DimensionedField.H"
#include "error.H"
#include "pcBrooksCorey.H"

namespace Foam 
{
namespace twoPhaseCapPressModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class RockType>
pcBrooksCorey<RockType>::pcBrooksCorey
(
    const word& name,
    const dictionary& transportProperties,
    const RockType& rock
)
:
    capPressModel<RockType, 2>(name, transportProperties, rock),
    alpha_
    (
        rock.mesh().template lookupObject<phase>
            (this->canonicalPhases_[0]).alpha()
    ),
    oneCellMesh_
    (
        (
            pcSmin0_.value() != -1 or pcSmax0_.value() != -1 
            or pc00_.value() != -1 
        )
        ?  new singleCellFvMesh
           (
               IOobject
               (
                   name+".oneCellMesh",
                   rock.mesh().polyMesh::instance(),
                   rock.mesh().time(),
                   IOobject::NO_READ,
                   IOobject::NO_WRITE
               ),
               rock.mesh()
           )
        :  nullptr
    ),
    pcSmin0_
    (
        this->pcDict_.found(alpha_.name()+".PcMin")
        ? dimensionedScalar
        (
            alpha_.name()+".PcMin", dimless, this->pcDict_
        )
        : dimensionedScalar(alpha_.name()+".PcMin", dimless, -1)
    ),
    pcSmax0_
    (
        this->pcDict_.found(alpha_.name()+".PcMax")
        ? dimensionedScalar
        (
            alpha_.name()+".PcMax", dimless, this->pcDict_
        )
        : dimensionedScalar(alpha_.name()+".PcMax", dimless, -1)
    ),
    pc00_
    (
        this->pcDict_.found("pc0")
        ? dimensionedScalar
        (
            "pc0", dimPressure, this->pcDict_
        )
        : dimensionedScalar("pc0", dimPressure, -1)
    ),
    pcSmin_
    (
        IOobject
        (
            alpha_.name()+".PcMin",
            rock.mesh().time().timeName(),
            rock.mesh(),
            pcSmin0_.value() == -1
                ? IOobject::MUST_READ
                : IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        pcSmin0_.value() == -1 ? rock.mesh() : oneCellMesh_(),
        pcSmin0_
    ),
    pcSmax_
    (
        IOobject
        (
            alpha_.name()+".PcMax",
            rock.mesh().time().timeName(),
            rock.mesh(),
            pcSmax0_.value() == -1
                ? IOobject::MUST_READ
                : IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        pcSmax0_.value() == -1 ? rock.mesh() : oneCellMesh_(),
        pcSmax0_
    ),
    pc0_
    (
        IOobject
        (
            name + ".pc0",
            rock.mesh().time().timeName(),
            rock.mesh(),
            pc00_.value() == -1
                ? IOobject::MUST_READ
                : IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        pc00_.value() == -1 ? rock.mesh() : oneCellMesh_(),
        pc00_
    ),
    n0_
    (
        this->pcDict_.found("n")
        ? dimensionedScalar
        (
            "n", dimless, this->pcDict_
        )
        : dimensionedScalar("n", dimless, -1)
    ),
    n_
    (
        IOobject
        (
            name + ".n",
            rock.mesh().time().timeName(),
            rock.mesh(),
            n0_.value() == -1
                ? IOobject::MUST_READ
                : IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        n0_.value() == -1 ? rock.mesh() : oneCellMesh_(),
        n0_
    )
{
    // TODO: Error checking
    // Alpha field bounds, coefficient legitimacy, print coeffs in debug
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class RockType>
pcBrooksCorey<RockType>::~pcBrooksCorey() {}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class RockType>
void pcBrooksCorey<RockType>::correct()
{
    // Refs to kr fields
    auto& pc = this->pcTable_[this->pcName(this->canonicalPhases_[0])];
    auto& dpc = this->pcTable_
        [this->dpcName(this->canonicalPhases_[0], this->canonicalPhases_[0])];

    forAll(alpha_.internalField(), ci)
    {
        if (alpha_[ci] <= pcSmin_[ci])
        {
            FatalErrorInFunction
                << "Capillary pressure for phase "
                << this->canonicalPhases_[0] << " is not defined where "
                << "the saturation is equal or less than " << pcSmin_[ci]
                << " at cell " << ci
                << exit(FatalError);
        } else {
            scalar SLower =  pcSmax_[ci] - pcSmin_[ci];
            scalar Snorm = (alpha_[ci] - pcSmin_[ci])/SLower;
            pc[ci] = pc0_[ci] * pow(Snorm, -n_[ci]);
            dpc[ci] = -n_[ci] * pc0_[ci] * pow(Snorm, -n_[ci]-1) / SLower;
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoPhaseCapPressModels

} // End namespace Foam

// ************************************************************************* //
