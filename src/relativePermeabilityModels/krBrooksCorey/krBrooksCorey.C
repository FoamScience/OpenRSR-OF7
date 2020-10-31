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
#include "krBrooksCorey.H"
#include "oneField.H"
#include "relPermModel.H"

namespace Foam 
{
namespace twoPhaseRelPermModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class RockType>
krBrooksCorey<RockType>::krBrooksCorey
(
    const word& name,
    const dictionary& transportProperties,
    const RockType& rock
)
:
    relPermModel<RockType, 2>(name, transportProperties, rock),
    otherPhase_
    (
        this->phaseNames_[0] == this->canonicalPhases_[0]
        ? this->phaseNames_[1]
        : this->phaseNames_[0]
    ),
    alpha_
    (
        rock.mesh().template lookupObject<phase>
            (this->canonicalPhases_[0]).alpha()
    ),
    oneCellMesh_
    (
        (
            Scr0_.value() != -1 or Sor0_.value() != -1 or
            krcMax0_.value() != -1 or krcMax0_.value() != -1 or
            mc0_.value() != -1 or mo0_.value() != -1
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
    Scr0_
    (
        this->krDict_.found(this->canonicalPhases_[0]+".alphaIrr")
        ? dimensionedScalar
        (
            this->canonicalPhases_[0]+".alphaIrr", dimless, this->krDict_
        )
        : dimensionedScalar(this->canonicalPhases_[0]+".alphaIrr", dimless, -1)
    ),
    Sor0_
    (
        this->krDict_.found(otherPhase_+".alphaRes")
        ? dimensionedScalar
        (
            otherPhase_+".alphaRes", dimless, this->krDict_
        )
        : dimensionedScalar(otherPhase_+".alphaRes", dimless, -1)
    ),
    Scr_
    (
        IOobject
        (
            this->canonicalPhases_[0]+".alphaIrr",
            rock.mesh().time().timeName(),
            rock.mesh(),
            Scr0_.value() == -1
                ? IOobject::MUST_READ
                : IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Scr0_.value() == -1 ? rock.mesh() : oneCellMesh_(),
        Scr0_
    ),
    Sor_
    (
        IOobject
        (
            otherPhase_+".alphaRes",
            rock.mesh().time().timeName(),
            rock.mesh(),
            Sor0_.value() == -1
                ? IOobject::MUST_READ
                : IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sor0_.value() == -1 ? rock.mesh() : oneCellMesh_(),
        Sor0_
    ),
    mc0_
    (
        this->krDict_.found(this->canonicalPhases_[0]+".m")
        ? dimensionedScalar
        (
            this->canonicalPhases_[0]+".m", dimless, this->krDict_
        )
        : dimensionedScalar(this->canonicalPhases_[0]+".m", dimless, -1)
    ),
    mo0_
    (
        this->krDict_.found(otherPhase_+".m")
        ? dimensionedScalar
        (
        otherPhase_+".m", dimless, this->krDict_
        )
        : dimensionedScalar(otherPhase_+".m", dimless, -1)
    ),
    mc_
    (
        IOobject
        (
            this->canonicalPhases_[0]+".m",
            rock.mesh().time().timeName(),
            rock.mesh(),
            mc0_.value() == -1
                ? IOobject::MUST_READ
                : IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mc0_.value() == -1 ? rock.mesh() : oneCellMesh_(),
        mc0_
    ),
    mo_
    (
        IOobject
        (
            otherPhase_+".m",
            rock.mesh().time().timeName(),
            rock.mesh(),
            Sor0_.value() == -1
                ? IOobject::MUST_READ
                : IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mo0_.value() == -1 ? rock.mesh() : oneCellMesh_(),
        mo0_
    ),
    krcMax0_
    (
        this->krDict_.found(this->canonicalPhases_[0]+".krMax")
        ? dimensionedScalar
        (
            this->canonicalPhases_[0]+".krMax", dimless, this->krDict_
        )
        : dimensionedScalar(this->canonicalPhases_[0]+".krMax", dimless, -1)
    ),
    kroMax0_
    (
        this->krDict_.found(otherPhase_+".krMax")
        ? dimensionedScalar
        (
        otherPhase_+".krMax", dimless, this->krDict_
        )
        : dimensionedScalar(otherPhase_+".krMax", dimless, -1)
    ),
    krcMax_
    (
        IOobject
        (
            this->canonicalPhases_[0]+".krMax",
            rock.mesh().time().timeName(),
            rock.mesh(),
            kroMax0_.value() == -1
                ? IOobject::MUST_READ
                : IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        krcMax0_.value() == -1 ? rock.mesh() : oneCellMesh_(),
        krcMax0_
    ),
    kroMax_
    (
        IOobject
        (
            otherPhase_+".krMax",
            rock.mesh().time().timeName(),
            rock.mesh(),
            Sor0_.value() == -1
                ? IOobject::MUST_READ
                : IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        kroMax0_.value() == -1 ? rock.mesh() : oneCellMesh_(),
        kroMax0_
    )
{
    // TODO: Error checking
    // Alpha field bounds, coefficient legitimacy, print coeffs in debug
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class RockType>
krBrooksCorey<RockType>::~krBrooksCorey() {}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class RockType>
void krBrooksCorey<RockType>::correct()
{
    // Refs to kr fields
    auto& kr1 = this->operator[](this->krName(this->canonicalPhases_[0]));
    auto& kr2 = this->operator[](this->krName(otherPhase_));
    auto& dkr1 = this->operator[]
        (this->dkrName(this->canonicalPhases_[0], this->canonicalPhases_[0]));
    auto& dkr2 = this->operator[]
        (this->dkrName(otherPhase_, this->canonicalPhases_[0]));

    forAll(alpha_.internalField(), ci)
    {
        if( alpha_[ci] < Scr_[ci] )
        {
            kr1[ci] = VSMALL;
            kr2[ci] = kroMax_[ci];
            dkr1[ci] = 0;
            dkr2[ci] = 0;
        } else if( alpha_[ci] > (1-Sor_[ci]) ) {
            kr1[ci] = krcMax_[ci];
            kr2[ci] = VSMALL;
            dkr1[ci] = 0;
            dkr2[ci] = 0;
        } else {
            // Calculate Krs
            scalar SoeUpper = 1-alpha_[ci]-Sor_[ci];
            scalar SceUpper = alpha_[ci]-Scr_[ci];
            scalar SceLower = 1-Scr_[ci]-Sor_[ci];
            kr1[ci] = krcMax_[ci] * pow(SceUpper/SceLower, mc_[ci]);
            kr2[ci] = kroMax_[ci] * pow(SoeUpper/SceLower, mo_[ci]);

            // Calculate Kr Derivatives irt S_
            if (SoeUpper != 0 && SceUpper != 0)
            {
                dkr1[ci] = mc_[ci]*kr1[ci]/SceUpper;
                dkr2[ci] = -mo_[ci]*kr2[ci]/SoeUpper;
            } else {
                dkr1[ci] = mc_[ci] * krcMax_[ci] *
                            pow(SceUpper/SceLower, mc_[ci]-1)/SceLower;
                dkr2[ci] = -mo_[ci] * kroMax_[ci] *
                            pow(SoeUpper/SceLower, mo_[ci]-1)/SceLower;
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoPhaseRelPermModels

} // End namespace Foam

// ************************************************************************* //
