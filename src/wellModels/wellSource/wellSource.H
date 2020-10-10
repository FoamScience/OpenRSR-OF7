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
    Foam::wellSource

Description
    A template base class for well sources for blackoil-like simulations; 
    Describing the relationship q = f(p), mainly how to calculate coeffs
    in q = coeff0 p + coeff1 BHP + coeff2

    Where:
    \vartable
        q               | Phase flow rate into/out of the well
        p               | Phase pressure
        BHP             | Bottom Hole Pressure of the well
    \endtable

Note
    Kind of stateless object; planning on keeping it that way.

SourceFiles
    wellSource.C

\*---------------------------------------------------------------------------*/

#ifndef wellSource_H
#define wellSource_H

#include "phase.H"
#include "rock.H"
#include "sourceProperties.H"
#include "fvScalarMatrix.H"
#include "fvMesh.H"
#include "relPermModel.H"
#include "capPressModel.H"
#include "addToTemplatedRunTimeSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Declare rock types
using iRock = rock<Isotropic>;
using dRock = rock<DiagAnisotropic>;
using aRock = rock<Anisotropic>;

/*---------------------------------------------------------------------------*\
                         Class wellSource Declaration
\*---------------------------------------------------------------------------*/

template<class RockType, int nPhases>
class wellSource
{

protected:

    // Protected Data

        //- Source name
        word name_;

        //- Attached phase
        const phase& phase_;

        //- Ref to Well Source dictionary
        const dictionary& wellSourceDict_;

        //- Const access to the rock object
        const RockType& rock_;

        //- Canonical Pressure Field
        const volScalarField& p_;

        //- Const-access to Kr model which is directly responsible for phase kr
        const relPermModel<RockType, nPhases>& krModel_;

        //- Const-access to Pc model which is directly responsible for pc
        const capPressModel<RockType, nPhases>& pcModel_;
        
public:

    //- Runtime type information
    ClassName("wellSource");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        wellSource,
        dictionary,
        (
            const word& name,
            const phase& attachedPhase,
            const dictionary& wellSourceDict,
            const RockType& rock
        ),
        (name, attachedPhase, wellSourceDict, rock)
    );

    // Constructors

        //- Construct from components
        wellSource
        (
            const word& name,
            const phase& attachedPhase,
            const dictionary& wellSourceDict,
            const RockType& rock
        );

        //- Construct from copy
        wellSource(const wellSource&) = delete;

    //- Destructor
    virtual ~wellSource();

    // Selectors

        //- Return a reference to the selected wellSource
        static autoPtr<wellSource> New
        (
            const word& name,
            const phase& attachedPhase,
            const dictionary& wellSourceDict,
            const RockType& rock
        );

    // Member Functions

        //- Return well source name
        const word& name() const {
            return name_;
        }

        //- Return well source dict
        const word& dict() const {
            return wellSourceDict_;
        }

        //- Return phase name
        const word& phaseName() const {
            return phase_.name();
        }

        //- Return const-access to the rock
        const RockType& rock() const {
            return rock_;
        }

        //- How to calculate coeff0 in q = coeff0 p + coeff1 BHP + coeff2
        virtual void calculateCoeff0
        (
            scalarList& coeff0,
            sourceProperties& srcProps,
            const labelList& cellIDs
        ) = 0;

        //- How to calculate coeff1 in q = coeff0 p + coeff1 BHP + coeff2
        virtual void calculateCoeff1
        (
            scalarList& coeff1,
            sourceProperties& srcProps,
            const labelList& cellIDs
        ) = 0;

        //- How to calculate coeff2 in q = coeff0 p + coeff1 BHP + coeff2
        virtual void calculateCoeff2
        (
            scalarList& coeff2,
            sourceProperties& srcProps,
            const labelList& cellIDs
        ) = 0;

    // Member Operators

        //- Disallow default bitwise assignment
        wellSource& operator=(const wellSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#define maketwoPhaseTemplatedWellSource\
(baseType,childType,TargName,ns,childTArg)                                     \
                                                                               \
    /* Define child type name */                                               \
    template<>                                                                 \
    const ::Foam::word ns::childType<childTArg >::                             \
        typeName(ns::childType<childTArg >::typeName_());                      \
                                                                               \
    /* Define child debug switch */                                            \
    template<>                                                                 \
    int ns::childType<childTArg >::debug                                       \
    (                                                                          \
        ::Foam::debug::debugSwitch(ns::childType<childTArg >::typeName_(), 0)  \
    );                                                                         \
    /* Register the debug switch */                                            \
    template<>                                                                 \
    const Foam::RegisterDebugSwitch<ns::childType<childTArg >>                 \
    Foam::RegisterDebugSwitch<ns::childType<childTArg >>                       \
        ::registerDebugSwitch(ns::childType<childTArg >::typeName_());         \
                                                                               \
    /* Add child to the base's run-time selection tabe */                      \
    namespace ns {                                                             \
    baseType<childTArg, 2>::                                                   \
    adddictionaryConstructorToTable<childType<childTArg >>                     \
    add##phaseType##TargName##argNames##ConstructorTo##baseType##Targ##Table_; \
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "wellSource.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
