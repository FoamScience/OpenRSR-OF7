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
    Foam::wellModel

Description
    A template base class for well models.

Resposibilities
    - Organise managed wells into groups
    - Approximate Well parameters on a specific mesh. Eg, calculate equivalent 
      radius and well index
    - Return explicit and implicit parts of well sources.

SourceFiles
    wellModel.C

\*---------------------------------------------------------------------------*/

#ifndef wellModel_H
#define wellModel_H

#include "fvMatricesFwd.H"
#include "objectRegistry.H"
#include "phase.H"
#include "rock.H"
#include "addToTemplatedRunTimeSelection.H"
#include "volFieldsFwd.H"
#include "well.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class wellModel Declaration
\*---------------------------------------------------------------------------*/

template<class RockType, int nPhases>
class wellModel
{

protected:

    // Protected Data

        //- Model name
        word name_;

        //- Phase names
        wordList phaseNames_;

        //- Model Dict
        dictionary wellsProperties_;

        //- Const access to the rock object
        const RockType& rock_;

        //- Const access to main pressure field
        const volScalarField& p_;

        //- List of pointers to well groups
        PtrList<objectRegistry> groups_;

        //- Per-phase source describers for all wells
        HashTable<autoPtr<wellSource<RockType, nPhases>>> sources_;

        //- Per-phase sources for all wells
        HashTable<fvScalarMatrix> matTable_;
        
        //- List of poitners to wells
        PtrList<well<RockType, nPhases>> wells_;


    // Protected Member Functions

        //- Collect well group names and return wells as list of dict entries
        PtrList<entry> assembleWellGroups();

        //- Create well objects
        void createWells();

public:

    //- Runtime type information
    ClassName("wellModel");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        wellModel,
        dictionary,
        (
            const word& name,
            const dictionary& transportProperties,
            const dictionary& wellsProperties,
            const RockType& rock
        ),
        (name, transportProperties, wellsProperties, rock)
    );

    // Constructors

        //- Construct from components
        wellModel
        (
            const word& name,
            const dictionary& transportProperties,
            const dictionary& wellsProperties,
            const RockType& rock
        );

        //- Construct from copy
        wellModel(const wellModel&) = delete;

    //- Destructor
    virtual ~wellModel();

    // Selectors

        //- Return a reference to the selected wellModel
        static autoPtr<wellModel> New
        (
            const word& name,
            const dictionary& transportProperties,
            const dictionary& wellsProperties,
            const RockType& rock
        );

    // Member Functions

        //- Return well model name
        const word& name() const {
            return name_;
        }

        //- Return well Model dictionary
        const dictionary& dict() const {
            return wellsProperties_;
        }

        //- Return const-access to the rock
        const RockType& rock() const {
            return rock_;
        }

        //- Const-access to list of well pointers
        const PtrList<well<RockType, nPhases>>& wells() const {
            return wells_;
        }

        ////- Return explicit portion of the source for given phase
        //const scalarField& sourceExplicitPart(const word& phase) const;

        ////- Return implicit portion of the source for given phase
        //const scalarField& sourceImplicitPart(const word& phase) const;

        //- Return Full source matrix for given phase
        tmp<fvScalarMatrix> source(const word& phase) const
        {
            return matTable_[phase];
        }

        //- Update well model sources
        virtual void correct() = 0;

    // Member Operators

        //- Disallow default bitwise assignment
        wellModel& operator=(const wellModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#define makeTemplatedWellModel\
(baseType,childType,TargName,ns,nPhases,childTArg)                             \
                                                                               \
    /* Define child type name */                                               \
    template<>                                                                 \
    const ::Foam::word ns::childType<childTArg, nPhases>::                     \
        typeName(ns::childType<childTArg, nPhases>::typeName_());              \
                                                                               \
    /* Define child debug switch */                                            \
    template<>                                                                 \
    int ns::childType<childTArg, nPhases>::debug                               \
    (                                                                          \
        ::Foam::debug::debugSwitch                                             \
        (ns::childType<childTArg, nPhases>::typeName_(), 0)                    \
    );                                                                         \
    /* Register the debug switch */                                            \
    template<>                                                                 \
    const Foam::RegisterDebugSwitch<ns::childType<childTArg, nPhases>>         \
    Foam::RegisterDebugSwitch<ns::childType<childTArg, nPhases>>               \
        ::registerDebugSwitch(ns::childType<childTArg, nPhases>::typeName_()); \
                                                                               \
    /* Add child to the base's run-time selection tabe */                      \
    namespace ns {                                                             \
    baseType<childTArg, nPhases>::                                             \
    adddictionaryConstructorToTable<childType<childTArg, nPhases>>             \
    add##phaseType##TargName##argNames##ConstructorTo##baseType##Targ##Table_; \
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "wellModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
