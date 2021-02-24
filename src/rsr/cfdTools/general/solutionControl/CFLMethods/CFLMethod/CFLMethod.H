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
    Foam::CFLMethod

Description
    A template base class for CFL calculation methods

SourceFiles
    CFLMethod.C

\*---------------------------------------------------------------------------*/

#ifndef CFLMethod_H
#define CFLMethod_H

#include "rock.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class CFLMethod Declaration
\*---------------------------------------------------------------------------*/

template<class RockType, int nPhases>
class CFLMethod
{
protected:

    // Protected Data

        //- Method name
        word name_;

        //- Method Dict
        dictionary cflDict_;

        //- Const access to Flux field
        const surfaceScalarField& phi_;

        //- Const access to the rock object
        const RockType& rock_;

        //- CFL number
        scalarField CFLNo_;

public:

    //- Runtime type information
    ClassName("CFLMethod");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        CFLMethod,
        dictionary,
        (
            const word& name,
            const dictionary& algorithmProperties,
            const surfaceScalarField& phi,
            const wordList& phaseNames,
            const RockType& rock
        ),
        (name, algorithmProperties, phi, phaseNames, rock)
    );

    // Constructors

        //- Construct from components
        CFLMethod
        (
            const word& name,
            const dictionary& algorithmProperties,
            const surfaceScalarField& phi,
            const wordList& phaseNames,
            const RockType& rock
        );

        //- Construct from copy
        CFLMethod(const CFLMethod&) = delete;

    //- Destructor
    virtual ~CFLMethod();

    // Selectors

        //- Return a reference to the selected CFLMethod
        static autoPtr<CFLMethod> New
        (
            const word& name,
            const dictionary& algorithmProperties,
            const surfaceScalarField& phi,
            const wordList& phaseNames,
            const RockType& rock
        );

    // Member Functions

        //- Return model name
        const word& name() const {
            return name_;
        }

        //- Return Model dictionary
        const dictionary& dict() const {
            return cflDict_;
        }

        //- Return const-access to the rock
        const RockType& rock() const {
            return rock_;
        }

        //- Return CFL number
        const scalarField& CFLNo() const
        {
            return CFLNo_;
        }

        //- Update CFL number
        virtual void correct() = 0;

    // Member Operators

        //- Disallow default bitwise assignment
        CFLMethod& operator=(const CFLMethod&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#define maketwoPhaseCFLMethodTemplatedModel\
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
    #include "CFLMethod.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
