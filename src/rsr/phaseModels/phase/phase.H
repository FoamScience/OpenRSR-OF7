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
    Foam::phase

Description
    A template base class for fluid phases.

SourceFiles
    phase.C

\*---------------------------------------------------------------------------*/

#ifndef phase_H
#define phase_H

#include "phaseCore.H"
#include "FVFModel.H"
#include "UniformityTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class phase Declaration
\*---------------------------------------------------------------------------*/

template<class CompressibilityType, class ViscosityType>
class phase
:
    public phaseCore
{
protected:

    // Protected Data

        //- Phase density
        CompressibilityType rho_;

        //- Phase viscosity
        ViscosityType mu_;

        //- FVF Model
        autoPtr<FVFModel<CompressibilityType>> BModel_;

public:

    //- Runtime type information
    ClassName("phase");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        phase,
        dictionary,
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& transportProperties,
            const mixtureType& mT
        ),
        (name, mesh, transportProperties, mT)
    );

    // Constructors

        // Construct from components
        phase
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& transportProperties,
            const mixtureType& mT = mixtureType::multiPhase
        );

        // Construct from copy
        phase(const phase&);


    //- Destructor
    virtual ~phase();

    // Selectors

        //- Return a reference to the selected IncompressiblePhase
        static autoPtr<phase> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& transportProperties,
            const mixtureType& mT = mixtureType::multiPhase
        );


    // Member Functions

        //- Return const-ref to density field
        const CompressibilityType& rho() const {
            return rho_;
        }

        //- Return const-ref to viscosity field
        const CompressibilityType& mu() const {
            return mu_;
        }

        //- Return const-access to FVF Model
        const FVFModel<CompressibilityType>& BModel()
        {
            return BModel();
        }

        //- Correct phase fields
        virtual void correct() = 0;

    // Member Operators

        //- Disallow default bitwise assignment
        phase& operator=(const phase&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "phase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
