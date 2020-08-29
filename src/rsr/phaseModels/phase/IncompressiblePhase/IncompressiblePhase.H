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
    Foam::IncompressiblePhase

Description
    A template base class for fluid IncompressiblePhases.

SourceFiles
    IncompressiblePhase.C

\*---------------------------------------------------------------------------*/

#ifndef IncompressiblePhase_H
#define IncompressiblePhase_H

#include "UniformityTypes.H"
#include "phase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class IncompressiblePhase Declaration
\*---------------------------------------------------------------------------*/

template<class ViscosityType>
class IncompressiblePhase
:
    public phase<Incompressible, ViscosityType>
{
public:
    
    //- Runtime type information
    TypeName("IncompressiblePhase");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        IncompressiblePhase,
        dictionary,
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const mixtureType& mT
        ),
        (name, mesh, dict, mT)
    );

    // Constructors

        // Construct from components
        IncompressiblePhase
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const mixtureType& mT = mixtureType::multiPhase
        );

        // Construct from copy
        IncompressiblePhase (const IncompressiblePhase&);

    // Selectors

        //- Return a reference to the selected IncompressiblePhase
        static autoPtr<IncompressiblePhase> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const mixtureType& mT = mixtureType::multiPhase
        );


    //- Destructor
    virtual ~IncompressiblePhase();

    // Member Operators

        //- Disallow default bitwise assignment
        IncompressiblePhase& operator=(const IncompressiblePhase&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IncompressiblePhase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
