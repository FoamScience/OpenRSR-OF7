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

#include "CompressiblePhase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ViscosityType>
Foam::autoPtr<Foam::CompressiblePhase<ViscosityType>>
Foam::CompressiblePhase<ViscosityType>::New
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const mixtureType& mT
)
{
    const word phaseType = dict.subDict(name).lookupOrDefault<word>
    (
      "phaseType",
      "blackoil"
    );
                                                                                
    Info<< "Selecting phase type " << phaseType << " for " << name << endl;             
    
    typename dictionaryConstructorTable::iterator cstrIter =                    
        dictionaryConstructorTablePtr_->find(phaseType);                        
                                                                                
    if (cstrIter == dictionaryConstructorTablePtr_->end())                      
    {                                                                           
        FatalErrorInFunction                                                    
            << "Unknown Compressible Phase type "                                  
            << phaseType << nl << nl                                            
            << "Valid Compressible Phase types:" << endl                           
            << dictionaryConstructorTablePtr_->sortedToc()                      
            << exit(FatalError);                                                
    }                                                                           
                                                                                
    return autoPtr<CompressiblePhase>                                             
    ( cstrIter()(name, mesh, dict, mT)  );          
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ViscosityType>
Foam::CompressiblePhase<ViscosityType>::CompressiblePhase
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const mixtureType& mT
)
:
    phase<Compressible, ViscosityType>(name, mesh, dict, mT)
{
}


template<class ViscosityType>
Foam::CompressiblePhase<ViscosityType>::CompressiblePhase
(
    const CompressiblePhase& ph
)
:
    phase<Compressible, ViscosityType>(ph)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ViscosityType>
Foam::CompressiblePhase<ViscosityType>::~CompressiblePhase()
{}

// ************************************************************************* //
