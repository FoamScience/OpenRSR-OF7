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
#include "impesControl.H"
#include "fvcDiv.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::impesControl<RockType, nPhases>::impesControl
(
    const word& name,
    const surfaceScalarField& phi,
    const wordList& phaseNames,
    const RockType& rock,
    const word& algorithmName
)
:
    fluidSolutionControl(const_cast<fvMesh&>(rock.mesh()), algorithmName),
    singleRegionConvergenceControl
    (
        static_cast<singleRegionSolutionControl&>(*this)
    ),
    dSMax_(1),
    CFLMethod_
    (
        CFLMethod<RockType, nPhases>::New
        (
            name, fluidSolutionControl::dict(), phi, phaseNames, rock
        )
    )
{
    if (!read())
    {
        FatalErrorInFunction
            << "Could not read solution parameters."
            << exit(FatalError);
    }
    printResidualControls();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::impesControl<RockType, nPhases>::~impesControl() {}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
bool Foam::impesControl<RockType,nPhases>::read()
{
    const dictionary& solutionDict = fluidSolutionControl::dict();

    // Read max dS
    dSMax_ = readScalar(solutionDict.lookup("dSMax"));

    return fluidSolutionControl::read() && readResidualControls();
}


template<class RockType, int nPhases>
bool Foam::impesControl<RockType,nPhases>::run(Time& time)
{
    read();

    if (!endIfConverged(time))
    {
        storePrevIterFields();
    }

    return time.run();
}


template<class RockType, int nPhases>
bool Foam::impesControl<RockType,nPhases>::loop(Time& time)
{
    read();

    if (!endIfConverged(time))
    {
        storePrevIterFields();
    }

    return time.loop();
}

template<class RockType, int nPhases>
template<class CoeffFieldType>
Foam::scalar Foam::impesControl<RockType,nPhases>::deltaTFromAlphaEquation
(
    const surfaceScalarField& phi,
    const scalarField& explicitSource,
    const CoeffFieldType& coefft
) const
{
    // Approximate ddt(phase.alpha)
    scalarField dSdt = mag
    (
        (-fvc::div(phi))->internalField()
        + explicitSource 
    ) / coefft;
    return dSMax_/(gMax(dSdt)+vSmall);
}

template<class RockType, int nPhases>
Foam::scalar Foam::impesControl<RockType,nPhases>::deltaT
(
    scalar maxCo,
    scalar dtForAlphaEqn,
    scalar maxDeltaT
)
{
    scalar gMaxCFLNo = -1;
    scalar goodDeltaT = -1;

    CFLMethod_->correct();

    gMaxCFLNo = gMax(CFLMethod_->CFLNo());
    scalar maxCFLDeltaT = maxCo/(gMaxCFLNo + vSmall);
    goodDeltaT = min(min(maxCFLDeltaT, 1 + 0.1*maxCFLDeltaT), 1.2);
    Info << goodDeltaT << endl;
    return min
    (
        dtForAlphaEqn,
        min(goodDeltaT*mesh().time().deltaTValue(), maxDeltaT)
    );
}
// ************************************************************************* //
