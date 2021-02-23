/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    twoPhaseIsoImpesFoam

Description
    Transient solver for two-phase BlackOil flow in isotropic porous media using
    the IMPES method.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "relPermModel.H"
#include "capPressModel.H"
#include "wellModel.H"
#include "impesControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createTimeControls.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createSaturationFields.H"
    #include "readTimeControls.H"

    impesControl<RockType,2> impes
    (
        "IMPES",
        phi,
        wordList(phaseNames),
        rockPtr()
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (impes.loop(runTime))
    {
        // Setup well sources
        wModel->correct();

        // Adjust timestep
        if (adjustTimeStep)
        {
            runTime.setDeltaT
            (
                impes.deltaT
                (
                    maxCo,
                    impes.deltaTFromAlphaEquation
                    (
                        phic,
                        wModel->explicitSource(canPhasePtr->name()),
                        alphaStorage
                    ),
                    maxDeltaT
                )
            );
        }
        Info << "deltaT: " << runTime.deltaTValue() << endl;

        Info<< "Time = " << runTime.timeName() << nl << endl;


        //- Solve saturation equation and update saturation fields
        #include "updateFields.H"
        #include "alphaEqn.H"
        #include "updateSaturationFields.H"

        //- Solve pressure equation and update pressure-related fluxes
        #include "pEqn.H"
        #include "continuityErrs.H"
        

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
