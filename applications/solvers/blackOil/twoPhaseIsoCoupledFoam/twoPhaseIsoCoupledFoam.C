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
    twoPhaseIsoCoupledFoam

Description
    Transient solver for two-phase flow in isotropic porous media using
    an AIM-like method where we treat capillary and gravity fluxes explicitely
    and solve the equations simulataniously

\*---------------------------------------------------------------------------*/

// Include this first to avoid namespace clashes with Foam
#include "coupledSolver.H"

#include "fvCFD.H"
#include "relPermModel.H"
#include "capPressModel.H"
#include "wellModel.H"

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
    #include "createSaturationFields.H"
    #include "readTimeControls.H"

    //- Create the coupled solver
    word solverName = Sc.name() + "-" + p.name();
    auto cs = autoPtr<coupledSolver>
    (
        new coupledSolver
        (
            solverName,
            mesh.name(),
            runTime,
            mesh.solutionDict()
        )
    );

    //- Setup the coupled system
    cs->insertMesh(mesh);
    cs->insertField(Sc);
    cs->insertField(p);

    //- Initiate timeStep
    #include "initTimeStep.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    wModel->correct();
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Update fields and add matrices to coupled system
        #include "updateFields.H"
        #include "alphaEqn.H"
        #include "pEqn.H"

        // Solve the coupled system
        cs->solve();

        // Update saturation fields
        #include "updateSaturationFields.H"

        //- Solve pressure equation
        Info << Sc.name() << " Min = " << gMin(Sc) 
            << " Max = " << gMax(Sc) << endl;

        #include "configureDeltaT.H"
        
        wModel->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
