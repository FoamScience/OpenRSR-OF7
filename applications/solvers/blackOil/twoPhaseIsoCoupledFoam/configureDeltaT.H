scalar adjustedDeltaT = 100;

if (adjustTimeStep)
{
    if (runTime.deltaTValue() != runTime.startTime().value())
    {
        scalarField dpf = p.internalField() - p.oldTime();
        maxDp = gMax(cmptMag(dpf));
        scalarField dScf = Sc.internalField() - Sc.oldTime();
        maxDSc = gMax(cmptMag(dScf));
    
        if (maxDp > vSmall and maxDSc > vSmall)
        {
            scalar dtp = runTime.deltaTValue() * dpLim / maxDp;
            scalar dtSc = runTime.deltaTValue() * dScLim / maxDSc;
            adjustedDeltaT = min(dtp, dtSc);
            runTime.setDeltaT(adjustedDeltaT);
        } else {
            runTime.setDeltaT(100);
        }
    } else {
        runTime.setDeltaT(100);
    }
}

Info << "deltaT: " << runTime.deltaTValue() << endl;
