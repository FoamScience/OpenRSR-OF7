{
    // canPhase eqn + nonCanPhase eqn
    fvScalarMatrix pEqn
    (
        fvc::ddt(pStorage, Sc)
        + fvm::laplacian(-Mf, p) + fvc::div(phiG)
        + fvc::div(phiPc)
        + wModel->source(canPhasePtr->name())
        + wModel->source(nonCanPhasePtr->name())
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve();

    // Update fluxes and velocities

    phiP = pEqn.flux();

    phi = phiP+phiG+phiPc;

    phic = Fcf*phiP + (Lcf/Lf)*phiG + phiPc;
    phin = phi - phic;

    U = fvc::reconstruct(phi);
    U.correctBoundaryConditions();

    Uc = fvc::reconstruct(phic);
    Un = U-Uc;

    Uc.correctBoundaryConditions();  
    Un.correctBoundaryConditions();

}
