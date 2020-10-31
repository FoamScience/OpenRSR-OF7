{
    fvScalarMatrix pEqn
        (
            fvm::laplacian(-Mf, p) + fvc::div(phiG)
            + fvc::div(phiPc)
            + wModel->source(canPhasePtr->name())
            + wModel->source(nonCanPhasePtr->name())
        );

    pEqn.solve();

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
