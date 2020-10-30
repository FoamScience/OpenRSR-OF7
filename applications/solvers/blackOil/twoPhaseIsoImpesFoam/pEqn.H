{
    fvScalarMatrix pEqn
        (
            fvm::laplacian(-Mf, p) + fvc::div(phiG)
            // capillary term
            + fvc::div(phiPc)//*activateCapillarity
            + wModel->source(canPhasePtr->name())
            + wModel->source(nonCanPhasePtr->name())
            // event source terms
            //- sourceTerm
        );

    pEqn.solve();

    phiP = pEqn.flux();

    phi = phiP+phiG+phiPc;//*activateCapillarity;

    phic = Fcf*phiP + (Lcf/Lf)*phiG + phiPc;//*activateCapillarity;
    phin = phi - phic;

    U = fvc::reconstruct(phi);
    U.correctBoundaryConditions();

    Uc = fvc::reconstruct(phic);
    Un = U-Uc;

    Uc.correctBoundaryConditions();  
    Un.correctBoundaryConditions();

}
