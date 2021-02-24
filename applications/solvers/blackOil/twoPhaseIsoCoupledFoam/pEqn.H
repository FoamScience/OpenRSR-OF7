{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(-Mf, p) + fvc::div(phiG)
        + fvc::div(phiPc)
        + wModel->source(canPhasePtr->name())
        + wModel->source(nonCanPhasePtr->name())
    );

    fvScalarMatrix ScInPEqn
    (
        fvm::ddt(pStorage, Sc)
    );

    pEqn.setReference(pRefCell, pRefValue);

    //pEqn.solve();
    cs->insertEquation(p.name(), p.name(), pEqn);
    cs->insertEquation(p.name(), Sc.name(), ScInPEqn);

    phiP = pEqn.flux();

}
