{

    //- Saturation equation for canonical phase
    fvScalarMatrix ScEqn
    (
        fvm::ddt(alphaStorage, Sc) + fvc::div(phicG) + fvc::div(phiPc) 
    );

    //- Pressure coeffs in canonical phase eqn
    fvScalarMatrix pInScEqn
    (
        fvm::laplacian(-Mcf, p)
        + wModel->source(canPhasePtr->name())
    );

    //pInScEqn.setReference(pRefCell, pRefValue);

    // Get conservative canonical flux
    //phiP = pInScEqn.flux();

    // Insert equations into coupled system
    cs->insertEquation(Sc.name(), Sc.name(), ScEqn);
    cs->insertEquation(Sc.name(), p.name(), pInScEqn);

}
