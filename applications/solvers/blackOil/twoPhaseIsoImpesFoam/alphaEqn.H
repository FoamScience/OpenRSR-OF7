{

    //- Place-holder code for "water aquifer" modeling
    Un.correctBoundaryConditions();
    Uc.correctBoundaryConditions();
    forAll(mesh.boundary(),patchi)
    {
        if (isA< fixedValueFvPatchField<vector> >(Un.boundaryField()[patchi]))
        {
            phin.boundaryFieldRef()[patchi] = 
                Un.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        }
        if (isA< fixedValueFvPatchField<vector> >(Uc.boundaryField()[patchi]))
        {
            phic.boundaryFieldRef()[patchi] = 
                Uc.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        }
    }
  
    //- Saturation equation
    fvScalarMatrix ScEqn
    (
        fvm::ddt(alphaTemp,Sc) + fvc::div(phic) 
    );
    ScEqn.source() += wModel->explicitSource(canPhasePtr->name());


    //- The solver should be hard-coded to use diagonal
    ScEqn.solve();

    Info << Sc.name() << " Min = " << gMin(Sc) << " Max = " << gMax(Sc) << endl;

}
