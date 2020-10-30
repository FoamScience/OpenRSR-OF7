{

    Un.correctBoundaryConditions();
    Uc.correctBoundaryConditions();
    forAll(mesh.boundary(),patchi)
    {
        if (isA< fixedValueFvPatchField<vector> >(Un.boundaryField()[patchi]))
        {
            phin.boundaryFieldRef()[patchi] = Un.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        }
        if (isA< fixedValueFvPatchField<vector> >(Uc.boundaryField()[patchi]))
        {
            phic.boundaryFieldRef()[patchi] = Uc.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        }
    }
  
    fvScalarMatrix ScEqn
        (
            porosity*fvm::ddt(Sc) + fvc::div(phic) 
            // event source terms
            // - sourceTerm
        );
    ScEqn.source() += wModel->explicitSource(canPhasePtr->name());


    ScEqn.solve();

    Info << Sc.name() << " Min = " << gMin(Sc) << " Max = " << gMax(Sc) << endl;

}
