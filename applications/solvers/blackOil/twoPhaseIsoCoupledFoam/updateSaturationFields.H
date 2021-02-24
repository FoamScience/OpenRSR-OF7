    phi = phiP+phiG+phiPc;

    phic = Fcf*phiP + (Lcf/Lf)*phiG + phiPc;
    phin = phi - phic;

    U = fvc::reconstruct(phi);
    U.correctBoundaryConditions();

    Uc = fvc::reconstruct(phic);
    Un = U-Uc;

    Uc.correctBoundaryConditions();  
    Un.correctBoundaryConditions();

Sn = 1 - Sc;

//- Update Kr
krModel->correct();
krcf = fvc::interpolate(krc,krc.name());
krnf = fvc::interpolate(krn, krn.name());
dkrcfdS = fvc::interpolate(dkrcdS,dkrcdS.name());
dkrnfdS = fvc::interpolate(dkrndS,dkrndS.name());

munf = fvc::interpolate(mun,"mu");
mucf = fvc::interpolate(muc,"mu");
rhonf = fvc::interpolate(rhon,"rho");
rhocf = fvc::interpolate(rhoc,"rho");
Bnf = fvc::interpolate(Bn,"rFVF");
Bcf = fvc::interpolate(Bc,"rFVF");

//- Update Moblities and fractional flow
Mnf = Bnf*Kf*krnf/munf;
Lnf = rhonf*Mnf;
Mcf = Bcf*Kf*krcf/mucf;
Lcf = rhocf*Mcf;
Mf = Mnf+Mcf;
Lf = Lnf+Lcf;
Fcf = Mcf/Mf;

//- "Auxilary" fluxes
pcModel->correct();
phiPc = Mcf * fvc::interpolate(pcModel()[dpcName],"pc") * fvc::snGrad(Sc)
    * mesh.magSf();
phiG = (Lf * g) & mesh.Sf();
phicG = (Lcf * g) & mesh.Sf();
