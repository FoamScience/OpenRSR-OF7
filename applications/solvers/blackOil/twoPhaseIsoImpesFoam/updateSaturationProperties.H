Sn = 1 - Sc;

//- relative permeability computation
krModel->correct();
krcf = fvc::interpolate(krc,krc.name());
krnf = fvc::interpolate(krn, krn.name());
dkrcfdS = fvc::interpolate(dkrcdS,dkrcdS.name());
dkrnfdS = fvc::interpolate(dkrndS,dkrndS.name());

munf = fvc::interpolate(mun,"mu");
mucf = fvc::interpolate(muc,"mu");
rhonf = fvc::interpolate(rhon,"rho");
rhocf = fvc::interpolate(rhoc,"rho");

//- mobility computation 
Mnf = Kf*krnf/munf;
Lnf = rhonf*Kf*krnf/munf;	
Mcf = Kf*krcf/mucf;
Lcf = rhocf*Kf*krcf/mucf;
Mf = Mnf+Mcf;
Lf = Lnf+Lcf;
Fcf = Mcf/Mf;

//- compute fluxes
pcModel().correct();
phiPc = Mcf * fvc::interpolate(pcModel()[dpcName],"pc")* fvc::snGrad(Sc) * mesh.magSf();
phiG = (Lf * g) & mesh.Sf();
