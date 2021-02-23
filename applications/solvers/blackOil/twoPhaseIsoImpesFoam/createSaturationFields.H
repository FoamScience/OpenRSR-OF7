// Relative Permeability
krModel->correct();

word krcName = relPermModel<RockType,2>::krName(canPhasePtr->name());
word krnName = relPermModel<RockType,2>::krName(nonCanPhasePtr->name());
word dkrcName = relPermModel<RockType,2>::dkrName(canPhasePtr->name(), canPhasePtr->name());
word dkrnName = relPermModel<RockType,2>::dkrName(nonCanPhasePtr->name(), canPhasePtr->name());

const auto& krc = krModel()[krcName];
const auto& krn = krModel()[krnName];
const auto& dkrcdS = krModel()[dkrcName];
const auto& dkrndS = krModel()[dkrnName];

surfaceScalarField krcf ("krcf",fvc::interpolate(krc, krc.name()));
surfaceScalarField krnf ("krnf",fvc::interpolate(krn, krn.name()));
surfaceScalarField dkrcfdS ("dkrcfdS",fvc::interpolate(dkrcdS, dkrcdS.name()));
surfaceScalarField dkrnfdS ("dkrnfdS",fvc::interpolate(dkrndS, dkrndS.name()));

// Face-interpolations of fluid properties

surfaceScalarField munf ("munf",fvc::interpolate(mun,"mu"));
surfaceScalarField mucf ("mucf",fvc::interpolate(mun,"mu"));
surfaceScalarField rhonf ("rhonf",fvc::interpolate(rhon,"rho"));
surfaceScalarField rhocf ("rhocf",fvc::interpolate(rhon,"rho"));
surfaceScalarField Bcf ("Bcf",fvc::interpolate(Bc,"rFVF"));
surfaceScalarField Bnf ("Bnf",fvc::interpolate(Bn,"rFVF"));

// Phase Mobility and fractional flow
surfaceScalarField Mnf ("Mnf", Bnf*Kf*krnf/munf);
surfaceScalarField Lnf ("Lnf", rhonf*Mnf);	
surfaceScalarField Mcf ("Mcf", Bcf*Kf*krcf/mucf);
surfaceScalarField Lcf ("Lcf", rhocf*Mcf);
surfaceScalarField Mf ("Mf", Mnf+Mcf);
surfaceScalarField Lf ("Lf", Lnf+Lcf);
surfaceScalarField Fcf ("Fcf", Mcf/Mf);

// Capillary Pressure
pcModel->correct();
surfaceScalarField phiG("phiG",(Lf * g) & mesh.Sf());
word dpcName = capPressModel<RockType,2>::dpcName
(
    canPhasePtr->name(), canPhasePtr->name()
);
surfaceScalarField phiPc
(
    "phiPc",
    Mcf * fvc::interpolate(pcModel()[dpcName],"pc")
    * fvc::snGrad(Sc) * mesh.magSf()
);

// Well model (Needs to be last registered thing to do)
autoPtr<wellModel<RockType,2>> wModel = wellModel<RockType,2>::New
(
    "wModel", transportProperties, wellsProperties, rockPtr()
);
