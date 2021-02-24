// Porous Medium Setup
using PermType = Isotropic;
using SurfacePermType = surfaceScalarField;
using RockType = rock<PermType>;

Info << "Reading pressure field p" << endl;

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << nl << "Reading transportProperties" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Collect phase names
Pair<word> phaseNames(transportProperties.lookup("phases"));

// Canonical phase
autoPtr<phase> canPhasePtr = phase::New(phaseNames.first(), mesh,transportProperties);
volVectorField& Uc = canPhasePtr->U();
volScalarField& Sc = canPhasePtr->alpha();
surfaceScalarField& phic = canPhasePtr->phi();
const auto& rhoc = canPhasePtr->rho();
const auto& muc = canPhasePtr->mu();
const auto& Bc = canPhasePtr->BModel().rFVF();

// Non-canonical phase
autoPtr<phase> nonCanPhasePtr = phase::New(phaseNames.second(), mesh,transportProperties);
volVectorField& Un = nonCanPhasePtr->U();
volScalarField& Sn = nonCanPhasePtr->alpha();
surfaceScalarField& phin = nonCanPhasePtr->phi();
const auto& rhon = nonCanPhasePtr->rho();
const auto& mun = nonCanPhasePtr->mu();
const auto& Bn = nonCanPhasePtr->BModel().rFVF();

// Viscosities ratio
volScalarField Mmu(muc/mun);

// Rock
autoPtr<RockType> rockPtr = RockType::New("rock", mesh, transportProperties);
const PermType& K = rockPtr->K();
const volScalarField& porosity = rockPtr->porosity();
SurfacePermType Kf("Kf", fvc::interpolate(K,"K"));

// Relative permeability / Capillary pressure models 
word krModelName = "krModel<"+phaseNames.first()+","+phaseNames.second()+">";
autoPtr<relPermModel<RockType,2>> krModel = relPermModel<RockType,2>::New
(
    krModelName,
    transportProperties,
    rockPtr()
);

word pcModelName = "pcModel<"+phaseNames.first()+","+phaseNames.second()+">";
autoPtr<capPressModel<RockType,2>> pcModel = capPressModel<RockType,2>::New
(
    pcModelName,
    transportProperties,
    rockPtr()
);

Info << nl << "Reading field U" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Un + Uc
);

#include "createPhi.H"
surfaceScalarField phiP = phi;

Info << nl << "Reading wellsProperties" << endl;

IOdictionary wellsProperties
(
    IOobject
    (
        "wellsProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Equation terms coeffs
//- Storage term coeff for saturation eqn
volScalarField alphaStorage
(
    IOobject
    (
        "porosity*conPhase.rFVF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    porosity * Bc
);
//- Storage term coff for pressure eqn
volScalarField pStorage
(
    IOobject
    (
        "porosity*(conPhase.rFVF+nonconPhase.rFVF)",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    porosity * (Bc - Bn)
);

label  pRefCell  = 0;
scalar pRefValue = 0;
setRefCell(p, mesh.solutionDict(), pRefCell, pRefValue);
