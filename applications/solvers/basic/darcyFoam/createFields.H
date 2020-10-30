// Rock type to use
using PermType = Isotropic;
using RockType = rock<PermType>; // isotropic medium


Info<< "Reading field p\n" << endl;

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

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading rockProperties\n" << endl;

IOdictionary rockProperties
(
    IOobject
    (
        "rockProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading wellsProperties\n" << endl;

IOdictionary wellsProperties
(
    IOobject
    (
        "wellsProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
word phaseName = wordList(transportProperties.lookup("phases"))[0];

autoPtr<phase> phasePtr = phase::New
(
    phaseName,
    mesh,
    transportProperties,
    mixtureType::singlePhase
);
autoPtr<RockType>  rockPtr  = RockType::New("rock", mesh, rockProperties);

PermType Tr
(
    IOobject
    (
        "Tr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    ( "Tr", dimArea/dimViscosity/dimDensity, 1.0),
    "zeroGradient"
);

forAll(Tr, ci)
{
    Tr[ci] = rockPtr->K()[ci]/phasePtr->mu()[ci];
}

// The well Model
autoPtr<wellModel<RockType,1>> wModel = wellModel<RockType,1>::New
    ("wellModel", transportProperties, wellsProperties, rockPtr());

//#include "createFvOptions.H"
