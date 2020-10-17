#include "IsotropyTypes.H"
#include "catch.H"
#include "fvCFD.H"
#include "wellModel.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

SCENARIO("Peaceman model for well source calculation", "[Virtual]")
{
    GIVEN("Valid mesh, kr and pc models, and a wellsProperties dictionary")
    {
        #include "createTestTimeAndMesh.H"
        #include "createTestBlackoilPhase.H"
        #include "createTestIsoRock.H"
        #include "readGravitationalAcceleration.H"

        dictionary transportProperties;
        transportProperties.add<wordList>("phases", {"water", "oil"});
        dictionary rockProperties;

        volScalarField p
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("p", dimPressure, 0.0)
        );

        createTestBlackoilPhase(water, 1.0, 1e-3, multiPhase);
        createTestBlackoilPhase(oil, 1.0, 1e-5, multiPhase); 

        createTestIsoRock(rk, 1e-12, 0.2, 1e-6);

        dictionary krDict("krModel<water,oil>");
        krDict.add("type", "BrooksCorey");
        krDict.add<dimensionedScalar>
        (
            "water.alphaIrr",
            dimensionedScalar("water.alphaIrr", dimless, 0.2)
        );
        krDict.add<dimensionedScalar>
        (
            "oil.alphaRes",
            dimensionedScalar("oil.alphaRes", dimless, 0.1)
        );
        krDict.add<dimensionedScalar>
        (
            "water.m",
            dimensionedScalar("mater.m", dimless, 2)
        );
        krDict.add<dimensionedScalar>
        (
            "oil.m",
            dimensionedScalar("oil.m", dimless, 2)
        );
        krDict.add<dimensionedScalar>
        (
            "water.krMax",
            dimensionedScalar("water.krMax", dimless, 1.0)
        );
        krDict.add<dimensionedScalar>
        (
            "oil.krMax",
            dimensionedScalar("oil.krMax", dimless, 0.9)
        );
        
        transportProperties.add(word("krModel<water,oil>"), krDict);

        auto krModel = relPermModel<iRock, 2>::New
        (
            word("krModel<water,oil>"),
            transportProperties,
            rkPtr()
        );

        dictionary pcDict("pcModel<water,oil>");
        pcDict.add("type", "BrooksCorey");
        pcDict.add<dimensionedScalar>
        (
            "water.alpha.PcMin",
            dimensionedScalar("water.alpha.PcMin", dimless, 0.15)
        );
        pcDict.add<dimensionedScalar>
        (
            "water.alpha.PcMax",
            dimensionedScalar("water.alpha.PcMax", dimless, 0.85)
        );
        pcDict.add<dimensionedScalar>
        (
            "n",
            dimensionedScalar("n", dimless, 0.238)
        );
        pcDict.add<dimensionedScalar>
        (
            "pc0",
            dimensionedScalar("pc0", dimPressure, 30)
        );
        
        transportProperties.add(word("pcModel<water,oil>"), pcDict);

        auto pcModel = capPressModel<iRock, 2>::New
        (
            word("pcModel<water,oil>"),
            transportProperties,
            rkPtr()
        );

        //FatalError.dontThrowExceptions();

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
        WHEN("Peaceman's well model is instantiated")
        {
            auto wModel = wellModel<iRock, 2>::New
            (
                "wModel", transportProperties, wellsProperties, rkPtr()
            );

            THEN("Well Cells and cell faces are selected correctly")
            {
                forAll(wModel->wells(), wi)
                {
                    CHECK
                    (
                        wModel->wells()[wi].cellIDs() 
                        == labelList{1,2,3,4,5,6,7,8}
                    );
                    CHECK
                    (
                        wModel->wells()[wi].faceIDs()
                        == labelList{1,2,3,4,5,6,7}
                    );
                }
            }
        }
    }
}

// ************************************************************************* //
