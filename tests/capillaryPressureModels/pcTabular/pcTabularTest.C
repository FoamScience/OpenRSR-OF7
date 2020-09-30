#include "IsotropyTypes.H"
#include "catch.H"
#include "autoPtr.H"
#include "fvCFD.H"
#include "capPressModel.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


using namespace Foam;

SCENARIO("Capillary pressure curve in tabular form", "[Virtual]")
{
    GIVEN("Valid mesh, two phases, and a valid transportProperties")
    {
        #include "createTestTimeAndMesh.H"
        #include "createTestBlackoilPhase.H"
        #include "createTestIsoRock.H"
        FatalError.dontThrowExceptions();

        dictionary transportProperties;
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

        dictionary pcDict("pcModel<water,oil>");
        pcDict.add("type", "tabular");
        dictionary pcData("pcData");
        pcData.add<fileName>("file", "testData/pc.dat");
        pcData.add("interpolationType", "linear"); // Default
        pcDict.add("pcData", pcData);
        
        
        transportProperties.add(word("pcModel<water,oil>"), pcDict);

        auto pcModel = capPressModel<iRock, 2>::New
        (
            word("pcModel<water,oil>"),
            transportProperties,
            rkPtr()
        );

        WHEN("correct() is called on a tabular relPerm model")
        {
            waterPtr->alpha()[0] = 0.3;
            pcModel->correct();

            const auto& pc = 
                pcModel()[capPressModel<iRock,2>::pcName("water")];
            const auto& dpc = 
                pcModel()[capPressModel<iRock,2>::dpcName("water","water")];

            THEN("Pc values must be consistent with input data")
            {
                REQUIRE_THAT
                (
                    pc[0], Catch::Matchers::WithinAbs(47.6162231, 1e-6)
                );
            }

            waterPtr->alpha()[0] = 0.202;
            pcModel->correct();

            THEN("dPc values must be correct near end-points")
            {
                REQUIRE_THAT
                (
                    dpc[0], Catch::Matchers::WithinAbs(-32433.89984, 1e-4)
                );
            }
        }
    }
}

// ************************************************************************* //
