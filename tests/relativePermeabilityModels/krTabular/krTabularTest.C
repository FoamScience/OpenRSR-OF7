#include "IsotropyTypes.H"
#include "catch.H"
#include "autoPtr.H"
#include "fvCFD.H"
#include "relPermModel.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


using namespace Foam;

SCENARIO("Relative permeability curves in tabular form", "[Virtual]")
{
    GIVEN("Valid mesh, two phases, and a valid transportProperties")
    {
        #include "createTestTimeAndMesh.H"
        #include "createTestBlackoilPhase.H"
        #include "createTestIsoRock.H"

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

        dictionary krDict("krModel<water,oil>");
        krDict.add("type", "tabular");
        dictionary krData("krData");
        krData.add<fileName>("file", "testData/kr.dat");
        krData.add("interpolationType", "linear"); // Default
        krDict.add("krData", krData);
        
        
        transportProperties.add(word("krModel<water,oil>"), krDict);

        auto krModel = relPermModel<iRock, 2>::New
        (
            word("krModel<water,oil>"),
            transportProperties,
            rkPtr()
        );

        WHEN("correct() is called on a tabular relPerm model")
        {
            waterPtr->alpha()[0] = 0.3;
            krModel->correct();

            const auto& kr1 = krModel()[waterPtr->name()+".kr"];
            const auto& kr2 = krModel()[oilPtr->name()+".kr"];
            const auto& dkr1 = krModel()[waterPtr->name()+".dkrdS(water)"];
            const auto& dkr2 = krModel()[oilPtr->name()+".dkrdS(water)"];

            THEN("Kr values must be consistent with input data")
            {
                REQUIRE_THAT
                (
                    kr1[0], Catch::Matchers::WithinAbs(0.00411522, 1e-6)
                );
                REQUIRE_THAT
                (
                    kr2[0], Catch::Matchers::WithinAbs(0.57248677, 1e-6)
                );
            }

            waterPtr->alpha()[0] = 0.88;
            krModel->correct();

            THEN("dKr values must be correct near end-points")
            {
                REQUIRE_THAT
                (
                    dkr1[0], Catch::Matchers::WithinAbs(4.0544217, 1e-6)
                );
                REQUIRE_THAT
                (
                    dkr2[0], Catch::Matchers::WithinAbs(-0.0122448, 1e-6)
                );
            }
        }
    }
}

// ************************************************************************* //
