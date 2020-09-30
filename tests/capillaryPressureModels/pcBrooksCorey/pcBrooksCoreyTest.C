#include "IsotropyTypes.H"
#include "catch.H"
#include "autoPtr.H"
#include "fvCFD.H"
#include "capPressModel.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


using namespace Foam;

SCENARIO("Brooks Corey Calculation of capillary pressure", "[Virtual]")
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

        dictionary pcDict("pcModel<water,oil>");
        pcDict.add("type", "BrooksCorey");
        pcDict.add<dimensionedScalar>
        (
            "water.alpha.PcMin",
            dimensionedScalar("water.alpha.PcMin", dimless, 0.2)
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

        WHEN("correct() is called on a B-C capPress model")
        {
            forAll(mesh.C(), ci)
            {
                waterPtr->alpha()[ci] = 0.2 + ci*(0.85-0.2)/(mesh.nCells()-1);
            }
            waterPtr->alpha()[0] += 0.001;
            pcModel->correct();
            THEN("Pc values must be consistent with model expression")
            {
                const std::vector<scalar> pcManual = {
                    140.150494, 50.609374, 42.912699, 38.965128,
                    36.386534, 34.504529, 33.039310, 31.849134,
                    30.852868, 30.0
                };
                
                const std::vector<scalar> dpcManual = {
                    -33355.81762, -166.7773549, -70.70692446, -42.80169505,
                    -29.97690644, -22.74113919, -18.14620587, -14.99359267,
                    -12.70900872, -10.98461538
                };

                const auto& pc = 
                    pcModel()[capPressModel<iRock,2>::pcName("water")];
                const auto& dpc = 
                    pcModel()[capPressModel<iRock,2>::dpcName("water","water")];

                REQUIRE_THAT
                (
                    pcManual, 
                    Catch::Matchers::Approx
                    (
                        std::vector<scalar> ( pc.begin(), pc.end() )
                    )
                );
                REQUIRE_THAT
                (
                    dpcManual, 
                    Catch::Matchers::Approx
                    (
                        std::vector<scalar> ( dpc.begin(), dpc.end() )
                    )
                );
            }
        }
    }
}

// ************************************************************************* //
