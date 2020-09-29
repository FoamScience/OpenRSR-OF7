#include "IsotropyTypes.H"
#include "catch.H"
#include "autoPtr.H"
#include "fvCFD.H"
#include "relPermModel.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


using namespace Foam;

SCENARIO("Brooks Corey Calculation of relative permeability", "[Virtual]")
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

        dictionary krDict("krModel(oil,water)");
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

        WHEN("correct() is called on a B-C relPerm model")
        {
            forAll(mesh.C(), ci)
            {
                waterPtr->alpha()[ci] = 0.2 + ci*(1-0.1-0.2)/(mesh.nCells()-1);
            }
            krModel->correct();
            THEN("Kr values must be consistent with model expression")
            {
                const std::vector<scalar> KrwManual = {
                    0, 0.012345679, 0.049382716, 0.111111111,
                    0.197530864, 0.308641975, 0.444444444, 0.604938272,
                    0.790123457, 1
                };
                
                const std::vector<scalar> KrnManual = {
                    0.9, 0.711111111, 0.544444444, 0.4, 0.277777778,
                    0.177777778, 0.1, 0.044444444, 0.011111111, 0
                };

                const DimensionedField<scalar, volMesh>& krw = krModel()[relPermModel<iRock, 2>::krName("water")];
                REQUIRE_THAT
                (
                    KrwManual, 
                    Catch::Matchers::Approx
                    (
                        std::vector<scalar> ( krw.begin(), krw.end() )
                    )
                );
            }
        }
    }
}

// ************************************************************************* //
