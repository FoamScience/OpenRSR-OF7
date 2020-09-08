#include "catch.H"
#include "error.H"
#include "fvCFD.H"
#include "FVFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

SCENARIO("Tabulated FVF Model with before-and-after bubble point data")
{
    GIVEN("Pressure field and FVF data file")
    {
        #include "createTestTimeAndMesh.H"

        // The pressure field
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

        // The phase dictionary
        dictionary dict("tabular");
        dict.add<word>("FVFModel", "tabularFVFvsPressure");

        dictionary fvfData;
        fvfData.add<fileName>("file", "testData/FVF.dat");
        fvfData.add("interpolationType", "linear"); // Default
        dict.add("FVFData", fvfData);

        WHEN("Attempting construction in incompressible context")
        {
            dictionary incDict(dict);
            incDict.add<dimensionedScalar>
            (
                "rFVF",
                dimensionedScalar("rFVF", dimless, 1.0)
            );
            THEN("Construction errors out.")
            {
                REQUIRE_THROWS(FVFModel::New("fvf", incDict, mesh));
            }
        }
        WHEN("Supplied valid P values")
        {
            auto fvf = FVFModel::New("fvf", dict, mesh);
            std::vector<scalar> pCheckList = {
                7.50e4, 1.00e5, 3.00e5, 4.00e5, 8.00e5, 1.50e6, 7.00e6,
                1.00e7, 1.80e7, 2.30e7
            };
            std::vector<scalar> rFVFManual = {
                0.94641579, 0.944130711, 0.925850035, 0.916709697, 0.885988783,
                0.86582205, 0.791810943, 0.769912425, 0.7798852, 0.786176065
            };
            std::vector<scalar> drFVFdPManual = {
                -9.13394e-8, -8.94851e-8, -7.46509e-8, -6.72338e-8, -3.50993e-8,
                -1.93656e-8, -8.70729e-9, 1.22407e-9, 1.43534e-9, 1.46173e-9
            };

            // Initialize pressure
            forAll(p.internalField(), ci)
            {
                p[ci] = pCheckList[ci];
            }

            // Perform model calculations
            fvf->correct();

            auto& rfvf = fvf->rFVF();
            auto& drfvf = fvf->drFVFdP();

            THEN("Linear interpolation must be fairly accurate")
            {
                REQUIRE_THAT
                (
                    rFVFManual,
                    Catch::Matchers::Approx
                    (
                        std::vector<scalar>(rfvf.begin(), rfvf.end())
                    )
                );
                REQUIRE_THAT
                (
                    drFVFdPManual,
                    Catch::Matchers::Approx
                    (
                        std::vector<scalar>(drfvf.begin(), drfvf.end())
                    )
                );
            }
        }
    }
}

// ************************************************************************* //
