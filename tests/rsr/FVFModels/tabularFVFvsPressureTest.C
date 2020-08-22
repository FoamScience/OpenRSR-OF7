#include "catch.H"
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
        dictionary dict;
        dict.add(word("FVFModel"), word("tabularFVFvsPressure"));

        dictionary fvfData;
        fvfData.add<fileName>("file", "testData/FVF.dat");
        fvfData.add("interpolationType", "linear"); // Default
        dict.add("FVFData", fvfData);

        WHEN("Supplied valid P values")
        {
            autoPtr<FVFModel<Compressible>> fvf =
                FVFModel<Compressible>::New("fvf", dict, mesh);
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

            std::vector<scalar> rFVFCalculated(pCheckList.size());
            for(size_t ci = 0; ci < rFVFCalculated.size(); ++ci)
            {
                rFVFCalculated[ci] = fvf->rFVF().internalField()[ci];
            }
            std::vector<scalar> drFVFdPCalculated(pCheckList.size());
            for(size_t ci = 0; ci < drFVFdPCalculated.size(); ++ci)
            {
                drFVFdPCalculated[ci] = fvf->drFVFdP().internalField()[ci];
            }

            THEN("Linear interpolation must be fairly accurate")
            {
                REQUIRE_THAT
                (
                    rFVFCalculated,
                    Catch::Matchers::Approx(rFVFManual)
                );
                REQUIRE_THAT
                (
                    drFVFdPCalculated,
                    Catch::Matchers::Approx(drFVFdPManual)
                );
            }
        }
    }
}

// ************************************************************************* //
