#include "UniformityTypes.H"
#include "autoPtr.H"
#include "catch.H"
#include "fvCFD.H"
#include "error.H"
#include "phase.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


using namespace Foam;

SCENARIO("Blackoil Phase Object creation for Uniform Viscosity", "[Virtual]")
{
    GIVEN("Valid mesh, FVFModel configuration and transportProperties")
    {
        #include "createTestTimeAndMesh.H"

        word phaseName = "water";
        dictionary transportProperties;
        dictionary waterDict;
        waterDict.add<word>("phaseType", "blackoil");
        waterDict.add<word>("FVFModel", "tabularFVFvsPressure");
        waterDict.add<scalar>("rhoSc", 10.0);

        dictionary fvfData;
        fvfData.add<fileName>("file", "testData/FVF.dat");
        fvfData.add("interpolationType", "linear"); // Default
        waterDict.add("FVFData", fvfData);
        transportProperties.add(phaseName, waterDict);

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

        forAll(p.internalField(), ci)
        {
            p[ci] = 7.50e4 + 1.7e7 * ci/mesh.nCells();
        }

        WHEN("Constructing blackoil phase in a singlePhase setup with"
             " linear FVF interpolation")
        {
            FatalError.dontThrowExceptions();
            // Needs the presence of '0/water.U' dictionary
            auto waterPtr = phase<Compressible,UniformMu>::New
            (
                phaseName,
                mesh,
                transportProperties,
                mixtureType::singlePhase
            );

            THEN("Calling correct() updates density field")
            {
                waterPtr->correct();
                auto rho = waterPtr->rho().internalField();

                // Expected rho vals
                std::vector<scalar> expectedRho
                { 9.46416, 8.60152, 8.30824, 8.09611, 7.92945,
                  7.78276, 7.70249, 7.7278, 7.74854, 7.77165};
                
                REQUIRE_THAT
                (
                    expectedRho, 
                    Catch::Matchers::Approx
                    ( std::vector<scalar>(rho.begin(), rho.end()) )
                );
                
            }
            FatalError.throwExceptions();
        }

    }
}

// ************************************************************************* //
