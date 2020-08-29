#include "UniformityTypes.H"
#include "autoPtr.H"
#include "catch.H"
#include "fvCFD.H"
#include "error.H"
#include "phase.H"
#include "samplePhase.H"
#include "IncompressiblePhase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


using namespace Foam;
using namespace phases;

SCENARIO("Phase Object creation", "[Virtual]")
{
    GIVEN("Valid mesh and transportProperties")
    {
        #include "createTestTimeAndMesh.H"

        word phaseName = "water";
        dictionary transportProperties;
        dictionary waterDict;
        waterDict.add<word>("FVFModel", "incompressible");
        waterDict.add<scalar>("rhoSc", 1.0);
        transportProperties.add(phaseName, waterDict);

        WHEN("Constructing phase in a singlePhase setup")
        {
            // Needs the presence of '0/water.U' dictionary
            samplePhase water
            (
                phaseName,
                mesh,
                transportProperties,
                mixtureType::singlePhase
            );

            THEN("Requesting a ref to alpha field errors out")
            {
                REQUIRE_THROWS(water.alpha());
            }
        }

        WHEN("Constructing phase in a multiPhase setup")
        {
            samplePhase water
            (
                phaseName,
                mesh,
                transportProperties,
                mixtureType::multiPhase
            );

            THEN("Requesting a ref to alpha field errors out")
            {
                REQUIRE_NOTHROW(water.alpha());
            }
        }

        WHEN("Operating on the velocity of a stack-allocated phase object")
        {
            samplePhase water(phaseName, mesh, transportProperties);

            surfaceScalarField& phi = water.phi();
            volVectorField& U = water.U();

            // Imitate a linear U profile
            forAll(U, ci)
            {
                U[ci] = vector(mesh.nCells() - ci, 0, 0);
            }
            U.correctBoundaryConditions();

            // Needs to set div(water.phi, water.U) in system/fvSchemes
            water.correct();

            THEN("Velocity must be constructed properly from phi "
                    "using Weller's approximation")
            {
                // Particularly test if fvc::reconstruct has changed!
                // fvc::reconstruct may be dissipative near some BCs
                surfaceScalarField newPhi( linearInterpolate(U) & mesh.Sf() );
                REQUIRE( phi.internalField()==newPhi.internalField() );
            }
        }

    }
}

// ************************************************************************* //
