#include "UniformityTypes.H"
#include "autoPtr.H"
#include "catch.H"
#include "fvCFD.H"
#include "samplePhase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


using namespace Foam;
using namespace phases;

SCENARIO("Phase Object creation", "[Virtual]")
{
    GIVEN("Valid mesh and transportProperties")
    {
        #include "createTestTimeAndMesh.H"

        word phaseName = "phase";
        dictionary transportProperties;
        dictionary phaseDict(phaseName);
        phaseDict.add<word>("FVFModel", "incompressible");
        phaseDict.add<dimensionedScalar>
        (
            "rFVF",
            dimensionedScalar("rFVF", dimless, 1.0)
        );
        phaseDict.add<dimensionedScalar>
        (
            "rhoSc",
            dimensionedScalar("rhoSc", dimDensity, 1.0)
        );
        phaseDict.add<dimensionedScalar>
        (
            "rho0",
            dimensionedScalar("rho0", dimDensity, 1.0)
        );
        phaseDict.add<dimensionedScalar>
        (
            "mu0",
            dimensionedScalar("mu0", dimViscosity*dimDensity, 1.0)
        );
        transportProperties.add(phaseName, phaseDict);

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

            THEN("Requesting a ref to alpha field doesn't throw")
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
