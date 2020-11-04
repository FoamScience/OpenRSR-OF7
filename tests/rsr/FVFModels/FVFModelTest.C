#include "catch.H"
#include "dimensionedScalarFwd.H"
#include "error.H"
#include "fvCFD.H"
#include "FVFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

SCENARIO("FVF Model Selection with Default configuration", "[Virtual]")
{
    GIVEN("Valid mesh and RunTime-Selectable implementation of 'FVFModel'")
    {
        #include "createTestTimeAndMesh.H"

        dictionary dict("water");
        dict.add<dimensionedScalar>
        (
            "rFVF",
            dimensionedScalar("rFVF", dimless, 1.0)
        );
        dict.add<word>("FVFModel", "");

        WHEN("Passing invalid model identifier to the selector")
        {
            // Try to instantiate a model from an "InvalidModelName" class
            dict.set("FVFModel", "InvalidModelName");

            THEN("An exception is raised")
            {
                REQUIRE_THROWS
                (
                    FVFModel::New("invalidModel", dict, mesh)
                );
            }
        }

        WHEN("Constructing a child FVF model in compressible mode")
        {
            dict.set("FVFModel", "childFVFModel");
            dict.set("incompressible", false);

            auto fvf = FVFModel::New("fvf", dict, mesh);
            THEN("The selected FVF mesh is the full one")
            {
                CHECK(fvf->rFVF().size() == mesh.nCells());
                CHECK(fvf->drFVFdP().size() == mesh.nCells());

                REQUIRE
                (
                    std::vector<scalar>(fvf->rFVF().begin(), fvf->rFVF().end())
                    == std::vector<scalar>(mesh.nCells(), 1.5)
                );
                REQUIRE
                (
                    std::vector<scalar>(fvf->drFVFdP().begin(),fvf->drFVFdP().end())
                    == std::vector<scalar>(mesh.nCells(), 0.0)
                );
            }
        }
    }
}

// ************************************************************************* //
