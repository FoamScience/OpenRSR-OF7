#include "catch.H"
#include "dimensionedScalarFwd.H"
#include "error.H"
#include "fvCFD.H"
#include "FVFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

SCENARIO("FVF Model Selection with Default configuration", "[Virtual]")
{
    GIVEN("Pressure field and RunTime-Selectable implementation of 'FVFModel'")
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

        WHEN("Constructing a child FVF model in incompressible mode")
        {
            dict.set("FVFModel", "childFVFModel");
            auto fvf = FVFModel::New("fvf", dict, mesh);
            THEN("The selected FVF mesh is a single-cell one")
            {
                CHECK(fvf->rFVF().size() == 1);
                CHECK(fvf->drFVFdP().size() == 1);

                // REQUIRE THAT 0/water.rFVF is preferred over dimensionedScalar
                REQUIRE(fvf->rFVF()[0] == 1.5);
                REQUIRE(fvf->drFVFdP()[0] == 0.0);
            }
        }
        WHEN("Constructing a child FVF model in compressible mode")
        {
            dict.set("FVFModel", "childFVFModel");
            dict.set
            (
                "rFVF",
                dimensionedScalar("rFVF", dimless, -1.0)
            );
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
