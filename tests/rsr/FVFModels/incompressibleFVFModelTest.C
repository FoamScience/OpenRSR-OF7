#include "catch.H"
#include "fvCFD.H"
#include "FVFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

SCENARIO("Incompressible FVF Construct")
{
    GIVEN("A Dimensioned scalar and a dummy mesh")
    {
        // The interface of FVFModel requires a mesh
        // but does not use it in incompressible state
        #include "createTestTimeAndMesh.H"

        // The phase dictionary
        dictionary dict("incompressible");
        dict.add<word>("FVFModel", "incompressible");
        dict.add<dimensionedScalar>
        (
            "rFVF",
            dimensionedScalar("rFVF", dimless, 1.0)
        );

        WHEN("Quering for rFVF and drFVFdP values")
        {
            CHECK_NOTHROW(FVFModel::New("fvf", dict, mesh));
            auto fvf = FVFModel::New("fvf", dict, mesh);
            // Perform model calculations
            fvf->correct();

            THEN("No changes to rFVF should be made")
            {
                CHECK( fvf->rFVF().size() == 1);
                REQUIRE( fvf->rFVF()[0] == 1.0 );
                REQUIRE( fvf->drFVFdP()[0] == 0.0 );
            }
        }
    }
}

// ************************************************************************* //
