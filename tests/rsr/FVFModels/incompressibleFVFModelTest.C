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
        dictionary dict;
        dict.add(word("FVFModel"), word("incompressible"));

        WHEN("Quering for rFVF and drFVFdP value")
        {
            CHECK_NOTHROW(incompressibleFVFModel::New("fvf", dict, mesh));
            auto fvf = incompressibleFVFModel::New("fvf", dict, mesh);
            // Perform model calculations
            fvf->correct();

            THEN("Linear interpolation must be fairly accurate")
            {
                REQUIRE ( fvf->rFVF().value() == 1.0 );
                REQUIRE ( fvf->drFVFdP().value() == 0.0 );
            }
        }
    }
}

// ************************************************************************* //
