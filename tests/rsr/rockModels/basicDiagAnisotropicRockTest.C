#include "IsotropyTypes.H"
#include "autoPtr.H"
#include "catch.H"
#include "dimensionedScalarFwd.H"
#include "fvCFD.H"
#include "rock.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


using namespace Foam;

SCENARIO("Rock object creation for diagonalTensor-Permeability","[Virtual]")
{
    GIVEN("Valid mesh and rockProperties dictionary with porosity,"
            " permeability and compressibility values")
    {
        #include "createTestTimeAndMesh.H"

        word rockName = "rock";

        dictionary rockDict;
        rockDict.add("incompressible", true);
        rockDict.add<dimensionedScalar>
        (
            "porosity",
            dimensionedScalar(rockName+".porosity", dimless, 0.2)
        );
        rockDict.add<dimensionedVector>
        (
            "permeability",
            dimensioned<vector>
            (
                rockName+".permeability", dimArea,
                vector(1e-12, 2e-12, 5e-15)
            )
        );
        rockDict.add<dimensionedScalar>
        (
            "compressibility",
            dimensionedScalar
            (
                rockName+".compressibility", dimless/dimPressure, 1e-6
            )
        );

        dictionary rockProperties;
        rockProperties.add(rockName, rockDict);

        WHEN("Constructing the default rock object")
        {
            // Needs the presence of '0/water.U' dictionary
            auto rockPtr = rock<DiagAnisotropic>::New
            (
                rockName,
                mesh,
                rockProperties
            );

            THEN("Members must be initialized correctly")
            {
                CHECK(rockPtr->Cf()[0] == 1e-6);
                CHECK(rockPtr->porosity()[0] == 0.2);
                CHECK(rockPtr->K()[0] == vector(1e-12, 2e-12, 5e-15));
            }
            FatalError.throwExceptions();
        }

    }
}

// ************************************************************************* //
