#include "catch.H"
#include "fvCFD.H"
#include "FVFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

SCENARIO("FVF Model Selection with Default configuration", "[Virtual]")
{
    //Foam::FatalError.dontThrowExceptions();

    GIVEN("Pressure field and RunTime-Selectable implementation of 'FVFModel'")
    {
        #include "createTestTimeAndMesh.H"

        dictionary dict;
        dict.add(word("FVFModel"), word(""));

        WHEN("Passing invalid model identifier to the selector")
        {
            // Try to instantiate a model from an "InvalidModelName" class
            dict.set("FVFModel", "InvalidModelName");

            THEN("An exception is raised")
            {
                REQUIRE_THROWS
                (
                    FVFModel<Compressible>::New("invalidModel", dict, mesh)
                );
            }
        }

        WHEN("Calling correct() on a valid FVF model object")
        {
            // Try to instantiate a model from an "InvalidModelName" class
            dict.set("FVFModel", "childFVFModel");
            auto fvf = FVFModel<Compressible>::New("fvf", dict, mesh);
            fvf->correct();
            THEN("1|FVF BCs must reflect a zeroGradient situation")
            {
                const volScalarField& rfvf = fvf->rFVF();
                forAll(mesh.boundary(), pi)
                {
                    const word& patchType = mesh.boundary()[pi].type();
                    if (patchType != "empty")
                    forAll(mesh.boundary()[pi], fi)
                    {
                        const label& bcell = 
                            mesh.boundaryMesh()[pi].faceCells()[fi];
                        REQUIRE
                        (
                            rfvf.internalField()[bcell] 
                            == rfvf.boundaryField()[pi][fi]
                        );
                    }
                }
            }
            THEN("d(1|FVF)dP BCs must reflect a zeroGradient situation")
            {
                const volScalarField& drfvfdp = fvf->rFVF();
                forAll(mesh.boundary(), pi)
                {
                    const word& patchType = mesh.boundary()[pi].type();
                    if (patchType != "empty")
                    forAll(mesh.boundary()[pi], fi)
                    {
                        const label& bcell = 
                            mesh.boundaryMesh()[pi].faceCells()[fi];
                        REQUIRE
                        (
                            drfvfdp.internalField()[bcell] 
                            == drfvfdp.boundaryField()[pi][fi]
                        );
                    }
                }
            }
        }
    }
}

// ************************************************************************* //
