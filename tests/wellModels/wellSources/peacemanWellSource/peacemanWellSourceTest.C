#include "IsotropyTypes.H"
#include "capPressModel.H"
#include "catch.H"
#include "autoPtr.H"
#include "fvCFD.H"
#include "wellSource.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


using namespace Foam;

SCENARIO("Peaceman's description of well source", "[Virtual]")
{
    GIVEN("Valid mesh, kr and pc models, and a valid sourceProperties object")
    {
        FatalError.dontThrowExceptions();
        #include "createTestTimeAndMesh.H"
        #include "createTestBlackoilPhase.H"
        #include "createTestIsoRock.H"
        #include "readGravitationalAcceleration.H"

        dictionary transportProperties;
        dictionary rockProperties;

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

        createTestBlackoilPhase(water, 1.0, 1e-3, multiPhase);
        createTestBlackoilPhase(oil, 1.0, 1e-5, multiPhase); 

        createTestIsoRock(rk, 1e-12, 0.2, 1e-6);

        dictionary krDict("krModel<water,oil>");
        krDict.add("type", "BrooksCorey");
        krDict.add<dimensionedScalar>
        (
            "water.alphaIrr",
            dimensionedScalar("water.alphaIrr", dimless, 0.2)
        );
        krDict.add<dimensionedScalar>
        (
            "oil.alphaRes",
            dimensionedScalar("oil.alphaRes", dimless, 0.1)
        );
        krDict.add<dimensionedScalar>
        (
            "water.m",
            dimensionedScalar("mater.m", dimless, 2)
        );
        krDict.add<dimensionedScalar>
        (
            "oil.m",
            dimensionedScalar("oil.m", dimless, 2)
        );
        krDict.add<dimensionedScalar>
        (
            "water.krMax",
            dimensionedScalar("water.krMax", dimless, 1.0)
        );
        krDict.add<dimensionedScalar>
        (
            "oil.krMax",
            dimensionedScalar("oil.krMax", dimless, 0.9)
        );
        
        transportProperties.add(word("krModel<water,oil>"), krDict);

        auto krModel = relPermModel<iRock, 2>::New
        (
            word("krModel<water,oil>"),
            transportProperties,
            rkPtr()
        );

        dictionary pcDict("pcModel<water,oil>");
        pcDict.add("type", "BrooksCorey");
        pcDict.add<dimensionedScalar>
        (
            "water.alpha.PcMin",
            dimensionedScalar("water.alpha.PcMin", dimless, 0.15)
        );
        pcDict.add<dimensionedScalar>
        (
            "water.alpha.PcMax",
            dimensionedScalar("water.alpha.PcMax", dimless, 0.85)
        );
        pcDict.add<dimensionedScalar>
        (
            "n",
            dimensionedScalar("n", dimless, 0.238)
        );
        pcDict.add<dimensionedScalar>
        (
            "pc0",
            dimensionedScalar("pc0", dimPressure, 30)
        );
        
        transportProperties.add(word("pcModel<water,oil>"), pcDict);

        auto pcModel = capPressModel<iRock, 2>::New
        (
            word("pcModel<water,oil>"),
            transportProperties,
            rkPtr()
        );

        dictionary wSrcDict("wellSourceDict");
        wSrcDict.add<word>("wellSourceType", "Peaceman");

        auto wSource = wellSource<iRock, 2>::New
        (
            "source", waterPtr(), wSrcDict, rkPtr()
        );

        dictionary srcPropsDict;
        srcPropsDict.add<dimensionedScalar>
        (
            "radius",
            dimensionedScalar("radius", dimLength, 0.3)
        );
        srcPropsDict.add<scalar>("skin", 2);
        srcPropsDict.add<word>("orientation", "vertical");
        srcPropsDict.add<word>("operationMode", "production");

        cellSet cSet(mesh, "cSet", 0);
        faceSet fSet(mesh, "fSet", 0);

        sourceProperties wellProps(mesh, srcPropsDict, cSet, fSet);

        WHEN("A set of cell IDs is passed to calculateCoeff*()")
        {
            forAll(mesh.C(), ci)
            {
                waterPtr->alpha()[ci] = 0.2 + ci*(1-0.1-0.2)/(mesh.nCells()-1);
            }
            scalarField coeff0(mesh.nCells(), 0);

            krModel->correct();
            pcModel->correct();

            labelList cellIDs(mesh.nCells(), 0);
            std::generate
            (
                cellIDs.begin(), cellIDs.end(),
                [n = 0] () mutable { return n++; }
            );
            wSource->calculateCoeff0(coeff0, wellProps, cellIDs);

            THEN("coeff0 values must be consistent with peaceman's expression")
            {
                std::vector<scalar> expectedCoeff0 = {
                    0, -2.1879e-11, -8.75159e-11, -1.96911e-10,
                    -3.50064e-10, -5.46975e-10, -7.87643e-10, -1.07207e-09,
                    -1.40025e-09, -1.7722e-09
                };
                REQUIRE_THAT
                (
                    expectedCoeff0, 
                    Catch::Matchers::Approx
                    (
                        std::vector<scalar> ( coeff0.begin(), coeff0.end() )
                    )
                );
            }
        }
    }
}

// ************************************************************************* //
