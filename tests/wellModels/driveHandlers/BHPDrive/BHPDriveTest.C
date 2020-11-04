#include "IsotropyTypes.H"
#include "capPressModel.H"
#include "catch.H"
#include "autoPtr.H"
#include "fvCFD.H"
#include "wellSource.H"
#include "volFieldsFwd.H"
#include "driveHandler.H"
#include "IOmanip.H"
#include "relPermModel.H"
#include "capPressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

SCENARIO("Imposed BHP for a well", "[Virtual]")
{
    GIVEN("Valid mesh, kr and pc models, sourceProperties object and "
          "a HashTable of phase matrices")
    {
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

        HashTable<autoPtr<wellSource<iRock,2>>> sources;
        sources.insert
        (
            "water",
            wellSource<iRock, 2>::New
            (
                "sourceW", waterPtr(), wSrcDict, rkPtr()
            )
        );
        sources.insert
        (
            "oil",
            wellSource<iRock, 2>::New
            (
                "sourceO", oilPtr(), wSrcDict, rkPtr()
            )
        );

        dictionary srcPropsDict;
        srcPropsDict.add<dimensionedScalar>
        (
            "radius",
            dimensionedScalar("radius", dimLength, 0.3)
        );
        srcPropsDict.add<scalar>("skin", 2);
        srcPropsDict.add<word>("orientation", "vertical");
        srcPropsDict.add<word>("operationMode", "injection");
        srcPropsDict.add<word>("injectedPhase", "water");

        // The cell and face sets
        cellSet cSet(mesh, "cSet", 8);
        cSet.insert(labelList{1,2,3,4,5,6,7,8});
        cSet.write();
        faceSet fSet(mesh, "fSet", 7);
        cellToFace fSrc(mesh, "cSet", cellToFace::BOTH);
        fSrc.applyToSet(topoSetSource::ADD, fSet);

        sourceProperties wellProps(mesh, srcPropsDict, cSet, fSet);
        HashPtrTable<fvScalarMatrix> matTable;
        matTable.insert("water", new fvScalarMatrix(p, dimless));
        matTable.insert("oil", new fvScalarMatrix(p, dimless));

        scalar BHP = 1.563e6;
        WHEN("Phase flowRate driveHandler is constructed and calls correct()")
        {
            forAll(mesh.C(), ci)
            {
                waterPtr->alpha()[ci] = 0.2 + ci*(1-0.1-0.2)/(mesh.nCells()-1);
                p[ci] = (10+2*ci)*6874.76;
            }

            dictionary driveDict("BHP");
            driveDict.add<fileName>("file", "testData/BHP.dat");


            autoPtr<driveHandler<iRock, 2>> dH = driveHandler<iRock, 2>::New
            (
                "prod.dH", driveDict, sources, wellProps, matTable
            );
            krModel->correct();
            pcModel->correct();
            dH->correct();

            THEN("Well Matrix for a phase must be consistent with expected one")
            {
                // Construct the reference matrix
                std::vector<scalar> expectedDiag ( mesh.nCells(), 0 );
                std::vector<scalar> expectedSource ( mesh.nCells(), 0 );
                // qi = ai * pi + bi * BHP + ci
                // Get coefficients from well source
                scalarList ai(cSet.size()),bi(cSet.size()),ci(cSet.size());
                sources["water"]->calculateCoeff0(ai, wellProps, cSet.toc());
                sources["water"]->calculateCoeff1(bi, wellProps, cSet.toc());
                sources["water"]->calculateCoeff2(ci, wellProps, cSet.toc());
                scalar apSum = 0;
                forAll(ai, ci)
                {
                    apSum += ai[ci]*p[cSet.toc()[ci]];
                }
                for (label i = 0; i < cSet.toc().size(); ++i) {
                    const label cell = cSet.toc()[i];
                    // diagonal and source
                    expectedDiag[cell] = ai[i];
                    expectedSource[cell] = bi[i]*BHP + ci[i];
                }
                REQUIRE_THAT
                (
                    expectedDiag, 
                    Catch::Matchers::Approx
                    (
                        std::vector<scalar>
                        ( 
                            (*matTable["water"]).diag().begin(),
                            (*matTable["water"]).diag().end()
                        )
                    )
                );
                REQUIRE_THAT
                (
                    expectedSource, 
                    Catch::Matchers::Approx
                    (
                        std::vector<scalar>
                        ( 
                            (*matTable["water"]).source().begin(),
                            (*matTable["water"]).source().end()
                        )
                    )
                );
            }
        }
    }
}

// ************************************************************************* //
