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

namespace Foam {

using matrix = std::vector<std::vector<scalar>>;

matrix lduMatrixToSparse(const fvMesh& mesh, fvScalarMatrix& mat)
{
    std::vector<std::vector<scalar>> fmat
    (
        mesh.nCells(), std::vector<scalar>(mesh.nCells(), 0.0)
    );
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const labelUList& lowerAddr = mat.lduAddr().lowerAddr();
    const labelUList& upperAddr = mat.lduAddr().upperAddr();
    for(int ii=0;ii<mesh.nCells();ii++)
    {
        for(int jj=0;jj<mesh.nCells();jj++)
        {
            if (jj == ii) fmat[ii][ii] = mat.diag()[ii];
        }
    }
    forAll(owner, facei)
    {
        label ii = lowerAddr[facei];
        label jj = upperAddr[facei];
        fmat[ii][jj] = mat.lower()[facei];
    }
    forAll(neighbour, facei)
    {
        label ii = lowerAddr[facei];
        label jj = upperAddr[facei];
        fmat[jj][ii] = mat.upper()[facei];
    }
    //Info << setprecision(3);
    //for(int ii=0;ii<mesh.nCells();ii++)
    //{
    //    for(int jj=0;jj<mesh.nCells();jj++)
    //    {
    //        Info << fmat[ii][jj] << setw(10);
    //    }
    //    Info << nl ;
    //}
    //Info << nl;
    //for (unsigned i = 0; i < mesh.nCells(); ++i) {
    //Info << mat.source()[i] << setw(10);
    //}
    Info << nl << endl;
    return fmat;
}
};


using namespace Foam;

SCENARIO("Imposed Phase-flowrate for a well", "[Virtual]")
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
        HashTable<fvScalarMatrix> matTable;
        matTable.insert("water", fvScalarMatrix(p, dimless));
        matTable.insert("oil", fvScalarMatrix(p, dimless));

        WHEN("Phase flowRate driveHandler is constructed and calls correct()")
        {
            forAll(mesh.C(), ci)
            {
                waterPtr->alpha()[ci] = 0.2 + ci*(1-0.1-0.2)/(mesh.nCells()-1);
                p[ci] = (10+2*ci)*6874.76; // start at 10psi, + 2psi per cell
            }

            dictionary driveDict("flowRate");
            driveDict.add<word>("phase", "water");
            driveDict.add<fileName>("file", "testData/water.rate.dat");

            autoPtr<driveHandler<iRock, 2>> dH = driveHandler<iRock, 2>::New
            (
                "prod.dH", driveDict, wSource(), wellProps, matTable
            );
            krModel->correct();
            pcModel->correct();
            dH->correct();
            matrix mat = lduMatrixToSparse(mesh, matTable["water"]);

            THEN("Well Matrix for a phase must be consistent with expected one")
            {
                // Construct the reference matrix
                matrix expectedMat
                (
                    mesh.nCells(), std::vector<scalar>(mesh.nCells(), 0)
                );
                std::vector<scalar> expectedMatSource ( mesh.nCells(), 0 );
                // qi = ai * pi + bi * BHP + ci
                // Get coefficients from well source
                scalarList ai(cSet.size()),bi(cSet.size()),ci(cSet.size());
                wSource->calculateCoeff0(ai, wellProps, cSet.toc());
                wSource->calculateCoeff1(bi, wellProps, cSet.toc());
                wSource->calculateCoeff2(ci, wellProps, cSet.toc());
                scalar biSum = sum(bi);
                forAll(bi, i)
                {
                    bi[i] /= biSum;
                }
                scalar qt = 1e-6; // target flowRate
                for (label i = 0; i < cSet.toc().size(); ++i) {
                    const label cell = cSet.toc()[i];

                    // diagonal and source
                    expectedMat[cell][cell] = ai[i]*(1-bi[i]);
                    expectedMatSource[cell] = ci[i] + bi[i]*qt - bi[i]*sum(ci);

                    // off-diagonal coeffs (Code for 1D meshes, well cells
                    // are assumed to be consecutive ...)
                    for (label j = 0; j < cSet.toc().size(); ++j) {
                        const label cellj = cSet.toc()[j];
                        if (cellj == cell+1 or cellj == cell-1)
                        {
                            expectedMat[cell][cellj] += -bi[i]*ai[j];
                        } else if (cellj != cell ) {
                            expectedMatSource[cell] += -bi[i]*ai[j]*p[cellj];
                        }
                    }
                }
                //Info << setprecision(3) ;
                //for (unsigned i = 0; i < expectedMat.size(); ++i) {
                //for (unsigned j = 0; j < expectedMat.size(); ++j)
                //    Info << expectedMat[i][j] << setw(10);
                //Info << nl;
                //}
                //for (unsigned i = 0; i < expectedMat.size(); ++i) {
                //Info << expectedMatSource[i] << setw(10);
                //}
                Info << nl;

                // Test equality of expected and calculated matrices
                REQUIRE_THAT
                (
                    expectedMatSource, 
                    Catch::Matchers::Approx
                    (
                        std::vector<scalar>
                        ( 
                            matTable["water"].source().begin(),
                            matTable["water"].source().end()
                        )
                    )
                );
                REQUIRE_THAT
                (
                    mat,
                    Catch::Matchers::Predicate<matrix>
                    (
                        [expectedMat](matrix const& m) -> bool 
                        {
                            bool eq = true;
                            for (size_t i = 0; i < m.size(); ++i) {
                            for (size_t j = 0; j < m.size(); ++j) {
                                if (Approx(m[i][j]) != expectedMat[i][j])
                                {
                                    eq = false;
                                    break;
                                }
                            }
                            }
                            return eq;
                        },
                        "Matrix coeffs should match expected ones."
                    )
                );
            }
        }
    }
}

// ************************************************************************* //
