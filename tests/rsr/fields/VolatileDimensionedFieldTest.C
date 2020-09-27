#include <vector>
#include "catch.H"
#include "error.H"
#include "fvCFD.H"
#include "VolatileDimensionedField.H"
#include "DimensionedTensorField.H"
#include "scalar.H"
#include "singleCellFvMesh.H"
#include "UniformDimensionedField.H"
#include "vector.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

#define TEST_DIMENSIONED_FIELD_BINARY_OP(t1,df1,df2,op,opName)                 \
        WHEN("Two volatile fields are " #opName)                               \
        {                                                                      \
            auto vdf3 = v##df1 op v##df2;                                      \
            auto df3 = df1 op df2;                                             \
            THEN("Operation is consistent with standard Dimensioned Fields")   \
            {                                                                  \
                REQUIRE                                                        \
                (                                                              \
                    std::vector<t1>(vdf3->begin(), vdf3->end())                \
                    == std::vector<t1>(df3->begin(), df3->end())               \
                );                                                             \
            }                                                                  \
        }


#define TEST_DIMENSIONED_FIELD_BINARY_FUNCTION(t1,vdf1,vdf2,FuncName)          \
        WHEN(#FuncName "(Two volatile fields)")                                \
        {                                                                      \
            auto vdf3 = FuncName(vdf1, vdf2);                                  \
            auto df3 = FuncName(df1, df2);                                     \
            THEN("Operation is consistent with standard Dimensioned Fields")   \
            {                                                                  \
                REQUIRE                                                        \
                (                                                              \
                    std::vector<t1>(vdf3->begin(), vdf3->end())                \
                    == std::vector<t1>(df3->begin(), df3->end())               \
                );                                                             \
            }                                                                  \
        }

#define TEST_VOLATILE_FIELD_BINARY_OP(t1,vdf1,vdf2,op,opName)                  \
        WHEN("Two volatile fields are " #opName)                               \
        {                                                                      \
            auto vdf3 = vdf1 op vdf2;                                          \
            THEN("Operation result match expected vector")                     \
            {                                                                  \
                REQUIRE_THAT                                                   \
                (                                                              \
                    std::vector<t1>(vdf3->begin(), vdf3->end()),               \
                    Catch::Matchers::Approx( expected##opName ).margin(1e-6)   \
                );                                                             \
            }                                                                  \
        }

#define TEST_TMP_DIMENSIONED_FIELD_BINARY_OP(t1,vdf1,vdf2,op1,op2,opName)      \
        WHEN("Two temporary volatile fields are " #opName)                     \
        {                                                                      \
            auto vdf3 = (vdf1 op1 vdf2) op2 (vdf1 op2 vdf2);                   \
            THEN("Operation is consistent with standard Dimensioned Fields")   \
            {                                                                  \
                REQUIRE_THAT                                                   \
                (                                                              \
                    std::vector<t1>(vdf3->begin(), vdf3->end()),               \
                    Catch::Matchers::Approx( expectedTmp##opName ).margin(1e-6)\
                );                                                             \
            }                                                                  \
        }


SCENARIO("Compatibility of volatile fields with standard DimensionedFields"
        " when they have the same size" )
{
    GIVEN("A valid mesh")
    {
        #include "createTestTimeAndMesh.H"
        DimensionedField<scalar, volMesh> df1
        (
            IOobject
            (
                "df1",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("df", 0.0)
        );
        DimensionedField<scalar, volMesh> df2
        (
            IOobject
            (
                "df2",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("df", 0.0)
        );
        VolatileDimensionedField<scalar, volMesh> vdf1
        (
            IOobject
            (
                "vdf1",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("vdf", 0.0)
        );
        VolatileDimensionedField<scalar, volMesh> vdf2
        (
            IOobject
            (
                "vdf2",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("vdf", 0.0)
        );

        for (label i = 0; i < mesh.nCells(); ++i)
        {
            df1[i] = 0.1 + 15*i/mesh.nCells();
            vdf1[i] = df1[i];
            df2[i] = -1.1 + 20*i/mesh.nCells();
            vdf2[i] = df2[i];
        }

        TEST_DIMENSIONED_FIELD_BINARY_OP(scalar,df1,df2,+,added)
        TEST_DIMENSIONED_FIELD_BINARY_OP(scalar,df1,df2,+,added)
        TEST_DIMENSIONED_FIELD_BINARY_OP(scalar,df1,df2,-,subtracted)
        TEST_DIMENSIONED_FIELD_BINARY_OP(scalar,df1,df2,*,outer-multiplied)

        TEST_DIMENSIONED_FIELD_BINARY_FUNCTION(scalar,df1,df2,pow)
    }
}

#undef TEST_DIMENSIONED_FIELD_BINARY_OP
#undef TEST_DIMENSIONED_FIELD_BINARY_FUNCTION

SCENARIO("Operations involving Single-Cell volatile fields")
{
    GIVEN("A valid mesh and its single-cell version")
    {
        #include "createTestTimeAndMesh.H"
        FatalError.dontThrowExceptions();
        VolatileDimensionedField<scalar, volMesh> vdf1
        (
            IOobject
            (
                "vdf1",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("vdf", 0.0)
        );
        singleCellFvMesh singleCell
        (
            IOobject
            (
                "oneCell",
                mesh.polyMesh::instance(),
                mesh.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        VolatileDimensionedField<scalar, volMesh> vdf2
        (
            IOobject
            (
                "vdf2",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            singleCell,
            dimensionedScalar("vdf", 0.0)
        );

        for (label i = 0; i < mesh.nCells(); ++i)
        {
            vdf1[i] = 0.1 + 15*i/mesh.nCells();
        }
        vdf2[0] = 1.6;

        std::vector<scalar> expectedAdded
        {
            1.7, 2.7, 4.7, 5.7, 7.7, 8.7, 10.7, 11.7, 13.7, 14.7
        };
        std::vector<scalar> expectedSubtracted
        {
            -1.5, -0.5, 1.5, 2.5, 4.5, 5.5, 7.5, 8.5, 10.5, 11.5
        };
        std::vector<scalar> expectedSubtracted2
        {
            1.5, 0.5, -1.5, -2.5, -4.5, -5.5, -7.5, -8.5, -10.5, -11.5
        };
        std::vector<scalar> expectedMultiplied
        {
            0.16, 1.76, 4.96, 6.56, 9.76, 11.36, 14.56, 16.16, 19.36, 20.96
        };
        std::vector<scalar> expectedCrossMultiplied
        {
            0.16, 1.76, 4.96, 6.56, 9.76, 11.36, 14.56, 16.16, 19.36, 20.96
        };

        scalarList tmpM(expectedAdded.size());
        for (size_t i = 0; i < expectedAdded.size(); ++i) {
            tmpM[i] = expectedAdded[i]*expectedSubtracted[i];
        }

        TEST_VOLATILE_FIELD_BINARY_OP(scalar,vdf1,vdf2,+,Added)
        TEST_VOLATILE_FIELD_BINARY_OP(scalar,vdf2,vdf1,+,Added)
        TEST_VOLATILE_FIELD_BINARY_OP(scalar,vdf1,vdf2,-,Subtracted)
        TEST_VOLATILE_FIELD_BINARY_OP(scalar,vdf2,vdf1,-,Subtracted2)
        TEST_VOLATILE_FIELD_BINARY_OP(scalar,vdf1,vdf2,*,Multiplied)

        std::vector<scalar> expectedTmpAdded
        {
            3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2
        };
        TEST_TMP_DIMENSIONED_FIELD_BINARY_OP(scalar,vdf1,vdf2,+,-,Added)

        tensor initialTensor { 0.1, 0.2, 0.3, 0.5, 0.7, 0.1, 0.5, 0.6, 0.4};

        VolatileDimensionedField<tensor, volMesh> tvdf1
        (
            IOobject
            (
                "vdf2",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<tensor>("vdf", tensor::zero)
        );
        forAll(tvdf1, i)
        {
            tvdf1[i].xx() = 1 + 0.2*i;
            tvdf1[i].yy() = 1 + 0.3*i;
            tvdf1[i].zz() = 2 - 0.1*i;
        }
        VolatileDimensionedField<tensor, volMesh> tvdf2
        (
            IOobject
            (
                "vdf2",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            singleCell,
            dimensioned<tensor>("vdf", initialTensor)
        );
        DimensionedField<tensor, volMesh> tdf1("tdf1", tvdf1);
        DimensionedField<tensor, volMesh> tdf2
        (
            IOobject
            (
                "df2",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<tensor>("vdf", initialTensor)
        );
        WHEN("Two volatile tensor fields are scalar-multiplied" )
        {
            auto tvdf3 = tvdf1 && tvdf2;
            auto tdf3 = tdf1 && tdf2;
            THEN("Operation result must match equivalent DimensionedField")
            {
                //forAll(tvdf3(), ci)
                //{
                //    Info << tvdf3()[ci] << nl << expectedTensor[ci] << nl;
                //}
                //std::vector<tensor> vt(tvdf3->begin(), tvdf3->end());
                REQUIRE
                (
                    std::vector<tensor>(tvdf3->begin(), tvdf3->end())
                    == std::vector<tensor>(tdf3->begin(), tdf3->end())
                );
            }
        }
        WHEN("Two volatile tensor fields are inner-multiplied" )
        {
            auto tvdf3 = tvdf1 & tvdf2;
            auto tdf3 = tdf1 & tdf2;
            THEN("Operation result must match equivalent DimensionedField")
            {
                //forAll(tvdf3(), ci)
                //{
                //    Info << tvdf3()[ci] << nl << expectedTensor[ci] << nl;
                //}
                //std::vector<tensor> vt(tvdf3->begin(), tvdf3->end());
                REQUIRE
                (
                    std::vector<tensor>(tvdf3->begin(), tvdf3->end())
                    == std::vector<tensor>(tdf3->begin(), tdf3->end())
                );
            }
        }
        vector initialVector { 1.1, 5.2, -3.3};
        VolatileDimensionedField<vector, volMesh> vvdf1
        (
            IOobject
            (
                "vdf2",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<vector>("vdf", vector::zero)
        );
        forAll(vvdf1, i)
        {
            vvdf1[i].x() = 1 + 0.2*i;
            vvdf1[i].y() = 1 + 0.3*i;
            vvdf1[i].z() = 2 - 0.1*i;
        }
        VolatileDimensionedField<vector, volMesh> vvdf2
        (
            IOobject
            (
                "vdf2",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            singleCell,
            dimensioned<vector>("vdf", initialVector)
        );
        DimensionedField<vector, volMesh> vecdf1("vdf1", vvdf1);
        DimensionedField<vector, volMesh> vecdf2
        (
            IOobject
            (
                "df2",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<vector>("vdf", initialVector)
        );

        WHEN("Two volatile vector fields are outer-multiplied" )
        {
            auto vvdf3 = vvdf1 * vvdf2;
            auto vdf3 = vecdf1 * vecdf2;
            THEN("Operation result must match equivalent DimensionedField")
            {
                REQUIRE
                (
                    std::vector<tensor>(vvdf3->begin(), vvdf3->end())
                    == std::vector<tensor>(vdf3->begin(), vdf3->end())
                );
            }
        }
        WHEN("Two volatile vector fields are cross-multiplied" )
        {
            auto vvdf3 = vvdf1 ^ vvdf2;
            auto vdf3 = vecdf1 ^ vecdf2;
            THEN("Operation result must match equivalent DimensionedField")
            {
                REQUIRE
                (
                    std::vector<vector>(vvdf3->begin(), vvdf3->end())
                    == std::vector<vector>(vdf3->begin(), vdf3->end())
                );
            }
        }
    }
}

#undef TEST_VOLATILE_FIELD_BINARY_OP

// ************************************************************************* //
