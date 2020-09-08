#include "DimensionedField.H"
#include "IOobject.H"
#include "UList.H"
#include "dimensionedScalar.H"
#include "singleCellFvMesh.H"
#include "scalarField.H"
#include "dictionary.H"

using namespace Foam;


int main(int argc, char* argv[]) {
    
    #include "createTestTimeAndMesh.H"

    autoPtr<singleCellFvMesh> oneMesh_
    (
        new singleCellFvMesh
        (
            IOobject
            (
                "single.mesh",
                mesh.polyMesh::instance(),
                runTime,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    dictionary dct;
    //dct.add<scalar>("ss", 100);
    dimensionedScalar ss("ss", dimless, dct.lookupOrAddDefault("ss", -1));


    DimensionedField<scalar,volMesh> sF
    (
            IOobject
            (
                "switchabale",
                mesh.polyMesh::instance(),
                runTime,
                ss.value() == -1
                    ? IOobject::MUST_READ
                    : IOobject::READ_IF_PRESENT,
                //IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            ((ss.value() == -1) ? mesh : oneMesh_() ),
            dimensionedScalar("", dimless, 2.0)
    );
    
    DimensionedField<scalar, volMesh> rF ("copy", sF+15.0);
    forAll(sF, i) { Info << rF[i] << endl;}

    return 0;
}
