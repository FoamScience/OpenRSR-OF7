// Create Sample Time object
Foam::Time runTime
(
    Foam::Time::controlDictName,
    "./testData",
    "",
    "system",
    "constant",
    false
);

// Create Sample Mesh object
Foam::fvMesh mesh
(
    Foam::IOobject
    (
        Foam::fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);
