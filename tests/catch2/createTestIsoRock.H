#define createTestIsoRock(rockName, K, phi, comp)\
    word rockName = #rockName; \
    dictionary rockName##Dict; \
    rockName##Dict.add("incompressible", true); \
    rockName##Dict.add<dimensionedScalar> \
    ( \
        "porosity", \
        dimensionedScalar(rockName+".porosity", dimless, phi)  \
    ); \
    rockName##Dict.add<dimensionedScalar> \
    ( \
        "permeability", \
        dimensionedScalar \
        ( \
            rockName+".permeability", dimArea, K \
        ) \
    ); \
    rockName##Dict.add<dimensionedScalar> \
    ( \
        "compressibility", \
        dimensionedScalar \
        ( \
            rockName+".compressibility", dimless/dimPressure, comp \
        ) \
    ); \
    rockProperties.add(rockName, rockName##Dict); \
    auto rockName##Ptr = rock<Isotropic>::New(rockName, mesh, rockProperties)

