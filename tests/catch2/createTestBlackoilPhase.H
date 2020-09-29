#define createTestBlackoilPhase(phaseName, rhoSc, mu, context)\
    word phaseName = #phaseName;\
    dictionary phaseName##Dict;\
    phaseName##Dict.add<word>("phaseType", "blackoil");\
    phaseName##Dict.add<word>("incompressible", "true");\
    phaseName##Dict.add<word>("FVFModel", "tabularFVFvsPressure");\
    dictionary phaseName##FvfData; \
    phaseName##FvfData.add<fileName>("file", "testData/FVF.dat"); \
    phaseName##FvfData.add("interpolationType", "linear"); \
    phaseName##Dict.add("FVFData", phaseName##FvfData); \
    phaseName##Dict.add<dimensionedScalar>\
    (\
        "rhoSc",\
        dimensionedScalar("rhoSc", dimDensity, rhoSc)\
    );\
    phaseName##Dict.add<dimensionedScalar>\
    (\
        "mu0",\
        dimensionedScalar("mu0", dimDensity, mu)\
    );\
    transportProperties.add(phaseName, phaseName##Dict);\
    auto phaseName##Ptr = phase::New\
    (\
        phaseName,\
        mesh,\
        transportProperties,\
        mixtureType::context\
    )

