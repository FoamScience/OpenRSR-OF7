
const dictionary& solverDict = 
    mesh.solutionDict().subDict("coupledSolvers").subDict(solverName);

scalar maxDp = solverDict.lookupOrDefault<scalar>("initPMax", 1e2);
scalar maxDSc = solverDict.lookupOrDefault<scalar>("initScMax", 1e-2);

scalar dpLim = readScalar(runTime.controlDict().lookup("dpMax")); 
scalar dScLim = readScalar(runTime.controlDict().lookup("d("+Sc.name()+")Max")); 

p.storeOldTime();
Sc.storeOldTime();
