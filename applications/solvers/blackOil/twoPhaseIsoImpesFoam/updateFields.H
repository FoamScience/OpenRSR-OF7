// Update phase properties
canPhasePtr->correct();
nonCanPhasePtr->correct();

// Update eqn. terms
alphaStorage = porosity * Bc;
pStorage = porosity * (Bc - Bn);
