// Update phase properties
canPhasePtr->correct();
nonCanPhasePtr->correct();

// Update eqn. terms
alphaTemp = porosity * Bc;
pTemp = porosity * (Bc - Bn);
