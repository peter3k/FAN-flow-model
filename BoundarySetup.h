/*
  licence: GNU GPL v3 licence

  Full-Scale Airway Network (FAN) flow model: BoundarySetup.h  
  Copyright (C) 2019  Minsuok Kim

*/

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "Structure.h"

#ifndef BOUNDARYSETUP_H_
#define BOUNDARYSETUP_H_

using namespace std;

class BoundarySetup {

    public:

    void SetBreathing(string InputConditions,
                      int& n,
    		      int& sn,
                      int& maxbreathcount,
		      int& outputstartbreath,
	              int& outputendbreath,
                      int& vtkstartbreath,
                      int& deepbreathing,
                      double& deepbreathingfactor,
                      int& solvertype,
                      int& maxiter,
                      double& maxIcrit,
                      double& cyc,
                      double& omega,
                      double& wcon1,
                      double& wcon2,
                      double& Pplscale,
                      double& dtime,
                      double& sdtdtime);


    void PleuralPre(double time,
                    double& cyc,
    		    int& breathcount,
		    int& deepbreathing,
		    double& deepbreathingtime,
		    double& deepbreathingfactor,
    		    double& omega,
		    double& wcon1,
		    double& wcon2,
		    double& Pplold,
		    double& Ppl,
		    double& newdPpl);


    void SetScalarInjection(string InputConditions,
    		            int& injstartbreath,
    	                    int& injendbreath,
			    double& InGesConc,
			    double& Mixfactor);


    void LungVolume(string InputConditions,
    		    double& TLC,
		    double& RV,
		    double& lambdazero,
		    double& Seg_LungVolume,
		    double& TemporalLungVolume);


    void LobalNodes(string InputConditions,
  		    std::vector<std::vector<int> >& Lobenode);

};

#endif /* BOUNDARYSETUP_H_ */
