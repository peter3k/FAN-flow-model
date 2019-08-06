/*
 * BoundarySetup.h
 *
 *  Created on: 16 Aug 2018
 *      Author: Minsuok Kim

This file is part of Full-scale Airway Network (FAN) Flow Model-V0.0.9
	
FAN is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
	
FAN is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
License for more details. The offer of FAN under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.
	
You should have received a copy of the GNU Lesser General Public License
along with FAN. If not, see <http://www.gnu.org/licenses/>.

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
