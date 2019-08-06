/*
 * DataAnalysis.h
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
#include <sstream>
#include <iomanip>

#ifndef DATAANALYSIS_H_
#define DATAANALYSIS_H_
using namespace std;

class DataAnalysis {

    public:

	void LungResistance(int& fgen,
                            int& fbr,
		            int& facini,
		            int& breathcount,
			    int& outputstartbreath,
			    int& outputendbreath,
	                    int& svnum,
			    int& scount,
			    double& sdtdt,
			    double& time,
		            double& Ppl,
	                    double& newdPpl,
	                    std::vector<int>& tbr4g,
			    std::vector<std::vector<int> >& br,
		            std::vector<std::vector<std::vector<int> > >& brconn,
                            std::vector<std::vector<double> >& Pt,
	                    std::vector<std::vector<double> >& Iit,
	                    std::vector<std::vector<double> >& Raw,
	                    std::vector<std::vector<double> >& Raw_sum,
	                    std::vector<std::vector<double> >& AcinarVolume,
			    std::ofstream& LungResi,
			    std::ofstream& GenRes);


	void ComputedParameters(int& fgen,
			        int& breathcount,
			        int& outputstartbreath,
			        int& outputendbreath,
				double& dtime,
	                        std::vector<int>& tbr4g,
				std::vector<std::vector<double> >& Q,
				std::vector<std::vector<double> >& Qold,
				std::vector<std::vector<double> >& Qold2,
       			        std::vector<std::vector<double> >& InspiredIit,
				std::vector<std::vector<double> >& Ptold,
				std::vector<std::vector<double> >& Pt_PeakInhal);
};

#endif /* DATAANALYSIS_H_ */
