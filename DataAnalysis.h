/*
  licence: GNU GPL v3 licence

  Full-Scale Airway Network (FAN) flow model: DataAnalysis.h  
  Copyright (C) 2019  Minsuok Kim

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
