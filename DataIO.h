/*
  licence: GNU GPL v3 licence

  Full-Scale Airway Network (FAN) flow model: DataIO.h  
  Copyright (C) 2019  Minsuok Kim

*/

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "Structure.h"

#ifndef DATAIO_H_
#define DATAIO_H_
using namespace std;

class DataIO : public Structure {

    public:

    void WriteStructure (int& fgen,
                         int& fbr,
			 int& facini,
			 std::vector<int>& tbr4g,
			 std::ofstream& Strdata);

    void WriteVTK(string& results1D,
		  int& fnode,
		  int& felem,
		  int& fgen,
		  int& breathcount,
		  int& svnum,
		  int& vtkscount,
	          int& vtkstartbreath,
	          double& time,
		  double& sdtdt,
		  std::vector<int>& nge,
		  std::vector<int>& nbr,
	          std::vector<int>& npt,
		  std::vector<int>& tbr4g,
	          std::vector<std::vector<double> >&tempncoord,
		  std::vector<std::vector<std::vector<int> > >& br2node,
		  std::vector<std::vector<std::vector<int> > >& brconn,
		  std::vector<double>& Nodalradius,
		  std::vector<double>& NodalAcinarVolume,
		  std::vector<std::vector<int> >& LobeID,
		  std::vector<std::vector<int> >& econn,
	          std::vector<double>& NodalDyeDensity,
		  std::vector<double>& NodalAcinarDensity,
		  std::vector<std::vector<double> > Pt,
		  std::vector<std::vector<double> > Iit,
		  std::vector<double> InspiredAcinarFlow,
		  std::vector<double> ExpiredAcinarFlow,
		  std::vector<std::vector<int> >& HorsfieldNo);

    void Interpolation(int& point0,
    		       int& point1,
    		       int& point2,
                       double& value0,
    		       double& value1,
    		       double& value2,
                       std::vector<std::vector<double> >&ncoord);


    void Write1D(int& breathcount,
		 int& maxbreathcount,
		 int& outputstartbreath,
		 int& outputendbreath,
		 int& svnum,
		 int& scount,
		 int& fgen,
		 double& time,
		 double& dtime,
		 double& sdtdt,
		 double& Ppl,
	         std::vector<int>& tbr4g,
	         std::vector<std::vector<int> >& br,
		 std::vector<std::vector<std::vector<int> > >& br2node,
		 std::vector<std::vector<std::vector<int> > >& brconn,
		 std::vector<std::vector<double> >& brradius,
		 std::vector<std::vector<double> >& Pt,
		 std::vector<std::vector<double> >& Iit,
		 std::vector<std::vector<double> >& Raw,
		 std::vector<std::vector<double> >& AcinarVolume,
 		 std::ofstream& PQdata,
		 std::ofstream& PSampledata,
		 std::ofstream& QSampledata,
		 double& InspiredFlowVolume);


    void DistalBranchLobe(int& gen,
	                  int& acn,
	                  std::vector<std::vector<std::vector<int> > >& brconn,
		          std::vector<std::vector<int> >& LobeID,
			  int& DistLobeMax,
			  int& DistLobeMin);

};

#endif /* DATAIO_H_ */
