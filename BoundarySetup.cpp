/*
 * BoundarySetup.cpp
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




#include "BoundarySetup.h"

void BoundarySetup::SetBreathing(string InputConditions,
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
				 double& sdtdtime)
{

    int writeinterval=0;
    double sdt=0.0;
    double maxPpl,minPpl,halfdPpl;

    Structure st;   
	string filename,line,rdata;

    ifstream inputdatafile;
    filename=InputConditions+"Pressure_Breathing_Solver_DataSaving.txt";
    inputdatafile.open(filename.c_str());

    if(!inputdatafile)
    {
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }

    for (int i=1;i<5;i++)
    {
        getline (inputdatafile,line);
    }


    getline (inputdatafile,line);
    st.findnum(line,5,rdata);
    maxPpl=atof(rdata.c_str());

    getline (inputdatafile,line);
    st.findnum(line,5,rdata);
    minPpl=atof(rdata.c_str());

    getline (inputdatafile,line);
    st.findnum(line,5,rdata);
    Pplscale=atof(rdata.c_str());

    halfdPpl=Pplscale*0.5*(maxPpl-minPpl);
    maxPpl=0.5*(maxPpl+minPpl)+halfdPpl;
    minPpl=maxPpl-2.0*halfdPpl;

    wcon1=-0.5*(minPpl+maxPpl);
    wcon2=0.5*(maxPpl-minPpl);


    getline (inputdatafile,line);
    getline (inputdatafile,line);
    st.findnum(line,4,rdata);
    cyc=atof(rdata.c_str());

    getline (inputdatafile,line);
    st.findnum(line,5,rdata);
    maxbreathcount=atoi(rdata.c_str());


    getline (inputdatafile,line);
    getline (inputdatafile,line);
    st.findnum(line,3,rdata);
    solvertype=atof(rdata.c_str());


    getline (inputdatafile,line);
    getline (inputdatafile,line);
    maxiter=0;
    maxIcrit=0.0;

   
    getline (inputdatafile,line);
    st.findnum(line,5,rdata);
    dtime=atof(rdata.c_str());

    getline (inputdatafile,line);
    st.findnum(line,6,rdata);
    writeinterval=atoi(rdata.c_str());

    getline (inputdatafile,line);
    getline (inputdatafile,line);
    st.findnum(line,7,rdata);
    outputstartbreath=atoi(rdata.c_str());

    getline (inputdatafile,line);
    st.findnum(line,7,rdata);
    outputendbreath=atoi(rdata.c_str());

    getline (inputdatafile,line);
    getline (inputdatafile,line);
    st.findnum(line,8,rdata);
    vtkstartbreath=atoi(rdata.c_str());

    getline (inputdatafile,line);
    getline (inputdatafile,line);
    st.findnum(line,7,rdata);
    deepbreathing=atoi(rdata.c_str());

    getline (inputdatafile,line);
    st.findnum(line,6,rdata);
    deepbreathingfactor=atof(rdata.c_str());

    cout<<"Solver type number = "<<solvertype<<endl;
    if (solvertype==2)
    {
    	cout<<"Maximum iteration = "<<maxiter
		    <<", Maximum stop criteria = "<<maxIcrit<<endl;
    }

    cout<<"Deepbreathing = "<<deepbreathing
    	<<", Deepbreathingfactor = "<<deepbreathingfactor
		<<endl;

    getline (inputdatafile,line);
    cout<<line<<endl<<endl;

    inputdatafile.close();

    omega=2.0*M_PI/cyc; 
    n=(30*cyc)/dtime+0.5;
    sdt=dtime*writeinterval;
    sdtdtime=sdt/dtime;
    sn=(30*cyc)/sdt+0.5;
}

void BoundarySetup::PleuralPre(double time,
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
			       double& newdPpl)
{

    if ((breathcount==deepbreathing)&&(deepbreathingtime==0)&&(Ppl==(-(wcon1+wcon2*std::sin(0.5*M_PI)))))
    {
        omega=2.0*M_PI/(cyc*deepbreathingfactor);
        deepbreathingtime=time;
    }

    Ppl=-(wcon1+wcon2*std::sin(omega*(time-deepbreathingtime)+0.5*M_PI));

    newdPpl=Ppl-Pplold;
    Pplold=Ppl;
}


void BoundarySetup::SetScalarInjection(string InputConditions,
		                       int& injstartbreath,
		                       int& injendbreath,
			               double& InGasConc,
			               double& Mixfactor)
{
    Structure st;
    string filename,line,rdata;

    ifstream scalardatafile;
    filename=InputConditions+"Scalar_Injection.txt";
    scalardatafile.open(filename.c_str());

    if(!scalardatafile)
    {
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }

    
    for (int i=1;i<5;i++)
    {
        getline (scalardatafile,line);
    }

    getline (scalardatafile,line);
    st.findnum(line,5,rdata);
    injstartbreath=atoi(rdata.c_str());

    getline (scalardatafile,line);
    st.findnum(line,5,rdata);
    injendbreath=atoi(rdata.c_str());

    getline (scalardatafile,line);
    st.findnum(line,4,rdata);
    InGasConc=atof(rdata.c_str());

    getline (scalardatafile,line);
    st.findnum(line,4,rdata);
    Mixfactor=atof(rdata.c_str());
}


void BoundarySetup::LungVolume(string InputConditions,
		               double& TLC,
			       double& RV,
			       double& lambdazero,
			       double& Seg_LungVolume,
			       double& TemporalLungVolume)
{

    double FRC,RVTLC,FRCRV;

    Structure st;
    string filename,line,rdata;

    ifstream volumedatafile;
    filename=InputConditions+"CaseInfo_LungVolume.txt";
    volumedatafile.open(filename.c_str());


    if(!volumedatafile) 
    {
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }

    
    for (int i=1;i<5;i++)
    {
        getline (volumedatafile,line);
    }

    
    getline (volumedatafile,line);
    cout<<line<<endl;
    getline (volumedatafile,line);
    cout<<line<<endl;
    getline (volumedatafile,line);
    cout<<line<<endl;
    getline (volumedatafile,line);
    cout<<line<<endl;

    getline (volumedatafile,line);
    st.findnum(line,5,rdata);
    TLC=atof(rdata.c_str());

    getline (volumedatafile,line);
    st.findnum(line,4,rdata);
    RV=atof(rdata.c_str());

    getline (volumedatafile,line);
    st.findnum(line,3,rdata);
    RVTLC=atof(rdata.c_str());

    getline (volumedatafile,line);
    st.findnum(line,5,rdata);
    FRC=atof(rdata.c_str());

    FRCRV=FRC/RV;

    getline (volumedatafile,line);

    getline (volumedatafile,line);
    st.findnum(line,5,rdata);
    Seg_LungVolume=atof(rdata.c_str());
    TemporalLungVolume=Seg_LungVolume;          
    cout<<"Segmented Lung Volume: "<<Seg_LungVolume<<"(L)"<<endl;
    getline (volumedatafile,line);
    getline (volumedatafile,line);
    getline (volumedatafile,line);
    getline (volumedatafile,line);
    getline (volumedatafile,line);

    getline (volumedatafile,line);
    cout<<line<<endl;

    lambdazero=powf((Seg_LungVolume/RV),0.333333333);
    volumedatafile.close();
}


void BoundarySetup::LobalNodes(string InputConditions,
		               std::vector<std::vector<int> >& Lobenode)
{

    string filename,line,rnumnode,rlobenode;
    std::vector<int> LobalBranchStartingNode;

    Structure st;

    ifstream lobalnodedata;
    filename=InputConditions+"LobalBranchStartingNodes.txt";
    lobalnodedata.open(filename.c_str());

    if(!lobalnodedata)
    {
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }


    for (int i=1;i<5;i++)
    {
        getline (lobalnodedata,line);
    }


    for (int i=0;i<5;i++)
    {
        getline (lobalnodedata,line);
        st.findnum(line,6,rnumnode);
        int numnode=atoi(rnumnode.c_str());

        for (int j=1;j<=numnode;j++)
        {
            st.findnum(line,6+j,rlobenode);
            LobalBranchStartingNode.push_back(atoi(rlobenode.c_str()));
        }
        Lobenode.push_back(LobalBranchStartingNode);
        LobalBranchStartingNode.clear();
    }


    for (int i=0;i<5;i++)
    {
    	int maxlobenodes=Lobenode[i].size();
    	for(int j=0;j<maxlobenodes;j++)
    	{
    		cout<<"Lobenode["<<i<<"]["<<j<<"]="<<Lobenode[i][j]<<endl;
    	}
    }

    getline (lobalnodedata,line);
    cout<<line<<endl;

    lobalnodedata.close();
}

