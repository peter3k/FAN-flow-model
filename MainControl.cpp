/*
  licence: GNU GPL v3 licence

  Full-Scale Airway Network (FAN) flow model: MainControl.cpp  
  Copyright (C) 2019  Minsuok Kim

*/


/*
* 
* Goal: To simlulate the dynamic lung ventilation in a full-scale conducting
         airways

* How it works: The FAN flow model takes the structure of airways and lobes 
                and lung tissue information extracted from thoracic CT images 
                as an input data set. The pulmonary function test data could 
                be applied to assume the initial and boundary condition. 
                The FAN model computes the dynamic air pressure, flow rate 
                and gas density to provide insights on lung ventilation.  

* Refs.: Kim M, Dognay O, Matin T, et al. CT-based airway flow model to assess 
         ventilation in COPD: A pilot study, Radiology. (accepted)
   
         Kim M, Bordas R, Vos W, et al. Dynamic flow characteristics in normal 
         and asthmatic lungs. Int J Numer Method Biomed Eng. 2015;31(12):e02730.
 
** Note: Currently available laboratory version of the code is for research 
         purpose only. The current file set does not include pre-processed 
         input files and functions for an operation. We are preparing a user-
         friendly public version of the software which includes examples and 
         all necessary functions. 
*
*/

#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <unistd.h> 

#include "Structure.h"
#include "BoundarySetup.h"
#include "PressureFlow.h"
#include "Refreshing.h"
#include "ScalarTransport.h"
#include "DataIO.h"
#include "DataAnalysis.h"

Structure Geometry;
AcinarDynamics ADynamics;
BoundarySetup BDsetup;
PressureFlow Pflow;
ScalarTransport ScalTrans;
DataIO DTIO;
Refreshing RefreshData;
DataAnalysis DTAnalysis;

using namespace std;

int main()
{
    int fnode,vfnode,felem,fgen,sacini,facini,fbr,tioid;
    double maxX,minX,maxY,minY,maxZ,minZ;
    double highac,lowac;
    double V0,V0TLC,sumAcf,TotalAcinarVolume;
    double Iit00,Iit00old;
    double TLC,RV,Seg_LungVolume,TemporalLungVolume;
    double lambdazero, scaling=1.0;

    string location1D,location3D,results1D,Geometry1D,InputConditions;
    string filename,filename1Ddata;
    string inletname="MOUTH";
    string headname="BD";

    char CurrentPath[256];
    GetCurrentDir(CurrentPath,255);
    strcat(CurrentPath, "/");
    location1D=CurrentPath;

    Geometry1D=location1D+"Geometry/";
    InputConditions=location1D+"Input_Conditions/";
    results1D=location1D+"Results/";

    BDsetup.LungVolume(InputConditions,TLC,RV,lambdazero,Seg_LungVolume,TemporalLungVolume);

    Geometry.nodenumber(Geometry1D,fnode);
    std::vector<std::vector<double> > ncoord;
    std::vector<std::vector<double> > tempncoord;
    std::vector<int> ioid(fnode);

    std::vector<int> couplebdnode;
    Geometry.nodedata(Geometry1D,ncoord,tempncoord,ioid,couplebdnode,tioid);

    std::vector<int> npt(fnode),nge(fnode),nbr(fnode);
    std::vector<std::vector<int> > econn;
    std::vector<std::vector <int> > nconn(fnode, std::vector<int> (2));
    std::vector<int> tbr4g;
    std::vector<int> accumbranch;

    Geometry.elementdata(Geometry1D,fgen,felem,fbr,npt,nge,nbr,tbr4g,accumbranch,nconn,econn);

    vfnode=fbr+1;
    int *vfnodept=&vfnode;

    std::vector<double> elength;
    Geometry.elementlength(econn,ncoord,elength);
    std::vector<double> angle;
    Geometry.gravity(econn,ncoord,angle);

    std::vector<std::vector<int> > br;
    std::vector<std::vector<int> > LobeID;
    std::vector<std::vector<int> > IDcoupledbds;
    std::vector<std::vector<int> > HorsfieldNo;
    std::vector<std::vector<double> > brlength;
    std::vector<std::vector<double> > brradius;
    std::vector<std::vector<double> > RoutTLC;
    std::vector<std::vector<double> > WallthickTLC;
    std::vector<std::vector<double> > brthickness;
    std::vector<std::vector<double> > r0;
    std::vector<std::vector<double> > brang;
    std::vector<std::vector<double> > Raw;
    std::vector<std::vector<double> > Raw_sum;
    std::vector<std::vector<double> > Rac;
    std::vector<std::vector<double> > AcinarElastance;
    std::vector<std::vector<double> > AirwayCompliance;

    std::vector<std::vector<double> > lambdaold;
    std::vector<std::vector<double> > lambda;
    std::vector<std::vector<std::vector<int> > > brconn;
    std::vector<std::vector<std::vector<int> > > br2node;

    std::vector<double> P1D(tioid);
    std::vector<double> P1Dold(tioid);
    std::vector<double> P3D(tioid);
    std::vector<double> P1D3Ddiff(tioid);
    std::vector<double> P1D3Ddiffold(tioid);
    std::vector<double> Q1D3D(tioid);
    std::vector<double> Kaccel(tioid);

    Geometry.vector2Dint(fgen,tbr4g,br);
    Geometry.vector2Dint(fgen,tbr4g,LobeID);
    Geometry.vector2Dint(fgen,tbr4g,IDcoupledbds);
    Geometry.vector2Dint(fgen,tbr4g,HorsfieldNo);
    Geometry.vector2Ddb(fgen,tbr4g,brlength);
    Geometry.vector2Ddb(fgen,tbr4g,brradius);
    Geometry.vector2Ddb(fgen,tbr4g,RoutTLC);
    Geometry.vector2Ddb(fgen,tbr4g,WallthickTLC);
    Geometry.vector2Ddb(fgen,tbr4g,brthickness);
    Geometry.vector2Ddb(fgen,tbr4g,r0);
    Geometry.vector2Ddb(fgen,tbr4g,brang);
    Geometry.vector2Ddb(fgen,tbr4g,Raw);
    Geometry.vector2Ddb(fgen,tbr4g,Raw_sum);
    Geometry.vector2Ddb(fgen,tbr4g,Rac);
    Geometry.vector2Ddb(fgen,tbr4g,AcinarElastance);
    Geometry.vector2Ddb(fgen,tbr4g,AirwayCompliance);

    Geometry.vector2Ddb(fgen,tbr4g,lambdaold);
    Geometry.vector2Ddb(fgen,tbr4g,lambda);
    Geometry.vector3D(fgen,tbr4g,3,brconn);
    Geometry.vector3D(fgen,tbr4g,2,br2node);

    Geometry.branchlength(elength,nge,nbr,econn,scaling,brlength);
    Geometry.maxmincoord(ncoord,maxX,minX,maxY,minY,maxZ,minZ);
    Geometry.circuit(npt,nge,nbr,econn,br,angle,brang,br2node);
    Geometry.branchconn(nge,npt,nbr,econn,nconn,brconn);
    Geometry.aciniheight(fgen,tbr4g,br,br2node,ncoord,highac,lowac);

    ADynamics.AcinivolumeIni(fgen,tbr4g,br,sacini,facini,RV,TLC,Seg_LungVolume,V0,V0TLC);

    std::vector<std::vector<double> > Iitold;

    std::vector<double> NodalDyeDensity(fnode);
    std::vector<double> InspiredAcinarFlow(fnode);
    std::vector<double> ExpiredAcinarFlow(fnode);
    std::vector<std::vector<double> > DyeDensity;
    std::vector<std::vector<double> > DyeDensityold;
    std::vector<std::vector<double> > DyeinAcini;

    std::vector<double> NodalAcinarVolume (fnode);
    std::vector<double> NodalAcinarDensity (fnode);
    std::vector<std::vector<double> > AcinarVolumeold2;
    std::vector<std::vector<double> > AcinarVolumeold;
    std::vector<std::vector<double> > AcinarVolume;
    std::vector<std::vector<double> > AcinarDensityold;
    std::vector<std::vector<double> > AirwayVolume;

    Pflow.InIvector2Ddb(fgen,tbr4g,Iitold);
    Pflow.InIvector2Ddb(fgen,tbr4g,DyeDensity);
    Pflow.InIvector2Ddb(fgen,tbr4g,DyeDensityold);
    Pflow.InIvector2Ddb(fgen,tbr4g,DyeinAcini);
    Pflow.InIvector2Ddb(fgen,tbr4g,AcinarVolumeold2);
    Pflow.InIvector2Ddb(fgen,tbr4g,AcinarVolumeold);
    Pflow.InIvector2Ddb(fgen,tbr4g,AcinarVolume);
    Pflow.InIvector2Ddb(fgen,tbr4g,AcinarDensityold);
    Pflow.InIvector2Ddb(fgen,tbr4g,AirwayVolume);

    Geometry.Horsfield(fgen,tbr4g,nge,nbr,br,brconn,HorsfieldNo);

    std::vector<std::vector <int> > Lobenode;
    BDsetup.LobalNodes(InputConditions,Lobenode);
    Geometry.LobeIndex(fnode,fgen,tbr4g,nge,nbr,br,br2node,brconn,Lobenode,LobeID);

    std::vector<double> Nodalradius(fnode);
    Geometry.fieldata(Geometry1D,fnode,fgen,tbr4g,scaling,brradius,r0,Nodalradius,
                      nge,nbr,econn,brthickness);

    ofstream Strdata;
    filename=results1D+"output-tbr4g.txt";
    Strdata.open(filename.c_str());
    DTIO.WriteStructure(fgen,fbr,facini,tbr4g,Strdata);

    std::vector<std::vector <double> > Ptold;
    std::vector<std::vector <double> > Ptold2;
    std::vector<std::vector <double> > Pt;
    std::vector<std::vector <double> > dPt;
    std::vector<std::vector <double> > Pt_PeakInhal;
    std::vector<std::vector <double> > Iit;
    std::vector<std::vector <double> > Q;
    std::vector<std::vector <double> > Qold;
    std::vector<std::vector <double> > Qold2;
    std::vector<std::vector <double> > InspiredIit;

    Pflow.InPvector2Ddb(fgen,tbr4g,Ptold);
    Pflow.InPvector2Ddb(fgen,tbr4g,Ptold2);
    Pflow.InPvector2Ddb(fgen,tbr4g,Pt);
    Pflow.InPvector2Ddb(fgen,tbr4g,dPt);
    Pflow.InPvector2Ddb(fgen,tbr4g,Pt_PeakInhal);
    Pflow.InIvector2Ddb(fgen,tbr4g,Iit);
    Pflow.InIvector2Ddb(fgen,tbr4g,Q);
    Pflow.InIvector2Ddb(fgen,tbr4g,Qold);
    Pflow.InIvector2Ddb(fgen,tbr4g,Qold2);
    Pflow.InIvector2Ddb(fgen,tbr4g,InspiredIit);

    int n,sn,timestep;
    int svnum,scount,vtkscount;
    int breathcount;
    int maxbreathcount;
    int outputstartbreath,outputendbreath;
    int vtkstartbreath;
    int deepbreathing;
    int coupleiter;

    double maxIcrit=0.01;
    double deepbreathingtime,deepbreathingfactor;
    double time,dtime,sdtdtime;
    double tpm,tmp;
    double cyc,omega,wcon1,wcon2;
    double Cmean,InspiredFlowVolume;
    double Pplscale,Ppl,Pplold,newdPpl;

    BDsetup.SetBreathing(InputConditions,n,sn,maxbreathcount,outputstartbreath,outputendbreath,
    		         vtkstartbreath,deepbreathing,deepbreathingfactor,cyc,omega,wcon1,wcon2,
                         Pplscale,dtime,sdtdtime);

    BDsetup.TissueDensity(InputConditions,fnode,facini,npt,nge,nbr,econn,RV,Seg_LungVolume,
                          AcinarVolume,lambda);

    int injstartbreath,injendbreath;
    double Mixfactor;
    double InGasConc;
    BDsetup.SetScalarInjection(InputConditions,injstartbreath,injendbreath,InGasConc,Mixfactor);

    ofstream PQdata;
    filename=results1D+"output-nodedata.txt";
    PQdata.open(filename.c_str());

    ofstream PSampledata;
    filename=results1D+"output-Psample.txt";
    PSampledata.open(filename.c_str());

    ofstream QSampledata;
    filename=results1D+"output-Qsample.txt";
    QSampledata.open(filename.c_str());
    filename=results1D+"test.txt";


    while (breathcount<maxbreathcount)
    {
        sumAcf=0.0;
        timestep++;

        if (time==0)
        {
            BDsetup.PleuralPre(time-dtime,cyc,breathcount,
            		       deepbreathing,deepbreathingtime,deepbreathingfactor,
            		       omega,wcon1,wcon2,Pplold,Ppl,newdPpl);
        }
       
        BDsetup.PleuralPre(time,cyc,breathcount,
        		   deepbreathing,deepbreathingtime,deepbreathingfactor,
        		   omega,wcon1,wcon2,Pplold,Ppl,newdPpl);

        Pflow.MatrixComputing(*vfnodept,fgen,facini,breathcount,maxbreathcount,
         		      time,dtime,cyc,newdPpl,
                              V0,lambdazero,Cmean,highac,lowac,
                              nge,nbr,tbr4g,accumbranch,Lobenode,LobeID,
                              ncoord,brradius,brlength,brang,
                              br,brconn,br2node,Rac,Raw,lambdaold,lambda,
                              AcinarElastance,AirwayCompliance,
                              Qold2,Qold,Q,Iit,Ptold,Pt,dPt,
			      emphysemalobe,emphystartbreath,emphyendbreath,emphyfactor,
			      EBVlobe,EBVneighbourlobe,EBVstartbreath,EBVendbreath,EBVfactor);

        DTAnalysis.ComputedParameters(fgen,breathcount,outputstartbreath,outputendbreath,dtime,
                                      tbr4g,Q,Qold,Qold2,InspiredIit,Ptold,Pt_PeakInhal);

        Pflow.BreathCount(breathcount,InspiredFlowVolume,time,dtime,tmp,tpm,Iit00old,Iit00,Iit);

        RefreshData.OneDPressureWallupdate(time,dtime,fgen,felem,sumAcf,TotalAcinarVolume,V0,V0TLC,
                                           scaling,tbr4g,
                                           br,brconn,br2node,brradius,Nodalradius,brlength,brthickness,
                                           AirwayCompliance,lambdaold,lambda,lambdazero,
                                           AcinarVolumeold2,AcinarVolumeold,AcinarVolume,
					   RoutTLC,WallthickTLC,ioid,couplebdnode,
					   Pt,dPt,Ptold,Ptold2,Iit,npt,nge,nbr,ncoord,econn);

        RefreshData.InspiredAcinarFlowVolume(fnode,dtime,npt,nge,nbr,Qold,Q,
       		                             Iit,InspiredAcinarFlow,ExpiredAcinarFlow;

        ScalTrans.ScalarFlow(fnode,fgen,breathcount,injstartbreath,injendbreath,InGasConc,Mixfactor,
        		     time,dtime,V0,nge,nbr,npt,tbr4g,br,ncoord,
			     br2node,brconn,brradius,brlength,lambda,Q,Qold,
			     NodalDyeDensity,DyeDensity,DyeDensityold,DyeinAcini,
			     NodalAcinarVolume,AcinarVolume,
			     NodalAcinarDensity,AcinarDensityold,AirwayVolume);

        DTIO.WriteVTK(results1D,fnode,felem,fgen,breathcount,
       		      svnum,vtkscount,vtkstartbreath,time,sdtdtime,
       		      nge,nbr,npt,tbr4g,tempncoord,br2node,brconn,Nodalradius,
                      NodalAcinarVolume,LobeID,econn,NodalDyeDensity,NodalAcinarDensity,Pt,Iit,
		      InspiredAcinarFlow,ExpiredAcinarFlow,HorsfieldNo);

        DTIO.Write1D(breathcount,maxbreathcount,outputstartbreath,outputendbreath,
       		     svnum,scount,fgen,time,dtime,sdtdtime,Ppl,
                     tbr4g,br,br2node,brconn,brradius,Pt,Iit,Raw,AcinarVolume,
                     PQdata,PSampledata,QSampledata,InspiredFlowVolume);

        time+=dtime;

        if ((time>100)&&(breathcount>=maxbreathcount+1))
        {
            break;
        }
    }

    Strdata.close();
    PQdata.close();
    PSampledata.close();
    QSampledata.close();

    cout<<"Mission completed"<<endl;

    return 0;
}
