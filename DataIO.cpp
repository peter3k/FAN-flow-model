/*
  licence: GNU GPL v3 licence

  Full-Scale Airway Network (FAN) flow model: DataIO.cpp
  Copyright (C) 2019  Minsuok Kim

*/

#include "DataIO.h"

void DataIO::WriteStructure(int& fgen,
		            int& fbr,
			    int& facini,
			    std::vector<int>& tbr4g,
		            std::ofstream& Strdata)
{
    Strdata<<fgen<<endl;
    Strdata<<fbr<<endl;
    Strdata<<facini<<endl;

    for (int count=0;count<fgen;count++)
    {
        Strdata<<tbr4g[count]<<endl;
    }

}

void DataIO::WriteVTK(string& results1D,
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
                      std::vector<std::vector<int> >& HorsfieldNo)
{
	int cline;
	double scale=1;

    if (((breathcount>=vtkstartbreath))&&(svnum==sdtdt))
    {
    	ostringstream filename;
    	filename<<results1D+"/VTK_files/output-vtk_";
    	filename<<setfill('0')<<setw(4)<<vtkscount;
    	filename<<".vtk";
    	ofstream vtkdata;
    	vtkdata.open(filename.str().c_str());

        vtkdata<<"# vtk DataFile Version 3.0"<<endl;
        vtkdata<<"vtk output"<<endl;
        vtkdata<<"ASCII"<<endl;
        vtkdata<<"DATASET POLYDATA"<<endl;
        vtkdata<<"FIELD FieldData 1"<<endl;
        vtkdata<<"TIME 1 1 double"<<endl;
	vtkdata<<time<<endl<<endl;


        std::vector<double> Ptnode(fnode);
        std::vector<double> Iitnode(fnode);
        std::vector<int> LobeIDnode(fnode);
        std::vector<int> Horsfieldnumber(fnode);

        for (int gen=fgen-1;gen>=0;gen--)
        {
            for (int acn=0;acn<tbr4g[gen];acn++)
            {
            	int corrnode0=br2node[gen][acn][0];
            	int corrnode1=br2node[gen][acn][1];
            	Ptnode[corrnode1]=Iit[gen][acn];
            	Iitnode[corrnode1]=Iit[gen][acn];
            	Horsfieldnumber[corrnode1]=HorsfieldNo[gen][acn];

            	if (gen>0)
            	{
            		LobeIDnode[corrnode1]=LobeID[gen][acn];
            	}

            }
        }
        Iitnode[0]=Iit[0][0];
        Horsfieldnumber[0]=HorsfieldNo[0][0];

        int beginnode,endnode;

        for (int ndnum=fnode-1;ndnum>0;ndnum--)
        {
            if (npt[ndnum]==2)
            {

            	if (nge[ndnum]==0)
                {
                    beginnode=0;
                    endnode=br2node[0][0][1];
                }
                else
                {
                    beginnode=br2node[nge[ndnum]-1][brconn[nge[ndnum]][nbr[ndnum]][0]][1];
                    endnode=br2node[nge[ndnum]][nbr[ndnum]][1];
                }

                Interpolation(beginnode,ndnum,endnode,
                                      Iitnode[beginnode],Iitnode[ndnum],Iitnode[endnode],
                                      tempncoord);
                Interpolation(beginnode,ndnum,endnode,
                                      Ptnode[beginnode],Ptnode[ndnum],Ptnode[endnode],
                                      tempncoord);

                LobeIDnode[ndnum]=LobeIDnode[beginnode];
                Horsfieldnumber[ndnum]=Horsfieldnumber[beginnode];
    		}
        }


        vtkdata<<"POINTS "<<fnode<<" double"<<endl;
        for (int ndnum=0;ndnum<fnode;ndnum++)
        {
            vtkdata<<scale*tempncoord[ndnum][0]<<" "
            	   <<scale*tempncoord[ndnum][1]<<" "
				   <<scale*tempncoord[ndnum][2]<<" "<<endl;
        }

        vtkdata<<"LINES "<<felem<<" "<<felem*3<<endl;
        for (int count=0;count<felem;count++)
        {
            vtkdata<<"2 "<<econn[count][0]-1<<" "<<econn[count][1]-1<<endl;
        }
        vtkdata<<endl;

        cline=0;
        vtkdata<<"POINT_DATA "<<fnode<<endl;
        vtkdata<<"SCALARS Pressure_Pa double"<<endl;
        vtkdata<<"LOOKUP_TABLE default"<<endl;

        for (int ndnum=0;ndnum<fnode;ndnum++)
        {
            cline++;
            vtkdata<<Ptnode[ndnum]<<" ";

            if (cline==5)
            {
                vtkdata<<endl;
                cline=0;
            }
        }
        vtkdata<<endl<<endl;


        cline=0;
        vtkdata<<"SCALARS FlowRate_mL double"<<endl;
        vtkdata<<"LOOKUP_TABLE default"<<endl;

        for (int ndnum=0;ndnum<fnode;ndnum++)
        {
            cline++;
            vtkdata<<Iitnode[ndnum]<<endl<<" ";

            if (cline==5)
            {
                vtkdata<<endl;
                cline=0;
            }
        }
        vtkdata<<endl<<endl;


        cline=0;
        vtkdata<<"SCALARS AwayRadius_mm double"<<endl;
        vtkdata<<"LOOKUP_TABLE default"<<endl;
        for (int ndnum=0;ndnum<fnode;ndnum++)
        {
            cline++;
            vtkdata<<Nodalradius[ndnum]*1000<<" ";

            if (cline==5)
            {
                vtkdata<<endl;
                cline=0;
            }
        }
        vtkdata<<endl<<endl;


        cline=0;
        vtkdata<<"SCALARS AcinarVolume_mL double"<<endl;
        vtkdata<<"LOOKUP_TABLE default"<<endl;
        for (int ndnum=0;ndnum<fnode;ndnum++)
        {
            cline++;
            vtkdata<<NodalAcinarVolume[ndnum]*1000000<<" ";

            if (cline==5)
            {
                vtkdata<<endl;
                cline=0;
            }
        }
        vtkdata<<endl<<endl;


        cline=0;
        vtkdata<<"SCALARS InspiredAcinarFlow_ml double"<<endl;
        vtkdata<<"LOOKUP_TABLE default"<<endl;
        for (int ndnum=0;ndnum<fnode;ndnum++)
        {
            cline++;
            vtkdata<<InspiredAcinarFlow[ndnum]<<" ";

            if (cline==5)
            {
                vtkdata<<endl;
                cline=0;
            }
        }
        vtkdata<<endl<<endl;


        cline=0;
        vtkdata<<"SCALARS ExpiredAcinarFlow_ml double"<<endl;
        vtkdata<<"LOOKUP_TABLE default"<<endl;
        for (int ndnum=0;ndnum<fnode;ndnum++)
        {
            cline++;
            vtkdata<<ExpiredAcinarFlow[ndnum]<<" ";

            if (cline==5)
            {
                vtkdata<<endl;
                cline=0;
            }
        }
        vtkdata<<endl<<endl;


        cline=0;
        vtkdata<<"SCALARS AcinarDye double"<<endl;
        vtkdata<<"LOOKUP_TABLE default"<<endl;

        for (int ndnum=0;ndnum<fnode;ndnum++)
        {
            cline++;
    	vtkdata<<NodalAcinarDensity[ndnum]<<" ";
            if (cline==5)
            {
                vtkdata<<endl;
                cline=0;
            }
        }
        vtkdata<<endl<<endl;


        cline=0;
        vtkdata<<"SCALARS Dye double"<<endl;
        vtkdata<<"LOOKUP_TABLE default"<<endl;

        for (int ndnum=0;ndnum<fnode;ndnum++)
        {
            cline++;
    	    vtkdata<<NodalDyeDensity[ndnum]<<" ";
            if (cline==5)
            {
                vtkdata<<endl;
                cline=0;
            }
        }
        vtkdata<<endl<<endl;


        cline=0;
        vtkdata<<"SCALARS LobeID int"<<endl;
        vtkdata<<"LOOKUP_TABLE default"<<endl;

        for (int ndnum=0;ndnum<fnode;ndnum++)
        {
            cline++;
            vtkdata<<LobeIDnode[ndnum]<<" ";

            if (cline==5)
            {
                vtkdata<<endl;
                cline=0;
            }
        }
        vtkdata<<endl<<endl;


        cline=0;
        vtkdata<<"SCALARS HorsfieldNumber int"<<endl;
        vtkdata<<"LOOKUP_TABLE default"<<endl;

        for (int ndnum=0;ndnum<fnode;ndnum++)
        {
            cline++;
            vtkdata<<Horsfieldnumber[ndnum]<<" ";

            if (cline==5)
            {
                vtkdata<<endl;
                cline=0;
            }
        }
        vtkdata<<endl<<endl;


	    vtkdata.close();
	    vtkscount++;
    }

}


void DataIO::Interpolation(int& point0,
                           int& point1,
                           int& point2,
                           double& value0,
    			   double& value1,
    			   double& value2,
                           vector<vector<double> >&ncoord)
{
        double x0,x1,x2,y0,y1,y2,z0,z1,z2;
        double distance1,distance2;

        x0=ncoord[point0][0];x1=ncoord[point1][0];x2=ncoord[point2][0];
        y0=ncoord[point0][1];y1=ncoord[point1][1];y2=ncoord[point2][1];
        z0=ncoord[point0][2];z1=ncoord[point1][2];z2=ncoord[point2][2];

        distance1=std::sqrt((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0)+(z2-z0)*(z2-z0));
        distance2=std::sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));

        value1=value2-(value2-value0)*distance2/distance1;
}

void DataIO::Write1D(int& breathcount,
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
		     double& InspiredFlowVolume)
{

	if (time==0)
	{
            PQdata<<"Time, NodeNumber, Radius, Pressure, Q, "
    		  <<"BranchPressure, Resistance, AcinarVolume"<<endl;
            PSampledata<<"time, BreathingCount, Ppl, Pt10, Pt20, Pt30, Pt40"<<endl;
            QSampledata<<"time, BreathingCount, Iit00, Iit10, Iit20, Iit30, InFlowVol"<<endl;
	}


    if ((breathcount>=outputstartbreath)&&(breathcount<outputendbreath))
    {
	    for (int gen=fgen-1;gen>=0;gen--)
	    {
	        for (int acn=0;acn<tbr4g[gen];acn++)
	        {
	            if ((time==0)||(svnum==sdtdt))
	            {
	                PQdata<<time
	                <<", "<<br2node[gen][acn][1]
	                <<", "<<brradius[gen][acn]
	                <<", "<<Pt[gen+1][acn]
	                <<", "<<Iit[gen][acn]
	                <<", "<<0.5*(Pt[gen+1][acn]+Pt[gen][brconn[gen][acn][0]])
	                <<", "<<Raw[gen][acn]*0.0000101972
                        <<", "<<AcinarVolume[gen][acn]
                        <<endl;

	                scount++;
	            }
            }
	    }
    }


    if ((time==0)||(svnum==sdtdt))
    {
        PSampledata<<time<<", "<<breathcount<<", "<<Ppl<<", "
		   <<Pt[1][0]<<", "<<Pt[2][0]
		   <<", "<<Pt[3][0]<<", "<<Pt[4][0]<<endl;

        QSampledata<<time<<", "<<breathcount<<", "
        	   <<Iit[0][0]<<", "<<Iit[1][0]<<", "
                   <<Iit[2][0]<<", "<<Iit[3][0]<<", "
		   <<InspiredFlowVolume<<endl;

        svnum=0;
    }

    svnum++;
}


void DataIO::DistalBranchLobe(int& gen,
		              int& acn,
			      std::vector<std::vector<std::vector<int> > >& brconn,
		              std::vector<std::vector<int> >& LobeID,
			      int& DistLobeMax,
			      int& DistLobeMin)
{
	int branch=acn;

    for (int i=1;i<=2;i++)
    {
        int branch1=brconn[gen][acn][i];
    	for (int j=1;j<=2;j++)
    	{
    		int branch2=brconn[gen+1][branch1][j];
            for (int k=1;k<=2;k++)
            {
            	int branch3=brconn[gen+2][branch2][k];
            	for (int m=1;m<=2;m++)
            	{
                    int branch4=brconn[gen+3][branch3][m];
                    for (int n=1;n<=2;n++)
                    {
                    	int branch5=brconn[gen+4][branch4][n];
            	        if (LobeID[gen+5][branch5]>DistLobeMax)
    		            {
    			            DistLobeMax=LobeID[gen+5][branch5];
    		            }
            	        if ((LobeID[gen+5][branch5]>0)&&(LobeID[gen+5][branch5]<DistLobeMin))
		                {
            	            DistLobeMin=LobeID[gen+5][branch5];
		                }
                    }
            	}

            }

    	}
    }
}

