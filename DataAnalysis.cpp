/*
  licence: GNU GPL v3 licence

  Full-Scale Airway Network (FAN) flow model: DataAnalysis.cpp
  Copyright (C) 2019  Minsuok Kim

*/

#include "DataAnalysis.h"

void DataAnalysis::LungResistance(int& fgen,
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
                                  std::ofstream& GenResi)
{
	double TotalAcinarVolume,sumRaw;

	if (time==0)
	{
	    GenResi<<"time, Generation, Avg_Resi, Tot_Resi"<<endl;
	    LungResi<<"time, Lung_Volume, Airway_Resi, Lung_Resi"<<endl;

   }

    if (((breathcount>=outputstartbreath)&&(breathcount<outputendbreath))&&((time==0)||(svnum==sdtdt)))
    {
    	TotalAcinarVolume=0.0;

	    for (int gen=fgen-1;gen>=0;gen--)
	    {
	    	sumRaw=0.0;
	        for (int acn=0;acn<tbr4g[gen];acn++)
	        {
	        	if (br[gen][acn]==1)
	        	{
	        		Raw_sum[gen][acn]=Raw[gen][acn];
	        		TotalAcinarVolume+=AcinarVolume[gen][acn];
	        	}
	        	else if ((br[gen][acn]==2)&&(fabs(Raw_sum[gen+1][brconn[gen][acn][1]])>0)
	        			                  &&(fabs(Raw_sum[gen+1][brconn[gen][acn][2]])>0))
                {
                    Raw_sum[gen][acn]=Raw[gen][acn]
								    +1/(1/fabs(Raw_sum[gen+1][brconn[gen][acn][1]])
                    		           +1/fabs(Raw_sum[gen+1][brconn[gen][acn][2]]));
                }
	        	sumRaw+=fabs(Raw[gen][acn]);
	        }

            GenResi<<time<<", "<<gen<<", "<<sumRaw*0.000010197/tbr4g[gen]<<sumRaw*0.000010197<<endl;
	    }
        LungResi<<time<<", "<<TotalAcinarVolume*1000<<", "<<Raw_sum[0][0]*0.000010197<<", "
        		<<newdPpl*0.0101972/(Iit[0][0]/1000)<<endl;
    }
}

void DataAnalysis::ComputedParameters(int& fgen,
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
		                      std::vector<std::vector<double> >& Pt_PeakInhal)
{
    if ((breathcount>=outputstartbreath)&&(breathcount<outputendbreath))
    {
    	if(Q[0][0]>0.0)
    	{
            for (int gen=fgen-1;gen>=0;gen--)
            {
    	        for (int acn=0;acn<tbr4g[gen];acn++)
                {
    	        	InspiredIit[gen][acn]+=Q[gen][acn]*1000000*dtime;
                }
            }

            if(((Q[0][0]-Qold[0][0])*(Qold[0][0]-Qold2[0][0]))<0)
            {
            	for (int gen=fgen-1;gen>=0;gen--)
                {
                    for (int acn=0;acn<tbr4g[gen];acn++)
                    {
                	    Pt_PeakInhal[gen+1][acn]=Ptold[gen+1][acn];
                    }
                }
            }

        }
    	else if ((Q[0][0]<0)&&(Qold[0][0]>0))
		{
            for (int gen=fgen-1;gen>=0;gen--)
            {
    	        for (int acn=0;acn<tbr4g[gen];acn++)
                {
    	        	InspiredIit[gen][acn]+=0.0;
                }
            }
		}

    }

}

