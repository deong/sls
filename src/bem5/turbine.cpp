#include <iostream>	
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>
#include <cstdlib>
#include "keywords.h"
#include "kvparse/kvparse.h"
#include "turbine.h"
#include "RotorProfiles.h"
#include "BEMLoop.h"
#include "WeibullAEP.h"
#include "TurbineCost.h"
#include "Burfell.h"
#include "LCOE.h"

using namespace std;

turbine_problem::turbine_problem() :
	turbine_dim(7)  // replace with number of GA variables
{
}

turbine_problem::~turbine_problem()
{
}

unsigned int turbine_problem::dimensions() const
{
	return turbine_dim;
}

unsigned int turbine_problem::objectives() const
{
	return 2;
}

std::pair<double,double> turbine_problem::parameter_range(unsigned int param_number) const
{
	switch(param_number)
	{
	case 0:
		return make_pair(0.5, 77.0); // Min-Max range for dR (i.e. rotors can be 0.5m to 77m)
		break;
	case 1:
		return make_pair(100000.0, 7580000.0); // Min-Max range for dGenCap (i.e. range from 1kW to 7.58MW)
		break;
	case 2:
		return make_pair(15, 80); // Min-Max range for iHH (i.e. range from 1m to 160m)
		break;
	case 3:
		return make_pair(1, 15); // Min-Max range for iRPMmin (i.e. range from 1 to 15 rpm)
		break;
	case 4:
		return make_pair(0, 20); // Min-Max range for iRPMRange (i.e. range from 0 to 20 rpm)
		break;
	case 5:
		return make_pair(0, 15); // Min-Max range for iPitchMin (i.e. range from 0 to 15 degrees)
		break;
	case 6:
		return make_pair(0, 15); // Min-Max range for iPitchRange (i.e. range from 0 to 15 degrees)
		break;
	default:
		// shut up the compiler warning
		cerr << "invalid parameter number passed to parameter range" << endl;
		exit(1);
		return make_pair(0, 0);
	}
}

void turbine_problem::initialize()
{
	kvparse::parameter_value(keywords::DB, dB, false);
	kvparse::parameter_value(keywords::DRHO, dRho, false);
	kvparse::parameter_value(keywords::DEFFICIENCY, dEfficiency, false);
	kvparse::parameter_value(keywords::IWINDCUTOUT, iWindCutOut, false);
	kvparse::parameter_value(keywords::IWINDCUTIN, iWindCutIn, false);
	kvparse::parameter_value(keywords::ALPHA, alpha, false);
	kvparse::parameter_value(keywords::IINVESTPERIOD, iInvestPeriod, false);
	kvparse::parameter_value(keywords::DDISCOUNTRATE, dDiscountRate, false);
	kvparse::parameter_value(keywords::IRPMFIX, iRPMFIX, false);
	kvparse::parameter_value(keywords::IPITCHFIX, iPitchFIX, false);
	kvparse::parameter_value(keywords::DINFLATION, dInflation, false);

	// read in the confidential wind data from the Burfell site
	string windData;
	kvparse::parameter_value(keywords::TURBINE_WIND_DATA, windData);
	read_wind_data(windData);
	
	// read in the confidential temperature data from the Burfell site
	string tempData;
	kvparse::parameter_value(keywords::TURBINE_TEMP_DATA, tempData);
	read_temperature_data(tempData);
}

bool turbine_problem::evaluate(const vector<double>& p, vector<double>& fit) const
{

	//double dB=3;				// Number of Blades
	//double dRho=1.225;			// Assumed density of air (kg/m3)
	//double dEfficiency=1;		// Efficiency Assumed for gearbox and generator
	//int	iWindCutOut=20;			//Cutout speed for wind turbine
	//int iWindCutIn=4;			//Cutin speed for wind turbine
	//double alpha=0.14005;		//Wind Shear Factor
	//int iInvestPeriod=20;		//Investment period - years (for LCoE)
	//double dDiscountRate=0.05;	// Discount rate (for LCoE calculation)
	//int iRPMFIX=0;				//if 0 RPM is variable, if 1 RPM is fixed
	//int iPitchFIX=0;			//#if 0 Pitch is variable, if 1 pitch is fixed //add non-GA paramters

	// pull variables out of p and set them to dR, dGenCap, etc.
	double	dR = p[0];
	double	dGenCap = p[1];
	int		iHH = static_cast<int>(floor(p[2]+0.5));
	int		iRPMmin = static_cast<int>(floor(p[3]+0.5)) ;			
	int		iRPMRange = static_cast<int>(floor(p[4]+0.5));
	int		iPitchMin = static_cast<int>(floor(p[5]+0.5));
	int		iPitchRange = static_cast<int>(floor(p[6]+0.5));
	using namespace std;

	//Initialize program stopwatch
	//clock_t t1,t2;
   // t1=clock();

	int		iRPMmax			= iRPMmin+iRPMRange;	// Calculated maximum rotor RPM for BEM theory loop
	int		iPitchMax		= iPitchMin+iPitchRange;// Calculated maximum rotor pitch angle for BEM theory loop

	if(iRPMFIX==1)
	{
		iRPMmin= static_cast<int>(0.5 * (iRPMmin+iRPMmax));
		iRPMmax= iRPMmin+4;
	}

	if(iPitchFIX==1)
		iPitchMax=iPitchMin;

	//cout << "------------------------------------------------------------" << endl;
	//cout << "GA PARAMETERS" << endl;
	//cout << "--------------------" << endl;
	//cout << "blade radius: " << dR << endl;
	//cout << "generator capacity: " << dGenCap/1000 << " kW" << endl;
	//cout << "Hub height: " << iHH << endl;
	//cout << "min and max RPM:" << iRPMmin << " & " << iRPMmax << endl;
	//cout << "min and max Pitch:" << iPitchMin << " & " << iPitchMax << endl << endl;
 
	//cout << "CONSTANTS FROM CFG" << endl;
	//cout << "number of blades: " << dB << endl;
	//cout << "air density: " << dRho << endl;
	//cout << "Efficiency: " << dEfficiency << endl;
	//cout << "wind speed range: " << iWindCutOut << " to " << iWindCutIn << endl;
	//cout << "Wind shear factor: " << alpha  << endl;
	//cout << "investment period: " << iInvestPeriod << endl;
	//cout << "discount rate: " << dDiscountRate << endl;
	//cout << "RPM and Pitch Fix: " << iRPMFIX << " & " << iPitchFIX << endl << endl;

	//STEP 1: BEM THEORY CALCULATIONS
	// Initialize Power Curve array, 20 elements for velocity, each with [dPower,dThrust,iRPM,iPitch]
	double adPowerCurve[20][4] = { 0 }; 

	//Runs BEM Theory modules and outputs resultant Power Curve to PCurve.csv
	BEMLoop(dB,dR,dGenCap,dRho,dEfficiency,adPowerCurve,iRPMmin,iRPMmax,iPitchMin,iPitchMax);

	//Prints Power Curve to Screen
	//cout << "Wind Speed (m/s), Power (W), Thrust (N), RPM, Pitch (deg) " << endl;
	//for( int i = 1; i<21; i++)	
	//cout << i <<"  "<< adPowerCurve[i-1][0] <<"  "<< adPowerCurve[i-1][1] <<"  "<< adPowerCurve[i-1][2] <<"  "<< adPowerCurve[i-1][3] <<"  "<< endl;

	//STEP 2: WIND DISTRIBUTION AND AEP CALCULATIONS
	//Read in Weibull Curve or Wind Data
	//double	dWeibullk		= 1.944;	//Declaration of Weibull shape function
	//double	dWeibullA		= 9.66;		//Declaration of Weibull scale function
	//double	dWeibullWidth	= 0.1;		//Declaration of bin width for wind speed histogram (m/s)
	
	//Enters 'WeibullAEP.h' to calculate the Annual Energy Production in kiloWatt-hours
	//double dAEPWeib = WeibullAEP(dWeibullk,dWeibullA,dWeibullWidth,iWindCutOut,adPowerCurve);

	//cout << "PRODUCTION ESTIMATES" << endl;
	//cout << "--------------------" << endl;
	//cout << "Weibull k-value = " << dWeibullk << endl;
	//cout << "Weibull A-value = " << dWeibullA << endl;
	//cout << "Weibull Estimated AEP = " << dAEPWeib/1000000 << " GWh/year" << endl << endl;

	//Verify AEP with estimation
	//double dCapFact = 0.35;		//Assumed capacity factor
	//double dAEPCheck = dCapFact*(dGenCap/1000)*8766;	//Capacity Factor * Generator Capacity * Hours/year = kWh/year

	//cout << "Capacity Factor Guess = " << dCapFact<< endl;
	//cout << "Capacity Factor Estimated AEP = " << dAEPCheck/1000000 << " GWh/year" << endl;
	//cout << "Difference from Weibull = " << 100*(dAEPWeib-dAEPCheck)/dAEPCheck << " %" << endl << endl; 

	//Enter 'Burfell.h' to check wind output based on 10min wind interval data and hub height, in kWh
	double dAEPBurf = Burfell(iHH,iWindCutOut,iWindCutIn,alpha,adPowerCurve, dRho);

	//cout << "Shear Factor (alpha) = " << alpha << endl;
	//cout << "10-min Time Series Estimate of AEP = " << dAEPBurf/1000000 << " GWh/year" << endl << endl;
	//cout << "Difference from Weibull = " << 100*(dAEPBurf-dAEPWeib)/dAEPWeib << "%" << endl;
	//cout << "Difference from Cap Factor Estimate = " << 100*(dAEPBurf-dAEPCheck)/dAEPCheck << "%" << endl << endl << endl;
	
	//STEP 3: COST CALCULATIONS
	double dCostInitial;
	double dCostFixed;
	double dCostVariable;

	TurbineCost(dB,dR,dGenCap,iHH,dAEPBurf,iPitchFIX, iRPMFIX, dCostInitial,dCostFixed,dCostVariable, dInflation);

	//cout << "COST ESTIMATES" << endl;
	//cout << "--------------------" << endl;

	//cout << "Initial Cost  =      $ " << static_cast<int>(dCostInitial)/1000000 << " million USD" << endl;
	//cout << "Fixed Cost    =      $ " << static_cast<int>(dCostFixed) << " USD per year" << endl;
	//cout << "Variable Cost Rate =   " << (100*(dCostVariable/dAEPBurf)) << " UScents/kWh" << endl;
	//cout << "Variable Cost =      $ " << static_cast<int>(dCostVariable) << " USD per year" << endl << endl;

	//STEP 4: KEY RESULT CALCULATIONS (i.e. AEP, $, LCoE)

	//cout << "KEY RESULTS" << endl;
	//cout << "--------------------" << endl;
	//AEP, Total Cost, LCoE, Fitness Cost, Capacity Factor

	//cout << "AEP = " << dAEPBurf/1000000 << " GWh/year" << endl;
	//cout << "Initial Cost = " << static_cast<int>(dCostInitial)/1000000 << "."<< (static_cast<int>(dCostInitial)/100000)-(10.0*(static_cast<int>(dCostInitial)/1000000)) << " million USD" << endl;
	//cout << "LCoE = " << LCOE(dCostInitial,dCostFixed,dCostVariable,iInvestPeriod,dDiscountRate,dAEPBurf) << " USD/kWh" << endl;
	//cout << "Capacity Factor = " << dAEPBurf / ((dGenCap/1000)*8766) << endl << endl;

	//Finish Clock stopwatch
	//t2=clock();
    //float diff ((float)t2-(float)t1);
	//float seconds = diff / CLOCKS_PER_SEC;
    //cout <<"BEM runtime: "<<seconds<<" seconds"<<endl;
	//cout <<"------------------------------------------------------------"<< endl << endl;
    
	//Send results of BEM model calculations to GA fitness vector (Fit 1 = LCOE in USD/MWh, Fit 2 = GWh/year)
	fit[0] = LCOE(dCostInitial,dCostFixed,dCostVariable,iInvestPeriod,dDiscountRate,dAEPBurf)*1000;
	fit[1] = (dAEPBurf/1000000);

	return true;
}
