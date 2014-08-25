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
	case 7:
		return make_pair(0, 500); // Min-Max range for dTurbDist (i.e. range from 0 to 500m)
		break;
	case 8:
		return make_pair(0, 360); // Min-Max range for dTurbAngle (i.e. range from 0 to 15 degrees)
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
	kvparse::parameter_value(keywords::IBLADE, iBlade, false);
	kvparse::parameter_value(keywords::ITWOTURBINES, iTwoTurbines, false);

	// read in the confidential wind data from the Burfell site
	string windData;
	kvparse::parameter_value(keywords::TURBINE_WIND_DATA, windData);
	read_wind_data(windData);
	
	// read in the confidential temperature data from the Burfell site
	string tempData;
	kvparse::parameter_value(keywords::TURBINE_TEMP_DATA, tempData);
	read_temperature_data(tempData);

	// read in the confidential wind direction data from the Burfell site
	string dirData;
	kvparse::parameter_value(keywords::TURBINE_DIR_DATA, dirData);
	read_direction_data(dirData);
}

bool turbine_problem::evaluate(const vector<double>& p, vector<double>& fit) const
{
	// pull variables out of p and set them to dR, dGenCap, etc.
	double	dR = p[0];
	double	dGenCap = p[1];
	int		iHH = static_cast<int>(floor(p[2]+0.5));
	int		iRPMmin = static_cast<int>(floor(p[3]+0.5)) ;			
	int		iRPMRange = static_cast<int>(floor(p[4]+0.5));
	int		iPitchMin = static_cast<int>(floor(p[5]+0.5));
	int		iPitchRange = static_cast<int>(floor(p[6]+0.5));
	double  dTurbDist = p[7];
	double  dTurbAngle = p[8];
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

	//STEP 1: BEM THEORY CALCULATIONS
	// Initialize Power Curve array, 20 elements for velocity, each with [dPower,dThrust,iRPM,iPitch,da]
	double adPowerCurve[20][5] = { 0 }; 

	//Runs BEM Theory modules and outputs resultant Power Curve to PCurve.csv
	BEMLoop(dB, dR, dGenCap, dRho, dEfficiency, adPowerCurve, iRPMmin, iRPMmax, 
			iPitchMin, iPitchMax, iBlade);


	//Enter 'Burfell.h' to check wind output based on 10min wind interval data and hub height, in kWh
	double dAEPBurf = Burfell(iHH, iWindCutOut, iWindCutIn, alpha, adPowerCurve, dRho, 
							  iTwoTurbines, dTurbDist, dTurbAngle, dR);
	
	//STEP 3: COST CALCULATIONS
	double dCostInitial;
	double dCostFixed;
	double dCostVariable;

	TurbineCost(dB, dR, dGenCap, iHH, dAEPBurf, iPitchFIX, iRPMFIX, dCostInitial,
				dCostFixed, dCostVariable, dInflation, iTwoTurbines);

	//Send results of BEM model calculations to GA fitness vector (Fit 1 = LCOE in USD/MWh, Fit 2 = GWh/year)
	fit[0] = LCOE(dCostInitial,dCostFixed,dCostVariable,iInvestPeriod,dDiscountRate,dAEPBurf)*1000;
	fit[1] = (dAEPBurf/1000000);

	return true;
}
