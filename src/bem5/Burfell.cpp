/**
 * @author Samuel Perkin <samuelp12@ru.is>
 * @date 21/01/2014
 * 
 * Copyright (c) 2014 Samuel Perkin
 */

/**
 * @file Burfell.cpp
 *
 * Module that recieves the angle of attack and returns a vector with
 * the interpolated Lift Coefficient in the first value, and the Drag
 * Coefficient in the second value
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <cstdlib>
#include "Burfell.h"

using namespace std;

// Burfell Wind Data 10min intervals
static double BurfellWind[52704];

// Burfell Wind Direction Data in 10min intervals
static double BurfellDir[52704];

// temperature at a height of 10m above ground level at Burfell
static double BurfellTemp[52704];

// The data used for this work is proprietary, obtained from Burfell, Iceland.
// We have to read the data from an external file.
void read_wind_data(const string& windFile) 
{
	ifstream in(windFile.c_str());
	if(!in) {
		cerr << "could not open burfell wind data file '" << windFile << "'" << endl;
		exit(1);
	}
	
	for(unsigned int i=0; i<sizeof(BurfellWind)/sizeof(double); ++i) {
		in >> BurfellWind[i];
	}
}

// The data used for this work is proprietary, obtained from Burfell, Iceland.
// We have to read the data from an external file.
void read_direction_data(const string& dirFile) 
{
	ifstream in(dirFile.c_str());
	if(!in) {
		cerr << "could not open burfell wind data file '" << dirFile << "'" << endl;
		exit(1);
	}
	
	for(unsigned int i=0; i<sizeof(BurfellDir)/sizeof(double); ++i) {
		in >> BurfellDir[i];
	}
}

// The data used for this work is proprietary, obtained from Burfell, Iceland.
// We have to read the data from an external file.
void read_temperature_data(const string& tempFile) 
{
	ifstream in(tempFile.c_str());
	if(!in) {
		cerr << "could not open burfell temperature data file '" << tempFile << "'" << endl;
		exit(1);
	}
	
	for(unsigned int i=0; i<sizeof(BurfellTemp)/sizeof(double); ++i) {
		in >> BurfellTemp[i];
	}
}
		
double Burfell(int iHH, int iWindCutOut, int iWindCutIn, double alpha, 
			   double adPowerCurve[20][5], double dRho, int iTwoTurbines, 
			   double dTurbDist, double dTurbAngle, double dR)
{
	double	dV = 0.0;		//Wind Speed read from Burfell data, used in iteration loop below
	double	dP = 0.0;		//Power at dV, used in iteration loop below
	double	dTime = 10.0/60.0;	//Time step for each data point (i.e. 10/60 hours)
	double	dAEP = 0.0;		//Annual Energy Production, returned by function
	double  dElev = 100;		//Assumed elevation at Burfell
	double  dPress0 = 100600;	//Air pressure at sea level [Pa] (Nawri, 2012)
	double	dTemp0 = 278.5;		//Air temperature at sea level [K] (Nawri, 2012)
	double  dLt = .0063;			//Terrain following temperature lapse rate [K/m] (Nawri, 2012)
	double  dL = .0057;			//Atmospheric temperature lapse rate [K/m] (Nawri, 2012)
	const int iRgas = 287;		//Specific gas constant of dry air [J/K.kg]
	const double g = 9.81;		//Gravitational accelleration constant [m/s2]
	double dTempHH = 0.0;		//Temperature at hub height, initialized for loops below
	double dRhoAdj = 0.0;		//Adjusted air density at hub height [kg/m3] for loops below
	double dShadowAlpha = 0.04; //Assumed rate of shadow radius growth relative to distance from turbine (m/m) [Gonzalez-Longatt, 2011]

	//Air Density Correction: pressure calculation
	double dPress = dPress0 * (pow(dTemp0 / (dTemp0 + (dLt * dElev)),g/(dLt*iRgas))) * (pow((dTemp0 + (dLt * dElev)) / (dTemp0 + (dLt * dElev)+(dL*iHH)),g/(dL*iRgas)));

	//Run normal algorithm if iTwoTurbines is 0
	if(iTwoTurbines == 0)
	{
		for(int iii=0; iii <52703; iii++)
			{

			dV = BurfellWind[iii];
			dTempHH = BurfellTemp[iii] + dL * (iHH - 10.0) + 273.15; //Temperature at Hub Height [K]
			dRhoAdj = dPress / (iRgas*dTempHH);

			//Wind shear correction, using alpha value and hub height
			dV *= pow((static_cast<double>(iHH)/10.0),alpha);

			if(dV > iWindCutOut || dV < iWindCutIn)
			{
				dAEP += 0;
				continue;
			}
			else if(dV > 20)
			{
				dP = adPowerCurve[19][0]; 
				dAEP += dP * dTime * (dRhoAdj/dRho);
				continue;
			}

			int iVDown = static_cast<int>(dV);
			int iVUp = static_cast<int>(dV+1);

			dP = ((adPowerCurve[iVUp-1][0] - adPowerCurve[iVDown-1][0])/(iVUp-iVDown))*(dV-iVDown) + adPowerCurve[iVDown-1][0];
			dAEP += dP * dTime * (dRhoAdj/dRho);
		} //end 'for 0 to 27408' loop
	} //end if two turbines = 0

	//If iTwoTurbines = 1 then use the modified algorithm to double the poweroutput and calculate windshadow
	// x Read in 2nd turbine relative position (distance and angle)
	// - Import 'a' from BEM Loop into Burfell.h by connecting it to the power curve as an additional dimension
	// - Check if either turbine is in wind shadow (i.e. if not, simply double power output)
	// - Calculate extent of windshadow (as a reduced % output of turbine)
	// - Calculate speed in wind shadow
	// - Calculate power output as the combination of the two individual power outputs (different velocities)
	// x Calculate cost as double

	if(iTwoTurbines == 1)
	{
		//Define two turbine specific variables to calculate shadow based on [Gonzalez-Longatt, 2011]:
		double dTurbD; //'D' variable, measuring distance between centre of turbine swept area and wind shadow [m]
		double dRshadow; // radius of the wind shadow [m]
		double dTurbL; //'L' variable, measures distance from centre of wind shadow to chord overlap between shadow and turbine swept area [m]
		double dTurbZ; //'Z' variable, measures maximum height of wind shadow overlap
		double dVshadow; // 'V shadow' describes the wind speed in the shadowed area
		double dAshadow; // 'A shadow' describes the area of the shadow that overlaps the turbine swept area
		double dShadowRatio; // Ratio of what proportion of the turbines swept area is shadowed
		double dVapp; // 'V apparent' is the weighted ratio of shaded and unshaded wind speeds on the turbine
		//double dTurbCT; // Turbine thrust coefficient

		for(int iii=0; iii <52703; iii++)
			{

			dV = BurfellWind[iii];
			dTempHH = BurfellTemp[iii] + dL * (iHH - 10.0) + 273.15; //Temperature at Hub Height [K]
			dRhoAdj = dPress / (iRgas*dTempHH);

			//Wind shear correction, using alpha value and hub height
			dV *= pow((static_cast<double>(iHH)/10.0),alpha);

		//Unshaded turbine power calculations
			if(dV > iWindCutOut || dV < iWindCutIn) //Power = 0 if outside of operating range
			{
				dAEP += 0;
			}
			else if(dV > 20) // Power = maximum if within capacity range
			{
				dP = adPowerCurve[19][0]; 
				dAEP += dP * dTime * (dRhoAdj/dRho);
			}
			else 			//Linearly interpolate power production for unshaded wind turbine and add to dAEP
			{
			int iVDown = static_cast<int>(dV);
			int iVUp = static_cast<int>(dV+1);
			dP = ((adPowerCurve[iVUp-1][0] - adPowerCurve[iVDown-1][0])/(iVUp-iVDown))*(dV-iVDown) + adPowerCurve[iVDown-1][0];
			dAEP += dP * dTime * (dRhoAdj/dRho);
			}

		//Shaded turbine equivalent velocity calculation
			dTurbD = fabs(dTurbDist*sin(fabs(BurfellDir[iii]-dTurbAngle)*(M_PI/180.0)));
			dRshadow = dR + dShadowAlpha*dTurbDist;
			dTurbL = 0.5*(dRshadow-dR+dTurbD);

			//check if 2nd turbine is in shadow
			if(dTurbD-dTurbL < dR) //if turbine is shadowed
			{
				//calculate shadow ratio
				dTurbZ = 2.0*pow(pow(dRshadow,2.0)-pow(dTurbL,2.0),0.5);
				dAshadow = pow(dRshadow,2.0)*acos(dTurbL/dRshadow)+pow(dR,2.0)*acos((dTurbD-dTurbL)/dR)-0.5*dTurbD*dTurbZ;
				dShadowRatio = dAshadow / (M_PI*pow(dR,2.0));
				if(dShadowRatio > 1.0)
					dShadowRatio = 1.0;

				//calculate CT for first turbine
			//int iVDown = static_cast<int>(dV);
			//int iVUp = static_cast<int>(dV+1);
			//dTurbCT = ((adPowerCurve[iVUp-1][4] - adPowerCurve[iVDown-1][4])/(iVUp-iVDown))*(dV-iVDown) + adPowerCurve[iVDown-1][4];

				//calculate vapparent based on shadow area
			dVshadow = dV + dV*(pow(1-0.4,0.5)-1.0)*pow(dR/dRshadow,2.0);
			dVapp = dShadowRatio*dVshadow + (1.0-dShadowRatio)*dV;
			}
			else
			{
				dVapp = dV;
			}

		//Unshaded turbine power calculations
			if(dVapp > iWindCutOut || dVapp < iWindCutIn) //Power = 0 if outside of operating range
			{
				dAEP += 0;
			}
			else if(dVapp > 20) // Power = maximum if within capacity range
			{
				dP = adPowerCurve[19][0]; 
				dAEP += dP * dTime * (dRhoAdj/dRho);
			}
			else 		//Linearly interpolate power production for unshaded wind turbine and add to dAEP
			{
			int iVDown = static_cast<int>(dVapp);
			int iVUp = static_cast<int>(dVapp+1);
			dP = ((adPowerCurve[iVUp-1][0] - adPowerCurve[iVDown-1][0])/(iVUp-iVDown))*(dV-iVDown) + adPowerCurve[iVDown-1][0];
			dAEP += dP * dTime * (dRhoAdj/dRho);
			}
			
		}
	}


	//Convert dAEP from Wh to kWh
	dAEP /= 1000;

	return dAEP;
}
