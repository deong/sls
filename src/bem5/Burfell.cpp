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
		
double Burfell(int iHH, int iWindCutOut, int iWindCutIn, double alpha, double adPowerCurve[20][4], double dRho)
{
	double	dV = 0.0;		    // Wind Speed read from Burfell data, used in iteration loop below
	double	dP = 0.0;		    // Power at dV, used in iteration loop below
	double	dTime = 10.0/60.0;	// Time step for each data point (i.e. 10/60 hours)
	double	dAEP = 0.0;		    // Annual Energy Production, returned by function
	double  dElev = 100;		// Assumed elevation at Burfell
	double  dPress0 = 100600;	// Air pressure at sea level [Pa] (Nawri, 2012)
	double	dTemp0 = 278.5;		// Air temperature at sea level [K] (Nawri, 2012)
	double  dLt = .0063;		// Terrain following temperature lapse rate [K/m] (Nawri, 2012)
	double  dL = .0057;			// Atmospheric temperature lapse rate [K/m] (Nawri, 2012)
	const int iRgas = 287;		// Specific gas constant of dry air [J/K.kg]
	const double g = 9.81;		// Gravitational accelleration constant [m/s2]
	double dTempHH = 0.0;		// Temperature at hub height, initialized for loops below
	double dRhoAdj = 0.0;		// Adjusted air density at hub height [kg/m3] for loops below

	// Air Density Correction: pressure calculation
	double dPress = dPress0 * (pow(dTemp0 / (dTemp0 + (dLt * dElev)),g/(dLt*iRgas))) *
	                (pow((dTemp0 + (dLt * dElev)) / (dTemp0 + (dLt * dElev)+(dL*iHH)),g/(dL*iRgas)));

	for(int iii=0; iii <52703; iii++) {
		dV = BurfellWind[iii];
		dTempHH = BurfellTemp[iii] + dL * (iHH - 10.0) + 273.15; // Temperature at Hub Height [K]
		dRhoAdj = dPress / (iRgas*dTempHH);

		// Wind shear correction, using alpha value and hub height
		dV *= pow((static_cast<double>(iHH)/10.0),alpha);

		if(dV > iWindCutOut || dV < iWindCutIn) {
			dAEP += 0;
			continue;
		} else if(dV > 20) {
			dP = adPowerCurve[19][0];
			dAEP += dP * dTime * (dRhoAdj/dRho);
			continue;
		}

		int iVDown = static_cast<int>(dV);
		int iVUp = static_cast<int>(dV+1);

		dP = ((adPowerCurve[iVUp-1][0] - adPowerCurve[iVDown-1][0])/(iVUp-iVDown))*(dV-iVDown) + adPowerCurve[iVDown-1][0];
		dAEP += dP * dTime * (dRhoAdj/dRho);
	}

	// Convert dAEP from Wh to kWh
	dAEP /= 1000;

	return dAEP;
}
