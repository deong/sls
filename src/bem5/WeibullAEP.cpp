/**
 * @author Samuel Perkin <samuelp12@ru.is>
 * @date 21/01/2014
 *
 * Copyright (c) 2014 Samuel Perkin
 */

/**
 * @file WeibullAEP.cpp
 */

#include <cmath>

using namespace std;

// Given some double Velocity 'dVel', interpolate the power curve to return P(dVel)
double PowInterp(double dVel,double adPowerCurve[20][4])
{
	double dPowerDown = 0.0;

	//Find integer wind speeds above and below dVel
	int iVDown = static_cast<int>(dVel);	// Find integer windspeed less than iVel
	int iVUp = static_cast<int>(dVel+1);	// Find integer windspeed above iVel

	// Return power values at integer wind speeds
	double dPowerUp = adPowerCurve[iVUp-1][0];		// Find power at iVLowUp

	if (iVDown > 0) {
		dPowerDown = adPowerCurve[iVDown-1][0];    // Find power at iVLowUp
	}

	double dPower = ((dPowerUp-dPowerDown)/(iVUp-iVDown))*(dVel-iVDown)+dPowerDown;

	return dPower;
}

// Given Weibull parameters k and A, histogram resolution (i.e. bin width), and Power Curve
// returns the AEP (annual energy production) in kiloWatt-hours
double WeibullAEP(double dWeibullk, double dWeibullA, double dWeibullWidth, int dWindCutOut, double adPowerCurve[20][4])
{
	// Initialize parameters for Weibull-AEP calculation loop
	int		iWindSteps		= (20)/dWeibullWidth;	// # of steps for Weibull-AEP calculation (up to 20 m/s)
	double	dVLow			= 0.0;					// Lower bound of Wind Speed bin
	double	dVHigh			= 0.0;					// Upper bound of Wind Speed bin
	double	dCDFHigh		= 0.0;					// Weibull CDF of upper bound wind speed
	double	dCDFLow			= 0.0;					// Weibull CDF of lower bound wind speed
	double	dCDF			= 0.0;					// Weibull CDF of bin
	double	dHoursPerYear	= 8766.0;				// Hours in a year (hrs/year)
	double	dHours			= 0.0;					// Hours of operation at windspeed dV (hrs)
	double	dPowerLow		= 0.0;					// Power output at 'dVLow'
	double	dPowerHigh		= 0.0;					// Power output at 'dVHigh'
	double	dPower			= 0.0;					// Incremental Power Output (W)
	double  dAEP			= 0.0;					// Annual Energy Production (kWh)

	// AEP calculation loop
	for (int iWindLoop = 1; iWindLoop<iWindSteps+1; iWindLoop++) {
		dVLow = (iWindLoop - 1.0) * dWeibullWidth;	// Calculate lower limit to bin
		dVHigh = iWindLoop * dWeibullWidth;			// Calculate upper limit to bin

		// Calculate Cumulative Distribution of Upper and Lower Bounds, then frequency of bin
		dCDFHigh = 1.0 - exp((-1.0)*(pow((dVHigh/dWeibullA),dWeibullk)));
		dCDFLow = 1.0 - exp((-1.0)*(pow((dVLow/dWeibullA),dWeibullk)));
		dCDF = dCDFHigh-dCDFLow;

		// Calculate number of hours in a year that wind speed is in bin
		dHours = dHoursPerYear * dCDF;

		// Find P(dVLow) and P(dVHigh) then average
		dPowerLow = PowInterp(dVLow,adPowerCurve);
		dPowerHigh = PowInterp(dVHigh,adPowerCurve);
		dPower = 0.5*(dPowerLow+dPowerHigh);

		// Add incremental AEP (in kWh) produced at 'dVLow < Wind Speed < dVHigh' to total AEP
		dAEP += (dPower*dHours)/1000.0;
	}

	// Calculate production between wind speed of 20 m/s and cut-out speed (defined above)
	dCDFHigh = 1.0 - exp((-1.0)*(pow((dWindCutOut/dWeibullA),dWeibullk)));
	dCDFLow = 1.0 - exp((-1.0)*(pow((20.0/dWeibullA),dWeibullk)));
	dCDF = dCDFHigh-dCDFLow;

	// Calculate number of hours in a year that wind speed is in bin
	dHours = dHoursPerYear * dCDF;
	dPower = adPowerCurve[19][0];
	dAEP += (dPower*dHours)/1000.0;

	return dAEP;
}

