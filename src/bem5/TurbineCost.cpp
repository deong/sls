/**
 * @author Samuel Perkin <samuelp12@ru.is>
 * @date 21/01/2014
 *
 * Copyright (c) 2014 Samuel Perkin
 */

/**
 * @file TurbineCost.cpp
 *
 * calculates the cost of the defined turbine
 */

#include <cmath>
#include "TurbineCost.h"

using namespace std;

// Given number of blades (dB), Radius (dR), Generator Capacity (dGenCap), Hub-Height (iHH) returns costs
int TurbineCost(double dB, double dR, double dGenCap, int iHH, double dAEP, int iPitchFIX,
                int iRPMFIX, double &dCostInitial, double &dCostFixed, double &dCostVariable,
                double dInflation, int iTwoTurbines)
{
	// Convert Generator Capacity from Watts to kiloWatts
	dGenCap /= 1000.0;

	// Annual Replacement Cost
	double dAnnualRepFix = 10.7 * dGenCap;

	// Operations and Maintenance
	double dOMVar = 0.007 * dAEP;

	// Land Lease Costs
	double dLLVar = 0.00108 * dAEP;


	dCostInitial = 4.1184 * pow(dR,3) + 72.864 * pow(dR,2) + 445.05 * dR + (2.0E-5) * pow(dGenCap,3) -
	               0.0741 * pow(dGenCap,2) + 507.3 * dGenCap + 481.378 * pow(iHH,0.4037) * pow(dR,0.8074) +
	               1.87223 * iHH * pow(dR,2) + 1.965 * pow(iHH*dR,1.1736) + 55539.6;
	dCostFixed = dAnnualRepFix;
	dCostVariable = dOMVar + dLLVar;

	//Adjustment factor to bring cost estimates in line with reality (based on real costs at Burfell)
	double dInitialAdj = 2.4022;

	//Adjust for inflation
	dCostInitial *= (1+dInflation)*dInitialAdj;
	dCostFixed *= (1+dInflation);
	dCostVariable *= (1+dInflation);

	//Adjust costs for 2nd turbine, if it exists
	if(iTwoTurbines == 1) {
		dCostInitial *= 2;
		dCostFixed *= 2;
		dCostVariable *= 2;
	}

	return 0;
}
