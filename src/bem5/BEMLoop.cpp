/**
 * @author Samuel Perkin <samuelp12@ru.is>
 * @date 21/01/2014
 * 
 * Copyright (c) 2014 Samuel Perkin
 */

/**
 * @file BEMLoop.cpp
 *
 * Module that recieves wind turbine characteristics and returns a
 * power curve in the form of an array
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include "BEMLoop.h"
#include "RotorProfiles.h"

#define _USE_MATH_DEFINES

using namespace std;

void BEMLoop(double dB, double dR, double dGenCap, double dRho, double dEfficiency, double adPowerCurve[20][4],
             int iRPMmin, int iRPMmax, double iPitchMin, double iPitchMax)
{
	const int	 iii_max    = 20;		// Maximum Iteration Count (for convergence of a and a')
	const int    iSections	= 25;		// Number of blade sections (i.e. 'dR' elements), should be data points for r/R, C/R and twist minus 1
	iRPMmin = 3;						// Allows model to calculate power output for low wind speeds!

	//Initialize parameters (where required)
	double dRadius			= 0.0; // Element Radius [r] (as per Hansen, 2008)
	double dRadius_old		= 0.0; // Stores the radius of the previous element (m)
	double dC				= 0.0; // Element Chord Width [C] (as per Hansen, 2008)
	double dTwist			= 0.0; // Element Twist angle (in degrees)
	double dElementWidth	= 0.0; // Element width [dr] (as per Hansen, 2008)
	double dCL				= 0.0; // Profile Lift Coefficient based on table (Martin, 2006)
	double dCD				= 0.0; // Profile Drag Coefficient based on table (Martin, 2006)
	double da				= 0.0; // Axial induction factor
	double da_dash			= 0.0; // Rotational induction factor
	double dPhiRads			= 0.0; // Flow Angle [radians]
	double dCn				= 0.0; // Normal Blade Coefficient
	double dCt				= 0.0; // Tangential Blade Coefficient
	double dMoment			= 0.0; // Total Blade Moment (N.m)
	double dPower			= 0.0; // Total Rotor power for given wind speed and settings (W)
	double dThrust			= 0.0; // Total thrust force on the rotor (N)
	int	iRPM				= iRPMmin; // RPM counter
	int iPitch				= 0; // Pitch counter
	int iterations			= 0; // Iteration counter

	// Begin loop for integer value of Wind Speed (1 to 25, as curve is flat or cut-out after 25 m/s)
	for(int iV = 1 ; iV < 21 ; iV++) {
		// Begin loop for RPM (min value to max value)
		for(iRPM = iRPMmin ; iRPM < iRPMmax+1 ; iRPM++) {
			double dGamma = (iRPM * 2.0 * M_PI) / 60.0; // Angular velocity of blade tip
			double dTSR = (dGamma * dR) / iV;
			if(dTSR > 9) {
				break;
			}

			// Begin loop for Pitch angle (min value to max value)
			for(iPitch = iPitchMin ; iPitch < iPitchMax+1 ; iPitch++) {
				iterations++;
				double dPT[iSections];

				// Reset key result parameters before new iteration of BEM algorithm
				dMoment  = 0.0;
				dPower   = 0.0;
				dThrust  = 0.0;

				// Begin loop for blade elements from element 1 to 19
				for(int iElement=1 ; iElement<iSections+1; iElement++) {

					// Step 1: Initialize a and a', typically a = a' = 0
					da = 0.0; //a
					da_dash = 0.0; //a'
					// Store previous radius value
					dRadius_old = dRadius;
					// Initialize and define Prandtl Tip Loss coefficients
					double dPrandtl_f = 0.0;
					double dPrandtl_F = 0.0;

					// Begin loop for convergence of a and a' for given blade section
					for(int iii=1 ; iii< iii_max; iii++) {
						// Call formula from 'RotorProfiles.h' to update blade parameters from table
						ChordTwist(iElement, dR, dRadius, dC, dTwist, dElementWidth);

						// Step 2: Compute the flow angle (phi) using equation 6.7
						dPhiRads = atan( ((1.0-da)*iV) / ((1.0+da_dash)*dGamma*dRadius) ); //Flow Angle [radians]
						// if dPhiRads is NaN, break the loop
						if(dPhiRads != dPhiRads) {
							break;
						}

						// Step 3: Local Angle of Attack using equation 6.6
						double dTheta = static_cast<double>(iPitch) + dTwist;	// Local pitch [degrees]
						double dAlpha = (dPhiRads * 180.0) / M_PI - dTheta;		// Local angle of attack [degrees]

						// Step 4: Determine Lift and Drag Coefficients from table
						// Look up Lift and Drag data values from table in 'RotorProfiles.h', and assign to params
						LiftDrag(dAlpha, dCL, dCD);

						// Chaviaropolous and Hansen correction (2000)
						double dChaviF = 0.0;
						if(dAlpha<15.0) {
							dChaviF = 1.0;
						} else if(dAlpha<25.0) {
							dChaviF = 0.5 * (cos(M_PI*((dAlpha - 15.0)/10))+1);
						} else {
							dChaviF = 0.0;
						}

						dCL += dChaviF*(2.2*pow(dC/dRadius,1.3)*pow(cos(dTwist),4.0)*(1.231 - dCL));
						dCD += dChaviF*(2.2*pow(dC/dRadius,1.3)*pow(cos(dTwist),4.0)*(dCD - 0.005));

						// Step 5: Compute the Normal and Tangential Coefficients (Equations 6.12 and 6.13)

						dCn = dCL * cos(dPhiRads) + dCD * sin(dPhiRads); // Normal coefficient
						dCt = dCL * sin(dPhiRads) - dCD * cos(dPhiRads); // Tangential coefficient

						// Step 6: Calculation of new a and a' (Equations 6.23 and 6.24)

						double dSigma = (dC * dB)/(2.0 * M_PI * dRadius);

						// Prandtl Correction 'f' and 'F' coefficients
						dPrandtl_f = (dB/2.0)*((dR - dRadius)/(dRadius*sin(dPhiRads)));
						// dPrandtl_F = (2.0/M_PI)*(acos(exp(-1.0*dPrandtl_f)));
						dPrandtl_F = (cos(exp((-1.0*dPrandtl_f)/2.5)));

						// Calculate new a and a' values, using Prandtl Tip Loss and Glauret corrections
						double da_new = 0.0;
						if(da<0.2) {
							da_new	= 1.0/(4.0 * dPrandtl_F * ((sin(dPhiRads)*sin(dPhiRads))/(dSigma * dCn)) + 1.0);
						} else {
							double dK = (4.0 * dPrandtl_F * (sin(dPhiRads)*sin(dPhiRads)))/(dSigma*dCn);
							da_new	= 0.5*(2.0 + dK * (1.0 - 2*0.2) - pow(pow(dK*(1.0 - 2.0 * 0.2) + 2,2) + 4.0 * (dK * 0.2 * 0.2 - 1),0.5));
						}

						double da_dash_new	= 1.0/(4.0 * dPrandtl_F * ((sin(dPhiRads)*cos(dPhiRads))/(dSigma * dCt)) - 1.0);

						// Exit conditions for loop
						// calculate error between new and old a values, to check loop
						double dError_a = da - da_new;
						double dError_a_dash = da_dash - da_dash_new;

						// Ensure error values are positive (i.e. remove sign)
						if(dError_a < 0.0) {
							dError_a *= -1.0;
						}
						if(dError_a_dash < 0.0) {
							dError_a_dash *= -1.0;
						}

						// sum errors in a and a' as a total error
						double dTotError = dError_a + dError_a_dash;

						// If errors are less than tolerance, exit loop, else update a and a'
						if(dTotError <0.00001) {
							break;
						} else {
							da = da_new;
							da_dash = da_dash_new;
						}
					}

					// Calculate the incremental thrust (dT) and tangential (dPT) forces of the element
					double dT = 0.5 * dRho * dB * ((iV*iV*(1-da)*(1-da))/(sin(dPhiRads)*sin(dPhiRads)))*dC*dCn*dElementWidth;
					dPT[iElement-1] = 0.5 * dRho * ((iV*(1-da)*dGamma*dRadius*(1+da_dash))/(sin(dPhiRads)*cos(dPhiRads)))*dC*dCt;

					dThrust += dT;

					// Interpolate the incremental tangential forces to calculate the moment on the blade
					if(iElement>1) {
						// Calculate the slope 'M' and intercept 'X' for linear integration of moments
						double dCoeffM = (dPT[iElement-1] - dPT[iElement-2])/(dRadius-dRadius_old);
						double dCoeffX = (dPT[iElement-2] * dRadius - dPT[iElement-1] * dRadius_old) / (dRadius-dRadius_old);
						dMoment += (1.0/3.0) * dCoeffM * ( pow(dRadius,3.0) - pow(dRadius_old,3.0)) + 0.5 * dCoeffX * (pow(dRadius,2.0) - pow(dRadius_old,2.0));
					}
					// End of blade element loop
				}

				// Multiple moment by number of rotor blades to find total rotor moment
				dMoment *= dB;

				// Calculate rotor power in Watts
				dPower = dMoment * dGamma * dEfficiency;

				if(dPower>dGenCap) {
					dPower=dGenCap;
				}

				if(dPower>adPowerCurve[iV-1][0]) {
					adPowerCurve[iV-1][0] = dPower;
					adPowerCurve[iV-1][1] = dThrust;
					adPowerCurve[iV-1][2] = iRPM;
					adPowerCurve[iV-1][3] = iPitch;
				}
				// End of RPM Loop
			}
			// End of Pitch Angle loop
		}

		// Raise the minimum RPM and Pitch to values of previous operating RPM and Pitch, to reduce wasteful calculations
		if(adPowerCurve[iV-1][2] > iRPMmin) {
			iRPMmin = adPowerCurve[iV-1][2]-1;
		}
		if(adPowerCurve[iV-1][3] > iPitchMin) {
			iPitchMin = adPowerCurve[iV-1][3]-1;
		}
	}
}
