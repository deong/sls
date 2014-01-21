/**
 * @author Samuel Perkin <samuelp12@ru.is>
 * @date 21/01/2014
 * 
 * Copyright (c) 2014 Samuel Perkin
 */

/**
 * @file RotorProfiles.cpp
 *
 * Module that recieves the angle of attack and returns a vector with
 * the interpolated Lift Coefficient in the first value, and the Drag
 * Coefficient in the second value
 */

#include <iostream>
#include <cassert>
#include "RotorProfiles.h"

using namespace std;

//Formula to look up Lift and Drag coefficients from a table based on a given alpha angle (in degrees)
void LiftDrag(double alpha, double &dLift, double &dDrag)
{
	//2-Dimensional array where entry 0 is alpha = -20 degrees; also column A = CL, column B = CD (Ramsav,1996 p.B-5)
	static const double adLiftDrag[97][2] = {
		{ -0.604, 0.016 },
		{ -0.499, 0.014 },
		{ -0.394, 0.012 },
		{ -0.288, 0.011 },
		{ -0.182, 0.009 },
		{ -0.077, 0.007 },
		{ 0.029, 0.005 },
		{ 0.146, 0.005 },
		{ 0.264, 0.006 },
		{ 0.381, 0.006 },
		{ 0.493, 0.008 },
		{ 0.606, 0.009 },
		{ 0.719, 0.01 },
		{ 0.787, 0.014 },
		{ 0.855, 0.018 },
		{ 0.923, 0.022 },
		{ 0.944, 0.032 },
		{ 0.963, 0.042 },
		{ 0.983, 0.052 },
		{ 1.003, 0.062 },
		{ 1.023, 0.072 },
		{ 1.043, 0.083 },
		{ 1.012, 0.104 },
		{ 0.982, 0.125 },
		{ 0.951, 0.146 },
		{ 0.776, 0.222 },
		{ 0.606, 0.296 },
		{ 0.625, 0.316 },
		{ 0.64, 0.338 },
		{ 0.656, 0.359 },
		{ 0.671, 0.381 },
		{ 0.715, 0.418 },
		{ 0.759, 0.456 },
		{ 0.803, 0.494 },
		{ 0.859, 0.543 },
		{ 0.915, 0.593 },
		{ 0.971, 0.643 },
		{ 1, 0.687 },
		{ 1.028, 0.73 },
		{ 1.056, 0.774 },
		{ 1.084, 0.817 },
		{ 1.112, 0.86 },
		{ 1.131, 0.906 },
		{ 1.151, 0.952 },
		{ 1.171, 0.998 },
		{ 1.191, 1.043 },
		{ 1.211, 1.089 },
		{ 1.214, 1.129 },
		{ 1.218, 1.17 },
		{ 1.223, 1.21 },
		{ 1.227, 1.251 },
		{ 1.231, 1.291 },
		{ 1.222, 1.33 },
		{ 1.213, 1.368 },
		{ 1.205, 1.407 },
		{ 1.196, 1.445 },
		{ 1.187, 1.484 },
		{ 1.178, 1.522 },
		{ 1.169, 1.561 },
		{ 1.16, 1.599 },
		{ 1.152, 1.638 },
		{ 1.143, 1.676 },
		{ 1.134, 1.715 },
		{ 1.125, 1.753 },
		{ 1.116, 1.792 },
		{ 1.108, 1.83 },
		{ 1.099, 1.869 },
		{ 1.075, 1.897 },
		{ 1.05, 1.925 },
		{ 1.025, 1.953 },
		{ 1.001, 1.981 },
		{ 0.976, 2.008 },
		{ 0.951, 2.036 },
		{ 0.926, 2.064 },
		{ 0.901, 2.092 },
		{ 0.876, 2.119 },
		{ 0.851, 2.147 },
		{ 0.82, 2.17 },
		{ 0.787, 2.193 },
		{ 0.754, 2.216 },
		{ 0.722, 2.238 },
		{ 0.689, 2.261 },
		{ 0.657, 2.284 },
		{ 0.624, 2.306 },
		{ 0.591, 2.329 },
		{ 0.559, 2.352 },
		{ 0.526, 2.375 },
		{ 0.489, 2.378 },
		{ 0.451, 2.379 },
		{ 0.414, 2.381 },
		{ 0.376, 2.382 },
		{ 0.338, 2.384 },
		{ 0.301, 2.385 },
		{ 0.263, 2.387 },
		{ 0.225, 2.388 },
		{ 0.187, 2.39 },
		{ 0.15, 2.391 }
	};

   // The adLiftDrag array lists the Lift and Drag coefficients for
   // each integer value of alpha within the range of -6 to 90
   // degrees.
   // 
   // Therefore adLiftDrag[0][1] is the Lift coefficient for an alpha
   // value of -6.
   // 
   // However, the alpha value is passed into the LiftDrag module as a
   // double value, so I use the 'down' and 'up' variables to round
   // the double alpha value to integer values, so that I can then
   // read the Lift/Drag values above and below the actual alpha value
   // so that it can be linearly interpolated.
   // 
   // Of course, if the alpha value is less than -6 (i.e. -7 as you
   // pointed out) you will end up with a negative array index which
   // is physically possible but is beyond the scope of the data I
   // have. If the alpha value is going below -6 it is normally
   // diverging away from where it should be, or the current scenario
   // being tested in the BEM loop is a nonsense scenario.
	if(alpha < -6.0) {
		dLift = -0.604;
		dDrag = 0.016;
		return;
	}
	
	// Find integer value above and below specified alpha, as int values required for look-up table
	int down = static_cast<int>(alpha);
	int up = static_cast<int>(alpha+1);
	
	//Return values for Lift & Drag for integer alpha values bounding actual alpha value
	double dLiftDown = adLiftDrag[6+down][0];
	double dDragDown = adLiftDrag[6+down][1];
	double dLiftUp = adLiftDrag[6+up][0];
	double dDragUp = adLiftDrag[6+up][1];

	//Linearly interpolate Lift and Drag values for alpha based on the values returned from the table
	dLift = ((dLiftUp - dLiftDown) / (up - down))*(alpha-down) + dLiftDown;
	dDrag = ((dDragUp - dDragDown) / (up - down))*(alpha-down) + dDragDown;
}

// Formula that looks up the radius ratio (r/R), chord ratio (C/R) and
// Twist at elements along the rotor blade Based on values listed in a
// pre-defined table
void ChordTwist(int iElement, double dR, double &dRadius, double &dC, double &dTwist, double &dElementWidth)
{
	// 3-Dimensional array where entry 0 is r = 0 metres, entry 20 is r=R metres;
	// Column A = Radius Ratio, column B = Chord Ratio, column C = twist angle (deg) (NREL,2000)
	static double RadChordTwist[26][3] = {
		{ 0,0,0, },
		{ 0.16,0.033,0, },
		{ 0.182,0.063,6.7, },
		{ 0.193,0.08,9.9, },
		{ 0.205,0.098,13.4, },
		{ 0.227,0.135,20.04, },
		{ 0.243,0.132,18.074, },
		{ 0.273,0.129,14.292, },
		{ 0.298,0.126,11.909, },
		{ 0.353,0.12,7.979, },
		{ 0.408,0.115,5.308, },
		{ 0.424,0.113,4.715, },
		{ 0.463,0.109,3.425, },
		{ 0.518,0.104,2.083, },
		{ 0.573,0.098,1.15, },
		{ 0.576,0.098,1.115, },
		{ 0.628,0.093,0.494, },
		{ 0.683,0.087,-0.015, },
		{ 0.727,0.083,-0.381, },
		{ 0.739,0.082,-0.475, },
		{ 0.794,0.076,-0.92, },
		{ 0.849,0.07,-1.352, },
		{ 0.864,0.069,-1.469, },
		{ 0.904,0.065,-1.775, },
		{ 0.959,0.059,-2.191, },
		{ 1,0.055,-2.5, },
	};
	
	dRadius = RadChordTwist[iElement][0]*dR;
	dC = dR*RadChordTwist[iElement][1];
	dTwist = (RadChordTwist[iElement][2]+ RadChordTwist[iElement-1][2])/2.0;
	dElementWidth = dR*RadChordTwist[iElement][0] - dR*RadChordTwist[iElement-1][0]; 
}
