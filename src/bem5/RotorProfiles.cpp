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

static const double adLiftDrag1[97][2] = {
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

static const double adLiftDrag2[97][2] = {
	{0,0, },
	{0,0, },
	{0,0, },
	{0,0, },
	{0,0, },
	{0,0, },
	{0,0, },
	{0,0, },
	{0,0, },
	{0,0, },
	{0,0, },
	{0,0, },
	{0,0, },
	{0,0, },
	{0,0, },
	{0,0, },
	{0,0, },
	{0,0, },
	{0,0, },
	{0,0, },
	{0,0.0337, },
	{0.83,0.0338, },
	{0.1534,0.0343, },
	{0.2009,0.0351, },
	{0.2003,0.0359, },
	{0.0328,0.0351, },
	{-0.1413,0.046, },
	{-0.1142,0.058, },
	{-0.0703,0.072, },
	{-0.0215,0.086, },
	{0.0311,0.101, },
	{0.0848,0.117, },
	{0.1387,0.134, },
	{0.1928,0.152, },
	{0.2468,0.171, },
	{0.3008,0.19, },
	{0.3548,0.21, },
	{0.4079,0.231, },
	{0.4606,0.252, },
	{0.5121,0.274, },
	{0.5838,0.297, },
	{0.6161,0.32, },
	{0.6687,0.344, },
	{0.7216,0.369, },
	{0.7744,0.394, },
	{0.8276,0.42, },
	{0.881,0.446, },
	{0.9345,0.473, },
	{0.928,0.505333333333333, },
	{0.9215,0.537666666666667, },
	{0.915,0.57, },
	{0.936,0.605, },
	{0.957,0.64, },
	{0.978,0.675, },
	{0.999,0.71, },
	{1.02,0.745, },
	{1.031,0.78, },
	{1.042,0.815, },
	{1.053,0.85, },
	{1.064,0.885, },
	{1.075,0.92, },
	{1.077,0.951, },
	{1.079,0.982, },
	{1.081,1.013, },
	{1.083,1.044, },
	{1.085,1.075, },
	{1.076,1.103, },
	{1.067,1.131, },
	{1.058,1.159, },
	{1.049,1.187, },
	{1.04,1.215, },
	{1.025,1.241, },
	{1.01,1.267, },
	{0.995,1.293, },
	{0.98,1.319, },
	{0.965,1.345, },
	{0.947,1.37, },
	{0.929,1.395, },
	{0.911,1.42, },
	{0.893,1.445, },
	{0.875,1.47, },
	{0.853,1.2075, },
	{0.831,0.945, },
	{0.809,0.6825, },
	{0.787,0.42, },
	{0.765,0.1575, },
	{0.742,0.459, },
	{0.719,0.7605, },
	{0.696,1.062, },
	{0.673,1.3635, },
	{0.65,1.665, },
	{0.623,1.679, },
	{0.596,1.693, },
	{0.569,1.707, },
	{0.542,1.721, },
	{0.515,1.735, },
	{0.486,1.4236, },
};

static const double RadChordTwist1[26][3] = {
	{0,0,0, },
	{0.16,0.033,0, },
	{0.182,0.063,6.7, },
	{0.193,0.08,9.9, },
	{0.205,0.098,13.4, },
	{0.227,0.135,20.04, },
	{0.243,0.132,18.074, },
	{0.273,0.129,14.292, },
	{0.298,0.126,11.909, },
	{0.353,0.12,7.979, },
	{0.408,0.115,5.308, },
	{0.424,0.113,4.715, },
	{0.463,0.109,3.425, },
	{0.518,0.104,2.083, },
	{0.573,0.098,1.15, },
	{0.576,0.098,1.115, },
	{0.628,0.093,0.494, },
	{0.683,0.087,-0.015, },
	{0.727,0.083,-0.381, },
	{0.739,0.082,-0.475, },
	{0.794,0.076,-0.92, },
	{0.849,0.07,-1.352, },
	{0.864,0.069,-1.469, },
	{0.904,0.065,-1.775, },
	{0.959,0.059,-2.191, },
	{1,0.055,-2.5, },
};

static const double RadChordTwist2[26][3] = {
	{0,0,0, },
	{0.16,0.119352,0, },
	{0.182,0.1209624,0, },
	{0.193,0.1214376,0, },
	{0.205,0.121792,0, },
	{0.227,0.1219744,0, },
	{0.243,0.1217696,0, },
	{0.273,0.1210176,0, },
	{0.298,0.120078333333333,0, },
	{0.353,0.1187912,0, },
	{0.408,0.115368,0, },
	{0.424,0.11396,0, },
	{0.463,0.110145882352941,0, },
	{0.518,0.1034176,0, },
	{0.573,0.0954612,0, },
	{0.576,0.0949944,0, },
	{0.628,0.086444,0, },
	{0.683,0.0764692,0, },
	{0.727,0.0678368,0, },
	{0.739,0.0654176,0, },
	{0.794,0.0537488,0, },
	{0.849,0.041288,0, },
	{0.864,0.037672,0, },
	{0.904,0.0279344,0, },
	{0.959,0.0136884,0, },
	{1,0.00252,0, },
};


//Formula to look up Lift and Drag coefficients from a table based on a given alpha angle (in degrees)
void LiftDrag(double alpha, double &dLift, double &dDrag, int iBlade)
{
	const double (*adLiftDrag)[2];

	if(alpha < -6.0) {
		dLift = -0.604;
		dDrag = 0.016;
		return;
	}

	//2-Dimensional array where entry 0 is alpha = -20 degrees; also column A = CL,
	// column B = CD (Ramsav,1996 p.B-5)
	adLiftDrag = (iBlade == 0) ? adLiftDrag1 : adLiftDrag2;

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
void ChordTwist(int iElement, double dR, double &dRadius, double &dC,
                double &dTwist, double &dElementWidth, int iBlade)
{
	const double (*RadChordTwist)[3] = (iBlade == 0) ? RadChordTwist1 : RadChordTwist2;

	dRadius = RadChordTwist[iElement][0]*dR;
	dC = dR*RadChordTwist[iElement][1];
	dTwist = (RadChordTwist[iElement][2]+ RadChordTwist[iElement-1][2])/2.0;
	dElementWidth = dR*RadChordTwist[iElement][0] - dR*RadChordTwist[iElement-1][0];
}
