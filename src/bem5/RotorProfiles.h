/**
 * @author Samuel Perkin <samuelp12@ru.is>
 * @date 21/01/2014
 * 
 * Copyright (c) 2014 Samuel Perkin
 */

#ifndef _ROTORPROFILES_H_
#define _ROTORPROFILES_H_

/**
 * @file RotorProfiles.h
 *
 * Module that recieves the angle of attack and returns a vector with
 * the interpolated Lift Coefficient in the first value, and the Drag
 * Coefficient in the second value
 */

// Formula to look up Lift and Drag coefficients from a table based on a given alpha angle (in degrees)
void LiftDrag(double alpha, double &dLift, double &dDrag);

// Formula that looks up the radius ratio (r/R), chord ratio (C/R) and Twist at elements along the rotor blade
// Based on values listed in a pre-defined table 
void ChordTwist(int iElement, double dR, double &dRadius, double &dC, double &dTwist, double &dElementWidth);

#endif
