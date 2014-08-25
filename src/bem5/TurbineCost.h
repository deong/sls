/**
 * @author Samuel Perkin <samuelp12@ru.is>
 * @date 21/01/2014
 * 
 * Copyright (c) 2014 Samuel Perkin
 */

#ifndef _TURBINECOST_H_
#define _TURBINECOST_H_

/**
 * @file TurbineCost.h
 *
 * calculates the cost of the defined turbine
 */

//Given number of blades (dB), Radius (dR), Generator Capacity (dGenCap), Hub-Height (iHH) returns costs
int TurbineCost(double dB, double dR, double dGenCap, int iHH, double dAEP, int iPitchFIX, 
				int iRPMFIX, double &dCostInitial, double &dCostFixed, double &dCostVariable, 
				double dInflation, int iTwoTurbines);

#endif
