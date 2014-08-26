/**
 * @author Samuel Perkin <samuelp12@ru.is>
 * @date 21/01/2014
 *
 * Copyright (c) 2014 Samuel Perkin
 */

#ifndef _LCOE_H_
#define _LCOE_H_

/**
 * @file LCOE.h
 *
 * Module that recieves wind turbine characteristics and returns a
 * power curve in the form of an array.
 */

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double LCOE(double dCostInitial,double dCostFixed,double dCostVariable,int iInvestPeriod, double dDiscountRate, double dAEP)
{
	// Assume uniform AEP (i.e. not a stochastic generation of energy output)
	double dAnnualSum = 0.0;	//Sum of levelized annual costs
	double dAEPSum =0.0;		//Sum of levelized annual electricity production

	for(int iii = 1 ; iii < iInvestPeriod ; iii++) {
		dAnnualSum += (dCostFixed+dCostVariable) / pow((1.0+dDiscountRate),iii);
		dAEPSum += dAEP / pow((1.0+dDiscountRate),iii);
	}

	return ((dCostInitial + dAnnualSum) / dAEPSum);
}

#endif
