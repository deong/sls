/**
 * @author Samuel Perkin <samuelp12@ru.is>
 * @date 21/01/2014
 * 
 * Copyright (c) 2014 Samuel Perkin
 */

#ifndef _WEIBULLAEP_H_
#define _WEIBULLAEP_H_

/**
 * @file WeibullAEP.h
 */

double PowInterp(double dVel,double adPowerCurve[20][4]);
double WeibullAEP(double dWeibullk, double dWeibullA, double dWeibullWidth, int dWindCutOut, double adPowerCurve[20][4]);

#endif
