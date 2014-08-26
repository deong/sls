/**
 * @author Samuel Perkin <samuelp12@ru.is>
 * @date 21/01/2014
 *
 * Copyright (c) 2014 Samuel Perkin
 */

#ifndef _BEMLOOP_H_
#define _BEMLOOP_H_

/**
 * @file BEMLoop.h
 *
 * Module that recieves wind turbine characteristics and returns a power curve
 * in the form of an array
 */

void BEMLoop(double dB, double dR, double dGenCap, double dRho, double dEfficiency,
             double adPowerCurve[20][5], int iRPMmin, int iRPMmax, double iPitchMin,
             double iPitchMax, int iBlade);

#endif
