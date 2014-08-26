/**
 * @author Samuel Perkin <samuelp12@ru.is>
 * @date 21/01/2014
 *
 * Copyright (c) 2014 Samuel Perkin
 */

#ifndef _BURFELL_H_
#define _BURFELL_H_

/**
 * @file Burfell.h
 *
 * Module that recieves the angle of attack and returns a vector with
 * the interpolated Lift Coefficient in the first value, and the Drag
 * Coefficient in the second value
 */

#include <string>

void read_wind_data(const std::string& windData);
void read_temperature_data(const std::string& tempData);
void read_direction_data(const std::string& dirData);
double Burfell(int iHH, int iWindCutOut, int iWindCutIn, double alpha,
               double adPowerCurve[20][5], double dRho, int iTwoTurbines,
               double dTurbDist, double dTurbAngle, double dR);

#endif
