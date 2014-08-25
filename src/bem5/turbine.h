#ifndef _TURBINE_H_
#define _TURBINE_H_

#include <iostream>
#include <vector>
#include "problems.h"

using namespace std;

class turbine_problem : public numeric_problem
{
private:
	unsigned int turbine_dim;	// number of GA variables
	double dB;				// Number of Blades
	double dRho;			// Assumed density of air (kg/m3)
	double dEfficiency;		// Efficiency Assumed for gearbox and generator
	int	iWindCutOut;		//Cutout speed for wind turbine
	int iWindCutIn;			//Cutin speed for wind turbine
	double alpha;			//Wind Shear Factor
	int iInvestPeriod;		//Investment period - years (for LCoE)
	double dDiscountRate;	// Discount rate (for LCoE calculation)
	int iRPMFIX;			//if 0 RPM is variable, if 1 RPM is fixed
	int iPitchFIX;			//#if 0 Pitch is variable, if 1 pitch is fixed //add non-GA paramters
	double dInflation;		//Inflation rate since 2002 (29.8% in 2003)
	int iBlade;				//Choose blade type (0 = S809 blade, 1 = NACA 0012H blade)
	int iTwoTurbines;		//Choose to include 2nd turbine (0 = no, 1 = yes)

public:
	turbine_problem();
	virtual ~turbine_problem();

	virtual void initialize();
	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual pair<double,double> parameter_range(unsigned int index) const;
	virtual bool evaluate(const vector<double>& p, vector<double>& fit) const;
};

#endif
