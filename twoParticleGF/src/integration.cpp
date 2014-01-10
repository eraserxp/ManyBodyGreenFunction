/*
 * integration.cpp
 *
 *  Created on: Jan 9, 2014
 *      Author: pxiang
 */

#include "integration.h"

double integrate(double ( * func)(double, void *), void * params, double x_min, double x_max,
               double epsabs, double epsrel, int maxInterval) {

	  gsl_integration_workspace * w = gsl_integration_workspace_alloc (maxInterval);
	  gsl_function F;
	  F.function = func;
	  //set up the additional parameters
	  F.params = params;
	  // integration
	  double result;
	  double error;
	  gsl_integration_qags (&F, x_min, x_max, epsabs, epsrel, maxInterval,
	                        w, &result, &error);
	  gsl_integration_workspace_free (w);
	  return result;
}
