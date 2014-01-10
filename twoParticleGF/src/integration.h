/*
 * integration.h
 *
 *  Created on: Jan 9, 2014
 *      Author: pxiang
 */

#ifndef INTEGRATION_H_
#define INTEGRATION_H_

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

double integrate(double ( * func)(double, void *), void * params, double x_min, double x_max,
               double epsabs=1.e-10, double epsrel=1e-7, int maxInterval=100000);

#endif /* INTEGRATION_H_ */
