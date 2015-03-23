/**********************************************************
 * This file is generated by 20-sim ANSI-C Code Generator
 *
 *  file:  %FILE_NAME%
 *  model: %MODEL_NAME%
 *  expmt: %EXPERIMENT_NAME%
 *  date:  %GENERATION_DATE%
 *  time:  %GENERATION_TIME%
 *  user:  %USER_NAME%
 *  from:  %COMPANY_NAME%
 *  build: %GENERATION_BUILD%
 **********************************************************/

/* This file describes the integration methods
that are supplied for computation.

   Currently only Euler, RungeKutta2 and RungeKutta4 are supplied,
   but it is easy for the user to add their own
   integration methods with these two as an example.
*/

#ifndef XX_INTEG_H
#define XX_INTEG_H

/* 20-sim include files */
#include "xxtypes.h"

/* the chosen integration method */
#define %INTEGRATION_METHOD_NAME%_METHOD

/* the integration methods */
#ifdef Discrete_METHOD
void %FUNCTIONPREFIX%DiscreteInitialize (void);
void %FUNCTIONPREFIX%DiscreteTerminate (void);
void %FUNCTIONPREFIX%DiscreteStep (void);
#endif 

#ifdef Euler_METHOD
void %FUNCTIONPREFIX%EulerInitialize (void);
void %FUNCTIONPREFIX%EulerTerminate (void);
void %FUNCTIONPREFIX%EulerStep (void);
#endif 

#ifdef RungeKutta2_METHOD
void %FUNCTIONPREFIX%RungeKutta2Initialize (void);
void %FUNCTIONPREFIX%RungeKutta2Terminate (void);
void %FUNCTIONPREFIX%RungeKutta2Step (void);
#endif

#ifdef RungeKutta4_METHOD
void %FUNCTIONPREFIX%RungeKutta4Initialize (void);
void %FUNCTIONPREFIX%RungeKutta4Terminate (void);
void %FUNCTIONPREFIX%RungeKutta4Step (void);
#endif

extern XXBoolean %VARPREFIX%major;

#endif
