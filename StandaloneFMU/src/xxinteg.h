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

   Currently the following methods are supported:
   * Discrete
   * Euler
   * RungeKutta2
   * RungeKutta4
   but it is easy for the user to add their own
   integration methods with these two as an example.
*/

#ifndef XX_INTEG_H
#define XX_INTEG_H

/* 20-sim include files */
#include "xxtypes.h"
#include "xxmodel.h" /* For the %VARPREFIX%ModelInstance typedef */

/* the integration methods */
#ifdef Discrete_METHOD
XXBoolean %FUNCTIONPREFIX%DiscreteInitialize (%VARPREFIX%ModelInstance* model_instance);
XXBoolean %FUNCTIONPREFIX%DiscreteTerminate (%VARPREFIX%ModelInstance* model_instance);
XXBoolean %FUNCTIONPREFIX%DiscreteStep (%VARPREFIX%ModelInstance* model_instance, XXDouble outputTime);
#endif 

#ifdef Euler_METHOD
XXBoolean %FUNCTIONPREFIX%EulerInitialize (%VARPREFIX%ModelInstance* model_instance);
XXBoolean %FUNCTIONPREFIX%EulerTerminate (%VARPREFIX%ModelInstance* model_instance);
XXBoolean %FUNCTIONPREFIX%EulerStep (%VARPREFIX%ModelInstance* model_instance, XXDouble outputTime);
#endif 

#ifdef RungeKutta2_METHOD
XXBoolean %FUNCTIONPREFIX%RungeKutta2Initialize (%VARPREFIX%ModelInstance* model_instance);
XXBoolean %FUNCTIONPREFIX%RungeKutta2Terminate (%VARPREFIX%ModelInstance* model_instance);
XXBoolean %FUNCTIONPREFIX%RungeKutta2Step (%VARPREFIX%ModelInstance* model_instance, XXDouble outputTime);
#endif

#ifdef RungeKutta4_METHOD
XXBoolean %FUNCTIONPREFIX%RungeKutta4Initialize (%VARPREFIX%ModelInstance* model_instance);
XXBoolean %FUNCTIONPREFIX%RungeKutta4Terminate (%VARPREFIX%ModelInstance* model_instance);
XXBoolean %FUNCTIONPREFIX%RungeKutta4Step (%VARPREFIX%ModelInstance* model_instance, XXDouble outputTime);
#endif

#endif
