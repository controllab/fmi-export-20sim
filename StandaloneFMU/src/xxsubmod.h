/**********************************************************
 * This file is generated by 20-sim ANSI-C Code Generator
 *
 *  file:  %FILE_NAME%
 *  subm:  %SUBMODEL_NAME%
 *  model: %MODEL_NAME%
 *  expmt: %EXPERIMENT_NAME%
 *  date:  %GENERATION_DATE%
 *  time:  %GENERATION_TIME%
 *  user:  %USER_NAME%
 *  from:  %COMPANY_NAME%
 *  build: %GENERATION_BUILD%
 **********************************************************/

/* This file describes the model functions
   that are supplied for computation.

   The model itself is the xxmodel.c file
*/

#ifndef XX_SUBMOD_H
#define XX_SUBMOD_H

/* Our own include files */
#include "xxmodel.h"



/* The submodel functions */
void %FUNCTIONPREFIX%InitializeSubmodel (XXModelInstance* %VARPREFIX%model_instance, XXDouble t);
void %FUNCTIONPREFIX%CalculateSubmodel (XXModelInstance* %VARPREFIX%model_instance, XXDouble t);
void %FUNCTIONPREFIX%TerminateSubmodel (XXModelInstance* %VARPREFIX%model_instance, XXDouble t);

#endif
