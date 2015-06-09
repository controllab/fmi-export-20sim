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

/* This file contains the implementation of the FMI functions
   Please check the fmiSunctions.h file for more details
*/

/* The FMI related headers */
#include "%FMI_PREFIX%Functions.h"
#include "fmiGUID.h"

/* Our own include files */
#include "xxsubmod.h"

/* The system include files */
#include <string.h>

/* make a FMI_Dll_Export that is for both FMI1 and FMI2 */
%IF%%FMI1%
#define FMI_Dll_Export DllExport
%ENDIF%
%IF%%FMI2%
#define FMI_Dll_Export FMI2_Export
%ENDIF%

/* our own component identifier 
   this can become a pointer to a full structure that is filled when instantiating
   for now a simple non-NULL string pointer is returned
*/
#define xxComponent "%SUBMODEL_NAME%"

/* Inquire version numbers of header files */
FMI_Dll_Export const char* %FMI_PREFIX%GetTypesPlatform()
{
%IF%%FMI1%
    return fmiPlatform;
%ENDIF%
%IF%%FMI2%
    return fmi2TypesPlatform;
%ENDIF%
}
FMI_Dll_Export const char* %FMI_PREFIX%GetVersion()
{
    return %FMI_PREFIX%Version;
}
%IF%%FMI1%
FMI_Dll_Export fmiStatus fmiSetDebugLogging  (fmiComponent c, fmiBoolean loggingOn)
{
    return fmiOK;       /* not yet */
}
%ENDIF%
%IF%%FMI2%
FMI_Dll_Export fmi2Status fmi2SetDebugLogging  (fmi2Component c, fmi2Boolean loggingOn,
											size_t nCategories,
											const fmi2String categories[])
{
    return fmi2OK;       /* not yet */
}
%ENDIF%
/* Data Exchange Functions*/
FMI_Dll_Export %FMI_PREFIX%Status %FMI_PREFIX%GetReal(%FMI_PREFIX%Component c,
									 const %FMI_PREFIX%ValueReference vr[],
									 size_t nvr, %FMI_PREFIX%Real value[])
{
    unsigned int i;
    for (i = 0; i < nvr; i++)
    {
        value [i] = %VARPREFIX%MEMORY[ vr[i] ];
    }
    return %FMI_PREFIX%OK;
}
FMI_Dll_Export %FMI_PREFIX%Status %FMI_PREFIX%GetInteger(%FMI_PREFIX%Component c,
								const %FMI_PREFIX%ValueReference vr[],
								size_t nvr,
								%FMI_PREFIX%Integer value[])
{
    return %FMI_PREFIX%Error;    /* not yet */
}
FMI_Dll_Export %FMI_PREFIX%Status %FMI_PREFIX%GetBoolean(%FMI_PREFIX%Component c,
								const %FMI_PREFIX%ValueReference vr[],
								size_t nvr,
								%FMI_PREFIX%Boolean value[])
{
    return %FMI_PREFIX%Error;    /* not yet */
}
FMI_Dll_Export %FMI_PREFIX%Status %FMI_PREFIX%GetString(%FMI_PREFIX%Component c,
								const %FMI_PREFIX%ValueReference vr[],
								size_t nvr,
								%FMI_PREFIX%String value[])
{
    return %FMI_PREFIX%Error;    /* not yet */
}

FMI_Dll_Export %FMI_PREFIX%Status %FMI_PREFIX%SetReal(%FMI_PREFIX%Component c,
								const %FMI_PREFIX%ValueReference vr[],
								size_t nvr,
								const %FMI_PREFIX%Real value[])
{
    unsigned int i;
    for (i = 0; i < nvr; i++)
    {
        %VARPREFIX%MEMORY[ vr[i] ] = value [i];
    }
    return %FMI_PREFIX%OK;
}
FMI_Dll_Export %FMI_PREFIX%Status %FMI_PREFIX%SetInteger (%FMI_PREFIX%Component c,
								const %FMI_PREFIX%ValueReference vr[],
								size_t nvr,
								const %FMI_PREFIX%Integer value[])
{
    return %FMI_PREFIX%Error;    /* not yet */
}
FMI_Dll_Export %FMI_PREFIX%Status %FMI_PREFIX%SetBoolean(%FMI_PREFIX%Component c,
								const %FMI_PREFIX%ValueReference vr[],
								size_t nvr,
								const %FMI_PREFIX%Boolean value[])
{
	unsigned int i;
	for (i = 0; i < nvr; i++)
	{
		xx_MEMORY[vr[i]] = value[i] ? 1.1 : 0.0;
	}
	return %FMI_PREFIX%OK;    /* temp implementation allowing boolean to real conversion, until proper boolean support is added. */
}
FMI_Dll_Export %FMI_PREFIX%Status %FMI_PREFIX%SetString(%FMI_PREFIX%Component c,
								const %FMI_PREFIX%ValueReference vr[],
								size_t nvr,
								const %FMI_PREFIX%String  value[])
{
    return %FMI_PREFIX%Error;   /* not yet */
}
%IF%%FMI1%
/* FMI functions for Co-Simulation 1.0 */
fmiComponent fmiInstantiateSlave(fmiString instanceName,
								 fmiString GUID,
								 fmiString fmuLocation,
								 fmiString mimeType,
								 fmiReal timeout,
								 fmiBoolean visible,
								fmiBoolean interactive,
								fmiCallbackFunctions functions,
								fmiBoolean loggingOn) 
{
 	/* we should remember the functions pointer in order to make callback functions */

	if (!functions.logger) 
        return NULL; // we cannot even log this problem
    if (!instanceName || strlen(instanceName)==0) { 
        functions.logger(NULL, "?", fmiError, "error",
                "Missing instance name.");
        return NULL;
    }

	/* Check whether the given GUID equals our GUID */
	if( strncmp(GUID, FMI_GUID, strlen(GUID)) != 0 )
	{
		functions.logger(NULL, instanceName, fmiError, "error",
			"Wrong GUID %s. Expected %s.", GUID, FMI_GUID);
		return NULL;
	}

    /* only one static instance for now */
    return (fmiComponent) xxComponent;
}
%ENDIF%

%IF%%FMI2%
fmi2Component fmi2Instantiate(fmi2String instanceName,
								fmi2Type fmuType,
								fmi2String fmuGUID,
								fmi2String fmuResourceLocation,
								const fmi2CallbackFunctions* functions,
								fmi2Boolean visible,
								fmi2Boolean loggingOn)
{
	/* we should remember the functions pointer in order to make callback functions */

    if (!functions) 
        return NULL; // we cannot even log this problem
    if (!functions->logger) 
        return NULL; // we cannot even log this problem
    if (!instanceName || strlen(instanceName)==0) { 
        functions->logger(NULL, "?", fmi2Error, "error",
                "Missing instance name.");
        return NULL;
    }
	/* Check whether the given GUID equals our GUID */
	if( strncmp(fmuGUID, FMI_GUID, strlen(fmuGUID)) != 0 )
	{
		functions->logger(NULL, instanceName, fmi2Error, "error",
			"Wrong GUID %s. Expected %s.", fmuGUID, FMI_GUID);
		return NULL;
	}
	/* check if we are setup for co-simulation, that's the only possible option for now */
	if( fmuType != fmi2CoSimulation )
	{
		functions->logger(NULL, instanceName, fmi2Error, "error",
			"FMU can only be used for Co-Simulation, not for Model Exchange");
		return NULL;
	}

    /* only one static instance for now */
    return (fmi2Component) xxComponent;
}
%ENDIF%
%IF%%FMI1%
fmiStatus fmiInitializeSlave(fmiComponent c,
							 fmiReal tStart,
							 fmiBoolean StopTimeDefined,
							 fmiReal tStop) 
{
    /* copy the arguments */
    %VARPREFIX%start_time = tStart;
    if (StopTimeDefined == fmiTrue) 
    {
        %VARPREFIX%finish_time = tStop;
    }

    /* initialize the submodel itself */
    %FUNCTIONPREFIX%InitializeSubmodel (tStart);

    /* all done */
    return fmiOK;
}
%ENDIF%
%IF%%FMI2%
fmi2Status fmi2SetupExperiment(fmi2Component c,
							fmi2Boolean toleranceDefined,
							fmi2Real tolerance,
							fmi2Real startTime,
							fmi2Boolean stopTimeDefined,
							fmi2Real stopTime)
{
    /* copy the arguments */
    %VARPREFIX%start_time = startTime;
    if (stopTimeDefined == fmi2True) 
    {
        %VARPREFIX%finish_time = stopTime;
    }

    /* initialize the submodel itself */
    %FUNCTIONPREFIX%InitializeSubmodel (startTime);

    /* all done */
    return fmi2OK;
}
fmi2Status fmi2EnterInitializationMode(fmi2Component c)
{
   /* nothing to do for now */
   return fmi2OK;
}
fmi2Status fmi2ExitInitializationMode(fmi2Component c)
{
   /* nothing to do for now */
   return fmi2OK;
}
%ENDIF%
%IF%%FMI1%
fmiStatus fmiTerminateSlave(fmiComponent c) 
%ENDIF%
%IF%%FMI2%
fmi2Status fmi2Terminate(fmi2Component c) 
%ENDIF%
{
	/* Perform the final calculations */
	%FUNCTIONPREFIX%TerminateSubmodel (%VARPREFIX%%XX_TIME%);

    /* all done */
    return %FMI_PREFIX%OK;
}
%IF%%FMI1%
fmiStatus fmiResetSlave(fmiComponent c) 
%ENDIF%
%IF%%FMI2%
fmi2Status fmi2Reset(fmi2Component c) 
%ENDIF%
{
    /* initialize the submodel itself */
    %FUNCTIONPREFIX%InitializeSubmodel (%VARPREFIX%start_time);

    /* all done */
    return %FMI_PREFIX%OK;
}

%IF%%FMI1%
void fmiFreeSlaveInstance(fmiComponent c) 
%ENDIF%
%IF%%FMI2%
void fmi2FreeInstance(fmi2Component c) 
%ENDIF%
{
    /* only one static instance (done automatically) */
}

%FMI_PREFIX%Status %FMI_PREFIX%SetRealInputDerivatives(%FMI_PREFIX%Component c,
									const %FMI_PREFIX%ValueReference vr[], size_t nvr,
									const %FMI_PREFIX%Integer order[],
									const %FMI_PREFIX%Real value[]) 
{
    /* not yet */
    return %FMI_PREFIX%Error;
}

%FMI_PREFIX%Status %FMI_PREFIX%GetRealOutputDerivatives(%FMI_PREFIX%Component c,
									const %FMI_PREFIX%ValueReference vr[],
									size_t nvr,
									const %FMI_PREFIX%Integer order[],
									%FMI_PREFIX%Real value[]) 
{
    /* not yet */
    return %FMI_PREFIX%Error;
}

%IF%%FMI1%
fmiStatus fmiCancelStep(fmiComponent c) 
%ENDIF%
%IF%%FMI2%
fmi2Status fmi2CancelStep(fmi2Component c) 
%ENDIF%
{
    /* not yet */
    return %FMI_PREFIX%Error;
}

%IF%%FMI1%
fmiStatus fmiDoStep(fmiComponent c,
					fmiReal currentCommunicationPoint,
					fmiReal communicationStepSize,
					fmiBoolean newStep) 
%ENDIF%
%IF%%FMI2%
fmi2Status fmi2DoStep(fmi2Component c,
					fmi2Real currentCommunicationPoint,
					fmi2Real communicationStepSize,
					fmi2Boolean noSetFMUStatePriorToCurrentPoint)
%ENDIF%
{
    /* Treat also case of zero step, i.e. during an event iteration */
    if (communicationStepSize == 0) 
    {
        return %FMI_PREFIX%OK;
    }

    /* as long as we are not passed our communication point */
    while (%VARPREFIX%%XX_TIME% < (currentCommunicationPoint + communicationStepSize))
    {
        /* check for termination first */
        if ((%VARPREFIX%%XX_TIME% > %VARPREFIX%finish_time) || (%VARPREFIX%stop_simulation == XXTRUE))
        {
            /* we're done */
            return %FMI_PREFIX%Error;
        }
        
		/* Call the submodel to calculate the output, and increase the time as well */
		%FUNCTIONPREFIX%CalculateSubmodel (%VARPREFIX%%XX_TIME%);
    }
    
    /* for now */
    return %FMI_PREFIX%OK;
}

%FMI_PREFIX%Status %FMI_PREFIX%GetStatus(%FMI_PREFIX%Component c, const %FMI_PREFIX%StatusKind s, %FMI_PREFIX%Status* value) 
{
    /* all fine? */
    return %FMI_PREFIX%OK;
}

%FMI_PREFIX%Status %FMI_PREFIX%GetRealStatus(%FMI_PREFIX%Component c, const %FMI_PREFIX%StatusKind s, %FMI_PREFIX%Real* value)
{
    /* all fine? */
    return %FMI_PREFIX%OK;
}

%FMI_PREFIX%Status %FMI_PREFIX%GetIntegerStatus(%FMI_PREFIX%Component c, const %FMI_PREFIX%StatusKind s, %FMI_PREFIX%Integer* value)
{
    /* not yet */
    return %FMI_PREFIX%Discard;
}

%FMI_PREFIX%Status %FMI_PREFIX%GetBooleanStatus(%FMI_PREFIX%Component c, const %FMI_PREFIX%StatusKind s, %FMI_PREFIX%Boolean* value)
{
    /* not yet */
    return %FMI_PREFIX%Discard;
}

%FMI_PREFIX%Status %FMI_PREFIX%GetStringStatus(%FMI_PREFIX%Component c, const %FMI_PREFIX%StatusKind s, %FMI_PREFIX%String*  value)
{
    /* not yet */
    return %FMI_PREFIX%Discard;
}

%IF%%FMI2%
fmi2Status fmi2GetContinuousStates(fmi2Component c, fmi2Real x[], size_t nx)
{
	unsigned i;
	if( %NUMBER_STATES% != nx )
	{
		return fmi2Error;
	}
	for( i = 0; i < %NUMBER_STATES%; ++i)
	{
		x[i] = %VARPREFIX%%XX_STATE_ARRAY_NAME%[i];
	}
	return fmi2OK;
}
%ENDIF%

%IF%%FMI2%
fmi2Status fmi2GetFMUstate (fmi2Component c, fmi2FMUstate* FMUstate) {
    
/* not yet */
    return fmi2Discard;
}
fmi2Status fmi2SetFMUstate (fmi2Component c, fmi2FMUstate FMUstate) {
    /* not yet */
    return fmi2Discard;
}
fmi2Status fmi2FreeFMUstate(fmi2Component c, fmi2FMUstate* FMUstate) {
    /* not yet */
    return fmi2Discard;
}
fmi2Status fmi2SerializedFMUstateSize(fmi2Component c, fmi2FMUstate FMUstate, size_t *size) {
    /* not yet */
    return fmi2Discard;
}
fmi2Status fmi2SerializeFMUstate (fmi2Component c, fmi2FMUstate FMUstate, fmi2Byte serializedState[], size_t size) {
    /* not yet */
    return fmi2Discard;
}
fmi2Status fmi2DeSerializeFMUstate (fmi2Component c, const fmi2Byte serializedState[], size_t size,
                                    fmi2FMUstate* FMUstate) {
    /* not yet */
    return fmi2Discard;
}

fmi2Status fmi2GetDirectionalDerivative(fmi2Component c, const fmi2ValueReference vUnknown_ref[], size_t nUnknown,
                                        const fmi2ValueReference vKnown_ref[] , size_t nKnown,
                                        const fmi2Real dvKnown[], fmi2Real dvUnknown[]) {
    /* not yet */
    return fmi2Discard;
}
%ENDIF%