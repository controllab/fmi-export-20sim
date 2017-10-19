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

%IF%%FMI2%
#if !defined FMI2_FUNCTION_PREFIX && !defined NO_FUNCTION_PREFIX
#define FMI2_FUNCTION_PREFIX %SUBMODEL_NAME%_
#endif
%ENDIF%

/* The FMI related headers */
#include "%FMI_PREFIX%Functions.h"
#include "fmiGUID.h"

/* Our own include files */
#include "xxmodel.h"
#include "xxsubmod.h"

/* The system include files */
#include <string.h>
#include <ctype.h>

/* make a FMI_Dll_Export that is for both FMI1 and FMI2 */
%IF%%FMI1%
#define FMI_Dll_Export DllExport
%ENDIF%
%IF%%FMI2%
#define FMI_Dll_Export FMI2_Export
%ENDIF%

#if defined WIN32 || defined WIN64
static const char native_path_separator = '\\';
static const char foreign_path_separator = '/';
#else
static const char native_path_separator = '/';
static const char foreign_path_separator = '\\';
#endif

/**
 * Convert an uri provided by the the co-simulator to a native path
 * @param uri Input path. native and file:/ and file:/// uri's right now
 * @return Natve path that ends with a path separator.
 * Note that the caller is responsible for free-ing the allocated string buffer
 * memory
 */
const char* URIToNativePath(%VARPREFIX%ModelInstance* model_instance, const char* uri)
{
	unsigned int path_start = 0;
	char* path = NULL;
	unsigned int path_len = 0;
	unsigned int uri_len = 0;
	unsigned int i = 0;
	unsigned int j = 0;
	char buf[3] = "00";

	if (!uri || model_instance == NULL)
	{
		return NULL;
	}

	uri_len = (unsigned int) strlen(uri);

	if (uri_len == 0)
	{
		return NULL;
	}

	/* Check if we got a file:/// uri */
	if (strncmp(uri, "file:///", 8) == 0)
	{
		if (uri[9] == ':')
		{
			/* Windows drive letter in the URI (e.g. file:///c:/ uri */
			/* Remove the file:/// */
			path_start = 8;
		}
		else
		{
			/* Remove the file:// but keep the third / */
			path_start = 7;
		}
	}
#if defined WIN32 || defined WIN64
	/* Check if we got a file://hostname/path uri */
	else if (strncmp(uri, "file://", 7) == 0)
	{
		/* Convert to a network share path: //hostname/path */
		path_start = 5;
	}
#endif
	/* Check if we got a file:/ uri */
	else if (strncmp(uri, "file:/", 6) == 0)
	{
		if (uri[7] == ':')
		{
			/* Windows drive letter in the URI (e.g. file:/c:/ uri */
			/* Remove the file:/ */
			path_start = 6;
		}
		else
		{
			/* Remove the file: but keep the / */
			path_start = 5;
		}
	}
	/* Assume that it is a native path */
	else
	{
		path_start = 0;
	}

	/* Check the length of the remaining string */
	path_len = (int)strlen(&uri[path_start]);
	if (path_len == 0)
	{
		return NULL;
	}

	/* Allocate memory for the return value including terminating \0 and extra path separator */
%IF%%FMI1%
	if (model_instance->fmiCallbackFunctions.allocateMemory != NULL)
	{
		path = (char*) model_instance->fmiCallbackFunctions.allocateMemory(path_len + 2, sizeof(char));
	}
%ENDIF%
%IF%%FMI2%
	if ((model_instance->fmiCallbackFunctions) &&( model_instance->fmiCallbackFunctions->allocateMemory != NULL))
	{
		path = (char*) model_instance->fmiCallbackFunctions->allocateMemory(path_len + 2, sizeof(char));
	}
%ENDIF%
	else
	{
		path = (char*) malloc(path_len + 2);
	}

	/* Copy the remainder of the uri and replace all percent encoded character
	* by their ASCII character and translate slashes to backslashes on Windows
	* and backslashes to slashes on other OSses
	*/
	for (i = path_start, j = 0; i < uri_len; i++, j++)
	{
		if (uri[i] == '%')
		{
			/* Replace the precent-encoded hexadecimal digits by its US-ASCII
			* representation */
			if (i < uri_len - 2)
			{
				if ((isxdigit(uri[i + 1])) && (isxdigit(uri[i + 2])))
				{
					strncpy(buf, uri + i + 1, 2);
					path[j] = (unsigned char)strtol(buf, NULL, 16);
					i += 2;
					path_len -= 2;
					if (path[j] == foreign_path_separator)
					{
						/* Translate slashes to backslashes on Windows and backslashes to slashes on other OSses */
						path[j] = native_path_separator;
					}
				}
				else
				{
					/* Not percent encoded, keep the % */
					path[j] = uri[i];
				}
			}
			else
			{
				/* Not percent encoded, keep the % */
				path[j] = uri[i];
			}
		}
		else if (uri[i] == foreign_path_separator)
		{
			/* Translate slashes to backslashes on Windows and backslashes to slashes on other OSses */
			path[j] = native_path_separator;
		}
		else
		{
			/* Just copy the character */
			path[j] = uri[i];
		}
	}

	/* Check if we need to add a path separator at the end */
	if (path[path_len - 1] == native_path_separator)
	{
		path[path_len] = '\0';
	}
	else
	{
		path[path_len] = native_path_separator;
	}
	/* Make sure that the string is always NULL terminated */
	path[path_len + 1] = '\0';

	return path;
}

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
	size_t i;
	%VARPREFIX%ModelInstance* model_instance = (%VARPREFIX%ModelInstance*) c;

	for (i = 0; i < nvr; i++)
	{
		value [i] = model_instance->MEMORY[ vr[i] ];
	}
	return %FMI_PREFIX%OK;
}
FMI_Dll_Export %FMI_PREFIX%Status %FMI_PREFIX%GetInteger(%FMI_PREFIX%Component c,
								const %FMI_PREFIX%ValueReference vr[],
								size_t nvr,
								%FMI_PREFIX%Integer value[])
{
	%VARPREFIX%ModelInstance* model_instance = (%VARPREFIX%ModelInstance*) c;

	/* 20-sim generated C-code uses doubles for the model equation calculations.
	   The variable type (i.e. integer or boolean) are not supported as data type,
	   but are transfered to doubles. In the FMI interface however, the double 
	   can be converted to its implicit type  */
	size_t i;
	for (i = 0; i < nvr; i++)
	{
		value [i] = (%FMI_PREFIX%Integer) model_instance->MEMORY[ vr[i] ];
	}
	return %FMI_PREFIX%OK;
}
FMI_Dll_Export %FMI_PREFIX%Status %FMI_PREFIX%GetBoolean(%FMI_PREFIX%Component c,
								const %FMI_PREFIX%ValueReference vr[],
								size_t nvr,
								%FMI_PREFIX%Boolean value[])
{
	%VARPREFIX%ModelInstance* model_instance = (%VARPREFIX%ModelInstance*) c;

	/* 20-sim generated C-code uses doubles for the model equation calculations.
	   The variable type (i.e. integer or boolean) are not supported as data type,
	   but are transfered to doubles. In the FMI interface however, the double
	   can be converted to its implicit type  */

	size_t i;
	for (i = 0; i < nvr; i++)
	{
		value [i] = (model_instance->MEMORY[ vr[i] ] == 1.0);
	}
	return %FMI_PREFIX%OK;
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
	%VARPREFIX%ModelInstance* model_instance = (%VARPREFIX%ModelInstance*) c;
	size_t i;
	for (i = 0; i < nvr; i++)
	{
		model_instance->MEMORY[ vr[i] ] = value [i];
	}
	return %FMI_PREFIX%OK;
}

FMI_Dll_Export %FMI_PREFIX%Status %FMI_PREFIX%SetInteger (%FMI_PREFIX%Component c,
								const %FMI_PREFIX%ValueReference vr[],
								size_t nvr,
								const %FMI_PREFIX%Integer value[])
{
	%VARPREFIX%ModelInstance* model_instance = (%VARPREFIX%ModelInstance*) c;

	/* 20-sim generated C-code uses doubles for the model equation calculations.
	   The variable type (i.e. integer or boolean) are not supported as data type,
	   but are transfered to doubles. In the FMI interface however, the double
	   can be converted to its implicit type  */
	size_t i;
	for (i = 0; i < nvr; i++)
	{
		model_instance->MEMORY[ vr[i] ] = (XXDouble) value [i];
	}
	return %FMI_PREFIX%OK;
}
FMI_Dll_Export %FMI_PREFIX%Status %FMI_PREFIX%SetBoolean(%FMI_PREFIX%Component c,
								const %FMI_PREFIX%ValueReference vr[],
								size_t nvr,
								const %FMI_PREFIX%Boolean value[])
{
	%VARPREFIX%ModelInstance* model_instance = (%VARPREFIX%ModelInstance*) c;
	/* temp implementation allowing boolean to real conversion, until proper boolean support is added. */
	size_t i;
	for (i = 0; i < nvr; i++)
	{
		model_instance->MEMORY[vr[i]] = value[i] ? 1.0 : 0.0;
	}
	return %FMI_PREFIX%OK;    
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
	%VARPREFIX%ModelInstance* model_instance;
	int offset = 0;
	
 	/* we should remember the functions pointer in order to make callback functions */
	if (!functions.logger) 
		return NULL; // we cannot even log this problem
	if (!instanceName || strlen(instanceName)==0)
	{
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
	if (!functions.allocateMemory || !functions.freeMemory){
		functions.logger(NULL, instanceName, fmiError, "error",
				"Missing memory callback function.");
		return NULL;
	}

	model_instance = (%VARPREFIX%ModelInstance *)functions.allocateMemory(1, sizeof(%VARPREFIX%ModelInstance));

	if(!model_instance)
	{
		functions.logger(NULL, instanceName, fmiError, "error",
			"Out of memory while allocating model instance");
		return NULL;
	}

	memset(model_instance, 0, sizeof(%VARPREFIX%ModelInstance));
	
	model_instance->instanceName = (%FMI_PREFIX%String) functions.allocateMemory(1 + strlen(instanceName), sizeof(char));
	
	if (!model_instance->instanceName)
	{
		functions.logger(NULL, instanceName, fmiError, "error",
			"fmiInstantiateSlave: Out of memory while allocating instance name");
		return NULL;
	}
	strcpy((char *)model_instance->instanceName, (char *)instanceName);

	/* Prepare the model instance struct */
	model_instance->start_time = %START_TIME%;
	model_instance->finish_time = 0.0;
	model_instance->m_use_finish_time = XXFALSE;
	model_instance->step_size = %TIME_STEP_SIZE%;
	model_instance->time = 0.0;
	model_instance->steps = 0;
	model_instance->%XX_INITIALIZE% = XXTRUE;
	model_instance->major = XXTRUE;
	model_instance->stop_simulation = XXFALSE;

	/* Set the offsets within the model_instance->MEMORY array */
%IF%%NUMBER_CONSTANTS%
	model_instance->%XX_CONSTANT_ARRAY_NAME% = &model_instance->MEMORY[offset]; /* constants offset */
	offset = offset + %VARPREFIX%constants_count;
%ENDIF%
%IF%%NUMBER_PARAMETERS%
	model_instance->%XX_PARAMETER_ARRAY_NAME% = &model_instance->MEMORY[offset];	/* parameters offset */
	offset = offset + %VARPREFIX%parameter_count;
%ENDIF%
%IF%%NUMBER_INITIAL_VALUES%
	model_instance->%XX_INITIAL_VALUE_ARRAY_NAME% = &model_instance->MEMORY[offset];		/* initial values offset */
	offset = offset + %VARPREFIX%initialvalue_count;
%ENDIF%
%IF%%NUMBER_VARIABLES%
	model_instance->%XX_VARIABLE_ARRAY_NAME% = &model_instance->MEMORY[offset];		/* variables offset */
	offset = offset + %VARPREFIX%variable_count;
%ENDIF%
%IF%%NUMBER_STATES%
	model_instance->%XX_STATE_ARRAY_NAME% = &model_instance->MEMORY[offset];		/* states offset */
	offset = offset + %VARPREFIX%state_count;
	model_instance->%XX_RATE_ARRAY_NAME% = &model_instance->MEMORY[offset];		/* rates offset */
	offset = offset + %VARPREFIX%state_count;
%ENDIF%

	/* Register the callback */
	model_instance->fmiCallbackFunctions = functions;

	/* Remember the resource folder location */
	model_instance->resourceLocation = URIToNativePath(model_instance, fmuLocation);

	/* Initialize our data arrays */
	%FUNCTIONPREFIX%ModelInitialize(model_instance);

	return (fmiComponent) model_instance;
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
	%VARPREFIX%ModelInstance* model_instance = NULL;
	int offset = 0;

	/* we should remember the functions pointer in order to make callback functions */
	if (!functions)
	{
		return NULL; // we cannot even log this problem
	}
	
	if (!functions->logger)
	{
		return NULL; // we cannot even log this problem
	}
	if (!instanceName || strlen(instanceName)==0)
	{
		functions->logger(NULL, "?", fmi2Error, "error",
				"Missing instance name.");
		return NULL;
	}
	if (!functions->allocateMemory || !functions->freeMemory) {
		functions->logger(functions->componentEnvironment, instanceName, fmi2Error, "error",
				"fmi2Instantiate: Missing memory callback function.");
		return NULL;
	}
	/* Check whether the given GUID equals our GUID */
	if (!fmuGUID || strlen(fmuGUID) == 0) {
		functions->logger(functions->componentEnvironment, instanceName, fmi2Error, "error",
				"fmi2Instantiate: Missing GUID.");
		return NULL;
	}
	if( strncmp(fmuGUID, FMI_GUID, strlen(fmuGUID)) != 0 )
	{
		functions->logger(NULL, instanceName, fmi2Error, "error",
			"fmi2Instantiate: Wrong GUID %s. Expected %s.", fmuGUID, FMI_GUID);
		return NULL;
	}
	
	model_instance = (%VARPREFIX%ModelInstance *)functions->allocateMemory(1, sizeof(%VARPREFIX%ModelInstance));

	if(!model_instance)
	{
		functions->logger(functions->componentEnvironment, instanceName, fmi2Error, "error",
			"fmi2Instantiate: Out of memory while allocating model instance");
		return NULL;
	}

	memset(model_instance, 0, sizeof(%VARPREFIX%ModelInstance));
	
	model_instance->instanceName = (%FMI_PREFIX%String) functions->allocateMemory(1 + strlen(instanceName), sizeof(char));

	if (!model_instance->instanceName)
	{
		functions->logger(functions->componentEnvironment, instanceName, fmi2Error, "error",
			"fmi2Instantiate: Out of memory while allocating instance name");
		return NULL;
	}
	strcpy((char *)model_instance->instanceName, (char *)instanceName);
	
	/* Prepare the model instance struct */
	model_instance->start_time = %START_TIME%;
	model_instance->finish_time = 0.0;
	model_instance->m_use_finish_time = XXFALSE;
	model_instance->step_size = %TIME_STEP_SIZE%;
	model_instance->time = 0.0;
	model_instance->steps = 0;
	model_instance->%XX_INITIALIZE% = XXTRUE;
	model_instance->major = XXTRUE;
	model_instance->stop_simulation = XXFALSE;

	/* Set the offsets within the model_instance->MEMORY array */
%IF%%NUMBER_CONSTANTS%
	model_instance->%XX_CONSTANT_ARRAY_NAME% = &model_instance->MEMORY[offset]; /* constants offset */
	offset = offset + %VARPREFIX%constants_count;
%ENDIF%
%IF%%NUMBER_PARAMETERS%
	model_instance->%XX_PARAMETER_ARRAY_NAME% = &model_instance->MEMORY[offset];	/* parameters offset */
	offset = offset + %VARPREFIX%parameter_count;
%ENDIF%
%IF%%NUMBER_INITIAL_VALUES%
	model_instance->%XX_INITIAL_VALUE_ARRAY_NAME% = &model_instance->MEMORY[offset];		/* initial values offset */
	offset = offset + %VARPREFIX%initialvalue_count;
%ENDIF%
%IF%%NUMBER_VARIABLES%
	model_instance->%XX_VARIABLE_ARRAY_NAME% = &model_instance->MEMORY[offset];		/* variables offset */
	offset = offset + %VARPREFIX%variable_count;
%ENDIF%
%IF%%NUMBER_STATES%
	model_instance->%XX_STATE_ARRAY_NAME% = &model_instance->MEMORY[offset];		/* states offset */
	offset = offset + %VARPREFIX%state_count;
	model_instance->%XX_RATE_ARRAY_NAME% = &model_instance->MEMORY[offset];		/* rates offset */
	offset = offset + %VARPREFIX%state_count;
%ENDIF%
	
	/* Register the callback */
	model_instance->fmiCallbackFunctions = functions;
	/* Remember the resource folder location */
	model_instance->resourceLocation = URIToNativePath(model_instance, fmuResourceLocation);
	
	/* check if we are setup for co-simulation, that's the only possible option for now */
	if( fmuType != fmi2CoSimulation )
	{
		model_instance->fmiCallbackFunctions->logger(NULL, instanceName, fmi2Error, "error",
			"FMU can only be used for Co-Simulation, not for Model Exchange");
		return NULL;
	}

	return (fmi2Component) model_instance;
}
%ENDIF%
%IF%%FMI1%
fmiStatus fmiInitializeSlave(fmiComponent c,
							 fmiReal tStart,
							 fmiBoolean StopTimeDefined,
							 fmiReal tStop) 
{
	%VARPREFIX%ModelInstance* model_instance = (%VARPREFIX%ModelInstance*) c;

	/* copy the arguments */
	model_instance->start_time = tStart;
	if (StopTimeDefined == fmiTrue)
	{
		model_instance->finish_time = tStop;
		model_instance->m_use_finish_time = XXTRUE;
	}

	/* initialize the submodel itself */
	%FUNCTIONPREFIX%InitializeSubmodel (model_instance);

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
	%VARPREFIX%ModelInstance* model_instance = (%VARPREFIX%ModelInstance*) c;

	/* copy the arguments */
	model_instance->start_time = startTime;
	if (stopTimeDefined == fmi2True)
	{
		model_instance->finish_time = stopTime;
		model_instance->m_use_finish_time = XXTRUE;
	}

	/* Initialize our data arrays */
	%FUNCTIONPREFIX%ModelInitialize(model_instance);

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
	%VARPREFIX%ModelInstance* model_instance = (%VARPREFIX%ModelInstance*) c;

	/* Initialize the submodel itself */
	%FUNCTIONPREFIX%InitializeSubmodel (model_instance);

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
	%VARPREFIX%ModelInstance* model_instance = (%VARPREFIX%ModelInstance*) c;

	/* Perform the final calculations */
	%FUNCTIONPREFIX%TerminateSubmodel (model_instance, model_instance->time);

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
	%VARPREFIX%ModelInstance* model_instance = (%VARPREFIX%ModelInstance*) c;

	/* initialize the submodel itself */
	%FUNCTIONPREFIX%InitializeSubmodel (model_instance);

	/* all done */
	return %FMI_PREFIX%OK;
}

%IF%%FMI1%
void fmiFreeSlaveInstance(fmiComponent c)
{
	%VARPREFIX%ModelInstance* model_instance = (%VARPREFIX%ModelInstance*) c;
	fmiCallbackFunctions fmiCallbackFunctions = model_instance->fmiCallbackFunctions;

	if (model_instance->resourceLocation != NULL)
	{
		model_instance->fmiCallbackFunctions.freeMemory((void*)model_instance->resourceLocation);
		model_instance->resourceLocation = NULL;
	}

	if(model_instance->instanceName != NULL)
	{
		model_instance->fmiCallbackFunctions.freeMemory((void *)model_instance->instanceName);
		model_instance->instanceName = NULL;
	}

	/* Copy the callback functions before freeing the model_instance */
	fmiCallbackFunctions.freeMemory((void *) model_instance);
	model_instance =  NULL;
}
%ENDIF%
%IF%%FMI2%
void fmi2FreeInstance(fmi2Component c)
{
	%VARPREFIX%ModelInstance* model_instance = (%VARPREFIX%ModelInstance*) c;
	const fmi2CallbackFunctions* fmiCallbackFunctions = model_instance->fmiCallbackFunctions;

	if (model_instance->resourceLocation != NULL)
	{
		model_instance->fmiCallbackFunctions->freeMemory((void*)model_instance->resourceLocation);
		model_instance->resourceLocation = NULL;
	}

	if(model_instance->instanceName != NULL)
	{
		model_instance->fmiCallbackFunctions->freeMemory((void *)model_instance->instanceName);
		model_instance->instanceName = NULL;
	}

	/* Copy the callback functions before freeing the model_instance */
	fmiCallbackFunctions->freeMemory((void *) model_instance);
	model_instance =  NULL;
}
%ENDIF%

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
	%VARPREFIX%ModelInstance* model_instance = (%VARPREFIX%ModelInstance*) c;

	/* Treat also case of zero step, i.e. during an event iteration */
	if (communicationStepSize == 0)
	{
		return %FMI_PREFIX%OK;
	}

	/* as long as we are not passed our communication point */
	while (model_instance->time < (currentCommunicationPoint + communicationStepSize))
	{
		/* check for termination first */
		if ( model_instance->m_use_finish_time && (model_instance->time > model_instance->finish_time) )
		{
%IF%%FMI2%
			if(model_instance->fmiCallbackFunctions != NULL && model_instance->fmiCallbackFunctions->logger != NULL)
			{
				model_instance->fmiCallbackFunctions->logger(NULL, "%SUBMODEL_NAME%", fmi2Error, "error",
					"Exceeded model finish time: %g > %g\n", model_instance->time, model_instance->finish_time);
			}
%ENDIF%
			
			/* we're done */
			return %FMI_PREFIX%Error;
		}
		/* Check for stop simulation */
		if (model_instance->stop_simulation == XXTRUE)
		{
			return %FMI_PREFIX%Error;
		}

		/* Call the submodel to calculate the output, and increase the time as well */
		%FUNCTIONPREFIX%CalculateSubmodel (model_instance, model_instance->time);
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
	/* all fine? */
	return %FMI_PREFIX%OK;
}

%FMI_PREFIX%Status %FMI_PREFIX%GetBooleanStatus(%FMI_PREFIX%Component c, const %FMI_PREFIX%StatusKind s, %FMI_PREFIX%Boolean* value)
{
	/* all fine? */
	return %FMI_PREFIX%OK;
}

%FMI_PREFIX%Status %FMI_PREFIX%GetStringStatus(%FMI_PREFIX%Component c, const %FMI_PREFIX%StatusKind s, %FMI_PREFIX%String*  value)
{
	/* not yet */
	return %FMI_PREFIX%Discard;
}

%IF%%FMI2%
fmi2Status fmi2GetContinuousStates(fmi2Component c, fmi2Real x[], size_t nx)
{
%IF%%NUMBER_STATES%
	size_t i;
	%VARPREFIX%ModelInstance* model_instance = (%VARPREFIX%ModelInstance*) c;

%ENDIF%
	if( %NUMBER_STATES% != nx )
	{
		return fmi2Error;
	}
%IF%%NUMBER_STATES%

	for( i = 0; i < %NUMBER_STATES%; ++i)
	{
		x[i] = model_instance->%XX_STATE_ARRAY_NAME%[i];
	}
%ELSE%
	/* the exported submodel has no states */
%ENDIF%
	return fmi2OK;
}
%ENDIF%

%IF%%FMI2%
fmi2Status fmi2GetFMUstate (fmi2Component c, fmi2FMUstate* FMUstate)
{
	/* not yet */
	return fmi2Discard;
}

fmi2Status fmi2SetFMUstate (fmi2Component c, fmi2FMUstate FMUstate)
{
	/* not yet */
	return fmi2Discard;
}

fmi2Status fmi2FreeFMUstate(fmi2Component c, fmi2FMUstate* FMUstate)
{
	/* not yet */
	return fmi2Discard;
}

fmi2Status fmi2SerializedFMUstateSize(fmi2Component c, fmi2FMUstate FMUstate, size_t *size)
{
	/* not yet */
	return fmi2Discard;
}

fmi2Status fmi2SerializeFMUstate (fmi2Component c, fmi2FMUstate FMUstate, fmi2Byte serializedState[], size_t size)
{
	/* not yet */
	return fmi2Discard;
}

fmi2Status fmi2DeSerializeFMUstate (fmi2Component c, const fmi2Byte serializedState[], size_t size,
                                    fmi2FMUstate* FMUstate)
{
	/* not yet */
	return fmi2Discard;
}

fmi2Status fmi2GetDirectionalDerivative(fmi2Component c, const fmi2ValueReference vUnknown_ref[], size_t nUnknown,
                                        const fmi2ValueReference vKnown_ref[] , size_t nKnown,
                                        const fmi2Real dvKnown[], fmi2Real dvUnknown[])
{
	/* not yet */
	return fmi2Discard;
}
%ENDIF%