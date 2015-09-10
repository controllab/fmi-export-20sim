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

/* This file contains support functions for several SIDOPS functions

   For flexibility, ANSI-C is created, and typedefs are used
   for integers and doubles, see the xxfuncs.h file for more
   information on these types.

   This means that all used functions follow the ANSI definition

   Please check the math.h file of your particular compiler
   to see if this is indeed the case. Otherwise, you might have
   to adapt the used functions below to obtain the same behavior.

*/

/* The system include files */
#include <stdlib.h>
#include <math.h>
%IF%%NUMBEROF_DELAYFUNCTION%
/* Include the header for memcpy
 * You may need to change this into <memory.h> for older compilers
 */
#include <string.h>
%ENDIF%

/* Our own include files */
#include "xxfuncs.h"

/* Constants that are used in our functions below */
XXDouble xx_logarithm_2 =  0.6931471805599;
XXDouble xx_logarithm_10 = 2.3025850929940;

/* The 20-sim SIDOPS support functions */
XXDouble XXAbsolute (XXDouble argument)
{
	return (XXDouble) fabs (argument);
}

XXDouble XXArcCosineHyperbolic (XXDouble argument)
{
	return (XXDouble) log (argument + sqrt(argument * argument - 1.0));
}

XXDouble XXArcSineHyperbolic (XXDouble argument)
{
	return (XXDouble) log (argument + sqrt(argument * argument + 1.0));
}

XXDouble XXArcTangentHyperbolic (XXDouble argument)
{
	return (XXDouble) 0.5 * log ((1.0 + argument) / (1.0 - argument));;
}

XXDouble XXExponent2 (XXDouble argument)
{
	return (XXDouble) exp (argument * xx_logarithm_2);
}

XXDouble XXExponent10 (XXDouble argument)
{
	return (XXDouble) exp (argument * xx_logarithm_10);
}

XXDouble XXIntegerDivide (XXDouble argument1, XXDouble argument2)
{
	XXInteger value;

	value = (XXInteger) (argument1 / argument2);
	return (XXDouble) value;
}

XXDouble XXIntegerModulo (XXDouble argument1, XXDouble argument2)
{
	XXInteger value;

	value = (XXInteger) (argument1 / argument2);
	return (XXDouble) argument1 - (value * argument2);
}

XXDouble XXLogarithm2 (XXDouble argument)
{
	return (XXDouble) log (argument) / xx_logarithm_2;
}

XXDouble XXLogarithm10 (XXDouble argument)
{
	return (XXDouble) log (argument) / xx_logarithm_10;
}

XXDouble XXPow2 (XXDouble argument)
{
	return argument * argument;
}

XXDouble XXPower (XXDouble argument1, XXDouble argument2)
{
	return (XXDouble) pow (argument1, argument2);
}

XXDouble XXRandom (XXDouble argument)
{
	XXDouble value;

	value = (XXDouble) rand() / (XXDouble) RAND_MAX - 0.5;
	return argument * 2.0 * value;
}

XXDouble XXRamp (XXDouble argument, XXDouble time)
{
	XXDouble value;

	if (time < argument)
		value = 0.0;
	else
		value = time - argument;
	return value;
}

XXDouble XXSign (XXDouble argument)
{
	XXDouble value;
	if (argument < 0.0)
		value = -1.0;
	else
		if (argument == 0.0)
			value = 0.0;
		else
			value = 1.0;
	return value;
}

XXDouble XXStep (XXDouble argument, XXDouble time)
{
	XXDouble value;

	if (time < argument)
		value = 0.0;
	else
		value = 1.0;
	return value;
}

XXDouble XXImpulse (XXDouble impulsestarttime, XXDouble impulseduration, XXDouble currenttime, XXDouble stepsize)
{
	XXDouble value;

	if (stepsize <= 0.0 || impulseduration <= 0.0)
		value = 0.0;
	else
	{
		if ((currenttime < impulsestarttime) || (currenttime > (impulsestarttime + impulseduration)))
			value = 0.0;
		else
		{
			if (stepsize < impulseduration)
				value = (1.0 / impulseduration);
			else
				value = (1.0 / stepsize);
		}
	}
	return value;
}

XXDouble XXXor(XXDouble argument1, XXDouble argument2)
{
	return (argument1 || argument2) && !(argument1 && argument2);
}

%IF%%NUMBEROF_DELAYFUNCTION%
XXDouble %VARPREFIX%delay_update_array[%NUMBEROF_DELAYFUNCTION%];
XXDouble %VARPREFIX%delay_last_values[%NUMBEROF_DELAYFUNCTION%];
void XXDelayUpdate()
{
	memcpy(%VARPREFIX%delay_update_array, %VARPREFIX%delay_last_values, %NUMBEROF_DELAYFUNCTION% * sizeof(XXDouble));
}
XXDouble XXDelay (XXDouble argument1, XXDouble argument2, XXInteger id)
{
	XXDouble value;

	if (%VARPREFIX%%XX_INITIALIZE%)
	{
		value = argument2;
	}
	else
	{
		value = %VARPREFIX%delay_update_array[id];
	}

	if (%VARPREFIX%major)
	{
		%VARPREFIX%delay_last_values[id] = argument1;
	}

	return value;
}
%ENDIF%


XXDouble XXRound (XXDouble argument)
{
	XXDouble leftOver, result;

	leftOver = argument - (XXInteger) argument;
	if (fabs (leftOver) < 0.5)
	{
		result = (XXDouble) ((XXInteger) argument);
	}
	else
	{
		if (argument >= 0)
			result = (XXDouble) ceil (argument);
		else
		{
			result = (XXDouble) floor (argument);
		}
	}
	return result;
}

%IF%%NUMBEROF_INITIALFUNCTION%
XXDouble %VARPREFIX%initial_value_array[%NUMBEROF_INITIALFUNCTION%];

XXDouble XXInitialValue (XXDouble argument, XXInteger identifier)
{
	XXDouble value;

	if (%VARPREFIX%%XX_INITIALIZE%)
	{
		value = argument;
		%VARPREFIX%initial_value_array[identifier] = value;
	}
	else
	{
		value = %VARPREFIX%initial_value_array[identifier];
	}
	return value;
}
%ENDIF%

XXInteger XXBitAnd(XXInteger argument1, XXInteger argument2)
{
	/* bitwise and */
	return (argument1 & argument2);
}

XXInteger XXBitOr(XXInteger argument1, XXInteger argument2)
{
	/* bitwise or */
	return  (argument1 | argument2);
}

XXInteger XXBitXor(XXInteger argument1, XXInteger argument2)
{
	/* bitwise xor */
	return (argument1 ^ argument2);
}

XXInteger XXBitCmp(XXInteger argument, XXInteger nrBits)
{
	XXInteger bits = (XXInteger)(pow(2.0, (double)nrBits)) - 1;

	/* only do the last number of asked bits (the and operator) and invert */
	return (bits - argument);
}

XXInteger XXBitGet(XXInteger argument, XXInteger bitPos)
{
	return (argument & (1 << (bitPos - 1)));
}

XXInteger XXBitSet(XXInteger argument, XXInteger bitPos)
{
	/* Set the bit to 1 */
	return (argument | (1 << (bitPos - 1)));
}

XXInteger XXBitClear(XXInteger argument, XXInteger bitPos)
{
	/* Clear the bit */
	return (argument & ~(1 << (bitPos - 1)));
}

XXInteger XXBitShift(XXInteger argument, XXInteger bitsToShift)
{
	if ( bitsToShift > 0 )
	{
		return (argument << bitsToShift);
	}
	else
	{
		return (argument >> bitsToShift);
	}
}

XXInteger XXSwapBytes(XXInteger argument)
{
	int arg1;
	int byte1;
	int byte2;
	int byte3;
	int byte4;
	int result;
	
	arg1 = (int) argument;	/* the argument to swap (only 32-bits swap) */

	byte1 =  0x000000ff & arg1;
	byte2 = (0x0000ff00 & arg1) >> 8;
	byte3 = (0x00ff0000 & arg1) >> 16;
	byte4 = (0xff000000 & arg1) >> 24;

	/* do the explicit 32-bit swap */
	result = (byte1 << 24) + (byte2 << 16) + (byte3 << 8) + byte4;

	return (XXInteger) result;
}

/* 20-sim stubs. Implement them yourself if needed */
XXDouble XXData (XXString name, XXInteger column, XXInteger id)
{
	return 0;
}

XXDouble XXTable (XXString name, XXDouble argument, XXInteger id)
{
	return 0;
}

XXBoolean XXEvent (XXDouble argument, XXInteger id)
{
	return 0;
}

XXBoolean XXEventUp (XXDouble argument, XXInteger id)
{
	return 0;
}

XXBoolean XXEventDown (XXDouble argument, XXInteger id)
{
	return 0;
}

XXBoolean XXFrequencyEvent (XXDouble argument, XXInteger id)
{
	return 0;
}

XXBoolean XXFrequencyEvent1 (XXDouble argument1, XXDouble argument2, XXInteger id)
{
	return 0;
}

XXBoolean XXTimeEvent (XXDouble argument, XXInteger id)
{
	return 0;
}

XXDouble XXTimeDelay (XXDouble argument, XXDouble time, XXInteger id)
{
	return 0;
}

XXBoolean XXWarning (XXString message, XXInteger id)
{
	return 0;
}

XXBoolean XXStopSimulation (XXString message, XXInteger id)
{
	%VARPREFIX%stop_simulation = XXTRUE;
	return 0;
}
