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

   This means that all used functions follow the ANSI definition.

   Please check the math.h file of your particular compiler
   to see if this is indeed the case. Otherwise, you might have
   to adapt the used functions below to obtain the same behavior.

*/

/* The system include files */
#include <stdlib.h>
#include <math.h>
%IF%%NUMBEROF_REALTIME%
#include <time.h>
%ENDIF%
/* Our own include files */
#include "xxfuncs.h"

/* Constants that are used in our functions below */
%IF%%OR(NUMBEROF_LOG2FUNCTION, NUMBEROF_EXP2FUNCTION)%
static const XXDouble xx_logarithm_2 =  0.6931471805599453;
%IF%%NUMBEROF_LOG2FUNCTION%
static const XXDouble xx_invLog2 = 1.4426950408889634;
%ENDIF%
%ENDIF%
%IF%%OR(NUMBEROF_LOG10FUNCTION, NUMBEROF_EXP10FUNCTION)%
static const XXDouble xx_logarithm_10 = 2.3025850929940457;
%IF%%NUMBEROF_LOG10FUNCTION%
#if __STDC_VERSION__ < 199901L
static const XXDouble xx_invLog10 = 0.4342944819032518;
#endif
%ENDIF%
%ENDIF%

typedef union
{
	double m_double;
	const char* m_char;
}str2dbl;

XXDouble XXString2Double(const char* argument)
{
	str2dbl myConversion;
	myConversion.m_char = argument;
	return myConversion.m_double;

}

const char* XXDouble2String(XXDouble argument)
{
	str2dbl myConversion;
	myConversion.m_double = argument;
	return myConversion.m_char;
}

/* The 20-sim SIDOPS functions */
%IF%%NUMBEROF_ABSFUNCTION%
XXDouble XXAbsolute (XXDouble argument)
{
	return (XXDouble) fabs (argument);
}

%ENDIF%
%IF%%NUMBEROF_ARCCOSHYPERBOLICFUNCTION%
XXDouble XXArcCosineHyperbolic (XXDouble argument)
{
	return (XXDouble) log (argument + sqrt(argument * argument - 1.0));
}

%ENDIF%
%IF%%NUMBEROF_ARCSINHYPERBOLICFUNCTION%
XXDouble XXArcSineHyperbolic (XXDouble argument)
{
	return (XXDouble) log (argument + sqrt(argument * argument + 1.0));
}

%ENDIF%
%IF%%NUMBEROF_ARCTANHYPERBOLICFUNCTION%
XXDouble XXArcTangentHyperbolic (XXDouble argument)
{
	return (XXDouble) 0.5 * log ((1.0 + argument) / (1.0 - argument));
}

%ENDIF%
%IF%%NUMBEROF_EXP2FUNCTION%
XXDouble XXExponent2 (XXDouble argument)
{
	return (XXDouble) exp (argument * xx_logarithm_2);
}

%ENDIF%
%IF%%NUMBEROF_EXP10FUNCTION%
XXDouble XXExponent10 (XXDouble argument)
{
	return (XXDouble) exp (argument * xx_logarithm_10);
}

%ENDIF%
%IF%%NUMBEROF_DIVISION%
XXDouble XXIntegerDivide (XXDouble argument1, XXDouble argument2)
{
	XXInteger value;

	value = (XXInteger) (argument1 / argument2);
	return (XXDouble) value;
}

%ENDIF%
%IF%%NUMBEROF_MODULO%
XXDouble XXIntegerModulo (XXDouble argument1, XXDouble argument2)
{
	XXInteger value;

	value = (XXInteger) (argument1 / argument2);
	return (XXDouble) argument1 - (value * argument2);
}

%ENDIF%
%IF%%NUMBEROF_LOG2FUNCTION%
XXDouble XXLogarithm2 (XXDouble argument)
{
	return (XXDouble) (log (argument) * xx_invLog2);
}

%ENDIF%
%IF%%NUMBEROF_LOG10FUNCTION%
XXDouble XXLogarithm10 (XXDouble argument)
{
#if (__STDC_VERSION__ >= 199901L) || (__cplusplus)
	/* C99 / C++ */
	return (XXDouble) log10 (argument);
#else
	/* Not C99 */
	return (XXDouble) (log (argument) * xx_invLog10);
#endif
}

%ENDIF%
%IF%%NUMBEROF_POWER%
XXDouble XXPow2 (XXDouble argument)
{
	return argument * argument;
}

XXDouble XXPower (XXDouble argument1, XXDouble argument2)
{
	return (XXDouble) pow (argument1, argument2);
}

%ENDIF%
%IF%%OR(NUMBEROF_GAUSSFUNCTION,NUMBEROF_RANDOMFUNCTION)%
XXDouble XXRandom (XXDouble argument)
{
	XXDouble value;

	value = (XXDouble) rand() / (XXDouble) RAND_MAX - 0.5;
	return argument * 2.0 * value;
}

%ENDIF%
%IF%%NUMBEROF_RAMPFUNCTION%
XXDouble XXRamp (XXDouble argument, XXDouble time)
{
	XXDouble value;

	if (time < argument)
		value = 0.0;
	else
		value = time - argument;
	return value;
}

%ENDIF%
%IF%%NUMBEROF_SIGNFUNCTION%
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

%ENDIF%
%IF%%NUMBEROF_STEPFUNCTION%
XXDouble XXStep (XXDouble steptime, XXDouble time)
{
	XXDouble value;

	if (time < steptime)
		value = 0.0;
	else
		value = 1.0;
	return value;
}

%ENDIF%
%IF%%NUMBEROF_IMPULSEFUNCTION%
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

%ENDIF%
%IF%%NUMBEROF_XOR%
XXDouble XXXor(XXDouble argument1, XXDouble argument2)
{
	return (argument1 || argument2) && !(argument1 && argument2);
}

%ENDIF%
%IF%%NUMBEROF_ROUNDFUNCTION%
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

%ENDIF%
%IF%%NUMBEROF_BITAND%
XXInteger XXBitAnd(XXInteger argument1, XXInteger argument2)
{
	/* bitwise and */
	return (argument1 & argument2);
}
%ENDIF%
%IF%%NUMBEROF_BITOR%
XXInteger XXBitOr(XXInteger argument1, XXInteger argument2)
{
	/* bitwise or */
	return  (argument1 | argument2);
}

%ENDIF%
%IF%%NUMBEROF_BITXOR%
XXInteger XXBitXor(XXInteger argument1, XXInteger argument2)
{
	/* bitwise xor */
	return (argument1 ^ argument2);
}

%ENDIF%
%IF%%NUMBEROF_BITCMP%
XXInteger XXBitCmp(XXInteger argument, XXInteger nrBits)
{
	XXInteger maxBits = (XXInteger) sizeof(XXInteger) << 3;

	/* calculate the maximum unsigned value for nrBits */
	if (nrBits < maxBits)
	{
		XXInteger bits = (XXInteger) (1 << nrBits) - 1;
		/* invert and only return the number of asked bits */
		return (~argument & bits);
	}

	return(~argument);
}

%ENDIF%
%IF%%NUMBEROF_BITGET%
XXInteger XXBitGet(XXInteger argument, XXInteger bitPos)
{
	/* get the bit itself (prevent double shifting) */
	return ((argument >> (bitPos - 1)) & 1);
}

%ENDIF%
%IF%%NUMBEROF_BITINV%
XXInteger XXBitInv(XXInteger argument)
{
	return ~argument;
}

%ENDIF%
%IF%%NUMBEROF_BITSET%
XXInteger XXBitSet(XXInteger argument, XXInteger bitPos)
{
	/* set the bit to 1 */
	return (argument | (1 << (bitPos - 1)));
}

%ENDIF%
%IF%%NUMBEROF_BITCLEAR%
XXInteger XXBitClear(XXInteger argument, XXInteger bitPos)
{
	/* reset the bit to 0 */
	return (argument & ~(1 << (bitPos - 1)));
}

%ENDIF%
%IF%%NUMBEROF_BITSHIFT%
XXInteger XXBitShift(XXInteger argument, XXInteger bitsToShift)
{
	if ( bitsToShift > 0 )
	{
		return (argument << bitsToShift);
	}
	else
	{
		return (argument >> (-bitsToShift));
	}
}
%ENDIF%
%IF%%NUMBEROF_BITSHIFTRIGHT%
XXInteger XXBitShiftRight(XXInteger argument, XXInteger bitsToShift)
{
	if ( bitsToShift > 0 )
	{
		return (argument >> bitsToShift);
	}
	else
	{
		return (argument << (-bitsToShift));
	}
}

%ENDIF%
%IF%%NUMBEROF_BITSWAPBYTES%
XXInteger XXSwapBytes(XXInteger argument)
{
	/* this function swaps the 4 bytes of a 32-bit integer (little-big endian) */
	XXCharacter byte1;
	XXCharacter byte2;
	XXCharacter byte3;
	XXCharacter byte4;
	int result;
	int arg1 = (int) argument;

	/* get the separate bytes */
	byte1 = (XXCharacter)(arg1 & 0xFF);
	byte2 = (XXCharacter)((arg1 >> 8) & 0xFF);
	byte3 = (XXCharacter)((arg1 >> 16) & 0xFF);
	byte4 = (XXCharacter)((arg1 >> 24) & 0xFF);

	/* do the explicit 32-bit swap */
	result = (byte1 << 24) | (byte2 << 16) | (byte3 << 8) | byte4;

	return (XXInteger) result;
}

%ENDIF%
/* 20-sim stubs. Implement them yourself if needed */
%IF%%NUMBEROF_DATAFUNCTION%
XXDouble XXData (XXString name, XXInteger column, XXInteger id)
{
#if defined _MSC_VER
#pragma message("warning: The 20-sim 'data' function is not yet implemented in this code generation template")
#elif defined __GNUC__
#warning The 20-sim 'data' function is not yet implemented in this code generation template
#endif
	return 0;
}

%ENDIF%
%IF%%NUMBEROF_TABLEFUNCTION%
XXDouble XXTable (XXString name, XXDouble argument, XXInteger id)
{
#if defined _MSC_VER
#pragma message("warning: The 20-sim 'table' function is not yet implemented in this code generation template")
#elif defined __GNUC__
#warning The 20-sim 'table' function is not yet implemented in this code generation template
#endif
	return 0;
}

%ENDIF%
%IF%%NUMBEROF_FREQUENCYEVENTFUNCTION%
#if defined _MSC_VER
#pragma message("warning: The 20-sim 'frequencyevent' function is not yet implemented in this code generation template")
#elif defined __GNUC__
#warning The 20-sim 'frequencyevent' function is not yet implemented in this code generation template
#endif
XXBoolean XXFrequencyEvent (XXDouble argument, XXInteger id)
{
	return 0;
}

XXBoolean XXFrequencyEvent1 (XXDouble argument1, XXDouble argument2, XXInteger id)
{
	return 0;
}

%ENDIF%
%IF%%NUMBEROF_TDELAYFUNCTION%
XXDouble XXTimeDelay (XXDouble argument, XXDouble time, XXInteger id)
{
#if defined _MSC_VER
#pragma message("warning: The 20-sim 'tdelay' function is not yet implemented in this code generation template")
#elif defined __GNUC__
#warning The 20-sim 'tdelay' function is not yet implemented in this code generation template
#endif
	return 0;
}

%ENDIF%
%IF%%NUMBEROF_REALTIME%
static time_t xx_start_run_time = 0;

/* Return the elapsed amount of seconds since the start of this program
 * This reference implementation has an accuracy of 1 {s}
 */
XXDouble XXRealTime(void)
{
	XXDouble seconds = 0.0;
	
	if (xx_start_run_time == 0)
	{
		time(&xx_start_run_time);
	}
	else
	{
		seconds = (XXDouble) difftime(time(NULL), xx_start_run_time);
	}
	
	return seconds;
}

%ENDIF%
