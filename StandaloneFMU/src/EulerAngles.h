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

/*****************************************************************************
--------------------------------------------------------------------------

 Copyright � 1995-2007 Controllab Products B.V. All Rights Reserved.

             www.controllab.nl
               www.20sim.com

 Authors: Frank Groen
--------------------------------------------------------------------------
 Description:

 Contains code for converting between different angle representations
 all 24 different kind of Euler angles
 Quaternions and H-Matrices and Rotation Matrices
--------------------------------------------------------------------------
******************************************************************************/

%IF%%NUMBEROF_DLL_EulerAngles%

/**** EulerAngles.h - Support for 24 angle schemes ****/

#ifndef _H_EulerAngles
#define _H_EulerAngles

/* 20-sim include files */
#include "xxtypes.h"

/*** Definitions ***/
typedef struct
{
	XXDouble x,
		   y,
		   z;
	XXDouble w;
} Quat; /* Quaternion */

/* separate typedef for EulerAngles */
typedef struct
{
	XXDouble x;
	XXDouble y;
	XXDouble z;
	XXInteger order;
} EulerAngles;

enum QuatPart {
	Q_X = 0,
	Q_Y,
	Q_Z,
	Q_W
};

//typedef XXDouble HMatrix[4][4]; /* Right-handed, for column vectors */
typedef XXDouble MyHMatrix[16]; /* now just linear space memory */
/* and a macro to calculate the index */
#define INDEX4x4(i,j) ((j<<2)+i)


/*** Order type constants, constructors, extractors ***/
    /* There are 24 possible conventions, designated by:    */
    /*	  o EulAxI = axis used initially		    */
    /*	  o EulPar = parity of axis permutation		    */
    /*	  o EulRep = repetition of initial axis as last	    */
    /*	  o EulFrm = frame from which axes are taken	    */
    /* Axes I,J,K will be a permutation of X,Y,Z.	    */
    /* Axis H will be either I or K, depending on EulRep.   */
    /* Frame S takes axes from initial static frame.	    */
    /* If ord = (AxI=X, Par=Even, Rep=No, Frm=S), then	    */
    /* {a,b,c,ord} means Rz(c)Ry(b)Rx(a), where Rz(c)v	    */
    /* rotates v around Z by c radians.			    */
#define EulFrmS	     0
#define EulFrmR	     1
#define EulRepNo     0
#define EulRepYes    1
#define EulParEven   0
#define EulParOdd    1

/* these are precalculated values (so no overhead) */
	/* Static axes */
#define EulOrdXYZs    0
#define EulOrdXYXs    2
#define EulOrdXZYs    4
#define EulOrdXZXs    6
#define EulOrdYZXs    8
#define EulOrdYZYs    10
#define EulOrdYXZs    12
#define EulOrdYXYs    14
#define EulOrdZXYs    16
#define EulOrdZXZs    18
#define EulOrdZYXs    20
#define EulOrdZYZs    22
    /* Rotating axes */
#define EulOrdZYXr    1
#define EulOrdXYXr    3
#define EulOrdYZXr    5
#define EulOrdXZXr    7
#define EulOrdXZYr    9
#define EulOrdYZYr    11
#define EulOrdZXYr    13
#define EulOrdYXYr    15
#define EulOrdYXZr    17
#define EulOrdZXZr    19
#define EulOrdXYZr    21
#define EulOrdZYZr    23

/* initialize and terminate function, which does nothing */
void EulerAngles_Initialize(void);
void EulerAngles_Terminate(void);

/*
void Eul_ToQuat(EulerAngles *ea, Quat *qu);
void Eul_ToHMatrix(EulerAngles *ea, MyHMatrix *M);
void Eul_FromHMatrix(MyHMatrix *M, XXInteger order, EulerAngles *ea);
void Eul_FromQuat(Quat *q, XXInteger order, EulerAngles *ea);
void Quat_ToHMatrix(Quat *q, MyHMatrix *M);
void Quat_FromHMatrix(MyHMatrix *M, Quat *q);
*/

/* relative from MyHMatrix */
XXInteger EulerAngles_EulXYXrFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulXYZrFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulXZXrFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulXZYrFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulYXYrFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulYXZrFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulYZXrFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulYZYrFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulZXYrFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulZXZrFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulZYXrFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulZYZrFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);

/* relative from RotationMatrix */
XXInteger EulerAngles_EulXYXrFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulXYZrFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulXZXrFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulXZYrFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulYXYrFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulYXZrFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulYZXrFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulYZYrFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulZXYrFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulZXZrFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulZYXrFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulZYZrFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);


/* static from MyHMatrix */
XXInteger EulerAngles_EulXYXsFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulXYZsFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulXZXsFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulXZYsFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulYXYsFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulYXZsFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulYZXsFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulYZYsFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulZXYsFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulZXZsFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulZYXsFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulZYZsFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);

/* static from RotationMatrix */
XXInteger EulerAngles_EulXYXsFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulXYZsFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulXZXsFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulXZYsFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulYXYsFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulYXZsFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulYZXsFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulYZYsFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulZXYsFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulZXZsFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulZYXsFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_EulZYZsFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);

/************************************************************************/
/* the other way round                                                  */
/************************************************************************/

/* relative from MyHMatrix */
XXInteger EulerAngles_HMatrixFromEulXYXr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulXYZr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulXZXr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulXZYr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulYXYr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulYXZr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulYZXr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulYZYr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulZXYr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulZXZr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulZYXr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulZYZr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);

/* static from MyHMatrix */
XXInteger EulerAngles_HMatrixFromEulXYXs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulXYZs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulXZXs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulXZYs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulYXYs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulYXZs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulYZXs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulYZYs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulZXYs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulZXZs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulZYXs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_HMatrixFromEulZYZs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);

/* relative from RotationMatrix */
XXInteger EulerAngles_RotationMatrixFromEulXYXr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulXYZr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulXZXr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulXZYr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulYXYr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulYXZr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulYZXr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulYZYr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulZXYr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulZXZr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulZYXr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulZYZr (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);

/* static from RotationMatrix */
XXInteger EulerAngles_RotationMatrixFromEulXYXs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulXYZs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulXZXs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulXZYs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulYXYs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulYXZs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulYZXs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulYZYs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulZXYs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulZXZs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulZYXs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);
XXInteger EulerAngles_RotationMatrixFromEulZYZs (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);

/* Quaternion to MyHMatrix */
XXInteger EulerAngles_HMatrixFromQuaternion (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);

/* MyHMatrix to Quaternion */
XXInteger EulerAngles_QuaternionFromHMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);

/* Quaternion to RotationMatrix */
XXInteger EulerAngles_RotationMatrixFromQuaternion (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);

/* RotationMatrix to Quaternion */
XXInteger EulerAngles_QuaternionFromRotationMatrix (XXDouble *inarr, XXInteger inputs, XXDouble *outarr, XXInteger outputs, XXInteger major);

#endif
%ENDIF%