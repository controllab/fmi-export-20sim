@ECHO OFF
rem ----Usage----
rem build [clean|noclean]
rem vs2010 for compiling with visual studio 2010
rem clean to force a full rebuild
rem noclean to force a build without clean
rem noprompt to avoid all prompts
CLS
COLOR 1B
TITLE 20-sim FMU %FMUVERSION% export - Build script
rem -------------------------------------------------------------
rem Config
rem If you get an error that Visual studio was not found, SET your path for VSNET main executable.
rem -------------------------------------------------------------
rem	CONFIG START
SET CURPATH=%~dp0
SET comp=vs2015
SET promptlevel=prompt
SET exitcode=0
SET buildmode=clean

cd ..\
SET ROOTPATH=%CD%
cd %CURPATH%

SET FMU=%ROOTPATH%\%SUBMODEL_NAME%.fmu
SET DLL=%SUBMODEL_NAME%.dll
SET ZIPTOOL=%TEMPLATE_DIR%\bin\7z.exe
rem The Github hosted template includes these two tools in order to support older 20-sim versions
SET XSLTTOOL=%TEMPLATE_DIR%\bin\msxsl.exe
SET GUIDTOOL=%TEMPLATE_DIR%\bin\GenerateGuid.exe

SET PYTHON="C:\ProgramData\Controllab Products B.V\Python34\python.exe"
SET RESOURCES_SCRIPT="%TEMPLATE_DIR%\bin\include_resources.py"

set FMU_DIR=%CURPATH%fmu
set BIN32_DIR=%FMU_DIR%\binaries\win32
set BIN64_DIR=%FMU_DIR%\binaries\win64
set SRC_DIR=%FMU_DIR%\sources
set DOC_DIR=%FMU_DIR%\documentation
set RES_DIR=%FMU_DIR%\resources

ECHO ------------------------------------------------------------
ECHO 20-sim standalone co-simulation FMU export for '%SUBMODEL_NAME%'
ECHO ------------------------------------------------------------
ECHO Creating an empty FMU
if not exist "%FMU_DIR%" mkdir "%FMU_DIR%"
if not exist "%SRC_DIR%" mkdir "%SRC_DIR%"
if not exist "%DOC_DIR%" mkdir "%DOC_DIR%"
if not exist "%RES_DIR%" mkdir "%RES_DIR%"

REM Generate a new GUID in XML and C format
REM -------------------
ECHO |SET /p=Generating a GUID: 
"%GUIDTOOL%" -c "%ROOTPATH%\src\fmiGUID.h" -x "%ROOTPATH%\src\GUID.xml"
ECHO.

ECHO Generating the modelDescription.xml
"%XSLTTOOL%"  "%ROOTPATH%\src\ModelConfiguration.xml" "%ROOTPATH%\template\mcf2modelDescription.xsl" -o "%FMU_DIR%\modelDescription.xml" SOURCEDIRECTORY="%ROOTPATH%\src"

ECHO Collecting resources
IF EXIST %PYTHON% (
	%PYTHON% %RESOURCES_SCRIPT% "%ROOTPATH%\src\ModelConfiguration.xml" "%RES_DIR%"
) ELSE (
	ECHO Unable to collect resources. Could not find Python in %PYTHON%
	ECHO Please re-install 20-sim 4.6 with Python support enabled.
)

ECHO ------------------------------------------------------------
ECHO Searching for Visual C++ compiler...

FOR %%b in (%1, %2, %3, %4, %5) DO (
	IF %%b==vs2010 SET comp=vs2010
	IF %%b==vs2013 SET comp=vs2013
	IF %%b==vs2015 SET comp=vs2015
	IF %%b==clean SET buildmode=clean
	IF %%b==noclean SET buildmode=noclean
	IF %%b==noprompt SET promptlevel=noprompt
)

SET buildconfig=Release
SET DEVENV=""
SET VSVARS32=""
SET BUILD_X64=1

rem Seach for VS 2015 / VS 2015 Express
IF %comp%==vs2015 (
	set PROJ_DIR=VS2015
	IF EXIST "%VS140COMNTOOLS%\vsvars32.bat" (
		set VSVARS32="%VS140COMNTOOLS%\vsvars32.bat"
		ECHO Found Visual C++ 2015
	) ELSE IF EXIST "%ProgramFiles%\Microsoft Visual Studio 14.0\Common7\Tools\vsvars32.bat" (
		set VSVARS32="%ProgramFiles%\Microsoft Visual Studio 14.0\Common7\Tools\vsvars32.bat"
		ECHO Found Visual C++ 2015
	) ELSE (
		rem Try an older compiler
		set comp=vs2013
	)
)

rem Seach for VS 2013 / VS 2013 Express / VS 2013 Community edition
IF %comp%==vs2013 (
	set PROJ_DIR=VS2013
	IF EXIST "%VS120COMNTOOLS%\vsvars32.bat" (
		set VSVARS32="%VS120COMNTOOLS%\vsvars32.bat"
		ECHO Found Visual C++ 2013
	) ELSE IF EXIST "%ProgramFiles%\Microsoft Visual Studio 12.0\Common7\Tools\vsvars32.bat" (
		set VSVARS32="%ProgramFiles%\Microsoft Visual Studio 12.0\Common7\Tools\vsvars32.bat"
		ECHO Found Visual C++ 2013
	) ELSE (
		rem Try an older compiler
		set comp=vs2010
	)
)

rem Seach for VS 2010 / VS 2010 Express
IF %comp%==vs2010 (
	set PROJ_DIR=VS2010
	IF EXIST "%VS100COMNTOOLS%\vsvars32.bat" (
		set VSVARS32="%VS100COMNTOOLS%\vsvars32.bat"
		IF EXIST "%VS100COMNTOOLS%\..\IDE\devenv.exe" (
			ECHO Found Visual C++ 2010
		) ELSE (
			REM The VC2010 express edition does not support x64 compilation
			ECHO Found Visual C++ 2010 express
			set BUILD_X64=0
		)
		
	) ELSE IF EXIST "%ProgramFiles%\Microsoft Visual Studio 10.0\Common7\Tools\vsvars32.bat" (
		set VSVARS32="%ProgramFiles%\Microsoft Visual Studio 10.0\Common7\Tools\vsvars32.bat"
		ECHO Found Visual C++ 2010
	) ELSE (
		echo "Could not find a supported Visual C++ (Express) compiler. Currently supported versions in this template: 2010, 2013 and 2015"
		pause
		goto END
	)
)

IF %DEVENV% NEQ "" (
	IF NOT EXIST %DEVENV% (
		echo "Could not find a suitable Visual C++ (Express) compiler. Supported versions: 2010, 2013, 2015 (including Community and Express editions)."
		echo "You can use the generated Makefile to compile the FMU using the MinGW compiler"
		pause
	)
)

REM Remove potential " from the path, otherwise loading vsvar32 fails
set PATH=%PATH:"=%

IF DEFINED VSVARS32 (
	echo Loading vsvar32
	call %VSVARS32%
)

:FMU_COMPILE
  IF EXIST %PROJ_DIR%\%buildconfig%\buildlog.html del %PROJ_DIR%\%buildconfig%\buildlog.html /q
  IF %buildmode%==clean goto COMPILE_FMU
  IF %buildmode%==noclean goto COMPILE_NO_CLEAN_FMU
  rem ---------------------------------------------
  rem	check for existing dll
  rem ---------------------------------------------

  IF EXIST "%FMU%" (
    echo "Found existing FMU named %FMU%"
    goto FMU_EXIST
  )
  goto COMPILE_FMU

:FMU_EXIST
  IF %promptlevel%==noprompt goto COMPILE_FMU
  ECHO Found an earlier generated FMU
  ECHO [1] a new FMU will be compiled.
  ECHO [2] existing FMU will be updated (quick mode compile).
  ECHO ------------------------------------------------------------
  set /P COMPILE_ANSWER=Compile a new FMU? [1/2]:
  if /I %COMPILE_ANSWER% EQU 1 goto COMPILE_FMU
  if /I %COMPILE_ANSWER% EQU 2 goto COMPILE_NO_CLEAN_FMU

:COMPILE_FMU
  IF EXIST "%FMU%" del /Q "%FMU%"
  msbuild.exe "%PROJ_DIR%\%SUBMODEL_NAME%.sln" /p:Configuration=%buildconfig%;Platform=win32 /t:Build /verbosity:minimal
  IF %BUILD_X64%==1 (
    msbuild.exe "%PROJ_DIR%\%SUBMODEL_NAME%.sln" /p:Configuration=%buildconfig%;Platform=x64 /t:Build /verbosity:minimal
  )

  REM Check whether the 32-bits dll has been generated.
  REM note: 64-bits will not for all visual studio installations be possible
  IF NOT EXIST %PROJ_DIR%\win32\%buildconfig%\%DLL% (
    set DIETEXT="%DLL% failed to build!  See ..\%PROJ_DIR%\%buildconfig%\BuildLog.htm for details."
    goto DIE
  )
  ECHO Build successful.
  GOTO MAKE_FMI

  
:COMPILE_NO_CLEAN_FMU
  ECHO Compiling 20-sim FMU for submodel "%SUBMODEL_NAME%"
  msbuild.exe "%PROJ_DIR%\%SUBMODEL_NAME%.vcxproj" /p:Configuration=%buildconfig%;Platform=win32 /t:Build /verbosity:minimal
  IF %BUILD_X64%==1 (
    msbuild.exe "%PROJ_DIR%\%SUBMODEL_NAME%.vcxproj" /p:Configuration=%buildconfig%;Platform=x64 /t:Build /verbosity:minimal
  )

  IF NOT EXIST "%FMU%" (
    set DIETEXT="%FMU% failed to build!  See ..\%PROJ_DIR%\%buildconfig%\BuildLog.htm for details."
    goto DIE
  )
  ECHO Build successful.
  ECHO ------------------------------------------------------------
  GOTO MAKE_FMI

:MAKE_FMI
  cd %CURPATH%
  if not exist "%BIN32_DIR%" mkdir "%BIN32_DIR%"
  ECHO Copy the compiled DLL %PROJ_DIR%\Win32\%buildconfig%\%DLL% to %BIN32_DIR%
  copy "%PROJ_DIR%\Win32\%buildconfig%\%DLL%" "%BIN32_DIR%"
  IF %BUILD_X64%==1 (
    if not exist "%BIN64_DIR%" mkdir "%BIN64_DIR%"
    ECHO Copy the compiled DLL %PROJ_DIR%\x64\%buildconfig%\%DLL% to %BIN64_DIR%
    copy "%PROJ_DIR%\x64\%buildconfig%\%DLL%" "%BIN64_DIR%"
  )

  ECHO copy the generated sources to %SRC_DIR%
  copy "%ROOTPATH%\src\*.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\*.h" "%SRC_DIR%"
%IF%%EQ(INTEGRATION_METHOD_NAME,VodeAdams)%
  ECHO copy the CVode sources to %SRC_DIR%
  mkdir "%SRC_DIR%\include"
  mkdir "%SRC_DIR%\include\cvode"
  mkdir "%SRC_DIR%\include\nvector"
  mkdir "%SRC_DIR%\include\sundials"
  copy "%ROOTPATH%\src\includes.txt" "%SRC_DIR%"
  copy "%ROOTPATH%\src\cvode\src\cvode\cvode.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\cvode\src\cvode\cvode_dense.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\cvode\src\cvode\cvode_direct.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\cvode\src\cvode\cvode_io.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\cvode\src\nvec_ser\nvector_serial.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\cvode\src\sundials\sundials_dense.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\cvode\src\sundials\sundials_direct.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\cvode\src\sundials\sundials_math.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\cvode\src\sundials\sundials_nvector.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\cvode\include\cvode\cvode.h" "%SRC_DIR%\include\cvode\"
  copy "%ROOTPATH%\src\cvode\include\cvode\cvode_direct.h" "%SRC_DIR%\include\cvode\"
  copy "%ROOTPATH%\src\cvode\include\cvode\cvode_dense.h" "%SRC_DIR%\include\cvode\"
  copy "%ROOTPATH%\src\cvode\src\cvode\cvode_direct_impl.h" "%SRC_DIR%"
  copy "%ROOTPATH%\src\cvode\src\cvode\cvode_impl.h" "%SRC_DIR%"
  copy "%ROOTPATH%\src\cvode\include\nvector\nvector_serial.h" "%SRC_DIR%\include\nvector"
  copy "%ROOTPATH%\src\cvode\include\sundials\sundials_math.h" "%SRC_DIR%\include\sundials"
  copy "%ROOTPATH%\src\cvode\include\sundials\sundials_nvector.h" "%SRC_DIR%\include\sundials"
  copy "%ROOTPATH%\src\cvode\include\sundials\sundials_config.h" "%SRC_DIR%\include\sundials"
  copy "%ROOTPATH%\src\cvode\include\sundials\sundials_dense.h" "%SRC_DIR%\include\sundials"
  copy "%ROOTPATH%\src\cvode\include\sundials\sundials_direct.h" "%SRC_DIR%\include\sundials"
  copy "%ROOTPATH%\src\cvode\include\sundials\sundials_types.h" "%SRC_DIR%\include\sundials"
%ENDIF%
%IF%%EQ(INTEGRATION_METHOD_NAME,MeBDFi)%
  ECHO copy the MeBDFi sources to %SRC_DIR%
  mkdir "%SRC_DIR%\include"
  copy "%ROOTPATH%\src\includes.txt" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\MeBDFi.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\MeBDFi.h" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\include\f2c.h" "%SRC_DIR%\include"
  copy "%ROOTPATH%\src\MeBDFi\include\fio.h" "%SRC_DIR%\include"
  copy "%ROOTPATH%\src\MeBDFi\include\fmt.h" "%SRC_DIR%\include"
  copy "%ROOTPATH%\src\MeBDFi\include\fp.h" "%SRC_DIR%\include"
  copy "%ROOTPATH%\src\MeBDFi\include\lio.h" "%SRC_DIR%\include"
  copy "%ROOTPATH%\src\MeBDFi\include\sysdep1.h" "%SRC_DIR%\include"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\close.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\dolio.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\endfile.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\err.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\fmt.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\fmtlib.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\lread.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\lwrite.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\open.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\pow_dd.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\pow_di.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\rdfmt.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\rsfe.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\sfe.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\sig_die.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\util.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\wref.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\wrtfmt.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\wsfe.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\libf2c\wsle.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\netlib\DOGLEG.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\netlib\DPMPAR.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\netlib\ENORM.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\netlib\FDJAC1.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\netlib\HYBRD.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\netlib\HYBRD1.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\netlib\mebdfiCore.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\netlib\QFORM.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\netlib\QRFAC.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\netlib\R1MPYQ.c" "%SRC_DIR%"
  copy "%ROOTPATH%\src\MeBDFi\netlib\R1UPDT.c" "%SRC_DIR%"
%ENDIF%

  ECHO Generate the FMU
  cd "%FMU_DIR%"
  if exist "%FMU%" del /Q "%FMU%"
  "%ZIPTOOL%" a -tzip -xr!.svn "%FMU%" *
  cd "%CURPATH%"

  IF NOT EXIST "%FMU%" (
    set DIETEXT="%FMU% failed to build!  See ..\%PROJ_DIR%\%buildconfig%\BuildLog.htm for details."
    goto DIE
  )
  ECHO ------------------------------------------------------------
  GOTO FMU_READY


:FMU_READY
  cd ..
  ECHO Your %SUBMODEL_NAME% FMU is ready and can be found in:
  ECHO   %CD%
  ECHO The name of your FMU is:
  ECHO   %FMU%
  ECHO ------------------------------------------------------------
  SET exitcode=0
  cd %CURPATH%
  GOTO END

:DIE
  set DIETEXT=Error: %DIETEXT%
  echo %DIETEXT%
  SET exitcode=1
  ECHO ------------------------------------------------------------

:END
  IF %promptlevel% NEQ noprompt (
  ECHO Press any key to exit...
  pause > NUL
  )
  IF %exitcode% NEQ 0 EXIT /B %exitcode%
  EXIT

