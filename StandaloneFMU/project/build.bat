@ECHO OFF
rem ----Usage----
rem build [clean|noclean]
rem vs2008 for compiling with visual studio 2008
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
SET comp=vs2013
SET promptlevel=prompt
SET exitcode=0
SET buildmode=clean

cd ..\
SET ROOTPATH=%CD%
cd %CURPATH%

SET FMU=%ROOTPATH%\%SUBMODEL_NAME%.fmu
SET DLL=%SUBMODEL_NAME%.dll
SET ZIPTOOL=%TEMPLATE_DIR%\bin\7z.exe
SET XSLTTOOL=%TEMPLATE_DIR%\bin\msxsl.exe
SET GUIDTOOL=%TEMPLATE_DIR%\bin\GenerateGuid.exe

set FMU_DIR=%CURPATH%fmu
set BIN_DIR=%FMU_DIR%\binaries\win32
set SRC_DIR=%FMU_DIR%\sources
set DOC_DIR=%FMU_DIR%\documentation

FOR %%b in (%1, %2, %3, %4, %5) DO (
	IF %%b==vs2010 SET comp=vs2010
	IF %%b==vs2013 SET comp=vs2013
	IF %%b==clean SET buildmode=clean
	IF %%b==noclean SET buildmode=noclean
	IF %%b==noprompt SET promptlevel=noprompt
)

SET buildconfig=Release
SET DEVENV=""
SET VSVARS32=""

ECHO ------------------------------------------------------------
ECHO 20-sim standalone co-simulation FMU export for '%SUBMODEL_NAME%'
ECHO ------------------------------------------------------------
ECHO Searching for Visual C++ compiler...

rem Seach for VS 2013 / VS 2013 Express / VS 2013 Community edition
IF %comp%==vs2013 (
	set PROJ_DIR=VS2013
	IF EXIST "%VS120COMNTOOLS%\vsvars32.bat" (
		set VSVARS32="%VS120COMNTOOLS%\vsvars32.bat"
		ECHO Found Visual C++ for Desktop 2013
	) ELSE IF EXIST "%ProgramFiles%\Microsoft Visual Studio 12.0\Common7\Tools\vsvars32.bat" (
		set VSVARS32="%VS120COMNTOOLS%\vsvars32.bat"
		ECHO Found Visual C++ for Desktop 2013
	) ELSE (
		rem Try an older compiler
		set comp=vs2010
	)
)

rem Seach for VS 2010 / VS 2010 Express
IF %comp%==vs2010 (
	set PROJ_DIR=VS2010
	IF EXIST "%VS100COMNTOOLS%\..\IDE\devenv.exe" (
		set DEVENV="%VS100COMNTOOLS%\..\IDE\devenv.exe"
		ECHO Found Visual C++ 2010
	) ELSE IF EXIST "%VS100COMNTOOLS%\..\IDE\VCExpress.exe" (
		set DEVENV="%VS100COMNTOOLS%\..\IDE\VCExpress.exe"
		ECHO Found Visual C++ Express 2010
	) ELSE IF EXIST "%ProgramFiles%\Microsoft Visual Studio 10.0\Common7\IDE\VCExpress.exe" (
		set DEVENV="%ProgramFiles%\Microsoft Visual Studio 10.0\Common7\IDE\VCExpress.exe"
		ECHO Found Visual C++ Express 2010
	) ELSE (
		echo "Could not find a supported Visual C++ (Express) compiler. Currently supported versions in this template: 2010 and 2013"
		pause
		goto END
	)
)


IF %DEVENV% NEQ "" (
	IF NOT EXIST %DEVENV% (
		echo "Could not find a suitable Visual C++ (Express) compiler. Supported versions: 2010, 2013 Community/Desktop (including Express editions)."
		echo "You can use the generated Makefile to compile the DLL using the MinGW compiler"
		pause
	)
)

IF %VSVARS32% NEQ "" (
	echo Loading vsvar32
	call %VSVARS32%
)

set OPTS_FMU="%PROJ_DIR%\%SUBMODEL_NAME%.sln" /build "%buildconfig%"
set CLEAN_FMU="%PROJ_DIR%\%SUBMODEL_NAME%.sln" /clean "%buildconfig%"

rem ECHO %PROJ_DIR%
rem ECHO %OPTS_DLL%
rem ECHO %CLEAN_DLL%

ECHO ------------------------------------------------------------
rem	CONFIG END
rem -------------------------------------------------------------

"%GUIDTOOL%" > "%ROOTPATH%\src\guid.txt"
rem generate GUID header
FOR /f "tokens=*" %%g IN (%ROOTPATH%\src\guid.txt) DO (
	rem generate a header with a define for the GUID
	echo #define FMI_GUID "{%%g}" > "%ROOTPATH%\src\fmiGUID.h"

	REM Generate a xml file (token format) with the generated GUID
	(  
		echo ^<?xml version="1.0" encoding="UTF-8"?^>
		echo ^<tokens^>
		echo ^<token name="GUID"^>^<![CDATA[%%g]]^>^</token^>
		echo ^</tokens^>
	) > "%ROOTPATH%\src\GUID.xml"
)

:FMU_COMPILE
  IF EXIST %PROJ_DIR%\%buildconfig%\buildlog.html del %PROJ_DIR%\%buildconfig%\buildlog.html /q
  IF %buildmode%==clean goto COMPILE_FMU
  IF %buildmode%==noclean goto COMPILE_NO_CLEAN_FMU
  rem ---------------------------------------------
  rem	check for existing dll
  rem ---------------------------------------------

  IF EXIST ..\%FMU% (
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
  IF EXIST ..\%FMU% del /Q ..\%FMU%
  IF %DEVENV% NEQ "" (
    ECHO Cleaning Solution...
    %DEVENV% "%PROJ_DIR%\%SUBMODEL_NAME%.sln" /clean "%buildconfig%"
    ECHO Compiling 20-sim FMU for submodel "%SUBMODEL_NAME%"
    %DEVENV% "%PROJ_DIR%\%SUBMODEL_NAME%.sln" /build "%buildconfig%"
  ) ELSE (
    ECHO Compiling 20-sim FMU for submodel "%SUBMODEL_NAME%"
    msbuild.exe "%PROJ_DIR%\%SUBMODEL_NAME%.vcxproj" /p:Configuration=%buildconfig% /t:Rebuild /verbosity:minimal
  )

  IF NOT EXIST %PROJ_DIR%\%buildconfig%\%DLL% (
    set DIETEXT="%DLL% failed to build!  See ..\%PROJ_DIR%\%buildconfig%\BuildLog.htm for details."
    goto DIE
  )
  ECHO Build successful.
  GOTO MAKE_FMI

  
:COMPILE_NO_CLEAN_FMU
  IF %DEVENV% NEQ "" (
    ECHO Compiling 20-sim FMU for submodel "%SUBMODEL_NAME%"
    %DEVENV% "%PROJ_DIR%\%SUBMODEL_NAME%.sln" /build "%buildconfig%"
  ) ELSE (
    ECHO Compiling 20-sim FMU for submodel "%SUBMODEL_NAME%"
    msbuild.exe "%PROJ_DIR%\%SUBMODEL_NAME%.vcxproj" /p:Configuration=Release /t:Build /verbosity:minimal
  )

  IF NOT EXIST ..\%FMU% (
    set DIETEXT="%FMU% failed to build!  See ..\%PROJ_DIR%\%buildconfig%\BuildLog.htm for details."
    goto DIE
  )
  ECHO Build successful.
  ECHO ------------------------------------------------------------
  GOTO MAKE_FMI

:MAKE_FMI
  cd %CURPATH%
  ECHO Creating an empty FMU
  if not exist %FMU_DIR% mkdir %FMU_DIR%
  if not exist %BIN_DIR% mkdir %BIN_DIR%
  if not exist %SRC_DIR% mkdir %SRC_DIR%
  if not exist %DOC_DIR% mkdir %DOC_DIR%
  ECHO Copy the compiled DLL %PROJ_DIR%\%buildconfig%\%DLL% to %BIN_DIR%
  copy "%PROJ_DIR%\%buildconfig%\%DLL%" "%BIN_DIR%"

  ECHO copy the generated sources to %SRC_DIR%
  copy "%ROOTPATH%\src\*.*" "%SRC_DIR%"

  ECHO Generate the modelDescription.xml
  "%XSLTTOOL%"  "%ROOTPATH%\src\ModelConfiguration.xml" "%ROOTPATH%\template\mcf2modelDescription.xsl" -o "%FMU_DIR%\modelDescription.xml" SOURCEDIRECTORY="%ROOTPATH%\src"

  ECHO Generate the FMU
  cd %FMU_DIR%
  if exist %FMU% del /Q %FMU%
  "%ZIPTOOL%" a -tzip -xr!.svn %FMU% *
  cd %CURPATH%

  IF NOT EXIST %FMU% (
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
