@echo off
cls
SET "curDirBefore=%CD%"
pushd "%~dp0"
SET "curDirAfter=%CD%"
IF "%curDirBefore%" NEQ "%curDirAfter%" (
	:: This means that PUSHD did map a temporary drive letter
	:: We need to wait until build.bat is finished otherwise POPD will remove the drive letter again
	set "WAIT=/wait"
)
IF "%XXSIM_SCRIPT_MODE%" == "1" (
	set "WAIT=/wait"
)
cd project
start %WAIT% build.bat
popd
