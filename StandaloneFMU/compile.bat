@echo off
cd project
IF "%XXSIM_SCRIPT_MODE%" == "1" (
	start /wait build.bat
) ELSE (
	start build.bat
)
cd ..