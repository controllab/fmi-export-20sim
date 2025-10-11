@echo on
cd "%GENERATION_DIR%"
"%20SIM_DIR%\bin\tokenparser.exe" /c "%GENERATION_DIR%\src\TokenParser.xml"
cd "%GENERATION_DIR%\project"
start build.bat
cd ..