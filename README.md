# fmi-export-20sim
FMI export code generation template for 20-sim

## Requirements
This template requires at least 20-sim version 4.6.0.

## Installation
Download the code generation template using the "Download ZIP" button on Github or using a "git clone"

1. Extract the zip file
2. Copy the folder StandaloneFMU to a location on your 20-sim PC.
3. Open 20-sim
4. Choose from the Menu: Tools | Options
5. Choose Folders
6. Choose C-Code Folders
7. Add the folder you copied in step 2.
8. To export an FMU, choose Tools | Real Time Toolbox | C-Code Generation
    The FMU export should now be visible in the list.

## External dependencies
1. For automated compilation of the .fmu, we currently rely on the availability of the Visual Studio 2010, 2013 or 2015 compiler.
Any version Express, Community, Professional or Enterprise will work. However, for a x64 FMU DLL, you cannot use the Express edition since it does not contain a x64 compiler.

Note that the FMU can be compiled with gcc (including MinGW) too on Windows, Linux and MacOSX. You can find the Makefile in the root directory of the generated code.

2. [Microsoft Visual C++ Redistributable Packages for Visual Studio 2013](http://www.microsoft.com/en-us/download/details.aspx?id=40784)
