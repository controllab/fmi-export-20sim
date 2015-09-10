# fmi-export-20sim
FMI export code generation template for 20-sim

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
1. For automated compilation of the .fmu, we currently rely on the availability of the Visual Studio 2010 or Visual Studio 2013 compiler.
(Any version Express, Professional or Community will do).

It can be compiled with gcc/MinGW too, but the addition of a proper Makefile is still a TODO.

2. [Microsoft Visual C++ Redistributable Packages for Visual Studio 2013](http://www.microsoft.com/en-us/download/details.aspx?id=40784)
