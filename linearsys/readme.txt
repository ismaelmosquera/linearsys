Building the linearsys shared library:

Maybe the library is already compiled, look at the lib folder to confirm it. Anyway, the library would be compiled as a *.dll, that is, linearsys.dll since the author uses the MinGW ( Minimalist GNU for Windows ), which is the Microsoft Windows GNU compilers suite.
If you are, for example, under Linux, please change in the Makefile the following: .dll to .so in the cleanup macro and in the shared library target; then, you will compile the desired library since this implementation uses only standard headers and libraries which are launched with any C compiler.
To compile the library just type:

make

And the shared library will be placed inside the lib folder.
 