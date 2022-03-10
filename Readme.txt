This is a short summary on how to run this software.

The code is run from the Matlab command line. For this purpose we have two files in the home folder:
MainSerial and MainParallel. MainSerial is used if you want everything to run in a serialized mode,
while MainParallel is used if you want to run the frequency sweeps in parallel (faster).

To run these files you just have to type:

mainSerial(arg1, arg2,..., argN)

or 

mainParallel(arg1, arg2,..., argN)

where arg1, agr2,..., argN are the function arguments. MainSerial takes 10 arguments, while mainParallel takes 11.
The description of what each argument means is included in the software user manual.

Several problems files are already defined in the problemFiles folder, and the meshes used are also included
in the meshes folder. The problems are described in my thesis and papers, and the structure of the problem files 
is also described in the user manual.

Furthermore, further to the function arguments, there are several switches defined at the beginning of the MainSerial
and MainParallel scripts that are used to turn on/off different options (coupling, POD, etc). Again, these switches
are described in the user manual.