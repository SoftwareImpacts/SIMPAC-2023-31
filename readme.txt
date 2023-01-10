*# INSTRUCTIONS - read.me
StochHyProp v. 2.01

*# This software will generate a number [Nstoch] of realizations of parameters according to Mean, Standard Error and correlation matrix following Cholesky decomposition.
*# Running parameters will be read from this file, default name [StochHyProp.stc] - another name may be specified in the command line. Default extension (if none specified): [.stc]
*# Output will be written to a comma separated values (.csv) file, and to a space separated (.out) file. Default name is the same as the input file, unless overridden by [OutFile].

*# Parameter names may be changed in the input .stc file and will be copied to the output file.
*#
*# The program may be run in batch mode, in which case it may be convenient to specify [Lkey]=0, avoiding the program to wait for a key-press at the end.
*#
*# Special feature 1: TAIL exclusion
*# In order to avoid extreme values, a tail fraction may be specified as [tau]. The tail will be applied on both sides (e.g., Tau=0.01 will eliminate a total of 2% of values).
*# A consequence of tail exclusion is a reduction of the standard error. This can be corrected for by specifying [Lsterr] = 1.
*# Note that the use of a tail will corrupt the assumed normal distribution and should be avoided or limited to the smallest reasonabel value.
*#
*# Special feature 2: Read statistics directly from RETC or Hydrus (inverse modeling) output files
*# Specify the [Lext] flag to 1 and the program will look for Mean, Standard Error and correlation matrix in the RETC or Hydrus ouptput files specified as [ExtFile]
*# instead of in the subsequent lines in this file.