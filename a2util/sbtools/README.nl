
The nl_* series of routines implements a fortran namelist.

A namelist starts with an "*" in the first column followed immediately
by the name of the namelist.  This is followed by any number of fields
of the form KEY=VAL separated by spaces or newlines.  Any number of
fields may appear on a single line provided the line is no more than 80
characters long.  The end of the namelist is signaled by the end of the
file or a line with an "*" in the first column.

Keys do not have to be a rigidly defined word.  Rather, they can be
abbreviated.  A sample keyword might be "CALC*ULATION".  This means that
the first four characters are required, and any remaining characters may
be omitted, but if they are given, must match the keyword appropriately.
In other words, "CALC", "CALCUL" and "CALCULATION" all work but "CALCCC"
won't.

