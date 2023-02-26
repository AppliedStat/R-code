R Codes for Example in the paper (https://doi.org/10.1016/j.jspi.2004.12.016)


Example.r  : file for Example in the draft
             (I just listed a part of the results for brevity)
Example.out : Output after running Example.r 

Put these three files at the same directory. 
Just look at "Example.r" file.

The following might be convenient for Linux/Unix users.

% R --slave --no-save < Example.r > Example.out 
% more Example.out 
