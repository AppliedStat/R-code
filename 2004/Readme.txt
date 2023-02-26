R Codes for Examples in the paper (https://doi.org/10.1109/TR.2003.821946)


Example.r  : file for Example in the draft
             (I just listed a part of the results for brevity)
Example.out : Output after running Example.r 

Put these three files at the same directory. 
Just look at "Example.r" file.

The following might be convenient for Linux/Unix users.

% R --slave --no-save < Example.r > Example.out 
% more Example.out 
