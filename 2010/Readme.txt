Example-2.r  : file for Example 2 in the draft
             (I just listed a part of the results for brevity)
Example-2.out : Output after running Example-2.r 

Put these three files at the same directory. 
Just look at "Example-2.r" file.

The following might be convenient for Linux/Unix users.

% R --slave --no-save < Example-2.r > Example-2.out 
% more Example-2.out 
