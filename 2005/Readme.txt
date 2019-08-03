Example-A.r : file for Example in the draft
Example-B.r : file for Example in the draft


The following might be convenient for Linux/Unix users.

% R --slave --no-save < Example-A.r > Example-A.out 
% more Example-A.out 

% R --slave --no-save < Example-B.r > Example-B.out 
% more Example-B.out 
