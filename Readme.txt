
Update 2021.02.18 the new version is the NILS

-----------------------
This work is already accepted by KBS.

Compile file: 
g++ BCP_main.cpp BCP_c_LS.cpp BCP_mem.cpp BCP_X.cpp -O3 -lm -Wall -o BCP_1

Run file:
BCP_1 -i tree10x2.rnd --seed 0 -rep 0 -alb 28 -L1 100 -L2 20 -L3 2000 -alpha 0.84

The benchmark instances are in the folder ITPS.

The option of "-alb" is the lower bound of the instance. It is important.
