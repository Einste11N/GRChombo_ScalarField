# GRChombo_ScalarField


cd GRChombo-main/Examples/ScalarFieldnexttoPBH


make realclean

make all -j 8

mpirun -np 8 XXX.ex params.txt
