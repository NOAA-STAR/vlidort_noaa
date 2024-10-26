rm *.o *.mod TNEW TOLD *~
exit

gfortran -c -Wall vfzmat_Numerical.f90
gfortran -c -Wall vfzmat_Rotation.f90
gfortran -c -Wall vfzmat_PhasMat.f90
gfortran -c -Wall vfzmat_DevelopCoeffs_New.f90
gfortran -c -Wall vfzmat_Pre_Master.f90
gfortran -c -Wall vfzmat_Post_Master.f90

gfortran -c -Wall vfzmat_DevelopCoeffs.f90
gfortran -c -Wall vfzmat_ExpandCoeffs.f90
gfortran -c -Wall vfzmat_Master.f90


gfortran -o TNEW -Wall Test_vfzmat_reeng.f90 vfzmat_Numerical.f90 vfzmat_Rotation.f90 vfzmat_PhasMat.f90 vfzmat_DevelopCoeffs_New.f90 vfzmat_Pre_Master.f90 vfzmat_Post_Master.f90
gfortran -o TOLD -Wall Test_vfzmat_original.f90 vfzmat_Numerical.f90 vfzmat_Rotation.f90 vfzmat_PhasMat.f90 vfzmat_DevelopCoeffs.f90 vfzmat_ExpandCoeffs.f90 vfzmat_Master.f90

exit


