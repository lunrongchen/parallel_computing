# parallel_computing
Multithread programming for Intro to Parallel Computing


ssh lct595@358smp.eecs.northwestern.edu


source /etc/profile.d/modules.csh
source /etc/profile.d/mpich-x86_64.csh
cd hw3/
mpicc -o executable recursive_bisection.c

mpiexec -launcher rsh -f machine_list_cluster -n 128 ./executable 256


setenv PATH /usr/lib64/mpich/bin:$PATH
setenv MPDIR /usr/lib64/mpich/bin 