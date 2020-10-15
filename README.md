# ZETA_analysis_repository
Some matlab code to run the benchmarks for ZETA

Do the following to run the ZETA benchmarks as shown in this paper:
1) Download the data files from Dryad and add the ZETA repository (https://github.com/JorritMontijn/ZETA) to your path
2) Edit the value of the variable strDataMasterPath in "runResponsivenessBenchmark.m" to match the location of the data files
3) Optionally change some parameters in the first cell, and then run the script.
4) Wait! A waitbar should pop up, but as the data consists of thousands of trials, it will take a while to finish. You can change boolUseSubset to true to use only a subset of the data and reduce the computation time. You can also reduce the number of resamplings, but this will make the results less accurate.

