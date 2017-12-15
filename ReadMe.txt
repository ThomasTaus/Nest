 ------
| Nest |
 ------

It is an R-package that allows you to simulate allele frequency trajectories (under neutrality and with selection) and estimate the effective population size (Ne). You can also load sync-files in R and polarize allele frequencies (e.g. for the rising allele).

PLEASE NOTE: All functions of the Nest package are also available in the poolSeq package (https://github.com/ThomasTaus/poolSeq), which provides additional functionality. It is therefor recommended to use poolSeq instead of Nest.


 ---------------
| Prerequisites |
 ---------------
It has the following dependencies, which need to be installed prior to the Nest package:
	•	R (>= 3.2.0)
	•	data.table (>= 1.9.4)
	•	foreach (>= 1.4.2)
	•	stringi (>= 0.4-1)
	•	matrixStats (>= 0.14.2)

If we can later submit it to CRAN, all the dependent packages will be installed automatically. But for now, the user has to do it manually.


 --------------
| Installation |
 --------------
Click on ‘Clone or download’ to download the whole project (e.g. as ZIP). You will find the most recent version in the folder 'Release’. To install Nest run the following line in an R-session:

install.packages("/Path/To/Nest_<version>.tar.gz", repos=NULL, type="source")


 -------
| Usage |
 -------
After installation you can load it like any other R-package:

library(Nest)

The main functions are:
	•	wf.traj
	•	estimateNe
	•	estimateWndNe
	•	read.sync

You can find detailed descriptions for each function including code examples by preceding the function name with '?' (as usual in R):

?wf.traj


-----------------------------------------------------------------------------------------------------------------------------------------------------
Enjoy using it and feel free to contact me if you have any questions or suggestions.
(Thomas.Taus@vetmeduni.ac.at)
