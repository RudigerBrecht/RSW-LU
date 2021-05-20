# Readme

With this repository we make the code available used to generate the data presented in the article:
**Rotating shallow water flow under location uncertainty with a structure-preserving discretization** 

Reproducing the data from this article can be done as follows:

 - Testcase on the plane: 
	 run `runstochasticDiffu_homgNoise(runNr,DT,128)`
	  for *runNr*=1..10 and *DT*=1, 5, 10, 50, 100
	 
	This creates in the folders *plane/128/DT/runNr*  a file *final.mat* containing the simulation data and diagnostics after 2 days.
- Testcase on the plane:
	- run `runLES(10242)` 
		-> This creates in the folder *sphere/LES/* files *1.mat* ... *80.mat*
	- run `runPIC_POD(runNR,10242)` for *runNR*=1..20 
		-> This creates in the folders *sphere/PIC-POD/run-runNR/* files *1.mat* ... *80.mat*
	- run `runPIC_SVD(runNR,10242)` for *runNR*=1..20 
		-> This creates in the folders *sphere/PIC-SVD/run-runNR/* files *1.mat* ... *80.mat*
	- run `runstochastic(runNR,10242)` for *runNR*=1..20 
		-> This creates in the folders *sphere/SVD/run-runNR/* files *1.mat* ... *80.mat*
	- run `runstochasticPOD(runNR,10242)` for *runNR*=1..20 
		-> This creates in the folders *sphere/POD/run-runNR/* files *1.mat* ... *80.mat*
	
	The fields are saved every 6 hour until day 20, this corresponds to the output files *1.mat* ...*80.mat*
