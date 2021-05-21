# Readme

With this repository we make the code available used to generate the data presented in the article:

**Rotating shallow water flow under location uncertainty with a structure-preserving discretization** 

Reproducing similar data (since the simulations are stochastic) from this article can be done as follows:

 - First, create the output folders with:

	`mkdir -p sphere/LES`

	`mkdir -p sphere/{PIC-POD,PIC-SVD,SVD,POD}/run-{1..20}`

	`mkdir -p plane/det/{1,5,10,50,100}`

	`mkdir -p plane/homNoise/{1,5,10,50,100}/run-{1..10}`

 - Testcase on the plane: 
	 - run `runDetPlane(DT)` for *DT*=1, 5, 10, 50, 100 

		 -> This creates in the folders *plane/det/DT/* a file *final.mat* containing the simulation data after 2 days.
	- run `runstochasticDiffu_homgNoise(runNr,DT)`  for *runNr*=1..10 and *DT*=1, 5, 10, 50, 100 
	
		->This creates in the folders *plane/homNoise/DT/runNr/*  a file *final.mat* containing the simulation data after 2 days.
- Testcase on the plane:
	- run `runLES()` 
	
		-> This creates in the folder *sphere/LES/* files *1.mat* ... *80.mat*
	- run `runPIC_POD(runNR)` for *runNR*=1..20 
	
		-> This creates in the folders *sphere/PIC-POD/run-runNR/* files *4.mat* ... *80.mat*
	- run `runPIC_SVD(runNR)` for *runNR*=1..20 
	
		-> This creates in the folders *sphere/PIC-SVD/run-runNR/* files *4.mat* ... *80.mat*
	- run `runstochastic(runNR)` for *runNR*=1..20 
	
		-> This creates in the folders *sphere/SVD/run-runNR/* files *1.mat* ... *80.mat*
	- run `runstochasticPOD(runNR)` for *runNR*=1..20 
	
		-> This creates in the folders *sphere/POD/run-runNR/* files *1.mat* ... *80.mat*
	
	The fields are saved every 6 hour until day 20, this corresponds to the output files *1.mat* ...*80.mat*
