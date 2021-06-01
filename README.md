# Rotating shallow water flow under location uncertainty with a structure-preserving discretization

The motion of geophysical fluids on the globe needs to be modelled to get some insights of tomorrow's climate. These forecasts must be precise enough while remaining computationally affordable. An ideal system should also deliver, across time, an accurate measurement of the uncertainties introduced through physical or numerical approximations. 
To address these issues, we use the rotating shallow water equations, which provide a simplified version of the dynamics, and a stochastic representation of the unresolved small-scale processes. The former is approximated with a structure preserving numerical model enabling the conservation of physical quantities such as mass and energy and 
the latter is modelled by the location uncertainty framework that relies on stochastic transport and has the great advantage to be energy conserving. 
Our method can directly be used in existing dynamical cores of global numerical weather prediction and climate models. Numerical results illustrate the energy conservation of the numerical model. Simulating a barotropically unstable jet on the sphere, we demonstrate that the random model better captures the structure of a large-scale flow than a comparable deterministic model. The random dynamical system is also shown to be associated with good uncertainty representation.


# Reproducibility

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
- Testcase on the sphere:
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
	


# How to cite 
- Article 


```
@article{brecht2021rotating,
  title={Rotating shallow water flow under location uncertainty with a structure-preserving discretization},
  author={Brecht, R{\"u}diger and Li, Long and Bauer, Werner and M{\'e}min, Etienne},
  journal={arXiv preprint arXiv:2102.03783},
  year={2021}
}
```

- Code:


```
@software{rudigerbrecht_2021_4884919,
  title        = {RudigerBrecht/RSW-LU: First release},
  month        = may,
  year         = 2021,
  publisher    = {Zenodo},
  version      = {v.1.0},
  doi          = {10.5281/zenodo.4884919},
  url          = {https://doi.org/10.5281/zenodo.4884919}
}
```

