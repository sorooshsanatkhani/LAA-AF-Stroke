This solver is consisted of a main function defined in ScalarAdvection.C and a boundary condition for PV inlets, DistributedVelocityInlet. Continuity equation, momentum equations and transport equations are solved in ScalarAdvection.C. A Piso loop is being used to solve momentum equations. Transport equation is being solved at each time step after pressure and velocity calculation in Piso loop. Further, PV flow is read from a CSV file at each time step and velocity at PV inlets is calculated based on the total surface area of PV inlets.

Directory structure for setup of a sample subject simulation using ScalarAdvection is presented below:
![Directory Structure](ScalarAdvection/DirectoryStructure.jpg "Structure")

where LAA surface geometry is presented in LAA_ID.stl. PV flow throughout a cardiac length is stored in PV_Flow.csv in two columns separated by a comma, where time is in the first column and PV flow (L min-1) in second column. Folder 0 is consisted of boundary conditions and initial values for pressure, velocity, and tracer concentration. CFD Mesh information is stored in polyMesh folder within the constant folder. Viscosity and density of blood are set in transportProperties file. Discretization schemes and model tolerances are defined in fvSchemes and fvSolution respectively. Number of domains that the CFD mesh is decomposed into for parallel processing is set in decomposeParDict. setFieldsDict is used to set the initial tracer concentration inside the LAA to 1. topoSetDict is used to define LAA as a zone so that it can later be used in a function in controlDict to calculate volumetric average of tracer concentration inside the LAA at each time step.


### Boundary conditions used with ScalarAdvection solver
				U	                    p	                    s
PVs		DistributedVelocityInlet	    zeroGradient	           fixedValue                                                                                                   uniform 0
Mitral	      zeroGradient	          fixedValue
                                      uniform 0	                zeroGradient
Wall	            noSlip	            zeroGradient	            zeroGradient


