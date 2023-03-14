# Subject-Specific Factors Affecting Particle Residence Time Distribution of Left Atrial Appendage in Atrial Fibrillation: A Computational Model-Based Study

## Abstract
Background: Atrial fibrillation (AF) is a prevalent arrhythmia, that causes thrombus formation, ordinarily in the left atrial appendage (LAA). The conventional metric of stroke risk stratification, CHA2DS2-VASc score, does not account for LAA morphology or hemodynamics. We showed in our previous study that residence time distribution (RTD) of blood-borne particles in the LAA and its associated calculated variables (i.e., mean residence time, tm, and asymptotic concentration, C∞) have the potential to improve CHA2DS2-VASc score. The purpose of this research was to investigate the effects of the following potential confounding factors on LAA tm and C∞: (1) pulmonary vein flow waveform pulsatility, (2) non-Newtonian blood rheology and hematocrit level, and (3) length of the simulation.
Methods: Subject-Specific data including left atrial (LA) and LAA cardiac computed tomography, cardiac output (CO), heart rate, and hematocrit level were gathered from 25 AF subjects. We calculate LAA tm and C∞ based on series of computational fluid dynamics (CFD) analyses.
Results: Both LAA tm and C∞ are significantly affected by the CO, but not by temporal pattern of the inlet flow. Both LAA tm and C∞ increase with increasing hematocrit level and both calculated indices are higher for non-Newtonian blood rheology for a given hematocrit level. Further, at least 20,000 s of CFD simulation is needed to calculate LAA tm and C∞ values reliably.
Conclusions: Subject-specific LA and LAA geometries, CO, and hematocrit level are essential to quantify the subject-specific proclivity of blood cell tarrying inside LAA in terms of the RTD function.
## OpenFOAM Solvers
Several OpenFOAM (version 8, The OpenFOAM Foundation Ltd, Inc., UK.) solvers and boundary conditions were developed for this study based on standard OpenFOAM codes. To investigate the effects of patient-specific PV flow waveform, we developed the **ScalarAdvection** solver which is based on original icoFoam solver of OpenFOAM. In this solver transport equation is coupled with the momentum equations. To investigate the effects of patient-specific hematocrit levels and non-Newtonian vs. Newtonian fluid modeling, we developed **nonNewtonianDistVel**, **icoFoamDistVel**, and **passiveScalarAdvection** solvers which are based on original nonNewtonianIcoFoam, icoFoam, and ScalarTransportFoam solvers of OpenFOAM, respectively. In our final simulations, nonNewtonianDistVel and passiveScalarAdvection solvers were used.
## Code of Conduct
Please cite our paper if you find this repository usefull:</br>
Sanatkhani, S., Nedios, S., Menon, P. G., Saba, S. F., Jain, S. K., Federspiel, W. J., & Shroff, S. G. (2023). Subject-specific factors affecting particle residence time distribution of left atrial appendage in atrial fibrillation: A computational model-based study. <i>Front Cardiovasc Med</i>, 10(1070498), 1-13. <a href = "https://doi.org/10.3389/fcvm.2023.1070498" target="_blank">doi: 10.3389/fcvm.2023.1070498</a></li>
