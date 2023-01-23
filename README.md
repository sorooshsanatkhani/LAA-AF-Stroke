# Hemodynamic Indices and Shape-Based Models of Left Atrial Appendage to Enhance Stroke Prediction in Atrial Fibrillation
## Abstract
Atrial fibrillation (AF) is the most common arrhythmia that leads to thrombus formation, mostly in the left atrial appendage (LAA). The current standard of stratifying stroke risk, based on the CHA2DS2-VASc score, does not consider LAA morphology/hemodynamics. The aim of this study was to determine whether LAA morphology and hemodynamics-based indices can stratify stroke risk independent of CHA2DS2-VASc score, left atrium size, and AF type. In a retrospective matched case-control study, patient-specific measurements in 128 AF patients included left atrial (LA) and LAA 3D geometry obtained by cardiac computed tomography, heart rate, cardiac output, and hematocrit. We quantified patient-specific 3D LAA morphology in terms of a novel LAA appearance complexity index (LAA-ACI) and employed computational fluid dynamics (CFD) analysis to quantify LAA mean residence time, tm and asymptotic concentration, C∞ of blood-borne particles.
Effects of confounding variables were examined to optimize the CFD analysis. cardiac output, but not by the temporal pattern of pulmonary vein inlet flow, significantly affected LAA tm. Both the hematocrit level and the blood rheology model (Newtonian vs. non-Newtonian) also significantly affected LAA tm. Finally, 10,000 s was found to be a sufficient length of CFD simulation to calculate LAA tm in a consistent and reliable manner.
LAA tm varied significantly within a given LAA morphology as defined by the current subjective method, and it was not simply a reflection of LAA geometry/appearance. In addition, LAA-ACI and tm varied significantly for a given CHA2DS2-VASc score, indicating that these two indices of stasis are not simply a reflection of the subjects’ clinical status. Using multiple logistic regression, we observed that ACI, tm, and C∞ had a modest, but statistically insignificant performance in predicting stroke (area under the ROC curve = 0.56–0.61). The temporal dissociation between adverse changes in LAA shape and hemodynamics-based indices and the actual stroke event can contribute to the negative result; a longitudinal study is necessary to address this issue. In addition, it is possible that a multiscale model that combines CFD-based hemodynamics simulation and biology-based thrombus formation can yield indices that can better stratify stroke risk in AF patients.

## OpenFOAM Solvers
Several OpenFOAM solvers and boundary conditions were developed for this dissertation based on original OpenFOAM codes. To investigate the effects of patient-specific PV flow waveform, we developed the ScalarAdvection solver which is based on original icoFoam solver of OpenFOAM. In this solver transport equation is coupled with the momentum equations. To investigate the effects of patient-specific hematocrit levels and non-Newtonian vs. Newtonian fluid modeling, we developed nonNewtonianDistVel, icoFoamDistVel, and passiveScalarAdvection solvers which are based on original nonNewtonianIcoFoam, icoFoam, and ScalarTransportFoam solvers of OpenFOAM, respectively. In our final simulations, nonNewtonianDistVel and passiveScalarAdvection solvers were used.
