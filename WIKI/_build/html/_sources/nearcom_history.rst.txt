*************************
**What is NearCoM**
*************************

* History of NearCoM model development

NearCoM (or Nearshore Community Model) is initially an extensible, user-configurable model system for nearshore wave, circulation and sediment processes developed during the National Oceanographic Partnership Program (NOPP). The model consists of a "backbone", i.e., the master program, handling data input and output as well as internal storage, together with a suite of "modules", each of which handles a focused subset of the physical processes being studied. A wave module will model wave transformation over arbitrary coastal bathymetry and predict radiation stresses and wave induced mass fluxes. A circulation module will model the slowly varying current field driven by waves, wind and buoyancy forcing, and will provide information about the bottom boundary layer structure. A seabed module will model sediment transport, determine the bedform geometry, parameterize the bedform effect on bottom friction, and compute morphological evolution resulting from spatial variations in local sediment transport rates. The model package includes more than 10 modules developed by a group of institutions involving this project. The initially package can be downloaded from :download:`here <download/initial_nearcom.zip>` 

Modules in the original NearCoM system are developed specifically for predicting nearshore waves and wave-induced nearshore processes, basically the region between the shoreline and about 10 m water depth. Hence, the applications of NearCoM to ocean-exposed coastal regions are limited. There is a growing demand recently to use NearCoM in more general coastal applications such as storm-induced coastal inundation, beach and dune erosion, and wave–current interaction in inlet systems. Therefore, a new model coupling system called NearCoM-TVD was developed based on one of the combinations of the original NearCoM system. NearCoM-TVD couples a newly developed nearshore circulation model, SHORECIRC, using a hybrid finite-difference finite-volume TVD-type scheme, the wave model SWAN and several sediment transport modules, such as Kobayshi et al. (2008), Soulsby (1997) and van Rijn et al. (2011).

Detailed documentation of NearCoM-TVD can be found in Chen et al. (2014) or this Wiki page. 

Chen J.-L., Shi F., Hsu T.-J., and Kirby, J. T., 2014, NearCoM-TVD - A quasi-3D nearshore circulation and sediment transport model, *Coastal Engineering*, **91**, 200-212. 

* New features

 * A conservative form of SHORECIRC equations in non-orthogonal curvilinear coordinates; 
 * MUSCL-TVD solver with adaptive Runge-Kutta time stepping;  
 * Wetting-drying moving boundary condition with incorporation of HLL construction method into the scheme; 
 * Fully parallel computation with equal CPU load for each processor.  
 * Concurrent Correction Method (Shi et al., 2015);

  Shi, F., Vittori, G. and Kirby, J. T., 2015, “Concurrent correction method for modeling morphological response to dredging an offshore sandpit”, Coastal Engineering , 97,1-10.

 