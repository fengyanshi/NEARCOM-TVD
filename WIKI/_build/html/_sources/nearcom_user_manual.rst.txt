*******************
**Users' Manual**
*******************

Program outline and flow chart
##################################

The code was written using Fortran 90 with the c preprocessor (cpp) statements for separation of the source code. Arrays are dynamically allocated at runtime. Precision is selected using the *selected_real_kind* Fortran intrinsic function defined in the makefile.  The default precision is single. 

The present version of NearCoM-TVD includes a number of options including (1) choice of serial or parallel code (2) Cartesian or curvilinear coordinate, (3) samples.

The flow chart is shown in :numref:`chart`. 

.. figure:: figures/chart.jpg
   :name: chart
   :scale: 50%
   :align: center

   Flow chart of the main program.





Permanent variables associated with coupling
##############################################

 Depth(): still water depth *h* at element point

 DepthX(): still water depth *h* at x-interface

 DepthY(): still water depth *h* at y-interface

 Eta():   surface elevation, for dry point, Eta() = MinDepth - Depth(), MinDepth is specified in input.txt. 

 Eta0(): :math:`\eta` at previous time level

 MASK(): 1 - wet,       0 - dry

 MASK\_STRUC(): 0 - permanent dry point

 U():  depth-averaged *u*
 
 V():  depth-averaged *v* 
 
 Uc(): contravariant component of depth-averaged velocity in :math:`\xi_1` direction 

 Vc(): contravariant component of depth-averaged velocity in :math:`\xi_2` direction 

 HU(): :math:`(h+\eta)u` at element

 HV(): :math:`(h+\eta)v` at element

 P(): :math:`(h+\eta)U`   at x- interface for Cartesian, and :math:`(h+\eta)Uc`  at :math:`\xi_1` -interface curvilinear

 Q(): :math:`(h+\eta)V`   at y- interface for Cartesian, and :math:`(h+\eta)Vc`  at :math:`\xi_2` -interface curvilinear

 Fx(): numerical flux F at x-interface

 Fy(): numerical flux F at y-interface

 Gx(): numerical flux G at x-interface

 Gy(): numerical flux G at y-interface

 Ubar(): :math:`HU` for Cartesian and :math:`JHU` for curvilinear

 Vbar(): :math:`HV` for Cartesian and :math:`JHV` for curvilinear

 EtaRxL(): :math:`\eta` Left value at x-interface

 EtaRxR(): :math:`\eta` Right value at x-interface

 EtaRyL(): :math:`\eta` Left value at y-interface

 EtaRyR(): :math:`\eta` Right value at y-interface

 HxL():   total depth  Left value at x-interface

 HxR():   total depth  Right value at x-interface

 HyL():   total depth  Left value at y-interface

 HyR():   total depth  Right value at y-interface

 HUxL(): :math:`(h+\eta)u` Left value at x-interface

 HUxR(): :math:`(h+\eta)u` Right value at x-interface

 HVyL(): :math:`(h+\eta)v` Left value at y-interface

 HVyR(): :math:`(h+\eta)v` Right value at y-interface

 PL(): Left P value at x-interface

 PR(): Right P value at x-interface

 QL(): Left Q value at y-interface

 QR(): Right Q value at y-interface




Installation and compilation
####################################

NearCoM-TVD is distributed in a compressed fie. To install the programs, first, uncompress the package. Then use 

  :math:`>` tar xvf :math:`*`.tar 

to extract files from the uncompressed package. The exacted files will be distributed in two new directories: **/CIRC\_SWAN** and **/work**.

To compile the program, go to **/CIRC\_SWAN** and modify Makefile if needed. There are several necessary flags in Makefile needed to specify below.

 -DDOUBLE_PRECISION: use double precision, default is single precision.

 -DPARALLEL: use parallel mode, default is serial mode.

 -DSAMPLES: include all samples, default is no sample included.

 -DCURVILINEAR: curvilinear version, otherwise Cartesian.  \  
   :math:`\bf NOTE:`  setting curvilinear is a must for SWAN and SHORECIRC coupled model.

 -DSEDIMENT: include sediment and seabed modules.

 -DINTEL: INTEL compiler.

 -DRESIDUAL: include tidal residual calculation.

 -DSTATIONARY: stationary mode for SHORECIRC

 CPP: path to CPP directory.

 FC: Fortran compiler. 


Then execute 

 :math:`>` make clean

 :math:`>` make


The executable file 'nearcom' will be generated and  copied from **/CIRC\_SWAN** to **/work/**. Note: use 'make clean' after any modification of Makefile.  

To run the model, go to **/work**. Modify INPUT if needed and run. 






Input
#############

Following are descriptions of parameters in **input.txt**  (:math:`\bf  NOTE:`  all parameter names are capital sensitive).

- **SWAN INPUT:**  refer to SWAN manual. Model run time is set in SWAN model. For example,

  COMPUTE NONSTAT 20081114.160000 1 MI 20081114.230000 

  The above setting means model run start from 2008 11 14 16:00 to 2008 11 14 23:00. The model call swan at :math:`DT_{\mbox{swan}}` = 1 minute. The loop number for SHORECIRC and SEDIMENT is estimated by  :math:`DT_{\mbox{swan}}` and the time step of SHORECIRC (time varying).

  **IMPORTANT SETTING IN SWAN:**

  1. in SET, always set CARTESIAN in order to make a grid orientation consistent with SHORECIRC
    
  2. in SET, always set [inrhog] as 1 to get a true wave energy dissipation.
    
  3. in COMPUTE, always set NONSTAT mode. 



- **WAVE CURRENT INTERACTION**

SWAN\_RUN: logical parameter to run SWAN

SHORECIRC\_RUN: logical parameter to run SHORECIRC 

WC\_BOUND\_WEST:  west bound region  (number of grid point) in which  wave-current is inactive. 

WC\_BOUND\_EAST : east bound region  (number of grid point) in which  wave-current is inactive.

WC\_BOUND\_SOUTH : south bound region  (number of grid point) in which  wave-current is inactive.

WC\_BOUND\_NORTH: north bound region  (number of grid point) in which  wave-current is inactive.

WC\_LAG :  time delay for wave-current interaction


- **TITLE**: 

 title for SHORECIRC log file

      
- **SPECIFICATION OF MULTI-PROCESSORS**

 PX:  processor numbers in X

 PY:  processor numbers in Y    \
    :math:`\bf NOTE:` PX and PY must be consistency with number of processors defined in mpirun command, e.g., mpirun -np n (where n = px :math:`\times` py). 

 
- **SPECIFICATION OF WATER DEPTH**
 
DEPTH\_TYPE: depth input type. 

  .. line-block::

       DEPTH\_TYPE=DATA: from a depth file.    
       The program includes several simple bathymetry configurations such as 
       DEPTH\_TYPE=FLAT:  flat bottom, need DEPTH\_FLAT     
       DEPTH\_TYPE=SLOPE:  plane beach along :math:`x` direction. It needs three parameters: slope,SLP,  slope starting point, Xslp and flat part of depth, DEPTH\_FLAT



DEPTH\_FILE: bathymetry file if  DEPTH\_TYPE=DATA, file dimension should be Mglob x Nglob    \
    with the first point as the south-west corner.  The read format in the code is shown below.

  ::

       DO J=1,Nglob       
        READ(1,*)(Depth(I,J),I=1,Mglob)
       ENDDO
 
DEPTH\_FLAT: water depth of flat bottom if DEPTH\_TYPE=FLAT or DEPTH\_TYPE=SLOPE    \
    (flat part of a plane beach).
 
SLP: slope if DEPTH\_TYPE=SLOPE

Xslp: starting :math:`x` (m) of a slope, if DEPTH\_TYPE=SLOPE


- **SPECIFICATION OF RESULT FOLDER**  
  
RESULT\_FOLDER: result folder name, e.g., RESULT\_FOLDER = /Users/fengyanshi/tmp/

- **SPECIFICATION OF DIMENSION**

Mglob: global dimension in :math:`x` direction.

Nglob: global dimension in :math:`y` direction.    \
     :math:`\bf NOTE:` For parallel runs, Mglob and Nglob can be divided by PX and PY, respectively. MAX(Mglob,Nglob) can be divided by PX :math:`\times` PY.


- **SPECIFICATION OF STATIONARY MODE**

N\_ITERATION: the iteration number for stationary mode of SHORECIRC   \
    (set -DSTATIONARY in Makefile).

WATER\_LEVEL\_FILE:  the file name of water level file containing time and water level, for   \
     stationary mode. The following example shows the format. 

  .. line-block::

        water levels for stationary mode
        5   - number of water level data
        0.0       0.0 ! Time (s), Level (m)
        3600.0    0.5
        7200.0    0.8660 
        10800.0   1.0
        14400.0   0.866 
        18000.0   0.5

- **SPECIFICATION OF TIME**
 
PLOT\_INTV: output interval in seconds (Note, output time is not exact because adaptive dt is used.)

SCREEN\_INTV: time interval (s) of screen print. 

PLOT\_INTV\_STATION: time interval (s) of gauge output


- **SPECIFICATION OF GRID**

DX: grid size(m) in :math:`x` direction, for Cartesian mode

DY: grid size(m) in :math:`y` direction, for Cartesian mode

X\_FILE: name of file to store x for curvilinear mode

Y\_FILE: name of file to store y for curvilinear mode   \
    :math:`\bf NOTE:` data format is the same as the depth data shown above. 

CORI\_CONSTANT: logical parameter for constant Coriolis parameter

LATITUDE: latitude if constant Coriolis parameter is used

LATITUDE\_FILE: name of file to store latitude at every grid point if not constant Coriolis   \
    :math:`\bf NOTE:` data format is the same as the depth data shown above. 


- **BOUNDARY CONDITIONS**

ETA\_CLAMPED: logical parameter for surface elevation clamped condition  

V\_CLAMPED: logical parameter for velocity clamped condition  

FLUX\_CLAMPED: logical parameter for flux clamped condition  

TIDE\_FILE: name of file to store tidal constituents 

     **DATA FORMAT:** please refer to **mk\_tide.f90**. 
     The formula of surface elevation at a tidal boundary can be expressed by


.. math:: \eta_0 (t) =  \sum_{n=1}^Na_{0}({\bf x}, n) f_c (n) \cos \left(\frac{2\pi}{T(n)} t - \phi({\bf x}, n)  + (V_0 +u_0)(n) \right)

where :math:`a_0`  and :math:`\phi` represent amplitude and phase lag, respectively, for a harmonic constituent at location :math:`\bf x`. :math:`T` is tidal period. :math:`f_c` and :math:`(V_0+u_0)` are the lunar node factor and the equilibrium argument, respectively, for a constituent. 

  .. line-block::

      The following is an example of M2 + O1.

  .. line-block::

      tidal boundary conditions
      150 --- number of days from Jan 1,  to simulation date 
      2 ---  number of constituents 
      1.000       0.000  --- :math:`f_c` and :math:`(V_0+u_0)` for M2
      0.980       0.000  --- :math:`f_c` and :math:`(V_0+u_0)` for O1
      80   ---  number of tidal boundary points
      1 , 1   ---  (i,j) grid location of tidal boundary  
      12.420   1.200   21.000 --- :math:`T`, amplitude :math:`a_0` and phase lag :math:`\phi` for M2 
      24.000   0.3  30.100 -- :math:`T`, amplitude :math:`a_0` and phase lag :math:`\phi` for O1 
      2 , 1  ---  (i,j) grid location of tidal boundary 
      12.420   1.200  21.000   --- :math:`T`, amplitude :math:`a_0` and phase lag :math:`\phi` for M2
      24.000   0.3   30.100 -- :math:`T`, amplitude :math:`a_0` and phase lag :math:`\phi` for O1
      3 ,  1
      ...
   

FLUX\_FILE: name of file to store time series of flux (e.g., unit width river flux) 

  .. line-block::

     **DATA FORMAT:**
     title
     Number of data, Number of flux point
     I, J, River orientation 
     Time, Flux, Angle in Cartesian 
     ...  
     where  (I,J) represent grid points of river location. River orientation represents the direction which a river flows from  in the  IMAGE domain (for curvilinear coordinates). Use W,E,S and N for the orientation.  For example, 'W' represents a river flowing into the domain from the west boundary (in IMAGE domain for curvilinear coordinates).

     Please refer to **mk\_river.f90**. The following is an example.

  .. line-block::

     river flux boundary condition \
     5       2     ! NumTimeData, NumFluxPoint \
     1  38  W      ! I, J, River\_Orientation\
     0.000       0.200       0.000\
     360000.000       0.200       0.000\
     720000.000       0.200       0.000\
     1080000.000       0.200       0.000\
     1440000.000       0.200       0.000\
     1  39  W      ! I, J, River\_Orientation\
     0.000       0.200       0.000\
     360000.000       0.200       0.000\
     720000.000       0.200       0.000\
     1080000.000       0.200       0.000\
     1440000.000       0.200       0.000\
     end of file\


- **WIND CONDITION**

  Spatially uniform wind field is assumed in this version.  

WindForce: logical parameter for wind condition, T or F. 

WIND\_FILE: name of file for a time series of wind speed.

  **DATA FORMAT:** the following is an example of wind data.

  .. line-block::

     wind data
     100  - number of data
     0.0 ,    -10.0 0.0   ---  time(s), wu, wv (m/s)
     2000.0,   -10.0,  0.0
     8000.0,  -10.0,   0.0
     ... 

Cdw: wind stress coefficient for the quadratic formula. 


- **SPECIFICATION OF INITIAL CONDITION**
 
INT\_UVZ : logical parameter for initial condition, default is FALSE
 
 
ETA\_FILE: name of file for initial :math:`\eta`, e.g., ETA\_FILE= /Users/fengyanshi/work/input/CVV\_H.grd, data format is the same as depth data.

U\_FILE:  name of file for initial :math:`u`, e.g.,U\_FILE= /Users/fengyanshi/work/input/CVV\_U.grd, data format is the same as depth data.

V\_FILE:  name of file for initial :math:`v`, e.g., V\_FILE= /Users/fengyanshi/work/input/CVV\_V.grd, data format is the same as depth data.


- **SPECIFICATION OF WAVEMAKER**
 
 There is no wavemaker implemented in SHORECIRC.


- **SPECIFICATION OF PERIODIC BOUNDARY CONDITION**

  (Note: only south-north periodic condition was implemented)

PERIODIC\_X: logical parameter for periodic boundary condition in x direction, T - periodic, F - wall boundary condition.

PERIODIC\_Y: logical parameter for periodic boundary condition in x direction.

Num\_Transit: grid numbers needed to make periodic condition for SWAN. The reason to set this parameter is that SWAN doesn't have an option for periodic boundary condition. In this implementation, a periodic boundary condition is implemented by making a transition from  a left array ( count to Num\_Transit from left boundary) to a right array. 


- **SPECIFICATION OF SPONGE LAYER**
 
SPONGE\_ON: logical parameter, T - sponge layer, F - no sponge layer.
 
Sponge\_west\_width: width (m) of sponge layer at west boundary.

Sponge\_east\_width:   width (m) of sponge layer at east boundary.

Sponge\_south\_width: width (m) of sponge layer at south boundary.

Sponge\_north\_width width (m) of sponge layer at north boundary

R\_sponge: decay rate in sponge layer. Its values are between 0.85 :math:`\sim` 0.95.

A\_sponge: maximum damping magnitude. The value is :math:`\sim` 5.0. 


- **SPECIFICATION OF OBSTACLES**

OBSTACLE\_FILE: name of obstacle file. 1 - water point, 0 - permanent dry point. Data dimension is (Mglob :math:`\times` Nglob). Data format is the same as the depth data. 
 

- **SPECIFICATION OF PHYSICS** 
  
Cd: quadratic bottom friction coefficient 

nu\_bkgd : background eddy viscosity parameter. 
  

- **SPECIFICATION OF NUMERICS**

Time\_Scheme: stepping option,  Runge\_Kutta or Predictor\_Corrector (not suggested for this version).

HIGH\_ORDER: spatial scheme option,  FOURTH for the fourth-order, THIRD for the third-order, and SECOND for the second-order (not suggested for Boussinesq modeling).

CONSTRUCTION: construction method,  HLL for HLL scheme, otherwise for averaging scheme.

CFL: CFL number, CFL :math:`\sim` 0.5.

FroudeCap: cap for Froude number in velocity calculation for efficiency. The value could be 5 :math:`\sim` 10.0.

MinDepth: minimum water depth (m) for wetting and drying scheme. Suggestion: MinDepth = 0.001 for lab scale and 0.01 for field scale. 

MinDepthFrc: minimum water depth (m) to limit bottom friction value. Suggestion: MinDepthFrc = 0.01 for lab scale and 0.1 for field scale. 


- **SPECIFICATION OF TIDAL RESIDUAL**

T\_INTV\_mean: time-averaging interval for Eulerian mean current and elevation. Note: use -DRESIDUAL in Makefile to make this option active.


- **SPECIFICATION OF SEDIMENT CALCULATION**

  Note: set -DSEDIMENT in Makefile to make sediment module active

T\_INTV\_sed: time interval to call sediment module

Factor\_Morpho: morphology factor.

D\_50 : :math:`D_{50}`

D\_90 : :math:`D_{90}`

por:  sediment porosity

RHO: water density

nu\_water: water eddy viscosity

S\_sed: specific gravity

SOULSBY: logical parameter for Soulsby (1997) total load formula, T = true, F = false 

z0: :math:`z_0`, bed roughness length. 

KOBAYASHI: logical parameter for KOBAYASHI's formula, T = true, F = false

angle\_x\_beach: coordinate rotation angle defined in Figure 1.

eB: :math:`e_B`,  suspension efficiency for energy dissipation rate due to wave breaking

ef: :math:`e_f`, suspension efficiency for energy dissipation rate due to bottom friction

a\_k: :math:`a`, empirical suspended load parameter. 

b\_k: :math:`b`, empirical bedload parameter.

TanPhi: :math:`\tan \phi`, where :math:`\phi` is the angle of internal friction of the sediment. 

Gm: :math:`G_m` for slope function (:math:`G_m=10`).

frc: friction coefficient in Kobayashi.

Si\_c: a coefficient in calculating :math:`P_b`.


- **SPECIFICATION OF OUTPUT VARIABLES**

NumberStations: number of station for output. If NumberStations :math:`> 0`, need input i,j in STATION\_FILE
 
DEPTH\_OUT: logical parameter for output depth. T or F. 

U: logical parameter for output :math:`u`. T or F. 

V: logical parameter for output :math:`v`. T or F. 

ETA: logical parameter for output :math:`\eta`. T or F. 

HS: logical parameter for output of significant wave height :math:`H_s`. T or F. 

WFC: logical parameter for output of wave force. T or F. 

WDIR: logical parameter for output of peak wave direction. T or F. 

WBV: logical parameter for output of wave orbital velocity. T or F. 

MASK: logical parameter for output wetting-drying MASK. T or F. 

SourceX: logical parameter for output source terms in :math:`x` direction. T or F. 

SourceY:  logical parameter for output source terms in :math:`y` direction. T or F. 

UV3D: logical parameter for output 3D structure. T or F.

Qstk:   logical parameter for output Stokes mass flux. T or F.

DepDt:  logical parameter for output depth variation rate. T or F.

Qsed:  logical parameter for output sediment transport rate. T or F.



Output
######################
The output files are saved in the result directory defined by RESULT\_FOLDER in INPUT. For outputs in ASCII,  a file name is a combination of variable name and an output series number such eta\_0001, eta\_0002, .... The format  and read/write algorithm are  consistent with a depth file.  Output for stations is a series of numbered files such as sta\_0001, sta\_0002 .... 




