*****************************************
**Numerical schemes for SHORECIRC**
*****************************************

Compact form of governing equations
#########################################

The compact form of governing equations can be written as

.. math:: \frac{\partial {\bf \Psi}}{\partial t} + \nabla \cdot {\bf \Theta}({\bf \Psi}) = {\bf S}
   :label: tvd

where :math:`{\bf \Psi}` and :math:`{\bf \Theta} ({\bf \Psi})` are  the vector of conserved variables and  the flux vector function, respectively, and are given by

.. math:: {\bf \Psi} = \left( \begin{array}{c} J\eta \\
          JP \\
          JQ \end{array} \right) 

.. math:: {\bf \Theta} = \left( \begin{array}{c} JP {\bf i} + JQ {\bf j} \\
          \left[JPu +\frac{1}{2}g (\eta^2 + 2 \eta h)JL^1_1 \right ] {\bf i} + \left[JQu +\frac{1}{2}g (\eta^2 + 2 \eta h)JL^2_1 \right ] {\bf j} \\
          \left[JPv +\frac{1}{2}g (\eta^2 + 2 \eta h)JL^1_2 \right ]  {\bf i} + \left[JQv +\frac {1}{2}g (\eta^2 + 2 \eta h)JL^2_2 \right ] {\bf j} \end{array}\right) .

.. math:: {\bf S} = \left( \begin{array}{c} 0\\
          \frac{1}{\rho} \left(F_{\mbox{wx}}  -\tau^b_x + \tau^s_x -f_x,  - ROT_x \right) \\
          \frac{1}{\rho} \left(F_{\mbox{wy}}  -\tau^b_y + \tau^s_y  -f_y - ROT_y \right) \end{array} \right)

where :math:`(u,v)` represent Cartesian components of velocity, :math:`(P,Q) = (Hu^1, Hu^2)`, the contravariant components of flow flux, (:math:`\tau^b_x, \tau^b_y`), (:math:`\tau^b_x, \tau^b_y`), (:math:`f_x, f_y`) and (:math:`ROT_x, ROT_y`) represent the bottom friction, wind stress,  Coriolis forcing and the rest of terms in :math:`x` and :math:`y` directions, respectively. 
(:math:`F_{\mbox{wx}}, F_{\mbox{wy}}`) are the gradient of radiation stresses expressed by

.. math:: F_{\mbox{wx}} = -  \left [ \frac{\partial }{\partial \xi^1} \left (S_{xx} J L^1_1+S_{xy} J L^1_2 \right) +  [ \frac{\partial }{\partial \xi^2} \left (S_{xx} J L^2_1+S_{xy} J L^2_2 \right)
 \right ]

.. math:: F_{\mbox{wy}} = -  \left [ \frac{\partial }{\partial \xi^1} \left (S_{xy} J L^1_1+S_{yy} J L^1_2 \right) + [ \frac{\partial }{\partial \xi^2} \left (S_{xy} J L^2_1+S_{yy} J L^2_2 \right) \right ]

The reason to put the wave forcing as a source term rather than a flux term is that the forcing  is calculated directly from SWAN in a curvilinear grid. 





Spatial discretization
###############################

A combined finite-volume and finite-difference method was applied to the spatial discretization. 
For the flux terms and the first-order derivative terms, MUSCL-TVD schemes from the 2nd-order to 4th-order  are implemented in the present model. The multi-order MUSCL-TVD scheme can be written in a compact form  according to Erduran et al. (2005) who modified Yamamoto et al.'s (1995)  fourth-order approach. 
In :math:`x` -direction, for example, the combined form of the interface construction can be written as follows:

.. math:: \phi^L_{i+1/2} = \phi_i +\frac{1}{4} \left[ (1-\kappa_1)\chi(r) \Delta^* \phi_{i-1/2} + (1+\kappa_1) \chi(1/r) \Delta^* \phi_{i+1/2} \right ]
   :label: left

.. math:: \phi^R_{i-1/2} = \phi_i -\frac{1}{4} \left[ (1+\kappa_1)\chi(r) \Delta^* \phi_{i-1/2} + (1-\kappa_1) \chi(1/r) \Delta^* \phi_{i+1/2} \right ]
   :label: right

where :math:`\phi^L_{i+1/2}` is the constructed value at the left-hand side of the interface :math:`i+\frac{1}{2}` and :math:`\phi^R_{i-1/2}` is the value at the right-hand side of the interface :math:`i - \frac{1}{2}`.  
The values of :math:`\Delta^* \phi` are evaluated as follows:

.. math::
   \begin{array}{l l l}
      \Delta^* \phi_{i+1/2} = \Delta \phi_{i+1/2} - \kappa_2 \Delta^3 \bar{\phi}_{i+1/2} /6, \\
      \Delta \phi_{i+1/2} = \phi_{i+1} - \phi_{i}, \\
      \Delta^3 \bar{\phi}_{i+1/2} = \Delta \bar{\phi}_{i+3/2} - 2 \Delta \bar{\phi}_{i+1/2} + \Delta \bar{\phi}_{i-1/2}, \\
      \Delta \bar{\phi}_{i-1/2} =\mbox{minmod} (\Delta \phi_{i-1/2}, \Delta \phi_{i+1/2}, \Delta \phi_{i+3/2}), \\
      \Delta \bar{\phi}_{i+1/2} =\mbox{minmod} (\Delta \phi_{i+1/2}, \Delta \phi_{i+3/2}, \Delta \phi_{i-1/2}), \\
      \Delta \bar{\phi}_{i+3/2} =\mbox{minmod} (\Delta \phi_{i+3/2}, \Delta \phi_{i-1/2}, \Delta \phi_{i+1/2}) 
   \end{array}
   :label: minmod

In :eq:`minmod`, minmod represents the minmod limiter and is given by

.. math:: \mbox{minmod} (j,k,l) = \mbox{sign} (j) \mbox{max} \{ 0, \mbox{min} [|j|, 2 \mbox{sign} (j) k, 2 \mbox{sign} (j) l ] \} .

:math:`\kappa_1` and :math:`\kappa_2` in :eq:`left` and :eq:`right` are  control parameters for orders of the scheme in the compact form.
The complete form with :math:`(\kappa_1, \kappa_2) = (1/3, 1)` is the fourth-order scheme given by Yamamoto et al. (1995). 
:math:`(\kappa_1, \kappa_2) = (1/3, 0)` yields a third-order scheme, while the second-order scheme can be retrieved using :math:`(\kappa_1, \kappa_2) = (-1, 0)`. 

:math:`\chi(r)` in :eq:`left` and :eq:`right` is the limiter function. 
The original scheme introduced by Yamamoto et al. (1998) uses the Minmod limiter as used in :eq:`minmod`. Erduran et al. (2005) found that the use of the van-Leer limiter for the third-order scheme gives more accurate results. 
Their finding was confirmed by  using the present model in the benchmark tests for wave runup conducted by Tehranirad et al. (2011).  
The van-Leer limiter can be expressed as

.. math:: \chi(r) = \frac{r+|r|}{1+r}

where

.. math:: r = \frac{\Delta^* \phi_{i+1/2}}{\Delta ^* \phi_{i-1/2}}.

The numerical fluxes are computed using a HLL approximate Riemann solver

.. math:: {\bf \Theta} ({\bf \Psi}^L, {\bf \Psi}^R)= \left \{ \begin{array}{ll} {\bf \Theta} ({\bf \Psi}^L) & \mbox{if} \ \ \ s_L \ge 0 \\
          {\bf \Theta}^* ({\bf \Psi}^L, {\bf \Psi}^R) & \mbox{if} \ \ \ s_L < 0 < s_R \\
          {\bf \Theta}({\bf \Psi}^R) & \mbox{if} \ \ \ s_R \le 0,\end{array}  \right.

where

.. math:: {\bf \Theta^*} ({\bf \Psi}^L, {\bf \Psi}^R) = \frac{s_R {\bf \Theta} ({\bf \Psi}^L) -s_L {\bf \Theta}({\bf \Psi}^R) + s_L s_R ({\bf \Psi}^R - {\bf \Psi}^L)}{s_R - s_L}

The wave speeds of the Riemann solver are given by

.. math:: s_L= \mbox{min} ({\bf V}^L \cdot {\bf n} - \sqrt{g (h+\eta)^L} , u_s -\sqrt{\varphi_s}),  
   :label: sl

.. math:: s_R= \mbox{max} ({\bf V}^R \cdot {\bf n} + \sqrt{g (h+\eta)^R},  u_s +\sqrt{\varphi_s}),  
   :label: sr

in which :math:`u_s` and :math:`\varphi_s` are estimated as

.. math:: u_s =\frac{1}{2} ({\bf V}^L + {\bf V}^R)\cdot {\bf n} + \sqrt{g (\eta + h)^L} - \sqrt{g(\eta+h)^R}

.. math:: \sqrt{\varphi_s} = \frac{\sqrt{g (\eta + h)^L}+\sqrt{g (\eta + h)^R}}{2} +\frac{({\bf V}^L - {\bf V}^R)\cdot {\bf n}}{4}

and :math:`{\bf n}` is the normalized side vector for a cell face.





Time stepping
##########################

The third-order Strong Stability-Preserving (SSP) Runge-Kutta scheme for nonlinear spatial discretization (Gottlieb et al., 2001) was adopted for time stepping. The scheme is given by

.. math::
   \begin{array}{l l l}
      {\bf \Psi}^{(1)} = {\bf \Psi}^{n}  + \Delta t (- \nabla \cdot {\bf \Theta} ({\bf \Psi}^n) + {\bf S}^{(1)} ) \\
      {\bf \Psi}^{(2)} = \frac{3}{4}{\bf \Psi}^{n}  + \frac{1}{4} \left[   {\bf \Psi}^{(1)} +  \Delta t \left (- \nabla \cdot {\bf \Theta} ({\bf \Psi}^{(1)} ) + {\bf S}^{(2)} \right) \right] \\
      {\bf \Psi}^{n+1} = \frac{1}{3}{\bf \Psi}^{n}  + \frac{2}{3} \left[   {\bf \Psi}^{(2)} +  \Delta t \left (- \nabla \cdot {\bf \Theta} ({\bf \Psi}^{(2)} ) + {\bf S}^{n+1} \right) \right]
   \end{array}
   :label: runge

in which :math:`{\bf \Psi}^{n}` denotes :math:`{\bf \Psi}` at time level :math:`n`. 
:math:`{\bf \Psi}^{(1)}` and :math:`{\bf \Psi}^{(2)}` are values at intermediate stages in the Runge-Kutta integration. 

  An adaptive time step is chosen, following the Courant-Friedrichs-Lewy (CFL) criterion:

.. math:: \Delta t = C  \mbox{min} \left ( \mbox{min} \frac{\Delta x}{|u_{i,j}| + \sqrt{g (h_{i,j} +\eta_{i,j})}},  \mbox{min} \frac{\Delta y}{|v_{i,j}| + \sqrt{g (h_{i,j} +\eta_{i,j})}} \right )
   :label: cfl

where :math:`C` is the Courant number and :math:`C=0.5` was used in the following examples.  





Wetting-drying schemes for shallow water
##############################################

The wetting-drying scheme for modeling a moving boundary is straightforward. 
The normal flux :math:`{\bf n} \cdot {\bf M}` at the cell interface of a dry cell is set to zero. 
A mirror boundary condition is applied to the high-order MUSCL-TVD scheme. 
It may be noted that the wave speeds of the Riemann solver  :eq:`sl` and :eq:`sr` for a dry cell are modified as 


.. math:: s_L= {\bf V}^L \cdot {\bf n} - \sqrt{g (h+\eta)^L} ,  \ \ \ \    s_R= {\bf V}^L \cdot {\bf n} +2  \sqrt{g (h+\eta)^L}  \  \  \  \mbox{(right dry cell)}

and

.. math:: s_L= {\bf V}^R \cdot {\bf n} - \sqrt{g (h+\eta)^R} ,  \ \ \ \    s_R= {\bf V}^R \cdot {\bf n} +2  \sqrt{g (h+\eta)^R}  \  \  \  \mbox{(left dry cell)}





Boundary conditions
##############################

To incorporate tide and river inflow into the circulation module, two types of open boundary conditions are implemented. 
One is the surface clamped boundary condition that reads at tidal open boundaries            

.. math:: \eta= \eta_0    \ \ \ \ \ \ \mbox{at tidal open boundaries}

where :math:`\eta_0` is measured or predicted surface elevations at open boundaries. 
The other is the specified flux boundary condition that is usually used at river boundaries:

.. math:: P^\alpha = Q^{\mbox{flux}}_\beta L^\alpha_\beta

where :math:`Q^{\mbox{flux}}` is the specified volume flux per unit width. 




