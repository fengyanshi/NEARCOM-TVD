*************************
**Formulations**
*************************

SHORECIRC equations
#########################

SHORECIRC is a quasi-3D nearshore circulation model for prediction of wave-induced nearshore circulation.
It is a 2D horizontal model which incorporates the effect of the vertical structure of horizontal flows.
In Putrevu and Svendsen (1999), the instantaneous horizontal velocity is split as

.. math:: u^{ins}_{\alpha}=u^{\prime}_{\alpha}+u_{w\alpha}+u_{\alpha}+u_{1\alpha} 

where :math:`u^{\prime}_{\alpha}, u_{w\alpha}, u_{\alpha}` and :math:`u_{1 \alpha}` are, respectively, the turbulence component, the wave component, the component of depth-averaged and short-wave-averaged velocity, and the vertical variation of the short-wave-averaged velocity. 
The complete SHORECIRC equations can be expressed in the Cartesian form as in Svendsen et al. (2003) or in the curvilinear form in Shi et al. (2003). In this section, we derive a conservative form of the equations in generalized curvilinear coordinates in order to use the TVD numerical scheme. 

SHORECIRC equations in Cartesian coordinates :math:`(x_1, x_2)` can be expressed as

.. math:: \frac{\partial \eta}{\partial t}+\frac{\partial Hu_\alpha}{\partial x_\alpha} =0 

.. math:: \frac{\partial H u_\alpha}{\partial t} + \frac{Hu_\alpha u_\beta}{\partial x_\beta} +f_{\alpha}+ g H \frac{\partial \eta}{\partial x_\alpha} + \frac{1}{\rho} \frac{\partial T_{\alpha \beta}}{\partial x_\beta} + \frac{1}{\rho} \frac{\partial S_{\alpha \beta}}{\partial x_\beta}  + \frac{\tau^b_\alpha}{\rho} - \frac{\tau^s_\alpha}{\rho} + \mbox{ROT}= 0  
   :label: euler

where :math:`\eta` represents the wave-averaged surface elevation, :math:`H = \eta + h`, in which :math:`h` is still water depth. The depth-averaged and wave-averaged velocity :math:`u_\alpha` is defined by the 'Lagrangian' averaging as

.. math:: u_{\alpha} = \frac{1}{H} \overline{\int_{-h}^\zeta u^{ins}_\alpha}

where :math:`\zeta` is the instantaneous surface elevation. This split is different from Haas et al. (2003) who used the 'Eulerian' concept for their split. For the Lagrangian averaging method, it is assumed that

.. math:: {\int}_{-h}^\eta u_{1\alpha} dz = - Q_{w \alpha}

where :math:`Q_{w \alpha}` is the short wave flux. 

In :eq:`euler`, :math:`f_\alpha` represents the Coriolis force which was not considered in the previous SHORECIRC equations (Putrevu and Svendsen, 1999). 
:math:`T_{\alpha \beta}, S_{\alpha \beta}, {\tau}^{s}_{\alpha}` and  :math:`{\tau}^{b}_{\alpha}` are the depth-integrated Reynolds' stress, the wave-induced radiation stress (Longuet-Higgins and Stewart, 1962, 1964), the surface shear stress and the bottom shear stress, respectively.   
ROT represents the rest of terms associated with 3D dispersion terms which are not presented here because they do not involve the coordinate transformation conducted below. 
Interested readers are referred to Putrevu and Svendsen (1999).

A curvilinear coordinate transformation is introduced in the general form

.. math:: xi^1 = \xi^1 (x_1,x_2), \ \ \ \   \xi^2 = \xi^2 (x_1,x_2)

where :math:`(\xi^1, \xi^2)` are the curvilinear coordinates. 
We use superscript indices, i.e., :math:`()^\alpha`, to represent the contravariant component and subscript indices for the Cartesian component  of a vector.
The relation between the Cartesian component :math:`u_\alpha` and the contravarient component :math:`u^\alpha` can be written, by definition, as

.. math:: u^\alpha = u_\beta L^\alpha_\beta

where

.. math:: L^\alpha_\beta = \frac{\partial \xi^\alpha}{\partial x_\beta}

Using the chain rule, the derivative of a function :math:`F` with respect to :math:`x_\alpha` in the Cartesian coordinates can be expressed in the curvilinear coordinates :math:`\xi^\alpha` by

.. math:: \frac{\partial F}{\partial x_\alpha} = L^\beta_\alpha \frac{\partial F}{\partial \xi^\beta}
   :label: der

Using the metric identity law (Thompson et al., 1985):

.. math:: \frac{\partial}{\partial \xi^\alpha}(J L^\alpha_\beta) \equiv 0

:eq:`der` can be rewritten in a different form as

.. math:: \frac{\partial F}{\partial x_\alpha} = \frac{1}{J} \frac{\partial F J L^\beta_\alpha}{\partial \xi^\beta}
   :label: der2

Using :eq:`der2`, the conservative form of SHORECIRC equations in curvilinear coordinates can be derived as

.. math:: frac{\partial \eta}{\partial t} + \frac{1}{J}\frac{\partial J P^\alpha}{\partial \xi^\alpha} = 0
   :label: mass1

.. math:: & \frac{\partial Hu_\alpha}{\partial t} + \frac{1}{J} \frac{\partial }{\partial \xi^\beta} \left[ JP^\beta u_\alpha  +\frac{1}{2} g (\eta^2 +2 \eta h)  J L^\beta_\alpha \right] +f_\alpha -  g \eta \frac{1}{J}  \frac{\partial }{\partial \xi^\beta} (h J L^\beta_\alpha) \\
          & + \frac{1}{\rho} \frac{1}{J}  \frac{\partial }{\partial \xi^\gamma} (S_{\alpha \beta} J L^\gamma_\beta) + \frac{1}{J}  \frac{\partial }{\partial \xi^\gamma} (\tau_{\alpha \beta} JHL^\gamma_\beta)+ \frac{\tau^b_\alpha}{\rho} - \frac{\tau^s_\alpha}{\rho} + \mbox{ROT} = 0 
   :label: conserv1

Note that  equations :eq:`mass1` and :eq:`conserv1` contain both the Cartesian and contravariant variables, which are different from the contravariant-only form in Shi et al., (2003). 
A difficulty in deriving the conservative form of the momentum equations using contravariant-only variables has be noticed by several authors (e.g., xxx) in shallow water applications.   
Here, we follow Shi and Sun (1995) who used both Cartesian and contravariant variables in the derivation of the momentum equations. 

In :eq:`mass1` and :eq:`conserv1`, :math:`P^\alpha = H u^\alpha`, denotes the contravariant component of volume flux, the Coriolis force :math:`f_\alpha` uses the Cartesian components, i.e., :math:`(-f H v,  fH u)`, :math:`S_{\alpha \beta}` represents the Cartesian component of radiation stress. 
In the present application, the gradient of radiation stresses (the fifth term of :eq:`conserv1`) is obtained directly from the wave module SWAN.
:math:`\tau_{\alpha \beta}` represents the Cartesian component of turbulent shear stress, :math:`\tau^b_\alpha` and :math:`\tau^s_\alpha` are the Cartesian components of bottom stress and wind stress.

In the derivation of the momentum equations, the surface gradient term is treated following Shi et al. (2011) who reorganized this term in order to make it well-balanced in a MUSCL-TVD scheme for a general order.
The  expression in curvilinear coordinates can be written as

.. math:: -gHJ\frac{\partial \eta}{\partial x_\alpha} = -\frac{\partial }{\partial \xi^\alpha} \left [\frac{1}{2} g (\eta^2 +2 \eta h) J L^\alpha_\beta \right] + g \eta \frac{\partial }{\partial \xi^\alpha} (h J L^\alpha_\beta)

There are several advantages in using equations :eq:`mass1` and :eq:`conserv1` versus the contravariant-only form in Shi et al. (2003). 
First, they are a conservative form which can be implemented  using a hybrid numerical scheme. Second, all forcing terms resume to a vector form  in Cartesian coordinates in contrast to the contravariant form in Shi et al. (2003). 
Third, the radiation stress term :math:`S_{\alpha \beta}` uses the original form defined in Cartesian coordinates, thus that  there is no need to make a transformation for the second-order tensor.  
The disadvantage is that :eq:`conserv1` contains both the Cartesian component :math:`u_\alpha` and the contravariant component :math:`P^\alpha`. 
However, in terms of the hybrid numerical  scheme used in the study, it is convenient to solve both of the two variables using the explicit numerical scheme, rather than  the implicit numerical scheme used  in Shi et al., (2007).

Wind stress in SHORECIRC was computed using Van Dorn's (1953) formula:
  
.. math:: \tau^s_\alpha = f_a \rho_a |{\bf W}| W_{\alpha} 
    
where :math:`{\bf W}` is wind speed at a 10 :math:`m` elevation above the water surface, :math:`f_a` is the drag coefficient which can be found in Dean and Dalrymple (1991). :math:`\rho_a` represents air density. 
For the bottom stress in SHORECRC, we adopted the formulation for wave-averaged bottom stress for combined currents and waves given by Svendsen and Putrevu (1990) which is written as

.. math:: \tau^b_\alpha = f_{cw} \rho u_0 (\beta_1 u_{b \alpha} +\beta_2 U_{w \alpha}) 

where :math:`U_{w \alpha}` is the amplitude of short-wave particle velocity evaluated at the bottom using wave bulk parameters, :math:`f_{cw}` is the friction factor, :math:`u_{b\alpha}` is the current velocity at the bottom obtained from the theoretical solution of the equation for the vertical variation, :math:`\beta_1` and :math:`\beta_2` are  the weight factors for the current and wave motion given by Svendsen and Putrevu (1990) and evaluated using  linear wave theory. 

The equation governing the vertical structure of horizontal velocity can be solved analytically using the lowest order of the equation for vertical variation of current:

.. math:: \frac{\partial u_{1\alpha}}{\partial t} -\frac{\partial }{\partial z} \left( \nu_t \frac{\partial u_{1\alpha}}{\partial z}\right) = F_\alpha
   :label: vertical

where :math:`\nu_t` is the eddy viscosity coefficient; and :math:`F_\alpha` is a general form of the local forcing described as

.. math:: F_\alpha = \frac{1}{\rho h} f_{w\alpha} - f^{lrad}_\alpha + \frac{\tau^b_\alpha}{\rho h} - \frac{\tau^s_\alpha}{\rho h}
   :label: f1

where :math:`f_\alpha^{lrad}` represents the local radiation stress defined by Putrevu and Svendsen (1999). In the offshore domain without wave forcing, :eq:`f1` reduces to

.. math:: F^\alpha =  \frac{\tau^b_\alpha}{\rho h} - \frac{\tau^s_\alpha}{\rho h}
   :label: f2
            
The solution of equation :eq:`vertical` can be solved analytically following Putrevu and Svendsen (1999). 
The bottom current velocity :math:`u_{b \alpha}` can be evaluated using :math:`u_\alpha` and :math:`u_{1\alpha}` at bottom:

.. math:: u_{b \alpha} = u_\alpha + u_{1\alpha} (z=-h)





SWAN equations
#########################

In this section, we briefly summarize equations used in SWAN. Detailed documentation is referred to SWAN Users' Manual at http://swanmodel.sourceforge.net/. 

The spectral wave model SWAN (Booij et al., 1999, Ris et al., 1999) solves the wave action balance equation. In generalized curvilinear coordinates, the governing equation reads

.. math:: \frac{\partial N}{\partial t} + \frac{1}{J}\frac{\partial (JC_g^\alpha)}{\partial \xi_\alpha} + \frac{\partial (C_{g\sigma}N)}{\partial \sigma} + \frac{\partial (C_{g\theta}N)}{\partial \theta} = \frac{S}{\sigma}
   :label: eq_action

where :math:`\xi_\alpha` represents curvilinear coordinates defined as the same as in the curvilinear SHORECIRC equations; 
:math:`\sigma` is the relative angular frequency; :math:`\theta` is propagation direction of each wave component; :math:`C_g^\alpha` represents contravariant component of the energy propagation speed  which can be obtained using coordinate transformation:

.. math:: C_g^\alpha = C_{g\beta} L^\alpha_\beta

in which :math:`C_{g_\beta} = (C_{gx}, C_{gy})` in the rectangular Cartesian coordinates.  
Equation :eq:`eq_action` is written in a tensor-invariant form in order to make it consistent with the circulation equations. The extended numerical form  can be found in Booij et al. (1997).
:math:`C_{g\sigma}` and :math:`C_{g \theta}` denote energy propagation speeds in the  :math:`\sigma` and :math:`\theta`-spaces, respectively; :math:`S` represents source and sink term in terms of energy density representing the effects of wave generation, dissipation and nonlinear wave-wave interactions; :math:`N` is wave action defined by

.. math:: N = E (\xi^\alpha, \sigma,\theta, t) / \sigma

in which :math:`E` is wave energy density. 

Current effects on wave transformation are presented by the following calculations:

1.total group velocity including the current component

.. math:: C_{g\alpha} = \frac{1}{2}\left( 1+ \frac{2kd}{\sinh2kd}\right)\frac{\sigma k_\alpha}{|{\bf k}|^2} + u_{E\alpha}

where :math:`k_\alpha` or :math:`{\bf k}` represent wave number, :math:`d` is short-wave-averaged water depth, and :math:`d = h + \eta`.  
:math:`u_{E\alpha} = u_\alpha-Q_{w\alpha}/H` in terms of undertow or Eulerian mean velocity. 

2.changes to relative frequency

.. math:: C_\sigma = \frac{\partial \sigma}{\partial d}\left( \frac{\partial d}{\partial t} + {\bf u_E} \cdot \nabla d\right) -C_g \bf{k}\cdot \frac{\partial {\bf u_E}} {\partial {\it s}}

3.wave refraction by current included in

.. math:: C_\theta = - \frac{1}{k} \left( \frac{\partial \sigma}{\partial d} \frac{\partial d}{\partial m} +\bf{k}\cdot \frac{\partial {\bf u_E} }{\partial  {\it m} } \right)

where :math:`s` is the space coordinate in direction :math:`\theta` and :math:`m` is a coordinate normal to :math:`s`.





Soulsby's sediment transport formula
##############################################

The total load transport by waves plus currents derived by  Soulsby's (1997) can be written as

.. math:: q_\alpha = A_s z u_\alpha \left[  \left( u_\alpha^2 + \frac{0.018}{C_d} u_{rms}^2 \right)^{1/2} - u_\alpha^{cr} \right]^{2.4} (1-1.6\tan \beta)

where   

.. math:: A_{sb} = \frac{0.005 h (d_{50}/d)^{1.2}}{\left[ (s-1)g d_{50}\right]^{1.2}}

.. math:: A_{ss} = \frac{0.012 d_{50}D_*^{-0.6}}{\left[ (s-1)g d_{50}\right]^{1.2}}

.. math:: A_s = A_{sb} + A_{ss}

:math:`u_{rms}` is root-mean-square wave orbital velocity; :math:`C_d` is drag coefficient due to current alone; :math:`u_\alpha^{cr}` represents threshold current velocity given by Van Rijn (1984); :math:`\beta` is the slope of bed in streamwise direction; :math:`s` is the specific gravity; and

.. math:: D_* = \left[ \frac{g(s-1)}{\nu} \right]^{1/3} d_{50}






Kobayashi's sediment transport formula
##############################################

Kobayashi et al. (2007) recently proposed sediment transport formulas for the cross-shore and longshore transport rates of suspended sand and bedload on beaches based on laboratory experiment and field data. 
For suspended sand, the suspended sediment volume :math:`V_s` per unit horizontal bottom area was introduced as

.. math:: V_s = P_s \frac{e_B D_r + e_f D_f}{\rho g (s-1) w_f} (1+S_b^2)^{0.5}

where :math:`S_b` is cross-shore bottom slope; :math:`e_B` and :math:`e_f` are suspension efficiencies for energy dissipation rates :math:`D_r` and :math:`D_f` due to wave breaking and bottom friction, respectively. 
:math:`w_f` is the sediment fall velocity; :math:`P_s` is the probability of sediment suspension formulated in Kobayashi et al. (2007). 
The corresponding cross-shore and alongshore suspended sediment transport rates may be expressed as

.. math:: q_{sc} = a \bar{U}_c V_s; \ \ \ \ \ \ q_{sl} = \bar{V}_l V_s
   :label: qs

where :math:`a` is empirical suspended load parameter; :math:`\bar{U}_c` is the Eulerian-mean cross-shore current which consists of wave-induced return flow (undertow) and tidal current in the cross-shore direction; :math:`\bar{V}_c` is alongshore current which includes wave-induced alongshore current and tidal current in the alongshore direction. 
Subscriptions :math:`c` and :math:`l` represent cross-shore and longshore, respectively (hereafter).  

In Kobayashi et al., the vector (:math:`\bar{U}_c, \bar{V}_l`) is defined specifically in cross-shore/alongshore directions which may not be appropriate for general 2-D applications, especially for complex coastal geometry. 
Here, we define :math:`\vec{\bar{V}} = (\bar{U},\bar{V})` as the Eulerian-mean current  vector in general Cartesian coordinates  which may not be restricted to the cross-shore/alongshore orientation.
Thus formula :eq:`qs` can be written in a general form using coordinate rotation as

.. math:: q_{sx} =( a \bar{U} \cos^2 \alpha + a \bar{V} \sin \alpha \cos \alpha + \bar{U} \sin^2 \alpha -\bar{V}\sin \alpha \cos \alpha) V_s

.. math:: q_{sy} = (a \bar{V} \sin^2 \alpha + a \bar{U} \sin \alpha \cos \alpha + \bar{V} \cos^2 \alpha -\bar{U}\sin \alpha \cos \alpha) V_s

in which :math:`\alpha` is an angle between :math:`x-` axis and beach-normal direction as shown in :numref:`beach_angle`.

.. figure:: figures/beach_direction.jpg
   :name: beach_angle
   :scale: 40 %
   :align: center
   
   Definition of an angle between :math:`x-` axis and beach-normal direction.
   

The bedload transport rates can be expressed as

.. math:: q_{bc} = \frac{b P_b}{g (s-1)} \sigma^3_T (1+ U_* V_*^2 +2 F_m \sin \theta) G_s

.. math:: q_{bl} = \frac{b P_b}{g (s-1)} \sigma^3_T \left[ V_* (1+U^2_* + V^2_*) -2 r_m \sin \theta \right ]

where :math:`b` is empirical bedload parameter; :math:`\sigma_T` is standard deviation of the oscillatory depth-averaged velocity with zero mean; :math:`(U_*, V_*) = (\bar{U}_c/\sigma_T, \bar{V}_l/\sigma_T)`; :math:`\theta` is wave angle; and 

.. math:: r_m = -(U_*\cos \theta + V_* \sin \theta)

.. math:: F_m = V_* \cos \theta - U_* \sin \theta

:math:`G_s` is the bottom slope function expressed by

.. math::
   \begin{array}{l l l}
      G_s = \tan \phi /(\tan \phi + S_b)  &  \mbox{for}  &  -\tan \phi < S_b < 0 \\
      G_s = (\tan \phi -2 S_b) /(\tan \phi -S_b)  & \mbox{for} & 0 < S_b < \tan \phi
    \end{array} 

in which :math:`\phi` is the angle of internal friction of the sediment and :math:`\tan \phi \simeq 0.63` for sand.
  
  In a general Cartesian system, (:math:`q_{bc}, q_{bl}`) are converted to (:math:`q_{bx}, q_{bl}`) using

.. math:: q_{bx} = q_{bc} \cos \alpha - q_{bl} \sin \alpha,  \ \ \ \ \ \  q_{by} = q_{bc} \sin \alpha + q_{bl} \cos \alpha





Seabed evolution equation
##################################

The sea bed evolution equation can be described using sediment transport flux in generalized curvilinear coordinates:

.. math:: (1-p) \frac{\partial h}{\partial t} + \frac{1}{J} \frac{\partial J f_{mor} q_\beta L^\alpha_\beta}{\partial \xi^\alpha} = 0
   :\label: bed

where :math:`p` is the bed porosity, :math:`f_{mor}` represents a morphology factor. :eq:`bed` is solved using the same TVD scheme as in SHORECIRC.  






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






*************************************************
**SWAN-SHORECIRC coupling and parallelization**
*************************************************

For the SHORECIRC code, we used a domain decomposition technique to subdivide the problem into multiple regions and assign each subdomain to a separate processor core. 
Each subdomain region contains an overlapping area of ghost cells, three-row deep, as required by the highest-order ( the fourth order) MUSCL-TVD scheme in the model. 
The Message Passing Interface (MPI) with non-blocking communication is used to exchange  data in the overlapping region between neighboring processors. 

For the SWAN code, we used the existing domain decomposition scheme which is slightly different from the SHORECIRC code.
The version 40.51AB of SWAN uses a single-direction domain decomposition and takes into account number of dry points in the grid splitting. Because of the different domain decomposition between SHORECIRC and SWAN, it is difficult to directly pass variables between the decomposed domains. 
We now used a two-step mapping method which is  to first gather a passing variable into the global domain, and then to distribute it into each sub-domain. 
Dry points in SWAN are set to be wet with a small water depth (1 mm) in order to make an equal split of a computational domain. 

:numref:`domain` shows an example of domain decomposition in SWAN and SHORECIRC. 

To investigate performance of the parallel program, numerical simulations of an 
idealized  case are tested with different numbers of processors  on a Linux cluster located at University of Delaware. 
The test case is set up in a numerical grid  of 1800 :math:`\times` 1800 cells.  
:numref:`fig1` shows the model speedup versus number of processors. 
It can be seen that performance scales nearly proportional to the number of processors, with some delay caused by inefficiencies in parallelization, such as inter-processor communication time.


.. figure:: figures/parallel.jpg
   :name: domain
   :scale: 40%
   :align: center

   An example of domain decomposition in SWAN and SHORECIRC/SEDIMENT.


.. figure:: figures/speedup.jpg
   :name: fig1
   :scale: 35%
   :align: center

   Variation in model performance with number of processors for a 1800 x 1800 domain.  Straight line indicates arithmetic speedup. Actual performance is shown in the curved line **(the figure needs to be updated using the coupled model !!!)**.







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



