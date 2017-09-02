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

   
