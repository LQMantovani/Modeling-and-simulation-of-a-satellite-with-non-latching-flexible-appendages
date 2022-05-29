This simulation environment was developed as a part of the MSc research published with title "Modeling and simulation of a satellite with non-latching flexible appendages", from the Aeronautics Institute of Technology (ITA - www.ita.br).
It aims to simulate space systems with flexible non-latching appendages, allowing an arbitrary number of rigid and flexible bodies places anywhere in the main body. It employs a multibody formulation based on
the floating frame of reference approach, which allows large rotations of the bodies, but just small deformations. 
The simulation environment uses modal shapes to described the flexible bodies. It also allows concentrated masses and inertias at the tip of the flexible beam (which are considered on both the modal shape and 
the inertia of the body), as well as concentrated masses and inertias along the boom (in this case the modal shape is not obtained considering the concentrated mass and inertia).

Two constraint options are available:
 - Fully contrained: Constrains three translational and three rotational degrees of freedom.
 - Revolute joint constraint: Constraints three translational degrees of freedom and two rotational degrees of freedom.

The code also has implemented constraint controls to avoid numerical instabilities associated with the constraints (due to numerical errors during the integration process). The gain of the constraints control can be selected (or set to zero to disconsider them).

The code also allows for the input of the non-latching mechanism parameter (spring deployment parameters and collision parameters). For more information, see references [1] and [2]. 
Furthermore, it is possible to introduce a control system in the simulation environment. The function called control has as an output the control force and torque to be applied at the first body's reference system.

The results will be displayed in graphics for each body. Next, a post-processing is applied to obtain the relative position and rotation of the bodies regarding the first body.

All data input is through a .txt file, excpet for the simulation time that is asked in the command window. The implemented codes allow for the generic slender flexible bodies (subjected to small deformations), 
obtaining the shape matrices each time a new body is introduced or any previous body's data is updated. The implemented dynamical models do not account for the dynamical stiffening, so flexible bodies with high angular velocities
(near their natural frequencies) will lead to large deformations and wrong results.

The gravitational torque is considering acting on the bodies, as well as the gravity of an ideal sphere (does not account for the Earth's oblateness and mass distribution). No drag nor solar radiation pressure is considered.

The simulation environment was developed with space systems in mind. Hence, the initial conditions inputs will be the body altitute, declination, and elevation. The simulation environment will then obtain the velocity of the body 
so it will be in a circular orbit. The initial conditions of the other bodies attached to the first body are obtained in a way to not violate the constraints. 

The simulation environment also can be used togheter with the Simulink World 3D Toolbox to visualize the behavior of the bodies. In this case, the model in Simulink needs to be prepare for each system, since the implemented code does not import the data automatically from the simulation environment. 
The function plot_relative_data uses parfor to obtain the relative angular positons of the bodies, which is a parallel computing method. It also uses the gpuArray function that allows to transfer vectors to the GPU for faster computation (is some cases). It is more of a example of possible approaches to 
improve the code efficiency, but it can be altered to a tradicional processing method (without a GPU) manually.
The code to obtain the FFT of the appendages relative angular position and velocity is also implemented. 


References:
[1] - Mantovani, L.Q. Modeling and simulation of a satellite with non-latching flexible appendages, Master thesis, Aeronautics Institute of Technology, 2022.
[2] - Mantovani, L.Q.; Santos, W.G.; Cardoso-Ribeiro, F.L.; Costa, L.E. V. L. Impact of non-latching booms on a nanosatellite attitude control system, Aerosp. Sci. Technol. 120 (2022) 107298. https://doi.org/10.1016/j.ast.2021.107298.