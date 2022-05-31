function U = control(t, X, Bodies)
%Put yout control function here
%Input:
%t: time
%X: state vector
%Bodies: Bodies data structure
%Output:
%U: Data structure with two fields: U.toruqe is the torque applied at the
%first body written in the first body reference frame. U.force is the force 
%applied to the first body and is written in the first body reference frame 


U.force = [0;0;0];
U.torque = [0;0;0];