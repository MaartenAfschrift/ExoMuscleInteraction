# Simulate influence exoskeleton on muscle dynamics



The goal of this project is to simulate and optimize the influence of exoskeleton assistive torques on muscle dynamics during activities of daily living.



## Dependencies

Several software packages are needed to run the program

- The OpenSim MATLAB interface is used to generate the inputs to the optimal control problem based on a scaled OpenSim model and the solution of inverse kinematics (providing the solution of inverse dynamics is optional). To this aim, install OpenSim and set up the OpenSim MATLAB interface (OpenSim: [https://simtk.org/frs/?group_id=91](https://simtk.org/frs/?group_id=91), OpenSim API: http://simtk-confluence.stanford.edu:8080/display/OpenSim/Scripting+with+Matlab.
- GPOPS-II is used to solve the optimal control problem using direct collocation (\url{http://www.gpops2.com/}). A one-time 30-day trial license is avaiable for all users who register. (Note, we are currently in the transition to Casadi.  Casadi is an open-source tool for nonlinear optimization and algorithmic differentiation (<https://web.casadi.org/>))
- ADiGator is used for automatic differentiation https://sourceforge.net/projects/adigator/.(not needed in Casadi version)
- Matlab code to solve muscle redundancy (https://github.com/antoinefalisse/solvemuscleredundancy_dev)
- Matlab code to compute metabolic energy consumption (https://github.com/MaartenAfschrift/MetabolicEnergy_Simulation)



## Example

 ### Optimize AFO

In this example, the actuation profile of an ankle foot exoskeleton is optimized to minimize metabolic energy consumption during walking. Actuation constratins are implemented based on Zhang2017 for comparision with the optimal profile determined using human-in-the-loop optimization (Zhang2017: http://dx.doi.org/10.1126/science.aal5054). The influence of achilles tendon stiffnes is evaluated to document the importance of muscle-tendon interaction during walking.

OpenSim gait2392 model is used to simulate muscle and skeleton dynamics. Joint kinematics and kinetics are computed based on normal walking on the treadmill without exoskeleton. Note that in inverse approach with invariant kinematics is used (i.e. joint kinematics and kinetics are invariant/ not optimized). In addition, muscle activations and exoskeleton assistance is optimized for one leg accounting for the ankle, knee and hip degrees of freedom.

Script to run example:



Script to plot results example:



Results:

