# TPM-Vibration 
  Version: 2.2
  
  Release: 14/08/2022
  
  Last Update: 09/04/2023
   

## Introduction
  Hello, my name is Flávio Dias. I'm a Mechanical Engineer graduated at Federal University of Amazonas (Brazil) and this project was developed for my final paper. The goal is to minimize the amplitude of the vibration signal of a spur gear mesh system (with contact ratio between 1 and 2) by performing a TPM - Tooth Profile Modification (specifically, changing the maximum amount of tip relief -adendum modification-) for a certain working condition.
  
  Link to the paperwork: <a href= http://riu.ufam.edu.br/handle/prefix/6549 > DIAS, Flávio. Minimizing gear mesh vibration by the tooth profile modification of
cylindrical spur gears. 2022. 113f. Bachelor's thesis (Mechanical Engineering) - Federal University of Amazonas, Manaus-AM, 2022.

## Abstract
  Gears are among the oldest devices invented by man and, among the various forms of power transmission, are generally the most robust and durable, with efficiency in the order of 98%. The dynamic characteristic of a gear transmission system is an important parameter to define its quality. The variation in the number of teeth in contact during rotation of the gears has a large contribution to the vibrations of these systems. From this context derives the need to evaluate the dynamic characteristics of the gear while still in the design phase. This was a study conducted on the addendum modification in the tooth profile of gears. For this purpose, dynamic model simulations were performed in order to determine an optimal value of the modification parameter capable of minimizing the vibrations of an ideal cylindrical spur gear for the specified operating conditions, considering time-varying gear stiffness and damping, including an improvement in the equivalent mass calculation. From the comparison between the simulation results of the optimized and original system, it is concluded that the optimal modification brings the dynamic responses closer to the static ones, with change rates close to zero and nearly constant transmission error, providing minimal vibration to the system.
  
## Code Description  
  
  The algorithm is structured as follows:

   <p align="center">
    <img src="https://user-images.githubusercontent.com/44821460/229648230-598b24b0-9d1d-4a95-9a83-59e314eb715f.png" alt="Optimization Flowchart" width="400" />
  </p>
  
  OPTIMIZATION: Has the properties of the optimization function with the goal to evaluate the main code (SIMULATION) to get the optimum amount of modification;
  
  SIMULATION: The main code, where all the constants and properties of the gear, mesh and operation. Holding all the other functions mentioned below;
  
  GEARS: Calculates the gear dimentions given the entries in SIMULATION;
  
  OPERATION: Calculates the gear mesh relations given the entries in SIMULATION;
  
  MESH: Calculates the mesh parameters given the gear and operation entries in SIMULATION;
  
  MASS: Calculates de equivalent mass of the sistem given the dimentions and the profile modification amount;
  
  INTEGR: Calculates de integrals for the stiffness by the energy formulas;
  
  ENERGY: Calculates de stiffness by the energy formulas of the original and modified systems;
  
  INTERPOLATION: Interpolate the discrete functions to obtain the continuos representation of the dynamic coefficients;
  
  MODEL: Sets the space-states variables together with the acceleration of the original and modified system;
  
  DYNAMIC: Solves the original dynamic model;
  
  DYNAMICTPM: Solves the modified dynamic model;
  
  ACCEL: Calculates the acceleration of the original system;
  
  ACCELTPM: Calculates the acceleration of the modified system;
  
  SIGNAL: Quantifies the vibration signals;
  
  PLOTS: Generate the graphics for analysis.
  
  
  OUTPUT EXAMPLES:
  
  <p align="center">
  <img src="https://user-images.githubusercontent.com/44821460/230794069-064d7b4e-b72a-4928-a2f7-47bb70ca2b29.png" alt="Vibration over a tooth cycle" width="800" />
  </p>

  <p align="center">
  <img src="https://user-images.githubusercontent.com/44821460/230794334-d4810585-1a9b-4f79-8c62-2e1ce66c7869.png" alt="Contact relations" width="800" />
  </p>
  
  
  For more technical details, take a look in the paperwork (for now, only available in portuguese - EN article incoming by the end of the year).
  
  If you have any questions or ideas regarding the theory, models, code, or anything related, feel free to reach me out.
  
  Cheers!
