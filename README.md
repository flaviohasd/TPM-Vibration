# TPM-Vibration 
  Version: 2.1
  
  Release: 14/08/2022
  
  Last Update: 22/10/2022
   

## Introduction
  Hello, my name is Flávio Dias. I'm a Mechanical Engineer graduated at Federal University of Amazonas (Brazil) and this project was developed for my final paper. The goal is to minimize the amplitude of the vibration signal of a spur gear mesh system (with contact ratio between 1 and 2) by performing a TPM - Tooth Profile Modification (specifically, changing the maximum amount of tip relief -adendum modification-) for a certain working condition.
  
  Link to the paperwork: <a href= http://riu.ufam.edu.br/handle/prefix/6549 > DIAS, Flávio. Minimizing gear mesh vibration by the tooth profile modification of
cylindrical spur gears. 2022. 113f. Undergraduate thesis (Mechanical Engineering) - Federal University of Amazonas, Manaus-AM, 2022.

## Abstract
  Gears are among the oldest devices invented by man and, among the various forms of power transmission, are generally the most robust and durable, with efficiency in the order of 98%. The dynamic characteristic of a gear transmission system is an important parameter to define its quality. The variation in the number of teeth in contact during rotation of the gears has a large contribution to the vibrations of these systems. From this context derives the need to evaluate the dynamic characteristics of the gear while still in the design phase. This was a study conducted on the addendum modification in the tooth profile of gears. For this purpose, dynamic model simulations were performed in order to determine an optimal value of the modification parameter capable of minimizing the vibrations of an ideal cylindrical spur gear for the specified operating conditions, considering time-varying gear stiffness and damping, including an improvement in the equivalent mass calculation. From the comparison between the simulation results of the optimized and original system, it is concluded that the optimal modification brings the dynamic responses closer to the static ones, with change rates close to zero and nearly constant transmission error, providing minimal vibration to the system.
  
## Code Description  
  
  Please note: the code is writen in MATLAB in PT-BR language (might make an EN version someday) and is structured as follows:
  
  <p align="center">
    <img src="https://i.imgur.com/x1Hqzyq.jpg" alt="Optimization Flowchart" width="400" />
  </p>
  
  OTIMIZACAO: Has the properties of the optimization function with the goal to evaluate the main code (SIMULACAO) to get the optimum amount of modification;
  
  SIMULACAO: The main code, where all the constants and properties of the gear, mesh and operation. Holding all the other functions mentioned below;
  
  ENGRENAGEM: Calculates the gear dimentions given the entries in SIMULACAO;
  
  OPERACAO: Calculates the gear mesh relations given the entries in SIMULACAO;
  
  ENGRENAMENTO: Calculates the mesh parameters given the gear and operation entries in SIMULACAO;
  
  MASSA: Calculates de equivalent mass of the sistem given the dimentions and the profile modification amount;
  
  INTEG: Calculates de integrals for the stiffness by the energy formulas;
  
  ENERGIA: Calculates de stiffness by the energy formulas of the original and modified systems;
  
  INTERP: Interpolate the discrete functions to obtain the continuos representation of the dynamic coefficients;
  
  MODELO: Sets the space-states variables together with the acceleration of the original and modified system;
  
  DINAMICA: Solves the original dynamic model;
  
  DINAMICATPM: Solves the modified dynamic model;
  
  ACCELERACAO: Calculates the acceleration of the original and modified systems;
  
  SINAL: Quantifies the vibration signals;
  
  PLOTS: Generate the graphics for analysis.
  
  
  OUTPUT EXAMPLES:
  
  <p align="center">
  <img src="https://i.imgur.com/24b0nM1.jpg" alt="Vibration over a tooth cycle" width="800" />
  </p>

  <p align="center">
  <img src="https://i.imgur.com/Qlf5hNn.jpg" alt="Contact relations" width="800" />
  </p>
  
  
  For more technical details, take a look in the paperwork (for now, only available in portuguese).
  
  If you have any questions regarding the theory, models, code, or anything related, fell free to reach me out.
  
  Cheers!
