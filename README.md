# TPM-Vibration 
  Version: 2.1
  
  Release: 14/08/2022
  
  Last Update: 22/10/2022
   

## Introduction
  Hello, my name is Flávio Dias. I'm a Mechanical Engineer graduated at Federal University of Amazonas (Brazil) and this project was developed for my final paper. The goal is to minimize the amplitude of the vibration signal of a spur gear mesh system (with contact ratio between 1 and 2) by performing a TPM - Tooth Profile Modification (specifically, changing the maximum amount of tip relief -adendum modification-) for a certain working condition.
  
  Link to the paperwork: http://riu.ufam.edu.br/handle/prefix/6549

## Description
  The idea of this post is to contribute with the learners, so they don't waste time doing what others already did and then develop even further. You can reproduce this as long you keep it well referenced.
  
  Please note: the code is writen in MATLAB in PT-BR language (might make an EN version someday) and is structured as follows:
  
  <p align="center">
  <img src="https://i.imgur.com/x1Hqzyq.jpg" alt="Optimization Flowchart" width="400" />
  </p>
  
  OTIMIZACAO: Has the properties of the optimization function with the goal to evaluate the main code to get the optimum value of modification;
  
  SIMULACAO: The main code, where all the constants and properties of the gear mesh and operation. Holding all the other functions mentioned below;
  
  ENGRENAGEM: It calculates the dimentions given the entries in SIMULACAO;
  
  OPERACAO: Calculates the gear mesh relations given the entries in SIMULACAO;
  
  ENGRENAMENTO: Calculates the mesh parameters given the gear and operation entries in SIMULACAO;
  
  MASSA: Calculates de equivalent mass of the sistem given the dimentions and the profile modification;
  
  INTEG: Calculates de integrals for the stiffness by the energy formulas;
  
  ENERGIA: Calculates de stiffness by the energy formula of the original and modified systems;
  
  INTERP: Interpolate the discrete functions to obtain the continuos representation of the dynamic coefficients;
  
  MODELO: Sets the space states variables together with the acceleration of the original and modified system;
  
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
  
  
  For more technical details, take a look in the paperwork (only available in portuguese).
  
  Cheers.
