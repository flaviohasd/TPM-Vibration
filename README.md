# TPM-Vibration 
  Version: 2.0
  
  Release: 26/08/2022
   
  Note: This project is still under development! Changes may occour often.

## Introduction
  Hello, my name is Flávio Dias. I'm a Mechanical Engineering undergraduate student at Federal University of Amazonas (Brazil) and this project was developed for my final paper. The aim of this project is to minimize the amplitude of the vibration signal of a spur gear mesh system by changing the maximum amount of tip relief (adendum modification) for a certain condition of operation. The paper is currently beeing writen at the moment and the link to it will be available when the official version becomes public.

## Description
  The idea of this post is to contribute with the learners, so they don't waste time doing what others already did and then develop even further. You can reproduce this as long you keep it well referenced.
  
  Please note: the code is writen in MATLAB in PT-BR language (might upload an EN version soon) and is structured as follows:
  
  <p align="center">
  <img src="https://i.imgur.com/wAHjoAF.jpg" alt="Optimization Flowchart" width="400" />
  </p>
  
  OPTMIZATION: Has the properties of the optimization function with the goal to evaluate the main code to get the optimum value of modification;
  
  SIMULATION: The main code, where all the constants and properties of the gear mesh and operation. Holding all the other functions mentioned below;
  
  GEAR: It calculates the dimentions given the entries in SIMULATION;
  
  OPERATION: Calculates the gear mesh relations given the entries in SIMULATION;
  
  MESH: Calculates the mesh parameters given the gear and operation entries in SIMULATION;
  
  MASS: Calculates de equivalent mass of the sistem given the dimentions and the profile modification;
  
  INTEG: Calculates de integrals for the stiffness by the energy formulas;
  
  ENERGY: Calculates de stiffness by the energy formula of the original and modified systems;
  
  INTERP: Interpolate the discrete functions to obtain the continuos representation of the dynamic coefficients;
  
  MODEL: Sets the space states variables together with the acceleration of the original and modified system;
  
  DYNAMIC: Solves the original dynamic model;
  
  DYNAMICTPM: Solves the modified dynamic model;
  
  ACCELERATION: Calculates the acceleration of the original and modified systems;
  
  SIGNAL: Quantifies the vibration signals;
  
  PLOTS: Generate the graphics for analysis.
  

  
  Cheers.
