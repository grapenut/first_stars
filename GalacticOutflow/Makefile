
Simulation += Simulation_data.o sim_checkStatus.o sim_initNova.o sim_interface.o sim_setupNova.o sim_markInitial.o sim_markOverdensity.o sim_markZoom.o

Simulation_init.o : Simulation_data.o
Simulation_initBlock.o : Simulation_data.o
Particles_initPositions.o : Simulation_data.o Particles_data.o 
sim_checkStatus.o : sim_interface.o
sim_initNova.o : sim_interface.o
sim_markInitial.o : sim_interface.o
sim_markOverdensity.o : sim_interface.o
sim_markZoom.o : sim_interface.o
sim_setupNova.o : sim_interface.o
sim_interface.o : Simulation_data.o
