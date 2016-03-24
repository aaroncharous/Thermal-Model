Thermal Model for CubeSat
Version 1: Archaeon
Written Finn Van Krieken -- edited by Aaron Charous, Eddie Williams, and Cameron Reid

This is a nodal thermal model for our 1U CubeSat. It uses an odeint (Ordinary Differential Equation INTegrating function)
to obtain a numerical solution for the temperature of the various nodes as the satellite orbits. Additional nodes can be
added in the function dy_dt(y,t). Just make sure you add them in the unpacking of the vector and the packing of the solution,
as well as defining the actual derivative of the added node (see the way each other node is handled). Make sure that if you
add a connection from one node to another, that the connection is in each differential equation, or suddenly energy is going
to be vanishing from the satellite or being generated out of nowhere.

As for the layout, you can see that the first section is thermal constants, segmented by what part of the system they
pertain to. Then are the various functions (orbit, latitude, orientation, etc.) and the differential equations/solutions,
and then the plotting/graphing section.

If you have questions let me know (finnvankrieken@gmail.com)

This uses numpy, which is pretty well documented if you need other mathematical functions (just
remember you have to enter them as np.func(x) rather than just func(x))

We chose not to model the controller board because it emitted almost no energy.

Additionally, when adding a new component, make sure you add the thermal resistance connections.

rtol is the “tolerance” of odeint, meaning the error bound of the numerical approximation.  Going above 10^-5 is dangerous, as it sometimes misses flashes.  Going below 10^-5 isn’t dangerous, but you may or may not be alive for the program to finish running.  The default rtol for odeint is 10^-8.

Quick Guide to adding a new module:
1. Enter necessary constants including emissivity, area, mass, specific heat, thermal conductivity, and thickness.
2. Define thermal resistance by dividing the thickness of the component by the product of the surface area in contact with other components and the conductivity.
3.  Define Thermal Resistances of connections.  Like a circuit, thermal resistances add in series.
4.  Pack a new module into y
5.  Define the differential equation for the temperature of the component.  Usually the differential equation includes radiation outwards, solar radiation, Earth’s albedo, and the connections between modules.  
6.  Pack the differential equation into dydt.
7.  In run_simulation, make sure you add another initial condition for the module, which you can do by adding 1 to y0s.
8.  Define the module as the answer[:, #], where # is the number of the module.
9.  Create a new figure and axes to plot the temperature of the specific module.

Happy modeling!
