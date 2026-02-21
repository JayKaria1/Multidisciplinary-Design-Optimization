# MDO Wing Optimization Project

This repository contains the code, scripts, and analysis from a Multi-Disciplinary Design Optimization (MDO) project I worked on together with a project partner. The overall idea was to build a reasonably complete MDO framework and then use it to redesign a transport aircraft wing.

The case study we focused on was the McDonnell Douglas DC-10-30. The objective was to reduce the aircraft’s maximum take-off weight (MTOW) for a fixed mission range, without increasing its original CO₂ emissions. Engine parameters were kept fixed, so any efficiency gains had to come from aerodynamics and structures rather than propulsion tweaks.

## General Setup

The framework couples four different disciplines:
- Aerodynamics  
- Loads  
- Structures  
- Performance  

These were integrated using a Multi-Discipline Feasible (MDF) architecture. Each discipline is solved sequentially inside the optimization loop until consistency is reached, which makes the setup slower but easier to manage and debug.

For aerodynamics, we used Q3D Solver with viscous corrections. Structural sizing was handled using EMWET, which provided wing weight estimates based on loads and geometry. Performance calculations relied on standard fuel fraction methods and Breguet’s range equation. Since the turbine inlet conditions were fixed and limiting, engine performance was not part of the design space.

An XDSM of this optimisation can be found [here](MDF_XDSM_Final.drawio.pdf)

## Design Variables

The optimizer worked with a total of 20 design variables. These covered a mix of geometric, aerodynamic, and operational parameters, including:

- Wing planform variables such as span, wing area, taper ratio, and leading-edge sweep
- Airfoil shape changes applied to a supercritical Whitcomb 135 airfoil, using CST parameters on both upper and lower surfaces
- Cruise conditions, including Mach number and cruise altitude
- Internal wing layout variables, specifically the front and rear spar locations, which directly affected the available fuel tank volume

This made the design space fairly high-dimensional and strongly coupled, especially between aerodynamics, loads, and structures.

## Constraints

A few hard constraints were enforced to keep the design realistic:

- Wing loading could not exceed the baseline aircraft value  
- Total CO₂ emissions had to be less than or equal to the original aircraft  
- The internal volume between the spars had to be sufficient to physically store all mission fuel  

Without these constraints, the optimizer would quickly drift toward impractical or non-certifiable designs.

## Optimization and Results

The optimization was performed in MATLAB using `fmincon` with a Sequential Quadratic Programming (SQP) algorithm. Despite the complexity of the problem, the optimizer converged to a well-behaved and efficient solution.

Key results:
- MTOW was reduced by 11.5%, from 263,636 kg to 233,260 kg  
- Wing weight was reduced by 59.3%, driven mainly by a reduction in leading-edge sweep and overall wing area  
- Fuel burn and CO₂ emissions were reduced by 14.5%  
- The optimized aircraft cruises slightly faster and at a marginally lower altitude than the baseline
- [Isometric comparison of original and optimised wing](wing_iso.pdf)



Overall, the project showed that meaningful reductions in MTOW and emissions are possible through integrated aerodynamic and structural optimization, even when propulsion parameters are frozen. The framework itself is modular and can be extended to other aircraft configurations or additional disciplines if needed.

