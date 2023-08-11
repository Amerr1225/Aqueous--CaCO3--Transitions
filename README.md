# Aqueous--CaCO3--Transitions  
This model simulates the transition of calcite to aragonite in water.
It takes in parameters including:
- Masses for calcite and aragonite
- Temperature
-pH 
- Amount of water  
- Timespan  
- Number of steps  
Etc...

The model uses the rate equation for calcite and aragonite, R = k * A * (1- SI)^n
where A is the surface area of the solid SI is the saturation n is the order of the reaction and k is the eq constant.

The code uses phreeqpy which is a python library which utalizes phreeqc software. The same code written as a phreeqc input file is also provided.

The model allows one to customise their own solutions composition, pH, density, etc or one can use preset Seawater and Pure water which are defined in the phreeqc guide from USGS.

The data is saved in a file called data1.dat which one can use to replot the data if wanted. The program also plots the data itself.

The Phreeqc guide goes into the details of how phreeqc works and the inner workings of both the phreeqc input file and the phreeqpy codes.

The sample data folder has data for both the salt and pure presets and then fresh water with an alkalinity set to 5. All of these runs had 0.05g of Aragonite and Calcite reacting and the solutions are equlibrated with CO2 before the kinetics is run. The reaction runs for 150 hours and has one million steps.

Hopefully this code provides some value either as a jumping of point for more refined modeling or as a tool to learn and experiment with phreeqc software and phreeqpy implementation.

To run this code through python one requires the following python extensions,
- numpy       -- Installed by default in jupyter enviornments only
- matplotlib  -- Installed by default in jupyter enviornemets only
- tqdm        -- Requries seperate installation
- phreeqpy    -- Requires seprate installation
- Ipython     -- Installed by default in jupyter enviornments
- io          -- Installed by default in jupyter and Pyhton default libraries
- sys         -- Installed by default in jupyter and Pyhton default libraries
- os          -- Installed by default in jupyter and Pyhton default libraries
- timeit      -- Installed by default in jupyter and Pyhton default libraries

Most of which are likely installed by default in jupyter however these are needed for the entire code to run without problems.

