"""
The following python code is based of of a PHREEQC input file which I wrote. Using PHREEQPY it runs and plots in a python environment.
The Model below simulates the transition between Aragonite and Calcite in water.
"""

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from IPython.display import clear_output
import time as tim
import io
import sys
import os
import timeit

# Simple Python 3 compatibility adjustment.
if sys.version_info[0] == 2:
    range = xrange

MODE = 'dll'  # 'dll' or 'com'

if MODE == 'com':
    import phreeqpy.iphreeqc.phreeqc_com as phreeqc_mod
elif MODE == 'dll':
    import phreeqpy.iphreeqc.phreeqc_dll as phreeqc_mod
else:
    raise Exception('Mode "%s" is not defined use "com" or "dll".' % MODE)

    
def run_phreeqc(input_string):
    """
    Run PHREEQC simulation with the given input string.
    """
    phreeqc = phreeqc_mod.IPhreeqc()
    phreeqc.load_database(r"C:\Program Files (x86)\USGS\Phreeqc Interactive 3.7.3-15968\database\wateq4f.DAT")
    phreeqc.run_string(input_string)
    return phreeqc


def extract_data(phreeqc, junkrow=1, i = 0):
    """
    Extract data from PHREEQC simulation results and save it to 'data.dat'.
    """
    # Get all selected output data as separate lists
    time = phreeqc.get_selected_output_column(0)
    pH = phreeqc.get_selected_output_column(1)
    SI_ara = phreeqc.get_selected_output_column(2)
    SI_cal = phreeqc.get_selected_output_column(3)
    Ara = phreeqc.get_selected_output_column(4)
    dAra = phreeqc.get_selected_output_column(5)
    Cal = phreeqc.get_selected_output_column(6)
    dCal = phreeqc.get_selected_output_column(7)

    while i < junkrow:
        time.pop(0)
        pH.pop(0)
        SI_ara.pop(0)
        SI_cal.pop(0)
        Ara.pop(0)
        dAra.pop(0)
        Cal.pop(0)
        dCal.pop(0)
        i += 1
        
    # Set Column headers.
    header = "Time\tpH\tSI_ara\tSI_cal\tAra\tdAra\tCal\tdCal\n"

    # Save data to 'data.dat'
    with open("data1.dat", "w") as data_file:
        data_file.write(header)
        for i in range(len(time)):
            data_file.write(f"{float(time[i]):.16f}\t{float(pH[i]):.16f}\t{float(SI_ara[i]):.16f}\t{float(SI_cal[i]):.16f}\t"
                            f"{float(Ara[i]):.16f}\t{float(dAra[i]):.16f}\t{float(Cal[i]):.16f}\t{float(dCal[i]):.16f}\n")

    return time, pH, SI_ara, SI_cal, Ara, dAra, Cal, dCal

def make_plots(time, data1, data2, label1, label2, x_label, y_label, Title, legend_loc='lower right', fig_size=(10, 4.5), semi=False, 
               compare=True, Limits = False, ylU=0, ylL=0, xlU=0, xlL=0):
    """
    Create plots for the given data.

    Parameters:
        time (list): Time data for the x-axis.
        data1 (list): First dataset to plot.
        data2 (list): Second dataset to plot.
        label1 (str): Label for the first dataset.
        label2 (str): Label for the second dataset.
        x_label (str): Label for the x-axis.
        y_label (str): Label for the y-axis.
        legend_loc (str, optional): Location of the legend. Default is 'lower right'.
        fig_size (tuple, optional): Size of the figure. Default is (10, 4.5).
        ylimit (Optional): if true will set limits.
        ylU and ylL (Optional): Set the upper and lower limits of the plots y axis.
        xlimit (Optional): if true will set limits.
        xlU and xlL (Optional): Set the upper and lower limits of the plots x axis.
    """

    fig, ax = plt.subplots(figsize=fig_size)
    ax.set_facecolor('grey')
    ax.set_xlabel(x_label, color='black')
    ax.set_ylabel(y_label, color='black')
    ax.tick_params(axis='x', colors='black')
    ax.tick_params(axis='y', colors='black')
    if compare:
        ax.set_title(Title, color='black')
        if semi:
            ax.semilogy(time, data1, color='Maroon', label=label1)
            ax.semilogy(time, data2, color='Pink', label=label2)
        else:
            ax.plot(time, data1, color='Maroon', label=label1)
            ax.plot(time, data2, color='Pink', label=label2)

        legend = ax.legend(loc=legend_loc, facecolor='grey')
        for text in legend.get_texts():
            text.set_color('black')
    else:
        ax.set_title(Title, color='black')
        if semi:
            ax.semilogy(time, data1, color='Maroon', label=label1)
        else:
            ax.plot(time, data1, color='Maroon', label=label1)

        legend = ax.legend(loc=legend_loc, facecolor='grey')
        legend.get_texts()[0].set_color('black')
    if Limits[0] == True:
        ax.set_xlim(xlL,xlU)
    if Limits[1] == True:
        ax.set_ylim(ylL,ylU)

    return fig

def plot_data(time, SI_ara, SI_cal, Ara, dAra, Cal, dCal, pH, t_unit, xLimits, yLimits, Limits, Titles, Labels, TTime, TSteps):
    """
    Plot the data extracted from PHREEQC simulation.

    Parameters:
        time (list): Time data.
        SI_ara (list): Saturation Index data for aragonite.
        SI_cal (list): Saturation Index data for calcite.
        Ara (list): Aragonite data.
        dAra (list): Aragonite dissolution rate data.
        Cal (list): Calcite data.
        dCal (list): Calcite dissolution rate data.
        t_unit (str): Time unit for the x-axis.
        TTime = total time
        TSteps = Total steps
    """
    Null = [0] * len(time)
    totmol = Ara[0] + Cal[0]
    stepsize = TTime / TSteps

    χ_ara = [(a / totmol) for a in Ara]
    # χ_Cal = [(c / totmol) for c in Cal]
    χ_cal = [1 - χ for χ in χ_ara]

    xlims = xLimits
    ylims = yLimits
    Limits = Limits
    Tits = Titles
    Labs = Labels

    # Initialize the progress bar
    progress_bar = tqdm(total=len(dAra))

    for i in range(len(dAra)):
        time[i] = time[i] / 3600
        dAra[i] = dAra[i] / stepsize
        dCal[i] = dCal[i] / stepsize

        # Uncomment the following if you want to have a progress bar, will impact performance!
        # Update the progress bar
        progress_bar.set_description(f"Processing data point: {i + 1}/{len(dAra)}")
        progress_bar.update(1)

    # Close the progress bar
    progress_bar.close()

    # Create figures for each plot
    fig1 = make_plots(time, Ara, Cal, "Aragonite", "Calcite", "Time (" + t_unit + ")", "Amount of Solid (Mol)", "Aragonite and Calcite Levels Over Time", semi=True,
                       Limits = Limits, ylU=ylims[1], ylL=ylims[0], xlU=xlims[1], xlL=xlims[0])
    fig2 = make_plots(time, Ara, Cal, "Ara", "Calcite", "Time (" + t_unit + ")", "Amount of solid Aragonite", "Amount of Aragonite and Calcite", compare=True,
                       Limits = Limits, ylU=ylims[1], ylL=ylims[0], xlU=xlims[1], xlL=xlims[0])    
    fig3 = make_plots(time, dAra, dCal, "Aragonite", "Calcite", "Time (" + t_unit + ")", "Rate Dissolved (mol/h)", "Rate of Change Over Time", legend_loc='upper right',
                       Limits = Limits, ylU=ylims[1], ylL=ylims[0], xlU=xlims[1], xlL=xlims[0])
    fig4 = make_plots(time, χ_ara, χ_cal, "Aragonite", "Calcite", "Time (" + t_unit + ")", "Mole Fraction", "Relative Levels Over Time", legend_loc='upper right',
                       Limits = Limits, ylU=ylims[1], ylL=ylims[0], xlU=xlims[1], xlL=xlims[0])
    fig5 = make_plots(time, SI_ara, SI_cal, "Aragonite", "Calcite", "Time(" + t_unit + ")", "Saturation Index", "Saturations Over Time",
                       Limits = Limits, ylU=ylims[1], ylL=ylims[0], xlU=xlims[1], xlL=xlims[0])
    fig6 = make_plots(time, pH, Null, "pH", "Null", "Time (" + t_unit + ")", "pH", "pH of Solution", compare=False,
                       Limits = Limits, ylU=ylims[1], ylL=ylims[0], xlU=xlims[1], xlL=xlims[0])

    # Save figures as images in memory
    buffer1 = io.BytesIO()
    buffer2 = io.BytesIO()
    buffer3 = io.BytesIO()
    buffer4 = io.BytesIO()
    buffer5 = io.BytesIO()
    buffer6 = io.BytesIO()

    fig1.savefig(buffer1, format='png')
    fig2.savefig(buffer2, format='png')
    fig3.savefig(buffer3, format='png')
    fig4.savefig(buffer4, format='png')
    fig5.savefig(buffer5, format='png')
    fig6.savefig(buffer6, format='png') 

    # Rewind the buffer and read the image data
    buffer1.seek(0)
    buffer2.seek(0)
    buffer3.seek(0)
    buffer4.seek(0)
    buffer5.seek(0)
    buffer6.seek(0)

    image1 = plt.imread(buffer1)
    image2 = plt.imread(buffer2)
    image3 = plt.imread(buffer3)
    image4 = plt.imread(buffer4)
    image5 = plt.imread(buffer5)
    image6 = plt.imread(buffer6)

    # Create a 3x2 grid to display the plots side by side
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(15, 15))
    fig.suptitle("Plotted Data", fontsize=16)

    # Adjust individual subplot sizes to make the images larger
    axes[0, 0].imshow(image1)
    axes[0, 0].axis('off')
    axes[0, 0].set_aspect('auto')  # Set the aspect ratio to auto to fill the available space

    axes[0, 1].imshow(image2)
    axes[0, 1].axis('off')
    axes[0, 1].set_aspect('auto')  # Set the aspect ratio to auto to fill the available space

    axes[1, 0].imshow(image3)
    axes[1, 0].axis('off')
    axes[1, 0].set_aspect('auto')  # Set the aspect ratio to auto to fill the available space

    axes[1, 1].imshow(image4)
    axes[1, 1].axis('off')
    axes[1, 1].set_aspect('auto')  # Set the aspect ratio to auto to fill the available space

    axes[2, 0].imshow(image5)
    axes[2, 0].axis('off')
    axes[2, 0].set_aspect('auto')  # Set the aspect ratio to auto to fill the available space

    axes[2, 1].imshow(image6)
    axes[2, 1].axis('off')
    axes[2, 1].set_aspect('auto')  # Set the aspect ratio to auto to fill the available space

    # Adjust spacing
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    
    # Close the solo plots
    plt.close(fig1)
    plt.close(fig2)
    plt.close(fig3)
    plt.close(fig4)
    plt.close(fig5)
    plt.close(fig6)

    # Show the grid plot
    plt.show()

def Generate_Solutions(Name, temp = "20", pH = "7", units = "mmol/L", density = "1",alkalinity = 0, mass = "1", Components = [] ,preset = "Null"):
    """
    This function generates a solution which will be saved and returned to main.
    Name == Name of soltuion
    Temp == temperature of solution
    pH == pH of solution
    units == Units of the concentration
    Density == Density of solvent
    mass= = mass of solvent in kg
    """

    if preset.lower() == 'salt':
        Solution = (
            "SOLUTION " + Name + "\n" 
            + "    temp        "+ temp + "\n" 
            + "    pH          8.22 \n"
            + "    density     1.023 \n"
            + "    units       ppm \n"
            + "    Ca          412.3 \n"
            + "    Mg          1291.8 \n"
            + "    K           399.1 \n"
            + "    Si          4.28 \n"
            + "    Cl          19353.0 \n"
            + "    Alkalitnity 141.682 as HCO3 \n"
            + "    S(6)        2712.0 \n"
            + "    -water      " + mass  
            + "\n\n" )
        
    elif preset.lower() == 'pure':
        Solution = (
            "SOLUTION " + Name + "\n" 
            + "    temp       "+ temp + "\n" 
            + "    pH         " + "7.0" + "\n"  
            + "    -water     " + mass
            + "\n\n" )

    else:   
        Solution = (
            "SOLUTION " + Name + "\n" 
            + "    temp       "+ temp + "\n" 
            + "    pH         " + pH + "\n"  
            + "    units      "+ units + "\n"
            + "    density    " + density+ "\n" 
            + "    Alkalinity " + alkalinity + "\n")
            
        for i in range(len(Components)):
            Solution = (Solution 
                        + "    " + Components[i][0] + "         " + Components[i][1] + "\n")
           
        Solution =  (Solution + "    -water     " + mass
                    + "\n\n" )
    
    return Solution

def Mixer(Name, Ratio_1 = 0.5):
    """
    This function generates a mixed solution based on the added components and their ratios to one another.
    Name == Name of newly mixed solution
    Ratio_1 and Ratio_2 are the amounts of each solution in proportion to each other.
    S_ means that the variable is a string.

    """
    Ratio_2 = 1 - Ratio_1 
    S_Ratio_1 = str(Ratio_1)
    S_Ratio_2 = str(Ratio_2)

    Mixed_sol = (
        "MIX " + Name + "\n"
        + "    1 " + S_Ratio_1 + "\n"
        + "    2 " + S_Ratio_2 + "\n"
        + "SAVE Solution 3"
        + "\n\n"
    )
    return Mixed_sol

def Equilibrate(Name, Mineral, CO2 = ["CO2", "-3.5", "10"], Target_SI = "0.0", Addsolid = "0.0", CO2eq = False):
    """
    This function calls and sets the equilibrium_phases block.
    Name == The name of the final solution.
    Mineral == Name of the mineral
    Target_SI == The desired saturation of the mineral
    Addsolid == The ammount of solid added to acheive this

    """
    if CO2eq == False:
        Eq = (
            "EQUILIBRIUM_PHASES " + Name + "\n"
            + "    " + Mineral + " " + Target_SI + " " + Addsolid + "\n"
            + "SAVE Solution " + Name
            + "\n"
        )
    elif CO2eq == True:
        Eq = (
            "EQUILIBRIUM_PHASES " + Name + "\n"
            + "    " + Mineral + " " + Target_SI + " " + Addsolid + "\n"
            + "    " + CO2[0]  + " " + CO2[1] + " " + CO2[2] + "\n"
            + "SAVE Solution " + Name
            + "\n"
        )

    return Eq

def Set_Rates(Name, Amount = 2, Minerals = ["Aragonite" , "Calcite"], Surface_Area = [100000, 100000], Molar_mass = [100.09, 100.09]):
    """
    The following sets up the RATES block.
    Name == The name of the rates block. I just number them.
    Amount == Number of minerals which need rates. It is assumed they follow a similar type of equation
    Minerals == The names of the minerals
    Surface_Area == The Area per gram of Solid
    Molar_mass == The molar mass of solid

    """
    ###The following initializes the countner and sets the header for the block
    i = 0
    Rates1 = "RATES " + Name + "\n"
    Rates2 = ""

    ###Iterates to make each segment of the rates block
    while i < Amount:
        mineral = Minerals[i]
        molmass = Molar_mass[i]
        SA = Surface_Area[i]
        TempRate = ("    " + str(mineral) + "\n" 
                + "    -start""" + "\n"
                + "    10 IF (TOTAL_TIME <= parm(4)) THEN k = 10^parm(1)" + "\n"
                + "    20 IF (TOTAL_TIME > parm(4)) THEN k = 10^parm(2)" + "\n"
                + "    30 n = parm(3)" + "\n"
                + "    40 A = " + str(SA)  +  " * " + str(molmass)  + ' * KIN("' + str(mineral) + '")' + "\n"
                + '    50 R = (k *A *  (1 - SR("' + mineral + '"))^n)' + "\n"
                + "    60 SAVE R * time" + "\n"
                + "    -end" + "\n"
        )
        Rates2 += TempRate
        i += 1
    Rates = Rates1 + Rates2 + "\n\n"
    return Rates

def Set_Kinetics(Name, Solution, Amount = 2, Minerals = ["Aragonite" ,"Calcite"], Formulas = ["CaCO3 1", "CaCO3 1"], masses = ["0", "0"],
                 Orders = ["3","5"], EqConstants1 = ["-8.34","-8.48"], EqConstants2 = ["-10.19","-10.33"], ChngTime = "5040", Tolerance = "1e-08",
                 Time = "100", Steps = "1000", unit = "hour"):
    """
    This function sets up and initializes the kinetics block.
    Name == Name of this block, I use numbers
    Solution == The solution the kniteics occurs in (REQUIRED TO RUN MAKE SURE TO DEFINE)
    Amount == The number of minerals which will react
    Minerals == The names of the minerals
    Formulas == The Chemical formula for the compounds
    mass = mass of solid
    Orders == the orders of each reactions
    EqConstants1 == the Constants for the first part of the reaction
    EqConstants2 == The constants for the second part of the reaction
    ChngTime == The time which the first and second parts of the reaction change.
    Tolerance == The allowed error
    Time == Timespan of the reaction
    Steps == number of steps
    unit == units of time
    """
    
    ###Initialize the counter and the first line of the block
    i = 0
    Kinetics1 = "KINETICS " + Name + "\n"
    Kinetics2 = "" 

    ###Define the loop
    while i < Amount:
        mineral = Minerals[i]
        order = Orders[i]
        EqConstant1 = EqConstants1[i]
        EqConstant2 = EqConstants2[i]
        Formula = Formulas[i]
        mass = masses[i]

        tempkinetics = ("    " + str(mineral) + "\n"
                    + "        -formula  " + str(Formula) + "\n"
                    + "        -m0       " + str(mass) + "\n"
                    + "        -parms    " + str(EqConstant1) + " " + str(EqConstant2) + " " + str(order) + " " + ChngTime + "\n"
                    + "        -tol      " + Tolerance + "\n" 
        )
        Kinetics2 += tempkinetics
        i += 1

    Kinetics3 = ( "-steps    " + Time + " " + unit + " " + "in " + Steps + "\n"
             + "-step_divide 1" + "\n"
             + "-runge_kutta 6" + "\n"
             + "-bad_step_max 500" + "\n"
             + "INCREMENTAL_REACTIONS true" + "\n"
             + "\n"
             + "USE Solution " + Solution  + "\n"
    )

    Kinetics = Kinetics1 + Kinetics2 + Kinetics3

    return Kinetics

def Generate_Input(Solution_1 = "", Solution_2 = "", Solution_3 = "", Eq1 = "", Eq2 = "", Rates = "", Kinetics = ""):
    input_string = (
    Solution_1 + 
    Eq1 + 
    Solution_2 + 
    Eq2 +
    Solution_3 +
    """END \n\n"""
    + Rates
    + Kinetics +
    """

    SELECTED_OUTPUT 1
        -high_precision       true
        -reset                false
        -time                 true
        -saturation_indices   Aragonite  Calcite
        -kinetic_reactants    Aragonite  Calcite
        -pH                   true
    END
    
    """)

    return input_string

def main(Minerals, pHs, temps, Names, units, Wmasses, Smasses, CO2s, CO2eq, alkalinity,
         Eqconstants1, Eqconstants2, Orders, chngTime,  Surface_Areas, Components1, Components2, presets, Molar_masses, 
         Time, Steps, LUnit, Tolerance, xLimits, yLimits, Limits, Titles, Labels):
    
    start_time = tim.time()
    #Define the units of time
    t_unit = "hours"
    
    """
    Main function to execute the PHREEQC simulation and plot the results.
    """
    plt.style.use('default')
    ### Make and mix Solutions
    Solution_1 = Generate_Solutions(Names[0], temp = temps[0], pH = pHs[0], units = units[0], alkalinity = alkalinity[0], mass = Wmasses[0], preset = presets[0], Components = Components1)
    Solution_2 = Generate_Solutions(Names[1], temp = temps[1], pH = pHs[1], units = units[1], alkalinity = alkalinity[0], mass = Wmasses[1], preset = presets[1], Components = Components2)
    Solution_3 = Mixer(Names[0], Ratio_1 = 0.5)
    Eq1 = Equilibrate(Names[0], Minerals[0], Addsolid = Smasses[0], CO2 = CO2s[0], CO2eq = CO2eq[0])
    Eq2 = Equilibrate(Names[0], Minerals[1], Addsolid = Smasses[1], CO2 = CO2s[1], CO2eq = CO2eq[1])

    ### Set up the kinetics
    Rates = Set_Rates(Names[0], Amount = len(Minerals), Minerals = Minerals, Surface_Area = Surface_Areas, Molar_mass = Molar_masses)
    Kinetics = Set_Kinetics(Name = Names[0], Solution = "3", Amount = 2, Minerals = Minerals, masses = Kmasses,
                            Orders = Orders, EqConstants1 = Eqconstants1, EqConstants2 = Eqconstants2, ChngTime = chngTime, 
                            Tolerance = Tolerance, Time = Time, Steps = Steps, unit = LUnit)

    input_string = Generate_Input(Solution_1 = Solution_1, Solution_2 = Solution_2, Solution_3 = Solution_3, Eq1 = Eq1, Eq2 = Eq2, Rates = Rates, Kinetics = Kinetics)
    phreeqc = run_phreeqc(input_string)
    time, pH, SI_ara, SI_cal, Ara, dAra, Cal, dCal = extract_data(phreeqc, junkrow = 1)
    plot_data(time, SI_ara, SI_cal, Ara, dAra, Cal, dCal, pH, t_unit, xLimits, yLimits, Limits, Titles, Labels, TTime = float(Time), TSteps = float(Steps))

    # Show the plots
    plt.show()
    end_time = tim.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} seconds")
    return time, pH, SI_ara, SI_cal, Ara, dAra, Cal, dCal, input_string  

if __name__ == '__main__':
    
    ###Solution information
    Minerals = ["Aragonite", "Calcite"]
    pHs = ["6.5","6.5"]
    temps = ["18", "18"]
    Names = ["1", "2", "3"]
    units = ["mmol/L", "mmol/L"]
    Wmasses = ["0.01", "0.01"] #Water masses
    Smasses = ["0.0004995", "0.0004995"] #Masses of Solids for equlibrate function
    CO2s = [["CO2(g)", "-3.5", "10"], ["CO2(g)", "-3.5", "10"]]
    CO2eq = [True, True]
    alkalinity = ["5", "5"]
    Components1 = [] #List of minerals followed by their concentration
    Components2 = []   # Eg.  [F, 10] is flouride with a concentration of 10
    presets = ["", ""]


    ###Kinetics information
    Eqconstants1 = ["-8.34","-8.48"]
    Eqconstants2 = ["-10.19","-10.33"]
    Orders = ["3", "5"]
    chngTime = "5040" #Time when the Eq. Constants change
    Surface_Areas = [100000, 100000]
    Molar_masses = [100.09, 100.09]
    Kmasses = ["0.0004995", "0.0002575"]

    ###Loop information
    Time = "150"
    Steps = "10000"
    LUnit = "hour" #Units of time
    Tolerance = "1e-08"

    ###Plotting information
    xLimits = [0, float(Time) - 140]
    yLimits = [0, 0]
    Limits = [False, False]
    Titles = ["Aragonite and Clacite Over Time", "Aragonite and Clacite Over Time", "Saturations Over Time", "Relative Levels Over Time", "pH of Solution", "Rate of Change Over Time"]
    Labels = ["Amount of Solid (mol)","Amount of Solid (mol)", "Saturation Index", "Mole Fraction", "pH", "Rate of Dissolution (mol/h)"]

    (time, pH, SI_ara, SI_cal, Ara, dAra, Cal, dCal, input_string) = main(Minerals, pHs, temps, Names, units, Wmasses, Smasses, CO2s, CO2eq, alkalinity,
                                                                          Eqconstants1, Eqconstants2, Orders, chngTime, Surface_Areas, Components1, Components2, presets,
                                                                          Molar_masses, Time, Steps, LUnit, Tolerance,
                                                                          xLimits, yLimits, Limits, Titles, Labels )
