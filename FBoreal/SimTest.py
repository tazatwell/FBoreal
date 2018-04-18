#from __future__ import division
######################################################################################################################################
#
#                FireSimulator (Split) serial and parallel version 0.05 - April 2018 
#                Author: Cristobal Pais; modifications by David L. Woodruff
#                example: mpiexec -n X python Path\Simulator1Beta.py  where X is the number of parallel processes
#
######################################################################################################################################
"""
CP Apr 2018: Previous version before using core computations in C
    Code has been cleaned (all classes and main), numerical operations have been optimized,
    FI and FS are being modified using dictionaries/sparse arrays, not
    needed NOW. Python 3.X version (no longer 2.7).
    Paper draft has been updated. 
    Parallel versions MPI and CUDA GPU are being tested using pure Python, pure C, and
    C wrapped in Python. Massive upload coming at the end of this month.
    ReadData.py renamed to FBP2PY.py
    Heuristic class is being tested: 5 heuristics + RDLearning algorithm
    All Parallel mentions were eliminated from this version
DLW Feb 2018 Note re: reading names. Caution: some of the names are "backwards"
    from the perspective of a native speaker of English, but some are not.
    Unrelated: FI and FS for sceanios: We need to stop
    using so much memory by using indicators (mostly zero) and we need
    to support a scenario time step that is different from the fire time step.
DLW Jan 2018 Notes re: Burnt/Burned and Burning:
    A cell that starts burning will be burned, so it could be
    listed as burned right away. I think "get_burned" was
    intended to mean the the fire starts based on a message.
NOTE: The paper uses "fuel cell" while the code often uses "fire cell"
    for the same thing.
DLW March 2018: changing Weather file format and weather to be a
    data frame; access to weather will be through a single object
    that knows about time. Also changing time tracking vis a vis 
    hours, weeks, etc.
March 2018: TBD: random weather needs a lot of rework
"""

__version__ = 0.05

# Importations
import os
import shutil
import numpy.random as npr
import numpy as np
from itertools import repeat
# from mpi4py import MPI
import time
import ctypes
from argparse import ArgumentParser

#Class importations
import WeatherFBP
import CellsFBP
import Forest 
import Lightning
import Plot
import FBP2PY
import ReadDataPrometheus
import SpottingFBP
import Output_Grid

#Excel Importations
from openpyxl import *

"""
deleted by dlw Feb 2018
# Input and outputfile names from the working directory (where the python script is, 
  otherwise, insert the Path inside the os.chdir function)
cwd = os.path.dirname(os.path.realpath(__file__))
os.chdir(cwd)
"""

# Shared library .so (CP: Careful with Windows vs UNIX)
soname = "FBPfunc5NODEBUG.so"
try:
    lib = ctypes.cdll.LoadLibrary(soname)
except:
    raise RuntimeError("Could not load the library="+soname)

        
# Main class with all the simulation scheme
class Main:        
######################################################################################################################################
#    
#        Simulator: All the simulation scheme is developed here
#
######################################################################################################################################

###Step -1: FBP
    # FType coefficients and data from FBP library
    listlen = 18
    ListOfCoefs = listlen * CellsFBP.fuel_coeffs   
    coef_ptr = ListOfCoefs()
    lib.setup_const(coef_ptr)
    FTypes2 = {"m1":0, "m2":1,"m3":2,"m4":3,
               "c1":4,"c2":5,"c3":6,"c4":7,"c5":8,"c6":9,"c7":10,
               "d1":11,"s1":12,"s2":13,"s3":14,"o1a":15,"o1b":16,"d2":17}
    
    # Forest Data (read table using c functions in class ReadData)
    parser = ArgumentParser()
                    
    parser.add_argument("--input-forest",
                        help="The name of the csv file that contains the forest grid data",
                        dest="input_forest",
                        type=str,
                        default=None)                
                    
    parser.add_argument("--input-fpblookup",
                        help="The name of the csv file that contains the forest grid data",
                        dest="input_fbplookup",
                        type=str,
                        default=None) 
    
    parser.add_argument("--input-instance-folder",
                        help="The path to the folder contains all the files for the simulation",
                        dest="input_folder",
                        type=str,
                        default=None)                
                    
    parser.add_argument("--output-folder",
                        help="The path to the folder for simulation output files",
                        dest="output_folder",
                        type=str,
                        default=None)                
                    
    parser.add_argument("--seed",
                        help="Seed for random numbers (default is None, which randomizes)",
                        dest="seed",
                        type=int,
                        default=None)        
    
    parser.add_argument("--weather",
                        help="The 'type' of weather: constant, random, rows (default rows)",
                        dest="input_weather",
                        type=str,
                        default="rows")
        
    parser.add_argument("--spotting-parameter-data-file-name",
                        help="Spotting parameter data file name (default None)\n  As of Feb 2018 a JSON file and the dictionary should have 'SPTANGLE', 'SPOT0PROB', and 'SPOT10TIME'",
                        dest="spotting_parameter_data_file_name",
                        type=str,
                        default=None)                
                    
    parser.add_argument("--plot",
                        help="Plots are created after the simulation",
                        dest="input_plottrue",
                        default=False,
                        action='store_true')    
    
    parser.add_argument("--verbose",
                        help="Output all the simulation log",
                        dest="verbose_input",
                        default=False,
                        action='store_true')

    parser.add_argument("--ignition-point",
                        help="The name of the csv file that contains the ignition data",
                        dest="input_ignitions",
                        type=str,
                        default="")        

    parser.add_argument("--ignitions",
                        help="Activates the predefined ignition points when using the folder execution",
                        dest="act_ignitions",
                        default=False,
                        action="store_true")    

    parser.add_argument("--sim-years",
                        help="Number of years per simulation (default 1)",
                        dest="sim_years",
                        type=int,
                        default=1)    
    
    parser.add_argument("--nsims",
                        help="Total number of simulations (replications)",
                        dest="input_nsims",
                        type=int,
                        default=1)    
                    
    parser.add_argument("--Fire-Period-Length",
                        help="Fire Period length in minutes (needed for ROS computations). Default 60",
                        dest="input_PeriodLen",
                        type=float,
                        default=60)                    
    
    parser.add_argument("--Weather-Period-Length",
                        help="Weather Period length in minutes (needed weather update). Default 60",
                        dest="weather_period_len",
                        type=float,
                        default=60)                    
    
    parser.add_argument("--ROS-Threshold",
                        help="A fire will not start or continue to burn in a cell if the head ros is not above this value (m/min) default 0.1.",
                        dest="ROS_Threshold",
                        type=float,
                        default=0.1)                    
    
    parser.add_argument("--HFI-Threshold",
                        help="A fire will not start or continue to burn in a cell if the HFI is not above this value (Kw/m) default is 10.",
                        dest="HFI_Threshold",
                        type=float,
                        default=10.0)                    
    
    parser.add_argument("--ROS-CV",
                        help="Coefficient of Variation for normal random ROS (e.g. 0.13), but default is 0 (deteriministic)",
                        dest="ROS_CV",
                        type=float,
                        default=0.0)                    
    
    parser.add_argument("--heuristic",
                        help="Harvesting heuristic/plan used (number)",
                        dest="input_heur",
                        type=int,
                        default=0)        
    
    parser.add_argument("--statistics",
                        help="Excel file with statistics is created at the end of the simulation",
                        dest="input_excel",
                        default=False,
                        action="store_true")    
                    
    parser.add_argument("--combine",
                        help="Combine fire evolution diagrams with the forest background",
                        dest="input_combine",
                        default=False,
                        action="store_true")                    
                    
    parser.add_argument("--save-memory",
                        help="Activates memory monitoring version (useful for very large instances)",
                        dest="input_save",
                        default=False,
                        action="store_true")                

    parser.add_argument("--scenarios-out",
                        help="Activates scenario output files; requires --output-grid",
                        dest="input_scenario",
                        default=False,
                        action="store_true")        
    
    parser.add_argument("--no-output",
                        help="Activates no-output mode ",
                        dest="no_output",
                        default=False,
                        action="store_true")        
                    
    parser.add_argument("--output-grid",
                        help="Activates forest grid's output mode (TEMPORARILY DISABLED by DLW)",
                        dest="output_grid",
                        default=False,
                        action="store_true")        
                    
    parser.add_argument("--FBP-tester-cell",
                        help="If not 0, Causes special output, then termination after one ROS calculation for the cell number given (default 0)",
                        dest="FBP_tester_cell",
                        type=int,
                        default=0)        
                    
    parser.add_argument("--max-fire-periods",
                        help="Maximum fire periods per year (default 2016)",
                        dest="max_fire_periods",
                        type=int,
                        default=2016)        
                    
    args = parser.parse_args()

    ######### a little "local" report writer #################
    '''
    Sim            int
    AreaCells      array of doubles
    NCells         int
    FireP_Period   int
    ASet           int set
    BSet           int set
    NSet           int set
    HSet           int set
    '''
    def special_report(args, Sim, AreaCells, NCells,
                       Fire_Period, ASet, BSet, NSet, HSet = None):
        """
        Used for quick analysis of output and for regression tests.
        Writes a self-documenting csv file
        """
        def writeone(appfile, name, val):
            appfile.write(", "+ name+ ", "+ str(val))

        fname = os.path.join(args.output_folder, 
                             "SimReport.csv")
        
        with open(fname, "a") as appfile:
            appfile.write("infolder, "+ args.input_folder )
            writeone(appfile, "ReplicateNum", Sim)
            writeone(appfile, "seed", args.seed)
            writeone(appfile, "ROSCV", args.ROS_CV)
            writeone(appfile, "CellArea", AreaCells)
            writeone(appfile, "NCells", NCells)
            writeone(appfile, "Burntcells", len(BSet))
            writeone(appfile, "NonBurnable", len(NSet))
            writeone(appfile, "StillAvailableCells", len(ASet))
            
            if HSet is not None:
                writeone(appfile, "HarvestedCells", len(HSet))
            writeone(appfile, "SimYears", args.sim_years)
            
            for fi in range(len(Fire_Period)):
                writeone(appfile, "fireperiods"+str(fi+1), Fire_Period[fi])
            appfile.write("\n")

    InFolder = args.input_folder
    OutFolder = args.output_folder
    SaveMem = args.input_save
    scenarios = args.input_scenario
    fileweather = True  # delete fileweather and this line...
    nooutput = args.no_output
    OutputGrid = args.output_grid
    MinutesPerWP = args.weather_period_len

    # Sanity checks (inputs)
    if scenarios and not OutputGrid:
        raise RuntimeError("--scenarios requires --output-grid")
    if InFolder != None:
        filename = os.path.join(InFolder, "Data.dat")
    else:
        filename = "Data.dat"

    DF = ReadData.inputData(filename)
    ##print "debug df after first read"
    ##print DF    

    if args.FBP_tester_cell > 0:
        cellID = args.FBP_tester_cell
        print ("Entering FBP_tester; trying cell=", cellID)
        mainstruct, headstruct, flankstruct, backstruct = FBP2PY.CalculateOne(DF, coef_ptr,cellID, verbose=True)
        quit()
    
    if os.path.exists(OutFolder):
        print ("Creating a new, empty folder=",OutFolder)
        shutil.rmtree(OutFolder)
    os.makedirs(OutFolder)
    
    #Getting FType for each cell from data 
    FTypeCells2 = []
    for i in DF['fueltype']:
        FTypeCells2.append(i)

    verbose = args.verbose_input
    plottrue = args.input_plottrue
    max_weeks = 12 # weeks
    Max_Fire_Periods = args.max_fire_periods

    TotalYears= args.sim_years
    TotalSims = args.input_nsims
    Sim = 1

    if verbose:
        print ("Setting initial pseudo-random number stream seed to ", str(args.seed))
    if args.seed is not None:
        npr.seed(int(args.seed))
    
    ###################################################
    ''' Modification October 2017
        Period length: in seconds for the moment
        Feb 2018: minutes
    '''
    FirePeriodLen = args.input_PeriodLen
    print ("-----------------------------------------------------------------------------------")
    print ("Fire Period Length for ROS computations [min]: ", FirePeriodLen)
    print ("-----------------------------------------------------------------------------------")
    '''
        Period length: in whatever the ROS is from FBP (I think minutes)
    '''
    ###################################################
    
    #    Simulations loop
    while Sim <= TotalSims:
        ### Step 0: Initializing Instances, classes, global parameters, sets, etc.
        # No output dominates verbose
        if nooutput == True:
            verbose = False
        
        #Global Parameters
        plottrue = args.input_plottrue
        week_number = 1 
        TotalYears= args.sim_years
        TotalSims = args.input_nsims
        Year = 1
        
        # Current fire period in a year (also records last fire period)
        Fire_Period = np.zeros(TotalYears)  #[0 for i in range(TotalYears)]
        NoIgnition=None
        MessagesSent=None
        plotnumber = 1
        
        print( "------------------------------------------------- Simulation Number ",Sim,"-------------------------------------------------" )
        '''
        CP April 2018: Outputgrid data structure is being modified (optimized), currently disabled
        '''
        if OutputGrid:
            FI = [0]*Max_Fire_Periods
            HPeriod = [0]
            
        # BurntP is an array (by period) of the array of cells burned in the period
        BurntP = [[] for i in repeat(None, Max_Fire_Periods)]
        
        Initial_Time = time.clock()
        if verbose ==  True:
            print("Initial Time:", Initial_Time)
            
        #Read Forest
        if nooutput == False:
            print ("\n----------------- Forest Data -----------------")
        
        # Get the arguments
        ForestFile = args.input_forest
        FBPlookup = args.input_fbplookup
        Ignitions = args.input_ignitions
        Thresholds = {}
        
        # If we have a folder with an instance read files
        if InFolder != None:
            ForestFile = os.path.join(InFolder, "Forest.asc")
            FBPlookup = os.path.join(InFolder, "fbp_lookup_table.csv")
            Ign = args.act_ignitions
            
            if Ign == True:
                Ignitions = os.path.join(InFolder, "IgnitionPoints.csv")
                
            if args.spotting_parameter_data_file_name is not None:
                SpottingFile = os.path.join(InFolder, args.spotting_parameter_data_file_name)
                SpottingParams = ReadDataPrometheus.ReadSpotting(SpottingFile,nooutput)
                
            else:
                SpottingParams = None
                
        '''
        CP April 2018: ReadDataPrometheus has been optimized in new version coming at the end of the month
        '''
        FBPDict,ColorsDict =  ReadDataPrometheus.Dictionary(FBPlookup)
        CellsGrid3,CellsGrid4,Rows,Cols,AdjCells,CoordCells,CellSide = ReadDataPrometheus.ForestGrid(ForestFile,FBPDict)
            
        # Check ignitions args
        if Ignitions != "":
            Ignitions = ReadDataPrometheus.IgnitionPoints(Ignitions)
            if nooutput == False:
                print("We have specific ignition points")
                print("Ignitions:", Ignitions)
        if Ignitions == "" and nooutput == False:
            print( "No ignition points")
            
        if nooutput == False:    
            print ("Rows:",Rows,"Cols:",Cols, "NCells:",Rows*Cols)
            
        # Initializing empty arrays
		FTypeCells = []
        StatusCells = []
        Colors = []
        RealCells = []
        cellcounter=1
            
        for i in range(0,len(CellsGrid4)):
            if str.lower(CellsGrid4[i]) not in FTypes2.keys():
                FTypeCells.append(0)
                StatusCells.append(4)
                CellsGrid4[i] = "s1"
                RealCells.append(0)
            else:
                FTypeCells.append(2)
                StatusCells.append(0)
                RealCells.append(cellcounter)
                cellcounter+=1
        
            if str(CellsGrid3[i]) not in ColorsDict.keys():
                Colors.append((1.0,1.0,1.0,1.0))
                
            if str(CellsGrid3[i]) in ColorsDict.keys():
                Colors.append(ColorsDict[str(CellsGrid3[i])])
        
        #Releasing memory
        del CellsGrid3
        del ColorsDict
            
        if nooutput == False:        
            print ("------------ End read forest data -------------")
        
        #Cells list (DLW says: yes, but this is dictionary....)
        # DLW, Feb 2018: Not sure what this has in it.... Seems to be just active cells
        # CP, April 2018: Dictionary with active cells, we "don't care" about not active cells. Used to be a list...
        Cells_Obj = {}
                
        # cell instance data (explicit for the moment, final version may read it from a source file)
        # CP: Heuristics will populate VolCells, AgeCells, and monetary parameters from a file
        #     kept here for future reference
        VolCells = 100 # ???  
        AgeCells = 1 # ???
        TerrainCells = 0 # ??? CP: deprecated
        AreaCells = CellSide * CellSide
        PerimeterCells = CellSide * 4 # ???
        
        #Main parameters (Forest)
        NCells = Rows*Cols
        IDF = 1
        Location = "My_Mind"
        Coord = [1,1,100]
        Area = NCells * AreaCells
        Vol = NCells * VolCells # ??? CP April 2018: see previous comment
        Age = 12.0 # ???
        # HEY TBD DLW Feb 2018 we should not assume a square forest
        # CP April 218: MODIFIED, however, not very important...
        Perimeter = 2 * CellSide * (Rows + Cols) # ???
        FTypes = {0: "NonBurnable", 1: "Normal", 2: "Burnable"}
        
        # Weather options
        WeatherOpt = str.lower(args.input_weather)
        wchoices = ['constant', 'random', 'rows']
        if WeatherOpt not in wchoices:
            print ("Valid weather choices are:",str(wchoices))
            raise RuntimeError ("invalid choice for --weather") 

        weatherperiod = 0
        if nooutput==False:
            print( "Reading weather file")
        Weather_Obj = WeatherFBP.Weather(os.path.join(args.input_folder, "Weather.csv"))
        DF = Weather_Obj.update_Weather_FBP(DF, WeatherOpt, weatherperiod)

        if verbose == True:
            Weather_Obj.print_info(0)

        if verbose == True:
            print( "DF", DF[["ws","waz","ps", "saz"]]    )

        """
        CP April 2018: Forest data will be populated from instance in next version, 
        currently being implemented.
        DLW, Feb 2018, Forest1 is not currently used. TBD: it should be initialized
        from the cell data, or at least it should match cell data.
        #Initializing Forest Instance
        Forest1 = Forest.Forest(IDF,Location,Coord,NCells,Area,Vol,Age,Perimeter,FTypes)
        
        if verbose == True:
            print "\n" + "------------------------ Forests lists ------------------------"
            Forest1.print_info()
        """
        
        #Initializing List of Cells and weather (printing info, if verbose)
        if verbose == True:
            Weather_Obj.print_info(weatherperiod)
        
        #Initializing plot object and plot the initial forest
        if plottrue == True:
            PlotPath = os.path.join(OutFolder, "Plots")
            if not os.path.exists(PlotPath):
                print("creating",PlotPath)
                os.makedirs(PlotPath)
            Plotter = Plot.Plot()
            emptylist = [[] for i in range(0,NCells)]
            
            if os.path.isfile(os.path.join(OutFolder, "ForestInitial.png")) == True:
                if nooutput == False:
                    print("Forest already exists")
            else:    
                Plotter.PlotForestOnly(Colors,CoordCells,plotnumber,0,Year,False,Rows,Cols,OutFolder)
                        
        #Initializing Lambda Instance (random lightning)
        Lambda_Strike = Lightning.Lightning()
            
        #Sets (Available cells, burning cells, burnt cells and harvest cells: starting from array and cast later)
        avail=[]
        nonbur=[]
        
        # Cell status: populate arrays before cast sets (EP)
        for i in range(0, NCells):
            if StatusCells[i] != 4:
                avail.append(i+1)
            if StatusCells[i] == 4:
                nonbur.append(i+1)
                
        # Cast to sets
        AvailCells_Set = set(avail)
        del avail
        NonBurnableCells_Set = set(nonbur)
        del nonbur
        BurningCells_Set = set()
        BurntCells_Set = set()
        HarvestCells_Set = set()
        
        #Printing info about the cells' status
        if verbose == True:
            print("\n" +"Set information period (week"+str(week_number)+")")
            print("Available Cells: " + str(AvailCells_Set))
            print("Non Burnable Cells:" +str(NonBurnableCells_Set))
            print("Burning Cells: " + str(BurningCells_Set))
            print("Burnt Cells: " + str(BurntCells_Set))
            print("Harvest Cells: " + str(HarvestCells_Set))
                
        # Spotting information
        spotting = args.spotting_parameter_data_file_name is not None
        
        
        ###Years' loop (start)
        while Year <= TotalYears:
            ### Step 1: Lightning/Ignition loop in order to find the first week with fire
            #Loop for finding next fire week
            
            #Printing information if verbose is true
            if verbose == True:
                print("\n","------------------------------------ Current Year: ",Year, "------------------------------------")
                print( "---------------------- Step 1: Ignition ----------------------")
                
            if SaveMem == True:
                for c in BurntCells_Set:
                    if (c-1) in Cells_Obj.keys():
                        if verbose == True:
                            print ("Deleting burnt cells from list---------------------------------------------")
                            print ("Deleted:",c)
                        del Cells_Obj[c-1]
                        
            #Parameters and variables
            aux=0
            loops=0
            NoIgnition = False
                
            # Starting fire week 
            if Ignitions == "":
                while week_number <= max_weeks:
                    #If a lightning occurs in a week, we select that week
                    if Lambda_Strike.Lambda_NH(week_number,verbose) == True:
                        Sel_Week = week_number
                        if verbose == True:
                            print("Selected Week: " + str(Sel_Week))
                        break
                    else: 
                        week_number+=1
                        weatherperiod += 1
                    
            else: 
                week_number=1 
                weatherperiod=0
                Sel_Week = 1
                
            #Go to that period/week
            ### DLW Feb 2018: deal with weather "weeks" then delete this comment
            week_number = Sel_Week
            print( "week_number=", Sel_Week)
            
            if verbose == True:
                print("\nCurrent Period (Week): " + str(week_number))
                print("Current Fire Period (Hour): " + str(Fire_Period[Year-1]))
                            
            # Select the "burning" cell
            aux=0
            loops=0
            NoIgnition = False
                            
            # Select the "burning" cell
            if Ignitions == "":
                while True:
                    aux = int(npr.uniform(0,1) * NCells)    
                    if StatusCells[aux] != 4 and (aux+1) not in BurntCells_Set:
                        # Initialize cell if needed
                        if aux not in Cells_Obj.keys():
                            Cells_Obj[aux] = CellsFBP.Cells((aux+1),AreaCells, CoordCells[aux],
                                                            AgeCells, FTypeCells[aux],
                                                            coef_ptr[FTypes2[str.lower(CellsGrid4[aux])]],
                                                            TerrainCells, VolCells, PerimeterCells,
                                                            StatusCells[aux], AdjCells[aux], Colors[aux],
                                                            RealCells[aux], OutputGrid)
                            '''
                                
                            '''
                            # Avail set modification for initialization
                            Cells_Obj[aux].InitializeFireFields(CoordCells, AvailCells_Set)
                                
                            '''
                        
                            '''
                        
                        
                        #If the cell is available and burnable
                        if Cells_Obj[aux].get_Status() == "Available" and Cells_Obj[aux].FType != 0:
                            if Cells_Obj[aux].ignition(Fire_Period[Year-1],Year,Ignitions,DF,coef_ptr,args.Start_Threshold, args.HFI_Threshold):
                               
                                #If we have an ignition, change the status of the forest
                                if OutputGrid:
                                    FI[Fire_Period[Year-1]-1]=Cells_Obj[aux].ID
                                    BurntP[Fire_Period[Year-1]-1]=Cells_Obj[aux].ID
                                
                                for i in range(Fire_Period[Year-1],Max_Fire_Periods):
                                    BurntP[i] = Cells_Obj[aux].ID

                                # Printing info about ignitions        
                                if verbose == True:
                                    print("Cell "+str(Cells_Obj[aux].ID)+" Ignites")
                                    print("Cell "+str(Cells_Obj[aux].ID)+" Status: " + Cells_Obj[aux].get_Status()    )
                                
                                if plottrue == True and SaveMem == True:
                                        Plotter.forest_plotV3_FreeMem(Cells_Obj,emptylist,plotnumber,Fire_Period[Year-1],Year,False,Rows,Cols,PlotPath,CoordCells,BurntCells_Set,Sim)
                                        plotnumber+=1
                                        
                                if plottrue == True and SaveMem != True:
                                    Plotter.forest_plotV3(Cells_Obj,emptylist,plotnumber,Fire_Period[Year-1],Year,False,Rows,Cols,PlotPath,Sim)
                                    plotnumber+=1    
                                        
                                break
                                    
                    #Updating parameters inside the loop
                    loops+=1
                        
                    #Maximum number of iterations
                    if loops>NCells:
                        NoIgnition = True
                        break
                    
            if Ignitions != "":
                if Ignitions[Year] not in BurntCells_Set and StatusCells[Ignitions[Year]-1]!=4:
                    if (Ignitions[Year]-1) not in Cells_Obj.keys():
                        Cells_Obj[Ignitions[Year]-1] = CellsFBP.Cells((Ignitions[Year]),AreaCells,CoordCells[Ignitions[Year]-1],AgeCells,FTypeCells[Ignitions[Year]-1],coef_ptr[FTypes2[str.lower(CellsGrid4[Ignitions[Year]-1])]],TerrainCells,VolCells,PerimeterCells,StatusCells[Ignitions[Year]-1],AdjCells[Ignitions[Year]-1],Colors[Ignitions[Year]-1],RealCells[Ignitions[Year]-1], OutputGrid)
                            
                            '''
                            
                            '''
                        Cells_Obj[Ignitions[Year]-1].InitializeFireFields(CoordCells, AvailCells_Set)
                            #Cells_Obj[Ignitions[Year]-1].InitializeFireFields(CoordCells)
                        
                            '''
                            
                            
                            '''
                        
                        
                    if Cells_Obj[Ignitions[Year]-1].get_Status() != "Available" or Cells_Obj[Ignitions[Year]-1].FType == 0:
                        NoIgnition = True
                        
                    if Cells_Obj[Ignitions[Year]-1].get_Status() == "Available" and Cells_Obj[Ignitions[Year]-1].FType != 0:
                        if Cells_Obj[Ignitions[Year]-1].ignition(Fire_Period[Year-1],Year, Ignitions,DF,coef_ptr,args.ROS_Threshold, args.HFI_Threshold):
                            if OutputGrid:                                
                                FI[Fire_Period[Year-1]-1]=Cells_Obj[Ignitions[Year]-1].ID
                            BurntP[Fire_Period[Year-1]-1]=Cells_Obj[Ignitions[Year]-1].ID
                            
                            for i in range(Fire_Period[Year-1],Max_Fire_Periods):
                                BurntP[i] = Cells_Obj[Ignitions[Year]-1].ID
                                        
                            # Printing info about ignitions        
                            if verbose == True:
                                print( "Cell "+str(Cells_Obj[Ignitions[Year]-1].ID)+" Ignites")
                                print( "Cell "+str(Cells_Obj[Ignitions[Year]-1].ID)+" Status: " + Cells_Obj[Ignitions[Year]-1].get_Status()    )
                            if plottrue == True and SaveMem == True:
                                Plotter.forest_plotV3_FreeMem(Cells_Obj,emptylist,plotnumber,Fire_Period[Year-1],Year,False,Rows,Cols,PlotPath,CoordCells,BurntCells_Set,Sim)
                                plotnumber+=1
                            if plottrue == True and SaveMem != True:
                                Plotter.forest_plotV3(Cells_Obj,emptylist,plotnumber,Fire_Period[Year-1],Year,False,Rows,Cols,PlotPath,Sim)
                                plotnumber+=1    
                        
                else:
                    NoIgnition = True
                    if nooutput == False:
                        print( "No ignition during year "+str(Year)+", cell "+str(Ignitions[Year])+ " is already burnt or non-burnable type")
                            
            #If ignition occurs, we update the forest status
            if NoIgnition == False:
                #Updating AvailCells and BurningCells sets    
                if Ignitions == "":
                    NewID = Cells_Obj[aux].ID
                    Aux_set = set([NewID])
                    BurningCells_Set = BurningCells_Set.union(Aux_set)
                    AvailCells_Set = AvailCells_Set.difference(BurningCells_Set)
                    
                else:
                    NewID = Cells_Obj[Ignitions[Year]-1].ID
                    Aux_set = set([Ignitions[Year]])
                    BurningCells_Set = BurningCells_Set.union(Aux_set)
                    AvailCells_Set = AvailCells_Set.difference(BurningCells_Set)

                    BurntCells_Set = BurntCells_Set.union(Aux_set)

            #Printing information about the forest
            if verbose == True:
                print( "Available cells: " + str(AvailCells_Set))
                print( "Non Burnable Cells:" +str(NonBurnableCells_Set))
                print( "Burning cells: " + str(BurningCells_Set))
                print( "Harvest cells: " + str(HarvestCells_Set))
                print( "Burnt cells (including burning): " +str(BurntCells_Set))
                                
            #Next Period: t=t+1. Update Weather
            Fire_Period[Year-1]+=1
            print( "debug Fire_Period[Year-1]", Fire_Period[Year-1])
            
            if WeatherOpt != 'constant' and Fire_Period[Year-1] * FirePeriodLen / MinutesPerWP > weatherperiod:
                weatherperiod +=1
                DF = Weather_Obj.update_Weather_FBP(DF, WeatherOpt, weatherperiod)
                
            if verbose == True:
                print( "Current Week: "+str(week_number))
                print( "Fire Period Starts: " + str(Fire_Period[Year-1]))
                print( "")
                Weather_Obj.print_info(weatherperiod)
                if WeatherOpt == 'constant':
                    print( "(NOTE: current weather is not used for ROS with constant option)")
                            
                # End of the ignition step
                print( "")
                print( "Next Fire Period: " + str(Fire_Period[Year-1]))
                                            
            if week_number == 12 or Fire_Period[Year-1] == 2016:
                print( "-------------------------------------------------------------------------\n"+ "End of the fire year",Year,"-------------------------------------------------------------------------"    )
                    
			### End Step 1
        
			#If no ignition occurs, go to next year (no multiple ignitions per year, only one)
			if NoIgnition == True:
				if verbose == True:
					print( "No ignition in year",Year)
				#Next year, reset the week
				Year+=1    
				week_number=1
				
			#If an ignition happened, enter in fire dynamic loop
			if NoIgnition == False:
					
				### Steps 2 and 3: Send/Receive messages (fire dynamic)
				#Fire dynamic loop (fire periods)
				while Fire_Period[Year-1] < Max_Fire_Periods:
					if Fire_Period[Year-1] == Max_Fire_Periods-1:
						print( "\n **** WARNING!!!! About to hit Max_Fire_Periods=",Max_Fire_Periods," ***\n\n")
						
					### Step 2: Sending messages
					# Initial Parameters
					MessagesSent=False
					SendMessageList = [[] for i in repeat(None, NCells)]
						
					#Printing info
					if verbose == True:
						print( "\n---------------------- Step 2: Sending Messages from Ignition ----------------------")
						print( "Current Week: " + str(week_number))
						print( "Current Fire Period: "+str(Fire_Period[Year-1]))
						print( "Burning Cells: "+str(BurningCells_Set))
						print( "Burnt Cells (should include burning): "+str(BurntCells_Set))
							
					'''
					Cleaning ROSAngleDir dictionaries based on the current burning cells
					'''
						
					if verbose == True: 
						print( " -------------------- NEW CLEANING STEP ------------------------")
					
					# For all the cells already initialized
					for cell in Cells_Obj.keys():
						# Check those cells with at least one possible neighbor to send fire
						if len(Cells_Obj[cell].ROSAngleDir) > 0:
								
							# Delete adjacent cells that are not available 
							for angle in Cells_Obj[cell].angle_to_nb:
								nb = Cells_Obj[cell].angle_to_nb[angle]
									
								if nb not in AvailCells_Set and angle in Cells_Obj[cell].ROSAngleDir:
										
									if verbose == True: 
										print( "Cell", cell+1, ": clearing ROSAngleDir")
							
									Cells_Obj[cell].ROSAngleDir.pop(angle)
								
					if verbose == True: 
						print( "----------------------------------------------------------------------")
						
					   
					
					# RepeatFire
					RepeatFire = False    
					BurnedOutList = []
					for cell in BurningCells_Set:
						########################################################
						''' New burning logic  (EXTRA INPUT FirePeriodLen)
								
						'''
								
						if verbose == True:
							print("Cell object new fields")
							print("ID:",Cells_Obj[cell-1].ID)
							print("FireProgress:",Cells_Obj[cell-1].FireProgress)
							print("AngleDict:",Cells_Obj[cell-1].AngleDict)
							print("ROSAngleDir:",Cells_Obj[cell-1].ROSAngleDir)
							print("DistToCenter:",Cells_Obj[cell-1].DistToCenter)
							print("angle_to_nb:",Cells_Obj[cell-1].angle_to_nb)
							
						# Check if the burning cell can send more messages or not 
						if len(Cells_Obj[cell-1].ROSAngleDir) > 0:
							if verbose == True:
								print( "Cell", cell, "can still send messages to neighbors")
							Aux_List = Cells_Obj[cell-1].manageFire(Fire_Period[Year-1],AvailCells_Set, 
																	verbose,DF,coef_ptr,spotting,SpottingParams, 
																	CoordCells,Cells_Obj, args)
								
						else:
							if verbose == True:
								print( "Cell", cell, "does not have any neighbor available for receiving messages")
							Aux_List = []
							
								'''    
								'''
						 ########################################################
								
						# Debug
						if verbose == True:
							print( "Aux list:", Aux_List)
								
						# Original condition
						if len(Aux_List) > 0 and Aux_List[0] != "True" :
							if verbose == True:
								print( "List is not empty")
							MessagesSent=True
							SendMessageList[Cells_Obj[cell-1].ID-1] = Aux_List
							if verbose == True:
								print( "SendMessageList: ", SendMessageList)
									
															
						'''
						Major modifications in Oct 17.........
						'''
						if len(Aux_List) > 0 and Aux_List[0] == "True":
							if verbose == True:
								print( "Fire is still alive and we may repeat if no other messages......")
							RepeatFire = True
										
						if len(Aux_List) == 0:
							BurnedOutList.append(Cells_Obj[cell-1].ID)
							if verbose == True:
								print( "Message and Aux Lists are empty; adding to BurnedOutList")
					### End burn loop
					
				BurningCells_Set.difference(set(BurnedOutList))
				# Check the conditions for repeating, stopping, etc.... 
				# Main change......
				Global_Message_Aux = [val for sublist in SendMessageList for val in sublist]
				if verbose == True:
					print( "Global_Message_Aux:", Global_Message_Aux)
					print( "RepeatFire:", RepeatFire)
							
				# if we have at least one cell that neeeds repetition and no other messages exists....
				# We repeat!!!!
				if RepeatFire == True and len(Global_Message_Aux) == 0:
					# Update fire period 
					print( "\nFires are still alive and no message has been generated during this period")
					print( "Current Fire period: ",Fire_Period[Year-1])
					Fire_Period[Year-1] += 1 
					print("New Fire period: ", Fire_Period[Year-1],"\n")
									
					if WeatherOpt != 'constant' and Fire_Period[Year-1] * FirePeriodLen / MinutesPerWP > weatherperiod:
						weatherperiod +=1
						DF = Weather_Obj.update_Weather_FBP(DF, WeatherOpt, weatherperiod)
				
				if RepeatFire == True and len(Global_Message_Aux) > 0:
					print("Messages have been sent, next step",
						  "Current Fire period: ",Fire_Period[Year-1])
					RepeatFire = False 
								
				# Checking if the list is empty and no repeat flag, then if it is empty, end of the actual fire dynamic period, next year
				if MessagesSent==False and RepeatFire==False:
					if verbose == True:
						print( "\n","No messages during the fire period, end of year",Year)
						   
					#Next year, reset weeks and update burnt cells from burning cells
					Year+=1
					week_number=1
					
					# DLW jan 2018: ???
					# CP apr 2018: Burning cells are labeled as Burnt cells (no more messages), then
					# if save memory flag is True, we delete the cells objects saving memory...
					BurntCells_Set = BurntCells_Set.union(BurningCells_Set)
					BurningCells_Set = set()
					
					# if no savemem flag, keep the cell object and update status
					if SaveMem != True:
						for br in BurntCells_Set:
							Cells_Obj[br-1].Status=2
						
						# Otherwise, delete the inactive (burnt cells)
					if SaveMem == True:
						for c in BurntCells_Set:
							if (c-1) in Cells_Obj.keys():
								if verbose == True:
									print( "Deleting burnt cells from dictionary---------------------------------------------")
									print( "Deleted:",c)
								del Cells_Obj[c-1]
						
					break
															
				#Otherwise, go to next fire period and receiving messages loop
				if MessagesSent==True and RepeatFire==False:
					#Global list with messages (all messages)
					Global_Message_List = [val for sublist in SendMessageList for val in sublist]
					Global_Message_List.sort()
												
					if verbose == True:
						print( "Lists of messages per Cell: "+ str(SendMessageList))
						print( "Global Message Lists: " + str(Global_Message_List))
						print( "We have at least one message: "+str(MessagesSent))
							
						
					# Initialize cell (getting a message) if needed
					for bc in Global_Message_List:
						if (bc-1) not in Cells_Obj.keys() and bc not in BurntCells_Set:
							Cells_Obj[bc-1] = CellsFBP.Cells((bc), AreaCells, CoordCells[bc-1], AgeCells,
															 FTypeCells[bc-1], coef_ptr[FTypes2[str.lower(CellsGrid4[bc-1])]],
															 TerrainCells, VolCells, PerimeterCells, StatusCells[bc-1],
															 AdjCells[bc-1], Colors[bc-1], RealCells[bc-1], OutputGrid)
										
							'''
							New cell initialization
							'''
							Cells_Obj[bc-1].InitializeFireFields(CoordCells, AvailCells_Set)
							#Cells_Obj[bc-1].InitializeFireFields(CoordCells)
							'''
							'''
									
									
									
					#  Only active cells are being plotted (SaveMem mode)
					if plottrue == True and SaveMem == True:
						Plotter.forest_plotV3_FreeMem(Cells_Obj,SendMessageList,plotnumber,
													  Fire_Period[Year-1],Year,False,Rows,Cols,
													  PlotPath,CoordCells,BurntCells_Set,Sim)
						plotnumber+=1
									
					# Otherwise, plot everything
					if plottrue == True and SaveMem != True:
						Plotter.forest_plotV3(Cells_Obj,SendMessageList,plotnumber,
											  Fire_Period[Year-1],Year,False,Rows,Cols,
											  PlotPath,Sim)
						plotnumber+=1
						
					if OutputGrid:
						for c in Cells_Obj.keys():
							Cells_Obj[c].got_burnt_from_mem(Fire_Period[Year-1],
															SendMessageList,Year-1,verbose)

				### End Step 2
					# Releasing Memory 
					del SendMessageList

					###Step 3: Receiving messages
					#Printing information
					if verbose == True:
						print( "\n"+"---------------------- Step 3: Receiving and processing messages from Ignition ----------------------")
						print( "Size of the world: ",size)

					# Initializing Parameters
					BurntList = []
					NMessages = [0]*NCells
					GotMsg_Set = set(Global_Message_List)

					# Check which cells got messages and how many of them
					if verbose == True:
						print( "Cells receiving messages: "+str(GotMsg_Set))

					for i in range(1,NCells+1):
						NMessages[i-1]=Global_Message_List.count(i)
						if verbose == True and NMessages[i-1] != 0:
							print( "Cell "+str(i)+ " receives "+str(NMessages[i-1]) + " messages" )

					# Releasing Memory 
					del Global_Message_List

					if verbose == True:
						print( "\n" + "Cells status")
					for bc in GotMsg_Set:
						if bc not in BurntCells_Set:
							if (bc-1) not in Cells_Obj.keys():
								Cells_Obj[bc-1] = CellsFBP.Cells((bc),AreaCells,CoordCells[bc-1],AgeCells,
																 FTypeCells[bc-1],coef_ptr[FTypes2[str.lower(CellsGrid4[bc-1])]],
																 TerrainCells,VolCells,PerimeterCells,StatusCells[bc-1],
																 AdjCells[bc-1],Colors[bc-1],RealCells[bc-1])

								'''
								New cell initialization
								'''
								print( "-------- Initializing new cell ", bc, "--------")
								#Cells_Obj[bc-1].InitializeFireFields(CoordCells)
								Cells_Obj[bc-1].InitializeFireFields(CoordCells, AvailCells_Set)
								print( "-----------------------------------------------")
								'''
								'''

							if Cells_Obj[bc-1].FType != 0:
								Check_Burnt = Cells_Obj[bc-1].get_burned(Fire_Period[Year-1],NMessages[bc-1],Year,verbose,DF,coef_ptr,args.ROS_Threshold)
							else:
								Check_Burnt  = False
							if verbose == True:
								print( "Cell "+ str(Cells_Obj[bc-1].ID) +" got burnt: " + str(Check_Burnt))
							if Check_Burnt == True:
								BurntList.append(Cells_Obj[bc-1].ID)

					if verbose == True:
						print( "\nResults")
						print( "newly Burnt (and/or burning) List: " +str(BurntList))


				# Update cells status (burnt or not burnt), Update AvailCells and BurntCells sets
				# changed by dlw Jan 2018
				# tbd feb 2018: clean up the set assignments
				Aux_set = set(BurntList) # newly burning
				BurntCells_Set = BurntCells_Set.union(Aux_set)
				BurntCells_Set = BurntCells_Set.union(set(BurnedOutList))
				BurningCells_Set = BurningCells_Set.union(Aux_set)
				AvailCells_Set = AvailCells_Set.difference(Aux_set)

				# Releasing Memory 
				del Aux_set

				if verbose == True:
					print( "Available cells: " + str(AvailCells_Set))
					print( "Non Burnable Cells: " +str(NonBurnableCells_Set))
					print( "Burning cells: " + str(BurningCells_Set))
					print( "Harvest cells: " + str(HarvestCells_Set))
					print( "Burnt and Burning cells: " + str(BurntCells_Set))

				for t in range(Fire_Period[Year-1],Max_Fire_Periods+1):
					BurntP[t-1]=list(BurningCells_Set)

				#Plot
				if plottrue == True:
					if SaveMem == True:
						Plotter.forest_plotV3_FreeMem(Cells_Obj,emptylist,plotnumber,
														  Fire_Period[Year-1],Year,False,Rows,Cols,
														  PlotPath,CoordCells,BurntCells_Set,Sim)
						plotnumber+=1
					
					if SaveMem != True:
						Plotter.forest_plotV3(Cells_Obj,emptylist,plotnumber,Fire_Period[Year-1],
												  Year,False,Rows,Cols,PlotPath,Sim)
						plotnumber+=1

				# Next Period: t=t+1. Update Weather
				Fire_Period[Year-1]+=1
				if WeatherOpt != 'constant' and Fire_Period[Year-1] * FirePeriodLen / MinutesPerWP > weatherperiod:
					weatherperiod +=1
					DF = Weather_Obj.update_Weather_FBP(DF, WeatherOpt, weatherperiod=weatherperiod)
					if verbose:
						print( "Just did weather update")

				if verbose == True:
					Weather_Obj.print_info(weatherperiod)
					if WeatherOpt == 'constant':
						print( "\n(NOTE: current weather is not used for ROS with constant option)")
						print( "\nCurrent week: "+str(week_number))
						print( "Current Fire Period: "+str(Fire_Period[Year-1]))

				#Equivalence between fire periods and weeks (hours for the moment)
				# 168 hours per week
				if Fire_Period[Year-1] >= 168 * 60 / FirePeriodLen: 
					week_number+=1

				### End Step 3

                ### End of fire dynamic
            ### End of if ignition
        ###Year's loop (end)    


        ### Step 4: Results
        #End of the code for one sim, output files with statistics and plots
        for br in BurntCells_Set:
            if (br-1) in Cells_Obj.keys():
                Cells_Obj[br-1].Status=2
        for bn in BurningCells_Set:
            if (bn-1) in Cells_Obj.keys():    
                Cells_Obj[bn-1].Status=2

        if nooutput == False:
            print ("\n"+"----------------------------- Results -----------------------------")
            # General results
            print( "--------------------------- Solution without Heuristic --------------------------")
            print( "Total Available Cells:    ", len(AvailCells_Set),"- % of the Forest: ", np.round(len(AvailCells_Set)/NCells*100.0,3),"%")
            print( "Total Burnt Cells:        ", len(BurntCells_Set),"- % of the Forest: ", np.round(len(BurntCells_Set)/NCells*100.0,3),"%")
            print( "Total Non-Burnable Cells: ", len(NonBurnableCells_Set),"- % of the Forest: ", np.round(len(NonBurnableCells_Set)/NCells*100.0,3),"%")

            special_report(args, Sim, AreaCells, NCells, Fire_Period, AvailCells_Set,
                           BurntCells_Set, NonBurnableCells_Set)

            # Statistics: Cells' status, Fire periods, start, end.
            if SaveMem != True:
                print( "\n"+"Cells status")
                for i in range(0,NCells):
                    if i in Cells_Obj.keys():
                        if Cells_Obj[i].get_Status() == "Burnt":
                            print( "Cell "+str(i+1)+" status: "+ str(Cells_Obj[i].get_Status()) +", Year: "+str(Cells_Obj[i].FireStartsSeason) + ", Fire starts (fire period): "+ str(Cells_Obj[i].Firestarts))


            if SaveMem == True:
                for br in BurntCells_Set:
                    print( "Cell "+str(br)+" status: Burnt" )
                if verbose == True:
                    for av in AvailCells_Set:
                        print( "Cell "+str(av)+" status: Available" )

            # Total simulation time
            Final_Time = time.clock()
            print( "\nFinal Time:",Final_Time)
            print( "Total simulation time: ",np.round(Final_Time-Initial_Time,2)," [s]")

        # Plot 
        if plottrue == True:
            if SaveMem == True:
                Plotter.forest_plotV3_FreeMem(Cells_Obj,emptylist,plotnumber,1,
                                              Year,False,Rows,Cols,PlotPath,CoordCells,
                                              BurntCells_Set,Sim)
                plotnumber+=1

            if SaveMem != True:
                Plotter.forest_plotV3(Cells_Obj,emptylist,plotnumber,1,Year,
                                      False,Rows,Cols,PlotPath,Sim)
                plotnumber+=1

            # Check combine flag
            combine = args.input_combine
            if combine == True:
                for i in range(1,plotnumber):
                    Plotter.Mix(OutFolder,i,Sim)


        #Scenarios and Excel
        #Excel
        exceltrue = args.input_excel
        heuristic = args.input_heur

        if exceltrue == True and SaveMem != True:

            # Call the function
            Output_Grid.ExcelOutput(Cells_Obj,Sim,OutFolder,heuristic,TotalSims,TotalYears,
                                    NCells,BurntCells_Set,AvailCells_Set)

            print( "Excel file has been successfully created ")

        #Scenarios
        # If we are not in the save-memory version, check all the initialized cells and print scenarios.dat
        if SaveMem != True and scenarios==True:
            for cell in BurntCells_Set:
                Cells_Obj[cell-1].FS_definition()

            #Printing forest's data to a txt file
            Output_Grid.ScenarioOutput(TotalYears,Cells_Obj,NCells,
                                       BurntCells_Set,OutFolder,Sim,spotting,verbose)

            if nooutput==False:
                print( "Scenarios have been successfully created ")

            if SaveMem == True and scenarios == True and nooutput == False:
                print( "\nScenarios cannot be output in Save Memory mode")

            if SaveMem == True and exceltrue == True and nooutput == False:
                print( "\nStatistics cannot be output in Save Memory mode")

        # TBD, DLW? what do you do about weather in a new year?

        # Forest Grid
        if OutputGrid == True:
            Output_Grid.OutputGrid(OutFolder,Rows,Cols,BurntCells_Set,
                                   Sim,spotting, verbose, AreaCells)
            if nooutput ==False:
                print( "Forest Grid has been created")

        Sim+=1
    ### end of sim loop
    if nooutput == True:
        Final_Time = time.clock()
        print("\nFinal Time:",Final_Time)
