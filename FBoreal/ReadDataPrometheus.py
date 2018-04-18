
######################################################################################################################################
#
#                FireSimulator FBP serial and parallel version 1.0b - April 2018
#                Authors: Cristobal Pais, David L. Woodruff
#                example: mpiexec -n X python Path\Simulator1Beta.py  where X is the number of parallel processes
#
######################################################################################################################################

import json

#=============
# some stuff to make json work with python 2.7 taken from the web
def json_load_byteified(file_handle):
    return _byteify(
        json.load(file_handle, object_hook=_byteify),
        ignore_dicts=True
    )

def json_loads_byteified(json_text):
    return _byteify(
        json.loads(json_text, object_hook=_byteify),
        ignore_dicts=True
    )

def _byteify(data, ignore_dicts = False):
    # if this is a unicode string, return its string representation
    if isinstance(data, unicode):
        return data.encode('utf-8')
    # if this is a list of values, return list of byteified values
    if isinstance(data, list):
        return [ _byteify(item, ignore_dicts=True) for item in data ]
    # if this is a dictionary, return dictionary of byteified keys and values
    # but only if we haven't already byteified it
    if isinstance(data, dict) and not ignore_dicts:
        return {
            _byteify(key, ignore_dicts=True): _byteify(value, ignore_dicts=True)
            for key, value in data.iteritems()
        }
    # if it's anything else, return it in its original form
    return data
#================

'''
Inputs:
filename        str
'''
# Reads fbp_lookup_table.csv and creates dictionaries for the fuel types and cells' colors
def Dictionary(filename):
    aux = 1
    file = open(filename, "r") 
    row = {}
    colors = {} 
    all = {}
    
    # Read file and save colors and ftypes dictionaries
    for line in file: 
        if aux > 1:
            aux +=1
            line = line.replace("-","")
            line = line.replace("\n","")
            line = line.replace("No","NF")
            line = line.split(",")
            
            row[line[0]] = line[3][0:2]
            colors[line[0]] = (float(line[4])/ 255.0, 
                                    float(line[5]) / 255.0,
                                    float(line[6]) / 255.0,
                                    1.0)
            all[line[0]] = line
    
        if aux == 1:
            aux +=1
            
    return row, colors
    

'''
Inputs:
filename        str
Dictionary      Dictionary
'''   
# Reads the ASCII file with the forest grid structure and returns an array with all the cells and grid dimensions nxm
# Modified Feb 2018 by DLW to read the forest params (e.g. cell size) as well
def ForestGrid(filename, Dictionary):
    AdjCells = []
    North = "N"
    South = "S"
    East = "E"
    West = "W"
    NorthEast = "NE"
    NorthWest = "NW"
    SouthEast = "SE"
    SouthWest = "SW"
    
    CoordCells = []
    
    with open(filename, "r") as f:
        filelines = f.readlines()

    line = filelines[4].replace("\n","")
    parts = line.split()
    if parts[0] != "cellsize":
        print ("line=",line)
        raise RuntimeError("Expected cellsize on line 5 of "+ filename)
    cellsize = float(parts[1])
    
    cells = 0
    row = 1
    trows = 0 
    tcols = 0
    gridcell1 = []
    gridcell2 = []
    gridcell3 = []
    gridcell4 = []
    grid = []
    grid2 = []
    
    # Read the ASCII file with the grid structure
    for row in range(6, len(filelines)):
        line = filelines[row]
        line = line.replace("\n","")
        line = line.split(" ")

        for c in range(0,len(line)-1):
            if line[c] not in Dictionary.keys():
                gridcell1.append("NData")
                gridcell2.append("NData")
                gridcell3.append("NData")
                gridcell4.append("NData")
            else:
                gridcell1.append(line[c])
                gridcell2.append(Dictionary[line[c]])
                gridcell3.append(int(line[c]))
                gridcell4.append(Dictionary[line[c]])
            tcols = max(tcols,len(line))

        grid.append(gridcell1)
        grid2.append(gridcell2)
        gridcell1 = []
        gridcell2 = []
            
    # Adjacent list
    n = 1
    for r in range(0,len(grid)):
        for c in range(0,tcols-1):
            CoordCells.append([c,len(grid)-r-1])
            
            if len(grid) >1:
                
                if r == 0:
                    if c == 0:
                        AdjCells.append({North:None,NorthEast:None,NorthWest:None, 
                                         South:[n+tcols-1], SouthEast:[n+tcols], 
                                         SouthWest:None, East:[n+1],West:None})
                        n+=1
                    if c == tcols-2:
                        AdjCells.append({North:None,NorthEast:None,NorthWest:None,
                                         South:[n+tcols-1],SouthEast:None,SouthWest:[n+tcols-2], 
                                         East:None, West:[n-1]})
                        n+=1
                    if c>0 and c<tcols-2:    
                        AdjCells.append({North:None,NorthEast:None,NorthWest:None,
                                         South:[n+tcols-1],SouthEast:[n+tcols], 
                                         SouthWest:[n+tcols-2], East:[n+1],West:[n-1]})
                        n+=1
                
                if r > 0 and r < len(grid)-1:
                    if c == 0:
                        AdjCells.append({North:[n-tcols+1], NorthEast:[n-tcols+2],NorthWest:None,
                                         South:[n+tcols-1], SouthEast:[n+tcols], SouthWest:None,
                                         East:[n+1], West:None})
                        n+=1
                    if c == tcols-2:
                        AdjCells.append({North:[n-tcols+1], NorthEast:None, NorthWest:[n-tcols],
                                         South:[n+tcols-1], SouthEast:None, SouthWest:[n+tcols-2],
                                         East:None, West:[n-1]})
                        n+=1
                    if c>0 and c<tcols-2:    
                        AdjCells.append({North:[n-tcols+1], NorthEast:[n-tcols+2], NorthWest:[n-tcols],
                                         South:[n+tcols-1], SouthEast:[n+tcols], SouthWest:[n+tcols-2],
                                         East:[n+1], West:[n-1]})
                        n+=1        
                
                if r == len(grid)-1:
                    if c == 0:
                        AdjCells.append({North:[n-tcols+1], NorthEast:[n-tcols+2], NorthWest:None,
                                         South:None, SouthEast:None, SouthWest:None,
                                         East:[n+1], West:None})
                        n+=1
                        
                    if c == tcols-2:
                        AdjCells.append({North:[n-tcols+1], NorthEast:None, NorthWest:[n-tcols],
                                         South:None, SouthEast:None, SouthWest:None,
                                         East:None, West:[n-1]})
                        n+=1
                        
                    if c>0 and c<tcols-2:    
                        AdjCells.append({North:[n-tcols+1], NorthEast:[n-tcols+2], NorthWest:[n-tcols],
                                         South:None, SouthEast:None,SouthWest:None,
                                         East:[n+1], West:[n-1]})
                        n+=1
            
            if len(grid)==1:
                if c == 0:
                    AdjCells.append({North:None, NorthEast:None, NorthWest:None,
                                     South:None, SouthEast:None, SouthWest:None,
                                     East:[n+1], West:None})
                    n+=1
                if c == tcols-2:
                    AdjCells.append({North:None, NorthEast:None, NorthWest:None,
                                     South:None, SouthEast:None, SouthWest:None,
                                     East:None,West:[n-1]})
                    n+=1
                if c>0 and c<tcols-2:    
                    AdjCells.append({North:None, NorthEast:None, NorthWest:None,
                                     South:None, SouthEast:None, SouthWest:None,
                                     East:[n+1], West:[n-1]})
                    n+=1
    
    
    return gridcell3, gridcell4, len(grid), tcols-1, AdjCells, CoordCells, cellsize

'''
Inputs:
filename        str
'''
# Reads fbp_lookup_table.csv and creates dictionaries for the fuel types and cells' colors (Pandas' version - TBD)    
def Dictionary_PD(filename):
    pd.set_option('display.mpl_style','default')
    plt.rcParams['figure.figsize'] = (15, 5)
    fbplookuptableDF = pd.read_csv(filename,sep=",",index_col="grid_value")
    fbplookuptableDF["export_value"]
    fbplookuptableDF
    

'''
Inputs:
filename        str
'''    
# Reads IgnitionPoints.csv file and creates an array with them 
def IgnitionPoints(filename):
    #Ignitions is a dictionary with years = keys and ncell = values
    aux = 1
    file = open(filename, "r") 
    ignitions = {}
    for line in file:
        if aux > 1:
            line = line.replace("\n","")
            line = line.split(",")
            ignitions[int(line[0])] = int(line[1])
        if aux==1:
            aux+=1    
    return ignitions        
    

'''
Inputs
filename        str
nooutput        boolean
'''
# Reads spotting parameters json file and returns dictionary
def ReadSpotting(filename,nooutput):
    # As of Jan 2018 the dictionary should have "SPOT" and "SPTANGLE"
    with open(filename, 'r') as f:
        SpottingParams = json_load_byteified(f)
    ### Thresholds["SPOT"] = 10
    ### Thresholds["SPTANGLE"] = 30
    
    if nooutput == False:
        print("---- Spotting Parameters ----")
        for i in SpottingParams:
            print(i,":",SpottingParams[i])
        print("......................")    
    return SpottingParams


