
################################################################################################################
#
#               FireSimulator serial and parallel version 1.0 - April 2018 
#               Author: Cristobal Pais
#               example: mpiexec -n X python Path\Simulator1Beta.py  where X is the number of parallel processes
#
################################################################################################################

# Importations
import numpy as np

# Lightning Class: Fire ignition class, determines the week when fire starts, based on a poisson strike	
class Lightning:
    PoissonRate = None

    '''
    Inputs:
    period       int
    '''
    # Returns true if a fire starts in a particular week of the fire season
    def Lambda_Simple_Test(self, period):
        Selected_Week = np.random.randint(1, 13)
        return Selected_Week

    
    '''
    Inputs:
    period       int
    verbose      string
    '''
    # Returns the probability of having a fire in a week of the fire season
    def Lambda_NH(self, period, verbose):
        Fire_Rate = 0.5
        AlfaFact = 0.1
        Poisson_Mean = (Fire_Rate / 2.0) * (2.0 + AlfaFact * ((period ** 2) - (period-1) ** 2 - 2.0))
        ProbsNoFire = np.round(np.exp(-Poisson_Mean), 2)
        if verbose == True:
            print("Probs of not fire (week "+ str(period)+"): "+str(ProbsNoFire))
         
        if np.round(np.random.uniform(0,1), 2) > ProbsNoFire:
            return True
        else:
            return False
        
    '''
    Inputs:
    period       int
    verbose      string
    '''        
    # Returns the probability of having a fire in a week of the fire season
    def Lambda_H(self, period, verbose):
        Poisson_Mean = 0.5
        ProbsNoFire = np.round(np.exp(-Poisson_Mean), 2)
        if verbose == True:
            print("Probs of not fire (week " + str(period)+"): " + str(ProbsNoFire))
            
        if np.round(np.random.uniform(0,1), 2) > ProbsNoFire:
            return True
        else:
            return False

