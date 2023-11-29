# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 16:33:07 2019

@author: eus358
# commented sections of code reomoved by R. Barton 10.19
"""

import VBASim #another python file
import Basic_Classes #another python file
import numpy as np

# parameter settings
m = 0 # size of jackson network or the number of servers in stage 1 (as well as 2)
WarmUpN = 0.0 #predetermined warmup number of customers
RunLengthN = 0.0
replications = 0
NumStations = 0
NumServStation = 0

Arate = []
Srate = []
QueueCapacity = []
Pij = []

# Statistics
WaitTime = []
Utilization = []

def init_run_params(input_m, input_WarmUpN, input_RunLengthN, input_replications):
    global m, WarmUpN, RunLengthN, replications, NumStations, NumServStation, Arate, Srate, QueueCapacity, Pij
    m = input_m
    WarmUpN = input_WarmUpN 
    RunLengthN = input_RunLengthN 
    replications = input_replications
    NumStations =  2*m+3 # total number of stations including the last 2*m+3rd station as exit
    NumServStation = NumStations-1 # service providing stations
    Arate = [1.0]*m
    Srate = [1.0]*m
    QueueCapacity = [0]*NumServStation 
    Pij = [[0.0 for x in range(NumStations)] for y in range(NumStations)] 
    

Clock = 0.0
#ZRNG = RNG.InitializeRNSeed()
Calendar = Basic_Classes.EventCalendar()


def init_input_params(input_Arate,input_Srate, input_QueueCapacity, input_Pij):
    global Arate, Srate, QueueCapacity, Pij, WaitTime, Utilization
    Arate = input_Arate
    Srate = input_Srate
    QueueCapacity = input_QueueCapacity
    Pij = input_Pij
    WaitTime = np.array([0.0]*replications,dtype='f')
    Utilization = [np.array([0.0]*NumServStation,dtype='f')  for y in range(replications)] 
    global StationQs, StationR, TheQueues, TheResources
    StationQs = [0]*NumServStation
    StationR = [0]*NumServStation
    for i in range(NumServStation):
        StationQs[i] = Basic_Classes.FIFOQueue()
        StationR[i] = Basic_Classes.Resource()
        StationR[i].SetUnits(1)
    
    for i in range(NumServStation):
        TheQueues.append(StationQs[i])
        TheResources.append(StationR[i])
    

# initialization
StationQs = [] 
StationR = []

TheCTStats = []
TheDTStats = []
TheQueues = []
TheResources = []

AllStationQs = []
AllStationService = []
AllStationR = []
AllStationQsNum = []
AllClock = []
    

route_i_list = []
    

def Route(i):
    cum_probs = np.cumsum(Pij[i,:])
    return np.where(cum_probs>np.random.uniform(0.0,1.0))[0][0]

def Arrival(j): # j represents the station number 
    if  j < m:
        VBASim.SchedulePlus(Calendar,"Arrival",np.random.exponential(1.0/Arate[j]),j,Clock)
    Customer = Basic_Classes.Entity(Clock)
    Customer.StationID = j    
    # following ensures if queue capacity is full, abandon the customer
    if StationQs[j].NumQueue() < QueueCapacity[j]:
        StationQs[j].Add(Customer,Clock)    
        if StationR[j].Busy < StationR[j].NumberOfUnits:
            StationR[j].Seize(1,Clock)
            VBASim.SchedulePlus(Calendar,"StartService",0,j,Clock)

def InternalArrival(Customer):
    # if queue capacity is full, abandon the customer
    if StationQs[Customer.StationID].NumQueue() < QueueCapacity[Customer.StationID]:
        StationQs[Customer.StationID].Add(Customer,Clock)    
        if StationR[Customer.StationID].Busy < StationR[Customer.StationID].NumberOfUnits:
            StationR[Customer.StationID].Seize(1,Clock)
            VBASim.SchedulePlus(Calendar,"StartService",0,Customer.StationID,Clock)

def StartService(j):
    global StationService
    Customer = StationQs[j].ThisQueue[0]
    ServiceTTime = np.random.exponential(1.0/Srate[j])
    Customer.CumulativeServiceTime=Customer.CumulativeServiceTime+ServiceTTime #variable accumulating servicetime for a particular customer at all stations
    VBASim.SchedulePlus(Calendar,"EndService",ServiceTTime,j,Clock)

def EndService(j):
    global StationWait # cumulative waiting time spent by all customers
    global NumCust # count for number of customers succesfully exiting the system
    global MeanWait # mean waiting time across all customers
    global done # variable used to triger the reset preocess after warmup
    Customer = StationQs[j].Remove(Clock)
    Customer.StationID = Route(j)
    if Customer.StationID == NumServStation:    # if the customer is exiting the system
        Customer.Wait=Clock - Customer.CreateTime- Customer.CumulativeServiceTime #Total time spent by the customer waiting in some or the other queue
        #print(Clock)
        #print(Customer.CreateTime)
        #print(Customer.CumulativeServiceTime)
        StationWait = StationWait+Customer.Wait 
        NumCust = NumCust+1  
        MeanWait = StationWait/NumCust        
    if StationQs[j].NumQueue() > 0: 
        VBASim.SchedulePlus(Calendar,"StartService",0,j,Clock)
    else:  
        StationR[j].Free(1,Clock)    
    if Customer.StationID < NumServStation:
        VBASim.SchedulePlus(Calendar,"InternalArrival",0,Customer,Clock)
    if done < 1: #checks if the warmup is reached
        if NumCust == WarmUpN:
            VBASim.Schedule(Calendar,"ClearIt",0,Clock)
            NumCust=0
            done = 1
            MeanWait=0.0
            StationWait = 0.0
    else:
        if NumCust==RunLengthN: #checks if runlength is reached
            VBASim.Schedule(Calendar,"EndSimulation",0,Clock)

Clock = 0.0
MeanWait=0.0
StationWait = 0.0
NumCust=0
done = 0

## Main code that runs simulations
def run_replications():
    global StationWait, NumCust, MeanWait, done, WaitTime, Clock, replications
    WaitTime = np.array([0.0]*replications,dtype='f')
    for reps in range(0,replications,1):   
        Clock = 0.0
        MeanWait=0.0
        StationWait = 0.0
        NumCust=0
        done = 0
        VBASim.VBASimInit(Calendar,TheQueues,TheCTStats,TheDTStats,TheResources,Clock)
        for i in range(m):
            VBASim.SchedulePlus(Calendar, "Arrival", np.random.exponential(1.0/Arate[i]),i, Clock)
        #VBASim.Schedule(Calendar,"EndSimulation",RunLength,Clock)
        #VBASim.Schedule(Calendar,"ClearIt",WarmUp,Clock)
        NextEvent = Calendar.Remove()
        Clock = NextEvent.EventTime
        if NextEvent.EventType == "Arrival":
            Arrival(NextEvent.WhichObject)
        elif NextEvent.EventType == "InternalArrival":
            InternalArrival(NextEvent.WhichObject) 
        elif NextEvent.EventType == "StartService":
            StartService(NextEvent.WhichObject)
        elif NextEvent.EventType == "EndService":
            EndService(NextEvent.WhichObject)
        elif NextEvent.EventType == "ClearIt":
           VBASim.ClearStats(TheCTStats,TheDTStats,Clock)
        while NextEvent.EventType !="EndSimulation":
            NextEvent = Calendar.Remove()
            Clock = NextEvent.EventTime
            #print(Clock)
            #print(NumCust)
            if NextEvent.EventType == "Arrival":
                Arrival(NextEvent.WhichObject)
            elif NextEvent.EventType == "InternalArrival":
                InternalArrival(NextEvent.WhichObject) 
            elif NextEvent.EventType == "StartService":
                StartService(NextEvent.WhichObject)
            elif NextEvent.EventType == "EndService":
                EndService(NextEvent.WhichObject) 
            elif NextEvent.EventType == "ClearIt":
               VBASim.ClearStats(TheCTStats,TheDTStats,Clock)
        WaitTime[reps]= MeanWait
        for r in range(0, NumServStation):
            Utilization[reps][r] = StationR[r].Mean(Clock);
