#!/usr/bin/env python3

#Setup
import itertools
import operator
import copy
import numpy as np
import pandas as pd
import math
#from matplotlib import pyplot as plt
#import matplotlib as mpl

#from statsmodels.nonparametric.api import KernelReg
#import pickle
import random
#import matplotlib.cm as cm
import os


#function to check if adjacent cells are empty

def adjacent_cell(lattice,i,j):
    free_cells = []
    change = [1,-1]
    for c in change:
        checki = lattice[i+c,j]
        if checki == None:
            free_cells.append([i+c,j])
        checkj = lattice[i, j+c]
        if checkj == None:
            free_cells.append([i,j+c])
        for h in change:
            checkd = lattice[i+h,j+c]
            if checkd == None:
                free_cells.append([i+h,j+c])

    return free_cells

def cancer_sim(founder_cells = 1,deltat=(1/float(24)),CSC = True,max_cells = 1000,
               proliferationCSC = 1, proliferationCC = 2, motility = 15, pmax = 10,
               mortality = 0.1, mutation = 0.1, selection = False,ps = 0.05,
              free_diffusion = False, JC = False, s = 0.1, driverMutCell = 100):

    #create classes
    class Cell:
        def __init__(self):

            self.parent_index = 0 #index of cell it was derived from (for df conversion)

            self.locx = None #x location in lattice
            self.locy = None #y location in lattice
            self.birthdate = None
            self.deathdate = None
            self.pmax = pmax #proliferation potential
            self.cellnum = None #unique cell number to identify
            self.alpha = mortality #probability of sponaneous death
            self.mutation_rate = mutation #probability of mutation per cell division

            self.proliferation_rate = proliferationCC #average cell divisions per time unit

            self.mutations = [] #store index of mutations
            self.sequence = [] #store simulated nucleotide sequence under JC model
            self.ps = 0
            self.cells_free = 8 #add to keep track of spatial constraint at time of cell birth
            self.driver = 0


#function to convert results into dictionary, when can be converted to dataframe
        def as_dict(self):


            return {'index': self.cellnum,
                    'parent_index': self.parent_index,
                    'birthdate': self.birthdate,
                    'locx': self.locx,
                    'locy': self.locy,
                    'deathdate': self.deathdate,
                    'mutations': self.mutations,
                    'proliferation_rate': self.proliferation_rate,
                    'alpha': self.alpha,
                    'mutation_rate': self.mutation_rate,
                    'cells_free':self.cells_free,
                    'sequence': self.sequence,
                    'driver': self.driver}


    class StemCell(Cell): #define stem cell class
        def __init__(self):
            Cell.__init__(self)
            self.pmax = float('inf') #infinite proliferation potential
            self.alpha = 0 #immortal
            self.ps = ps #probability of symmetric division
            self.proliferation_rate = proliferationCSC


    #create lattice array
    N = max_cells # starting demensions of lattice
    lattice = np.empty( (N,N), dtype=object)

    #parameters
    proliferation_rate_CC = proliferationCC #proliferation rate of clonal cell
    proliferation_rate_CSC = proliferationCSC #proliferation rate of cancer stem cell
    pmax_CC = pmax

    alpha_CC = mortality
    cur_cellnum = 1
    cur_mutnum = 1

    mutation_rate = mutation
    time = 0

    #initialize list to keep track of cells (all time and alive)
    alive_cells = []
    cells = []

    nucleotides = ["A", "G", "C", "T"]

#to place founder cell in center of lattice
    i = int(N/2)
    modi = [0,0,1,1]
    modj = [0,1,0,1]

    if CSC == True:
        for x in range(founder_cells):
            Fcell = StemCell() #initiate founder cells, stem-cells
            Fcell.ps = ps
            Fcell.proliferation_rate = proliferation_rate_CSC #average cell divisions per day
            Fcell.locx = i+modi[x] #x location in lattice
            Fcell.locy = i+modj[x] #y location in lattice
            Fcell.birthdate = 0
            Fcell.cellnum = cur_cellnum #unique cell number to identify
            cur_cellnum += 1
            Fcell.mutation_rate = mutation_rate #probability of mutation per cell division


            cells.append(Fcell)
            alive_cells.append(Fcell)

            lattice[Fcell.locx,Fcell.locy] = Fcell
    else:
        for x in range(founder_cells):
            Fcell = StemCell() #initiate founder cells, non-stem cells
            Fcell.proliferation_rate = proliferation_rate_CC
            Fcell.locx = i+modi[x] #x location in lattice
            Fcell.locy = i+modj[x] #y location in lattice
            Fcell.alpha = mortality
            Fcell.ps = 0
            Fcell.birthdate = 0
            Fcell.cellnum = cur_cellnum  #unique cell number to identify
            cur_cellnum += 1
            Fcell.mutation_rate = mutation_rate #probability of mutation per cell division


            cells.append(Fcell)
            alive_cells.append(Fcell)
            lattice[Fcell.locx,Fcell.locy] = Fcell


    #time parameters
    dt = deltat # time is equilavent to 1/24 of a day or 1 hour

    while len(alive_cells) < max_cells:
        if len(alive_cells) < 1:
            print('no cells alive')
            #start overr simulation
            cells = []
            alive_cells = []
            time = 0
            cur_cellnum = 1
            cur_mutnum = 1

            lattice = np.empty( (N,N), dtype=object)
            for x in range(founder_cells):
                Fcell = StemCell() #initiate founder cells, non-stem cells
                Fcell.proliferation_rate = proliferation_rate_CC
                Fcell.locx = i+modi[x] #x location in lattice
                Fcell.locy = i+modj[x] #y location in lattice
                Fcell.alpha = mortality
                Fcell.ps = 0
                Fcell.birthdate = 0
                Fcell.cellnum = cur_cellnum  #unique cell number to identify
                cur_cellnum += 1
                Fcell.mutation_rate = mutation_rate #probability of mutation per cell division per time unit


                cells.append(Fcell)
                alive_cells.append(Fcell)
                lattice[Fcell.locx,Fcell.locy] = Fcell


        time += dt
        cell_stack = random.sample(alive_cells,len(alive_cells)) #random order of cells

        for cell in cell_stack:
            alive = True
            r = random.uniform(0, 1) #decision for prolferation
            r2 = random.uniform(0, 1) #decision for cell death
            pd = cell.proliferation_rate * dt #probability of proliferation in time dt

            if free_diffusion: #for free diffusion model cell is naive to neighbors in movement
                free_cells = [[cell.locx + 1, cell.locy + 1],
                             [cell.locx - 1, cell.locy + 1],
                             [cell.locx + 1, cell.locy - 1],
                             [cell.locx - 1, cell.locy - 1]]
            else:
                free_cells = adjacent_cell(lattice,cell.locx,cell.locy)

            p_die = cell.alpha * dt #probability of cell death in time dt

            if r < p_die:
                alive = False
                lattice[cell.locx,cell.locy] = None
                alive_cells.remove(cell)

                cell.deathdate = time


            elif r2 < pd: # Does cell attempt to divide?
                if len(free_cells) > 0: #is there any space to divide?

                    if cell.pmax > 0: #is cell proliferation capacity exhausted?

                        cell.pmax -= 1
                        r = random.uniform(0,1)
                        new_cell = copy.deepcopy(cell)
                        new_cell_2 = copy.deepcopy(cell)

                        if r <= cell.ps: #does cell divide asymmetrically?
                            new_cell.pmax = pmax
                            new_cell.proliferation_rate = proliferation_rate_CC #average cell divisions per day
                            new_cell.alpha = mortality
                            new_cell.ps = 0

                        newloc = random.choice(free_cells)
                        new_cell.locx = newloc[0] #x location in lattice
                        new_cell.locy = newloc[1] #y location in lattice
                        new_cell_2.locx = cell.locx #x location in lattice
                        new_cell_2.locy = cell.locy #y location in lattice
                        new_cell.birthdate = time
                        new_cell_2.birthdate = time
                        new_cell.cellnum = cur_cellnum #unique cell number to identify

                        #add selection, when driver mutation occurs is dependent on driverMutCell
                        if cur_cellnum == driverMutCell and selection:

                            new_cell.proliferation_rate = (1+s)* new_cell.proliferation_rate
                            new_cell.driver = 1

                        cur_cellnum += 1
                        new_cell_2.cellnum = cur_cellnum
                        new_cell.parent_index = cell.cellnum
                        new_cell_2.parent_index = cell.cellnum


                        if cur_cellnum == driverMutCell and selection:

                            new_cell_2.proliferation_rate = (1+s)* new_cell_2.proliferation_rate
                            new_cell_2.driver = 1

                        cur_cellnum += 1

                        lattice[new_cell.locx,new_cell.locy] = new_cell
                        lattice[cell.locx,cell.locy] = new_cell_2

                        if not free_diffusion:
                            new_cell.cells_free = len(adjacent_cell(lattice,new_cell.locx,new_cell.locy))
                            new_cell_2.cells_free = len(adjacent_cell(lattice,new_cell_2.locx,new_cell_2.locy))

                        r1 = random.uniform(0,1)
                        if r1 < cell.mutation_rate: #does cell gain a mutation?

                            if JC: #if Jukes-Cantor model

                                original_nucleotide = random.choice(["A", "G", "C", "T"])

                                #put original nucleotide in all sequences
                                for cell_0 in cells:
                                    cell_0.sequence.append(original_nucleotide)

                                #cell 2 has not been added yet
                                new_cell_2.sequence.append(original_nucleotide)

                                #replace mutation in this cells sequence

                                possible_nucs = ["A", "G", "C", "T"]
                                possible_nucs.remove(original_nucleotide)
                                new_nucleotide = random.choice(possible_nucs)
                                new_cell.sequence.append(new_nucleotide)

                            new_mut = cur_mutnum
                            cur_mutnum +=1
                            new_cell.mutations.append(new_mut)

                        #does second cell get a mutation
                        r2 = random.uniform(0,1)
                        if r2 < cell.mutation_rate: #does cell gain a mutation?

                            if JC: #if Jukes-Cantor model

                                original_nucleotide = random.choice(["A", "G", "C", "T"])

                                #put original nucleotide in all sequences
                                for cell_0 in cells:
                                    cell_0.sequence.append(original_nucleotide)

                                new_cell.sequence.append(original_nucleotide)

                                possible_nucs = ["A", "G", "C", "T"]
                                possible_nucs.remove(original_nucleotide)
                                new_nucleotide = random.choice(possible_nucs)
                                new_cell_2.sequence.append(new_nucleotide)

                            new_mut = cur_mutnum
                            cur_mutnum +=1
                            new_cell_2.mutations.append(new_mut)


                        alive_cells.append(new_cell)
                        alive_cells.append(new_cell_2)
                        cells.append(new_cell)
                        cells.append(new_cell_2)
                        alive_cells.remove(cell)

                        cell.deathdate = time


    for cell in alive_cells:

        cell.deathdate = time + dt #need to add one time step to differentiate cells that died in last gen


    return cells, alive_cells, lattice

    #Simulation with pushing
    ##Same as above except simulation creates new cell when pushed
    ##Also saves locations through time since they can change during a cells lifespan

def pushing_cancer_sim2(founder_cells = 1,deltat=(1/float(24)),CSC = True,max_cells = 1000,
                   proliferationCSC = 1, proliferationCC = 2, motility = 15, pmax = 10,
                   mortality = 0.1, mutation = 0.1, selection = False, JC = True):

        #create classes
        class Cell:
            def __init__(self):

                self.parent_index = 0 #index of cell it was derived from (for df conversion)

                self.locx = None #x location in lattice
                self.locy = None #y location in lattice
                self.birthdate = None
                self.deathdate = None
                self.pmax = pmax #proliferation potential
                self.cellnum = None #unique cell number to identify
                self.alpha = mortality #probability of sponaneous death
                self.mutation_rate = mutation #probability of mutation per cell division

                self.proliferation_rate = proliferationCC #average cell divisions per time unit

                self.mutations = [] #store index of mutations
                self.sequence = [] #store simulated nucleotide sequence under JC model
                #self.ps = 0
                self.cells_free = 8 #add to keep track of spatial constraint at time of cell birth


    #function to convert results into dictionary, when can be converted to dataframe
            def as_dict(self):


                return {'index': self.cellnum,
                        'parent_index': self.parent_index,
                        'birthdate': self.birthdate,
                        'locx': self.locx,
                        'locy': self.locy,
                        'deathdate': self.deathdate,
                        'mutations': self.mutations,
                        'proliferation_rate': self.proliferation_rate,
                        'alpha': self.alpha,
                        'mutation_rate': self.mutation_rate,
                        'cells_free':self.cells_free,
                        'sequence': self.sequence}

            def get_loc_dict(self):


                return {'index': self.cellnum,
                        'locx': self.locx,
                        'locy': self.locy}


        class StemCell(Cell): #define stem cell class
            def __init__(self):
                Cell.__init__(self)
                self.pmax = float('inf') #infinite proliferation potential
                self.alpha = 0 #immortal
                #self.ps = ps #probability of symmetric division
                self.proliferation_rate = proliferationCSC


        #create lattice array
        N = max_cells # starting demensions of lattice
        lattice = np.empty( (N,N), dtype=object)

        #parameters
        proliferation_rate_CC = proliferationCC #proliferation rate of clonal cell
        proliferation_rate_CSC = proliferationCSC #proliferation rate of cancer stem cell
        pmax_CC = pmax

        alpha_CC = mortality
        cur_cellnum = 1
        cur_mutnum = 1

        mutation_rate = mutation
        time = 0

        #initialize list to keep track of cells (all time and alive)
        alive_cells = []
        cells = []


        nucleotides = ["A", "G", "C", "T"]

        #create data frame for location time stamps
        #df_empty = pd.DataFrame({'A' : []})
        all_locations_df = pd.DataFrame(columns = ['index', 'locx', 'locy', 't'])


    #to place founder cell in center of lattice
        i = int(N/2)
        modi = [0,0,1,1]
        modj = [0,1,0,1]

        for x in range(founder_cells):
            Fcell = StemCell() #initiate founder cells, non-stem cells
            Fcell.proliferation_rate = proliferation_rate_CC
            Fcell.locx = i+modi[x] #x location in lattice
            Fcell.locy = i+modj[x] #y location in lattice
            Fcell.alpha = mortality

            Fcell.birthdate = 0
            Fcell.cellnum = cur_cellnum  #unique cell number to identify
            cur_cellnum += 1
            Fcell.mutation_rate = mutation_rate #probability of mutation per cell division


            cells.append(Fcell)
            alive_cells.append(Fcell)
            lattice[Fcell.locx,Fcell.locy] = Fcell


        #time parameters
        dt = deltat # time is equilavent to 1/24 of a day or 1 hour
    #     stop_time =stime#stop simulation after this many days


        #while time < stop_time:
        while len(alive_cells) < max_cells:

            if len(alive_cells) < 1:
                print('no cells alive')
                #start overr simulation
                cells = []
                alive_cells = []
                time = 0
                cur_cellnum = 1
                cur_mutnum = 1

                lattice = np.empty( (N,N), dtype=object)
                all_locations_df = pd.DataFrame(columns = ['index', 'locx', 'locy', 't'])

                for x in range(founder_cells):
                    Fcell = StemCell() #initiate founder cells
                    Fcell.proliferation_rate = proliferation_rate_CC
                    Fcell.locx = i+modi[x] #x location in lattice
                    Fcell.locy = i+modj[x] #y location in lattice
                    Fcell.alpha = mortality
                   # Fcell.ps = 0
                    Fcell.birthdate = 0
                    Fcell.cellnum = cur_cellnum  #unique cell number to identify
                    cur_cellnum += 1
                    Fcell.mutation_rate = mutation_rate #probability of mutation per cell division per time unit


                    cells.append(Fcell)
                    alive_cells.append(Fcell)
                    lattice[Fcell.locx,Fcell.locy] = Fcell


            time += dt
            cell_stack = random.sample(alive_cells,len(alive_cells)) #random order of cells



            for cell in cell_stack:

                #print(cell.cellnum)
                #print(cell in alive_cells)
                alive = True
                r = random.uniform(0, 1) #decision for prolferation
                r2 = random.uniform(0, 1) #decision for cell death
                p_d = cell.proliferation_rate * dt #probability of proliferation in time dt


                free_cells = adjacent_cell(lattice,cell.locx,cell.locy)



                p_die = cell.alpha * dt #probability of cell death in time dt

                if r < p_die:
                    alive = False
                    lattice[cell.locx,cell.locy] = None
                    alive_cells.remove(cell)

                    cell.deathdate = time


                elif r2 < p_d: # Does cell attempt to divide?

                    #print(cell in alive_cells)
                    cell.pmax -= 1
                    r = random.uniform(0,1)
                    new_cell = copy.deepcopy(cell)
                    new_cell_2 = copy.deepcopy(cell)

                    new_cell_2.locx = cell.locx #x location in lattice
                    new_cell_2.locy = cell.locy #y location in lattice

                    #alive_cells.remove(cell)
                    if len(free_cells) > 0:

                        newloc = random.choice(free_cells)
                        new_cell.locx = newloc[0] #x location in lattice
                        new_cell.locy = newloc[1] #y location in lattice
                        lattice[new_cell.locx, new_cell.locy] = new_cell

                    else:

                        #choose direction to push (x,y vector)


                        pushx = random.choice([-1,0,1])
                        pushy = random.choice([-1,0,1])

                        while pushx == 0 and pushy == 0: #redraw if no shift

                            pushx = random.choice([-1,0,1])
                            pushy = random.choice([-1,0,1])

                        assert pushx != 0 or pushy != 0, "Push direction is 0,0"

                        new_cell.locx = cell.locx + pushx #x location in lattice
                        new_cell.locy = cell.locy + pushy #y location in lattice

                        lattice[new_cell.locx, new_cell.locy] = new_cell

                        #then push all other alive cells
                        curr_overlap = [new_cell.locx, new_cell.locy]

                        curr_occupying_cell = new_cell


                        while curr_overlap:



                            overlapping_cell = next((x for x in alive_cells if ([x.locx,x.locy] == curr_overlap and x != curr_occupying_cell)), None)



                            if overlapping_cell:



                                curr_free_cells = adjacent_cell(lattice, overlapping_cell.locx, overlapping_cell.locy)


                                if len(curr_free_cells) > 0:


                                    push_newloc = random.choice(curr_free_cells)
                                    overlapping_cell.locx = push_newloc[0] #x location in lattice
                                    overlapping_cell.locy = push_newloc[1] #y location in lattice
                                    curr_overlap = None
                                    curr_occupying_cell = overlapping_cell


                                else:    #if there is no empty lattice spot then keep pushing in same direction


                                    overlapping_cell.locx = overlapping_cell.locx + pushx
                                    overlapping_cell.locy = overlapping_cell.locy + pushy
                                    curr_overlap = [overlapping_cell.locx, overlapping_cell.locy]


                                    curr_occupying_cell = overlapping_cell

                                #update lattice and occupying cell
                                lattice[overlapping_cell.locx, overlapping_cell.locy] = overlapping_cell

                            else:

                                curr_overlap = None

                    new_cell.move = False
                    new_cell_2.move = False

                    new_cell.birthdate = time
                    new_cell_2.birthdate = time
                    new_cell.cellnum = cur_cellnum #unique cell number to identify



                    cur_cellnum += 1

                    new_cell_2.cellnum = cur_cellnum
                    new_cell.parent_index = cell.cellnum
                    new_cell_2.parent_index = cell.cellnum


                    cur_cellnum += 1
                    lattice[cell.locx,cell.locy] = new_cell_2


                    new_cell.cells_free = len(adjacent_cell(lattice,new_cell.locx,new_cell.locy))
                    new_cell_2.cells_free = len(adjacent_cell(lattice,new_cell_2.locx,new_cell_2.locy))

                    r1 = random.uniform(0,1)
                    if r1 < cell.mutation_rate: #does cell gain a mutation?

                        if JC: #if Jukes-Cantor model

                            original_nucleotide = random.choice(["A", "G", "C", "T"])

                            #put original nucleotide in all sequences
                            for cell_0 in cells:
                                cell_0.sequence.append(original_nucleotide)

                            #cell 2 has not been added yet
                            new_cell_2.sequence.append(original_nucleotide)

                            #replace mutation in this cells sequence

                            possible_nucs = ["A", "G", "C", "T"]
                            possible_nucs.remove(original_nucleotide)
                            new_nucleotide = random.choice(possible_nucs)
                            new_cell.sequence.append(new_nucleotide)

                        new_mut = cur_mutnum
                        cur_mutnum +=1
                        new_cell.mutations.append(new_mut)

                    #does second cell get a mutation
                    r2 = random.uniform(0,1)
                    if r2 < cell.mutation_rate: #does cell gain a mutation?

                        if JC: #if Jukes-Cantor model

                            original_nucleotide = random.choice(["A", "G", "C", "T"])

                            #put original nucleotide in all sequences
                            for cell_0 in cells:
                                cell_0.sequence.append(original_nucleotide)

                            new_cell.sequence.append(original_nucleotide)

                            possible_nucs = ["A", "G", "C", "T"]
                            possible_nucs.remove(original_nucleotide)
                            new_nucleotide = random.choice(possible_nucs)
                            new_cell_2.sequence.append(new_nucleotide)

                        new_mut = cur_mutnum
                        cur_mutnum +=1
                        new_cell_2.mutations.append(new_mut)


                    alive_cells.append(new_cell)
                    alive_cells.append(new_cell_2)
                    cells.append(new_cell)
                    cells.append(new_cell_2)

                    alive_cells.remove(cell)

                    cell.deathdate = time

            #snapshot of tumor locations
            location_df = pd.DataFrame([x.get_loc_dict() for x in alive_cells])
            location_df['t'] = time
            #all_locations_df = all_locations_df.append(location_df, ignore_index = False)
            all_locations_df = pd.concat([all_locations_df, location_df], ignore_index = True)

        for cell in alive_cells:

            cell.deathdate = time + dt #need to add one time step to differentiate cells that died in last gen


        return cells, alive_cells, all_locations_df

##simulations for x-y plot of estimated diversification rates versus true
##this is modified by death rate, keeping proliferation rate constant
##
# random.seed(822156)
# n_cells = [10000]
# dr_list = np.arange(0, 0.43, 0.005).tolist()
# for N in n_cells:
#     for dr in dr_list:
#
#         #run simulated
#         cells_CC, alive_cells_CC, lattice_CC = cancer_sim(founder_cells = 1,max_cells = N, proliferationCSC = 0.5,
#         proliferationCC = 0.5, CSC = False, pmax = 10,mutation = 0.1,mortality = dr, JC = True)
#
#         #convert to dataframe
#         CC_df = pd.DataFrame([x.as_dict() for x in cells_CC])
#
#         file_name = f"cells_death_rate_validation_pop_{N}_dr_%0.3f.csv" %(dr)
#         CC_df.to_csv(file_name)


##multiple replicates of same death rates

# random.seed(31221)
# n_cells = [10000]
# dr_list = np.arange(0.0, 0.43, 0.025).tolist()
# i_list = np.arange(0,10,1).tolist()
# for N in n_cells:
#     for dr in dr_list:
#         for i in i_list:
#
#             #run simulated
#             cells_CC, alive_cells_CC, lattice_CC = cancer_sim(founder_cells = 1,max_cells = N, proliferationCSC = 0.5,
#             proliferationCC = 0.5, CSC = False, pmax = 10,mutation = 0.1,mortality = dr, JC = True)
#
#             #convert to dataframe
#             CC_df = pd.DataFrame([x.as_dict() for x in cells_CC])
#
#             file_name = f"cells_death_rate_validation_pop_{N}_dr_%0.2f_i_{i}.csv" %(dr)
#             CC_df.to_csv(file_name)


# #iterations for branching and clock rate signal analysis
# random.seed(2121)
#
# N=10000
#
# dr = 0.05
# i_list = np.arange(0, 100, 1).tolist()
#
# for itr in i_list:
#     #run simulated
#     cells_CC, alive_cells_CC, lattice_CC = cancer_sim(founder_cells = 1,max_cells = N, proliferationCSC = 0.5,
#     proliferationCC = 0.5, CSC = False, pmax = 10,mutation = 0.1, mortality = dr, JC = True)
#
#     #convert to dataframe
#     CC_df = pd.DataFrame([x.as_dict() for x in cells_CC])
#
#     file_name = f"cells_death_rate_validation_{N}_i_{itr}_dr_%0.2f.csv" %(dr)
#     CC_df.to_csv(file_name)


# #run multiple pushing simulations at dr = 0.1
# random.seed(9918)
#
# N=10000
#
# dr = 0.05
# i_list = np.arange(0, 10, 1).tolist()
#
# for itr in i_list:
#     #run simulated
#     cells_CC, alive_cells_CC, all_locations_df_CC = pushing_cancer_sim2(founder_cells = 1,max_cells = N, proliferationCSC = 0.5,
#     proliferationCC = 0.5, CSC = False, pmax = 10,mutation = 0.1, mortality = dr, JC = True)
#
#     #save locations
#     locs_file_name = f"cells_pushing_pop2_{N}_dr_%0.3f_i_{itr}_locs.csv" %(dr)
#     all_locations_df_CC.to_csv(locs_file_name)
#     #convert to dataframe
#     CC_df = pd.DataFrame([x.as_dict() for x in cells_CC])
#
#     file_name = f"cells_pushing_pop2_{N}_dr_%0.3f_i_{itr}.csv" %(dr)
#     CC_df.to_csv(file_name)
#
# #Run one iteration accross range of death rates
random.seed(8161)

N=10000

dr_list = np.arange(0, 0.43, 0.005).tolist()

for dr in dr_list:
    print(dr)
    #run simulated
    cells_CC, alive_cells_CC, all_locations_df_CC = pushing_cancer_sim2(founder_cells = 1,max_cells = N, proliferationCSC = 0.5,
    proliferationCC = 0.5, CSC = False, pmax = 10,mutation = 0.1, mortality = dr, JC = True)

    #save locations
    locs_file_name = f"cells_pushing_pop_{N}_dr_%0.3f_locs.csv" %(dr)
    all_locations_df_CC.to_csv(locs_file_name)
    #convert to dataframe

    #convert to dataframe
    CC_df = pd.DataFrame([x.as_dict() for x in cells_CC])

    file_name = f"cells_pushing_pop_{N}_dr_%0.3f.csv" %(dr)
    CC_df.to_csv(file_name)
#
# #run all pushing validations
# #interations
# #intentionally same seed as above
# random.seed(1211)
#
# N=10000
#
# dr_list = np.arange(0.025, 0.33, 0.0025).tolist()
# i_list = np.arange(0,10,1).tolist()
#
# for dr in dr_list:
#     for i in i_list:
#         #run simulated
#         cells_CC, alive_cells_CC, all_locations_df_CC = pushing_cancer_sim2(founder_cells = 1,max_cells = N, proliferationCSC = 0.5,
#         proliferationCC = 0.5, CSC = False, pmax = 10,mutation = 0.1, mortality = dr, JC = True)
#
#
#         #save locations
#         locs_file_name = f"cells_pushing_pop_{N}_dr_%0.3f_i_{i}_locs.csv" %(dr)
#         all_locations_df_CC.to_csv(locs_file_name)
#         #convert to dataframe
#         CC_df = pd.DataFrame([x.as_dict() for x in cells_CC])
#
#         file_name = f"cells_pushing_pop_{N}_dr_%0.3f_i_{i}.csv" %(dr)
#         CC_df.to_csv(file_name)
#
# #run all pushing validations at dr = 0.1
# random.seed(2121)
#
# N=10000
#
# dr = 0.05
# i_list = np.arange(0, 100, 1).tolist()
#
# for itr in i_list:
#     #run simulated
#     cells_CC, alive_cells_CC, lattice_CC = pushing_cancer_sim(founder_cells = 1,max_cells = N, proliferationCSC = 0.5,
#     proliferationCC = 0.5, CSC = False, pmax = 10,mutation = 0.1, mortality = dr, JC = True)
#
#     #convert to dataframe
#     CC_df = pd.DataFrame([x.as_dict() for x in cells_CC])
#
#     file_name = f"cells_pushing_pop_{N}_i_{itr}_dr_%0.2f.csv" %(dr)
#     CC_df.to_csv(file_name)
