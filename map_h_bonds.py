from ase.atom import Atom
from ase.atoms import Atoms
from ase.io import read, write
from ase import neighborlist

import numpy as np
from scipy import sparse



# We're looking for any bond from a H3PO4 H atom to an O atom in another H3PO4
def find_h_bonds(atoms):
    neighborList = neighborlist.NeighborList(h_cutoffs, self_interaction=False, bothways=True)
    neighborList.update(atoms)
    this_hbmatrix = np.zeros(shape=(59,59))
    nhbonds = 0
    for atom in atoms:
        if atom.index >= 288 and atom.symbol == 'H': #Look one way, make matrix symmetric
            indices, offsets = neighborList.get_neighbors(atom.index)
            for idx in indices:
                if atoms[idx].symbol == 'O':
                    #check different molelcule
                    if master_component_list[atom.index] != master_component_list[idx]:
                        nhbonds +=1
                        #print(master_component_list[atom.index],master_component_list[idx])
                        this_hbmatrix[master_component_list[atom.index],master_component_list[idx]] = 1
                        this_hbmatrix[master_component_list[idx],master_component_list[atom.index]] = 1
    return(this_hbmatrix,nhbonds)

# COF = atoms 1-288 (1-index)
# 57 * each H3PO4 successive 8 atoms
#fulltraj = read('geo_end_small.xyz',index=':')
fulltraj = read('geo_end.xyz',index='-3:-1') #testing

#xyz has no pbc:
pbc_info =  read('geo_end.cif')
for i in range(len(fulltraj)):
    fulltraj[i].pbc = pbc_info.pbc
    fulltraj[i].cell = pbc_info.cell

natoms = len(pbc_info)
#make master component_list
master_component_list = np.zeros(len(pbc_info), dtype=int)
master_component_list[144:288] = 1
for idx in range(57):
    master_component_list[288 + idx * 8 : 297 + idx * 8] = idx + 2 
master_n_components=59

h_cutoffs = neighborlist.natural_cutoffs(pbc_info)
#print(h_cutoffs)
for atom in fulltraj[-1]:
    if atom.symbol == 'H' or atom.symbol == 'O':
        h_cutoffs[atom.index] += 0.115# Yields H-O cutoff of 1.8A

#with np.printoptions(threshold=np.inf):
#    print(find_h_bonds(fulltraj[-1]))

total_nhbonds = np.zeros(len(fulltraj))
average_nhbonds = np.zeros(len(fulltraj))
max_chain = np.zeros(len(fulltraj))

#neighborList = neighborlist.NeighborList(neighborlist.natural_cutoffs(fulltraj[-1]), self_interaction=False, bothways=True)
neighborList = neighborlist.NeighborList(h_cutoffs, self_interaction=False, bothways=True)
for img in range(len(fulltraj)):
    print(img)
    neighborList.update(fulltraj[img])
    matrix = neighborList.get_connectivity_matrix()
    n_components, component_list = sparse.csgraph.connected_components(matrix)
    unique, counts = np.unique(component_list, return_counts=True)
    max_chain_length = counts[2:].max() / 8
    this_matrix, this_nhbonds = find_h_bonds(fulltraj[img])
    total_nhbonds[img] = this_nhbonds
    average_nhbonds[img] = this_nhbonds / 57
    max_chain[img] = max_chain_length
    print("Image: ", img)
    print(this_nhbonds, max_chain_length)
    #print(max_chain_length)
    #with np.printoptions(threshold=np.inf):
    #    print(n_components)
    #    print(component_list)
    #    print(unique, counts)


#with open("hbond_analysis.txt","w+") as outfile:
#    print("#\"Total number of H-bonds\tAvg number of H-bonds\tMax chain length",file=outfile)
#    for i in zip(total_nhbonds,average_nhbonds,max_chain):
#        print(i, file=outfile)
    

#import inspect
#print(inspect.getfile(neighborlist))

#this_matrix, this_nhbonds = find_h_bonds(fulltraj[-1])
#print("\n",this_nhbonds)
#print(this_matrix)
#with np.printoptions(threshold=np.inf):
#    print(this_matrix)


#n_components, component_list = sparse.csgraph.connected_components(matrix)
#with np.printoptions(threshold=np.inf):
#    print(n_components)
#    print(component_list)

#neighborList = neighborlist.NeighborList(h_cutoffs, self_interaction=False, bothways=True)
#neighborList.update(fulltraj[-1])
#matrix = neighborList.get_connectivity_matrix()
#n_components, component_list = sparse.csgraph.connected_components(matrix)







#for molno in range(2,master_n_components): 
#    molIdxs = [ i for i in range(len(master_component_list)) if master_component_list[i] == molIdx ]
#    print("The following atoms are part of molecule {}: {}".format(molno, molIdxs))



#hbonds = np.zeros(shape=(57,57,len(fulltraj)))

#for i in range(len(fulltraj)):
    #print(i)


#h3po4_indices = [atom.index for atom in fulltraj[-1] if atom.symbol=='P']
#h3po4_indices = np.array(h3po4_indices)
