import numpy as np
import matplotlib.pyplot as plt
import math
import copy
import string
import cmath
from scipy.linalg import block_diag
import itertools
from itertools import permutations, repeat
import random
from matplotlib import cm
import matplotlib.path as mpltPath
from termcolor import colored



#model parameters

onsite = 0.
dx_onsite = -0.1
dy_onsite = -0.1
dxsqrt_onsite = -1.2
dxy_onsite = -1.2
tdsqrt = -0.5
delta = 2.5
t = -2.0
td = 0.4
en_fermi = 0.071


'''
onsite = 0.
dx_onsite = -10
dy_onsite = -10
dxsqrt_onsite = -10
dxy_onsite = -10
tdsqrt = 0.
delta = 0.
t = -2.0
td = 0.
'''

#print('td-x2_y2 = ',tdsqrt)
#print('td-xz,yz = ',td)

#original graphene lattice vectors
lat = np.array([[3./2.,-np.sqrt(3.)/2.],[3./2.,np.sqrt(3.)/2.]])

#define coordinates of two pz orbitals in unit cell
orb = np.array([[1./3.,1./3.],[2./3.,2./3.]])


#defines the unit cell itself reasonable are of the form vec = [n, 3m+n]

vec = [-1,5]

#defines the number of unit cells in the lattice
length = 30





def super_hexagon(vec):
    vec2=[-vec[1],vec[1]+vec[0]]
    return [vec,vec2]

coef = super_hexagon(vec)
mega_lat = np.dot(coef,lat)
imp_loc = np.array([0.,0.])



#make a grid
expand = list(itertools.product(np.arange(-10,10,1), np.arange(-10,10,1)))
grid = []
for lat_vec in expand:
    for orbit in orb:
        tmp = orbit + lat_vec
        grid.append(np.array(tmp))

def cartez_coord(point,lat_vecs):
    x_val = np.array(point[0]*lat_vecs[0,0]+point[1]*lat_vecs[1,0])
    y_val = np.array(point[0]*lat_vecs[0,1]+point[1]*lat_vecs[1,1])
    return np.array([x_val,y_val])


#print('checks', np.dot([0.5,0.5],lat),cartez_coord([0.5,0.5],lat))

cartez_grid = [cartez_coord(point,lat) for point in grid]

unit_hex = np.array([[1/3,1/3],[-1/3,2/3],[-2/3,1/3],[-1/3,-1/3],[1/3,-2/3],[2/3,-1/3],[1/3,1/3]])

hex_loc  = [cartez_coord(i,mega_lat) for i in unit_hex]

path_hex = mpltPath.Path(hex_loc)


def insider_check(pnt, path , ends):
    if path.contains_points([pnt],radius = 1e-9) == True:
        return True
    else:
        return False

list_inside   = []
for gen_point in grid:
    if insider_check(cartez_coord(gen_point,lat), path_hex, hex_loc) == True:
        list_inside.append(gen_point)
list_inside.append(np.array([0.,0.]))
list_inside.append(np.array([0.,0.]))
list_inside.append(np.array([0.,0.]))
list_inside.append(np.array([0.,0.]))

inside_hexagon  = np.array(list_inside)
inside          = [cartez_coord(point,lat) for point in inside_hexagon]
mega_orbits     = np.array([np.dot(np.linalg.inv(np.array(mega_lat).T),orb) for orb in inside])


expand = [[0,i] for i in range(length)]
unit   = [ ]
for lat_vec in expand:
    for orbit in mega_orbits:
        tmp = orbit + lat_vec
        unit.append(np.array(tmp))


inside  = [cartez_coord(point,mega_lat) for point in unit]


def unit_cell_drawer(grid):
    fig,ax = plt.subplots(figsize=(10,10))
    cartez_grid = [cartez_coord(point,mega_lat) for point in unit]
    x_vals = np.array([cartez_grid[i][0] for i in range(len(cartez_grid))])
    y_vals = np.array([cartez_grid[i][1] for i in range(len(cartez_grid))])
    ax.scatter(x_vals, y_vals, marker='.',color = 'red')



    #start  = [r[0] for r in hex_loc]
    #finish = [r[1] for r in hex_loc]
    #plt.plot(start,finish,'ro-',color='b',markersize=2)
    #for i in range(len(hex_loc)):
    #    ax.annotate(i, (start[i], finish[i]))


    x_inside = np.array([inside[i][0] for i in range(len(inside))])
    y_inside = np.array([inside[i][1] for i in range(len(inside))])
    #ax.scatter(x_inside, y_inside,marker='o',color='r')
    start  = [r[0] for r in inside]
    finish = [r[1] for r in inside]
    for i in range(len(inside)):
        ax.annotate(i, (start[i], finish[i]))


    plt.ylabel('y')
    plt.xlabel('x')
    plt.xlim(-20,20)
    plt.ylim(-10,20)
    #plt.savefig('hexagon_obc_c6symmetric.pdf')
    return plt.show()
#unit_cell_drawer(grid)





def ribbon_drawer():
    fig,ax = plt.subplots(figsize=(10,10))

    hex_grid  = []
    grid_list = []
    for span in [[i,0] for i in range(20)]:
        for orbit in unit:
                tmp = orbit + span
                grid_list.append(np.array(tmp))
    hex_grid  = np.array(grid_list)

    grid_cartez = [cartez_coord(i,mega_lat) for i in hex_grid]
    x_vals = np.array([grid_cartez[i][0] for i in range(len(grid_cartez))])
    y_vals = np.array([grid_cartez[i][1] for i in range(len(grid_cartez))])
    ax.scatter(x_vals, y_vals, marker='.',color='b')



    ax.scatter(mega_lat[0][0],mega_lat[0][1],color='g',marker='.')
    ax.scatter(mega_lat[1][0],mega_lat[1][1],color='g',marker='.')
    ax.scatter(0.,0.,color='g',marker='.')



    plt.ylabel('y')
    plt.xlabel('x')
    plt.xlim(-50,50)
    plt.ylim(-10,90)
    #plt.savefig('hexagon_obc_c6symmetric.pdf')
    return plt.show()

#ribbon_drawer()





def hopping_list(coef):
    R_list   = [[1,0],[0,1],[-1,1]]
    hop_list = []

    for subset in itertools.combinations(mega_orbits[:-4], 2):
        s  = np.sqrt(((subset[1][0]-subset[0][0])*mega_lat[0][0]+(subset[1][1]-subset[0][1])*mega_lat[1][0])**2+((subset[1][0]-subset[0][0])*mega_lat[0][1]+(subset[1][1]-subset[0][1])*mega_lat[1][1])**2)
        if np.isclose(s,1.):
            tmp = []
            n = mega_orbits.tolist().index(list(subset[0]))
            p = mega_orbits.tolist().index(list(subset[1]))
            tmp.append(p)
            tmp.append(n)
            tmp.append([0,0])
            hop_list.append(tmp)
    for R in R_list:
        neigh_list = np.array([R+orb for orb in mega_orbits[:-4]])
        for orb in mega_orbits[:-4]:
            for neigh in neigh_list:
                l = np.sqrt(((neigh[0]-orb[0])*mega_lat[0][0]+(neigh[1]-orb[1])*mega_lat[1][0])**2+((neigh[0]-orb[0])*mega_lat[0][1]+(neigh[1]-orb[1])*mega_lat[1][1])**2)
                if np.isclose(l,1.):
                    tmp = []
                    g   = mega_orbits.tolist().index(list(orb))
                    z   = neigh_list.tolist().index(list(neigh))
                    tmp.append(g)
                    tmp.append(z)
                    tmp.append(R)
                    hop_list.append(tmp)
    return hop_list
def onsite_gen(coef,imp_loc):
    int_orb     = np.array([np.dot(np.linalg.inv(np.array(coef).T),vec) for vec in orb])
    dist        = np.sqrt(((int_orb[1][0]-int_orb[0][0])*mega_lat[0][0]+(int_orb[1][1]-int_orb[0][1])*mega_lat[1][0])**2+((int_orb[1][0]-int_orb[0][0])*mega_lat[0][1]+(int_orb[1][1]-int_orb[0][1])*mega_lat[1][1])**2)
    on_site     = []
    for site in mega_orbits:
        s = np.sqrt(((site[0]-imp_loc[0])*mega_lat[0][0]+(site[1]-imp_loc[1])*mega_lat[1][0])**2+((site[0]-imp_loc[0])*mega_lat[0][1]+(site[1]-imp_loc[1])*mega_lat[1][1])**2)
        if np.isclose(s,dist):
            n = mega_orbits.tolist().index(list(site))
            on_site.append(n)
    return on_site


def distance(a,b):
    return np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)
def rot_matrix(num_rot,vec):
    tmp = np.dot(vec,mega_lat)
    rot_matrix = np.array([[np.cos(num_rot*2*np.pi/6),-np.sin(num_rot*2*np.pi/6)],[np.sin(num_rot*2*np.pi/6),np.cos(num_rot*2*np.pi/6)]])
    rot_vec    = np.dot(rot_matrix,tmp)
    return np.dot(rot_vec,np.linalg.inv(mega_lat))
def d_hoppings(mega_orbits):
    unit_size = len(mega_orbits) - 4
    onsites   = onsite_gen(coef,imp_loc)
    xz_list   = []
    dhoppings = []
    vec       = mega_orbits[onsites[0]]
    for num_rot in range(6):
        check_vec  = rot_matrix(num_rot,vec)
        check_list = list([check_vec[0],check_vec[1]])
        index      = [mega_orbits.tolist().index(list(i)) for i in mega_orbits if np.all(np.isclose(check_list,i))]
        xz_list.append(index[0])
    dhoppings.append([unit_size,xz_list[0],2*td/np.sqrt(3)])
    dhoppings.append([unit_size,xz_list[1],td/np.sqrt(3)])
    dhoppings.append([unit_size,xz_list[2],-td/np.sqrt(3)])
    dhoppings.append([unit_size,xz_list[3],-2*td/np.sqrt(3)])
    dhoppings.append([unit_size,xz_list[4],-td/np.sqrt(3)])
    dhoppings.append([unit_size,xz_list[5],td/np.sqrt(3)])
    dhoppings.append([unit_size+1,xz_list[1],-td])
    dhoppings.append([unit_size+1,xz_list[2],-td])
    dhoppings.append([unit_size+1,xz_list[4],td])
    dhoppings.append([unit_size+1,xz_list[5],td])
    return dhoppings
def dsqrt_hoppings(mega_orbits):
    unit_size = len(mega_orbits) - 4
    onsites   = onsite_gen(coef,imp_loc)
    xz_list   = []
    dhoppings = []
    vec       = mega_orbits[onsites[0]]
    for num_rot in range(6):
        check_vec   = rot_matrix(num_rot,vec)
        check_list  = list([check_vec[0],check_vec[1]])
        index       = [mega_orbits.tolist().index(list(i)) for i in mega_orbits if np.all(np.isclose(check_list,i))]
        xz_list.append(index[0])
    dhoppings.append([unit_size+3,xz_list[0],tdsqrt])
    dhoppings.append([unit_size+3,xz_list[1],-tdsqrt/2])
    dhoppings.append([unit_size+3,xz_list[2],-tdsqrt/2])
    dhoppings.append([unit_size+3,xz_list[3],tdsqrt])
    dhoppings.append([unit_size+3,xz_list[4],-tdsqrt/2])
    dhoppings.append([unit_size+3,xz_list[5],-tdsqrt/2])
    dhoppings.append([unit_size+2,xz_list[1],tdsqrt*np.sqrt(3)/2])
    dhoppings.append([unit_size+2,xz_list[2],-tdsqrt*np.sqrt(3)/2])
    dhoppings.append([unit_size+2,xz_list[4],tdsqrt*np.sqrt(3)/2])
    dhoppings.append([unit_size+2,xz_list[5],-tdsqrt*np.sqrt(3)/2])
    return dhoppings



hoppings   = hopping_list(coef)
onsites    = onsite_gen(coef,imp_loc)
d_hops     = d_hoppings(mega_orbits)
dsqrt_hops = dsqrt_hoppings(mega_orbits)




def hex_tb_ham(k_vec,orb):
    unit_size  = len(orb)
    tb_ham     = np.zeros((unit_size,unit_size), dtype=complex)
    for i in range(unit_size):
        tb_ham[i,i] = onsite
    for s in onsites:
        tb_ham[s,s] += delta
    for hop in hoppings:
        i = int(hop[0])
        j = int(hop[1])
        r_vec = hop[2]
        #r_vec = r+orb[j]-orb[i]
        phase = np.exp(2.0j*np.pi*np.dot(k_vec,r_vec))
        amp   = t*phase
        tb_ham[i,j] += amp
        tb_ham[j,i] += amp.conjugate()
    for d in d_hops:
        i = int(d[0])
        j = int(d[1])
        r_vec = [0,0]
        #r_vec = r + orb[j]
        phase = np.exp(2.0j*np.pi*np.dot(k_vec,r_vec))
        amp   = d[2]*phase
        tb_ham[i,j] += amp
        tb_ham[j,i] += amp.conjugate()
    for dsqrt in dsqrt_hops:
        i = int(dsqrt[0])
        j = int(dsqrt[1])
        r_vec = [0,0]
        #r_vec = r+orb[j]
        phase = np.exp(2.0j*np.pi*np.dot(k_vec,r_vec))
        amp   = dsqrt[2]*phase
        tb_ham[i,j] += amp
        tb_ham[j,i] += amp.conjugate()
    tb_ham[unit_size-4,unit_size-4] = dx_onsite
    tb_ham[unit_size-3,unit_size-3] = dy_onsite
    tb_ham[unit_size-2,unit_size-2] = dxsqrt_onsite
    tb_ham[unit_size-1,unit_size-1] = dxy_onsite
    return tb_ham






hoppings_inside = []
hoppings_neigh  = []
hoppings_inbetween   = []

for hop in hoppings:
    if hop[2]  != [0,0] and hop[2]  != [0,1]:
        if hop[2] == [-1,1]:
            tmp = [hop[0],hop[1],[0,1]]
            hoppings_neigh.append(tmp)
        else:
            tmp = [hop[0],hop[1],[0,0]]
            hoppings_neigh.append(tmp)
    if hop[2]  == [0,1]:
        hoppings_inbetween.append(hop)
    if hop[2] == [0,0]:
        hoppings_inside.append(hop)

#print('inside-t',hoppings_inside)
#print('exp-outside',hoppings_neigh)
#print('between t also',hoppings_inbetween)

def unit_cell_tb(orb):
    unit_size  = len(orb)
    tb_ham     = np.zeros((unit_size,unit_size), dtype=complex)
    for i in range(unit_size):
        tb_ham[i,i] = onsite
    for s in onsites:
        tb_ham[s,s] += delta
    for hop in hoppings_inside:
        i = int(hop[0])
        j = int(hop[1])
        amp   = t
        tb_ham[i,j] += amp
        tb_ham[j,i] += amp.conjugate()
    for d in d_hops:
        i = int(d[0])
        j = int(d[1])
        amp   = d[2]
        tb_ham[i,j] += amp
        tb_ham[j,i] += amp.conjugate()
    for dsqrt in dsqrt_hops:
        i = int(dsqrt[0])
        j = int(dsqrt[1])
        amp   = dsqrt[2]
        tb_ham[i,j] += amp
        tb_ham[j,i] += amp.conjugate()
    tb_ham[unit_size-4,unit_size-4] = dx_onsite
    tb_ham[unit_size-3,unit_size-3] = dy_onsite
    tb_ham[unit_size-2,unit_size-2] = dxsqrt_onsite
    tb_ham[unit_size-1,unit_size-1] = dxy_onsite
    return tb_ham




def orbit_localizer(cell_loc,orbit_loc):
    number_inexpand  = expand.index(list(cell_loc))
    new_loc = number_inexpand * len(mega_orbits) + orbit_loc
    return new_loc


#print('HERE',orbit_localizer([0,1],38))


np.set_printoptions(linewidth = np.inf)

tb_unit_hexagon = unit_cell_tb(mega_orbits)

#print(expand)


def ribbon_tb(k):
    final_hex = np.zeros((len(mega_orbits)*len(expand),len(mega_orbits)*len(expand)),dtype = complex)
    final_hex       = block_diag(*([tb_unit_hexagon] * (len(expand))))
    for unit in expand:
        for hop_el in hoppings_inbetween:
            one_neighb     = [sum(x) for x in zip(list(unit), hop_el[2])]
            if any(one_neighb == local for local in expand):
                i_value = orbit_localizer(list(unit),hop_el[0])
                j_value = orbit_localizer(one_neighb,hop_el[1])
                final_hex[i_value,j_value] = t
                final_hex[j_value,i_value] = t
    for unit in expand:
        for hop_el in hoppings_neigh:
            one_neighb = [sum(x) for x in zip(list(unit), hop_el[2])]
            #print(unit, hop_el,one_neighb)
            if any(one_neighb == local for local in expand):
                i_value = orbit_localizer(list(unit),hop_el[0])
                j_value = orbit_localizer(list(one_neighb),hop_el[1])
                tmp1 = [i_value,j_value]
                #print('inside',unit, hop_el,one_neighb)
                #print('r from here',hop_el[2])
                if hop_el[2] ==[0,0]:
                    r = 1
                if hop_el[2] ==[0,1]:
                    r = -1
                #print('i-j',tmp1,'R',r)
                phase = np.exp(2.0j*np.pi*k*r)
                amp   = t*phase
                final_hex[i_value,j_value] += amp
                final_hex[j_value,i_value] += amp.conjugate()
                #print('goodones',i_value,j_value,amp)
                #print('goodones',j_value,i_value,amp.conjugate())
    return final_hex



def hermitian_check(matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if np.transpose(matrix)[i,j] == matrix[i,j]:
                print('hermitian')
            else:
                print('NOT')

#hermitian_check(ribbon_tb(0.2))

np.set_printoptions(linewidth = np.inf)
#print('matrix_size',len(ribbon_tb(0.0)))
#print('size of verticel unit cell',len(mega_orbits)*len(expand))


def chop(matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            matrix[i][j] = round(np.real(matrix[i][j]),2)+1j*round(np.imag(matrix[i][j]),2)
    return matrix

#print('resulting matrix after chop')
np.set_printoptions(linewidth = np.inf)
#print(chop(ribbon_tb(0.)))

#print(ribbon_tb(0.2)[0,84])
#print(ribbon_tb(0.2)[92,130],-2*np.exp(1j*2*np.pi*0.2))


k_graph  = np.linalg.inv(lat).T
k_mega_lat = np.dot(k_graph.T,np.linalg.inv(coef)).T



def ribbon_energy_bands(k_list):
    eval_list = []
    for k in k_list:
        eval  = np.linalg.eigvalsh(ribbon_tb(k))
        eval_list.append(eval)
    return eval_list




nk = 100

def ribbon_spectrum(nk):
    #basis vectors in k space
    #k_graph = np.linalg.inv(lat).T
    global mega_orbits
    k_mega_lat = np.dot(k_graph.T,np.linalg.inv(coef)).T

    kx_list = np.arange(-1/2,1/2,1/nk)
    k_dist = []
    for point in kx_list:
        k_el   = k_mega_lat[1]*point
        l = np.sqrt(k_el[0]**2+k_el[1]**2)
        k_dist.append(l)
    #print('dist',k_dist)

    evals = np.array(ribbon_energy_bands(kx_list))



    file = open("energy_list_ribbon.txt","w")
    for i in evals:
        for l in i:
            if l != i[-1]:
                file.write(str(l)+',')
            if l == i[-1]:
                file.write(str(l))
        file.write('\n')
    file.close()

    evals = evals.T
    #print(len(evals[0]))
    fig, ax = plt.subplots()
    #ax.set_xlabel('k1')
    ax.set_xlabel('tdxz,yz='+str(td)+',tdx2_y2,xy='+str(tdsqrt)+'.')

    mega_size = len(mega_orbits)*len(expand)
    for g in range(mega_size):
        ax.plot(kx_list,evals[g],color='darkblue')

    fig.tight_layout()
    plt.ylim(0.,1.)
    fig.savefig('ribbon_bands.pdf')
    plt.xlabel('k')
    plt.ylabel('energy')
    plt.show()

    return


ribbon_spectrum(nk)
