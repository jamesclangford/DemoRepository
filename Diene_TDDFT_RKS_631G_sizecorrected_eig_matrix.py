#!/usr/bin/env python
#
# Author: Qiming Sun <osirpt.sun@gmail.com>
#

'''
TDDFT analysis of trans-1,3-butadiene, using fast method
code by James Langford
Yang Research Group 
9-1-2020
'''

from pyscf import gto, scf, dft, tdscf
from pyscf.tdscf.rks import TDDFT
import xlwt
import numpy

numpy.set_printoptions(threshold = numpy.inf)

book = xlwt.Workbook(encoding="utf-8")
sheet1 = book.add_sheet("Raw Data")
sheet2 = book.add_sheet("Mol_info")
sheet3 = book.add_sheet("Dot Product")

f1 = open('Diene_fast_eigenvector_RKS_631G_SF_000_dotmatrix', 'w')
f2 = open('Diene_fast_eigenvector_RKS_631G_SF_100_dotmatrix', 'w')  
f3 = open('Diene_fast_eigenvector_RKS_631G_dotmatrix', 'w') 
g = open('Diene_fast_sizecorrected_MO_orbitals_TDDFT_RKS_631G_eig_dotmarix', 'w')

mol = gto.Mole()
mol.cart = True
#mol.charge = 0
#mol.spin = 0

atom = '''C 0.00000000 0.00000000 0.00000000
	C 0.00000000 -1.34527300 0.00000000
	C -1.20530000 -2.16958300 0.00000000
	C -1.20530000 -3.51485600 0.00000000
	H -2.12642000 -4.08847500 0.00000000
	H -0.27931500 -4.08543500 0.00000000
	H -2.15458500 -1.63280300 0.00000000
	H 0.94928500 -1.88205300 0.00000000
	H -0.92598500 0.57057900 0.00000000
	H 0.92112000 0.57361900 0.00000000''' # in Angstrom

basis = '631g'
symmetry = True

mol.build(atom = atom, basis = basis, symmetry = symmetry)
mf = dft.RKS(mol)
mf.xc = 'b3lyp'
mf.kernel()
e_tot = mf.energy_tot()   

sheet2.write(1,1, 'coordinates = ' +  atom)
sheet2.write(2,1, 'functional = ' + mf.xc) 
sheet2.write(3,1, 'basis set = ' + basis)

f1.write('coordinates = ' +  atom)
f1.write('\nfunctional = ' + mf.xc)
f1.write('\nbasis set = ' + basis)
f1.write('\n')

f2.write('coordinates = ' +  atom)
f2.write('\nfunctional = ' + mf.xc) 
f2.write('\nbasis set = ' + basis)
f2.write('\n')  

f3.write('coordinates = ' +  atom)
f3.write('\nfunctional = ' + mf.xc) 
f3.write('\nbasis set = ' + basis)
f3.write('\n')  

g.write('coordinates = ' +  atom)
g.write('\nfunctional = ' + mf.xc) 
g.write('\nbasis set = ' + basis)
g.write('\n')

sheet2.write(4,1, 'total energy = ' + str(e_tot))
f1.write('RKS energy at scale factor 1 = ' + str(e_tot) + '\n') 
f2.write('RKS energy at scale factor 1 = ' + str(e_tot) + '\n') 
f3.write('RKS energy at scale factor 1 = ' + str(e_tot) + '\n') 
g.write('RKS energy at scale factor 1 = ' + str(e_tot) + '\n') 

mytd = TDDFT(mf)
mytd.nstates =25
e, xy = mytd.kernel()

scale = 0.0
col_index = 1

#Scale Factor = 0
print('scale factor = ', scale)
f1.write('\n scale factor = ' +  str(scale) + '\n')
nstates_len = False
while(nstates_len == False):
    e_scale0, xy_scale0 = mytd.scaled_kernel(scale)
    if len(e_scale0) == mytd.nstates:
        nstates_len = True
print('len e_scale0 = ' + str(len(e_scale0)))
f1.write('\n' + str(xy_scale0))    
sheet1.write(1, col_index, 'scale = ' + str(scale))
sheet1.write(2, col_index, scale)
for i in range(len(e_scale0)):
    sheet1.write(i+3, col_index, e_scale0[i].real)
col_index = col_index + 1

#Scale Factor = 1
scale = 0 
print('scale factor = ', scale)
f2.write('\n scale factor = ' +  str(scale) + '\n')
nstates_len = False
while(nstates_len == False):
    e_scale1, xy_scale1 = mytd.scaled_kernel(scale)
    if len(e_scale1) == mytd.nstates:
        nstates_len = True  
print('len e_scale = ' + str(len(e_scale1)))
f2.write('\n' + str(xy_scale1))
sheet1.write(1, col_index, 'scale = ' + str(scale))
sheet1.write(2, col_index, scale)
for i in range(len(e_scale1)):
    sheet1.write(i+3, col_index, e_scale1[i].real)

#Calculate Dot Product
size0 = len(e_scale0)
size1 = len(e_scale1)
xy_scale0_array = numpy.asarray(xy_scale0)
xy_scale1_array = numpy.asarray(xy_scale1)
dot_matrix = numpy.zeros((size1,size0))
xy0_shape = numpy.shape(xy_scale0_array[0,0])
xy1_shape = numpy.shape(xy_scale1_array[0,0])
print('\nxy_i_shape = ' + str(xy0_shape) + '\n')
#iterate over scale factor = 1
for i in range(0, size1):
    #convert matrix representation of eigenvector to vector representation 
    vector1 = numpy.array([]) 
    for v1 in range(0, xy1_shape[0]):
        vector1 = numpy.concatenate([vector1, xy_scale1_array[i,0,v1]])

    #iterate over scale factor = 0
    for j in range(0, size0):
        #convert matrix representation of eigenvector to vector representation 
        vector0 = numpy.array([])
        for v2 in range(0, xy0_shape[0]):
            vector0 = numpy.concatenate([vector0, xy_scale0_array[j,0,v2]])
        #if number is trivially small, enter as "0" for ease of reading file
        if numpy.dot(vector0, vector1) < .01: 
            dot_matrix[i][j] = 0
        else:
            dot_matrix[i][j] = numpy.dot(vector0,vector1)
f3.write('\nconvention: rows are scale factor = 1, columns are scale factor = 0\n')
f3.write('\ndot product matrix = \n' + str(dot_matrix))

print('\norbital energies of RKS calculation')
mf.analyze()
print()
mytd.analyze()

#g.write('\n')
#g.write(a)

book.save('Scaled_TDDFT_RKS_631G_Diene_sizecorrected_result_eV_fast_eig_dotmatrix.xls')

f1.close()
f2.close()
f3.close()
g.close()
