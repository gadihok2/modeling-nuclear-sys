import numpy as np
from numpy.linalg import inv, eig, norm, eigh
import os
import csv
import xlsxwriter

# Read input file

os.chdir(r'/Users/nalingadihoke/Desktop/Spring 2018/npre 247/CP 2')

import pandas as pd
from pandas import ExcelWriter, ExcelFile
df = pd.read_excel('input_cp_2.xls',sheetnumber=0, header=2)

# making lists out of data. These will be used later via indexing. 
Sa = df['Sa']
nSf = df['nSf']
chi = df['c (chi)']
tg_1 = df['two_1']
tg_2 = df['two_2']
Sa1 = df['Sa.1']
nsf1 = df['nSf.1']
chi1 = df['c (chi).1']
eg_1 = df['eight_1']
eg_2 = df['eight_2']
eg_3 = df['eight_3']
eg_4= df['eight_4']
eg_5= df['eight_5']
eg_6= df['eight_6']
eg_7= df['eight_7']
eg_8= df['eight_8']

# 2 Group Eigen values and vectors 

a1 = np.matrix([[Sa[0], 0],[0 , Sa[1]]])
si1 = np.matrix([[0, 0],[tg_1[1] , 0]])
so1 = np.matrix([[tg_1[1], 0],[0 , 0]])
f1 = np.matrix([[nSf[0] ,nSf[1]],[0, 0]])
m1 = a1 + so1 - si1
b1 = inv(m1)*f1
x1,y1 = eig(b1)

# 2 Group Power Iteration Method

phi_1 = np.matrix('1;1')
k_1 = 1
counter = 1

for i in range(20):
	phi_new = (b1*phi_1)/norm((b1*phi_1),ord=2)
	k_new = (((b1*phi_1).transpose())*phi_new)/((phi_new.transpose())*phi_new)
	if k_new>k_1:	
		phi_1 = phi_new
		k_1 = k_new
		counter += 1
	else:
		break


# 8 Group 

# Function to make absorption matrix
def A_creator(x):
	matrix = np.zeros((8,8))
	a = 0
	for i in range(0,8):
		for j in range(0,8):
			if i == j:
				matrix[i,j] = x[a]
				a += 1
			else:
				continue
	return np.asmatrix(matrix)

# Function to make fission matrix
def F_creator(chi,n):
	matrix = np.zeros((8,8))
	a = 0
	b = 0
	for i in range(0,8):
		for j in range(0,8):
			matrix[i,j] = chi[a] * n[b]
		
			if b == 7:
				a += 1
				b = 0
			else:
				b += 1
	return np.asmatrix(matrix)

# Function to make Scatter-in matrix 
def Si_creator(x):
	matrix = np.zeros((8,8))
	a = 0
	for i in range(0,8):
		for j in range(0,8):
			if i == j:
				continue
			else:
				matrix[i,j] = x[a]
				a += 1
	return np.asmatrix(matrix)

# Function to make Scatter-out matrix
def So_creator(m):
	matrix = np.zeros((8,8))
	index = 0
	for i in range(0,8):
		for j in range(0,8):
			if i == j:
				matrix[i,j] = sum(m[index])
				index += 1
			else:
				continue
	return np.asmatrix(matrix)

# Creating Si
list_1 = [0,0,0,0,0,0,0]
list_2 = [eg_1[1],0,0,0,0,0,0]
list_3 = [eg_1[2],eg_2[2],0,0,0,0,0]
list_4 = [eg_1[3],eg_2[3],eg_3[3],0,0,0,0]
list_5 = [0,0,eg_3[4],eg_4[4],0,0,0]
list_6 = [0,0,0,eg_4[5],eg_5[5],eg_7[5],0]
list_7 = [0,0,0,eg_4[6],eg_5[6],eg_6[6],eg_8[6]]
list_8 = [0,0,0,eg_4[7],eg_5[7],eg_6[7],eg_7[7]]
list_final_Si = list_1+list_2+list_3+list_4+list_5+list_6+list_7+list_8
si2 = Si_creator(list_final_Si)

# Creating A
list_for_A = [Sa1[0],Sa1[1],Sa1[2],Sa1[3],Sa1[4],Sa1[5],Sa1[6],Sa1[7]]
list_for_so = [eg_1[0],eg_2[1],eg_3[2],eg_4[3],eg_5[4],eg_6[5],eg_7[6],eg_8[7]]
a2 = A_creator(list_for_A)

# Creating So
sigma_1 = [eg_1[1],eg_1[2],eg_1[3]]
sigma_2 = [eg_2[2],eg_2[3]]
sigma_3 = [eg_3[3],eg_3[4]]
sigma_4 = [eg_4[4],eg_4[5],eg_4[6],eg_4[7]]
sigma_5 = [eg_5[5],eg_5[6],eg_5[7]]
sigma_6 = [eg_6[6],eg_6[7]]
sigma_7 = [eg_7[7],eg_7[5]]
sigma_8 = [eg_8[6]]
sigma = [sigma_1,sigma_2,sigma_3,sigma_4,sigma_5,sigma_6,sigma_7,sigma_8]
so2 = So_creator(sigma)

# Creating f
chi_list = [chi1[0],chi1[1],chi1[2],0,0,0,0,0]
other_list = [nsf1[0],nsf1[1],nsf1[2],nsf1[3],nsf1[4],nsf1[5],nsf1[6],nsf1[7]]
f2 = F_creator(chi_list,other_list)

# Eigen values and vectors
m2 = a2 + so2 - si2
b2 = inv(m2)*f2
x2,y2 = eig(b2)

# 8 Group Power Iteration method
phi_1_8 = np.matrix('1;1;1;1;1;1;1;1')
k_1_8 = 1
counter2 = 1

for i in range(20):

	phi_new_8 = (b2*phi_1_8)/norm((b2*phi_1_8),ord=2)
	k_new_8 = (((b2*phi_1_8).transpose())*phi_new_8)/((phi_new_8.transpose())*phi_new_8)
	if k_new_8>k_1_8:
		phi_1_8 = phi_new_8
		k_1_8 = k_new_8
		counter2 += 1
	else:
		break

# Writing data to output file
workbook = xlsxwriter.Workbook('cp_2_output.xlsx')
worksheet = workbook.add_worksheet()

def writer(a,n,g,c=1):
	worksheet.write(c-1, n, g)
	row = c
	for i in a:
		worksheet.write(row, n,i)
		row+=1

	return

# 2 G data
writer(x1,0,'Eigen values for 2 Group')
writer(y1[:,0],1,'Eigen vector 1 for 2 Group')
writer(y1[:,1],2,'Eigen vector 2 for 2 Group')

# 2 G data power iteration
worksheet.write(0, 3, 'Largest Eigen value 2G')
worksheet.write(1, 3, k_new)
writer(phi_new,4,'Largest Eigen vector 2G')

# 8 G data
writer(x2,0,'Eigen values for 8 Group',5)
writer(np.real(y2[:,0]),1,'Eigen vector 1 for 8 Group',5)
writer(np.real(y2[:,1]),2,'Eigen vector 2 for 8 Group',5)
writer(np.real(y2[:,2]),3,'Eigen vector 3 for 8 Group',5)
writer(np.real(y2[:,3]),4,'Eigen vector 4 for 8 Group',5)
writer(np.real(y2[:,4]),5,'Eigen vector 5 for 8 Group',5)
writer(np.real(y2[:,5]),6,'Eigen vector 6 for 8 Group',5)
writer(np.real(y2[:,6]),7,'Eigen vector 7 for 8 Group',5)
writer(np.real(y2[:,7]),8,'Eigen vector 8 for 8 Group',5)

# 8 G data power iteration
worksheet.write(16, 0, 'Largest Eigen value 8G')
worksheet.write(17, 0, k_new_8)
writer(phi_new_8[:,0],1,'Largest Eigen vector 8G',16)

width= len("Eigen vector 2 for 2 Group")
worksheet.set_column(0, 22, width)

workbook.close()


