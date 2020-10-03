import numpy as np
from scipy.integrate import odeint
import pylab as plt
import os
import csv
import xlsxwriter

#some useful functions

def decay_constant(t): # to calculate decay constant fron half life which is the input
	constant = np.log(2) / t
	return constant


os.chdir(r'/Users/nalingadihoke/Desktop/Spring 2018/npre 247/CP 1')

def read_data( input_file ):
   d = { }
   for line in open( input_file ):
       line = line.rstrip( )
       parameter,value = line.split( "," )
       d[parameter] =  float(value)
   return d

d1 = read_data('input.csv')

#initial parameters extracted from input file 'input.csv'


t_half_a = d1['t_half_A (hours)'] # hours
t_half_b = d1['t_half_B (hours)']

N_0_a = d1['N(0)_A'] # initial conditions
N_0_b = d1['N(0)_B']
N_0_c = d1['N(0)_C']

t_final = d1['t_final'] # hours

# values of timestep (to be used in part 2, the numerical solutions)
t_del1 = 1 
t_del2 = 0.5
t_del3 = 0.25
t_del4 = 0.125

# decay constants
lambda_a = decay_constant(t_half_a) 
lambda_b = decay_constant(t_half_b)

N_0_list = [N_0_a, N_0_b,N_0_c] # a list of initial number of Na, Nb and Nc
lambda_list = [lambda_a,lambda_b] # a list of decay constants of a and b

t_max = (np.log(lambda_list[1]/lambda_list[0]))/(lambda_list[1] - lambda_list[0])

#part 1 - ANALYTICAL SOLLUTION. Here we solve the three odes and store the three array values of Na, Nb and Nc over [0 , t_final] in n as a list.

t = np.linspace(0,t_final, 1500) # 1500 data points considered for a smooth analytical solution in Figure 1

def model_b(N_0_list,t):
	na = N_0_list[0]
	nb = N_0_list[1]
	nc = N_0_list[2]
	dnadt = -lambda_a*na
	dnbdt = (-lambda_b*nb)+(lambda_a*na)
	dncdt = lambda_b*nb
	return [dnadt,dnbdt,dncdt]



n = odeint(model_b,N_0_list,t) # n = [Na,Nb,Nc] where each element is an array

'''plt.plot(t,n[:,0],'b-')
plt.plot(t,n[:,1],'r--')
plt.plot(t,n[:,2],'g')
plt.show()'''

#PART 2. NUMERICAL SOLUTION 

'''This function accepts lists of initial number of nuclides of the three elements and of their decay constants.
Also accepts a timestep (Δt) and returns a list [na,nb,nc,t] where:
na, nb and nc are the arrays containg the numerical solution using the finite difference method for Δt
t is an array with the time from 0 to t_final in Δt timesteps'''

def numerical(N_0_list,lambda_list,t_del):
	n = int((50/t_del)+1)
	na = np.zeros(n)
	nb = np.zeros(n)
	nc = np.zeros(n)
	t = np.linspace(0,50,((50/t_del)+1))
	na[0] = N_0_list[0]
	nb[0] = N_0_list[1]
	nc[0] = N_0_list[2] 
	for i in range(1,n):
		na[i] = (- lambda_list[0] * na[i-1] * t_del) + na[i-1]
		nb[i] = ((-lambda_list[1]*nb[i-1]) + (lambda_list[0]*na[i-1]))*t_del + nb[i-1]
		nc[i] = (lambda_list[1]*nb[i-1]*t_del) + nc[i-1]
	return [na,nb,nc,t]

#Using above function for different values of Δt and assigning necessary values to variables for plotting

sol1 = numerical(N_0_list,lambda_list,t_del1)
nb_num_1 = sol1[1]
t_num_1 = sol1[3]

sol2 = numerical(N_0_list,lambda_list,t_del2)
nb_num_2 = sol2[1]
t_num_2 = sol2[3]

sol3 = numerical(N_0_list,lambda_list,t_del3)
nb_num_3 = sol3[1]
t_num_3 = sol3[3]

sol4 = numerical(N_0_list,lambda_list,t_del4)
na_num_4 = sol4[0]
nb_num_4 = sol4[1]
nc_num_4 = sol4[2]
t_num_4 = sol4[3]
n_total_4 = na_num_4 + nb_num_4 + nc_num_4

# Plotting statements for Number of radionuclides of Nb vs time for three values of Δt

plt.figure(1)
plt.plot(t_num_1,nb_num_1,'b-', label='delta_t = 1 (coarse)')
plt.plot(t_num_2,nb_num_2,'r--', label='delta_t = 0.5 (medium)')
plt.plot(t_num_3,nb_num_3,'g', label='delta_t = 0.25 (fine)')
#plt.plot(t_num_3,nb_num_3,'y', label='Δt = 0.125')
plt.plot(t,n[:,1],'k', label='Analytical solution Nb(t)')
axes = plt.gca()
axes.set_xlim([0,50])
axes.set_ylim([0,85])
plt.ylabel('Number of radionuclides of Nb')
plt.xlabel('Time (hours)')
plt.legend(loc='upper right')
plt.title('Number of radionuclides of Nb vs time for three values of delta_t')
plt.savefig("/Users/nalingadihoke/Desktop/Spring 2018/npre 247/CP 1/Nb_vs_Time.png")

# Plotting statements of Number of radionuclides of Na,Nb and Nc vs time 


plt.figure(2)
plt.plot(t_num_4,na_num_4,'b-', label='Na(t)')
plt.plot(t_num_4,nb_num_4,'r--', label='Nb(t)')
plt.plot(t_num_4,nc_num_4,'g', label='Nc(t)')
plt.plot(t_num_4,n_total_4,'y', label='Na(t) + Nb(t) + Nc(t')
axes1 = plt.gca()
axes1.set_xlim([0,50])
axes1.set_ylim([-1,105])
plt.ylabel('Number of radionuclides')
plt.xlabel('Time (hours)')
plt.legend(loc='center right')
plt.title('Number of radionuclides of Na,Nb and Nc vs time (Δt = 0.125)')
plt.savefig("/Users/nalingadihoke/Desktop/Spring 2018/npre 247/CP 1/radionuclides_vs_time.png")


# For time of maximum nb vs. 1/t
#considering 310 values of 1/del_t between 1 and 310

one_over_delt = np.linspace(1,310,310)
#this means del_t ranges from 1 to 0.0032258 in 310 time-steps
timeformax = np.zeros(310)
for j in range(0,310):
	del_t = 1/one_over_delt[j]
	sol = numerical(N_0_list,lambda_list,del_t)
	time_array = sol[3]
	timeformax[j] = time_array[np.argmax(sol[1])]

time_for_max_nb_analytic = t[np.argmax(n[:,1])] * np.ones(len(one_over_delt))

#print(time_for_max_nb_analytic[-2], timeformax[-1])

# Plotting statements for Time to reach maximum Nb vs 1/Δt

plt.figure(3)
plt.plot(one_over_delt,timeformax, 'b-', label='Numerical solutions')
plt.plot(one_over_delt,time_for_max_nb_analytic,'r',label='Analytic solution = 3.835')
axes2 = plt.gca()
axes2.set_xlim([0,315])
axes2.set_ylim([2.8,3.9])
plt.ylabel('Time to reach maximum Nb (hours)')
plt.xlabel('1/Δt (hours^-1)')
plt.legend(loc='center right')
plt.title('Time to reach maximum Nb vs 1/Δt')
plt.savefig("/Users/nalingadihoke/Desktop/Spring 2018/npre 247/CP 1/max_Nb_vs_oneoverdel_time.png")


# Writing data into output file

workbook = xlsxwriter.Workbook('output.xlsx')
worksheet = workbook.add_worksheet()

def writer(a,n,g):
	worksheet.write(0, n, g)
	row = 1
	for i in a:
		worksheet.write(row, n,i)
		row+=1
	return

# data for replication of figure 1
writer(t,0,'Time steps for Anaylytical solution')
writer(n[:,1],1,'Analytical solutions of Nb')
worksheet.write(0,2, 'Initial Parameters')
worksheet.write(1, 2, 't_half_A')
worksheet.write(2, 2, 't_half_B')
worksheet.write(3, 2, 'N(0)A')
worksheet.write(4, 2, 'N(0)B')
worksheet.write(5, 2, 'N(0)C')
worksheet.write(6, 2, 't_final')
worksheet.write(1, 3, t_half_a)
worksheet.write(2, 3, t_half_b)
worksheet.write(3, 3, N_0_a)
worksheet.write(4, 3, N_0_b)
worksheet.write(5, 3, N_0_c)
worksheet.write(6, 3, t_final)
writer(t_num_1,4,'Time steps Δt = 1 hour')
writer(nb_num_1,5,'Nb for Δt = 1 hour')
writer(t_num_2,7,'Time steps Δt = 0.5 hour')
writer(nb_num_2,8,'Nb for Δt = 0.5 hour')
writer(t_num_3,10,'Time steps Δt = 0.25 hour')
writer(nb_num_3,11,'Nb for Δt = 0.25 hour')

#data for replication of figure 2
writer(t_num_4,13,'Time steps Δt = 0.125 hour')
writer(na_num_4,14,'Na for Δt = 0.125 hour')
writer(nb_num_4,15,'Nb for Δt = 0.125 hour')
writer(nc_num_4,16,'Nc for Δt = 0.125 hour')
writer(n_total_4,17,'N total for Δt = 0.125 hour')

#data for replication of figure 3
writer(one_over_delt ,19,'1/Δt')
writer(time_for_max_nb_analytic ,20,'time for max Nb analytical approach')
writer(timeformax ,21,'time for max Nb analytical approach numerical approach')

width= len("time for max Nb analytical approach numerical approach")
worksheet.set_column(0, 22, width)

workbook.close()

plt.show()






