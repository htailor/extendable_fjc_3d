#!/usr/bin/python
import os
import math
import time
import multiprocessing
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

from multiprocessing import Pool
from libmath import *

plt.rc('text',usetex=True)
plt.rc('font',**{'family':'serif','serif':['Computer Modern']})

ResultsDir = ['results','results_small_approx']

# Definition of constants in the FJC

link_length=1
extendable_length=0.1


def ConstantParameter(n_, R_):

    return pi*pow(-1,n_)*pow(2*pi,n_-2)/float(R_)

def Eta(k1_, k2_, k3_, k4_):

    return (link_length + half*extendable_length)*(k1_ - k4_) + (link_length - half*extendable_length)*(k2_ - k3_)

def Sj(j1_, j2_, j3_, j4_):

    return Multinomial(j1_, j2_, j3_, j4_)*pow(-1,j2_+j4_)

def Sl(l1_, l2_, l3_, l4_):

    return Multinomial(l1_, l2_, l3_, l4_)*pow(-link_length-half*extendable_length,l1_+l4_)*pow(link_length-half*extendable_length,l2_+l3_)

def Const(n_, alpha_):

    return Binomial(n_,alpha_)*Inverse(Factorial(2*n_+alpha_-2))

def NormalisingFunction(n_, R_):

    return 4*pi*Squared(R_)

def Kernal(n_, R_, j1_, j2_, j3_, j4_, l1_, l2_, l3_, l4_, alpha_):

    return pow(Eta(j1_, j2_, j3_, j4_) + Eta(l1_, l2_, l3_, l4_) - R_,2*n_+alpha_-2)*Sgn(Eta(j1_, j2_, j3_, j4_) + Eta(l1_, l2_, l3_, l4_) - R_)

def PartitionFunction(n_, R_):

    SumPF = 0

    for alpha in range(0,int(n_+1)):
        for l1 in range(0,int(n_-alpha+1)):
            for l2 in range(0,int(n_-alpha-l1+1)):
                for l3 in range(0,int(n_-alpha-l1-l2+1)):
                    l4 = (n_-alpha)-l1-l2-l3
                    for j1 in range(0,int(alpha+1)):
                        for j2 in range(0,int(alpha-j1+1)):
                            for j3 in range(0,int(alpha-j1-j2+1)):
                                j4 = alpha - j1 - j2 - j3
                                SumPF = SumPF + Const(n_,alpha)*Sj(j1, j2, j3, j4)*Sl(l1, l2, l3, l4)*Kernal(n_, R_, j1, j2, j3, j4, l1, l2, l3, l4, alpha)

    return SumPF*ConstantParameter(n_,R_)




def SAConstantParameter(n_, R_):

    return pow(4*pi*link_length*extendable_length,n_)*pi*Inverse(Squared(2*pi)*R_*pow(2,n_)*Factorial(n_-2))

def SAEta(n_, k_, R_):

    return link_length*(n_-2*k_)-R_

def SAKernal(n_, k_, R_):

    return Binomial(n_,k_)*pow(-1,k_)*pow(SAEta(n_, k_, R_),n_-2)*Sgn(SAEta(n_, k_, R_))

def SAPartitionFunction(n_, R_):

    SumPF = 0
    for k in range(0,int(n_+1)):
        SumPF = SumPF + SAKernal(n_, k, R_)

    return SumPF*SAConstantParameter(n_,R_)



def ZR_data(n_):

    nrange = n_ + 1

    output_filename1 = 'n'+str(n_)+'_ZR.txt'
    out1 = open(os.path.join(ResultsDir[0],output_filename1),'w')

    output_filename2 = 'n'+str(n_)+'_ZR_normalised.txt'
    out2 = open(os.path.join(ResultsDir[0],output_filename2),'w')

    
    rx = np.linspace(0.01,nrange,nrange*100)
    zr = []
    zr_norm = []


    for i in range(len(rx)):
        z = PartitionFunction(n_,rx[i])
        zn = z*NormalisingFunction(n_ ,rx[i])

        zr.append(z)
        zr_norm.append(zn)

        out1.write(str(rx[int(i)]) + ',' + str(zr[int(i)]) + '\n')
        out2.write(str(rx[int(i)]) + ',' + str(zr_norm[int(i)]) + '\n')


    output_filename_plot_PDF1 = os.path.join(ResultsDir[0],'n'+str(n_)+'_ZR.pdf')
    output_filename_plot_PDF2 = os.path.join(ResultsDir[0],'n'+str(n_)+'_ZR_normalised.pdf')

    PlotData(rx,zr,'R/a','Z(R)/C^{N}',output_filename_plot_PDF1,'#feb119','None')   
    PlotData(rx,zr_norm,'R/a','4 \pi R^{2} Z(R)/C^{N}',output_filename_plot_PDF2,'#feb119','#FFEFD1')

def ZR_dataSA(n_):
    if n_ > 2:
        nrange = n_ + 1

        output_filename1 = 'n'+str(n_)+'_ZR.txt'
        output_filename2 = 'n'+str(n_)+'_ZR_normalised.txt'

        SAout1 = open(os.path.join(ResultsDir[1],output_filename1),'w')
        SAout2 = open(os.path.join(ResultsDir[1],output_filename2),'w')

        rx = np.linspace(0.01,nrange,nrange*100)
        zr = []
        zr_norm = []

        for i in range(len(rx)):
            z = SAPartitionFunction(n_,rx[i])
            zn = z*NormalisingFunction(n_ ,rx[i])

            zr.append(z)
            zr_norm.append(zn)

        SAout1.write(str(rx[int(i)]) + ',' + str(zr[int(i)]) + '\n')
        SAout2.write(str(rx[int(i)]) + ',' + str(zr_norm[int(i)]) + '\n')

        SAoutput_filename_plot_PDF1 = os.path.join(ResultsDir[1],'n'+str(n_)+'_ZR.pdf')
        SAoutput_filename_plot_PDF2 = os.path.join(ResultsDir[1],'n'+str(n_)+'_ZR_normalised.pdf')

        PlotData(rx,zr,'R/a','Z(R)/C^{N}',SAoutput_filename_plot_PDF1,'#feb119','None')
        PlotData(rx,zr_norm,'R/a','4 \pi R^{2} Z(R)/C^{N}',SAoutput_filename_plot_PDF2,'#feb119','#FFEFD1')



def PlotData(xData_, yData_, xlabel_, ylabel_, outfilenamepdf_, plotcolor_, fillcolor_):

    plt.figure()
    ax = plt.gca() # gca() - get current axes
    plt.plot(xData_, yData_, color=plotcolor_, linewidth=3.0, zorder=3)
    
    for loc, spine in ax.spines.items():  #ax.spines is a dictionary
        spine.set_zorder(10)

    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')

    plt.fill(xData_, yData_, color=fillcolor_) 
    plt.xlabel(r'$' + xlabel_ + '$',fontsize=12)
    plt.ylabel(r'$' + ylabel_ + '$', fontsize=12)
    plt.grid(color='0.85', linestyle='--')
    
    plt.ylim([0,max(yData_)+0.04*max(yData_)])

    plt.savefig(outfilenamepdf_, bbox_inches='tight')



def main():

    os.system('cls' if os.name=='nt' else 'clear')
    
    for Dir in ResultsDir:
        if os.path.exists(Dir):
            for file in os.listdir(Dir):	
                file_path = os.path.join(Dir, file)
                try:
                    if os.path.isfile(file_path):
                        os.unlink(file_path)
                except Exception, e:
                    print e
        else:
            os.makedirs(Dir)

    cpu_cores=multiprocessing.cpu_count()
    start = time.clock()

    pool = Pool(processes=cpu_cores)   
    N = [1,2,3,4,6,8]
    
  
    print "\nExtendable Freely Jointed Chain 3D"
    print "----------------------------------\n"
    print "Number of Links\t\t:", N
    print "Link Length\t\t:", link_length
    print "Extendable Link Length\t:", extendable_length
    print  
    
    pool.map(ZR_data, N)

    N = [3,4,6,8,10,12]
    print "\nExtendable Freely Jointed Chain 3D - Small Approximation"
    print "--------------------------------------------------------\n"
    print "Number of Links\t\t:", N
    print "Link Length\t\t:", link_length
    print "Extendable Link Length\t:", extendable_length
    print

    pool.map(ZR_dataSA, N)

    print 
    print 'Runtime (' + str(cpu_cores) + ' cpu cores) : ' + str(time.clock()-start)
    print
    

if __name__ == '__main__':
    main()

  