import numpy as np
import matplotlib.pyplot as plt
import math
import simulation_data_import as sdi  # Simulation_data_import is where we create functions to import simulation data
                                      # from PSIM´s .txt files. 

psim_txtfile = 'simulation_data3.txt'

Da = sdi.create_signal( sdi.get_column('Da',psim_txtfile) ) # Imports column 'Da' from PSIM´s txt file into a list.
Db = sdi.create_signal( sdi.get_column('Db',psim_txtfile) ) 
Dc = sdi.create_signal( sdi.get_column('Dc',psim_txtfile) )

h=1e-6 # Time step
tMax = 0.08;  # Simulation time.
time = np.arange(0,tMax,h) # Time vector

# GROUP B:----------------------------------------------------------------------------------------------------------------

Dgb = sdi.create_signal( sdi.get_column('Dgb',psim_txtfile) ) # Imports column 'Dgb' from txt file.
 
# Parameters for group B:
Vb=300
Lb=1e-3
Cb=10e-3 
Rb=1


Xb=[]
Ub=[]
Xb.append(np.zeros( (2,1) )) 
Ub.append(np.matrix( [ [Vb] , [0] ])) # Vb is the internal input for group B.

Ib=np.identity(2)

#GROUP T: --------------------------------------------------------------------------------------------------------------

# Parameters for group T:
Lt=3e-3
Rt=5
Ct=400e-6

Xt = []
Ut = []
Xt.append(np.matrix( [ [0],[0],[0],[0] ] )) 
Ut.append( np.matrix(np.zeros( (4,1) ) )) 


#  ea eb and ec are the alternate sources in group T. They are the internal inputs for group T.
freq = 60   
w = 2*np.pi*freq
ea = math.sqrt(2)*110*np.sin(w*time)
eb = math.sqrt(2)*110*np.sin(w*time - 2*np.pi/3)
ec = math.sqrt(2)*110*np.sin(w*time + 2*np.pi/3)

It=np.identity(4) 

# Defining the nodal analysis function: ----------------------------------------


def nodal_analysis(i,A):      # i = Yhistb/Wb + Yhistt/Wt , A= 1/Wt  + 1/Wb,  
    V=np.linalg.inv(A)*i   # i is the vector with the current sources from the Norton equivalent and A the admittance matrix.
    return V

V=[]# All voltages computed in the nodal analysis. 

nPoints = len(time)



for n in range(nPoints):
    
    # We start calcuating W and Yhist for both groups:
    
    #GROUP T
    Sa=Da[n]
    Sb=Db[n]
    Sc=Dc[n]
    
    At = np.matrix([[-Rt/Lt,0,0, (1/Lt)*(2*Sa/3 - Sb/3 -Sc/3)] , [0,-Rt/Lt,0, (1/Lt)*(2*Sb/3 - Sa/3 -Sc/3)] , [0,0,-Rt/Lt, (1/Lt)*(2*Sc/3 - Sb/3 -Sa/3)] , [-Sa/Ct, -Sb/Ct, -Sc/Ct, 0] ])
    Bt = np.matrix([ [-1/Lt,0,0,0], [0,-1/Lt,0,0], [0,0,-1/Lt,0], [0,0,0,-1/Ct] ])
    Atc=np.linalg.inv(It-h*At)
    Btc=Atc*h*Bt                     # Where matrix Mc indicates matrix M in the discretized state space equation.
    Btc_i = Btc[:,0:3]
    Btc_n = Btc[:,3]
    Ctn = np.matrix([[0,0,0,1]]) #  Ct = Ct_n
    
    Wt=Ctn*Btc_n
    Yhistt = Ctn*(Atc*Xt[-1] + Btc_i*Ut[-1][0:3,0]) #Ut_i
        
    #GROUP B
    s = Dgb[n]
    
    Ab=np.matrix([[0 ,(1-s)*(-1/Lb)],[(1-s)*1/Cb, -1/(Rb*Cb)]])
    Bb=np.matrix([[1/Lb,0],[0,(1/Cb-1)*s - 1/Cb]]) 
    Abc=np.linalg.inv(Ib - h*Ab) 
    Bbc=Abc*h*Bb
    Bbc_i=Bbc[:,0]
    Bbc_n=Bbc[:,1]
    Cbn = np.matrix([[0, 1]]) # Cb = Cb_n
    
    Wb = Cbn*Bbc_n 
    Yhistb = Cbn*(Abc*Xb[-1] + Bbc_i*Ub[-1][0,0]) #Ub_i
    #----------------------------------------------------------
    
    #Nodal analysis is computed using Yhistb and Yhistt:
    
    V_new = nodal_analysis( Yhistb/Wb + Yhistt/Wt , 1/Wb + 1/Wt ) # Voltage in the interface Node. 
    V.append( V_new[0,0] )
    
    # Since group B and group T are both I-type groups, U_n holds currents.  We compute the current in the interface
    # of both groups so we can compute vector U. 
    
    # Thevenin equivalent is converted to Norton and the current is calculated with a simple circuit analysis: 
    i = Yhistb/Wb - V_new/Wb
    
    
    #Using the current calculated and the internal inputs we calculate Uk+1.
    Ut.append(np.matrix( [ [ea[n]],[eb[n]],[ec[n]],[i] ])) 
    Ub.append(np.matrix([[Vb],[-i]]))
    
    #Computing the new vectors X k+1.  =>   X k+1 = Ac * Xk  +  Bc * Uk+1
    
    Xt.append(Atc*Xt[-1] + Btc*Ut[-2]) 
    Xb.append(Abc*Xb[-1] + Bbc*Ub[-2]) 
    
figure1=plt.figure(1)
figure1.suptitle('Capacitor Voltages')
subplot1=figure1.add_subplot(111)
subplot1.plot(time, V, label = 'VDC SSN')

VDC = sdi.create_signal( sdi.get_column('VDC',psim_txtfile) ) #Imports column VDC form PSIM´s text file. Interface Node Voltage.
subplot1.plot(time, VDC[:len(time)], label = 'VDC PSIM' )
subplot1.legend()




#figure2=plt.figure(2)
#figure2.suptitle('Ea')
#e=figure2.add_subplot(111)
#e.plot(time, ea)
#e.plot(time, eb)
#e.plot(time, ec)
#e.grid()
#    
#
#figure3=plt.figure(3)
#figure3.suptitle('Duty cicles')
#e2=figure3.add_subplot(111)
#e2.plot(time, Dc[:len(time)])
#e2.plot(time, Db[:len(time)])
#e2.plot(time, Da[:len(time)])
#e2.grid()
#


    
