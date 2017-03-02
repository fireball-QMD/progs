#!/usr/bin/env python
"""
A server to control thermodynamic integration option of fireball.  This uses one massive nonequilibrium simulation.
"""
# user-defined vars---------------------------
numAtoms1=182 #254
numAtoms2=174 #250
numAtoms3=4
dir1='whole'
dir2='vacancy'
dir3='molecules'

# code---------------------------------------

import select,socket,sys,os,popen2
from my import *

timesteps=int(sys.argv[1]) # these are one way only!
forward=[x/float(timesteps) for x in range(timesteps)]
reverse=[x/float(timesteps) for x in range(timesteps)]
reverse.reverse()
lam=forward+[1]+reverse
energies=[]

# server initialization
host = 'localhost'
port = 900001 # hopefully this isn't used
backlog = 5
size = 100
server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
server.bind((host,port))
server.listen(backlog)

# start fireball
os.popen(dir1+'/cleanup.com 2> junk.dat')
os.popen(dir2+'/cleanup.com 2> junk.dat')
os.popen(dir3+'/cleanup.com 2> junk.dat')
os.chdir(dir1)
fire1=popen2.Popen4('./fireball.x > f'+str(timesteps)+'.log')
os.chdir('../'+dir2)
fire2=popen2.Popen4('./fireball.x > f'+str(timesteps)+'.log')
os.chdir('../'+dir3)
fire3=popen2.Popen4('./fireball.x > f'+str(timesteps)+'.log')
os.chdir('..')

# open energy file
f_en=file('energies_'+str(timesteps)+'.txt','w')

def lamFunc(tau):
 return 1-tau**5*(70*tau**4-315*tau**3+540*tau**2-420*tau+126)

# begin loop
for i in ind(lam):
    print 'lam is',lam[i]
    dataNum=0
    clientData=[0.,0.,0.]   
    input = [server]
    running=True
    while running:
        inputready,outputready,exceptready = select.select(input,[],[])
        for s in inputready:
            if s == server:
                # handle the server socket
                client, address = server.accept()
                input.append(client)
            else:
                # handle all other sockets
                data = s.recv(size)
                if data:
                    # print data
                    clientData[dataNum]=data
                    dataNum+=1
                    if dataNum>2:
                        # print 'data recieved is',clientData[0],'and',clientData[1]
                        returnSignal='done'
                        # read the energies of the current time step from both processes
                        # and average them
                        f1=file(dir1+'/energyValues.txt','r');
                        en1=float(f1.readline())
                        f2=file(dir2+'/energyValues.txt','r');
                        en2=float(f2.readline())
                        f3=file(dir3+'/energyValues.txt','r');
                        en3=float(f3.readline())
                        # print 'en1 is',en1
                        # print 'en2 is',en2
                        # print 'en3 is',en3
                        tot_en=-en1 + en2 + en3
                        print 'tot_en is',tot_en
                        print >>f_en, tot_en
                        f_en.flush()
                        energies.append(tot_en)
                        f1.close();f2.close();f3.close()
                        # read the forces and average them
                        # system 1
                        f1=file(dir1+'/forceValues.txt','r')
                        ftot1=[]
                        for line in f1:
                            ftot1.append([float(coor) for coor in line.split()])
                        # print 'ftot1',ftot1
                        # system 2 -- pad with three 0's so can add to system 1
                        f2=file(dir2+'/forceValues.txt','r')
                        ftot2=[]
                        for line in f2:
                            ftot2.append([float(coor) for coor in line.split()])
                        # print 'ftot2',ftot2
                        ftot2.append([0,0,0])
                        ftot2.append([0,0,0])
                        ftot2.append([0,0,0])
                        ftot2.append([0,0,0])
                        # system 3 -- fill with 0's before reading file so can add to system 1
                        f3=file(dir3+'/forceValues.txt','r')
                        ftot3=[[0,0,0] for iii in range(numAtoms2)]
                        for line in f3:
                            ftot3.append([float(coor) for coor in line.split()])
                        # print 'ftot3',ftot3
                        # add them all together
                        ftot=[[(1-lamFunc(lam[i]))*ftot1[ii][jj] + lamFunc(lam[i])*ftot2[ii][jj] \
                         + lamFunc(lam[i])*ftot3[ii][jj] for jj in range(3)] for ii in ind(ftot1)] 
                        f1.close();f2.close();f3.close()
                        
                        # write new forces back out
                        nf1=open(dir1+'/newforce.txt','w')
                        nf2=open(dir2+'/newforce.txt','w')
                        nf3=open(dir3+'/newforce.txt','w')
                        # system 1
                        for j in range(numAtoms1):
                            print >>nf1, ftot[j][0],ftot[j][1],ftot[j][2]
                        # system 2
                        for j in range(numAtoms2): 
                            print >>nf2, ftot[j][0],ftot[j][1],ftot[j][2]
                        # system 3
                        for j in range(numAtoms3):
                            print >>nf3, ftot[j+numAtoms2][0],ftot[j+numAtoms2][1],ftot[j+numAtoms2][2]
                        nf1.close();nf2.close();nf3.close()
                        # send signal to clients to continue
                        for s in input:
                            if s == server:
                                pass
                            else:
                                s.send(returnSignal)
                        dataNum=0
                        clientData=[0.,0.,0.]
                        running=False  # get out of loop to process data
                    break
#    os.system('cp')
server.close()
#print energies
## energiesStr=[str(energies[i])+'\n' for i in ind(energies)]
## f_en_back=file('energies.txt','w')
## f_en_back.writelines(energiesStr)
