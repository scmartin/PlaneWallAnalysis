"""
This program performs error analysis for the planar hard wall system.
It takes in the density profile trajectory and calculates rho_i and rho_bulk
and the associated error for all choices of interface volume cutoff. Excess
volume is calculated for each of these interface volumes, and the error is
propagated through using either the error in rho_bulk or rho_i
"""
import numpy as np
import argparse

""" Takes in command line arguments for the name of the density trajectory,
the distance between hard walls, the length of the hard walls, and the number
of histogram bins in the density profiles, respectively """

parser = argparse.ArgumentParser(description="VXS error propagation")
parser.add_argument('infile',type=str)
parser.add_argument('ydim',type=float)
parser.add_argument('xdim',type=float)
parser.add_argument('nbin',type=int)
args = parser.parse_args()
ydim = args.ydim
xdim = args.xdim
nbin = args.nbin

vol = ydim*xdim  # volume of the simulation box
rho = np.fromfile(args.infile,sep=" ")
steps = int(len(rho)/nbin) # number of profiles in the trajectory
rho.shape = (steps,nbin)
#print(rho[0:5,0:10])
rhotot = np.mean(rho)
N = rhotot*vol
ybin = ydim/float(nbin) # width of each histogram bin
z = np.arange(ybin,ydim,ybin,float) # interface cutoff values
#print(z)
dfile = open("./dens.txt","w")  # Holds the rho_i,rho_b, and the associated errors
vfile = open("./vxs.txt","w")  # Holds the excess volumes and the associated error
                               #   propagated from both rho_i and rho_b
""" The following three arrays hold block standard error estimate data for rho_i of each wall
and rho_b for each choice of interface volume. Each row contains estimates for 10:200
blocks for a given choice of interface volume
"""
ri1err = np.zeros((int(nbin/2)-5,15))  # error in rho_i on one side of sim box
ri2err = np.zeros((int(nbin/2)-5,15))  # error in rho_i of other hard wall
rberr = np.zeros((int(nbin/2)-5,15))   # error in rho_b
for i in range(int(nbin/2)-5):  # loops over choices of interface volumes
    yi = ybin*float(i+1)  # distance the interface extends away from the wall
#    efname = "./error"+str(yi)+".txt"
#    errfile = open(efname,"w")
    rhoi = (np.mean(rho[:,:i+1]) + np.mean(rho[:,nbin-i-1:nbin:1]))/(2.0) # total interface density 
    rhoi1 = np.mean(rho[:,:i+1])  # interface density for one wall
    rhoi2 = np.mean(rho[:,nbin-i-1:nbin:1]) # interface density for other wall
    rhob = (np.mean(rho[:,i+1:nbin-i-1])) # bulk density
    for l in range(15):  # block averaging loop
        block = 10+5*(l)
        m = int(steps/block)
        bdens = np.zeros(block) 
        f1dens = np.zeros(block)
        f2dens = np.zeros(block)
        for j in range(block):
            f1dens[j] = np.mean(rho[j*m:(j+1)*m,:i+1])
            f2dens[j] = np.mean(rho[j*m:(j+1)*m,nbin-i-1:nbin])
            bdens[j] = np.mean(rho[j*m:(j+1)*m,i+1:nbin-i-1])
        """ this section calculates the block standard error for the density of each interface independendently,
            as well as the BSE of the bulk density """
        varbdens = ((np.sum((bdens-rhob)**2)))/(float(block-1))  # variance of the blocks
        varf1dens = ((np.sum((f1dens-rhoi1)**2)))/(float(block-1)) 
        varf2dens = ((np.sum((f2dens-rhoi2)**2)))/(float(block-1))
        varbd_mean = varbdens/float(block)  # variance in the mean
        varf1d_mean = varf1dens/float(block)
        varf2d_mean = varf2dens/float(block)
        bdensbse = np.sqrt(varbd_mean)  # block standard error
        f1densbse = np.sqrt(varf1d_mean)
        f2densbse = np.sqrt(varf2d_mean)
        ri1err[i,l] = f1densbse  
        ri2err[i,l] = f2densbse
        rberr[i,l] = bdensbse
    """ using the error in the two interfaces, we can estimate the error in the interface density as BSE(rho_i)/sqrt(2)
        because the interfaces are independent. Therefore, I take the mean of the two estimates ri1err and ri2err as the 
        BSE of rho_i """
    rierr = (np.mean(ri1err[i,:])+np.mean(ri2err[i,:]))/(2.0*np.sqrt(2.0)) 
    rberr2 = (np.mean(rberr[i,:]))
    vxs = yi*(1.0-rhoi/rhob) # excess volume
    vi2err = rierr*yi*N/(rhob**2*(vol-yi*2.0*xdim)) # error in excess volume propagated from error in rho_i
    vberr = np.mean(rberr[i,:])*yi*N/(rhob**2*(yi*2.0*xdim)) # error in the excess volume propagated from the error in rho_b
    vfile.write("{:12.9f}".format(z[i]))
    vfile.write("\t")
    vfile.write("{:12.9f}".format(vxs))
    vfile.write("\t")
    vfile.write("{:12.9f}".format(vi2err))
    vfile.write("\t")
    vfile.write("{:12.9f}".format(vberr))
    vfile.write("\n")
    dfile.write("{:12.9f}".format(z[i]))
    dfile.write("\t")
    dfile.write("{:12.9f}	{:12.9f}".format(rhoi,rierr))
    dfile.write("\t")
    dfile.write("{:12.9f}	{:12.9f}".format(rhob,rberr2))
    dfile.write("\n")
vfile.close()
dfile.close()
np.savetxt("ri1err.txt",ri1err) #error of one interface
np.savetxt("ri2err.txt",ri2err) #error of the second interface
np.savetxt("rberr.txt",rberr) #BSE file
