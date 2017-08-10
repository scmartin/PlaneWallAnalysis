import matplotlib

matplotlib.use('tkagg')

import argparse
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description="Vxs error plotting")
parser.add_argument('xmin',type=float)
parser.add_argument('xmax',type=float)
parser.add_argument('ymin',type=float)
parser.add_argument('ymax',type=float)
args = parser.parse_args()
xmin = args.xmin
xmax = args.xmax
ymin = args.ymin
ymax = args.ymax
vxs = np.fromfile("vxs.txt",sep=" ")
vxs.shape = (int(len(vxs)/4),4)
#print(vxs)

plt.plot(vxs[:,0],vxs[:,1],vxs[:,0],vxs[:,1]+vxs[:,2],vxs[:,0],vxs[:,1]-vxs[:,2])
plt.axis([xmin,xmax,ymin,ymax])
plt.show()
