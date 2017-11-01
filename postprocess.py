# process .xyz file dumped by lammps
# Xia Xiuyang, Beijing Univ. Chem. Tech.
#


# History
    # 2017-10-30,  Xia Xiuyang: originatoml version

# To Do List


#Sample


# Variables



#import
import numpy as np
from scipy.spatial import distance
#Mean Square Displacement
def Msd(atom_type):

    msd_coord=Whichatoms(atom_type,coord,atominf)
    msd_coord=Coord_image(msd_coord, boxinf)
    Xcoord = msd_coord - np.tile(msd_coord[0],(msd_coord.shape[0],1,1))
    Xdistance = np.add.reduce(np.square(Xcoord), 2)
    msd = np.add.reduce(Xdistance, 1) / msd_coord.shape[1]
    Print('msd.txt',time,msd)



#Radial distribution function
def Rdf(atom_type1,atom_type2,time_of_snapshot,nhist=100):
    rdf_atominf, rdf_coord, rdf_boxinf = Whichsnapshot(time_of_snapshot, atominf, coord, boxinf)
    rdf_coord1=Coord_unimage(Whichatoms(atom_type1,rdf_coord,rdf_atominf))
    rdf_coord2=Coord_unimage(Whichatoms(atom_type2,rdf_coord,rdf_atominf))
    rdf_delta=np.stack([distance.cdist(rdf_coord1[:,i].reshape((-1,1)),rdf_coord2[:,i].reshape((-1,1)),'euclidean') for i in range(3)])
    length=rdf_boxinf.mean()
    rdf_delta[rdf_delta>length/2]=length-rdf_delta[rdf_delta>length/2]
    rdf_distance=np.sqrt(np.add.reduce(np.square(rdf_delta),0))
    thickness=length/2/nhist
    rho=rdf_coord2.shape[0]/(length**3)
    rdf,ni=[],[]
    r=np.arange(0,length/2,thickness)
    for hist in range(nhist):
        if hist==0:
            rdf.append(0)
        else:
            n=(rdf_distance[rdf_distance>hist*thickness])[(rdf_distance[rdf_distance>hist*thickness])<(hist+1)*thickness].shape[0]/rdf_coord1.shape[0]
            rdf.append(n/(4*rho*np.pi*((hist*thickness)**2)*thickness))
    Print('rdf.txt',r,rdf)
    return length,r,rdf,rho

#Static structure factor
def Ssf(atom_type1,atom_type2,time_of_snapshot,nhist=100):
    length,r,rdf,rho=Rdf(atom_type1,atom_type2,time_of_snapshot,nhist)
    k=np.linspace(0,np.pi/length,nhist,endpoint=False)
    ssf=[]
    k=np.delete(k,0)
    r=np.delete(r,0)
    rdf=np.delete(np.array(rdf),0)
    for ki in k:
        ssf.append(1+4*np.pi*rho*np.add.reduce((r**2)*rdf*np.sin(ki*r)/(ki*r)*(length/2/nhist)))
    Print('ssf.txt',k,ssf)

#import .xyz file
def Importfile(filename,starttime,endtime,timestep):
    global time,atominf,coord,boxinf,natom
    nsnapshot=int((endtime-starttime)/timestep+1)
    with open(filename,'r') as dump:
        time,bx,by,bz,atominf,mid,coord,allinf,boxinf=[],[],[],[],[],[],[],[],[]
        for snapshot in range(nsnapshot):
            dump.readline()
            time.append(dump.readline())
            dump.readline()
            natom=int(dump.readline())
            dump.readline()
            bx.append(list(dump.readline().split()))
            by.append(list(dump.readline().split()))
            bz.append(list(dump.readline().split()))
            dump.readline()
            mid=[list(dump.readline().split()) for line in range(natom)]
            allinf.append(mid)
            mid=[]
    time=np.array(time).astype(int)
    atominf=np.array(allinf)[:,:,:2].astype(int)
    coord=np.array(allinf)[:,:,2:].astype(float)
    bx,by,bz=np.array(bx).astype(float),np.array(by).astype(float),np.array(bz).astype(float)
    boxinf=np.column_stack((bx[:,1]-bx[:,0],by[:,1]-by[:,0],bz[:,1]-bz[:,0]))
    return time,atominf,coord,boxinf,natom




def Whichsnapshot(time_of_snapshot,pre_atominf,pre_coord,pre_boxinf):
    if time_of_snapshot=='end':
        one_atominf=atominf[-1]
        one_coord=pre_coord[-1]
        one_boxinf=boxinf[-1]
    elif isinstance(time_of_snapshot,int) is True:
        one_atominf=atominf[time==time_of_snapshot]
        one_coord=pre_coord[time==time_of_snapshot]
        one_boxinf=boxinf[time==time_of_snapshot]
    return one_atominf,one_coord,one_boxinf




def Whichatoms(atom_type,pre_coord,atominf):
    if pre_coord.ndim==3:
        if atom_type=='all':
            new_coord=pre_coord
        elif isinstance(atom_type,list) is True:
            new_coord=np.hstack([pre_coord[atominf[:,:,1]==type].reshape((pre_coord.shape[0],-1,6)) for type in atom_type])
        elif isinstance(atom_type,int) is True:
            new_coord=pre_coord[atominf[:,:,1]==atom_type].reshape((pre_coord.shape[0],-1,6))
        return new_coord
    elif pre_coord.ndim==2:
        if atom_type=='all':
            new_coord=pre_coord
        elif isinstance(atom_type,list) is True:
            new_coord=np.row_stack([pre_coord[atominf[:,1]==type].reshape((-1,6)) for type in atom_type])
        elif isinstance(atom_type,int) is True:
            new_coord=pre_coord[atominf[:,1]==atom_type]
        return new_coord


def Coord_image(pre_coord, boxinf):
    if pre_coord.ndim == 3:
        boxinf = np.tile(boxinf.reshape((boxinf.shape[0], 1, 3)), (1, pre_coord.shape[1], 1))
        coord_image = pre_coord[:, :, 0:3] + pre_coord[:, :, 3:6] * boxinf
        return coord_image
    elif pre_coord.ndim == 2:
        boxinf = np.tile(boxinf, (pre_coord.shape[1], 1))
        coord_image = pre_coord[:, 0:3] + pre_coord[:, 3:6] * boxinf
        return coord_image

def Coord_unimage(pre_coord):
    if pre_coord.ndim == 3:
        coord_unimage = pre_coord[:, :, 0:3]
        return coord_unimage
    elif pre_coord.ndim == 2:
        coord_unimage = pre_coord[:, 0:3]
        return coord_unimage


def Print(filename, a, b):
    with open(filename, 'w') as file:
        for i in range(a.shape[0]):
            file.write('%g\t%-10g\n' % (a[i], b[i]))
    return 0

#以下是测试
Importfile('as.lammpstrj', 110000, 1100000, 10000)
#print(time.shape)
#print(atominf.shape)
#print(coord.shape)
#print(boxinf.shape)
#
#coord_image=Coord_image(coord,boxinf)
#print(coord_image.shape)
#
#coord_unimage=Coord_unimage(coord,boxinf)
#print(coord_unimage.shape)
#one_coord,one_boxinf=Whichsnapshot('end')
#print(one_coord.shape)
#print(one_boxinf.shape)



Ssf(2,2,'end',50)









#output .txt file
#output .pkl file