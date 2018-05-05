'''
  Copyright 2018 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

   OPM is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with OPM.  If not, see <http://www.gnu.org/licenses/>.
'''

from Tkinter import *
import tkFileDialog
import ttk
import os
import sys
import matplotlib
import math
import imageio
matplotlib.use('TkAgg')


from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# implement the default mpl key bindings
from matplotlib.figure import Figure
import pydicom
import pylab
import numpy as np
from os import listdir
from os.path import join
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import csv

# Linear interpolation
def lerp(a, b, f): return (1 - f) * a + b * f
#Functions for getting voxels from a 3D texture
# texture - 3D array of values (3D texture)
# x, y, z - texture coordinates (from 0 to 1)
def get_pixel_simple(texture, x, y, z):
    # simple, no interpolation
    zmax = len(texture) - 1
    ymax = len(texture[0]) - 1
    xmax = len(texture[1]) - 1
    vZ = int(z * zmax)
    vY = int(y * ymax)
    vX = int(x * xmax)
    return texture[vZ][vY][vX]
def get_pixel_trilinear(texture, x, y, z):
    # triliniear interpolation
    zmax = len(texture) - 1
    ymax = len(texture[0]) - 1
    xmax = len(texture[1]) - 1
    vZ = z * zmax
    vY = y * ymax
    vX = x * xmax
    z0 = min(int(math.floor(vZ)), zmax - 1)
    z1 = z0 + 1
    y0 = min(int(math.floor(vY)), ymax - 1)
    y1 = y0 + 1
    x0 = min(int(math.floor(vX)), xmax - 1)
    x1 = x0 + 1
    xd = (vX - x0) / (x1 - x0)
    yd = (vY - y0) / (y1 - y0)
    zd = (vZ - z0) / (z1 - z0)
    c000 = texture[z0][y0][x0]
    c100 = texture[z0][y0][x1]
    c010 = texture[z0][y1][x0]
    c001 = texture[z1][y0][x0]
    c110 = texture[z0][y1][x1]
    c101 = texture[z1][y0][x1]
    c011 = texture[z1][y1][x0]
    c111 = texture[z1][y1][x1]
    c00 = lerp(c000, c100, xd)
    c01 = lerp(c001, c101, xd)
    c10 = lerp(c010, c110, xd)
    c11 = lerp(c011, c111, xd)
    c0 = lerp(c00, c10, yd)
    c1 = lerp(c01, c11, yd)
    c = lerp(c0, c1, zd)
    return c
def get_pixel_coordinates(picX, picY, imgWidth, imgHeight, radius = 1.0, offset = (0, 0), startAngle = 0.0):
    # this function maps x,y of a image to a cylinder coordinate in 3d texture 
    # center is at (0.5, 0.5, 0.5), all coordinates are between 0 and 1 (start and end of the texture)
    angleStep = 2 * math.pi / imgWidth
    volX = (0.5 + radius * 0.5 * math.cos(picX * angleStep + startAngle)) + offset[0]
    volY = (0.5 + radius * 0.5 * math.sin(picX * angleStep + startAngle)) + offset[1]
    volZ = 1.0 * picY / imgHeight
    return (volX, volY, volZ)
def create_cylinder(imageWidth, imageHeight, voxels, radius = 1.0, offset = (0, 0), startAngle = 0.0):
    imageData = np.zeros([imageHeight, imageWidth])
    for y in range(0, imageHeight):
        for x in range(0, imageWidth):
            voxelPos = get_pixel_coordinates(x, y, imageWidth, imageHeight, radius, offset, startAngle)
            imageData[y][x] = get_pixel_trilinear(voxels, voxelPos[0], voxelPos[1], voxelPos[2])
    return imageData


#Define a class to return the values slice by slice from the dicom

class Dicom:
    def __init__(self,Path,display_img,Pixel):

        # Initialisation of the arrays

        sliceList=[]
        middleslicex=[]
        middleslicey=[]
        voxels = []
        x=[]
        i=Padding_bottom

        # loop for each .dcm file
        length=len(listdir(Path))-1
        files=[f for f in listdir(Path)]
        files=sorted(files)

        pb_hD.maximum=length

        for i,f in enumerate(files):


            if StopRun ==True: break
            if (i<=length-Padding_top and i>=Padding_bottom):
                if int(i*100/length)-int((i-1)*100/length)==1:
                    pb_hD.step(1)
                    pb_hD.update()

                ds = pydicom.read_file(join(Path, f))
                x=np.array(ds.pixel_array)
                voxels.append(x)

                if (i==Padding_bottom):
                    a = x.shape[0]/2
                    Offsetr=2*a*Offsetx/Diameter
                    Offsetc=2*a*Offsety/Diameter

                # Filter the values

                x_masked=GetMaskedValues(x,Offsetr,Offsetc)

                sliceList.append(x_masked)
                middleslicex.insert(0,x[:,int(a+Offsetc)])
                middleslicey.insert(0,x[int(a+Offsetr),:])

        #gather slices into a single np-matrix
        self.CTPixel=np.zeros([len(sliceList),len(sliceList[0])])
        self.CTSlice=np.zeros([len(middleslicex),len(middleslicex[0])])
        self.CTSlicey=np.zeros([len(middleslicey),len(middleslicey[0])])
        self.CTVoxels=np.array(tuple(voxels))
        self.CTCylinder=create_cylinder(
            int(math.pi * len(self.CTVoxels[0]) * (1.0 - Crop_pct)),
            len(self.CTVoxels), self.CTVoxels, 1.0 - Crop_pct, (Offsetx/Diameter, Offsety/Diameter))

        for j in range(0,len(sliceList)-1):
            self.CTPixel[j]=sliceList[j]
            self.CTSlice[j]=middleslicex[j]
            self.CTSlicey[j]=middleslicey[j]

        pb_hD.stop()

        pathTuple = os.path.split(Path)
        imageio.imsave(os.path.join(pathTuple[0], pathTuple[1] + "_cylinder.jpg"), self.CTCylinder)
        imageio.imsave(os.path.join(pathTuple[0], pathTuple[1] + "_xz.jpg"), self.CTSlice)
        imageio.imsave(os.path.join(pathTuple[0], pathTuple[1] + "_yz.jpg"), self.CTSlicey)
        imageio.imsave(os.path.join(pathTuple[0], pathTuple[1] + "_xy.jpg"), ds.pixel_array)


    # Show the image

        if (display_img):

            fig1 = plt.figure()
            ax0 = fig1.add_subplot(234,aspect='equal')
            ax0.set_title("Edge Slice")
            ax0.set_xlabel('Angle')
            ax0.set_ylabel('Z')
            ax0.matshow(self.CTCylinder, cmap=plt.cm.gray,aspect='equal')

            ax1 = fig1.add_subplot(232,aspect='equal')
            ax1.set_title("XZ Slice")
            ax1.set_xlabel('X')
            ax1.set_ylabel('Z')
            ax1.matshow(self.CTSlice, cmap=plt.cm.gray,aspect='equal')
            ax1.add_patch(patches.Rectangle((a-Offsetr-a*Crop_pct, 0), 2*a*Crop_pct,    length, fill=False,edgecolor="red"))

            ax3 = fig1.add_subplot(233,aspect='equal')
            ax3.set_title("YZ Slice")
            ax3.set_xlabel('Y')
            ax3.set_ylabel('Z')
            ax3.matshow(self.CTSlicey, cmap=plt.cm.gray,aspect='equal')
            ax3.add_patch(patches.Rectangle((a-Offsetc-a*Crop_pct, 0), 2*a*Crop_pct,    length, fill=False,edgecolor="red"))

            ax2 = fig1.add_subplot(231,aspect='equal')
            ax2.set_title("XY Slice")
            ax2.set_xlabel('X')
            ax2.set_ylabel('Y')
            ax2.imshow(ds.pixel_array, cmap=pylab.cm.bone)
            circle=plt.Circle((a+Offsetr,a+Offsetc),a*Crop_pct,color='r',linewidth=1,fill=False)
            plt.gcf().gca().add_artist(circle)
            ax2.plot((a-10+Offsetr , a+10+Offsetr), (a+Offsetc, a+Offsetc), 'k')
            ax2.plot((a+Offsetr, a+Offsetr),(a-10+Offsetc , a+10+Offsetc), 'k')
            plt.tight_layout()

            canvas = FigureCanvasTkAgg(fig1, master=stepFour)
            canvas.show()
            canvas.get_tk_widget().grid(row=0, column=0, columnspan=3, pady=2, sticky='W')

            form.update()






# This function is used to return the min/max for the filtered circle area

def GetMaskedValues(x,Offsetr,Offsetc):

    # A circle shape is generated based on the diameter and the x/y offsets

    n=x.shape[0]

    a=n/2+Offsetc
    b=n/2+Offsetr
    r=n/2*Crop_pct
    ny,nx = np.ogrid[-a:n-a, -b:n-b]

    # The mask will apply to the values inside the circle

    mask = nx*nx + ny*ny <= r*r

    return x[mask]


# This function is used for the plotting

def Plot(lowenergy,a,b,c,depth):


    fig2=plt.figure(1)
    plt.subplot(141)
    plt.subplot(141).set_title("Density", fontweight='bold')
    plt.plot(a,depth,'b')
    plt.subplot(141).ticklabel_format(useOffset=False)


    plt.subplot(142)
    plt.subplot(142).set_title("Photoelectric factor", fontweight='bold')
    plt.plot(b,depth,'b')
    plt.subplot(142).ticklabel_format(useOffset=False)


    plt.subplot(143)
    plt.subplot(143).set_title("Porosity", fontweight='bold')
    plt.plot(c,depth,'b')
    plt.subplot(143).ticklabel_format(useOffset=False)

    plt.subplot(144).matshow(lowenergy.CTSlice, cmap=plt.cm.gray,aspect='auto')
    plt.subplot(144).axis('off')
    plt.tight_layout()

    canvas = FigureCanvasTkAgg(fig2, master=stepFour)
    canvas.show()
    canvas.get_tk_widget().grid(row=0, column=0, columnspan=3, pady=2, sticky='W')
    form.update()



def Plot2(lowenergy,a,b,c,mslices):

    fig2=plt.figure(1)
    plt.subplot(141)
    plt.subplot(141).set_title("Density", fontweight='bold')
    plt.plot(a,mslices,'b')
    plt.subplot(141).ticklabel_format(useOffset=False)


    plt.subplot(142)
    plt.subplot(142).set_title("Photoelectric factor", fontweight='bold')
    plt.plot(b,mslices,'b')
    plt.subplot(142).ticklabel_format(useOffset=False)


    plt.subplot(143)
    plt.subplot(143).set_title("Porosity", fontweight='bold')
    plt.plot(c,mslices,'b')
    plt.subplot(143).ticklabel_format(useOffset=False)

    plt.subplot(144).matshow(lowenergy.CTSlice, cmap=plt.cm.gray,aspect='auto')
    plt.subplot(144).axis('off')
    plt.tight_layout()

    canvas = FigureCanvasTkAgg(fig2, master=stepFour)
    canvas.show()
    canvas.get_tk_widget().grid(row=0, column=0, columnspan=3, pady=2, sticky='W')
    form.update()



# Function to check the validity of the inputs

def CheckInput():

    return True

def Calculate_Parameters(low,high,Mask_saturated):

    low.m=low.CTPixel #-1024 #uncomment -1024 when CT values are with offset
    high.m=high.CTPixel #-1024 #uncomment -1024 when CT values are with offset
    nnslices=len(low.m)
    Tbox.delete('1.0', END)
    Tbox.insert(END, "Low Energy CT mean value: "+str(low.m.mean())+"\n")
    Tbox.insert(END, "High Energy CT mean value: "+str(high.m.mean())+"\n")
    Tbox.insert(END, "Low Energy CT max value: "+str(low.m.max())+"\n")
    Tbox.insert(END, "High Energy CT max value: "+str(high.m.max())+"\n")


    # Calculate parameters
    if Mask_saturated==1:
        try:
            masked_low = np.ma.masked_where(low.m>3070,low.m)
            masked_high = np.ma.masked_where(low.m>3070,high.m)
            print "Applied the mask on saturated pixels"
        except:
            print "Did not manage to apply the mask on saturated pixels"
    else:
        masked_low=low.m
        masked_high=high.m

    masked_low_slices=np.zeros(nnslices)
    for i in range(0,nnslices-1):
        masked_low_slices[i]=masked_low[i].mean()
    #for i in range(0,np.shape(masked_high)[0]-1):
    masked_high_slices=np.zeros(nnslices)
    for i in range(0,nnslices-1):
        masked_high_slices[i]=masked_high[i].mean()

    Density=fM*masked_low_slices+fP*masked_high_slices+fQ
    Zeff=((fR*masked_low_slices+fS*masked_high_slices+fT)/(0.9342*(fM*masked_low_slices+fP*masked_high_slices+fQ)+0.1759))**(1.0/3.6)
    Phi=(rho_matrix-Density)/(rho_matrix-rho_fluid)
    AvgLowCT=masked_low_slices
    AvgHighCT=masked_high_slices
    Pe=(Zeff/10)**3.6

    return Density,Pe,Phi,nnslices,AvgLowCT,AvgHighCT,Zeff

def create_csv(path,depth,a,b,c,d,e,f,maincsv):

    with open(path+'.csv', 'wb') as myfile:
        csv.register_dialect('custom', delimiter=';')
        wr = csv.writer(myfile,dialect='custom')
        wr2 = csv.writer(maincsv,dialect='custom')
        rows = zip(depth,a,b,c,d,e,f)
        wr.writerow(["--depth","Density","Zeff","Pe","Phi","AvgLowCT","AvgHighCT"])
        wr2.writerow(["--depth","Density","Zeff","Pe","Phi","AvgLowCT","AvgHighCT"])
        for row in rows:
            wr.writerow(row)
            wr2.writerow(row)

#Check Inputs
if not CheckInput():
    print "Missing inputs"
    sys.exit(0)


def Calc_for_well(myfolder,maincsv):

    for subfolder in os.listdir(myfolder):
        path=os.path.join(myfolder,subfolder);print path
    #look for high and low energy folders (here 140kV and 80kV)
        if (os.path.exists(path+'/140kV') and os.path.exists(path+'/80kV')):
            split=myfolder.split("_")
            split2=split[-1].split("-")
            bottom=float(split2[0])
            top=float(split2[1])
            count=len([name for name in os.listdir(path+'/140kV') ])-0.5-Padding_top-Padding_bottom
            step=(top-bottom)/count
            depth=np.arange(bottom,top,step)

            for name in os.listdir(path):

                if (name=="140kV"):
                    highenergy=Dicom(os.path.join(path,name),Check,Mask_saturated)
                if (name=="80kV"):
                    lowenergy=Dicom(os.path.join(path,name),Check,Mask_saturated)


            (Density,Pe,Phi,nnslices,AvgHighCT,AvgLowCT,Zeff)=Calculate_Parameters(lowenergy,highenergy,Mask_saturated)

            Tbox.insert(END, "Density mean value: "+str(np.mean(Density))+"\n")
            Tbox.insert(END, "Photoelectric factor mean value: "+str(np.mean(Pe))+"\n");Tbox.insert(END, "Atomic number mean value: "+str(np.mean(Zeff))+"\n")
            Tbox.insert(END, "Porosity mean value: "+str(np.mean(Phi))+"\n")
            create_csv(path+"/Output_WFL",depth,Density,Zeff,Pe,Phi,AvgLowCT,AvgHighCT,maincsv)
            if Show_Plot:
                Plot(lowenergy,Density,Pe,Phi,depth)
        else:
            print "Skipped Folder "+myfolder

def Main():
    if Browse:
        rootdir=dirnameroot

        maincsv=tkFileDialog.asksaveasfile(title='export results to file', mode='w', defaultextension=".csv")
        for folder in os.listdir(rootdir):
            if os.path.isdir(os.path.join(rootdir, folder)) :
                Calc_for_well(os.path.join(rootdir, folder),maincsv)
                plt.close('all')
    else:
        rootdir=dirnamelow
        rootdir2=dirnamehigh
        highenergy=Dicom(rootdir,Check,Mask_saturated)
        lowenergy=Dicom(rootdir2,Check,Mask_saturated)

        (Density,Pe,Phi,nnslices,AvgLowCT,AvgHighCT,Zeff)=Calculate_Parameters(lowenergy,highenergy,Mask_saturated)
        Tbox.insert(END, "Density mean value: "+str(np.mean(Density))+"\n")
        Tbox.insert(END, "Photoelectric factor mean value: "+str(np.mean(Pe))+"\n")
        Tbox.insert(END, "Atomic number mean value: "+str(np.mean(Zeff))+"\n")
        Tbox.insert(END, "Porosity mean value: "+str(np.mean(Phi))+"\n")
        path=os.path.abspath(os.path.join(rootdir, os.pardir))
        mslices=np.arange(0,nnslices,1)
        if Show_Plot:
            Plot2(lowenergy,Density,Pe,Phi,mslices)



        sys.stdout.write("\n Script finished \n Skipped+Skippedfolders+Folders")

if __name__ == '__main__':
    form = Tk()

    getFld = IntVar()

    form.wm_title('OpenDECT')

    stepZero = LabelFrame(form, text=" 1. Select DECT Files: ")
    stepZero.grid(row=0, columnspan=7, sticky='W', \
                 padx=8, pady=8, ipadx=8, ipady=8)

    stepOne = LabelFrame(form, text=" 2. Enter DECT Parameters from Calibration: ")
    stepOne.grid(row=1, columnspan=7, sticky='W', \
                 padx=8, pady=8, ipadx=8, ipady=8)

    stepTwo = LabelFrame(form, text=" 3. Cropping Filter: ")
    stepTwo.grid(row=2, columnspan=7, sticky='W', \
                 padx=5, pady=5, ipadx=5, ipady=5)

    stepThree = LabelFrame(form, text=" 4. Calculation Parameters: ")
    stepThree.grid(row=3, columnspan=7, sticky='W', \
                   padx=5, pady=5, ipadx=5, ipady=5)

    stepFour = LabelFrame(form, text=" 5. Crop Inspection ")
    stepFour.grid(row=0, column=9, columnspan=7, rowspan=3, \
                  sticky='W', ipadx=5, ipady=5)

    stepFive = LabelFrame(form, text=" 6. Info : ")
    stepFive.grid(row=3,column=9,columnspan=7,  rowspan=3,sticky='W', \
                ipadx=5, ipady=5)

############
#Files Pane#
############
    var_TwoFoldChk = IntVar()
    var_RootFoldChk= IntVar()
    TwoFoldChk = Checkbutton(stepZero, \
                            text="Single sample ", \
                            onvalue=1, offvalue=0, variable=var_TwoFoldChk)
    TwoFoldChk.grid(row=0, column=0)

    def LoadDECTLow():
        global dirnamelow
        dirnamelow = tkFileDialog.askdirectory(parent=stepZero, initialdir='/home/', \
                                               title="Select Low Energy DECT folder")
    def LoadDECTHigh():
        global dirnamehigh
        dirnamehigh = tkFileDialog.askdirectory(parent=stepZero, initialdir='/home/', \
                                                title="Select High Energy DECT folder")


    loadlownrjbtn = Button(stepZero, text="Choose Low Energy Directory", command=LoadDECTLow)
    loadlownrjbtn.grid(row=1, column=0, sticky='W', padx=5, pady=2)
    loadhighnrjbtn = Button(stepZero, text="Choose High Energy Directory", command=LoadDECTHigh)
    loadhighnrjbtn.grid(row=2, column=0, sticky='W', padx=5, pady=2)

    RootFoldChk = Checkbutton(stepZero, \
                             text="Multiple samples ", \
                             onvalue=1, offvalue=0, variable=var_RootFoldChk)
    RootFoldChk.grid(row=0, column=3)

    def LoadDECTmutltiple():
        global dirnameroot
        dirnameroot = tkFileDialog.askdirectory(parent=stepZero, initialdir='/home/', title="Select Root Folder")



    loadRootbtn = Button(stepZero, text="Choose Root Directory", command=LoadDECTmutltiple)
    loadRootbtn.grid(row=1, column=3, sticky='W', padx=5, pady=2)


############
#Input Pane#
############

    label3 = Label(stepOne, text='M')
    label3.grid(row=2, column=0, sticky='E', padx=5, pady=2)
    M = Entry(stepOne, width=25)
    M.grid(row=2, column=1, sticky='E', padx=5, pady=2)

    label4 = Label(stepOne, text='P')
    label4.grid(row=3, column=0, sticky='E', padx=5, pady=2)
    P = Entry(stepOne, width=25)
    P.grid(row=3, column=1, sticky='E', padx=5, pady=2)

    label5 = Label(stepOne, text='Q')
    label5.grid(row=4, column=0, sticky='E', padx=5, pady=2)
    Q = Entry(stepOne, width=25)
    Q.grid(row=4, column=1, sticky='E', padx=5, pady=2)

    label6 = Label(stepOne, text='R')
    label6.grid(row=5, column=0, sticky='E', padx=5, pady=2)
    R = Entry(stepOne, width=25)
    R.grid(row=5, column=1, sticky='E', padx=5, pady=2)

    label7 = Label(stepOne, text='S')
    label7.grid(row=6, column=0, sticky='E', padx=5, pady=2)
    S = Entry(stepOne, width=25)
    S.grid(row=6, column=1, sticky='E', padx=5, pady=2)

    label8 = Label(stepOne, text='T')
    label8.grid(row=7, column=0, sticky='E', padx=5, pady=2)
    T = Entry(stepOne, width=25)
    T.grid(row=7, column=1, sticky='E', padx=5, pady=2)

    def saveCallback():
        with open("output.txt", "w") as f:
            f.write(RhoMat.get() + "\t")
            f.write(RhoFluid.get() + "\t")
            f.write(M.get() + "\t")
            f.write(P.get() + "\t")
            f.write(Q.get() + "\t")
            f.write(R.get() + "\t")
            f.write(S.get() + "\t")
            f.write(T.get() + "\t")

    savebtn = Button(stepOne, text="Save to File", command=saveCallback)
    savebtn.grid(row=8, column=3, sticky='W', padx=5, pady=2)


    def SetDefault ():
        M.delete(0, END)
        P.delete(0, END)
        Q.delete(0, END)
        R.delete(0, END)
        S.delete(0, END)
        T.delete(0, END)
        M.insert(0,"-0.77")
        P.insert(0,"1.98")
        Q.insert(0,"1007")
        R.insert(0,"36597.06")
        S.insert(0,"-35330.83")
        T.insert(0,"233946.02")

    defaultbtn = Button(stepOne, text="Default Parameters", command=SetDefault)
    defaultbtn.grid(row=8, column=0, sticky='W', padx=5, pady=2)

    def LoadFileInput ():
        from tkFileDialog import askopenfilename
        Tk().withdraw()
        filename = askopenfilename()
        f = open(filename, 'r')
        line = f.readline()
        M2,P2,Q2,R2,S2,T2=line.split()
        M.delete(0, END)
        M.insert(0,M2)
        P.delete(0, END)
        P.insert(0,P2)
        Q.delete(0, END)
        Q.insert(0,Q2)
        R.delete(0, END)
        R.insert(0,R2)
        S.delete(0, END)
        S.insert(0,S2)
        T.delete(0, END)
        T.insert(0,T2)
        f.close()

    loadbtn = Button(stepOne, text="Load from File", command=LoadFileInput)
    loadbtn.grid(row=8, column=1, sticky='W', padx=5, pady=2)

###############
#Cropping Pane#
###############
    var_getImaChk = IntVar()

    label10 = Label(stepTwo, text='Core Diameter (mm)')
    label10.grid(row=0, column=0, sticky='E', padx=5, pady=2)
    CoreDia = Entry(stepTwo, width=25)
    CoreDia.grid(row=0, column=1, sticky='E', padx=5, pady=2)

    label11 = Label(stepTwo, text='Crop_Frac')
    label11.grid(row=1, column=0, sticky='E', padx=5, pady=2)
    Crop_Frac = Entry(stepTwo, width=25)
    Crop_Frac.grid(row=1, column=1, sticky='E', padx=5, pady=2)

    label12 = Label(stepTwo, text='Offset X (mm)')
    label12.grid(row=2, column=0, sticky='E', padx=5, pady=2)
    offsetx = Entry(stepTwo, width=25)
    offsetx.grid(row=2, column=1, sticky='E', padx=5, pady=2)

    label13 = Label(stepTwo, text='Offset Y (mm)')
    label13.grid(row=3, column=0, sticky='E', padx=5, pady=2)
    offsety = Entry(stepTwo, width=25)
    offsety.grid(row=3, column=1, sticky='E', padx=5, pady=2)

    label14 = Label(stepTwo, text='Remove top slices(number)')
    label14.grid(row=4, column=0, sticky='E', padx=5, pady=2)
    paddingtop = Entry(stepTwo, width=25)
    paddingtop.grid(row=4, column=1, sticky='E', padx=5, pady=2)

    label15 = Label(stepTwo, text='Remove bottom slices(number)')
    label15.grid(row=5, column=0, sticky='E', padx=5, pady=2)
    paddingbottom = Entry(stepTwo, width=25)
    paddingbottom.grid(row=5, column=1, sticky='E', padx=5, pady=2)

    getImaChk = Checkbutton(stepTwo, \
                           text="Check cropping image",\
                           onvalue=1, offvalue=0,variable=var_getImaChk)
    getImaChk.grid(row=6, column=1, columnspan=3, pady=2, sticky='W')


##################
#Calculation Pane#
##################
    var_getMaskPixel = IntVar()

    label1 = Label(stepThree, text='Matrix Density (kg.m-3)')
    label1.grid(row=0, column=0, sticky='E', padx=5, pady=2)
    RhoMat = Entry(stepThree, width=25)
    RhoMat.grid(row=0, column=1, sticky='E', padx=5, pady=2)

    label2 = Label(stepThree, text='Fluid Density (kg.m-3)')
    label2.grid(row=1, column=0, sticky='E', padx=5, pady=2)
    RhoFluid = Entry(stepThree, width=25)
    RhoFluid.grid(row=1, column=1, sticky='E', padx=5, pady=2)

    def SetDefaultDens ():
        RhoMat.delete(0, END)
        RhoFluid.delete(0, END)
        RhoMat.insert(0, "2600")
        RhoFluid.insert(0, "1")


    defaultbtn = Button(stepThree, text="Default Parameters", command=SetDefaultDens)
    defaultbtn.grid(row=4, column=0, sticky='W', padx=5, pady=2)

    getMaskPixel = Checkbutton(stepThree, text="Remove Saturated Pixel from Calculations", onvalue=1, offvalue=0,variable=var_getMaskPixel)
    getMaskPixel.grid(row=3, column=0, columnspan=3, pady=2, sticky='W')

#############
#Figure Pane#
#############



def GetValues():
    global StopRun
    StopRun=False
    global Padding_top,Padding_bottom,Diameter,Crop_pct,Offsetx,Offsety,Check,Show_Plot,Mask_saturated,Browse,rho_matrix,rho_fluid,fM,fP,fQ,fR,fS,fT

    Padding_top=float(paddingtop.get())
    Padding_bottom=float(paddingbottom.get())
    Diameter=float(CoreDia.get())
    Crop_pct=float(Crop_Frac.get())
    Offsetx=float(offsetx.get())
    Offsety=float(offsety.get())
    Check=var_getImaChk.get()
    Show_Plot=1
    Mask_saturated=var_getMaskPixel.get()
    Browse=var_RootFoldChk.get()

    rho_matrix=float(RhoMat.get())
    rho_fluid=float(RhoFluid.get())
    fM=float(M.get())
    fP=float(P.get())
    fQ=float(Q.get())
    fR=float(R.get())
    fS=float(S.get())
    fT=float(T.get())
    Main()

def Stop():
    global StopRun
    StopRun=True


f = Figure(figsize=(5, 4), dpi=100)

# a tk.DrawingArea
canvas = FigureCanvasTkAgg(f, master=stepFour)
canvas.show()
canvas.get_tk_widget().grid(row=0, column=0, columnspan=3, pady=2, sticky='W')



#############
#Info Pane  #
#############

Tbox = Text(stepFive, height=10, width=70)
Tbox.pack()

pb_hD = ttk.Progressbar(stepFive, orient='horizontal', mode='determinate',length=500)
pb_hD.pack(expand=True, side=TOP)

Runbtn = Button(form, text="Run", command=lambda :GetValues())
Runbtn.grid(row=4, sticky='W', padx=5, pady=2)

Stopbtn = Button(form, text="Stop", command=lambda :Stop())
Stopbtn.grid(row=4, sticky='W', padx=50, pady=2)

form.mainloop()
