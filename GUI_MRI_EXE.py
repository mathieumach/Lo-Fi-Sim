import numpy as np
import scipy.fft as fft
import scipy.io
import os
import matplotlib.pyplot as plt 
import cv2
import math
import cmath
import datetime
import matplotlib.image
import sys
import webbrowser
import interactiveAxesPlot_SE as pltaxes_SE
import interactiveAxesPlot_GE as pltaxes_GE
import interactiveAxesPlot_IN as pltaxes_IN
import interactiveAxesPlot_DIN as pltaxes_DIN
import interactiveAxesPlot_FLAIR as pltaxes_FLAIR
import interactiveAxesPlot_SSFP as pltaxes_SSFP
import interactiveAxesPlot_Dif as pltaxes_Dif
import interactiveAxesPlot_TSE as pltaxes_TSE
import interactiveAxesPlot_Lowpass as pltaxes_low
import interactiveAxesPlot_Highpass as pltaxes_high
import interactiveAxesPlot_Gausspass as pltaxes_gauss
import interactiveAxesPlot_Nonlocalpass as pltaxes_nonlocal
import interactivePlot_SNR as pltsnr
import interactiveAxesPlot_viewdataset as pltviewdatasets
import nibabel as nib
import matplotlib.image as mpimg
from Sequences import *
from scipy import signal
from scipy import constants
from scipy.ndimage import zoom
from tkinter import *
from tkinter import filedialog
from tkinter import ttk
from PIL import ImageTk, Image
from scipy import ndimage
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from scipy.interpolate import RectBivariateSpline
from numpy import inf
from scipy.ndimage import gaussian_filter
from skimage.restoration import denoise_nl_means, estimate_sigma # estimate_sigma --> wavelet-based estimator of the (Gaussian) noise standard deviation.

root = Tk()
root.title("Low field MRI simulator")
root.state('zoomed') # To have a full screen display

#/////// Functions getting the data's path to create an executable /////////#  
def resource_path(relative_path):
    """ 
    Get absolute path to resource, for PyInstaller to create the .exe
    """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)

# ///////////////////////////////////////////////////////////////////////// #
################################ Simulator ##################################
# ///////////////////////////////////////////////////////////////////////// #

#///////////////////////////// FUNCTIONS ///////////////////////////////////#       
#//////////////////// Functions for the tutorials //////////////////////////#  
def forward(image_num):
    """
    Will change the slide in the tutorial to the next one
    Input:  image_num --> the number of the slide
    Output: none
    """
    # To be able to update the buttons, need global buttons
    global Image_label
    global forward_button
    global back_button

    Image_label.grid_forget()                                    # To remove the original image in my_label when we click forward
    Image_label = Label(Tutorial, image=image_list[image_num-1]) # Going from 1st to 2nd image means going from 0 to 1, we passed 2 so need a -1

    # Updates, so that forward doesn't have the '2' as input
    forward_button = Button(Tutorial, text=">>", command=lambda: forward(image_num+1)) # Want the next image so +1
    back_button = Button(Tutorial, text="<<", command=lambda: back(image_num-1))       # Want the previous image so -1

    if image_num == len(image_list):
        forward_button = Button(Tutorial, text=">>", state=DISABLED)

    Image_label.grid(row=0,column=0,columnspan=3)
    back_button.grid(row=1,column=0)
    forward_button.grid(row=1,column=2)
    
    # We want the status to be updated is we press forward button
    status = Label(Tutorial, text="Image " + str(image_num)  + " of " + str(len(image_list)), bd=1, relief=SUNKEN, anchor=E)
    status.grid(row=2,column=0, columnspan=3, sticky=W+E)

def back(image_num):
    """
    Will change the slide in the tutorial to the previous one
    Input:  image_num --> the number of the slide
    Output: none
    """
    # To be able to update the buttons, need global buttons
    global Image_label
    global forward_button
    global back_button

    Image_label.grid_forget()
    Image_label = Label(Tutorial, image=image_list[image_num-1])    # Going from 1st to 2nd image means going from 0 to 1, we passed 2 so need a -1

    # Updates
    forward_button = Button(Tutorial, text=">>", command=lambda: forward(image_num+1)) # Want the next image so +1
    back_button = Button(Tutorial, text="<<", command=lambda: back(image_num-1))       # Want the previous image so -1

    if image_num == 1:
        back_button = Button(Tutorial, text="<<", state=DISABLED)

    Image_label.grid(row=0,column=0,columnspan=3)
    back_button.grid(row=1,column=0)
    forward_button.grid(row=1,column=2)
    
    # We want the status to be updated is we press back button
    status = Label(Tutorial, text="Image " + str(image_num)  + " of " + str(len(image_list)), bd=1, relief=SUNKEN, anchor=E)
    status.grid(row=2,column=0, columnspan=3, sticky=W+E)

def tutorial_simu():
    """
    Will open a new window with the first slide of the tutorial and enable the user to change slides or exit the tutorial
    Input: none 
    Output: none 
    """
    global Tutorial
    global Image_label
    global image_list
    Tutorial = Toplevel() # New window
    Tutorial.title("Tutorial")
    
    # The images have to come from the correct folder
    image1 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide1.jpg")))
    image2 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide2.jpg")))
    image3 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide3.jpg")))
    image4 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide4.jpg")))
    image5 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide5.jpg")))
    image6 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide6.jpg")))
    image7 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide7.jpg")))
    image8 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide8.jpg")))
    image9 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide9.jpg")))
    image10 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide10.jpg")))
    image11 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide11.jpg")))
    image12 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide12.jpg")))
    image13 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide13.jpg")))
    image14 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide14.jpg")))
    image15 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide15.jpg")))
    image16 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide16.jpg")))
    image17 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide17.jpg")))
    image18 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide18.jpg")))
    image19 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide19.jpg")))
    image20 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide20.jpg")))
    image21 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide21.jpg")))
    image22 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide22.jpg")))
    image23 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide23.jpg")))
    image24 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide24.jpg")))
    image25 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide25.jpg")))
    image26 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial simulator\Slide26.jpg")))
    
    # Creating a list of all images
    image_list = [image1, image2, image3, image4, image5, image6, image7, image8, image9, image10, image11, image12, image13, image14, image15, image16, image17, image18, image19, image20, image21, image22, image23, image24, image25, image26] 
    
    # The status, bd for border, relief for sunnken style, anchor E -> so it is fixed at buttom right
    status = Label(Tutorial, text="Image 1 of " + str(len(image_list)), bd=1, relief=SUNKEN, anchor=E)
    # sticky for stretching, here from west to east
    status.grid(row=2,column=0, columnspan=3, sticky=W+E)
    
    Image_label = Label(Tutorial, image=image1) # Show the first image
    Image_label.grid(row=0,column=0,columnspan=3)
    
    # Buttons
    back_button = Button(Tutorial, text="<<", command=back, state = DISABLED) # We don't want to be able to go back when the first image arrives
    exit_button = Button(Tutorial, text="Exit tutorial", command=Tutorial.destroy)
    forward_button = Button(Tutorial, text=">>", command=lambda: forward(2)) # We pass '2' to go to the next image
    back_button.grid(row=1,column=0)
    exit_button.grid(row=1,column=1)
    forward_button.grid(row=1,column=2)   
    
def tutorial_gene():
    """
    Will open a new window with the first slide of the generator tutorial and enable the user to change slides or exit the tutorial
    Input: none 
    Output: none 
    """
    global Tutorial
    global Image_label
    global image_list
    Tutorial = Toplevel() # New window
    Tutorial.title("Tutorial")
    
    # The images have to come from the correct folder
    image1 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial generator\Slide1.jpg")))
    image2 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial generator\Slide2.jpg")))
    image3 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial generator\Slide3.jpg")))
    image4 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial generator\Slide4.jpg")))
    image5 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial generator\Slide5.jpg")))
    image6 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial generator\Slide6.jpg")))
    image7 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial generator\Slide7.jpg")))
    image8 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial generator\Slide8.jpg")))
    image9 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial generator\Slide9.jpg")))
    image10 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial generator\Slide10.jpg")))
    image11 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial generator\Slide11.jpg")))
    image12 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial generator\Slide12.jpg")))
    image13 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial generator\Slide13.jpg")))
    image14 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial generator\Slide14.jpg")))
    image15 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial generator\Slide15.jpg")))
    image16 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial generator\Slide16.jpg")))
    image17 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial generator\Slide17.jpg")))
    image18 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial generator\Slide18.jpg")))
    image19 = ImageTk.PhotoImage(Image.open(resource_path("Tutorial images\Tutorial generator\Slide19.jpg")))
    
    # Creating a list of all images
    image_list = [image1, image2, image3, image4, image5, image6, image7, image8, image9, image10, image11, image12, image13, image14, image15, image16, image17, image18, image19] 
    
    # The status, bd for border, relief for sunnken style, anchor E -> so it is fixed at buttom right
    status = Label(Tutorial, text="Image 1 of " + str(len(image_list)), bd=1, relief=SUNKEN, anchor=E)
    # sticky for stretching, here from west to east
    status.grid(row=2,column=0, columnspan=3, sticky=W+E)
    
    Image_label = Label(Tutorial, image=image1) # Show the first image
    Image_label.grid(row=0,column=0,columnspan=3)
    
    # Buttons
    back_button = Button(Tutorial, text="<<", command=back, state = DISABLED) # We don't want to be able to go back when the first image arrives
    exit_button = Button(Tutorial, text="Exit tutorial", command=Tutorial.destroy)
    forward_button = Button(Tutorial, text=">>", command=lambda: forward(2)) # We pass '2' to go to the next image
    back_button.grid(row=1,column=0)
    exit_button.grid(row=1,column=1)
    forward_button.grid(row=1,column=2)    
    
#//////////////////// Functions updating the GUI //////////////////////////#      
def UPDATE(event):
    """
    Updating function when keyboard buttons are pressed
    Input: event --> when a keyboard is pressed
    Output: None
    """
    if SIM == 1 and GEN == 0: # Simulator is open
        #the time of scan, bandwidth/pixel, the size of the data matrix, the minimum TE (or TR in SSFP) and the minimum TE and the b coefficient for the diffusion sequence are automatically updated
        global Pre_def_seq
        global TR
        global TE
        global TI
        global TI2
        global FOV
        global FOV_check
        global Data_mat
        global Data_matrix_check
        global Resolution
        global Res_check
        global Bandwidth
        global dif_TEmin
        global minimumTE
        global minimumTR    
        global Time_scan_num_label; global time  
        global s
        global Bd_by_pixel_label;   global Bd_by_pixel
        global SNR_num_label
        global minimumTE_num_label
        global minimumTR_num_label
        global TEeff_val_label
        global Post_TSE_TEeff_label_val
        global Bval_label
        global TEminval_label
        global B
        global dif_TEmin
        TR = TR_entry.get();        TR = int(TR);   TR = np.divide(TR,1000) # Divide by 1000 to have the values in milli
        TE = TE_entry.get();        TE = int(TE);   TE = np.divide(TE,1000)
        TI = TI_entry.get();        TI = int(TI);   TI = np.divide(TI,1000)
        TI2 = TI2_entry.get();      TI2 = int(TI2); TI2 = np.divide(TI2,1000)
        fov1 = FOV1_entry.get();    fov1 = int(fov1)
        fov2 = FOV2_entry.get();    fov2 = int(fov2)
        fov3 = FOV3_entry.get();    fov3 = int(fov3)
        res1 = Res1_entry.get();    res1 = float(res1) 
        res2 = Res2_entry.get();    res2 = float(res2)
        res3 = Res3_entry.get();    res3 = float(res3)
        dat1 = Data_mat1_entry.get(); dat1 = int(dat1)
        dat2 = Data_mat2_entry.get(); dat2 = int(dat2)
        dat3 = Data_mat3_entry.get(); dat3 = int(dat3)
        bd = Bandwidth_entry.get(); Bandwidth = int(bd)
        Alpha = Alpha_entry.get();  Alpha = int(Alpha)
        ETL = ETL_entry.get();      ETL = int(ETL)
        FOV = [fov1, fov2, fov3]
        Resolution = [res1, res2, res3]
        Data_mat = [dat1, dat2, dat3]

        if np.prod(Data_mat) != Data_matrix_check: # User changed the data matrix size

            Resolution[0] = np.round(np.divide(FOV[0],Data_mat[0]),2)
            Resolution[1] = np.round(np.divide(FOV[1],Data_mat[1]),2)
            Resolution[2] = np.round(np.divide(FOV[2],Data_mat[2]),2)
            Res1_entry.delete(0,END)
            Res1_entry.insert(0,float(Resolution[0]))
            Res2_entry.delete(0,END)
            Res2_entry.insert(0,float(Resolution[1]))
            Res3_entry.delete(0,END)
            Res3_entry.insert(0,float(Resolution[2]))

        elif (np.prod(FOV) != FOV_check) or (np.prod(Resolution) != Res_check): # User changed the FOV or the resolution

            Data_mat = [int(np.round(np.divide(FOV[0], Resolution[0]))), int(np.round(np.divide(FOV[1], Resolution[1]))), int(np.round(np.divide(FOV[2], Resolution[2])))]
            Data_mat1_entry.delete(0,END)
            Data_mat1_entry.insert(0,int(Data_mat[0]))
            Data_mat2_entry.delete(0,END)
            Data_mat2_entry.insert(0,int(Data_mat[1]))
            Data_mat3_entry.delete(0,END)
            Data_mat3_entry.insert(0,int(Data_mat[2]))

        # Updating the values to check for next update
        Data_matrix_check = np.prod(Data_mat)
        FOV_check = np.prod(FOV)
        Res_check = np.prod(Resolution)

        time_scan = TR * Data_mat[1] * Data_mat[2]

        if Pre_def_seq == 'TSE':
            time_scan = time_scan/ETL
            trajectory = traj.get()

            if trajectory == 'Linear':
                eff = 0.5 * TE * ETL
            elif trajectory == 'In-out':
                eff = TE
            elif trajectory == 'Out-in':
                eff = TE * ETL   

            TEeff_val_label.grid_forget()
            TEeff_val_label = Label(frame3, text = np.round(1000*eff,2), font = ("Helvetica", 12)); TEeff_val_label.grid(row = 8, column = 3) 

        if Pre_def_seq == 'IN' or Pre_def_seq == 'Double IN' or Pre_def_seq == 'FLAIR' or Pre_def_seq == 'Dif':
            post_TSE_TE = np.divide(int(post_TSE_TE_entry.get()),1000)
            post_TSE_ETL = int(post_TSE_ETL_entry.get())
            posttraj = post_traj.get()

            if posttraj == 'Linear':
                posteff = 0.5 * post_TSE_TE * post_TSE_ETL
            elif posttraj == 'In-out':
                posteff = post_TSE_TE
            elif posttraj == 'Out-in':
                posteff = post_TSE_TE * post_TSE_ETL   
            Post_TSE_TEeff_label_val.grid_forget()
            Post_TSE_TEeff_label_val = Label(frame3, text = np.round(1000*posteff,2), font = ("Helvetica", 12))
            Post_TSE_TEeff_label_val.grid(row = 15, column = 3)

        if Pre_def_seq == 'IN':
            time_scan = (TR+TI) * Data_mat[1] * Data_mat[2]
            time_scan = time_scan/ETL

        if Pre_def_seq == 'Double IN':
            time_scan = (TR+TI+TI2) * Data_mat[1] * Data_mat[2]
            time_scan = time_scan/ETL

        if Pre_def_seq == 'FLAIR':
            time_scan = (TR+TI) * Data_mat[1] * Data_mat[2]
            time_scan = time_scan/ETL

        time = datetime.timedelta(seconds=time_scan) # Converts amount of seconds to hours:minutes:seconds
        Bd_by_pixel = np.divide(Bandwidth,Data_mat[0])

        Time_scan_num_label.grid_forget()
        Bd_by_pixel_label.grid_forget()

        Time_scan_num_label = Label(frame11, text = str(time), font=("Helvetica", 12)); Time_scan_num_label.grid(row = 0, column = 1)   
        Bd_by_pixel_label = Label(frame11, text = str(round(Bd_by_pixel,2)), font=("Helvetica", 12)); Bd_by_pixel_label.grid(row = 1, column = 1)

        # minimum TE and TR
        minimumTE = Data_mat[0]/(2*Bandwidth)
        minimumTR = 2*minimumTE
        minimumTE_num_label.grid_forget()

        if Pre_def_seq == 'SSFP':
            minimumTE_label = Label(frame11, text = "TRmin (ms) ", font=("Helvetica", 12)).grid(row = 3, column = 0)
            minimumTE_num_label = Label(frame11, text = np.round(minimumTR * 1000,2), font=("Helvetica", 12))
            minimumTE_num_label.grid(row = 3, column = 1)
        else:
            minimumTE_label = Label(frame11, text = "TEmin (ms)", font=("Helvetica", 12)).grid(row = 3, column = 0)
            minimumTE_num_label = Label(frame11, text = np.round(minimumTE * 1000,2), font=("Helvetica", 12))
            minimumTE_num_label.grid(row = 3, column = 1)

        # B and minTE for diffusion sequence
        if Pre_def_seq == 'Dif':
            gamma =  42.58*(10**6)#*2*constants.pi      # gyromagnetic ratio for hydrogen 42.58 [MHz/T] 

            G = G_entry.get(); G = float(G); G = np.divide(G,1000)
            smalldelta = smalldelta_entry.get(); smalldelta = float(smalldelta); smalldelta = np.divide(smalldelta,1000)
            bigdelta = bigdelta_entry.get(); bigdelta = float(bigdelta); bigdelta = np.divide(bigdelta,1000)

            B = (bigdelta - smalldelta/3)*(gamma*G*smalldelta)**2; 
            dif_TEmin = smalldelta + bigdelta; 

            Bval_label.grid_forget()
            Bval_label = Label(frame3, text = np.round(B,2), font=("Helvetica", 12)); Bval_label.grid(row = 7, column = 4)
            TEminval_label.grid_forget()
            TEminval_label = Label(frame3, text = np.round(dif_TEmin*1000,2), font=("Helvetica", 12)); TEminval_label.grid(row = 8, column = 4)
            
    elif SIM == 0 and GEN == 1: # Generator is open
        
        #Updating the b coefficient for the diffusion sequence and effective TE are automatically 
        global TEeff_val_label_gene
        global Bval_label_gene
        global TEminval_label_gene
        global FOV_check_gene
        global Res_check_gene
        global Data_matrix_check_gene
        FOV_gene = [int(FOV1_entry_gene.get()), int(FOV2_entry_gene.get()), int(FOV3_entry_gene.get())]
        Resolution_gene = [float(Res1_entry_gene.get()), float(Res2_entry_gene.get()), float(Res3_entry_gene.get())]
        Data_mat_gene = [int(Data_mat1_entry_gene.get()), int(Data_mat2_entry_gene.get()), int(Data_mat3_entry_gene.get())]

        if np.prod(Data_mat_gene) != Data_matrix_check_gene: # User changed the data matrix size
            Resolution_gene[0] = np.round(np.divide(FOV_gene[0],Data_mat_gene[0]),2)
            Resolution_gene[1] = np.round(np.divide(FOV_gene[1],Data_mat_gene[1]),2)
            Resolution_gene[2] = np.round(np.divide(FOV_gene[2],Data_mat_gene[2]),2)
            Res1_entry_gene.delete(0,END)
            Res1_entry_gene.insert(0,float(Resolution_gene[0]))
            Res2_entry_gene.delete(0,END)
            Res2_entry_gene.insert(0,float(Resolution_gene[1]))
            Res3_entry_gene.delete(0,END)
            Res3_entry_gene.insert(0,float(Resolution_gene[2]))

        elif (np.prod(FOV_gene) != FOV_check_gene) or (np.prod(Resolution_gene) != Res_check_gene): # User changed the FOV or the resolution
            Data_mat_gene = [int(np.round(np.divide(FOV_gene[0], Resolution_gene[0]))), int(np.round(np.divide(FOV_gene[1], Resolution_gene[1]))), int(np.round(np.divide(FOV_gene[2], Resolution_gene[2])))]
            Data_mat1_entry_gene.delete(0,END)
            Data_mat1_entry_gene.insert(0,int(Data_mat_gene[0]))
            Data_mat2_entry_gene.delete(0,END)
            Data_mat2_entry_gene.insert(0,int(Data_mat_gene[1]))
            Data_mat3_entry_gene.delete(0,END)
            Data_mat3_entry_gene.insert(0,int(Data_mat_gene[2]))

        # Updating the values to check for next update
        Data_matrix_check_gene = np.prod(Data_mat_gene)
        FOV_check_gene = np.prod(FOV_gene)
        Res_check_gene = np.prod(Resolution_gene)

        if Pre_def_seq_gene == 'TSE':
            eff = 0
            TE = TE_entry_gene.get(); TE = int(TE); TE = np.divide(TE,1000)
            ETL = int(ETL_entry_gene.get()) 
            trajectory_gene = traj_gene.get()
            if trajectory_gene == 'Linear':
                eff = 0.5 * TE * ETL
            elif trajectory_gene == 'In-out':
                eff = TE
            elif trajectory_gene == 'Out-in':
                eff = TE * ETL
            TEeff_val_label_gene.grid_forget() 
            TEeff_val_label_gene = Label(frame7, text = str(int(1000*eff)), font = ("Helvetica", f))
            TEeff_val_label_gene.grid(row = 8, column = 3) 

        # B for diffusion sequence
        if Pre_def_seq_gene == 'Dif':
            gamma =  42.58*(10**6)#*2*constants.pi      # gyromagnetic ratio for hydrogen 42.58 [MHz/T] 

            G = G_entry_gene.get(); G = float(G); G = np.divide(G,1000)
            smalldelta = smalldelta_entry_gene.get(); smalldelta = float(smalldelta); smalldelta = np.divide(smalldelta,1000)
            bigdelta = bigdelta_entry_gene.get(); bigdelta = float(bigdelta); bigdelta = np.divide(bigdelta,1000)

            B = (bigdelta - smalldelta/3)*(gamma*G*smalldelta)**2; 
            dif_TEmin = smalldelta + bigdelta;

            Bval_label_gene.grid_forget()
            Bval_label_gene = Label(frame7, text = str(round(B,2)), font=("Helvetica", f)); Bval_label_gene.grid(row = 8, column = 4)
            TEminval_label_gene.grid_forget()
            TEminval_label_gene = Label(frame7, text = str(round(1000*dif_TEmin,2)), font=("Helvetica", f));  TEminval_label_gene.grid(row = 9, column = 4)

        #Update of effective TE for post TSE sequence
        if Pre_def_seq_gene == 'IN' or Pre_def_seq_gene == 'Double IN' or Pre_def_seq_gene == 'FLAIR' or Pre_def_seq_gene == 'Dif':
            post_showTEeff_gene()

#### Function resting the parameters #####        
def reset():
    """
    Reset all the parameters to their inital value and plot the initial proton density images
    Input: None 
    Output: None
    """
    global Time_scan_num_label
    global Bd_by_pixel_label
    global SNR_num_label
    global minimumTE_num_label
    global minimumTR_num_label
    global warning_label
    global Pre_def_seq
    global seq_label2
    global SNR_length
    global SNR_noise_box_center
    global SNR_mean_box_center
    global canvas1
    global canvas2
    global canvas3
    global FOV_check
    global Res_check
    global Data_matrix_check
    
    FOV_check = 20625000 # this value is 250*300*275 and will be used to see if the FOV was changed by the user
    Res_check = 9.81045 # this value is 1.95*2.34*2.15 and will be used to see if the resolution was changed by the user
    Data_matrix_check = 2097152 # this value is 128*128*128 and will be used to see if the Data matrix size was changed by the user
    
    canvas1.get_tk_widget().destroy()
    canvas2.get_tk_widget().destroy()
    canvas3.get_tk_widget().destroy()
    
    SNR_length = 5
    SNR_noise_box_center = [105,20]
    SNR_mean_box_center = [50,80]
    
    J = [1, 2, 3]
    num = [128, 150, 135]
    
    # Create axial plot
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.rot90(M0_3D[:,:,num[2]]), cmap='gray')  # Create image plot
    plt.axis('off')
    canvas1 = FigureCanvasTkAgg(fig, root)                  # Tkinter canvas which contains matplotlib figure
    canvas1.draw()
    canvas1.get_tk_widget().grid(row = 1, column = 3, rowspan = 4)       # Placing canvas on Tkinter window
    plt.close()

    # Create coronal plot
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.rot90(M0_3D[:,num[1],:]), cmap='gray') 
    plt.axis('off')
    canvas2 = FigureCanvasTkAgg(fig, root)         
    canvas2.draw()
    canvas2.get_tk_widget().grid(row = 1, column = 2, rowspan = 4)
    plt.close()

    # Create sagittal plot
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.flip(np.rot90(M0_3D[num[0],:,:]), axis = 1), cmap='gray')
    plt.axis('off')
    canvas3 = FigureCanvasTkAgg(fig, root)         
    canvas3.draw()
    canvas3.get_tk_widget().grid(row = 1, column = 1, rowspan = 4) 
    plt.close()

    TR_entry.delete(0, END); TR_entry.insert(0, '500')
    TE_entry.delete(0, END); TE_entry.insert(0, '20')
    TI_entry.delete(0, END); TI_entry.insert(0, '250')
    FOV1_entry.delete(0, END); FOV1_entry.insert(0, '250')
    FOV2_entry.delete(0, END); FOV2_entry.insert(0, '300')
    FOV3_entry.delete(0, END); FOV3_entry.insert(0, '275')
    Res1_entry.delete(0, END); Res1_entry.insert(0, '1.95')
    Res2_entry.delete(0, END); Res2_entry.insert(0, '2.34')
    Res3_entry.delete(0, END); Res3_entry.insert(0, '2.15')
    Data_mat1_entry.delete(0,END); Data_mat1_entry.insert(0,'128')
    Data_mat2_entry.delete(0,END); Data_mat2_entry.insert(0,'128')
    Data_mat3_entry.delete(0,END); Data_mat3_entry.insert(0,'128')
    Bandwidth_entry.delete(0, END); Bandwidth_entry.insert(0, '40000')
    Alpha_entry.delete(0, END); Alpha_entry.insert(0, '45')
    G_entry.delete(0,END); G_entry.insert(0, '10')
    smalldelta_entry.delete(0,END); smalldelta_entry.insert(0, '1')
    bigdelta_entry.delete(0,END);bigdelta_entry.insert(0, '2')

    Time_scan_num_label.grid_forget()
    Bd_by_pixel_label.grid_forget()
    SNR_num_label.grid_forget()
    minimumTE_num_label.grid_forget()
    minimumTR_num_label.grid_forget()
    warning_label.grid_forget()
    
    minimumTE = 0.00128 * 1000 # *1000 to have it in ms
    minimumTR = 2 * minimumTE
    Time_scan_num_label = Label(frame11, text = "2:16:32", font=("Helvetica", f));         Time_scan_num_label.grid(row = 0, column = 1)
    Bd_by_pixel_label = Label(frame11, text = "390.62", font=("Helvetica", f));            Bd_by_pixel_label.grid(row = 1, column = 1)
    SNR_num_label = Label(frame11, text = "      ", font=("Helvetica", f));                SNR_num_label.grid(row = 2, column = 1)
    minimumTE_num_label = Label(frame11, text = minimumTE, font=("Helvetica", f));         minimumTE_num_label.grid(row = 3, column = 1)
    minimumTR_num_label = Label(frame11, text = minimumTR, font=("Helvetica", f));         minimumTR_num_label.grid(row = 3, column = 1)
    minimumTR_num_label.grid_forget() 
    
    remove_widgets(Pre_def_seq)
    Pre_def_seq = " "
    
    seq_label2.grid_forget()
    seq_label2 = Label(frame1, text = " ", font=("Helvetica", 18))
    seq_label2.grid(row = 10, column = 0, rowspan = 2)
        
###### Functions regarding the noise and SNR ######
def noise_generation(Data_mat, Resolution, Bandwidth, B0_field):
    """
    Creates a random noise tensor of same size as 'Data_mat' (data shape of simulated data) 
    The std of the noise to add is based on the resolution, bandwidth and number of voxels in the final data
    
    Input: Data_mat   --> shape of the simulated data
           Resolution --> resolution of the simulated data
           Bandwidth  --> bandwidth
           B0_field   --> B0 field strength
                      
    Output: n --> noise tensor
    """
    tot_res = Resolution[0]*Resolution[1]*Resolution[2]
    tot_data = Data_mat[0]*Data_mat[1]*Data_mat[2]

    # Equation below links the SNR with the resolution, bandwidth and FOV
    snr_relation = tot_res / np.sqrt(np.divide(Bandwidth,tot_data))
    
    # snr_relation = 71,69 (for Yiming's scan parameters) and snr ~= 12,6 and std ~= 10,9 from Yiming's scan
    # --> 10,9 = const/71,69 => const = 781,4
    #print(snr_relation) 71.69
    
    if B0_field == 'B0_46mT':
        const = 781
    elif B0_field == 'B0_15T':
        const = 100
    elif B0_field == 'B0_3T':
        const = 1

    std = const / snr_relation
    n = np.abs(np.random.normal(0,std,Data_mat))  
    return n

#// Function computing SNR based on 2 boxes defined by the parameters (box "s" is in a noisy region of the image (to compute the std of noise) and box "m" is in a image region where their is signal (to compute the mean of the signal)) //#
def snr_homemade(im,s1,s2,s3,s4,m1,m2,m3,m4):
    """
    Computes the SNR of an image as the ration of the mean in the box defined by points {m1,m2,m3,m4} by the standard deviation in the bx defined by the points {s1,s2,s3,s4}
    Input: im --> image
           s1,s2 --> top and bottum of "noise box" (height)
           s3,s4 --> left and rigth of "noise box" (width)
           m1,m2 --> top and bottum of "signal box" (height)
           m3,m4 --> left and rigth of "signal box" (width)
           
    Output: snr --> signal to noise ration of im
    """
    m = np.mean(im[m1:m2,m3:m4])
    std = np.std(im[s1:s2,s3:s4])
    
    if std == 0:
        snr = 1000
    else:
        snr = m/std    
    return snr

def SNR_measure():
    """
    Creates an interactive plot to visualise and define the boxes used in the computation of the SNR
    Input: None       
    Output: None 
    """
    global SNR_length
    global SNR_noise_box_center
    global SNR_mean_box_center
    global B0_field
    global SNR_num_label
    global final_im_ax
    
    TR = TR_entry.get(); TR = int(TR); TR = np.divide(TR,1000) # Divide by 1000 to have the values in milli
    TE = TE_entry.get(); TE = int(TE); TE = np.divide(TE,1000)
    fov1 = int(FOV1_entry.get()); fov2 = int(FOV2_entry.get()); fov3 = int(FOV3_entry.get()); FOV = [fov1, fov2, fov3]
    res1 = float(Res1_entry.get()); res2 = float(Res2_entry.get()); res3 = float(Res3_entry.get()); Resolution = [res1, res2, res3]
    Bandwidth = int(Bandwidth_entry.get())
    Data_mat = np.divide(FOV, Resolution)
    Data_mat = [int(Data_mat[0]), int(Data_mat[1]), int(Data_mat[2])]  
    
    fig, ax = plt.subplots() # figsize=(4, 4)
    figSNR = pltsnr.interactivePlot(fig, ax, final_im_ax, Data_mat, snr_homemade, plotAxis = 2, fov = FOV)  
    plt.show() 
    
    # REMARK the coordinates in the boxes center vectors are defined on the data rotated 90 degree!
    SNR_length = figSNR.length
    SNR_noise_box_center = [int(figSNR.noise_box_center[0]), int(figSNR.noise_box_center[1])]
    SNR_mean_box_center = [int(figSNR.mean_box_center[0]), int(figSNR.mean_box_center[1])]
    
    if np.sum(SNR_noise_box_center) != 0 and np.sum(SNR_mean_box_center) != 0:
        SNR_num_label.grid_forget()
        SNR_num_label = Label(frame11, text = str(np.abs(round(figSNR.snr,2))), font=("Helvetica", 12))
        SNR_num_label.grid(row = 2, column = 1) 
        
def SNR_visu():
    """
    Creates an interactive plot to visualise how the SNR is currently defined
    Input: None       
    Output: None 
    """
    global SNR_length
    global SNR_noise_box_center
    global SNR_mean_box_center
    global B0_field
    global SNR_num_label
    global final_im_ax
    
    fov1 = int(FOV1_entry.get()); fov2 = int(FOV2_entry.get()); fov3 = int(FOV3_entry.get()); FOV = [fov1, fov2, fov3]
    res1 = float(Res1_entry.get()); res2 = float(Res2_entry.get()); res3 = float(Res3_entry.get()); Resolution = [res1, res2, res3]
    Data_mat = np.divide(FOV, Resolution)
    Data_mat = [int(Data_mat[0]), int(Data_mat[1]), int(Data_mat[2])]  
    
    arr2 = np.zeros((Data_mat[1], Data_mat[0])) # This array is only zeros with 2 white boxes
    m = SNR_length
    M = 2*m
    SNR_length = m
    arrstd = np.ones(M)
    arrm = np.ones(M)
    arr2 = np.zeros((Data_mat[1], Data_mat[0])) # This array is only zeros with 2 white boxes
    # std box
    A = int(SNR_noise_box_center[1])
    B = int(SNR_noise_box_center[0])
    arr2[A-m:A+m, B-m] = arrstd
    arr2[A-m:A+m, B+m-1] = arrstd
    arr2[A-m, B-m:B+m] = arrstd
    arr2[A+m-1, B-m:B+m] = arrstd
    # mean box
    a = int(SNR_mean_box_center[1])
    b = int(SNR_mean_box_center[0])
    arr2[a-m:a+m, b-m] = arrm
    arr2[a-m:a+m, b+m-1] = arrm
    arr2[a-m, b-m:b+m] = arrm
    arr2[a+m-1, b-m:b+m] = arrm
        
    plt.imshow(np.rot90(final_im_ax) + 150*arr2, cmap = 'gray')
    plt.title('SNR = ' + str(snr_homemade(np.rot90(final_im_ax),A-m,A+m,B-m,B+m,a-m,a+m,b-m,b+m)))
    plt.show()
        
#/////////// Functions computing a trapeze style 3D integral //////////////#        
def Trapeze_integral_3D(data):
    """
    Performs a 3D integral of each voxel with the bordres of the voxel as integral borders
    Input: data --> 3D array 
    Output: res --> the integrated values of each voxel (3D array)
    """
    s = data.shape
    res = np.zeros((s))
    for i in range(s[0]-1):
        for j in range(s[1]-1):
            for k in range(s[2]-1): 
                a = data[i,j,k]
                b = data[i,j,k+1]
                c = data[i,j+1,k]
                d = data[i,j+1,k+1]
                e = data[i+1,j,k]
                f = data[i+1,j,k+1]
                g = data[i+1,j+1,k]
                h = data[i+1,j+1,k+1]
                minim = np.min([a,b,c,d,e,f,g,h])
                res[i,j,k] = 0.125 * ((a-minim)+(b-minim)+(c-minim)+(d-minim)+(e-minim)+(f-minim)+(g-minim)+(h-minim))
    return res
        
#/////////// Functions updating the widgets for the specific sequences //////////////#  
def remove_widgets(seq):
    """
    Removes the widgets (labels and entries) that were needed for previously selected sequence
    Input : seq --> string of Pre_def_seq
    Output: None
    """
    if seq == 'SE':
        TE_entry.grid_forget();         TE_label.grid_forget()
    elif seq == 'GE':
        TE_entry.grid_forget();         TE_label.grid_forget()
        Alpha_entry.grid_forget();      Alpha_label.grid_forget()
    elif seq == 'IN':
        frame5.grid_forget()
        TE_entry.grid_forget();         TE_label.grid_forget()
        TI_entry.grid_forget();         TI_label.grid_forget()
        Post_TSE_TEeff_label.grid_forget()
        Post_TSE_TEeff_label_val.grid_forget()
    elif seq == 'Double IN':
        frame5.grid_forget()
        TE_entry.grid_forget();         TE_label.grid_forget()
        TI_entry.grid_forget();         TI_label.grid_forget()
        TI2_entry.grid_forget();        TI2_label.grid_forget()
        Post_TSE_TEeff_label.grid_forget()
        Post_TSE_TEeff_label_val.grid_forget()
    elif seq == 'FLAIR':
        frame5.grid_forget()
        TE_entry.grid_forget();         TE_label.grid_forget()
        TI_entry.grid_forget();         TI_label.grid_forget()
        Post_TSE_TEeff_label.grid_forget()
        Post_TSE_TEeff_label_val.grid_forget()
    elif seq == 'Dif':
        frame5.grid_forget()
        TE_entry.grid_forget();         TE_label.grid_forget()
        G_label.grid_forget();          G_entry.grid_forget();        
        smalldelta_label.grid_forget(); smalldelta_entry.grid_forget()
        bigdelta_label.grid_forget();   bigdelta_entry.grid_forget()
        B_label.grid_forget()
        Bval_label.grid_forget()
        TEmin_label.grid_forget()   
        TEminval_label.grid_forget()     
        Post_TSE_TEeff_label.grid_forget()
        Post_TSE_TEeff_label_val.grid_forget()
    elif seq == 'TSE':
        TE_entry.grid_forget();         TE_label.grid_forget()
        ETL_entry.grid_forget();        ETL_label.grid_forget()
        Kspace_traj_label.grid_forget();tse_drop.grid_forget()
        TEeff_label.grid_forget();      tse_read_drop.grid_forget()
        TEeff_val_label.grid_forget()
        tse_read_label.grid_forget()
    elif seq == 'SSFP':
        Alpha_entry.grid_forget();      Alpha_label.grid_forget()
        
def show_widgets(seq):
    """
    Show the widgets (labels and entries) and the minimum TE (or TR) that are needed for the selected sequence 
    Input : seq --> string of Pre_def_seq
    Output: None
    """    
    global TR
    global TE
    global FOV
    global Data_mat
    global Resolution
    global Bandwidth
    global minimumTE_num_label
    
    TR = TR_entry.get(); TR = int(TR); TR = np.divide(TR,1000) # Divide by 1000 to have the values in milli
    TE = TE_entry.get(); TE = int(TE); TE = np.divide(TE,1000)
    fov1 = FOV1_entry.get(); fov1 = int(fov1)
    fov2 = FOV2_entry.get(); fov2 = int(fov2)
    fov3 = FOV3_entry.get(); fov3 = int(fov3)
    res1 = Res1_entry.get(); res1 = float(res1) 
    res2 = Res2_entry.get(); res2 = float(res2)
    res3 = Res3_entry.get(); res3 = float(res3)
    bd = Bandwidth_entry.get(); Bandwidth = int(bd)
    FOV = [fov1, fov2, fov3]
    Resolution = [res1, res2, res3]
    Data_mat = np.divide(FOV, Resolution)
    Data_mat = [int(Data_mat[0]), int(Data_mat[1]), int(Data_mat[2])]
    
    minimumTE = Data_mat[0]/(2*Bandwidth)
    minimumTR = 2*minimumTE
    
    minimumTE_num_label.grid_forget()
    
    if seq == 'SE':
        TE_label.grid(row = 6, column = 0);         TE_entry.grid(row = 6, column = 1)
    elif seq == 'GE':
        TE_label.grid(row = 6, column = 0);         TE_entry.grid(row = 6, column = 1)
        Alpha_label.grid(row = 7, column = 0);      Alpha_entry.grid(row = 7, column = 1)
    elif seq == 'IN':
        frame5.grid(row = 4, column = 2, rowspan = 3 , columnspan = 4)
        TI_label.grid(row = 6, column = 0);         TI_entry.grid(row = 6, column = 1); TI_entry.delete(0,END); TI_entry.insert(0, '250')  
        Post_TSE_TEeff_label.grid(row = 15, column = 2)
        Post_TSE_TEeff_label_val.grid(row = 15, column = 3) 
    elif seq == 'Double IN':
        frame5.grid(row = 4, column = 2, rowspan = 3 , columnspan = 4)
        TI_label.grid(row = 6, column = 0);         TI_entry.grid(row = 6, column = 1); TI_entry.delete(0,END); TI_entry.insert(0, '250')
        TI2_label.grid(row = 7, column = 0);        TI2_entry.grid(row = 7, column = 1)
        Post_TSE_TEeff_label.grid(row = 15, column = 2)
        Post_TSE_TEeff_label_val.grid(row = 15, column = 3)
    elif seq == 'FLAIR':
        frame5.grid(row = 4, column = 2, rowspan = 3 , columnspan = 4)
        TI_label.grid(row = 6, column = 0);         TI_entry.grid(row = 6, column = 1); TI_entry.delete(0,END); TI_entry.insert(0, '1700') 
        Post_TSE_TEeff_label.grid(row = 15, column = 2)
        Post_TSE_TEeff_label_val.grid(row = 15, column = 3)
    elif seq == 'Dif':
        frame5.grid(row = 4, column = 2, rowspan = 3 , columnspan = 4)
        G_label.grid(row = 6, column = 0);          G_entry.grid(row = 6, column = 1);  
        smalldelta_label.grid(row = 7, column = 0); smalldelta_entry.grid(row = 7, column = 1)
        bigdelta_label.grid(row = 8, column = 0);   bigdelta_entry.grid(row = 8, column = 1)
        B_label.grid(row = 7, column = 3)
        Bval_label.grid(row = 7, column = 4)
        TEmin_label.grid(row = 8, column = 3)     
        TEminval_label.grid(row = 8, column = 4)
        Post_TSE_TEeff_label.grid(row = 15, column = 2)
        Post_TSE_TEeff_label_val.grid(row = 15, column = 3)            
    elif seq == 'TSE':
        TE_label.grid(row = 6, column = 0);         TE_entry.grid(row = 6, column = 1)
        ETL_label.grid(row = 7, column = 0);        ETL_entry.grid(row = 7, column = 1)
        Kspace_traj_label.grid(row = 8, column = 0);tse_drop.grid(row = 8, column = 1)
        TEeff_label.grid(row = 8, column = 2);      tse_read_drop.grid(row = 5, column = 5)
        TEeff_val_label.grid(row = 8, column = 3)
        tse_read_label.grid(row = 5, column = 4)
    elif seq == 'SSFP':
        Alpha_label.grid(row = 6, column = 0);      Alpha_entry.grid(row = 6, column = 1)      
        
    if seq == 'SSFP':
        minimumTE_label = Label(frame11, text = "TRmin (ms) ", font=("Helvetica", 12)).grid(row = 3, column = 0)
        minimumTE_num_label = Label(frame11, text = np.round(minimumTR * 1000,2), font=("Helvetica", 12))
        minimumTE_num_label.grid(row = 3, column = 1)
    else:
        minimumTE_label = Label(frame11, text = "TEmin (ms)", font=("Helvetica", 12)).grid(row = 3, column = 0)
        minimumTE_num_label = Label(frame11, text = np.round(minimumTE * 1000,2), font=("Helvetica", 12))
        minimumTE_num_label.grid(row = 3, column = 1)        
        
#///////////////////// FUNCTIONS FOR THE HISTORY FRAME /////////////////////#        
def bind():
    """
    Function linking the sequences in the history to action of right clicking on them
    Input: None
    Output: None
    """  
    global history_label1, history_label2, history_label3, history_label4, history_label5
    global h1, h2, h3, h4, h5
    
    history_label1.bind('<Button-3>', h1)
    history_label2.bind('<Button-3>', h2)
    history_label3.bind('<Button-3>', h3)
    history_label4.bind('<Button-3>', h4)
    history_label5.bind('<Button-3>', h5)

def h1(*args):
    """
    Function updating the history labels & dictionnary when right clicking on the first sequence in history calls the open history function
    Input: None
    Output: None
    """  
    global history_label1, history_label2, history_label3, history_label4, history_label5
    
    key = 'h1'
    h_open(key)
    
    l1 = history_label1.cget('text')
    d1 = history_dict['h1']
    
    if history_label2.cget('text') ==  " ... ":     
        history_label2 = Label(frame10, text = l1, font=("Helvetica", f)); history_label2.grid(row = 1, column = 0)
        history_dict['h2'] = d1
    elif history_label3.cget('text') ==  " ... ":     
        history_label3 = Label(frame10, text = l1, font=("Helvetica", f)); history_label3.grid(row = 2, column = 0)
        history_dict['h3'] = d1
    elif history_label4.cget('text') ==  " ... ":     
        history_label4 = Label(frame10, text = l1, font=("Helvetica", f)); history_label4.grid(row = 3, column = 0)
        history_dict['h4'] = d1
    elif history_label5.cget('text') ==  " ... ":     
        history_label5 = Label(frame10, text = l1, font=("Helvetica", f)); history_label5.grid(row = 4, column = 0)
        history_dict['h5'] = d1
    else:    
        history_label1.grid_forget()
        history_label1 = Label(frame10, text = history_label2.cget('text'), font=("Helvetica", f)); history_label1.grid(row=0, column=0)
        history_dict['h1'] = history_dict['h2']
        history_label2.grid_forget()
        history_label2 = Label(frame10, text = history_label3.cget('text'), font=("Helvetica", f)); history_label2.grid(row=1, column=0)
        history_dict['h2'] = history_dict['h3']
        history_label3.grid_forget()
        history_label3 = Label(frame10, text = history_label4.cget('text'), font=("Helvetica", f)); history_label3.grid(row=2, column=0)
        history_dict['h3'] = history_dict['h4']
        history_label4.grid_forget()
        history_label4 = Label(frame10, text = history_label5.cget('text'), font=("Helvetica", f)); history_label4.grid(row=3, column=0)
        history_dict['h4'] = history_dict['h5']
        history_label5.grid_forget()
        history_label5 = Label(frame10, text = l1, font=("Helvetica", f)); history_label5.grid(row=4, column=0)
        history_dict['h5'] = d1
        
    bind()
        
def h2(*args):
    """
    Function updating the history labels & dictionnary when right clicking on the second sequence in history calls the open history function
    Input: None
    Output: None
    """  
    global history_label1, history_label2, history_label3, history_label4, history_label5
    
    key = 'h2'
    h_open(key)
    
    l2 = history_label2.cget('text')
    d2 = history_dict['h2']
    
    if history_label3.cget('text') ==  " ... ":     
        history_label3 = Label(frame10, text = l2, font=("Helvetica", f)); history_label3.grid(row = 2, column = 0)
        history_dict['h3'] = d2
    elif history_label4.cget('text') ==  " ... ":     
        history_label4 = Label(frame10, text = l2, font=("Helvetica", f)); history_label4.grid(row = 3, column = 0)
        history_dict['h4'] = d2
    elif history_label5.cget('text') ==  " ... ":     
        history_label5 = Label(frame10, text = l2, font=("Helvetica", f)); history_label5.grid(row = 4, column = 0)
        history_dict['h5'] = d2
    else:    
        history_label1.grid_forget()
        history_label1 = Label(frame10, text = history_label2.cget('text'), font=("Helvetica", f)); history_label1.grid(row=0, column=0)
        history_dict['h1'] = history_dict['h2']
        history_label2.grid_forget()
        history_label2 = Label(frame10, text = history_label3.cget('text'), font=("Helvetica", f)); history_label2.grid(row=1, column=0)
        history_dict['h2'] = history_dict['h3']
        history_label3.grid_forget()
        history_label3 = Label(frame10, text = history_label4.cget('text'), font=("Helvetica", f)); history_label3.grid(row=2, column=0)
        history_dict['h3'] = history_dict['h4']
        history_label4.grid_forget()
        history_label4 = Label(frame10, text = history_label5.cget('text'), font=("Helvetica", f)); history_label4.grid(row=3, column=0)
        history_dict['h4'] = history_dict['h5']
        history_label5.grid_forget()
        history_label5 = Label(frame10, text = l2, font=("Helvetica", f)); history_label5.grid(row=4, column=0)
        history_dict['h5'] = d2
    bind()
    
def h3(*args):
    """
    Function updating the history labels & dictionnary when right clicking on the third sequence in history calls the open history function
    Input: None
    Output: None
    """  
    global history_label1, history_label2, history_label3, history_label4, history_label5
    
    key = 'h3'
    h_open(key)
    
    l3 = history_label3.cget('text')
    d3 = history_dict['h3']
    
    if history_label4.cget('text') ==  " ... ":     
        history_label4 = Label(frame10, text = l3, font=("Helvetica", f)); history_label4.grid(row = 3, column = 0)
        history_dict['h4'] = d3
    elif history_label5.cget('text') ==  " ... ":     
        history_label5 = Label(frame10, text = l3, font=("Helvetica", f)); history_label5.grid(row = 4, column = 0)
        history_dict['h5'] = d3
    else:    
        history_label1.grid_forget()
        history_label1 = Label(frame10, text = history_label2.cget('text'), font=("Helvetica", f)); history_label1.grid(row=0, column=0)
        history_dict['h1'] = history_dict['h2']
        history_label2.grid_forget()
        history_label2 = Label(frame10, text = history_label3.cget('text'), font=("Helvetica", f)); history_label2.grid(row=1, column=0)
        history_dict['h2'] = history_dict['h3']
        history_label3.grid_forget()
        history_label3 = Label(frame10, text = history_label4.cget('text'), font=("Helvetica", f)); history_label3.grid(row=2, column=0)
        history_dict['h3'] = history_dict['h4']
        history_label4.grid_forget()
        history_label4 = Label(frame10, text = history_label5.cget('text'), font=("Helvetica", f)); history_label4.grid(row=3, column=0)
        history_dict['h4'] = history_dict['h5']
        history_label5.grid_forget()
        history_label5 = Label(frame10, text = l3, font=("Helvetica", f)); history_label5.grid(row=4, column=0)
        history_dict['h5'] = d3
        
    bind()
    
def h4(*args):
    """
    Function updating the history labels & dictionnary when right clicking on the forth sequence in history calls the open history function
    Input: None
    Output: None
    """  
    global history_label1, history_label2, history_label3, history_label4, history_label5
    
    key = 'h4'
    h_open(key)
    
    l4 = history_label4.cget('text')
    d4 = history_dict['h4']
    
    if history_label5.cget('text') ==  " ... ":     
        history_label5 = Label(frame10, text = l4, font=("Helvetica", f)); history_label5.grid(row = 4, column = 0)
        history_dict['h5'] = d4
    else:    
        history_label1.grid_forget()
        history_label1 = Label(frame10, text = history_label2.cget('text'), font=("Helvetica", f)); history_label1.grid(row=0, column=0)
        history_dict['h1'] = history_dict['h2']
        history_label2.grid_forget()
        history_label2 = Label(frame10, text = history_label3.cget('text'), font=("Helvetica", f)); history_label2.grid(row=1, column=0)
        history_dict['h2'] = history_dict['h3']
        history_label3.grid_forget()
        history_label3 = Label(frame10, text = history_label4.cget('text'), font=("Helvetica", f)); history_label3.grid(row=2, column=0)
        history_dict['h3'] = history_dict['h4']
        history_label4.grid_forget()
        history_label4 = Label(frame10, text = history_label5.cget('text'), font=("Helvetica", f)); history_label4.grid(row=3, column=0)
        history_dict['h4'] = history_dict['h5']
        history_label5.grid_forget()
        history_label5 = Label(frame10, text = l4, font=("Helvetica", f)); history_label5.grid(row=4, column=0)
        history_dict['h5'] = d4
        
    bind()
    
def h5(*args):
    """
    Function updating the history labels & dictionnary when right clicking on the fifth sequence in history calls the open history function
    Input: None
    Output: None
    """  
    global history_label1, history_label2, history_label3, history_label4, history_label5
    
    key = 'h5'
    h_open(key)
    
    l5 = history_label5.cget('text')
    d5 = history_dict['h5']
    
    history_label1.grid_forget()
    history_label1 = Label(frame10, text = history_label2.cget('text'), font=("Helvetica", f)); history_label1.grid(row=0, column=0)
    history_dict['h1'] = history_dict['h2']
    history_label2.grid_forget()
    history_label2 = Label(frame10, text = history_label3.cget('text'), font=("Helvetica", f)); history_label2.grid(row=1, column=0)
    history_dict['h2'] = history_dict['h3']
    history_label3.grid_forget()
    history_label3 = Label(frame10, text = history_label4.cget('text'), font=("Helvetica", f)); history_label3.grid(row=2, column=0)
    history_dict['h3'] = history_dict['h4']
    history_label4.grid_forget()
    history_label4 = Label(frame10, text = history_label5.cget('text'), font=("Helvetica", f)); history_label4.grid(row=3, column=0)
    history_dict['h4'] = history_dict['h5']
    history_label5.grid_forget()
    history_label5 = Label(frame10, text = l5, font=("Helvetica", f)); history_label5.grid(row=4, column=0)
    history_dict['h5'] = d5
    
    bind()

def h_open(key):
    """
    Function opening a sequence in the history based on the key provided
    Input: key --> string {h1, h2, h3, h4, h5} representing which sequence the user wants to open from the history
    Output: None
    """  
    global Pre_def_seq
    global final_seq
    global final_im_ax
    global final_im_cor
    global final_im_sag
    global final_im_ax_original
    global final_im_cor_original
    global final_im_sag_original
    global TE_entry;         global TE_label
    global TI_entry;         global TI_label
    global TI2_entry;        global TI2_label
    global Alpha_entry;      global Alpha_label
    global G_entry;          global G_label
    global smalldelta_entry; global smalldelta_label
    global bigdelta_entry;   global bigdelta_label
    global post_TSE_ETL_entry
    global post_TSE_TE_entry
    global B_label
    global Bval_label
    global TEmin_label;      global minimumTE_label
    global TEminval_label;   
    global ETL_entry;        global ETL_label
    global Kspace_traj_label
    global traj
    global tse_drop
    global post_tse_drop
    global TEeff_label
    global TEeff_val_label
    global frame5 
    global seq_label2
    global warning_label
    global Time_scan_num_label 
    global Bd_by_pixel_label  
    global SNR_num_label
    global minimumTE_num_label
    global minimumTR_num_label
    global minimumTR
    global minimumTE
    global canvas1
    global canvas2
    global canvas3
    
    canvas1.get_tk_widget().destroy()
    canvas2.get_tk_widget().destroy()
    canvas3.get_tk_widget().destroy()

    X = history_dict[key]
    final_seq = X[0]
    Data_mat = X[1]
    FOV = X[2]
    seq_used = X[6] # The sequences that was used for those images
    Reso = X[7]
    
    Pre_def_seq = seq_used
    
    # The following lines until "S = [int..." are just showing thecorrect widgets and labels
    if seq_used == 'SE':
        name = 'Spin echo'
    elif seq_used == 'GE':
        name = 'Gradient echo'
    elif seq_used == 'IN':
        name = 'Inversion recovery'
    elif seq_used == 'Double IN':
        name = 'Double inversion recovery'
    elif seq_used == 'FLAIR':
        name = 'FLAIR'
    elif seq_used == 'SSFP':
        name = 'Steady-state'
    elif seq_used == 'Dif':
        name = 'Diffusion'
    elif seq_used == 'TSE':
        name = 'Turbo spin echo'
        
    seq_label2.grid_forget()
    seq_label2 = Label(frame1, text = name, font=("Helvetica", 18))
    seq_label2.grid(row = 10, column = 0, rowspan = 2)
        
    TE_entry.grid_forget();         TE_label.grid_forget()
    Alpha_entry.grid_forget();      Alpha_label.grid_forget()
    TI_entry.grid_forget();         TI_label.grid_forget()
    TI2_entry.grid_forget();        TI2_label.grid_forget()
    G_label.grid_forget();          G_entry.grid_forget()     
    smalldelta_label.grid_forget(); smalldelta_entry.grid_forget()
    bigdelta_label.grid_forget();   bigdelta_entry.grid_forget()    
    B_label.grid_forget()
    Bval_label.grid_forget()
    TEmin_label.grid_forget()   
    TEminval_label.grid_forget()    
    ETL_entry.grid_forget();        ETL_label.grid_forget()
    Kspace_traj_label.grid_forget();tse_drop.grid_forget()
    TEeff_label.grid_forget()
    TEeff_val_label.grid_forget()
    frame5.grid_forget()  
    
    show_widgets(seq_used)
    
    # show parameters used for an old sequence in the history
    Bandwidth_entry.delete(0, END)
    Bandwidth_entry.insert(0, str(X[14]))
    FOV1_entry.delete(0, END); FOV1_entry.insert(0, int(FOV[0]))
    FOV2_entry.delete(0, END); FOV2_entry.insert(0, int(FOV[1]))
    FOV3_entry.delete(0, END); FOV3_entry.insert(0, int(FOV[2]))
    Res1_entry.delete(0, END); Res1_entry.insert(0, float(Reso[0])) 
    Res2_entry.delete(0, END); Res2_entry.insert(0, float(Reso[1])) 
    Res3_entry.delete(0, END); Res3_entry.insert(0, float(Reso[2])) 
    Data_mat1_entry.delete(0,END); Data_mat1_entry.insert(0,int(Data_mat[0]))
    Data_mat2_entry.delete(0,END); Data_mat2_entry.insert(0,int(Data_mat[1]))
    Data_mat3_entry.delete(0,END); Data_mat3_entry.insert(0,int(Data_mat[2]))
    
    warning_label.grid_forget()
    Time_scan_num_label.grid_forget()
    Bd_by_pixel_label.grid_forget()
    minimumTE_num_label.grid_forget()
    minimumTR_num_label.grid_forget()
    
    Time_scan_num_label = Label(frame11, text = str(X[22]), font=("Helvetica", 12));        Time_scan_num_label.grid(row = 0, column = 1) 
    Bd_by_pixel_label = Label(frame11, text = str(round(X[23],2)), font=("Helvetica", 12)); Bd_by_pixel_label.grid(row = 1, column = 1)
    
    if seq_used == 'SSFP':
        minimumTE_label = Label(frame11, text = "TRmin (ms) ", font=("Helvetica", 12)).grid(row = 3, column = 0)
        minimumTR_num_label = Label(frame11, text = np.round(X[24] * 1000,2), font=("Helvetica", 12))
        minimumTR_num_label.grid(row = 3, column = 1)
    else:
        minimumTE_label = Label(frame11, text = "TEmin (ms)", font=("Helvetica", 12)).grid(row = 3, column = 0)
        minimumTE_num_label = Label(frame11, text = np.round(X[25] * 1000,2), font=("Helvetica", 12))
        minimumTE_num_label.grid(row = 3, column = 1)
    
    if seq_used == 'SE':
        TR_entry.delete(0,END);TR_entry.insert(0, str(X[8]))
        TE_entry.delete(0,END);TE_entry.insert(0, str(X[9]))
    elif seq_used == 'GE':
        TE_entry.delete(0,END);TE_entry.insert(0, str(X[9]))
        Alpha_entry.delete(0,END);Alpha_entry.insert(0, str(X[12]))
    elif seq_used == 'IN':
        post_traj.set(X[18])
        post_TSE_TE_entry.delete(0,END);post_TSE_TE_entry.insert(0,X[15])
        post_TSE_ETL_entry.delete(0,END);post_TSE_ETL_entry.insert(0, X[16])
        TE_entry.delete(0,END);TE_entry.insert(0, str(X[9]))
        TI_entry.delete(0,END);TI_entry.insert(0, str(X[10]))
    elif seq_used == 'Double IN':
        post_traj.set(X[18])
        post_TSE_TE_entry.delete(0,END);post_TSE_TE_entry.insert(0,X[15])
        post_TSE_ETL_entry.delete(0,END);post_TSE_ETL_entry.insert(0, X[16])
        TE_entry.delete(0,END);TE_entry.insert(0, str(X[9]))
        TI_entry.delete(0,END);TI_entry.insert(0, str(X[10]))
        TI2_entry.delete(0,END);TI2_entry.insert(0, str(X[11]))
    elif seq_used == 'FLAIR':
        post_traj.set(X[18])
        post_TSE_TE_entry.delete(0,END);post_TSE_TE_entry.insert(0,X[15])
        post_TSE_ETL_entry.delete(0,END);post_TSE_ETL_entry.insert(0, X[16])
        TR_entry.delete(0,END);TR_entry.insert(0, str(X[8]))
        TE_entry.delete(0,END);TE_entry.insert(0, str(X[9]))
        TI_entry.delete(0,END);TI_entry.insert(0, str(X[10]))
    elif seq_used == 'Dif':
        post_traj.set(X[18])
        post_TSE_TE_entry.delete(0,END);post_TSE_TE_entry.insert(0,X[15])
        post_TSE_ETL_entry.delete(0,END);post_TSE_ETL_entry.insert(0, X[16])
        TE_entry.delete(0,END);TE_entry.insert(0, str(X[9]))
        G_entry.delete(0,END);G_entry.insert(0, X[19])  
        smalldelta_entry.delete(0,END);smalldelta_entry.insert(0, X[20])
        bigdelta_entry.delete(0,END);bigdelta_entry.insert(0, X[21])
        Bval_label.grid_forget()
        Bval_label = Label(frame3, text = np.round(X[26],2), font=("Helvetica", 12));     Bval_label.grid(row = 8, column = 4)
        TEminval_label.grid_forget()
        TEminval_label = Label(frame3, text = np.round(X[27],2), font=("Helvetica", 12)); TEminval_label.grid(row = 9, column = 4)
    elif seq_used == 'TSE':
        TE_entry.delete(0,END);TE_entry.insert(0, str(X[9]))
        ETL_entry.delete(0,END);ETL_entry.insert(0, str(X[13]))
        traj.set(X[17])
    elif seq_used == 'SSFP':
        Alpha_entry.delete(0,END);Alpha_entry.insert(0, str(X[12]))  
    
    S = [int(Data_mat[0]/2),int(Data_mat[1]/2),int(Data_mat[2]/2)]   # S is a vector with the center of the Data matrix
    
    if type(final_seq) != list:
        final_im_sag = final_seq[S[0],:,:]
        final_im_cor = final_seq[:,S[1],:]
        final_im_ax = final_seq[:,:,S[2]]
    elif type(final_seq) == list: # non local mean algo takes to much time to run on full 3D data, so only on 2D
        final_im_sag = X[3]
        final_im_cor = X[4]
        final_im_ax = X[5]

    final_im_sag_original = final_im_sag
    final_im_cor_original = final_im_cor
    final_im_ax_original = final_im_ax
    
    # To change the axis value to mm (the fov)
    ax_pos = [[0, Data_mat[0]/2, Data_mat[0]-1], [0, Data_mat[1]/2, Data_mat[1]-1], [0, Data_mat[2]/2, Data_mat[2]-1]]
    ax_label = [[0, FOV[0]/2, FOV[0]], [0, FOV[1]/2, FOV[1]], [0, FOV[2]/2, FOV[2]]]
    Labels = ['Left/Right [mm]', 'Anterior/Posterior [mm]', 'Foot/Head [mm]']

    # Create axial plot
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.rot90(final_im_ax, k = 1), cmap='gray')  # Create image plot
    plt.xlabel(Labels[0]); plt.xticks(ax_pos[0], ax_label[0])
    plt.ylabel(Labels[1]); plt.yticks(ax_pos[1], ax_label[1])
    canvas1 = FigureCanvasTkAgg(fig, root)                  # Tkinter canvas which contains matplotlib figure
    canvas1.draw()
    canvas1.get_tk_widget().grid(row = 1, column = 3, rowspan = 4)       # Placing canvas on Tkinter window
    plt.close()

    # Create coronal plot
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.rot90(final_im_cor, k = 1), cmap='gray')   
    plt.xlabel(Labels[0]); plt.xticks(ax_pos[0], ax_label[0])
    plt.ylabel(Labels[2]); plt.yticks(ax_pos[2], ax_label[2])
    canvas2 = FigureCanvasTkAgg(fig, root)         
    canvas2.draw()
    canvas2.get_tk_widget().grid(row = 1, column = 2, rowspan = 4)
    plt.close()

    # Create sagittal plot
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.flip(np.rot90(final_im_sag, k = 1), axis = 1), cmap='gray')   
    plt.xlabel(Labels[1]); plt.xticks(ax_pos[1], ax_label[1])
    plt.ylabel(Labels[2]); plt.yticks(ax_pos[2], ax_label[2])
    canvas3 = FigureCanvasTkAgg(fig, root)         
    canvas3.draw()
    canvas3.get_tk_widget().grid(row = 1, column = 1, rowspan = 4) 
    plt.close() 
    
    m = SNR_length
    # std box
    A = int(SNR_noise_box_center[1])
    B = int(SNR_noise_box_center[0])
    # mean box
    a = int(SNR_mean_box_center[1])
    b = int(SNR_mean_box_center[0])

    snr = snr_homemade(np.rot90(final_im_ax_original),A-m,A+m,B-m,B+m,a-m,a+m,b-m,b+m)
    SNR_num_label.grid_forget()
    SNR_num_label = Label(frame11, text = str(np.abs(round(snr,2))), font=("Helvetica", 12))
    SNR_num_label.grid(row = 2, column = 1)

def history(se):
    """
    Function that will update the history frame and dictonnary when the user has simulated a new sequence or used a filter
    Input: se --> string {'seq','low','high','gauss','nonloc','undo'} designing if the function was called when a sequence ('seq') or a filter ('low','high','gauss','nonloc') was run or when undoing the filters actions ('undo') 
    Output: None
    """    
    global history_label1, history_label2, history_label3, history_label4, history_label5
    global Pre_def_seq
    global TR, TE, TI, TI2
    global ETL, Bandwidth, Alpha
    global ETL_entry
    global traj
    global post_TSE_ETL, post_TSE_TE
    global G        
    global smalldelta
    global bigdelta
    global B
    global dif_TEmin
    global low_pass_entry
    global high_pass_entry
    global Gauss_entry
    global Non_local_h_entry
    global Non_local_psize_entry
    global Non_local_pdist_entry
    global history_dict
    global history_list
    global final_seq
    global final_im_sag
    global final_im_cor
    global final_im_ax
    global FOV
    global Resolution
    global Data_mat 
    global Time_scan_num_label; global time     
    global Bd_by_pixel_label;   global Bd_by_pixel
    global SNR_num_label
    global minimumTE_num_label
    global minimumTR_num_label
    global minimumTR
    global minimumTE
    global Bval_label
    global TEminval_label
    
    # Defining the label message
    if se == 'seq' or se == 'undo':
        t = '{} TR/TE/TI/TI2 -- {}/{}/{}/{} '.format(Pre_def_seq, int(TR*1000), int(TE*1000), int(TI*1000), int(TI2*1000))
    elif se == 'low':
        
        if type(final_seq) != list: # The previous final seq doesn't come from the non local mean
            final_seq = cv2.blur(final_seq,(int(low_pass_entry.get()),int(low_pass_entry.get())))
            
        if type(final_seq) == list:# The previous final seq doesn't come from the non local mean
            final_seq[0] = cv2.blur(final_seq[0],(int(low_pass_entry.get()),int(low_pass_entry.get())))
            final_seq[1] = cv2.blur(final_seq[1],(int(low_pass_entry.get()),int(low_pass_entry.get())))
            final_seq[2] = cv2.blur(final_seq[2],(int(low_pass_entry.get()),int(low_pass_entry.get())))
        t = '{} Lowpass kernel size -- {}'.format(Pre_def_seq, int(low_pass_entry.get()))
        
    elif se == 'high':
        
        if type(final_seq) != list:
            final_seq = cv2.Laplacian(final_seq, -1, ksize = int(high_pass_entry.get()))

        if type(final_seq) == list:
            final_seq[0] = cv2.Laplacian(final_seq[0], -1, ksize = int(high_pass_entry.get()))
            final_seq[1] = cv2.Laplacian(final_seq[1], -1, ksize = int(high_pass_entry.get()))
            final_seq[2] = cv2.Laplacian(final_seq[2], -1, ksize = int(high_pass_entry.get()))
        t = '{} Highpass kernel size -- {} '.format(Pre_def_seq, int(high_pass_entry.get()))
        
    elif se == 'gauss':
        
        if type(final_seq) != list:
            final_seq = gaussian_filter(final_seq, sigma=float(Gauss_entry.get()))
            
        if type(final_seq) == list:
            final_seq[0] = gaussian_filter(final_seq[0], sigma=float(Gauss_entry.get()))
            final_seq[1] = gaussian_filter(final_seq[1], sigma=float(Gauss_entry.get()))
            final_seq[2] = gaussian_filter(final_seq[2], sigma=float(Gauss_entry.get()))
        t = '{} Gaussian pass std -- {}'.format(Pre_def_seq, float(Gauss_entry.get()))
        
    elif se == 'nonloc':
        
        H = float(Non_local_h_entry.get())
        S = int(Non_local_psize_entry.get())
        d = int(Non_local_pdist_entry.get())
        
        final_seq = [final_im_sag, final_im_cor, final_im_ax]
        
        t = '{} Non local mean h/s/d -- {}/{}/{}'.format(Pre_def_seq, H, S, d)
        
    if history_label1.cget('text') ==  " ... ":
        history_label1 = Label(frame10, text = t, font=("Helvetica", f)); history_label1.grid(row = 0, column = 0)
        history_dict['h1'] = [final_seq, Data_mat, FOV, final_im_sag, final_im_cor, final_im_ax, Pre_def_seq,Resolution,int(1000*TR),int(1000*TE),int(1000*TI),int(1000*TI2),int(Alpha),int(ETL),int(Bandwidth),int(1000*post_TSE_TE),int(post_TSE_ETL),traj.get(),post_traj.get(),int(1000*G),int(1000*smalldelta),int(1000*bigdelta),time,Bd_by_pixel,minimumTR,minimumTE,B,1000*dif_TEmin]
    elif history_label2.cget('text') ==  " ... ":     
        history_label2 = Label(frame10, text = t, font=("Helvetica", f)); history_label2.grid(row = 1, column = 0)
        history_dict['h2'] = [final_seq, Data_mat, FOV, final_im_sag, final_im_cor, final_im_ax, Pre_def_seq,Resolution,int(1000*TR),int(1000*TE),int(1000*TI),int(1000*TI2),int(Alpha),int(ETL),int(Bandwidth),int(1000*post_TSE_TE),int(post_TSE_ETL),traj.get(),post_traj.get(),int(1000*G),int(1000*smalldelta),int(1000*bigdelta),time,Bd_by_pixel,minimumTR,minimumTE,B,1000*dif_TEmin]
    elif history_label3.cget('text') ==  " ... ":     
        history_label3 = Label(frame10, text = t, font=("Helvetica", f)); history_label3.grid(row = 2, column = 0)
        history_dict['h3'] = [final_seq, Data_mat, FOV, final_im_sag, final_im_cor, final_im_ax, Pre_def_seq,Resolution,int(1000*TR),int(1000*TE),int(1000*TI),int(1000*TI2),int(Alpha),int(ETL),int(Bandwidth),int(1000*post_TSE_TE),int(post_TSE_ETL),traj.get(),post_traj.get(),int(1000*G),int(1000*smalldelta),int(1000*bigdelta),time,Bd_by_pixel,minimumTR,minimumTE,B,1000*dif_TEmin]
    elif history_label4.cget('text') ==  " ... ":     
        history_label4 = Label(frame10, text = t, font=("Helvetica", f)); history_label4.grid(row = 3, column = 0)
        history_dict['h4'] = [final_seq, Data_mat, FOV, final_im_sag, final_im_cor, final_im_ax, Pre_def_seq,Resolution,int(1000*TR),int(1000*TE),int(1000*TI),int(1000*TI2),int(Alpha),int(ETL),int(Bandwidth),int(1000*post_TSE_TE),int(post_TSE_ETL),traj.get(),post_traj.get(),int(1000*G),int(1000*smalldelta),int(1000*bigdelta),time,Bd_by_pixel,minimumTR,minimumTE,B,1000*dif_TEmin]
    elif history_label5.cget('text') ==  " ... ":     
        history_label5 = Label(frame10, text = t, font=("Helvetica", f)); history_label5.grid(row = 4, column = 0)
        history_dict['h5'] = [final_seq, Data_mat, FOV, final_im_sag, final_im_cor, final_im_ax, Pre_def_seq,Resolution,int(1000*TR),int(1000*TE),int(1000*TI),int(1000*TI2),int(Alpha),int(ETL),int(Bandwidth),int(1000*post_TSE_TE),int(post_TSE_ETL),traj.get(),post_traj.get(),int(1000*G),int(1000*smalldelta),int(1000*bigdelta),time,Bd_by_pixel,minimumTR,minimumTE,B,1000*dif_TEmin]
    else:
        history_label1.grid_forget()
        history_label1 = Label(frame10, text = history_label2.cget('text'), font=("Helvetica", f)); history_label1.grid(row=0, column=0)
        history_dict['h1'] = history_dict['h2']
        history_label2.grid_forget()
        history_label2 = Label(frame10, text = history_label3.cget('text'), font=("Helvetica", f)); history_label2.grid(row=1, column=0)
        history_dict['h2'] = history_dict['h3']
        history_label3.grid_forget()
        history_label3 = Label(frame10, text = history_label4.cget('text'), font=("Helvetica", f)); history_label3.grid(row=2, column=0)
        history_dict['h3'] = history_dict['h4']
        history_label4.grid_forget()
        history_label4 = Label(frame10, text = history_label5.cget('text'), font=("Helvetica", f)); history_label4.grid(row=3, column=0)
        history_dict['h4'] = history_dict['h5']
        history_label5.grid_forget()
        history_label5 = Label(frame10, text = t, font=("Helvetica", f)); history_label5.grid(row=4, column=0)
        history_dict['h5'] = [final_seq, Data_mat, FOV, final_im_sag, final_im_cor, final_im_ax, Pre_def_seq,Resolution,int(1000*TR),int(1000*TE),int(1000*TI),int(1000*TI2),int(Alpha),int(ETL),int(Bandwidth),int(1000*post_TSE_TE),int(post_TSE_ETL),traj.get(),post_traj.get(),int(1000*G),int(1000*smalldelta),int(1000*bigdelta),time,Bd_by_pixel,minimumTR,minimumTE,B,1000*dif_TEmin]
    bind()        

#/////////// Functions showing the effective TE for TSE and post-TSE sequences //////////////#      
def showTEeff(*args):
    """
    Computes and prints the effective TE for the TSE sequence when the k-space trajectory is changed
    Input : None
    Output: None
    """
    global ETL
    global TE
    global TEeff_val_label
    ETL = ETL_entry.get(); ETL = int(ETL) 
    TE = TE_entry.get(); TE = int(TE)
    trajectory = traj.get()
    eff = 0
    
    if trajectory == 'Linear':
        eff = 0.5 * TE * ETL
    elif trajectory == 'In-out':
        eff = TE
    elif trajectory == 'Out-in':
        eff = TE * ETL   
        
    TEeff_val_label.grid_forget()
    TEeff_val_label = Label(frame3, text = np.round(eff,2), font = ("Helvetica", 12)); TEeff_val_label.grid(row = 8, column = 3) 

def post_showTEeff(*args):
    """
    Computes and prints the effective TE for the post TSE sequence when the k-space trajectory is changed
    Input : None
    Output: None
    """
    global Post_TSE_TEeff_label_val
    global post_TSE_TE_entry
    global post_TSE_ETL_entry
    post_ETL = int(post_TSE_ETL_entry.get()) 
    post_TE = int(post_TSE_TE_entry.get())
    posttraj = post_traj.get()
    eff = 0
    
    if posttraj == 'Linear':
        eff = 0.5 * post_TE * post_ETL
    elif posttraj == 'In-out':
        eff = post_TE
    elif posttraj == 'Out-in':
        eff = post_TE * post_ETL   
        
    Post_TSE_TEeff_label_val.grid_forget()
    Post_TSE_TEeff_label_val = Label(frame3, text = np.round(eff,2), font = ("Helvetica", 12))
    Post_TSE_TEeff_label_val.grid(row = 15, column = 3) 

#/////////// Functions updating the sequence variable //////////////# 
def seq_simu(name):
    """
    Function informing the user of selected sequence and updates the widgets (labels and entries) for the specific sequence 
    Input : seq --> string of Pre_def_seq
    Output: None
    """
    global Pre_def_seq
    global seq_label2
    
    remove_widgets(Pre_def_seq)
    
    if name == 'Spin echo':
        Pre_def_seq = 'SE'
    elif name == 'Gradient echo':
        Pre_def_seq = 'GE'
    elif name == 'Inversion recovery':
        Pre_def_seq = 'IN'
    elif name == 'Double inversion recovery':
        Pre_def_seq = 'Double IN'
    elif name == 'Fluid-attenuated inversion recovery':
        Pre_def_seq = 'FLAIR'
    elif name == 'Steady-state':
        Pre_def_seq = 'SSFP'
    elif name == 'Diffusion':
        Pre_def_seq = 'Dif'
    elif name == 'Turbo spin echo':
        Pre_def_seq = 'TSE'
        
    seq_label2.grid_forget()
    if Pre_def_seq == 'FLAIR':
        seq_label2 = Label(frame1, text = 'FLAIR', font=("Helvetica", 18))
    else:
        seq_label2 = Label(frame1, text = name, font=("Helvetica", 18))
        
    seq_label2.grid(row = 10, column = 0, rowspan = 2)
    show_widgets(Pre_def_seq)
    
#/////////// Function computing the TE effective base don the kspace trajectory and multiplying the data by the exponential of it //////////////#     
def post_TSE_TEeff(data, postTE, postETL, T2, traj):
    """
    Function performing an exponential multiplication of the data based on the appropriate kspace trajectory
    Inputs:
    data    -> simulated sequence
    postTE  -> TE used for the post sequence TSE acquisition
    postETL -> ETL used for the post sequence TSE acquisition
    T2      -> T2 relaxation map
    traj    -> kspace trajectory
    Ouputs:
    data    -> simulated sequence multiplied by the exponentional
    """
    if traj == 'Out-in':    
        TEeff = postTE * postETL    
        TE_div = np.divide(-TEeff, T2, out=np.zeros_like(T2), where=T2!=0)
        data = data*np.exp(TE_div)
    elif traj == 'In-out':
        TEeff = postTE
        TE_div = np.divide(-TEeff, T2, out=np.zeros_like(T2), where=T2!=0)
        data = data*np.exp(TE_div)
    elif traj == 'Linear':
        TEeff = 0.5 * postTE * postETL
        TE_div = np.divide(-TEeff, T2, out=np.zeros_like(T2), where=T2!=0)
        data = data*np.exp(TE_div)     
    return data

#/////////// Function resizing the data to the correct dimension //////////////#   
def resize(data, T1, data_size):
    """
    Function resizing the data to the appropriate final shape defined by the user
    Inputs:
    data       -> simulated sequence
    T1         -> T1 relaxation map
    data_size  -> shape of the output data 
    Ouputs:
    data_out   -> simulated sequence reshaped
    """
    # Resizing the data
    n_seq = np.zeros((T1.shape[0], data_size[1], data_size[2]))
    data_out = np.zeros((data_size))
    for x in range(T1.shape[0]):
        n_seq[x,:,:] = cv2.resize(data[x,:,:], dsize=(data_size[2], data_size[1]))
    for x in range(data_size[1]):
        data_out[:,x,:] = cv2.resize(n_seq[:,x,:], dsize=(data_size[2], data_size[0]))   
    del n_seq
    return data_out
        
#/////////// Function running the simulation of the sequences //////////////#   
def run():   
    """
    Compute a 3D simulation of the sequence defined be the user
    Input : None
    Output: None
    """
    global imag_size
    global Pre_def_seq
    global TR, TE, TI, TI2, ETL, Alpha
    global post_TSE_ETL, post_TSE_TE
    global FOV
    global Data_mat
    global Resolution
    global Bandwidth
    global Labels
    global ax_pos
    global ax_label
    global final_seq
    global S
    global G        
    global smalldelta
    global bigdelta
    global B
    global dif_TEmin
    global final_im_sag
    global final_im_cor
    global final_im_ax
    global final_im_sag_original
    global final_im_cor_original
    global final_im_ax_original
    global minimumTE
    global minimumTR
    global Readout_axis
    global B0_3D
    global B0_3D_ori
    global Time_scan_num_label
    global Bd_by_pixel_label
    global SNR_num_label
    global minimumTE_num_label
    global minimumTR_num_label
    global warning_label
    global history_label1 
    global history_label2 
    global history_label3 
    global history_label4 
    global history_label5 
    global B0_field
    global canvas1
    global canvas2
    global canvas3
    global time 
    global Bd_by_pixel
    global SNR_noise_box_center
    global SNR_mean_box_center
    
    canvas1.get_tk_widget().destroy()
    canvas2.get_tk_widget().destroy()
    canvas3.get_tk_widget().destroy()
    
    TR = TR_entry.get(); TR = int(TR); TR = np.divide(TR,1000) # Divide by 1000 to have the values in milli
    TE = TE_entry.get(); TE = int(TE); TE = np.divide(TE,1000)
    TI = TI_entry.get(); TI = int(TI); TI = np.divide(TI,1000)
    TI2 = TI2_entry.get(); TI2 = int(TI2); TI2 = np.divide(TI2,1000)
    fov1 = FOV1_entry.get(); fov1 = int(fov1)
    fov2 = FOV2_entry.get(); fov2 = int(fov2)
    fov3 = FOV3_entry.get(); fov3 = int(fov3)
    res1 = Res1_entry.get(); res1 = float(res1) 
    res2 = Res2_entry.get(); res2 = float(res2)
    res3 = Res3_entry.get(); res3 = float(res3)
    bd = Bandwidth_entry.get(); Bandwidth = int(bd)
    Alpha = Alpha_entry.get(); Alpha = int(Alpha)
    ETL = ETL_entry.get(); ETL = int(ETL)
    FOV = [fov1, fov2, fov3]
    Resolution = [res1, res2, res3]
    Data_mat = np.divide(FOV, Resolution)
    Data_mat = [int(Data_mat[0]), int(Data_mat[1]), int(Data_mat[2])]
    S = [int(Data_mat[0]/2),int(Data_mat[1]/2),int(Data_mat[2]/2)]   # S is a vector with the center of the Data matrix
    
    # For diffusion
    G = G_entry.get(); G = float(G); G = np.divide(G,1000)
    smalldelta = smalldelta_entry.get(); smalldelta = float(smalldelta); smalldelta = np.divide(smalldelta,1000)
    bigdelta = bigdelta_entry.get(); bigdelta = float(bigdelta); bigdelta = np.divide(bigdelta,1000)
    
    # Post sequence TSE entries
    post_TSE_TE = post_TSE_TE_entry.get(); post_TSE_TE = int(post_TSE_TE); post_TSE_TE = np.divide(post_TSE_TE,1000)
    post_TSE_ETL = post_TSE_ETL_entry.get(); post_TSE_ETL = int(post_TSE_ETL)
    
    Time_scan_num_label.grid_forget()
    Bd_by_pixel_label.grid_forget()
    SNR_num_label.grid_forget()
    minimumTE_num_label.grid_forget()
    minimumTR_num_label.grid_forget()
    
    minimumTE = Data_mat[0]/(2*Bandwidth)
    minimumTR = 2*minimumTE
        
    final_zeros = np.zeros((Data_mat))
    
    # If ok --> 'yes' the sequence will be simulated, if ok --> 'no', there is a problem with the parameters and the images will be completly black
    ok = 'yes'
    
    # Checking if the parameters make physical sense
    if Pre_def_seq == 'SSFP':
        TE = 0 # Because the original (inserted) value of TE can be higher than a desired TR, but TE isn't used in the SSFP seq
     
        if (TR) < minimumTR:  
            print('TR can not be smaller than the minimum TR!!!')
            ok = 'no' 
            w1 = "Warning!!"
            w2 = "Error, TR < minTR!!!"
            
    elif Pre_def_seq == 'IN':
        if TR < (TE + TI):
            print('TR can not be smaller than TE + TI!!!')
            ok = 'no'
            w1 = "Warning!!"
            w2 = "TR smaller than TE + TI!!!"

    elif Pre_def_seq == 'Double IN':
        if TR < (TE + TI + TI2):
            print('TR can not be smaller than TE + TI + TI2!!!')
            ok = 'no'
            w1 = "Warning!!"
            w2 = "TR smaller than TE + TI + TI2!!!"
    
    else:
        if TR < TE:
            print('TE can not be larger than TR!!!')
            ok = 'no'
            w1 = "Warning!!"
            w2 = "Error, TR < TE!!!"

        if (TE) < minimumTE:
            print('TE can not be smaller than the minimum TE!!!')
            ok = 'no'
            w1 = "Warning!!"
            w2 = "Error, TE < TEmin!!!"

        if B != 0:
            if Pre_def_seq == 'Dif':
                if dif_TEmin > TE:
                    print('Error, TE must be grater than diffusion TEmin!!!')
                    ok = 'no'
                    w1 = "Warning!!"
                    w2 = "Error, TE < diffusion TEmin !!!"
                
    if ok == 'no': # There is a problem with the parameters and the images will be completly black
        
        warning_label.grid_forget()
        warning_label = Label(root, text = w2, font=("Helvetica", 15), fg='#f00'); warning_label.grid(row = 8, column = 0)
        
        # Create axial plot
        fig = plt.figure(figsize=(imag_size,imag_size))
        plt.imshow(np.rot90(final_zeros[:,:,0], k = 1), cmap='gray')  # Create image plot
        plt.axis('off')
        canvas1 = FigureCanvasTkAgg(fig, root)                         # Tkinter canvas which contains matplotlib figure
        canvas1.draw()
        canvas1.get_tk_widget().grid(row = 1, column = 3, rowspan = 4) # Placing canvas on Tkinter window
        plt.close()
        
        # Create coronal plot
        fig = plt.figure(figsize=(imag_size,imag_size))
        plt.imshow(np.rot90(final_zeros[:,0,:], k = 1), cmap='gray')   
        plt.axis('off')
        canvas2 = FigureCanvasTkAgg(fig, root)         
        canvas2.draw()
        canvas2.get_tk_widget().grid(row = 1, column = 2, rowspan = 4)
        plt.close()
        
        # Create sagittal plot
        fig = plt.figure(figsize=(imag_size,imag_size))
        plt.imshow(np.flip(np.rot90(final_zeros[0,:,:], k = 1), axis = 1), cmap='gray')   
        plt.axis('off')
        canvas3 = FigureCanvasTkAgg(fig, root)         
        canvas3.draw()
        canvas3.get_tk_widget().grid(row = 1, column = 1, rowspan = 4) 
        plt.close()
                
    else: # The sequence will be simulated
        
        warning_label.grid_forget()
        
        if Pre_def_seq == 'SSFP':
            minimumTE_label = Label(frame11, text = "TRmin (ms) ", font=("Helvetica", 12)).grid(row = 3, column = 0)
            minimumTR_num_label = Label(frame11, text = np.round(minimumTR * 1000,2), font=("Helvetica", 12))
            minimumTR_num_label.grid(row = 3, column = 1)
        else:
            minimumTE_label = Label(frame11, text = "TEmin (ms)", font=("Helvetica", 12)).grid(row = 3, column = 0)
            minimumTE_num_label = Label(frame11, text = np.round(minimumTE * 1000,2), font=("Helvetica", 12))
            minimumTE_num_label.grid(row = 3, column = 1)
        
        time_scan = TR * Data_mat[1] * Data_mat[2]
        
        if Pre_def_seq == 'TSE':
            ETL = ETL_entry.get(); ETL = int(ETL)
            time_scan = time_scan/ETL

        time = datetime.timedelta(seconds=time_scan) # Converts amount of seconds to hours:minutes:seconds
        Bd_by_pixel = np.divide(Bandwidth,Data_mat[0])
        
        # Plotting the parameters
        Time_scan_num_label = Label(frame11, text = str(time),font=("Helvetica",12))
        Time_scan_num_label.grid(row = 0, column = 1)
        Bd_by_pixel_label = Label(frame11, text = str(round(Bd_by_pixel,2)), font=("Helvetica", 12))
        Bd_by_pixel_label.grid(row = 1, column = 1)

        b_field = np.copy(B0_3D)
        
        # loading the 3D distorted data from the appropriate B0 strength
        if B0_field == 'B0_46mT':       
            T1_3D_grad = np.load(resource_path("Data\T1_3D_cer_lip_grad.npy"))
            T2_3D_grad = np.load(resource_path("Data\T2_3D_cer_lip_grad.npy"))
            M0_3D_grad = np.load(resource_path("Data\M0_3D_cer_lip_grad.npy"))
            B1map_3D_grad = np.load(resource_path("Data\B1_3D_cer_lip_grad.npy"))
            flipAngleMaprescale_3D_grad = np.load(resource_path("Data\FlipAngleMaprescale_3D_cer_lip_grad.npy"))
            t2_star_3D_grad = np.load(resource_path("Data\T2_star_tensor_3D_cer_lip_grad.npy"))
            ADC_3D_grad = np.load(resource_path("Data\ADC_3D_cer_lip_grad.npy"))     
            B0_3D = B0_3D_ori
        elif B0_field == 'B0_15T':
            T1_3D_grad = np.load(resource_path("Data\T1_15.npy"))
            T2_3D_grad = np.load(resource_path("Data\T2_15.npy"))
            M0_3D_grad = np.load(resource_path("Data\M0_3D_cerb_lip.npy"))
            B1map_3D_grad = np.load(resource_path("Data\B1_3D.npy"))
            flipAngleMaprescale_3D_grad = np.load(resource_path("Data\FlipAngleMaprescale_3D.npy"))
            t2_star_3D_grad = np.load(resource_path("Data\T2_star_15.npy"))
            ADC_3D_grad = np.load(resource_path("Data\ADC_3D.npy"))
            B0_3D = b_field * 32.608 # 1.5T / 46mT ~= 32.608
            B0_3D[B0_3D>0] = 1.5
        elif B0_field == 'B0_3T':
            T1_3D_grad = np.load(resource_path("Data\T1_3.npy"))
            T2_3D_grad = np.load(resource_path("Data\T2_3.npy"))
            M0_3D_grad = np.load(resource_path("Data\M0_3D_cerb_lip.npy"))
            B1map_3D_grad = np.load(resource_path("Data\B1_3D.npy"))
            flipAngleMaprescale_3D_grad = np.load(resource_path("Data\FlipAngleMaprescale_3D.npy"))
            t2_star_3D_grad = np.load(resource_path("Data\T2_star_3.npy"))
            ADC_3D_grad = np.load(resource_path("Data\ADC_3D.npy"))
            B0_3D = b_field * 65.217  # 3T / 46mT ~= 65.217
            B0_3D[B0_3D>0] = 3
            
        # Check if maps have been imported (If n_import is a list/array --> return True, if n_import is an integer --> return False)
        if isinstance(T1_import,np.ndarray): 
            T1_3D_grad = T1_import
        if isinstance(T2_import,np.ndarray): 
            T2_3D_grad = T2_import
        if isinstance(T2_star_import,np.ndarray): 
            t2_star_3D_grad = T2_star_import
        if isinstance(B0_import,np.ndarray): 
            B0_3D = B0_import
          
        ## change the data to match the FOV specified by the user ##
        # The number or 2D arrays to delete from each axis (at the beginning and end, for each axis)
        to_delete = [int((T1_3D_grad.shape[0] - FOV[0])/2), int((T1_3D_grad.shape[1] - FOV[1])/2), int((T1_3D_grad.shape[2] - FOV[2])/2)]
        # take a sub sample of the data corresponding to the FOV (from the center)
        sh = T1_3D_grad.shape
        T1_3D_grad = T1_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
        T2_3D_grad = T2_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
        M0_3D_grad = M0_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
        B1map_3D_grad = B1map_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
        flipAngleMaprescale_3D_grad = flipAngleMaprescale_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
        t2_star_3D_grad = t2_star_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
        ADC_3D_grad = ADC_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
            
        final_seq = np.zeros((Data_mat))
        
        if Pre_def_seq == 'SE':
            SE_3D = spin_echo_seq(TR, TE, T1_3D_grad, T2_3D_grad, M0_3D_grad)
            SE_3D = np.multiply(SE_3D, B1map_3D_grad) 
            final_seq = resize(SE_3D, T1_3D_grad, Data_mat)

        elif Pre_def_seq == 'GE':         
            angle = flipAngleMaprescale_3D_grad/Alpha
            GE_3D = Gradient_seq(TR, TE, T1_3D_grad, t2_star_3D_grad, M0_3D_grad, angle)
            GE_3D = np.multiply(GE_3D, B1map_3D_grad)
            final_seq = resize(GE_3D, T1_3D_grad, Data_mat)

        elif Pre_def_seq == 'IN': 
            IN_3D = IN_seq(TR, TI, T1_3D_grad, T2_3D_grad, M0_3D_grad)
            IN_3D = np.multiply(IN_3D, B1map_3D_grad)      
            trajectory = post_traj.get() 
            IN_3D = post_TSE_TEeff(IN_3D, post_TSE_TE, post_TSE_ETL, T2_3D_grad, trajectory)
            final_seq = resize(IN_3D, T1_3D_grad, Data_mat)
            final_seq = TSE_seq_R(post_TSE_TE, post_TSE_ETL, T2_3D_grad, trajectory, Readout_axis, final_seq)
                
        elif Pre_def_seq == 'Double IN':
            DIN_3D = DoubleInversion_seq(TR, post_TSE_TE, TI, TI2, T1_3D_grad, T2_3D_grad, M0_3D_grad)
            DIN_3D = np.multiply(DIN_3D, B1map_3D_grad)  
            trajectory = post_traj.get()     
            DIN_3D = post_TSE_TEeff(DIN_3D, post_TSE_TE, post_TSE_ETL, T2_3D_grad, trajectory)                       
            final_seq = resize(DIN_3D, T1_3D_grad, Data_mat)                
            final_seq = TSE_seq_R(post_TSE_TE, post_TSE_ETL, T2_3D_grad, trajectory, Readout_axis, final_seq)
            
        elif Pre_def_seq == 'FLAIR':
            FLAIR_3D = IN_seq(TR, TI, T1_3D_grad, T2_3D_grad, M0_3D_grad)
            FLAIR_3D = np.multiply(FLAIR_3D, B1map_3D_grad)
            trajectory = post_traj.get()                 
            FLAIR_3D = post_TSE_TEeff(FLAIR_3D, post_TSE_TE, post_TSE_ETL, T2_3D_grad, trajectory)            
            final_seq = resize(FLAIR_3D, T1_3D_grad, Data_mat)                
            final_seq = TSE_seq_R(post_TSE_TE, post_TSE_ETL, T2_3D_grad, trajectory, Readout_axis, final_seq)
            
        elif Pre_def_seq == 'Dif':
            DIF_3D = Diffusion_seq(TR, T1_3D_grad, T2_3D_grad, M0_3D_grad, B, ADC_3D_grad)
            DIF_3D = np.multiply(DIF_3D, B1map_3D_grad)   
            trajectory = post_traj.get()               
            DIF_3D = post_TSE_TEeff(DIF_3D, post_TSE_TE, post_TSE_ETL, T2_3D_grad, trajectory)           
            final_seq = resize(DIF_3D, T1_3D_grad, Data_mat)                
            final_seq = TSE_seq_R(post_TSE_TE, post_TSE_ETL, T2_3D_grad, trajectory, Readout_axis, final_seq)
            
        elif Pre_def_seq == 'TSE':
            trajectory = traj.get()
            TSE_3D = 0
            TEeff = 0
            if trajectory == 'Out-in':    
                TEeff = TE * ETL    
                TSE_3D = spin_echo_seq(TR, TEeff, T1_3D_grad, T2_3D_grad, M0_3D_grad) # Computing a spin echo
            elif trajectory == 'In-out':
                TEeff = TE
                TE_div = np.divide(TEeff, T2_3D_grad, out=np.zeros_like(T2_3D_grad), where=T2_3D_grad!=0)
                TSE_3D = np.abs(M0_3D_grad * np.exp(-TE_div))                         # Computing a spin echo with the (1 - exp(-TR/T1)) set to 1
            elif trajectory == 'Linear':
                TEeff = 0.5 * TE * ETL
                TSE_3D = spin_echo_seq(TR, TEeff, T1_3D_grad, T2_3D_grad, M0_3D_grad) # Computing a spin echo
             
            TSE_3D = np.multiply(TSE_3D, B1map_3D_grad)          
            final_seq = resize(TSE_3D, T1_3D_grad, Data_mat)            
            final_seq = TSE_seq_R(TE, ETL, T2_3D_grad, trajectory, Readout_axis, final_seq)
            
        elif Pre_def_seq == 'SSFP':        
            gamma = 42.58*(10**6)*2*constants.pi # 267538030.37970677 [rad/sT]
            omega = B0_3D * gamma
            center = np.divide(omega.shape,2).astype(int)
            center_freq_value = omega[center[0],center[1],center[2]]
            offset = omega - center_freq_value
            phi = offset * TR
            
            if B0_field == 'B0_46mT':
                SSFP_3D = SSFP_Echo_seq(T1_3D_grad, T2_3D_grad, M0_3D_grad, Alpha, phi)
            elif B0_field == 'B0_3T' or 'B0_15T':
                SSFP_3D = SSFP_new(TR, T1_3D_grad, T2_3D_grad, M0_3D_grad, Alpha, phi)          
            
            SSFP_3D = np.multiply(SSFP_3D, B1map_3D_grad)       
            final_seq = resize(SSFP_3D, T1_3D_grad, Data_mat)
        
        # Noise computation     
        n = noise_generation(Data_mat, Resolution, Bandwidth, B0_field)
        final_seq = final_seq + n        
        final_im_sag = final_seq[S[0],:,:]
        final_im_cor = final_seq[:,S[1],:]
        final_im_ax = final_seq[:,:,S[2]]
        final_im_sag_original = final_im_sag
        final_im_cor_original = final_im_cor
        final_im_ax_original = final_im_ax  

        # To change the axis value to mm (the fov)
        ax_pos = [[0, Data_mat[0]/2, Data_mat[0]-1], [0, Data_mat[1]/2, Data_mat[1]-1], [0, Data_mat[2]/2, Data_mat[2]-1]]
        ax_label = [[0, FOV[0]/2, FOV[0]], [0, FOV[1]/2, FOV[1]], [0, FOV[2]/2, FOV[2]]]
        Labels = ['Left/Right [mm]', 'Anterior/Posterior [mm]', 'Foot/Head [mm]']

        # Create axial plot
        fig = plt.figure(figsize=(imag_size,imag_size))
        plt.imshow(np.rot90(final_im_ax, k = 1), cmap='gray')  # Create image plot
        plt.xlabel(Labels[0]); plt.xticks(ax_pos[0], ax_label[0])
        plt.ylabel(Labels[1]); plt.yticks(ax_pos[1], ax_label[1])
        canvas1 = FigureCanvasTkAgg(fig, root)                  # Tkinter canvas which contains matplotlib figure
        canvas1.draw()
        canvas1.get_tk_widget().grid(row = 1, column = 3, rowspan = 4)       # Placing canvas on Tkinter window
        plt.close()
        
        # Create coronal plot
        fig = plt.figure(figsize=(imag_size,imag_size))
        plt.imshow(np.rot90(final_im_cor, k = 1), cmap='gray')   
        plt.xlabel(Labels[0]); plt.xticks(ax_pos[0], ax_label[0])
        plt.ylabel(Labels[2]); plt.yticks(ax_pos[2], ax_label[2])
        canvas2 = FigureCanvasTkAgg(fig, root)         
        canvas2.draw()
        canvas2.get_tk_widget().grid(row = 1, column = 2, rowspan = 4)
        plt.close()
        
        # Create sagittal plot
        fig = plt.figure(figsize=(imag_size,imag_size))
        plt.imshow(np.flip(np.rot90(final_im_sag, k = 1), axis = 1), cmap='gray')   
        plt.xlabel(Labels[1]); plt.xticks(ax_pos[1], ax_label[1])
        plt.ylabel(Labels[2]); plt.yticks(ax_pos[2], ax_label[2])
        canvas3 = FigureCanvasTkAgg(fig, root)         
        canvas3.draw()
        canvas3.get_tk_widget().grid(row = 1, column = 1, rowspan = 4) 
        plt.close()
        
        # Check if the boxes used for the SNR are in the image, if one is outside the position is redefined
        if SNR_noise_box_center[0] > Data_mat[0] or SNR_noise_box_center[1] > Data_mat[1]:
            SNR_noise_box_center = [(Data_mat[0]-15),15]
        if SNR_mean_box_center[0] > Data_mat[0] or SNR_mean_box_center[1] > Data_mat[1]:
            SNR_mean_box_center = [(Data_mat[0]//5)*2,(Data_mat[1]//6)*4] # To have the mean box inside the brain
                     
        m = SNR_length
        # std box
        A = int(SNR_noise_box_center[1])
        B = int(SNR_noise_box_center[0])
        # mean box
        a = int(SNR_mean_box_center[1])
        b = int(SNR_mean_box_center[0])
        snr = snr_homemade(np.rot90(final_im_ax),A-m,A+m,B-m,B+m,a-m,a+m,b-m,b+m)
        SNR_num_label.grid_forget()
        SNR_num_label = Label(frame11, text = str(np.abs(round(snr,2))), font=("Helvetica", 12))
        SNR_num_label.grid(row = 2, column = 1)
            
        # Update the history frame
        history('seq')
         
#/////////// Function opening a parameter interactive window //////////////#  
def param_vis():
    """
    Opens a new window where the user can interact with the simulated data in 3D and change it's parameters to see their effect in real time
    Input : None
    Output: None
    """
    global Pre_def_seq
    global TR
    global TE
    global TI
    global TI2
    global FOV
    global Data_mat
    global Resolution
    global Bandwidth
    global B
    global ETL
    global B0_field
    global Readout_axis
    global B0_3D
    global B0_import
    global T2_star_import
    
    TR = TR_entry.get(); TR = int(TR); TR = np.divide(TR,1000) # Divide by 1000 to have the values in milli
    TE = TE_entry.get(); TE = int(TE); TE = np.divide(TE,1000)
    TI = TI_entry.get(); TI = int(TI); TI = np.divide(TI,1000)
    TI2 = TI2_entry.get(); TI2 = int(TI2); TI2 = np.divide(TI2,1000)
    fov1 = FOV1_entry.get(); fov1 = int(fov1)
    fov2 = FOV2_entry.get(); fov2 = int(fov2)
    fov3 = FOV3_entry.get(); fov3 = int(fov3)
    res1 = Res1_entry.get(); res1 = float(res1) 
    res2 = Res2_entry.get(); res2 = float(res2)
    res3 = Res3_entry.get(); res3 = float(res3)
    bd = Bandwidth_entry.get(); Bandwidth = int(bd)
    Alpha = Alpha_entry.get(); Alpha = int(Alpha)
    ETL = ETL_entry.get(); ETL = int(ETL) 
    FOV = [fov1, fov2, fov3]
    Resolution = [res1, res2, res3]
    Data_mat = np.divide(FOV, Resolution)
    Data_mat = [int(Data_mat[0]), int(Data_mat[1]), int(Data_mat[2])]    
    
    # loading the data from the appropriate B0 strength       
    if B0_field == 'B0_46mT':       
        T1_3D_grad = np.load(resource_path("Data\T1_3D_cer_lip_grad.npy"))
        T2_3D_grad = np.load(resource_path("Data\T2_3D_cer_lip_grad.npy"))
        M0_3D_grad = np.load(resource_path("Data\M0_3D_cer_lip_grad.npy"))
        B1map_3D_grad = np.load(resource_path("Data\B1_3D_cer_lip_grad.npy"))
        flipAngleMaprescale_3D_grad = np.load(resource_path("Data\FlipAngleMaprescale_3D_cer_lip_grad.npy"))
        t2_star_3D_grad = np.load(resource_path("Data\T2_star_tensor_3D_cer_lip_grad.npy"))
        ADC_3D_grad = np.load(resource_path("Data\ADC_3D_cer_lip_grad.npy"))
        B0_3D = B0_3D_ori
    elif B0_field == 'B0_15T':
        T1_3D_grad = np.load(resource_path("Data\T1_15.npy"))
        T2_3D_grad = np.load(resource_path("Data\T2_15.npy"))
        M0_3D_grad = np.load(resource_path("Data\M0_3D_cerb_lip.npy"))
        B1map_3D_grad = np.load(resource_path("Data\B1_3D.npy"))
        flipAngleMaprescale_3D_grad = np.load(resource_path("Data\FlipAngleMaprescale_3D.npy"))
        t2_star_3D_grad = np.load(resource_path("Data\T2_star_15.npy"))
        ADC_3D_grad = np.load(resource_path("Data\ADC_3D.npy"))
        B0_3D = B0_3D * 32.608 # 1.5T / 46mT ~= 32.608
        #B0_3D[B0_3D>0] = 1.5
    elif B0_field == 'B0_3T':
        T1_3D_grad = np.load(resource_path("Data\T1_3.npy"))
        T2_3D_grad = np.load(resource_path("Data\T2_3.npy"))
        M0_3D_grad = np.load(resource_path("Data\M0_3D_cerb_lip.npy"))
        B1map_3D_grad = np.load(resource_path("Data\B1_3D.npy"))
        flipAngleMaprescale_3D_grad = np.load(resource_path("Data\FlipAngleMaprescale_3D.npy"))
        t2_star_3D_grad = np.load(resource_path("Data\T2_star_3.npy"))
        ADC_3D_grad = np.load(resource_path("Data\ADC_3D.npy"))
        B0_3D = B0_3D * 65.217  # 3T / 46mT ~= 65.217
        #B0_3D[B0_3D>0] = 3
        
    if isinstance(T1_import,np.ndarray): 
        T1_3D_grad = T1_import
    if isinstance(T2_import,np.ndarray): 
        T2_3D_grad = T2_import
    if isinstance(T2_star_import,np.ndarray): 
        t2_star_3D_grad = T2_star_import
    if isinstance(B0_import,np.ndarray): 
        B0_3D = B0_import
        
    ## Change the data to match the FOV specified by the user ##
    # The number or 2D arrays to delete from each axis (at the beginning and end, for each axis)
    to_delete = [int((T1_3D_grad.shape[0] - FOV[0])/2), int((T1_3D_grad.shape[1] - FOV[1])/2), int((T1_3D_grad.shape[2] - FOV[2])/2)]
    # take a sub sample of the data corresponding to the FOV (from the center)
    sh = T1_3D_grad.shape
    T1_3D_grad = T1_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
    T2_3D_grad = T2_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
    M0_3D_grad = M0_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
    B1map_3D_grad = B1map_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
    flipAngleMaprescale_3D_grad = flipAngleMaprescale_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
    t2_star_3D_grad = t2_star_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
    ADC_3D_grad = ADC_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
 
    # Resizing
    n_t1 = np.zeros((T1_3D_grad.shape[0], Data_mat[1], Data_mat[2])); final_t1 = np.zeros((Data_mat))
    n_t2 = np.zeros((T1_3D_grad.shape[0], Data_mat[1], Data_mat[2])); final_t2 = np.zeros((Data_mat))
    n_m0 = np.zeros((T1_3D_grad.shape[0], Data_mat[1], Data_mat[2])); final_m0 = np.zeros((Data_mat))
    n_B1 = np.zeros((T1_3D_grad.shape[0], Data_mat[1], Data_mat[2])); final_B1 = np.zeros((Data_mat))
    n_t2star = np.zeros((T1_3D_grad.shape[0], Data_mat[1], Data_mat[2]));  final_t2star = np.zeros((Data_mat))
    n_flipmap = np.zeros((T1_3D_grad.shape[0], Data_mat[1], Data_mat[2])); final_flipmap = np.zeros((Data_mat))
    n_adc = np.zeros((T1_3D_grad.shape[0], Data_mat[1], Data_mat[2])); final_adc = np.zeros((Data_mat))
    n_seq = np.zeros((T1_3D_grad.shape[0], Data_mat[1], Data_mat[2])); final_seq = np.zeros((Data_mat))
    n_phi = np.zeros((T1_3D_grad.shape[0], Data_mat[1], Data_mat[2])); final_phi = np.zeros((Data_mat))
    n_offset = np.zeros((T1_3D_grad.shape[0], Data_mat[1], Data_mat[2])); final_offset = np.zeros((Data_mat))

    for x in range(T1_3D_grad.shape[0]):
        n_t1[x,:,:] = cv2.resize(T1_3D_grad[x,:,:], dsize=(Data_mat[2], Data_mat[1]))
        n_t2[x,:,:] = cv2.resize(T2_3D_grad[x,:,:], dsize=(Data_mat[2], Data_mat[1]))
        n_m0[x,:,:] = cv2.resize(M0_3D_grad[x,:,:], dsize=(Data_mat[2], Data_mat[1]))
        n_B1[x,:,:] = cv2.resize(B1map_3D_grad[x,:,:], dsize=(Data_mat[2], Data_mat[1]))
        n_t2star[x,:,:] = cv2.resize(t2_star_3D_grad[x,:,:], dsize=(Data_mat[2], Data_mat[1]))
        n_flipmap[x,:,:] = cv2.resize(flipAngleMaprescale_3D_grad[x,:,:], dsize=(Data_mat[2], Data_mat[1]))
        n_adc[x,:,:] = cv2.resize(ADC_3D_grad[x,:,:], dsize=(Data_mat[2], Data_mat[1]))
    
    for x in range(Data_mat[1]):
        final_t1[:,x,:] = cv2.resize(n_t1[:,x,:], dsize=(Data_mat[2], Data_mat[0]))
        final_t2[:,x,:] = cv2.resize(n_t2[:,x,:], dsize=(Data_mat[2], Data_mat[0]))
        final_m0[:,x,:] = cv2.resize(n_m0[:,x,:], dsize=(Data_mat[2], Data_mat[0]))
        final_B1[:,x,:] = cv2.resize(n_B1[:,x,:], dsize=(Data_mat[2], Data_mat[0]))
        final_t2star[:,x,:] = cv2.resize(n_t2star[:,x,:], dsize=(Data_mat[2], Data_mat[0]))
        final_flipmap[:,x,:] = cv2.resize(n_flipmap[:,x,:], dsize=(Data_mat[2], Data_mat[0]))
        final_adc[:,x,:] = cv2.resize(n_adc[:,x,:], dsize=(Data_mat[2], Data_mat[0]))
 
    if Pre_def_seq == 'SE':
        SE_3D = spin_echo_seq(TR, TE, T1_3D_grad, T2_3D_grad, M0_3D_grad)
        SE_3D = np.multiply(SE_3D, B1map_3D_grad)
        final_seq = resize(SE_3D, T1_3D_grad, Data_mat)
        n = noise_generation(Data_mat, Resolution, Bandwidth, B0_field)
        fig, ax = plt.subplots(1,3)
        fig3D = pltaxes_SE.interactivePlot(fig, ax, final_seq + n, TR, TE, final_t1, final_t2, final_m0, final_B1, n, Data_mat, plotAxis = 2, fov = FOV)
        plt.show()
        
    elif Pre_def_seq == 'GE':         
        angle = flipAngleMaprescale_3D_grad/Alpha
        GE_3D = Gradient_seq(TR, TE, T1_3D_grad, t2_star_3D_grad, M0_3D_grad, angle)
        GE_3D = np.multiply(GE_3D, B1map_3D_grad)
        final_seq = resize(GE_3D, T1_3D_grad, Data_mat)
        n = noise_generation(Data_mat, Resolution, Bandwidth, B0_field)
        fig, ax = plt.subplots(1,3)
        fig3D = pltaxes_GE.interactivePlot(fig, ax, final_seq + n, TR, TE, final_t1, final_t2star, final_m0, final_flipmap, Alpha, final_B1, n, Data_mat, plotAxis = 2, fov = FOV)
        plt.show()
        
    elif Pre_def_seq == 'IN': 
        IN_3D = IN_seq(TR, TI, T1_3D_grad, T2_3D_grad, M0_3D_grad)
        IN_3D = np.multiply(IN_3D, B1map_3D_grad)
        trajectory = post_traj.get()              
        IN_3D = post_TSE_TEeff(IN_3D, post_TSE_TE, post_TSE_ETL, T2_3D_grad, trajectory)            
        final_seq = resize(IN_3D, T1_3D_grad, Data_mat)
        final_seq = TSE_seq_R(post_TSE_TE, post_TSE_ETL, T2_3D_grad, trajectory, Readout_axis, final_seq)
        n = noise_generation(Data_mat, Resolution, Bandwidth, B0_field)      
        fig, ax = plt.subplots(1,3)
        fig3D = pltaxes_IN.interactivePlot(fig, ax, final_seq + n, TR, post_TSE_TE, TI, final_t1, final_t2, final_m0, final_B1, trajectory, Readout_axis, post_TSE_ETL, n, Data_mat, plotAxis = 2, fov = FOV)
        plt.show()
            
    elif Pre_def_seq == 'Double IN':         
        DIN_3D = DoubleInversion_seq(TR, post_TSE_TE, TI, TI2, T1_3D_grad, T2_3D_grad, M0_3D_grad)
        DIN_3D = np.multiply(DIN_3D, B1map_3D_grad)   
        trajectory = post_traj.get()           
        DIN_3D = post_TSE_TEeff(DIN_3D, post_TSE_TE, post_TSE_ETL, T2_3D_grad, trajectory)
        final_seq = resize(DIN_3D, T1_3D_grad, Data_mat)
            
        # This code is because in double inversion recovery their is an exp with a positive value (usualy negative)
        # and the resizing code use a linear interpolation that can give some very small values for T1 --> TR/T1 quite big
        # resulting in masive (overflow) values. So the code bellow puts all values of resized T1 in the correct T1 bins (0 / 0.157 / 0.272 / 0.33 / 4)
        t = final_t1
        t[(0.4 < t)] = 4
        t[t <= 0.157/2] = 0
        t[((0.157/2 < t) & (t <= (0.272-0.157)/2 + 0.157))] = 0.157
        t[(((0.272-0.157)/2 + 0.157 < t) & (t <= (0.33-0.272)/2 + 0.272))] = 0.272
        t[(((0.33-0.272)/2 + 0.272 < t) & (t <= 0.4))] = 0.33
        
        final_seq = TSE_seq_R(post_TSE_TE, post_TSE_ETL, T2_3D_grad, trajectory, Readout_axis, final_seq) 
        n = noise_generation(Data_mat, Resolution, Bandwidth, B0_field)       
        fig, ax = plt.subplots(1,3)
        fig3D = pltaxes_DIN.interactivePlot(fig, ax, final_seq + n, TR, post_TSE_TE, TI, TI2, t, final_t2, final_m0, final_B1, trajectory, Readout_axis, ETL, n, Data_mat, plotAxis = 2, fov = FOV)        
        plt.show()
                
    elif Pre_def_seq == 'FLAIR':       
        FLAIR_3D = IN_seq(TR, TI, T1_3D_grad, T2_3D_grad, M0_3D_grad)
        FLAIR_3D = np.multiply(FLAIR_3D, B1map_3D_grad)
        trajectory = post_traj.get()          
        FLAIR_3D = post_TSE_TEeff(FLAIR_3D, post_TSE_TE, post_TSE_ETL, T2_3D_grad, trajectory)
        final_seq = resize(FLAIR_3D, T1_3D_grad, Data_mat)
        final_seq = TSE_seq_R(post_TSE_TE, post_TSE_ETL, T2_3D_grad, trajectory, Readout_axis, final_seq)
        n = noise_generation(Data_mat, Resolution, Bandwidth, B0_field)        
        fig, ax = plt.subplots(1,3)
        fig3D = pltaxes_FLAIR.interactivePlot(fig, ax, final_seq + n, TR, post_TSE_TE, TI, final_t1, final_t2, final_m0, final_B1, trajectory, Readout_axis, ETL, n, Data_mat, plotAxis = 2, fov = FOV)
        plt.show()
        
    elif Pre_def_seq == 'Dif':
        DIF_3D = Diffusion_seq(TR, T1_3D_grad, T2_3D_grad, M0_3D_grad, B, ADC_3D_grad)
        DIF_3D = np.multiply(DIF_3D, B1map_3D_grad)
        trajectory = post_traj.get()  
        DIF_3D = post_TSE_TEeff(DIF_3D, post_TSE_TE, post_TSE_ETL, T2_3D_grad, trajectory)           
        final_seq = resize(DIF_3D, T1_3D_grad, Data_mat)
        final_seq = TSE_seq_R(post_TSE_TE, post_TSE_ETL, T2_3D_grad, trajectory, Readout_axis, final_seq)
        n = noise_generation(Data_mat, Resolution, Bandwidth, B0_field)        
        fig, ax = plt.subplots(1,3)
        fig3D = pltaxes_Dif.interactivePlot(fig, ax, final_seq + n, TR, TE, final_t1, final_t2, final_m0, B, final_adc, final_B1, trajectory, Readout_axis, ETL, n, Data_mat, plotAxis = 2, fov = FOV)
        plt.show()

    elif Pre_def_seq == 'TSE':
        trajectory = traj.get()
        TSE_3D = 0
        TEeff = 0
        if trajectory == 'Out-in':    
            TEeff = TE * ETL    
            TSE_3D = spin_echo_seq(TR, TEeff, T1_3D_grad, T2_3D_grad, M0_3D_grad) # Computing a spin echo
        elif trajectory == 'In-out':
            TEeff = TE
            TE_div = np.divide(TEeff, T2_3D_grad, out=np.zeros_like(T2_3D_grad), where=T2_3D_grad!=0)
            TSE_3D = np.abs(M0_3D_grad * np.exp(-TE_div))                         # Computing a spin echo with the (1 - exp(-TR/T1)) set to 1
        elif trajectory == 'Linear':
            TEeff = 0.5 * TE * ETL
            TSE_3D = spin_echo_seq(TR, TEeff, T1_3D_grad, T2_3D_grad, M0_3D_grad) # Computing a spin echo
        
        TSE_3D = np.multiply(TSE_3D, B1map_3D_grad)
        final_seq = resize(TSE_3D, T1_3D_grad, Data_mat)            
        final_seq = TSE_seq_R(TE, ETL, T2_3D_grad, trajectory, Readout_axis, final_seq)
        n = noise_generation(Data_mat, Resolution, Bandwidth, B0_field)              
        fig, ax = plt.subplots(1,3)
        fig3D = pltaxes_TSE.interactivePlot(fig, ax, final_seq + n, TR, TE, final_t1, final_t2, ETL, trajectory, Readout_axis, final_m0, final_B1, n, Data_mat, plotAxis = 2, fov = FOV)
        plt.show()
        
    elif Pre_def_seq == 'SSFP': 
        gamma = 42.58*(10**6)*2*constants.pi # 267538030.37970677 [rad/sT]
        omega = B0_3D * gamma
        center = np.divide(omega.shape,2).astype(int)
        center_freq_value = omega[center[0],center[1],center[2]]
        offset = omega - center_freq_value
        phi = offset * TR
        SSFP_3D = SSFP_new(TR, T1_3D_grad, T2_3D_grad, M0_3D_grad, Alpha, phi)
        SSFP_3D = np.multiply(SSFP_3D, B1map_3D_grad)
                
        # Resizing
        for x in range(T1_3D_grad.shape[0]):
            n_seq[x,:,:] = cv2.resize(SSFP_3D[x,:,:], dsize=(Data_mat[2], Data_mat[1]))
            n_phi[x,:,:] = cv2.resize(phi[x,:,:], dsize=(Data_mat[2], Data_mat[1]))
            n_offset[x,:,:] = cv2.resize(offset[x,:,:], dsize=(Data_mat[2], Data_mat[1]))
        for x in range(Data_mat[1]):
            final_seq[:,x,:] = cv2.resize(n_seq[:,x,:], dsize=(Data_mat[2], Data_mat[0]))
            final_phi[:,x,:] = cv2.resize(n_phi[:,x,:], dsize=(Data_mat[2], Data_mat[0]))
            final_offset[:,x,:] = cv2.resize(n_offset[:,x,:], dsize=(Data_mat[2], Data_mat[0]))
            
        n = noise_generation(Data_mat, Resolution, Bandwidth, B0_field) 
        fig, ax = plt.subplots(1,3)
        fig3D = pltaxes_SSFP.interactivePlot(fig, ax, final_seq + n, TR, final_t1, final_t2, final_m0, Alpha, final_phi, final_offset, final_B1, n, Data_mat, plotAxis = 2, fov = FOV)
        plt.show()
        
#/////////// Function opening a filter interactive window //////////////# 
def filter_vis(*args):
    """
    Opens a new window where the user can interact with the simulated filtered data in 3D and change the filter's parameters to see their effect in real time
    Input : None
    Output: None
    """
    global final_seq
    global low_pass_entry
    global high_pass_entry
    global Gauss_entry
    global Non_local_h_entry
    global Non_local_psize_entry
    global Non_local_pdist_entry
    
    if fi.get() == 'low':
        ks = low_pass_entry.get(); ks = int(ks)
        fig, ax = plt.subplots(1,3)
        fig3D = pltaxes_low.interactivePlot(fig, ax, final_seq, ks, Data_mat, plotAxis = 2, fov = FOV)
        plt.show()
        
    elif fi.get() == 'high':
        ks = high_pass_entry.get(); ks = int(ks)
        fig, ax = plt.subplots(1,3)
        fig3D = pltaxes_high.interactivePlot(fig, ax, final_seq, ks, Data_mat, plotAxis = 2, fov = FOV)
        plt.show()
        
    elif fi.get() == 'gauss':
        sig = Gauss_entry.get()
        fig, ax = plt.subplots(1,3)
        fig3D = pltaxes_gauss.interactivePlot(fig, ax, final_seq, sig, Data_mat, plotAxis = 2, fov = FOV)
        plt.show()
        
#/////////////////// Functions for post-processing ///////////////////
def lowpass(s):
    """
    Applys a low pass filter
    Input : s --> kernel size
    Output: None
    """
    global final_im_sag
    global final_im_cor
    global final_im_ax
    global SNR_num_label
    global Data_mat
    global FOV
    global undo
    global final_im_sag_original
    global final_im_cor_original
    global final_im_ax_original
    global canvas1
    global canvas2
    global canvas3
    
    canvas1.get_tk_widget().destroy()
    canvas2.get_tk_widget().destroy()
    canvas3.get_tk_widget().destroy()
    
    # To change the axis value to mm (the fov)
    ax_pos = [[0, Data_mat[0]/2, Data_mat[0]-1], [0, Data_mat[1]/2, Data_mat[1]-1], [0, Data_mat[2]/2, Data_mat[2]-1]]
    ax_label = [[0, FOV[0]/2, FOV[0]], [0, FOV[1]/2, FOV[1]], [0, FOV[2]/2, FOV[2]]]
    Labels = ['Left/Right [mm]', 'Anterior/Posterior [mm]', 'Foot/Head [mm]']
    
    s = int(s)
        
    if undo == 'no': # will compute the filter on the images that could have already being filtered once

        final_im_sag = cv2.blur(final_im_sag,(s,s))
        final_im_cor = cv2.blur(final_im_cor,(s,s))
        final_im_ax = cv2.blur(final_im_ax,(s,s))
        
    elif undo == 'yes': # will compute the filter on the original simulate images 
    
        final_im_sag = cv2.blur(final_im_sag_original,(s,s))
        final_im_cor = cv2.blur(final_im_cor_original,(s,s))
        final_im_ax = cv2.blur(final_im_ax_original,(s,s))
    
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.flip(np.rot90(final_im_sag, k = 1), axis = 1), cmap='gray')  
    plt.xlabel(Labels[1]); plt.xticks(ax_pos[1], ax_label[1])
    plt.ylabel(Labels[2]); plt.yticks(ax_pos[2], ax_label[2])
    canvas1 = FigureCanvasTkAgg(fig, root)                  
    canvas1.draw()
    canvas1.get_tk_widget().grid(row = 1, column = 1, rowspan = 4)      
    plt.close()
    
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.rot90(final_im_cor, k = 1), cmap='gray')  
    plt.xlabel(Labels[0]); plt.xticks(ax_pos[0], ax_label[0])
    plt.ylabel(Labels[2]); plt.yticks(ax_pos[2], ax_label[2])
    canvas2 = FigureCanvasTkAgg(fig, root)                  
    canvas2.draw()
    canvas2.get_tk_widget().grid(row = 1, column = 2, rowspan = 4)      
    plt.close()
    
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.rot90(final_im_ax, k = 1), cmap='gray')  
    plt.xlabel(Labels[0]); plt.xticks(ax_pos[0], ax_label[0])
    plt.ylabel(Labels[1]); plt.yticks(ax_pos[1], ax_label[1])
    canvas3 = FigureCanvasTkAgg(fig, root)                  
    canvas3.draw()
    canvas3.get_tk_widget().grid(row = 1, column = 3, rowspan = 4)      
    plt.close()
    
    m = SNR_length
    # std box
    A = int(SNR_noise_box_center[1])
    B = int(SNR_noise_box_center[0])
    # mean box
    a = int(SNR_mean_box_center[1])
    b = int(SNR_mean_box_center[0])

    snr = snr_homemade(np.rot90(final_im_ax),A-m,A+m,B-m,B+m,a-m,a+m,b-m,b+m)
    SNR_num_label.grid_forget()
    SNR_num_label = Label(frame11, text = str(np.abs(round(snr,2))), font=("Helvetica", 12))
    SNR_num_label.grid(row = 2, column = 1)
        
    # Update the history frame
    history('low')
    
def highpass(s):
    """
    Applys a high pass filter
    Input : s --> kernel size
    Output: None
    """
    global final_im_sag
    global final_im_cor
    global final_im_ax
    global SNR_num_label
    global Data_mat
    global FOV
    global undo
    global final_im_sag_original
    global final_im_cor_original
    global final_im_ax_original
    global canvas1
    global canvas2
    global canvas3
    
    canvas1.get_tk_widget().destroy()
    canvas2.get_tk_widget().destroy()
    canvas3.get_tk_widget().destroy()
    
    ax_pos = [[0, Data_mat[0]/2, Data_mat[0]-1], [0, Data_mat[1]/2, Data_mat[1]-1], [0, Data_mat[2]/2, Data_mat[2]-1]]
    ax_label = [[0, FOV[0]/2, FOV[0]], [0, FOV[1]/2, FOV[1]], [0, FOV[2]/2, FOV[2]]]
    Labels = ['Left/Right [mm]', 'Anterior/Posterior [mm]', 'Foot/Head [mm]']
    
    s = int(s)  
    
    if undo == 'no': # will compute the filter on the images that could have already being filtered once

        final_im_sag = cv2.Laplacian(final_im_sag, -1, ksize = s) # The -1 is so that the final image as the same depth as the original one
        final_im_cor = cv2.Laplacian(final_im_cor, -1, ksize = s)
        final_im_ax = cv2.Laplacian(final_im_ax, -1, ksize = s)
        
    elif undo == 'yes': # will compute the filter on the original simulate images 
    
        final_im_sag = cv2.Laplacian(final_im_sag_original, -1, ksize = s) 
        final_im_cor = cv2.Laplacian(final_im_cor_original, -1, ksize = s)
        final_im_ax = cv2.Laplacian(final_im_ax_original, -1, ksize = s)
    
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.flip(np.rot90(final_im_sag, k = 1), axis = 1), cmap='gray')  
    plt.xlabel(Labels[1]); plt.xticks(ax_pos[1], ax_label[1])
    plt.ylabel(Labels[2]); plt.yticks(ax_pos[2], ax_label[2])
    canvas1 = FigureCanvasTkAgg(fig, root)                  
    canvas1.draw()
    canvas1.get_tk_widget().grid(row = 1, column = 1, rowspan = 4)      
    plt.close()
    
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.rot90(final_im_cor, k = 1), cmap='gray')  
    plt.xlabel(Labels[0]); plt.xticks(ax_pos[0], ax_label[0])
    plt.ylabel(Labels[2]); plt.yticks(ax_pos[2], ax_label[2])
    canvas2 = FigureCanvasTkAgg(fig, root)                  
    canvas2.draw()
    canvas2.get_tk_widget().grid(row = 1, column = 2, rowspan = 4)      
    plt.close()
    
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.rot90(final_im_ax, k = 1), cmap='gray')  
    plt.xlabel(Labels[0]); plt.xticks(ax_pos[0], ax_label[0])
    plt.ylabel(Labels[1]); plt.yticks(ax_pos[1], ax_label[1])
    canvas3 = FigureCanvasTkAgg(fig, root)                  
    canvas3.draw()
    canvas3.get_tk_widget().grid(row = 1, column = 3, rowspan = 4)      
    plt.close()
    
    m = SNR_length
    # std box
    A = int(SNR_noise_box_center[1])
    B = int(SNR_noise_box_center[0])
    # mean box
    a = int(SNR_mean_box_center[1])
    b = int(SNR_mean_box_center[0])

    snr = snr_homemade(np.rot90(final_im_ax),A-m,A+m,B-m,B+m,a-m,a+m,b-m,b+m)
    SNR_num_label.grid_forget()
    SNR_num_label = Label(frame11, text = str(np.abs(round(snr,2))), font=("Helvetica", 12))
    SNR_num_label.grid(row = 2, column = 1)

    # Update the history frame
    history('high')
        
def gauss(sig):
    """
    Applys a Gaussian filter
    Input : s --> standard deviation of Gaussian kernel
    Output: None
    """
    global final_im_sag
    global final_im_cor
    global final_im_ax
    global SNR_num_label
    global Data_mat
    global FOV
    global undo
    global final_im_sag_original
    global final_im_cor_original
    global final_im_ax_original
    global canvas1
    global canvas2
    global canvas3
    
    canvas1.get_tk_widget().destroy()
    canvas2.get_tk_widget().destroy()
    canvas3.get_tk_widget().destroy()
    
    ax_pos = [[0, Data_mat[0]/2, Data_mat[0]-1], [0, Data_mat[1]/2, Data_mat[1]-1], [0, Data_mat[2]/2, Data_mat[2]-1]]
    ax_label = [[0, FOV[0]/2, FOV[0]], [0, FOV[1]/2, FOV[1]], [0, FOV[2]/2, FOV[2]]]
    Labels = ['Left/Right [mm]', 'Anterior/Posterior [mm]', 'Foot/Head [mm]']
    
    if undo == 'no': # will compute the filter on the images that could have already being filtered once

        final_im_sag = gaussian_filter(final_im_sag, sigma=float(sig))
        final_im_cor = gaussian_filter(final_im_cor, sigma=float(sig))
        final_im_ax = gaussian_filter(final_im_ax, sigma=float(sig))
        
    elif undo == 'yes': # will compute the filter on the original simulate images 
    
        final_im_sag = gaussian_filter(final_im_sag_original, sigma=float(sig))
        final_im_cor = gaussian_filter(final_im_cor_original, sigma=float(sig))
        final_im_ax = gaussian_filter(final_im_ax_original, sigma=float(sig))

    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.flip(np.rot90(final_im_sag, k = 1), axis = 1), cmap='gray')  
    plt.xlabel(Labels[1]); plt.xticks(ax_pos[1], ax_label[1])
    plt.ylabel(Labels[2]); plt.yticks(ax_pos[2], ax_label[2])
    canvas1 = FigureCanvasTkAgg(fig, root)                  
    canvas1.draw()
    canvas1.get_tk_widget().grid(row = 1, column = 1, rowspan = 4)      
    plt.close()
    
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.rot90(final_im_cor, k = 1), cmap='gray')  
    plt.xlabel(Labels[0]); plt.xticks(ax_pos[0], ax_label[0])
    plt.ylabel(Labels[2]); plt.yticks(ax_pos[2], ax_label[2])
    canvas2 = FigureCanvasTkAgg(fig, root)                  
    canvas2.draw()
    canvas2.get_tk_widget().grid(row = 1, column = 2, rowspan = 4)      
    plt.close()
    
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.rot90(final_im_ax, k = 1), cmap='gray')  
    plt.xlabel(Labels[0]); plt.xticks(ax_pos[0], ax_label[0])
    plt.ylabel(Labels[1]); plt.yticks(ax_pos[1], ax_label[1])
    canvas3 = FigureCanvasTkAgg(fig, root)                  
    canvas3.draw()
    canvas3.get_tk_widget().grid(row = 1, column = 3, rowspan = 4)      
    plt.close()
    
    m = SNR_length
    # std box
    A = int(SNR_noise_box_center[1])
    B = int(SNR_noise_box_center[0])
    # mean box
    a = int(SNR_mean_box_center[1])
    b = int(SNR_mean_box_center[0])

    snr = snr_homemade(np.rot90(final_im_ax),A-m,A+m,B-m,B+m,a-m,a+m,b-m,b+m)
    SNR_num_label.grid_forget()
    SNR_num_label = Label(frame11, text = str(np.abs(round(snr,2))), font=("Helvetica", 12))
    SNR_num_label.grid(row = 2, column = 1)
        
    # Update the history frame
    history('gauss')
    
def non_local(H, s, d):
    """
    Applys a non-local filter
    Input : H --> Cut-off distance, higher h results in a more permissive filter in accepting patches
            s --> patch size
            d --> patch distance
    Output: None
    """
    global final_im_sag
    global final_im_cor
    global final_im_ax
    global SNR_num_label
    global Data_mat
    global FOV
    global undo
    global final_im_sag_original
    global final_im_cor_original
    global final_im_ax_original
    global canvas1
    global canvas2
    global canvas3
    
    canvas1.get_tk_widget().destroy()
    canvas2.get_tk_widget().destroy()
    canvas3.get_tk_widget().destroy()
        
    # To change the axis value to mm (the fov)
    ax_pos = [[0, Data_mat[0]/2, Data_mat[0]-1], [0, Data_mat[1]/2, Data_mat[1]-1], [0, Data_mat[2]/2, Data_mat[2]-1]]
    ax_label = [[0, FOV[0]/2, FOV[0]], [0, FOV[1]/2, FOV[1]], [0, FOV[2]/2, FOV[2]]]
    Labels = ['Left/Right [mm]', 'Anterior/Posterior [mm]', 'Foot/Head [mm]']
    
    H = float(H); s = int(s); d = int(d)

    if undo == 'no': # will compute the filter on the images that could have already being filtered once
    
        # estimate the noise standard deviation from the noisy image
        sigma_est_sag = estimate_sigma(final_im_sag, channel_axis=None)
        sigma_est_cor = estimate_sigma(final_im_cor, channel_axis=None)
        sigma_est_ax = estimate_sigma(final_im_ax, channel_axis=None)

        final_im_sag = denoise_nl_means(final_im_sag, h=H * sigma_est_sag, sigma=sigma_est_sag, fast_mode=False, patch_size=s, patch_distance=d)
        final_im_cor = denoise_nl_means(final_im_cor, h=H * sigma_est_cor, sigma=sigma_est_sag, fast_mode=False, patch_size=s, patch_distance=d)
        final_im_ax = denoise_nl_means(final_im_ax, h=H * sigma_est_ax, sigma=sigma_est_sag, fast_mode=False, patch_size=s, patch_distance=d)
        
    elif undo == 'yes': # will compute the filter on the original simulate images 
    
        # estimate the noise standard deviation from the noisy image
        sigma_est_sag = estimate_sigma(final_im_sag_original, channel_axis=None)
        sigma_est_cor = estimate_sigma(final_im_cor_original, channel_axis=None)
        sigma_est_ax = estimate_sigma(final_im_ax_original, channel_axis=None)

        final_im_sag = denoise_nl_means(final_im_sag_original, h=H * sigma_est_sag, sigma=sigma_est_sag, fast_mode=False, patch_size=s, patch_distance=d)
        final_im_cor = denoise_nl_means(final_im_cor_original, h=H * sigma_est_cor, sigma=sigma_est_sag, fast_mode=False, patch_size=s, patch_distance=d)
        final_im_ax = denoise_nl_means(final_im_ax_original, h=H * sigma_est_ax, sigma=sigma_est_sag, fast_mode=False, patch_size=s, patch_distance=d)
        
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.flip(np.rot90(final_im_sag, k = 1), axis = 1), cmap='gray')  
    plt.xlabel(Labels[1]); plt.xticks(ax_pos[1], ax_label[1])
    plt.ylabel(Labels[2]); plt.yticks(ax_pos[2], ax_label[2])
    canvas1 = FigureCanvasTkAgg(fig, root)                  
    canvas1.draw()
    canvas1.get_tk_widget().grid(row = 1, column = 1, rowspan = 4)      
    plt.close()
    
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.rot90(final_im_cor, k = 1), cmap='gray')  
    plt.xlabel(Labels[0]); plt.xticks(ax_pos[0], ax_label[0])
    plt.ylabel(Labels[2]); plt.yticks(ax_pos[2], ax_label[2])
    canvas2 = FigureCanvasTkAgg(fig, root)                  
    canvas2.draw()
    canvas2.get_tk_widget().grid(row = 1, column = 2, rowspan = 4)      
    plt.close()
    
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.rot90(final_im_ax, k = 1), cmap='gray')  
    plt.xlabel(Labels[0]); plt.xticks(ax_pos[0], ax_label[0])
    plt.ylabel(Labels[1]); plt.yticks(ax_pos[1], ax_label[1])
    canvas3 = FigureCanvasTkAgg(fig, root)                  
    canvas3.draw()
    canvas3.get_tk_widget().grid(row = 1, column = 3, rowspan = 4)      
    plt.close()
    
    m = SNR_length
    # std box
    A = int(SNR_noise_box_center[1])
    B = int(SNR_noise_box_center[0])
    # mean box
    a = int(SNR_mean_box_center[1])
    b = int(SNR_mean_box_center[0])

    snr = snr_homemade(np.rot90(final_im_ax),A-m,A+m,B-m,B+m,a-m,a+m,b-m,b+m)
    SNR_num_label.grid_forget()
    SNR_num_label = Label(frame11, text = str(np.abs(round(snr,2))), font=("Helvetica", 12))
    SNR_num_label.grid(row = 2, column = 1)
        
    # Update the history frame
    history('nonloc')
    
def filter_undo():
    """
    Undo any previous filter application
    Input : None
    Output: None
    """
    global final_im_sag_original 
    global final_im_cor_original 
    global final_im_ax_original 
    global Labels
    global SNR_num_label
    global ax_pos
    global ax_label
    global undo
    global canvas1
    global canvas2
    global canvas3
    
    canvas1.get_tk_widget().destroy()
    canvas2.get_tk_widget().destroy()
    canvas3.get_tk_widget().destroy()
    
    # Create axial plot
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.rot90(final_im_ax_original, k = 1), cmap='gray')
    plt.xlabel(Labels[0]); plt.xticks(ax_pos[0], ax_label[0])
    plt.ylabel(Labels[1]); plt.yticks(ax_pos[1], ax_label[1])
    canvas1 = FigureCanvasTkAgg(fig, root)                
    canvas1.draw()
    canvas1.get_tk_widget().grid(row = 1, column = 3, rowspan = 4)     
    plt.close()

    # Create coronal plot
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.rot90(final_im_cor_original, k = 1), cmap='gray')   
    plt.xlabel(Labels[0]); plt.xticks(ax_pos[0], ax_label[0])
    plt.ylabel(Labels[2]); plt.yticks(ax_pos[2], ax_label[2])
    canvas2 = FigureCanvasTkAgg(fig, root)         
    canvas2.draw()
    canvas2.get_tk_widget().grid(row = 1, column = 2, rowspan = 4)
    plt.close()

    # Create sagittal plot
    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.flip(np.rot90(final_im_sag_original, k = 1), axis = 1), cmap='gray')   
    plt.xlabel(Labels[1]); plt.xticks(ax_pos[1], ax_label[1])
    plt.ylabel(Labels[2]); plt.yticks(ax_pos[2], ax_label[2])
    canvas3 = FigureCanvasTkAgg(fig, root)         
    canvas3.draw()
    canvas3.get_tk_widget().grid(row = 1, column = 1, rowspan = 4) 
    plt.close()
    
    m = SNR_length
    # std box
    A = int(SNR_noise_box_center[1])
    B = int(SNR_noise_box_center[0])
    # mean box
    a = int(SNR_mean_box_center[1])
    b = int(SNR_mean_box_center[0])

    snr = snr_homemade(np.rot90(final_im_ax_original),A-m,A+m,B-m,B+m,a-m,a+m,b-m,b+m)
    SNR_num_label.grid_forget()
    SNR_num_label = Label(frame11, text = str(np.abs(round(snr,2))), font=("Helvetica", 12))
    SNR_num_label.grid(row = 2, column = 1)
    
    history('undo')
    undo = 'yes'
    
#/////////////////// Functions regarding the change in field strength ///////////////////
def pr_b0_46(*args):
    """
    Changes the field strength variable and label to 46mT 
    Input : None
    Output: None
    """
    global B0_field
    global field_strength_label_value
    B0_field = 'B0_46mT'
    field_strength_label_value.grid_forget()
    field_strength_label_value = Label(root, text = "46 mT", font=("Helvetica", 17));field_strength_label_value.grid(row = 6, column = 0)
    
def pr_b0_15(*args):
    """
    Changes the field strength variable and label to 1.5T 
    Input : None
    Output: None
    """
    global B0_field
    global field_strength_label_value
    B0_field = 'B0_15T'
    field_strength_label_value.grid_forget()
    field_strength_label_value = Label(root, text = "1.5 T", font=("Helvetica", 17));field_strength_label_value.grid(row = 6, column = 0)
    
def pr_b0_3(*args):
    """
    Changes the field strength variable and label to 3T 
    Input : None
    Output: None
    """
    global B0_field
    global field_strength_label_value
    B0_field = 'B0_3T'
    field_strength_label_value.grid_forget()
    field_strength_label_value = Label(root, text = "3 T", font=("Helvetica", 17));field_strength_label_value.grid(row = 6, column = 0)
    
def settings():
    """
    Function showing a new window with the different fields strength the user can choose from and which maps are affected (B0, T1, and T2) 
    Input: None
    Output: 
    """
    global B0_3D
    global top
    global B0_selected_label
    global FOV
    global Resolution

    T1_3D_46 = np.load(resource_path("Data\T1_3D_cer_lip_grad.npy"))
    T2_3D_46 = np.load(resource_path("Data\T2_3D_cer_lip_grad.npy"))
    T1_3D_15 = np.load(resource_path("Data\T1_15.npy"))
    T2_3D_15 = np.load(resource_path("Data\T2_15.npy")) 
    T1_3D_3 = np.load(resource_path("Data\T1_3.npy"))
    T2_3D_3 = np.load(resource_path("Data\T2_3.npy"))

    FOV = [int(FOV1_entry.get()), int(FOV2_entry.get()), int(FOV3_entry.get())]
    Resolution = [float(Res1_entry.get()), float(Res2_entry.get()), float(Res3_entry.get())]
    Data_mat = np.divide(FOV, Resolution)
    Data_mat = [int(Data_mat[0]), int(Data_mat[1]), int(Data_mat[2])]   
    
    # To change the axis value to mm (the fov)
    ax_pos = [[0, FOV[0]/2, FOV[0]], [0, FOV[1]/2, FOV[1]], [0, FOV[2]/2, FOV[2]]]
    ax_label = ax_pos
    Labels = ['Left/Right [mm]', 'Anterior/Posterior [mm]', 'Foot/Head [mm]']
    
    top = Toplevel()
    top.title("Settings")
    
    # Labels
    B0_label = Label(top, text="B0 fields", font=("Helvetica", 16)).grid(row = 1, column = 0)
    B01_label = Label(top, text="46 mT", font=("Helvetica", 16)).grid(row = 0, column = 1)
    B02_label = Label(top, text="1.5 T", font=("Helvetica", 16)).grid(row = 0, column = 2)
    B03_label = Label(top, text="3 T", font=("Helvetica", 16)).grid(row = 0, column = 3)
    T1_label = Label(top, text="T1 values", font=("Helvetica", 16)).grid(row = 2, column = 0)
    T2_label = Label(top, text="T2", font=("Helvetica", 16)).grid(row = 3, column = 0)
    
    s = B0_3D.shape
    B0_46mT = B0_3D[:,:,134]
    B0_15T = B0_46mT * 32.608 # 1.5T / 46mT ~= 32.608
    B0_3T = B0_46mT * 65.217  # 3T / 46mT ~= 65.217
    
    B0_46mT[B0_46mT==0] = np.nan
    B0_15T[B0_15T==0] = np.nan
    B0_3T[B0_3T==0] = np.nan
    
    B0_15T[B0_15T>0] = 1.5
    B0_3T[B0_3T>0] = 3
    
    proton_grad = np.copy(M0_3D_grad[:,:,135])
    proton_grad[B0_46mT == np.nan] = np.nan
    
    proton = np.copy(M0_3D[:,:,135])
    proton[B0_46mT == np.nan] = np.nan

    imag_size = 3
       
    fig = plt.figure(figsize=(imag_size,imag_size))                
    plt.imshow(np.rot90(B0_46mT + proton_grad * 0.00003229247027748)) # 0.00003229247027748 is just a scaling factor
    plt.colorbar()
    canvas = FigureCanvasTkAgg(fig, top)                          
    canvas.draw()
    canvas.get_tk_widget().grid(row = 1, column = 1)     
    canvas.get_tk_widget().bind('<Button-1>', pr_b0_46)
    plt.close()
    
    fig = plt.figure(figsize=(imag_size,imag_size))                
    plt.imshow(np.rot90(B0_15T + (proton * 0.001)))   
    plt.colorbar()
    canvas = FigureCanvasTkAgg(fig, top)                          
    canvas.draw()
    canvas.get_tk_widget().grid(row = 1, column = 2) 
    canvas.get_tk_widget().bind('<Button-1>', pr_b0_15)
    plt.close()
    
    fig = plt.figure(figsize=(imag_size,imag_size))                
    plt.imshow(np.rot90(B0_3T + (proton * 0.001)))          
    plt.colorbar()
    canvas = FigureCanvasTkAgg(fig, top)                          
    canvas.draw()
    canvas.get_tk_widget().grid(row = 1, column = 3)     
    canvas.get_tk_widget().bind('<Button-1>', pr_b0_3)
    plt.close()
    
    ### T1 maps 
    fig = plt.figure(figsize=(imag_size,imag_size))                
    plt.imshow(np.rot90(T1_3D_46[:,:,134]))          
    plt.colorbar()
    canvas = FigureCanvasTkAgg(fig, top)                          
    canvas.draw()
    canvas.get_tk_widget().grid(row = 2, column = 1)     
    canvas.get_tk_widget().bind('<Button-1>', pr_b0_46)
    plt.close()
    
    fig = plt.figure(figsize=(imag_size,imag_size))                
    plt.imshow(np.rot90(np.divide(T1_3D_15[:,:,134],1000)))
    plt.colorbar()
    canvas = FigureCanvasTkAgg(fig, top)                          
    canvas.draw()
    canvas.get_tk_widget().grid(row = 2, column = 2)     
    canvas.get_tk_widget().bind('<Button-1>', pr_b0_15)
    plt.close()
    
    fig = plt.figure(figsize=(imag_size,imag_size))                
    plt.imshow(np.rot90(np.divide(T1_3D_3[:,:,134],1000)))
    plt.colorbar()
    canvas = FigureCanvasTkAgg(fig, top)                          
    canvas.draw()
    canvas.get_tk_widget().grid(row = 2, column = 3)     
    canvas.get_tk_widget().bind('<Button-1>', pr_b0_3)
    plt.close()
    
    ### T2 maps 
    fig = plt.figure(figsize=(imag_size,imag_size))                
    plt.imshow(np.rot90(T2_3D_46[:,:,134]))    
    plt.colorbar()
    canvas = FigureCanvasTkAgg(fig, top)                          
    canvas.draw()
    canvas.get_tk_widget().grid(row = 3, column = 1)     
    canvas.get_tk_widget().bind('<Button-1>', pr_b0_46)
    plt.close()
    
    fig = plt.figure(figsize=(imag_size,imag_size))                
    plt.imshow(np.rot90(np.divide(T2_3D_15[:,:,134],1000))) 
    plt.colorbar()
    canvas = FigureCanvasTkAgg(fig, top)                          
    canvas.draw()
    canvas.get_tk_widget().grid(row = 3, column = 2)     
    canvas.get_tk_widget().bind('<Button-1>', pr_b0_15)
    plt.close()
    
    fig = plt.figure(figsize=(imag_size,imag_size))                
    plt.imshow(np.rot90(np.divide(T2_3D_3[:,:,134],1000))) 
    plt.colorbar()
    canvas = FigureCanvasTkAgg(fig, top)                          
    canvas.draw()
    canvas.get_tk_widget().grid(row = 3, column = 3)     
    canvas.get_tk_widget().bind('<Button-1>', pr_b0_3)
    plt.close()
        
    # To make the buttons/labels adjust to window size
    top.columnconfigure(0, weight=1)
    top.columnconfigure(1, weight=1)
    top.columnconfigure(2, weight=1)
    top.columnconfigure(3, weight=1)
    top.rowconfigure(0, weight=1)
    top.rowconfigure(1, weight=1)
    top.rowconfigure(2, weight=1)
    top.rowconfigure(3, weight=1)
        
#/////////////////// Functions updating the readout variable for TSE and post TSE sequences  ///////////////////
def readout_selection(*args):
    """
    Assign the selected readout axis to the TSE axis readout variable
    Input : None
    Output: None
    """
    global Readout_axis
    R = Readout.get()

    if R == 'Along Sagittal Foot/head':
        Readout_axis = 'FH'
    elif R == 'Along Sagittal Anterior/posterior':
        Readout_axis = 'AP'
    elif R == 'Along Coronal Foot/head':
        Readout_axis = 'FH'
    elif R == 'Along Coronal Left/right':
        Readout_axis = 'LR'
    elif R == 'Along Axial Anterior/posterior':
        Readout_axis = 'AP'
    elif R == 'Along Axial Left/right':
        Readout_axis = 'LR'
        
def post_readout_selection(*args):
    """
    Assign the selected readout axis to the post TSE axis readout variable
    Input : None
    Output: None
    """
    global Readout_axis
    R = Post_Readout.get()

    if R == 'Along Sagittal Foot/head':
        Readout_axis = 'FH'
    elif R == 'Along Sagittal Anterior/posterior':
        Readout_axis = 'AP'
    elif R == 'Along Coronal Foot/head':
        Readout_axis = 'FH'
    elif R == 'Along Coronal Left/right':
        Readout_axis = 'LR'
    elif R == 'Along Axial Anterior/posterior':
        Readout_axis = 'AP'
    elif R == 'Along Axial Left/right':
        Readout_axis = 'LR'

#/////////////////// Functions displaying the different frames ///////////////////
def display_param_frame():
    """
    Function displaying the parameters frame
    Input: None
    Output: None
    """
    global param_frame_on
    global post_process_frame_on
    global param_visu_frame_on
    
    param_frame_on = 'yes'
    
    frame3.grid(row = 7, column = 1, columnspan = 2, rowspan = 2)
    Display_param_frame_button.grid_forget()
    
    if post_process_frame_on == 'no':
        Display_postpro_frame_button.grid(row = 7, column = 3)
    elif post_process_frame_on == 'yes':
        frame4.grid(row = 7, column = 3, rowspan = 2)
        
    if param_visu_frame_on == 'no':
        Display_parvisu_frame_button.grid(row = 9, column = 3)
    elif param_visu_frame_on == 'yes':
        frame12.grid(row = 9, column = 3, rowspan = 2)

def display_postpro_frame():
    """
    Function displaying the post-processing options frame
    Input: None
    Output: None
    """
    global param_frame_on
    global post_process_frame_on
    global param_visu_frame_on
    
    post_process_frame_on = 'yes'
    Display_postpro_frame_button.grid_forget()
    
    if param_frame_on == 'yes':
        frame4.grid(row = 7, column = 3, rowspan = 2)
    elif param_frame_on == 'no':
        frame4.grid(row = 7, column = 2, rowspan = 2)

def display_visu_frame():
    """
    Function displaying the parameters visualisation option frame
    Input: None
    Output: None
    """
    global param_frame_on
    global post_process_frame_on
    global param_visu_frame_on
    
    param_visu_frame_on = 'yes'
    Display_parvisu_frame_button.grid_forget()
    
    if param_frame_on == 'yes':
        frame12.grid(row = 9, column = 3, rowspan = 2)
    elif param_frame_on == 'no':
        frame12.grid(row = 9, column = 2, rowspan = 2)
        
### Function closing the parameters and post processing frames on right click ###
def close_parvisu_display(*args):
    """
    Function closing the parameters frame
    Input: None
    Output: None
    """
    global param_frame_on
    global post_process_frame_on
    global param_visu_frame_on
    
    param_visu_frame_on = 'no'
    frame12.grid_forget() 
    
    if param_frame_on == 'yes':
        Display_parvisu_frame_button.grid(row = 9, column = 3)
    elif param_frame_on == 'no':
        Display_parvisu_frame_button.grid(row = 9, column = 2)   
        
def close_postpro_display(*args):
    """
    Function closing the post-processing options frame
    Input: None
    Output: None
    """
    global param_frame_on
    global post_process_frame_on
    global param_visu_frame_on
    
    post_process_frame_on = 'no'
    frame4.grid_forget() 
    
    if param_frame_on == 'yes':
        Display_postpro_frame_button.grid(row = 7, column = 3)
    elif param_frame_on == 'no':
        Display_postpro_frame_button.grid(row = 7, column = 2)   
        
def close_param_display(*args):
    """
    Function closing the parameters visualisation option frame
    Input: None
    Output: None
    """
    global param_frame_on
    global post_process_frame_on
    global param_visu_frame_on
    
    param_frame_on = 'no'
    frame3.grid_forget()
    Display_param_frame_button.grid(row = 7, column = 1)
    
    if post_process_frame_on == 'no':
        Display_postpro_frame_button.grid(row = 7, column = 2)
    elif post_process_frame_on == 'yes':
        frame4.grid(row = 7, column = 2, rowspan = 2)
        
    if param_visu_frame_on == 'no':
        Display_parvisu_frame_button.grid(row = 9, column = 2)
    elif param_visu_frame_on == 'yes':
        frame12.grid(row = 9, column = 2, rowspan = 2)
        
#/////////// Function regarding the saving of 2D images from the simulator //////////////#  
def save_img(axe):
    """
    Saves a 2D numpy array as a .png grayscale image
    Input: None
    Output: None
    """
    global final_im_sag
    global final_im_cor
    global final_im_ax
    global Pre_def_seq
    
    er_label.grid_forget(); 
    all_axes = ['sagittal', 'coronal', 'axial']
    
    if Pre_def_seq == " ":
        er_label.grid(row = 2, column = 0)
    else:
        if axe in all_axes:       
            filename = filedialog.asksaveasfilename()
            if not filename:
                return
            if axe == 'sagittal':
                matplotlib.image.imsave((filename + str('.png')), np.flip(np.rot90(final_im_sag),axis = 1),cmap='gray')
            elif axe == 'coronal':
                matplotlib.image.imsave((filename + str('.png')), np.rot90(final_im_cor),cmap='gray')
            elif axe == 'axial':
                matplotlib.image.imsave((filename + str('.png')), np.rot90(final_im_ax),cmap='gray')
        else:
            er_label.grid(row = 2, column = 0)
        
#//// Functions opening webpages giving the user information for each sequences or label ///
def info_SE(*args):
    webbrowser.open("https://mriquestions.com/spin-echo1.html")
def info_GE(*args):
    webbrowser.open("https://mriquestions.com/gradient-echo.html")
def info_IN(*args):
    webbrowser.open("https://mriquestions.com/what-is-ir.html")
def info_double_IN(*args):
    webbrowser.open("https://www.mriquestions.com/double-ir.html")
def info_FLAIR(*args):
    webbrowser.open("https://en.wikipedia.org/wiki/Fluid-attenuated_inversion_recovery")
def info_TSE(*args):
    webbrowser.open("https://mriquestions.com/what-is-fsetse.html")
def info_Dif(*args):
    webbrowser.open("https://en.wikipedia.org/wiki/Diffusion_MRI")
def info_SSFP(*args):
    webbrowser.open("https://mriquestions.com/what-is-ssfp.html")
def info_tr_te(*args):
    webbrowser.open("https://mriquestions.com/tr-and-te.html")
def info_ti(*args):
    webbrowser.open("https://mriquestions.com/ti-to-null-a-tissue.html")
def info_ti2(*args):
    webbrowser.open("https://www.mriquestions.com/double-ir.html")
def info_alpha(*args):
    webbrowser.open("https://mriquestions.com/what-is-flip-angle.html")
def info_etl(*args):
    webbrowser.open("https://mriquestions.com/fse-parameters.html")
def info_Kspace_traj(*args):
    webbrowser.open("https://mriquestions.com/k-space-trajectories.html")
def info_teeff(*args):
    webbrowser.open("https://mriquestions.com/fse-parameters.html")
def info_G(*args):
    webbrowser.open("https://mriquestions.com/what-is-the-b-value.html")
def info_smalldelta(*args):
    webbrowser.open("https://mriquestions.com/what-is-the-b-value.html")
def info_bigdelta(*args):
    webbrowser.open("https://mriquestions.com/what-is-the-b-value.html")
def info_B(*args):
    webbrowser.open("https://mriquestions.com/tr-and-te.html")
def info_bd(*args):   
    webbrowser.open('https://mriquestions.com/receiver-bandwidth.html')
def info_fov(*args):   
    webbrowser.open('https://mriquestions.com/field-of-view-fov.html')
def info_res(*args):   
    webbrowser.open('https://mrimaster.com/index.4.html')
def info_snr(*args):   
    webbrowser.open('https://mriquestions.com/signal-to-noise.html')
def info_lowpass(*args):
    webbrowser.open("https://en.wikipedia.org/wiki/Low-pass_filter#:~:text=A%20low%2Dpass%20filter%20is,depends%20on%20the%20filter%20design.")
def info_highpass(*args):
    webbrowser.open("https://en.wikipedia.org/wiki/High-pass_filter")
def info_gausspass(*args):
    webbrowser.open("https://en.wikipedia.org/wiki/Gaussian_filter#:~:text=In%20electronics%20and%20signal%20processing,would%20have%20infinite%20impulse%20response).")    
def info_nonlocal(*args):
    webbrowser.open("https://en.wikipedia.org/wiki/Non-local_means")
    
#//// Functions opening and closing the simulator ///
def open_simu():
    """
    Opens the simulator, display all labels/entries/buttons/frames important for the simulator
    If maps are imported, this function will compute related information (ex: a B0 map is imported --> a new T2* must be computed)
    Input: None
    Output: None
    """
    global canvas1
    global canvas2
    global canvas3
    global B0_import
    global T2_star_import
    global T2_import
    global SIM
    
    SIM = 1
    T2_check = isinstance(T2_import,np.ndarray)
    B0_check = isinstance(B0_import,np.ndarray)
    
    # Check if maps have been imported (If n_import is a list/array --> return True, if n_import is an integer --> return False)
    if (B0_check == True) and (T2_check == False):
        Delta_B0_import = 0
        s = B0_import.shape

        # Compute delta B0 (for T2*) and rescaling
        # If dimension sum of imported map > 300, first downrescale, compute deltaB0, and then rescale (for speed)
        if s[0]+s[1]+s[2] > 300:
            bb = zoom(B0_import, (100/s[0], 100/s[1], 100/s[2]), order=0)
            s = bb.shape
            Delta_B0_import = Trapeze_integral_3D(bb)                     
            Delta_B0_import = zoom(Delta_B0_import, (250/s[0], 300/s[1], 275/s[2]), order=0)

        else: # If dimension sum of imported map < 300 , first compute deltaB0 and then rescale
            Delta_B0_import = Trapeze_integral_3D(B0_import)
            Delta_B0_import = zoom(Delta_B0_import, (250/s[0], 300/s[1], 275/s[2]), order=0)

        T2 = np.load(resource_path("Data\T2_3D_cer_lip_grad.npy")) # unit of T2 map is in [s]
        t2_inverse = np.divide(1, T2, out=np.zeros_like(T2), where=T2!=0) # 1/t2
        gamma =  42.58*(10**6)   # gyromagnetic ratio for hydrogen 42.58 [MHz/T]
        D_gam = Delta_B0_import * gamma
        inv_t2star = t2_inverse + D_gam
        T2_star_import = np.divide(1, inv_t2star, out=np.zeros_like(inv_t2star), where=t2_inverse!=0)

        t2_star_3D_grad = np.load(resource_path("Data\T2_star_tensor_3D_cer_lip_grad.npy"))

        plt.subplot(131)
        plt.imshow(t2_star_3D_grad[:,:,134])
        plt.colorbar()
        plt.subplot(132)
        plt.imshow(T2_star_import[:,:,134])
        plt.colorbar()
        plt.subplot(133)
        plt.imshow(Delta_B0_import[:,:,134])
        plt.colorbar()
        plt.show()
        
        s = B0_import.shape
        B0_import_big = zoom(B0_import, (250/s[0], 300/s[1], 275/s[2]), order=0)
        
        plt.imshow(np.rot90(B0_import[:,:,s[2]//2]))
        plt.colorbar()
        plt.show()
        
    elif (B0_check == False) and (T2_check == True):

        D = np.load(resource_path("Data\DELTA_B0_3D.npy"))
        
        T2_ori = np.load(resource_path("Data\T2_3D_cer_lip_grad.npy")) # unit of T2 map is in [s]
        t2_inverse = np.divide(1, T2_import, out=np.zeros_like(T2_import), where=T2_import!=0) # 1/t2
        gamma =  42.58*(10**6)   # gyromagnetic ratio for hydrogen 42.58 [MHz/T]
        D_gam = D * gamma
        inv_t2star = t2_inverse + D_gam
        T2_star_import = np.divide(1, inv_t2star, out=np.zeros_like(inv_t2star), where=t2_inverse!=0)

        t2_star_3D_grad = np.load(resource_path("Data\T2_star_tensor_3D_cer_lip_grad.npy"))

        plt.subplot(121)
        plt.imshow(T2_ori[:,:,134])
        plt.colorbar()
        plt.subplot(122)
        plt.imshow(T2_import[:,:,134])
        plt.colorbar()
        plt.show()
        
        plt.subplot(131)
        plt.imshow(t2_star_3D_grad[:,:,134])
        plt.colorbar()
        plt.subplot(132)
        plt.imshow(T2_star_import[:,:,134])
        plt.colorbar()
        plt.subplot(133)
        plt.imshow(D[:,:,134])
        plt.colorbar()
        plt.show()   
        
    elif (B0_check == True) and (T2_check == True):
        Delta_B0_import = 0
        s = B0_import.shape

        # Compute delta B0 (for T2*) and rescaling
        # If dimension sum of imported map > 300, first downrescale, compute deltaB0, and then rescale (for speed)
        if s[0]+s[1]+s[2] > 300:
            bb = zoom(B0_import, (100/s[0], 100/s[1], 100/s[2]), order=0)
            s = bb.shape
            Delta_B0_import = Trapeze_integral_3D(bb)                     
            Delta_B0_import = zoom(Delta_B0_import, (250/s[0], 300/s[1], 275/s[2]), order=0)

        else: # If dimension sum of imported map < 300 , first compute deltaB0 and then rescale
            Delta_B0_import = Trapeze_integral_3D(B0_import)
            Delta_B0_import = zoom(Delta_B0_import, (250/s[0], 300/s[1], 275/s[2]), order=0)

        T2_ori = np.load(resource_path("Data\T2_3D_cer_lip_grad.npy")) # unit of T2 map is in [s]
        t2_inverse = np.divide(1, T2_import, out=np.zeros_like(T2_import), where=T2_import!=0) # 1/t2
        gamma =  42.58*(10**6)   # gyromagnetic ratio for hydrogen 42.58 [MHz/T]
        D_gam = Delta_B0_import * gamma
        inv_t2star = t2_inverse + D_gam
        T2_star_import = np.divide(1, inv_t2star, out=np.zeros_like(inv_t2star), where=t2_inverse!=0)

        t2_star_3D_grad = np.load(resource_path("Data\T2_star_tensor_3D_cer_lip_grad.npy"))

        plt.subplot(121)
        plt.imshow(T2_ori[:,:,134])
        plt.colorbar()
        plt.subplot(122)
        plt.imshow(T2_import[:,:,134])
        plt.colorbar()
        plt.show()
        plt.subplot(131)
        plt.imshow(t2_star_3D_grad[:,:,134])
        plt.colorbar()
        plt.subplot(132)
        plt.imshow(T2_star_import[:,:,134])
        plt.colorbar()
        plt.subplot(133)
        plt.imshow(Delta_B0_import[:,:,134])
        plt.colorbar()
        plt.show()
        
    canvas1.get_tk_widget().destroy()
    canvas2.get_tk_widget().destroy()
    
    # Opening the simulator tutorial
    tutorial_simu()
    
    frame1.grid(row = 0, column = 0, rowspan = 2)
    frame3.grid(row = 7, column = 1, columnspan = 2, rowspan = 2)
    Display_postpro_frame_button.grid(row = 7, column = 3)
    Display_parvisu_frame_button.grid(row = 9, column = 3)
    frame10.grid(row = 7, column = 0)
    frame11.grid(row = 3, column = 0)
    setting_button.grid(row = 4, column = 0)
    field_strength_label.grid(row = 5, column = 0)
    field_strength_label_value.grid(row = 6, column = 0)
    Sagittal_label.grid(row = 0, column = 1)
    Coronal_label.grid(row = 0, column = 2)
    Axial_label.grid(row = 0, column = 3)
    reset_button.grid(row = 2, column = 0)
    Exit_simu.grid(row = 9, column = 0) 
    
    imag_size = 5.5
    fig = plt.figure(figsize=(imag_size,imag_size))                
    plt.imshow(np.flip(np.rot90(M0_3D[128,:,:]), axis = 1), cmap='gray')          
    plt.axis('off')
    canvas1 = FigureCanvasTkAgg(fig, root)                          
    canvas1.draw()
    canvas1.get_tk_widget().grid(row = 1, column = 1, rowspan = 4)        
    plt.close()

    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.rot90(M0_3D[:,150,:]), cmap='gray') 
    plt.axis('off')
    canvas2 = FigureCanvasTkAgg(fig, root)         
    canvas2.draw()
    canvas2.get_tk_widget().grid(row = 1, column = 2, rowspan = 4)        
    plt.close()

    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(np.rot90(M0_3D[:,:,135]), cmap='gray')
    plt.axis('off')
    canvas3 = FigureCanvasTkAgg(fig, root)         
    canvas3.draw()
    canvas3.get_tk_widget().grid(row = 1, column = 3, rowspan = 4)          
    plt.close()
    
    open_simulator.grid_forget()
    open_generator.grid_forget()
    continue_to_simu.grid_forget()
    frame13.grid_forget()
    frame14.grid_forget()

def exit_simu():
    """
    Closes the simulator and all labels/entries/buttons/frames related to it
    Input: None
    Output: None
    """
    global history_label1, history_label2, history_label3, history_label4, history_label5
    global canvas1
    global canvas2
    global canvas3
    global T1_import
    global T2_import
    global T2_star_import
    global B0_import
    global Delta_B0_import
    global seq_label2
    global SIM
    global SNR_num_label
    global SNR_length
    global SNR_mean_box_center
    global SNR_noise_box_center
    global Time_scan_num_label
    global Bd_by_pixel_label
    global SNR_num_label
    global minimumTE_num_label
    global minimumTR_num_label
    global warning_label
    
    TR_entry.delete(0, END); TR_entry.insert(0, '500')
    TE_entry.delete(0, END); TE_entry.insert(0, '20')
    TI_entry.delete(0, END); TI_entry.insert(0, '250')
    FOV1_entry.delete(0, END); FOV1_entry.insert(0, '250')
    FOV2_entry.delete(0, END); FOV2_entry.insert(0, '300')
    FOV3_entry.delete(0, END); FOV3_entry.insert(0, '275')
    Res1_entry.delete(0, END); Res1_entry.insert(0, '1.95')
    Res2_entry.delete(0, END); Res2_entry.insert(0, '2.34')
    Res3_entry.delete(0, END); Res3_entry.insert(0, '2.15')
    Data_mat1_entry.delete(0,END); Data_mat1_entry.insert(0,'128')
    Data_mat2_entry.delete(0,END); Data_mat2_entry.insert(0,'128')
    Data_mat3_entry.delete(0,END); Data_mat3_entry.insert(0,'128')
    Bandwidth_entry.delete(0, END); Bandwidth_entry.insert(0, '40000')
    Alpha_entry.delete(0, END); Alpha_entry.insert(0, '45')
    G_entry.delete(0,END); G_entry.insert(0, '10')
    smalldelta_entry.delete(0,END); smalldelta_entry.insert(0, '1')
    bigdelta_entry.delete(0,END);bigdelta_entry.insert(0, '2')

    Time_scan_num_label.grid_forget()
    Bd_by_pixel_label.grid_forget()
    SNR_num_label.grid_forget()
    minimumTE_num_label.grid_forget()
    minimumTR_num_label.grid_forget()
    warning_label.grid_forget()
    
    minimumTE = 0.00128 * 1000 # *1000 to have it in ms
    minimumTR = 2 * minimumTE
    Time_scan_num_label = Label(frame11, text = "2:16:32", font=("Helvetica", f));         Time_scan_num_label.grid(row = 0, column = 1)
    Bd_by_pixel_label = Label(frame11, text = "390.62", font=("Helvetica", f));            Bd_by_pixel_label.grid(row = 1, column = 1)
    SNR_num_label = Label(frame11, text = "      ", font=("Helvetica", f));                SNR_num_label.grid(row = 2, column = 1)
    minimumTE_num_label = Label(frame11, text = minimumTE, font=("Helvetica", f));         minimumTE_num_label.grid(row = 3, column = 1)
    minimumTR_num_label = Label(frame11, text = minimumTR, font=("Helvetica", f));         minimumTR_num_label.grid(row = 3, column = 1)
    minimumTR_num_label.grid_forget() 
    
    SIM = 0
        
    seq_label2.grid_forget()
    seq_label2 = Label(frame1, text = ' '); seq_label2.grid(row = 10, column = 0)
    
    history_label1.grid_forget(); history_label1 = Label(frame10, text = " ... ", font=("Helvetica", f)); history_label1.grid(row=0, column=0)
    history_label2.grid_forget(); history_label2 = Label(frame10, text = " ... ", font=("Helvetica", f)); history_label2.grid(row=1, column=0)
    history_label3.grid_forget(); history_label3 = Label(frame10, text = " ... ", font=("Helvetica", f)); history_label3.grid(row=2, column=0)
    history_label4.grid_forget(); history_label4 = Label(frame10, text = " ... ", font=("Helvetica", f)); history_label4.grid(row=3, column=0)
    history_label5.grid_forget(); history_label5 = Label(frame10, text = " ... ", font=("Helvetica", f)); history_label5.grid(row=4, column=0)
    history_dict = dict()
    history_list = ['h1','h2','h3','h4','h5']
    
    canvas1.get_tk_widget().destroy()
    canvas2.get_tk_widget().destroy()
    canvas3.get_tk_widget().destroy()
    
    SNR_num_label.grid_forget();SNR_num_label = Label(frame11, text = "  ", font=("Helvetica", f));SNR_num_label.grid(row = 2, column = 1)
    
    T1_import = 0
    T2_import = 0
    T2_star_import = 0
    B0_import = 0
    Delta_B0_import = 0
    SNR_length = 5
    SNR_noise_box_center = [105,20]
    SNR_mean_box_center = [50,80]
    
    frame1.grid_forget()
    frame3.grid_forget()
    frame4.grid_forget()
    frame10.grid_forget()
    frame11.grid_forget()
    frame12.grid_forget()
    Display_postpro_frame_button.grid_forget()
    Display_param_frame_button.grid_forget()
    Display_parvisu_frame_button.grid_forget()
    setting_button.grid_forget()
    field_strength_label.grid_forget()
    field_strength_label_value.grid_forget()
    Sagittal_label.grid_forget()
    Coronal_label.grid_forget()
    Axial_label.grid_forget()
    reset_button.grid_forget()
    continue_to_simu.grid_forget()
    Exit_simu.grid_forget()
    
    imag_size = 5.5
    fig = plt.figure(figsize=(imag_size,imag_size))                
    plt.imshow(img1)           
    plt.axis('off')
    canvas1 = FigureCanvasTkAgg(fig, root)                          
    canvas1.draw()
    canvas1.get_tk_widget().grid(row = 1, column = 1, rowspan = 4)        
    plt.close()

    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(img2)   
    plt.axis('off')
    canvas2 = FigureCanvasTkAgg(fig, root)         
    canvas2.draw()
    canvas2.get_tk_widget().grid(row = 1, column = 2, rowspan = 4)        
    plt.close()
    
    open_simulator.grid(row = 0, column = 1)
    open_generator.grid(row = 0, column = 2)

#///////////// FRAMES, BUTTONS, ENTRIES, LABELS ///////////////////////# 
    
# //////////////////////////////////////////////////////////////////// #
############# Buttons, entry, labels not in any frame ##################
# //////////////////////////////////////////////////////////////////// #    
    
# Font size parameter
f = 12 # labels
e = 10 # entries

# Functions for the settings
global B0_field
B0_field = 'B0_46mT'
    
# Reset button
reset_button = Button(root, text = "Reset images and parameters", font=("Helvetica", f), command = reset)
reset_button.grid(row = 2, column = 0)
     
# Setting button
setting_button = Button(root, text = "Settings", font=("Helvetica", f), command = settings)
setting_button.grid(row = 4, column = 0)

# Field strength label
field_strength_label = Label(root, text = "Magnetic field strength: ", font=("Helvetica", 17))
field_strength_label.grid(row = 5, column = 0)
field_strength_label_value = Label(root, text = "46 mT", font=("Helvetica", 17))
field_strength_label_value.grid(row = 6, column = 0)

# Variables of imported maps
T1_import = 0
T2_import = 0
T2_star_import = 0
B0_import = 0

# //////////////////////////////////////////////////////////////////// #
#################### Frame regarding the sequences ####################
# //////////////////////////////////////////////////////////////////// #
frame1 = LabelFrame(root, text = "Sequences", font=("Helvetica", 15))
frame1.grid(row = 0, column = 0, rowspan = 2)

# Label informing user of predefined sequence choosen
lab = Label(frame1, text = "Select one sequence: ", font=("Helvetica", f)).grid(row = 0, column = 0)

# Sequences buttons
SE_button = Button(frame1, text = "Spin echo", font=("Helvetica", f), activebackground='#999999', command = lambda: seq_simu(SE_button.cget('text')))
SE_button.grid(row = 1, column = 0)
GE_button = Button(frame1, text = "Gradient echo", font=("Helvetica", f), activebackground='#999999', command = lambda: seq_simu(GE_button.cget('text')))
GE_button.grid(row = 2, column = 0)
IN_button = Button(frame1, text = "Inversion recovery", font=("Helvetica", f), activebackground='#999999', command = lambda: seq_simu(IN_button.cget('text')))
IN_button.grid(row = 3, column = 0)
double_IN_button = Button(frame1, text = "Double inversion recovery", font=("Helvetica", f), activebackground='#999999', command = lambda: seq_simu(double_IN_button.cget('text')))
double_IN_button.grid(row = 4, column = 0)
FLAIR_button = Button(frame1, text = "Fluid-attenuated inversion recovery", font=("Helvetica", f), activebackground='#999999', command = lambda: seq_simu(FLAIR_button.cget('text')))
FLAIR_button.grid(row = 5, column = 0)
TSE_button = Button(frame1, text = "Turbo spin echo", font=("Helvetica", f), activebackground='#999999', command = lambda: seq_simu(TSE_button.cget('text')))
TSE_button.grid(row = 6, column = 0)
Dif_button = Button(frame1, text = "Diffusion", font=("Helvetica", f), activebackground='#999999', command = lambda: seq_simu(Dif_button.cget('text')))
Dif_button.grid(row = 7, column = 0)
SSFP_button = Button(frame1, text = "Steady-state", font=("Helvetica", f), activebackground='#999999', command = lambda: seq_simu(SSFP_button.cget('text')))
SSFP_button.grid(row = 8, column = 0)

global Pre_def_seq
Pre_def_seq = " "
seq_label1 = Label(frame1, text = "You have chosen: ", font=("Helvetica", f)).grid(row = 9, column = 0)
seq_label2 = Label(frame1, text = Pre_def_seq); seq_label2.grid(row = 10, column = 0); seq_label2.grid_forget()

# Binding the sequence buttons to open webpages with information on the specific sequences (what is spin echo, TSE, etc...) when right clicking on them 
SE_button.bind('<Button-3>', info_SE)
GE_button.bind('<Button-3>', info_GE)
IN_button.bind('<Button-3>', info_IN)
double_IN_button.bind('<Button-3>', info_double_IN)
FLAIR_button.bind('<Button-3>', info_FLAIR)
TSE_button.bind('<Button-3>', info_TSE)
Dif_button.bind('<Button-3>', info_Dif)
SSFP_button.bind('<Button-3>', info_SSFP)

# //////////////////////////////////////////////////////////////////// #
####### Frame with time of scan, Bd/pixel, SNR, minimum TE or TR #######
# //////////////////////////////////////////////////////////////////// #

# This value is computed from the initial FOV, resolution, bandwidth values ( 128/(50000*2) = 0.00128)
minimumTE = 0.00128 * 1000 # *1000 to have it in ms
minimumTR = 2 * minimumTE

frame11 = LabelFrame(root)
frame11.grid(row = 3, column = 0)

Time_scan_label = Label(frame11, text = "Total scan duration ", font=("Helvetica", f));Time_scan_label.grid(row = 0, column = 0)
Bd_pix_label = Label(frame11, text = "Bandwidth / pixel ", font=("Helvetica", f));     Bd_pix_label.grid(row = 1, column = 0)
SNR_label = Label(frame11, text = "SNR ", font=("Helvetica", f));                      SNR_label.grid(row = 2, column = 0)
minimumTE_label = Label(frame11, text = "TEmin (ms) ", font=("Helvetica", f));         minimumTE_label.grid(row = 3, column = 0)
minimumTR_label = Label(frame11, text = "TRmin (ms) ", font=("Helvetica", f));         minimumTE_label.grid(row = 3, column = 0)
minimumTR_label.grid_forget()

global Time_scan_num_label
global Bd_by_pixel_label
global SNR_num_label
global minimumTE_num_label
global minimumTR_num_label
Time_scan_num_label = Label(frame11, text = "2:16:32", font=("Helvetica", f));         Time_scan_num_label.grid(row = 0, column = 1)
Bd_by_pixel_label = Label(frame11, text = "390.62", font=("Helvetica", f));            Bd_by_pixel_label.grid(row = 1, column = 1)
SNR_num_label = Label(frame11, text = "      ", font=("Helvetica", f));                SNR_num_label.grid(row = 2, column = 1)
minimumTE_num_label = Label(frame11, text = minimumTE, font=("Helvetica", f));         minimumTE_num_label.grid(row = 3, column = 1)
minimumTR_num_label = Label(frame11, text = minimumTR, font=("Helvetica", f));         minimumTR_num_label.grid(row = 3, column = 1)
minimumTR_num_label.grid_forget()

# //////////////////////////////////////////////////////////////////// #
###### Frame regarding the history of the sequences & filter ran #######
# //////////////////////////////////////////////////////////////////// #

frame10 = LabelFrame(root, text = "Sequence & filter history", font=("Helvetica", 15))
frame10.grid(row = 7, column = 0)

global history_label1, history_label2, history_label3, history_label4, history_label5

history_label1 = Label(frame10, text = " ... ", font=("Helvetica", f)); history_label1.grid(row = 0, column = 0)
history_label2 = Label(frame10, text = " ... ", font=("Helvetica", f)); history_label2.grid(row = 1, column = 0)
history_label3 = Label(frame10, text = " ... ", font=("Helvetica", f)); history_label3.grid(row = 2, column = 0)
history_label4 = Label(frame10, text = " ... ", font=("Helvetica", f)); history_label4.grid(row = 3, column = 0)
history_label5 = Label(frame10, text = " ... ", font=("Helvetica", f)); history_label5.grid(row = 4, column = 0)

history_dict = dict()
history_list = ['h1','h2','h3','h4','h5']

# Binding functions to history labels
history_label1.bind('<Button-3>', h1)
history_label2.bind('<Button-3>', h2)
history_label3.bind('<Button-3>', h3)
history_label4.bind('<Button-3>', h4)
history_label5.bind('<Button-3>', h5)

# //////////////////////////////////////////////////////////////////// #
#################### Frame regarding the parameters ####################
# //////////////////////////////////////////////////////////////////// #

# String used to keep track if the parameter frame is up or hidden
global param_frame_on
param_frame_on = 'yes'

# Button showing the parameter frame
Display_param_frame_button = Button(root, text = "Parameters", font=("Helvetica", 18), command = display_param_frame)
Display_param_frame_button.grid(row = 7, column = 1)
Display_param_frame_button.grid_forget()

frame3 = LabelFrame(root, text = "Parameters", font=("Helvetica", 15))
frame3.grid(row = 7, column = 1, columnspan = 2, rowspan = 2)

# Labels of the parameters
# Parameters that will always be shown
readout_label = Label(frame3, text = "Readout gradient ", font=("Helvetica", f));          readout_label.grid(row = 0, column = 1)
phase1_label = Label(frame3, text = "Phase gradient 1", font=("Helvetica", f));            phase1_label.grid(row = 0, column = 3)
phase2_label = Label(frame3, text = "Phase gradient 2", font=("Helvetica", f));            phase2_label.grid(row = 0, column = 5)
FOV1_label = Label(frame3, text = "Field of view (mm) ", font=("Helvetica", f));           FOV1_label.grid(row = 1, column = 0)
FOV2_label = Label(frame3, text = "x", font=("Helvetica", f));                             FOV2_label.grid(row = 1, column = 2)
FOV3_label = Label(frame3, text = "x", font=("Helvetica", f));                             FOV3_label.grid(row = 1, column = 4)
Resolution1_label = Label(frame3, text = "Voxel resolution (mm) ", font=("Helvetica", f)); Resolution1_label.grid(row = 2, column = 0)
Resolution2_label = Label(frame3, text = "x", font=("Helvetica", f));                      Resolution2_label.grid(row = 2, column = 2)
Resolution3_label = Label(frame3, text = "x", font=("Helvetica", f));                      Resolution3_label.grid(row = 2, column = 4)
Data_matrix1_label = Label(frame3, text = "Data matrix ", font=("Helvetica", f));          Data_matrix1_label.grid(row = 3, column = 0)
Data_matrix2_label = Label(frame3, text = "x", font=("Helvetica", f));                     Data_matrix2_label.grid(row = 3, column = 2)
Data_matrix3_label = Label(frame3, text = "x", font=("Helvetica", f));                     Data_matrix3_label.grid(row = 3, column = 4)
Bandwidth_label = Label(frame3, text = "Bandwidth (Hz) ", font=("Helvetica", f));          Bandwidth_label.grid(row = 4, column = 0, pady = 10)
TR_label = Label(frame3, text = "TR (ms) ", font=("Helvetica", f));                        TR_label.grid(row = 5, column = 0)

# Parameters dependent on the sequence
TE_label = Label(frame3, text = "TE (ms) ", font=("Helvetica", f));                        TE_label.grid(row = 6, column = 0)
TI_label = Label(frame3, text = "TI (ms) ", font=("Helvetica", f));                        TI_label.grid(row = 7, column = 0)
TI2_label = Label(frame3, text = "TI2 (ms) ", font=("Helvetica", f));                      TI2_label.grid(row = 7, column = 2)
Alpha_label = Label(frame3, text = "Alpha ", font=("Helvetica", f));                       Alpha_label.grid(row = 8, column = 0)
ETL_label = Label(frame3, text = "ETL ", font=("Helvetica", f));                           ETL_label.grid(row = 8, column = 0)
Kspace_traj_label = Label(frame3, text = "Kspace trajectory ", font = ("Helvetica", f));   Kspace_traj_label.grid(row = 8, column = 1)    
TEeff_label = Label(frame3, text = "TEeff (ms) ", font = ("Helvetica", f));                TEeff_label.grid(row = 8, column = 2) 
TEeff_val_label = Label(frame3, text = "80.0", font = ("Helvetica", f));                   TEeff_val_label.grid(row = 8, column = 3)  
G_label = Label(frame3, text = "G (mT/mm) ", font=("Helvetica", f));                       G_label.grid(row = 10, column = 0)
smalldelta_label = Label(frame3, text = "small delta (ms) ", font=("Helvetica", f));       smalldelta_label.grid(row = 11, column = 0)
bigdelta_label = Label(frame3, text = "big delta (ms) ", font=("Helvetica", f));           bigdelta_label.grid(row = 12, column = 0)
B_label = Label(frame3, text = "b (s/mm2) ", font=("Helvetica", f));                       B_label.grid(row = 10, column = 3)
TEmin_label = Label(frame3, text = "TEmin diffusion (ms) ", font=("Helvetica", f));        TEmin_label.grid(row = 11, column = 3)
Bval_label = Label(frame3, text = "302.18", font=("Helvetica", f));                        Bval_label.grid(row = 8, column = 4)
TEminval_label = Label(frame3, text = "3.0", font=("Helvetica", f));                       TEminval_label.grid(row = 9, column = 4)

# Entries of the parameters
global entry_width
entry_width = 7
# Parameters that will always be shown
FOV1_entry = Entry(frame3, width=entry_width, font=("Helvetica", e)); FOV1_entry.grid(row = 1, column = 1); FOV1_entry.insert(0, '250')
FOV2_entry = Entry(frame3, width=entry_width, font=("Helvetica", e)); FOV2_entry.grid(row = 1, column = 3); FOV2_entry.insert(0, '300')
FOV3_entry = Entry(frame3, width=entry_width, font=("Helvetica", e)); FOV3_entry.grid(row = 1, column = 5); FOV3_entry.insert(0, '275')
FOV_check = 20625000 # this value is 250*300*275 and will be used to see if the FOV was changed by the user
Res1_entry = Entry(frame3, width=entry_width, font=("Helvetica", e)); Res1_entry.grid(row = 2, column = 1); Res1_entry.insert(0, '1.95')
Res2_entry = Entry(frame3, width=entry_width, font=("Helvetica", e)); Res2_entry.grid(row = 2, column = 3); Res2_entry.insert(0, '2.34')
Res3_entry = Entry(frame3, width=entry_width, font=("Helvetica", e)); Res3_entry.grid(row = 2, column = 5); Res3_entry.insert(0, '2.15')
Res_check = 9.81045 # this value is 1.95*2.34*2.15 and will be used to see if the resolution was changed by the user
Data_mat1_entry = Entry(frame3, width=7, font=("Helvetica", e)); Data_mat1_entry.grid(row = 3, column = 1); Data_mat1_entry.insert(0, '128')
Data_mat2_entry = Entry(frame3, width=entry_width, font=("Helvetica", e)); Data_mat2_entry.grid(row = 3, column = 3); Data_mat2_entry.insert(0, '128')
Data_mat3_entry = Entry(frame3, width=entry_width, font=("Helvetica", e)); Data_mat3_entry.grid(row = 3, column = 5); Data_mat3_entry.insert(0, '128')
Data_matrix_check = 2097152 # this value is 128*128*128 and will be used to see if the Data matrix size was changed by the user
Bandwidth_entry=Entry(frame3,width=entry_width,font=("Helvetica", e));Bandwidth_entry.grid(row=4,column=1); Bandwidth_entry.insert(0, '40000')
TR_entry = Entry(frame3, width=entry_width, font=("Helvetica", e)); TR_entry.grid(row = 5, column = 1);   TR_entry.insert(0, '500')

# Parameters dependent on the sequence
TE_entry = Entry(frame3, width=entry_width, font=("Helvetica", e)); TE_entry.grid(row = 6, column = 1);         TE_entry.insert(0, '20')
TI_entry = Entry(frame3, width=entry_width, font=("Helvetica", e)); TI_entry.grid(row = 7, column = 1);         TI_entry.insert(0, '250')
TI2_entry = Entry(frame3, width=entry_width, font=("Helvetica", e));TI2_entry.grid(row = 7, column = 3);        TI2_entry.insert(0, '250')
Alpha_entry = Entry(frame3, width=entry_width, font=("Helvetica", e)); Alpha_entry.grid(row = 8, column = 1);      Alpha_entry.insert(0, '45')
ETL_entry = Entry(frame3, width=entry_width, font=("Helvetica", e)); ETL_entry.grid(row = 9, column = 1);        ETL_entry.insert(0, '8')
G_entry = Entry(frame3, width=entry_width, font=("Helvetica", e)); G_entry.grid(row = 10, column = 1);          G_entry.insert(0, '10')
smalldelta_entry = Entry(frame3,width=entry_width,font=("Helvetica",e));smalldelta_entry.grid(row=11,column=1);smalldelta_entry.insert(0, '1')
bigdelta_entry = Entry(frame3,width=entry_width,font=("Helvetica", e));bigdelta_entry.grid(row =12,column=1);  bigdelta_entry.insert(0, '2')
  
# Dropdown menu to select kspace trajectory in TSE sequence   
options = [
    "Linear",
    "In-out",
    "Out-in"
]
traj = StringVar()
traj.set(options[0])                                           # Default value, or could use; options[0]
tse_drop = OptionMenu(frame3, traj, *options, command = showTEeff) # Using a list, NEEDS a star in front
tse_drop.grid(row = 8, column = 1)

global Readout_axis
Readout_axis = 'FH'

# Dropdown menu to select readout axis in TSE sequence   
options = [
    "Along Sagittal Foot/head",
    "Along Sagittal Anterior/posterior",
    "Along Coronal Foot/head",
    "Along Coronal Left/right",
    "Along Axial Anterior/posterior",
    "Along Axial Left/right"
]
Readout = StringVar()
Readout.set(options[0])                                            # Default value, or could use; options[0]
tse_read_drop = OptionMenu(frame3, Readout, *options, command = readout_selection) # Using a list, NEEDS a star in front
tse_read_drop.grid(row = 5, column = 5)
tse_read_label = Label(frame3, text = "Readout axe", font=("Helvetica", f)); tse_read_label.grid(row = 5, column = 4)
Post_TSE_TEeff_label = Label(frame3, text = "TEeff = ", font=("Helvetica", f)); Post_TSE_TEeff_label.grid(row = 15, column = 2)
Post_TSE_TEeff_label_val = Label(frame3, text = "20.0", font=("Helvetica", f)); Post_TSE_TEeff_label_val.grid(row = 15, column = 3)

# Hide the labels and the entries of specific sequences, they will appear only if the sequence where they are involed is selected
# (Label ---- Entry ---- Dropdown)
TE_label.grid_forget();         TE_entry.grid_forget();           tse_drop.grid_forget()
TI_label.grid_forget();         TI_entry.grid_forget();           tse_read_drop.grid_forget()
TI2_label.grid_forget();        TI2_entry.grid_forget()
Alpha_label.grid_forget();      Alpha_entry.grid_forget()
ETL_label.grid_forget();        ETL_entry.grid_forget()
Kspace_traj_label.grid_forget()
TEeff_label.grid_forget()
TEeff_val_label.grid_forget()
G_label.grid_forget();          G_entry.grid_forget()
smalldelta_label.grid_forget(); smalldelta_entry.grid_forget()
bigdelta_label.grid_forget();   bigdelta_entry.grid_forget()
B_label.grid_forget()
Bval_label.grid_forget()
TEmin_label.grid_forget()       
TEminval_label.grid_forget()
tse_read_label.grid_forget()
Post_TSE_TEeff_label.grid_forget()
Post_TSE_TEeff_label_val.grid_forget()

# Those values are computed from the G/small delta/big delta values
B = 302.18
dif_TEmin = np.divide(3.0,1000)

# //////////////////////////////////////////////////////////////////// #
#################### Frame of the post sequence TSE ####################
# //////////////////////////////////////////////////////////////////// #
frame5 = LabelFrame(frame3, text = "TSE readout", font=("Helvetica", 15))
frame5.grid(row = 4, column = 2, rowspan = 3 , columnspan = 4)

# Labels
Post_seq_TSE_TE_label = Label(frame5, text = "TE (ms) ", font=("Helvetica", f))          
Post_seq_TSE_TE_label.grid(row = 0, column = 0)
Post_seq_TSE_ETL_label = Label(frame5, text = "ETL ", font=("Helvetica", f))            
Post_seq_TSE_ETL_label.grid(row = 1, column = 0)
post_seq_TSE_traj_label = Label(frame5, text = "Kspace trajectory ", font = ("Helvetica", f))    
post_seq_TSE_traj_label.grid(row = 0, column = 2) 

# Dropdown menu to select kspace trajectory in post sequence TSE   
options = [
    "Linear",
    "In-out",
    "Out-in"
]
post_traj = StringVar()
post_traj.set(options[1])                                           # Default value, or could use; options[0]
post_tse_drop = OptionMenu(frame5, post_traj, *options, command = post_showTEeff) # Using a list, NEEDS a star in front
post_tse_drop.grid(row = 1, column = 2)
# Entries
post_TSE_TE_entry = Entry(frame5,width=entry_width,font=("Helvetica",e));post_TSE_TE_entry.grid(row=0,column=1); post_TSE_TE_entry.insert(0, '20')
post_TSE_ETL_entry = Entry(frame5,width=entry_width,font=("Helvetica",e));post_TSE_ETL_entry.grid(row=1,column=1); post_TSE_ETL_entry.insert(0, '8')

options = [
    "Along Sagittal Foot/head",
    "Along Sagittal Anterior/posterior",
    "Along Coronal Foot/head",
    "Along Coronal Left/right",
    "Along Axial Anterior/posterior",
    "Along Axial Left/right"
]

post_tse_read_label = Label(frame5, text = "Readout axe", font=("Helvetica", f)); post_tse_read_label.grid(row = 2, column = 0, columnspan = 2)
Post_Readout = StringVar()
Post_Readout.set(options[0])                                            # Default value, or could use; options[0]
post_tse_read_drop = OptionMenu(frame5, Post_Readout, *options, command = post_readout_selection) # Using a list, NEEDS a star in front
post_tse_read_drop.grid(row = 2, column = 2)
frame5.grid_forget()

# Button for applying the parameters and simulate the chosen sequence
Run_seq_button = Button(frame3, text = "Run sequence", font=("Helvetica", 18), command = run).grid(row = 15, column = 5)

# Binding function to open webpage that offers the user information on the specific label (what is TR, TE, etc...)    
TR_label.bind('<Button-3>', info_tr_te)
TE_label.bind('<Button-3>', info_tr_te)
TI_label.bind('<Button-3>', info_ti)
TI2_label.bind('<Button-3>', info_ti2)
Alpha_label.bind('<Button-3>', info_alpha)
ETL_label.bind('<Button-3>', info_etl)
Kspace_traj_label.bind('<Button-3>', info_Kspace_traj)
TEeff_label.bind('<Button-3>', info_teeff)
G_label.bind('<Button-3>', info_G)
smalldelta_label.bind('<Button-3>', info_smalldelta)
bigdelta_label.bind('<Button-3>', info_bigdelta)
B_label.bind('<Button-3>', info_B)
Bandwidth_label.bind('<Button-3>', info_bd)
FOV1_label.bind('<Button-3>', info_fov)
Resolution1_label.bind('<Button-3>', info_res)
SNR_label.bind('<Button-3>', info_snr)

# ///////////////////////////////////////////////////////////////////////// #
#################### Frame regarding the post-processing ####################
# ///////////////////////////////////////////////////////////////////////// #

# String used to keep track if the post processing frame is up or hidden
global post_process_frame_on
post_process_frame_on = 'no'

# Button showing the parameter frame
Display_postpro_frame_button = Button(root, text = "Post-processing", font=("Helvetica", 18), command = display_postpro_frame)
Display_postpro_frame_button.grid(row = 7, column = 3)

frame4 = LabelFrame(root, text = "Post-processing", font=("Helvetica", 15))
frame4.grid(row = 7, column = 3, rowspan = 2)

# Labels and buttons of post-processing options
Postprocess_label = Label(frame4, text = "Post-processing options", font=("Helvetica", f)).grid(row = 0, column = 0)
# Labels
low_pass_label = Label(frame4, text = "Kernel size --> ", font=("Helvetica", f)).grid(row = 1, column = 0)
high_pass_label = Label(frame4, text = "Kernel size (must be odd)--> ", font=("Helvetica", f)).grid(row = 2, column = 0)
Gauss_label = Label(frame4, text = "Std of kernel filter --> ", font=("Helvetica", f)).grid(row = 3, column = 0)
non_local_label = Label(frame4, text = "h, patch size and distance --> ", font=("Helvetica", f)).grid(row = 4, column = 0)
# Entries
low_pass_entry = Entry(frame4,width=entry_width,font=("Helvetica", e)); low_pass_entry.grid(row = 1, column = 1, columnspan = 3); low_pass_entry.insert(0,'3')
high_pass_entry = Entry(frame4,width=entry_width,font=("Helvetica", e));high_pass_entry.grid(row = 2, column = 1, columnspan =3); high_pass_entry.insert(0,'3')
Gauss_entry = Entry(frame4,width=entry_width,font=("Helvetica", e));Gauss_entry.grid(row = 3, column = 1, columnspan = 3); Gauss_entry.insert(0, '0.5')
Non_local_h_entry = Entry(frame4, font=("Helvetica", e), width = 5);     Non_local_h_entry.grid(row = 4, column = 1); Non_local_h_entry.insert(0, '0.8')
Non_local_psize_entry = Entry(frame4, font=("Helvetica", e), width = 5); Non_local_psize_entry.grid(row = 4, column = 2); Non_local_psize_entry.insert(0, '7')
Non_local_pdist_entry = Entry(frame4, font=("Helvetica", e), width = 5); Non_local_pdist_entry.grid(row = 4, column = 3); Non_local_pdist_entry.insert(0, '11')
# Buttons (lambda is to use the parentheses and arguments)
Low_pass_button = Button(frame4, text = "Low-pass filter", font=("Helvetica", f), command = lambda: lowpass(low_pass_entry.get()))
Low_pass_button.grid(row = 1, column = 4)
High_pass_button = Button(frame4, text = "High-pass filter", font=("Helvetica", f), command = lambda: highpass(high_pass_entry.get()))
High_pass_button.grid(row = 2, column = 4)
Gaussian_button = Button(frame4, text = "Gaussian filter", font=("Helvetica", f), command = lambda: gauss(Gauss_entry.get()))
Gaussian_button.grid(row = 3, column = 4) 
Non_local_button = Button(frame4, text = "Non local filter", font=("Helvetica", f), command = lambda: non_local(Non_local_h_entry.get(), Non_local_psize_entry.get(), Non_local_pdist_entry.get()))
Non_local_button.grid(row = 4, column = 4) 

# string indicating if the 'undo filter function' was used or not
undo = 'no'

fil_undo_button = Button(frame4, text = "Reset sequence image unfiltered", font=("Helvetica", f), command = filter_undo)
fil_undo_button.grid(row = 0, column = 1, columnspan = 4)

# Button to see how the SNR is currently defined   
SNR_visu_button = Button(frame4, text = "Visualizing SNR", font=("Helvetica", f), command = SNR_visu).grid(row = 5, column = 0)

# Button to define the boxes used in the computation of the SNR
SNR_measure_button = Button(frame4, text = "Measuring SNR", font=("Helvetica", f), command = SNR_measure).grid(row = 5, column = 1, columnspan =3)

SNR_length = 5
SNR_noise_box_center = [105,20]
SNR_mean_box_center = [50,80]

# Binding function to open webpage that offers the user information on the specific label (what is TR, TE, etc...?)    
Low_pass_button.bind('<Button-3>', info_lowpass)
High_pass_button.bind('<Button-3>', info_highpass)
Gaussian_button.bind('<Button-3>', info_gausspass)
Non_local_button.bind('<Button-3>', info_nonlocal)

# ///////////////////////////////////////////////////////////////////////// #
############### Frame regarding the parameters visualisation ################
# ///////////////////////////////////////////////////////////////////////// #

# String used to keep track if the post processing frame is up or hidden
global param_visu_frame_on
param_visu_frame_on = 'no'

# Button showing the parameter frame
Display_parvisu_frame_button = Button(root, text = "Parameter visualisation", font=("Helvetica", 18), command = display_visu_frame)
Display_parvisu_frame_button.grid(row = 9, column = 3)

frame12 = LabelFrame(root, text = "Parameter visualisation", font=("Helvetica", 15))
frame12.grid(row = 9, column = 3)

Filter_vis_button = Label(frame12, text = "Filter visualization", font=("Helvetica", f)).grid(row = 0, column = 0)
# Dropdown menu to select filter to visualize   
options = [
    "low",
    "high",
    "gauss"
]
fi = StringVar()
fi.set("No selection")                                               # Default value, or could use; options[0]
filter_drop = OptionMenu(frame12, fi, *options, command = filter_vis) # Using a list, NEEDS a star in front
filter_drop.grid(row = 0, column = 1)

Param_vis_button = Button(frame12, text = "Parameters visualization", font=("Helvetica", f), command = param_vis).grid(row = 0, column = 2)

# Binding the frames to their closing function
frame3.bind('<Button-3>', close_param_display)    
frame4.bind('<Button-3>', close_postpro_display) 
frame12.bind('<Button-3>', close_parvisu_display)

# //////////////////////////////////////////////////////////////// #
#################### Frame regarding the saving ####################
# //////////////////////////////////////////////////////////////// #

frame2 = LabelFrame(frame3, text = "Saving images", font=("Helvetica", f))
frame2.grid(row = 15, column = 0, columnspan = 2)

Axe_save_label = Label(frame2, text = "sagittal, coronal, or axial ", font=("Helvetica", f)); Axe_save_label.grid(row = 0, column = 0)
Axe_save_entry = Entry(frame2, font=("Helvetica", e)); Axe_save_entry.grid(row = 1, column = 0); Axe_save_entry.insert(0,'sagittal')

Sag_save_button = Button(frame2, text="Save", command = lambda: save_img(Axe_save_entry.get()))
Sag_save_button.grid(row = 1, column = 1)

er_label = Label(frame2, text = "Misspelled axe or no sequence selected! ", font=("Helvetica", 10)); er_label.grid(row = 2, column = 0); 
er_label.grid_forget();

# Warning label if there is a paradox in the parameters
warning_label = Label(root, text = "  ", font=("Helvetica", 15)); warning_label.grid(row = 8, column = 0); warning_label.grid_forget()

# Exit button
Exit_button = Button(root, text = "Exit program", font=("Helvetica", f), command = root.destroy)
Exit_button.grid(row = 10, column = 0) 

# Simulator image labels
Sagittal_label = Label(root, text = "Sagittal", font=("Helvetica", 15)); Sagittal_label.grid(row = 0, column = 1)
Coronal_label = Label(root, text = "Coronal", font=("Helvetica", 15)); Coronal_label.grid(row = 0, column = 2)
Axial_label = Label(root, text = "Axial", font=("Helvetica", 15)); Axial_label.grid(row = 0, column = 3)

frame1.grid_forget()
frame3.grid_forget()
frame4.grid_forget()
frame10.grid_forget()
frame11.grid_forget()
frame12.grid_forget()
Display_postpro_frame_button.grid_forget()
Display_parvisu_frame_button.grid_forget()
setting_button.grid_forget()
Sagittal_label.grid_forget()
Coronal_label.grid_forget()
Axial_label.grid_forget()
reset_button.grid_forget()
field_strength_label.grid_forget()
field_strength_label_value.grid_forget()

#// Classe creating pop up windows when mouse is over button or label //#     
class ToolTip(object):
    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0
    def showtip(self, text):
        "Display text in tooltip window"
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 57
        y = y + cy + self.widget.winfo_rooty() +27
        self.tipwindow = tw = Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        label = Label(tw, text=self.text, justify=LEFT,
                      background="#ffffe0", relief=SOLID, borderwidth=1,
                      font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)
    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

def CreateToolTip(widget, text):
    toolTip = ToolTip(widget)
    def enter(event):
        toolTip.showtip(text)
    def leave(event):
        toolTip.hidetip()
    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)
    
#//// Functions that prints small window beneath the mouse, \n make a line break ///
# Sequence frame
CreateToolTip(SE_button, text = 'Spin echo sequence. Usual parameter range are:.\n'
                 'TR ≈ 3000ms, TE ≈ 160ms\n'
                 '(For more information right click)')
CreateToolTip(GE_button, text = 'Gradient echo sequence. Usual parameter range are:.\n'
                 'TR ≈ 50ms, TE ≈ 4ms, Alpha angle ≈ 30°\n'
                 '(For more information right click)')
CreateToolTip(IN_button, text = 'Inversion recovery sequence. Usual parameter range are:.\n'
                 'TR ≈ 5000ms, TE ≈ 50ms, TI ≈ 400ms\n'
                 '(For more information right click)')
CreateToolTip(double_IN_button, text = 'Double inversion recovery sequence. Usual parameter range are:.\n'
                 'TR ≈ 5000ms, TE ≈ 160ms, TI ≈ 650ms, TI2 ≈ 150ms\n'
                 '(For more information right click)')
CreateToolTip(FLAIR_button, text = 'FLAIR sequence. Usual parameter range are:.\n'
                 'TR ≈ 11000ms, TE ≈ 20ms\n'
                 '(For more information right click)')
CreateToolTip(TSE_button, text = 'Turbo spin echo sequence. Usual parameter range are:.\n'
                 'TR ≈ 2000ms, TE ≈ 20ms, ETL ≈ 6\n'
                 '(For more information right click)')
CreateToolTip(Dif_button, text = 'Diffusion sequence. Usual parameter range are:.\n'
                 'TR ≈ 3000ms, TE ≈ 100ms, b ≈ 483s/mm²\n'
                 '(For more information right click)')
CreateToolTip(SSFP_button, text = 'Steady-state free precession sequence. Usual parameter range are:.\n'
                 'TR ≈ 16ms, Alpha angle ≈ 45°\n'
                 '(For more information right click)')
# Parameter frame
CreateToolTip(TR_label, text = 'Repetition time.\n'
                 '(For more information right click)')
CreateToolTip(TE_label, text = 'Echo time.\n'
                 '(For more information right click)')
CreateToolTip(TI_label, text = 'Inversion time.\n'
                 '(For more information right click)')
CreateToolTip(TI2_label, text = 'Second inversion time.\n'
                 '(For more information right click)')
CreateToolTip(Alpha_label, text = 'Flip angle.\n'
                 '(For more information right click)')
CreateToolTip(ETL_label, text = 'Echo train length.\n'
                 '(For more information right click)')
CreateToolTip(Kspace_traj_label, text = 'Trajectory used to fill the kspace.\n'
                 '(For more information right click)')
CreateToolTip(TEeff_label, text = 'Effective echo time.\n'
                 '(For more information right click)')
CreateToolTip(G_label, text = 'Gradient strength used in the diffusion sequence.\n'
                 '(For more information right click)')
# Scan characteristics
CreateToolTip(Bandwidth_label, text = 'The bandwidth of an mri image usualy refers to either the range of frequencies\n'
                 'used in the RF pulses (transmission of signal), or to the range of frequencies used to recieve the signal.'
                 '(For more information right click)')
CreateToolTip(FOV1_label, text = 'The field of view refers to the distance over which an image is acquired.\n'
                 '(For more information right click)')
CreateToolTip(Resolution1_label, text = 'Resolution of an mri image refers to the size of a voxel,\n'
                 'or the amount of voxel in the field of view (FOV).\n'
                 '(For more information right click)')
CreateToolTip(SNR_label, text = 'Signal to noise ratio. The SNR is used to illustrate image quality.\n'
                 '(For more information right click)')
# Post-processing fram
CreateToolTip(Low_pass_button, text = 'Low pass filter cutting high frequencies, resulting in a smoother image with less details, less noise, but could be blurry.\n'
                 '(For more information right click)')
CreateToolTip(High_pass_button, text = 'Hig pass filter cutting low frequencies, resulting in an image with more details (very sharp), but losses a bit of its smoothness shape structure.\n'
                 '(For more information right click)')
CreateToolTip(Gaussian_button, text = 'Low pass filter based on the gaussian function.\n'
                 '(For more information right click)')
CreateToolTip(Non_local_button, text = 'Low pass denoising filter based on the mean of group of pixels/voxel not necessarly neighbors of the targeted pixel/voxel.\n'
                 '(For more information right click)')
   
# ///////////////////////////////////////////////////////////////////////// #
################################ Generator ##################################
# ///////////////////////////////////////////////////////////////////////// #

#///////////////////////////// FUNCTIONS ///////////////////////////////////#          
#/////////// Functions showing the effective TE for TSE and post-TSE sequences //////////////#  
def showTEeff_gene(*args):
    """
    Computes and prints the effective TE for the TSE sequence when the k-space trajectory is changed
    Input : None
    Output: None
    """
    global TEeff_val_label_gene
    
    eff = 0
    TE = TE_entry_gene.get(); TE = int(TE); TE = np.divide(TE,1000)
    ETL = int(ETL_entry_gene.get()) 
    trajectory_gene = traj_gene.get()
    
    if trajectory_gene == 'Linear':
        eff = 0.5 * TE * ETL
    elif trajectory_gene == 'In-out':
        eff = TE
    elif trajectory_gene == 'Out-in':
        eff = TE * ETL   
        
    TEeff_val_label_gene.grid_forget() 
    TEeff_val_label_gene = Label(frame7, text = str(int(1000*eff)), font = ("Helvetica", f))
    TEeff_val_label_gene.grid(row = 8, column = 3) 
    
def post_showTEeff_gene(*args):
    """
    Computes and prints the effective TE for the post TSE sequence when the k-space trajectory is changed
    Input : None
    Output: None
    """
    global TEeff_val_label_gene
    post_ETL = int(post_TSE_ETL_entry_gene.get()) 
    post_TE = np.divide(int(post_TSE_TE_entry_gene.get()),1000)
    posttraj = post_traj_gene.get()
    eff = 0
    
    if posttraj == 'Linear':
        eff = 0.5 * post_TE * post_ETL
    elif posttraj == 'In-out':
        eff = post_TE
    elif posttraj == 'Out-in':
        eff = post_TE * post_ETL   
        
    TEeff_val_label_gene.grid_forget() 
    TEeff_val_label_gene = Label(frame7, text = str(int(1000*eff)), font = ("Helvetica", f))
    if Pre_def_seq_gene == 'Dif':
        TEeff_val_label_gene.grid(row = 10, column = 3)
    elif Pre_def_seq_gene == 'Double IN':
        TEeff_val_label_gene.grid(row = 9, column = 3)
    else:
        TEeff_val_label_gene.grid(row = 8, column = 3) 

#/////// Functions updating the widgets for the specific sequences /////////#  
def remove_widgets_gene(seq):
    """
    Removes the widgets (labels and entries) that were needed for previously selected sequence
    Input : seq --> string of Pre_def_seq
    Output: None
    """
    if seq == 'SE':
        TE_entry_gene.grid_forget();         TE_label_gene.grid_forget()
    elif seq == 'GE':
        TE_entry_gene.grid_forget();         TE_label_gene.grid_forget()
        Alpha_entry_gene.grid_forget();      Alpha_label_gene.grid_forget()
    elif seq == 'IN':
        TI_entry_gene.grid_forget();         TI_label_gene.grid_forget()
        frame15.grid_forget()
        TEeff_label_gene.grid_forget()
        TEeff_val_label_gene.grid_forget()
    elif seq == 'Double IN':
        TI_entry_gene.grid_forget();         TI_label_gene.grid_forget()
        TI2_entry_gene.grid_forget();        TI2_label_gene.grid_forget()
        frame15.grid_forget()
        TEeff_label_gene.grid_forget()
        TEeff_val_label_gene.grid_forget()
    elif seq == 'FLAIR':
        TI_entry_gene.grid_forget();         TI_label_gene.grid_forget()
        frame15.grid_forget()
        TEeff_label_gene.grid_forget()
        TEeff_val_label_gene.grid_forget()
    elif seq == 'Dif':
        G_label_gene.grid_forget();          G_entry_gene.grid_forget();        
        smalldelta_label_gene.grid_forget(); smalldelta_entry_gene.grid_forget()
        bigdelta_label_gene.grid_forget();   bigdelta_entry_gene.grid_forget()
        B_label_gene.grid_forget()
        Bval_label_gene.grid_forget()
        TEmin_label_gene.grid_forget()   
        TEminval_label_gene.grid_forget()     
        frame15.grid_forget()
        TEeff_label_gene.grid_forget()
        TEeff_val_label_gene.grid_forget()
    elif seq == 'TSE':
        TE_entry_gene.grid_forget();         TE_label_gene.grid_forget()
        ETL_entry_gene.grid_forget();        ETL_label_gene.grid_forget()
        Kspace_traj_label_gene.grid_forget();tse_drop_gene.grid_forget()
        TEeff_label_gene.grid_forget()
        TEeff_val_label_gene.grid_forget()
        tse_read_label_gene.grid_forget();   tse_read_drop_gene.grid_forget()
    elif seq == 'SSFP':
        Alpha_entry_gene.grid_forget();      Alpha_label_gene.grid_forget()
        
def show_widgets_gene(seq):
    """
    Show the widgets (labels and entries) and the minimum TE (or TR) that are needed for the selected sequence 
    Input : seq --> string of Pre_def_seq
    Output: None
    """     
    global TR
    global TE
    global FOV
    global Data_mat
    global Resolution
    global Bandwidth
    global minimumTE_num_label
    
    TR = TR_entry_gene.get(); TR = int(TR); TR = np.divide(TR,1000) # Divide by 1000 to have the values in milli
    TE = TE_entry_gene.get(); TE = int(TE); TE = np.divide(TE,1000)
    fov1 = FOV1_entry_gene.get(); fov1 = int(fov1)
    fov2 = FOV2_entry_gene.get(); fov2 = int(fov2)
    fov3 = FOV3_entry_gene.get(); fov3 = int(fov3)
    res1 = Res1_entry_gene.get(); res1 = float(res1) 
    res2 = Res2_entry_gene.get(); res2 = float(res2)
    res3 = Res3_entry_gene.get(); res3 = float(res3)
    bd = Bandwidth_entry_gene.get(); Bandwidth = int(bd)
    FOV = [fov1, fov2, fov3]
    Resolution = [res1, res2, res3]
    Data_mat = np.divide(FOV, Resolution)
    Data_mat = [int(Data_mat[0]), int(Data_mat[1]), int(Data_mat[2])]
    
    minimumTE = Data_mat[0]/(2*Bandwidth)
    minimumTR = 2*minimumTE
    
    if seq == 'SE':
        TE_label_gene.grid(row = 6, column = 0);         TE_entry_gene.grid(row = 6, column = 1)
    elif seq == 'GE':
        TE_label_gene.grid(row = 6, column = 0);         TE_entry_gene.grid(row = 6, column = 1)
        Alpha_label_gene.grid(row = 7, column = 0);      Alpha_entry_gene.grid(row = 7, column = 1)
    elif seq == 'IN':
        TI_label_gene.grid(row = 6, column = 0);         TI_entry_gene.grid(row = 6, column = 1)
        frame15.grid(row = 8, column = 4, rowspan = 1, columnspan = 2)
        TEeff_label_gene.grid(row = 8, column = 2)
        TEeff_val_label_gene.grid(row = 8, column = 3)
        post_showTEeff_gene()
    elif seq == 'Double IN':
        TI_label_gene.grid(row = 7, column = 0);         TI_entry_gene.grid(row = 7, column = 1)
        TI2_label_gene.grid(row = 8, column = 0);        TI2_entry_gene.grid(row = 8, column = 1)
        frame15.grid(row = 9, column = 4, rowspan = 2, columnspan = 2)
        TEeff_label_gene.grid(row = 9, column = 2)
        TEeff_val_label_gene.grid(row = 9, column = 3)
        post_showTEeff_gene()
    elif seq == 'FLAIR':
        TI_label_gene.grid(row = 7, column = 0);         TI_entry_gene.grid(row = 7, column = 1) ; TI_entry_gene.delete(0,END); TI_entry_gene.insert(0, '1700')  
        frame15.grid(row = 8, column = 4, rowspan = 1, columnspan = 2)
        TEeff_label_gene.grid(row = 8, column = 2)
        TEeff_val_label_gene.grid(row = 8, column = 3)
        post_showTEeff_gene()
    elif seq == 'Dif':
        G_label_gene.grid(row = 7, column = 0);          G_entry_gene.grid(row = 7, column = 1);  
        smalldelta_label_gene.grid(row = 8, column = 0); smalldelta_entry_gene.grid(row = 8, column = 1)
        bigdelta_label_gene.grid(row = 9, column = 0);   bigdelta_entry_gene.grid(row = 9, column = 1)
        B_label_gene.grid(row = 8, column = 3)
        Bval_label_gene.grid(row = 8, column = 4)
        TEmin_label_gene.grid(row = 9, column = 3)     
        TEminval_label_gene.grid(row = 9, column = 4)
        frame15.grid(row = 10, column = 4, rowspan = 1, columnspan = 2)
        TEeff_label_gene.grid(row = 10, column = 2)
        TEeff_val_label_gene.grid(row = 10, column = 3)
        post_showTEeff_gene()
    elif seq == 'TSE':
        TE_label_gene.grid(row = 6, column = 0);         TE_entry_gene.grid(row = 6, column = 1)
        ETL_label_gene.grid(row = 7, column = 0);        ETL_entry_gene.grid(row = 7, column = 1)
        Kspace_traj_label_gene.grid(row = 8, column = 0);tse_drop_gene.grid(row = 8, column = 1)
        TEeff_label_gene.grid(row = 8, column = 2)
        TEeff_val_label_gene.grid(row = 8, column = 5)
        tse_read_label_gene.grid(row = 6, column = 4);   tse_read_drop_gene.grid(row = 6, column = 5)
        showTEeff_gene()
    elif seq == 'SSFP':
        Alpha_label_gene.grid(row = 6, column = 0);      Alpha_entry_gene.grid(row = 6, column = 1)      
        
    if seq == 'SSFP':
        minimumTE_label_gene = Label(frame11, text = "TRmin (ms) ", font=("Helvetica", 12)).grid(row = 3, column = 0)
        minimumTE_num_label_gene = Label(frame11, text = np.round(minimumTR * 1000,2), font=("Helvetica", 12))
        minimumTE_num_label_gene.grid(row = 3, column = 1)
    else:
        minimumTE_label_gene = Label(frame11, text = "TEmin (ms)", font=("Helvetica", 12)).grid(row = 3, column = 0)
        minimumTE_num_label_gene = Label(frame11, text = np.round(minimumTE * 1000,2), font=("Helvetica", 12))
        minimumTE_num_label_gene.grid(row = 3, column = 1)

#/////////// Functions updating the sequence variable //////////////#
def seq_gene(name):
    """
    Function informing the user of selected sequence and updates the widgets (labels and entries) for the specific sequence 
    Input : seq --> string of Pre_def_seq
    Output: None
    """
    global Pre_def_seq_gene
    global seq_name
    global seq_label2_gene
    remove_widgets_gene(Pre_def_seq_gene)
    
    seq_name = name
    if name == 'Spin echo':
        Pre_def_seq_gene = 'SE'
    elif name == 'Gradient echo':
        Pre_def_seq_gene = 'GE'
    elif name == 'Inversion recovery':
        Pre_def_seq_gene = 'IN'
    elif name == 'Double inversion recovery':
        Pre_def_seq_gene = 'Double IN'
    elif name == 'Fluid-attenuated inversion recovery':
        Pre_def_seq_gene = 'FLAIR'
    elif name == 'Steady-state ':
        Pre_def_seq_gene = 'SSFP'
    elif name == 'Diffusion':
        Pre_def_seq_gene = 'Dif'
    elif name == 'Turbo spin echo':
        Pre_def_seq_gene = 'TSE'
    elif name == 'Steady-state':
        Pre_def_seq_gene = 'SSFP'
    
    seq_label2_gene.grid_forget()
    seq_label2_gene = Label(frame6, text = name, font=("Helvetica", 18))
    seq_label2_gene.grid(row = 10, column = 0, rowspan = 2)
    show_widgets_gene(Pre_def_seq_gene)
            
#/////// Functions simulating a high and low SNR MRI sequence /////////#     
def create_3D_low_and_high_SNR_data(FOV, Resolution, Bandwidth, seq, TR, TE, TI, TI2, alpha, noise_factor\
                                   ,T1_3D, T2_3D, M0_3D, B1_3D, flipangle_3D, t2_star_3D, ADC_3D, b, c, met, ETL, phi\
                                   , Readout_axis, post_traj, post_ETL, post_TE, B0_gene):
    """
    Function that will create two 3D simulation of low-field MRI sequences, one without noise and one with noise
    
    Inputs:
    
    //////// parameters to be chosen for runing the sequence ////////
    FOV          -> 1x3 array of the field of view
    Resolution   -> 1x3 array of the resolution
    Bandwidth    -> Bandwidth of aquisition
    seq          -> string defining the sequence to simulate; {'SE','GE','IN','Double IN','FLAIR','Dif','TSE','SSFP'} 
    The seq strings correspond to these sequences; {spin echo, gradient echo, inversion recovery, double inversion recovery, FLAIR, diffusion}
    TR           -> repetition time, must be in milli-seconds (ex: 3000 for 3 seconds, 50 for 50 milli-seconds)
    TE           -> echo time, must be in milli-seconds (ex: 160 for 0.160 seconds, 50 for 50 milli-seconds)
    TI           -> inversion time, must be in milli-seconds (ex: 650 for 0.650 seconds, 50 for 50 milli-seconds)
    TI2          -> second inversion time, must be in milli-seconds (ex: 150 for 0.150 seconds, 50 for 50 milli-seconds)
    alpha        -> flip angle, between {0-90}
    noise_factor -> multiplying noise factor
    b            -> b coefficient in diffusion sequence
    c            -> the number of points in the time signal generated in TSE seq (10 works well and is not to slow in terms of computation)
    met          -> is a string specifying the kspace trajectory in TSE; {'Out-in','In-out','Linear'}
    ETl          -> Echo train length
    phi          -> map, create from B0 map, representing of the difference in frequencies from center frequency
    Readout_axis -> readout axis of TSE or post TSE sequence
    post_traj    -> kspace trajectory of post TSE seq
    post_ETL     -> ETL of post TSE seq
    post_TE      -> TE of post TSE seq
    B0_gene      -> B0 strength of the generator; {'46','15','3'}

    //////// low-field systems and physiological maps ////////
    T1_3D        -> T1 relaxation map
    T2_3D        -> T2 relaxation map
    M0_3D        -> Proton density map
    B1_3D        -> B1 map
    flipangle_3D -> flipangle tensor map
    t2_star_3D   -> t2 star tensor map
    ADC_3D       -> apparent diffusion coefficient (ADC) tensor map
    
    Outputs:
    clean_data -> the simulated tensor sequence without noise
    noisy_data -> the simulated tensor sequence with noise
    """
    # Final data matrix size (field of view / reolution) 
    Data_mat = np.divide(FOV, Resolution)
    Data_mat = [int(Data_mat[0]), int(Data_mat[1]), int(Data_mat[2])]
    
    B0_const = 'B0_46mT'
    if B0_gene == '46':
        B0_const = 'B0_46mT'
        T1_3D = np.load(resource_path("Data\T1_3D_cer_lip_grad.npy"))
        T2_3D = np.load(resource_path("Data\T2_3D_cer_lip_grad.npy"))
        M0_3D = np.load(resource_path("Data\M0_3D_cer_lip_grad.npy"))
        B1_3D = np.load(resource_path("Data\B1_3D_cer_lip_grad.npy"))
        flipAngleMaprescale_3D = np.load(resource_path("Data\FlipAngleMaprescale_3D_cer_lip_grad.npy"))
        t2_star_3D = np.load(resource_path("Data\T2_star_tensor_3D_cer_lip_grad.npy"))
        ADC_3D = np.load(resource_path("Data\ADC_3D_cer_lip_grad.npy"))
    elif B0_gene == '15':
        B0_const = 'B0_15T'
        T1_3D = np.load(resource_path("Data\T1_15.npy"))
        T2_3D = np.load(resource_path("Data\T2_15.npy"))
        M0_3D = np.load(resource_path("Data\M0_3D_cerb_lip.npy"))
        B1_3D = np.load(resource_path("Data\B1_3D.npy"))
        flipAngleMaprescale_3D = np.load(resource_path("Data\FlipAngleMaprescale_3D.npy"))
        t2_star_3D = np.load(resource_path("Data\T2_star_15.npy"))
        ADC_3D = np.load(resource_path("Data\ADC_3D.npy"))
        B0_3D = b_field * 32.608 # 1.5T / 46mT ~= 32.608
        B0_3D[B0_3D>0] = 1.5
    elif B0_gene == '3':
        B0_const = 'B0_3T'
        T1_3D = np.load(resource_path("Data\T1_3.npy"))
        T2_3D = np.load(resource_path("Data\T2_3.npy"))
        M0_3D = np.load(resource_path("Data\M0_3D_cerb_lip.npy"))
        B1_3D = np.load(resource_path("Data\B1_3D.npy"))
        flipAngleMaprescale_3D = np.load(resource_path("Data\FlipAngleMaprescale_3D.npy"))
        t2_star_3D = np.load(resource_path("Data\T2_star_3.npy"))
        ADC_3D = np.load(resource_path("Data\ADC_3D.npy"))
        B0_3D = b_field * 65.217  # 3T / 46mT ~= 65.217
        B0_3D[B0_3D>0] = 3
        
    ## change the data to match the FOV specified by the user ##
    # The number or 2D arrays to delete from each axis (at the beginning and end, for each axis)
    to_delete = [int((T1_3D.shape[0] - FOV[0])/2), int((T1_3D.shape[1] - FOV[1])/2), int((T1_3D.shape[2] - FOV[2])/2)]
    # take a sub sample of the data corresponding to the FOV (from the center)
    sh = T1_3D.shape
    T1_3D = T1_3D[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
    T2_3D = T2_3D[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
    M0_3D = M0_3D[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
    B1_3D = B1_3D[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
    flipAngleMaprescale_3D = flipAngleMaprescale_3D[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
    t2_star_3D = t2_star_3D[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
    ADC_3D = ADC_3D[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
    
    # Computing the 3D sequence and resizing
    if seq == 'SE':
        Data_3D = spin_echo_seq(TR, TE, T1_3D, T2_3D, M0_3D)
        Data_3D = np.multiply(Data_3D, B1_3D)
        clean_data = resize(Data_3D, T1_3D, Data_mat)
    elif seq == 'GE':
        angle = flipangle_3D/alpha
        Data_3D =  Gradient_seq(TR, TE, T1_3D, t2_star_3D, M0_3D, angle)
        Data_3D = np.multiply(Data_3D, B1_3D)
        clean_data = resize(Data_3D, T1_3D, Data_mat)
    elif seq == 'IN':
        Data_3D = IN_seq(TR, TI, T1_3D, T2_3D, M0_3D) 
        Data_3D = np.multiply(Data_3D, B1_3D)
        Data_3D = post_TSE_TEeff(Data_3D, post_TE, post_ETL, T2_3D, post_traj)
        clean_data = resize(Data_3D, T1_3D, Data_mat)
        clean_data = TSE_seq_R(post_TE, post_ETL, T2_3D, post_traj, Readout_axis, clean_data)
    elif seq == 'Double IN':
        Data_3D = DoubleInversion_seq(TR, post_TE, TI, TI2, T1_3D, T2_3D, M0_3D)
        Data_3D = np.multiply(Data_3D, B1_3D)
        Data_3D = post_TSE_TEeff(Data_3D, post_TE, post_ETL, T2_3D, post_traj)
        clean_data = resize(Data_3D, T1_3D, Data_mat)        
        clean_data = TSE_seq_R(post_TE, post_ETL, T2_3D, post_traj, Readout_axis, clean_data)
    elif seq == 'FLAIR':
        #TI = np.log(2) * 3.695 
        Data_3D = IN_seq(TR, TI, T1_3D, T2_3D, M0_3D)
        Data_3D = np.multiply(Data_3D, B1_3D)   
        Data_3D = post_TSE_TEeff(Data_3D, post_TE, post_ETL, T2_3D, post_traj)
        clean_data = resize(Data_3D, T1_3D, Data_mat)
        clean_data = TSE_seq_R(post_TE, post_ETL, T2_3D, post_traj, Readout_axis, clean_data)       
    elif seq == 'Dif':
        Data_3D = Diffusion_seq(TR, T1_3D, T2_3D, M0_3D, b, ADC_3D)
        Data_3D = np.multiply(Data_3D, B1_3D)    
        Data_3D = post_TSE_TEeff(Data_3D, post_TE, post_ETL, T2_3D, post_traj)
        clean_data = resize(Data_3D, T1_3D, Data_mat)
        clean_data = TSE_seq_R(post_TE, post_ETL, T2_3D, post_traj, Readout_axis, clean_data)       
    elif seq == 'TSE':
        TSE_3D = 0
        TEeff = 0
        if post_traj == 'Out-in':    
            TEeff = TE * ETL    
            Data_3D = spin_echo_seq(TR, TEeff, T1_3D_grad, T2_3D_grad, M0_3D_grad) # Computing a spin echo
        elif post_traj == 'In-out':
            TEeff = TE
            TE_div = np.divide(TEeff, T2_3D_grad, out=np.zeros_like(T2_3D_grad), where=T2_3D_grad!=0)
            Data_3D = np.abs(M0_3D_grad * np.exp(-TE_div))                         # Computing a spin echo with the (1 - exp(-TR/T1)) set to 1
        elif post_traj == 'Linear':
            TEeff = 0.5 * TE * ETL
            Data_3D = spin_echo_seq(TR, TEeff, T1_3D_grad, T2_3D_grad, M0_3D_grad) # Computing a spin echo

        Data_3D = np.multiply(Data_3D, B1_3D)
        clean_data = resize(Data_3D, T1_3D, Data_mat)        
        clean_data = TSE_seq_R(post_TE, post_ETL, T2_3D, post_traj, Readout_axis, clean_data)  
    elif seq == 'SSFP':
        Data_3D = SSFP_new(TR, T1_3D, T2_3D, M0_3D, alpha, phi)
        Data_3D = np.multiply(Data_3D, B1_3D)
        clean_data = resize(Data_3D, T1_3D, Data_mat)
        
    ##### NOISE #####
    noise = noise_generation(Data_mat, Resolution, Bandwidth, B0_const)
    noise = noise*noise_factor      # constant that multiply the amount of noise to add
    noisy_data = clean_data + noise # Adding the noise
    
    return np.abs(clean_data), np.abs(noisy_data)    

#/////// Functions regarding the generation/visualisation/saving of a single dataset defined by the user /////////# 
def create_dataset():
    """
    Generates a single dataset comprised of two 3D matrices (low and high SNR)
    Images will be ploted on the GUI
    Input: None
    Output: None
    """
    global Update_dataset_creation
    global seq_name
    global clean_data
    global noise_data
    global FOV_gene
    global Resolution_gene
    global bd_gene
    global Pre_def_seq_gene
    global TR_gene
    global TE_gene 
    global TI_gene
    global TI2_gene
    global Alpha_gene
    global Generator_noise_factor
    global T1_3D_grad
    global T2_3D_grad
    global M0_3D_grad
    global B1map_3D_grad
    global flipAngleMaprescale_3D_grad
    global t2_star_3D_grad 
    global ADC_3D_grad
    global c_gene
    global trajectory_gene
    global ETL_gene
    global phi_gene
    global B0_gene
    global canvas1
    global canvas2
    global canvas3
    
    canvas1.get_tk_widget().destroy()
    canvas2.get_tk_widget().destroy()
    canvas3.get_tk_widget().destroy()
    
    TR_gene = np.divide(int(TR_entry_gene.get()),1000) # Divide by 1000 to have the values in milli
    TE_gene = np.divide(int(TE_entry_gene.get()),1000)
    TI_gene = np.divide(int(TI_entry_gene.get()),1000)
    TI2_gene = np.divide(int(TI2_entry_gene.get()),1000)
    FOV_gene = [int(FOV1_entry_gene.get()), int(FOV2_entry_gene.get()), int(FOV3_entry_gene.get())]
    Resolution_gene = [float(Res1_entry_gene.get()), float(Res2_entry_gene.get()), float(Res3_entry_gene.get())]
    Data_mat_gene = [int(Data_mat1_entry_gene.get()), int(Data_mat2_entry_gene.get()), int(Data_mat3_entry_gene.get())]
    bd_gene = int(Bandwidth_entry_gene.get())
    Alpha_gene = int(Alpha_entry.get())
    ETL_gene = int(ETL_entry.get()) 
    Generator_noise_factor = float(noise_factor_entry_gene.get())
    c_gene = 10
    trajectory_gene = traj_gene.get()
    Read = Readout_axis_gene
    post_traj = post_traj_gene.get()
    post_ETL = int(post_TSE_ETL_entry_gene.get())
    post_TE = np.divide(int(post_TSE_TE_entry_gene.get()),1000)
    
    phi_gene = 1
    if Pre_def_seq_gene == 'SSFP':
        # Phi is specific to SSFP seq
        gamma = 42.58*(10**6)*2*constants.pi # 267538030.37970677 [rad/sT]
        omega = B0_3D * gamma
        center = np.divide(omega.shape,2).astype(int)
        center_freq_value = omega[center[0],center[1],center[2]]
        offset = omega - center_freq_value
        phi_gene = offset * np.divide(int(TR_entry_gene.get()),1000) 
        
    # B and minTE for diffusion sequence
    b_gene = 0
    if Pre_def_seq_gene == 'Dif':
        # Diffusion parameters
        global Bval_label_gene
        global TEminval_label_gene
        gamma =  42.58*(10**6)#*2*constants.pi      # gyromagnetic ratio for hydrogen 42.58 [MHz/T] 

        G = float(G_entry_gene.get()); G = np.divide(G,1000)
        smalldelta = float(smalldelta_entry_gene.get()); smalldelta = np.divide(smalldelta,1000)
        bigdelta = float(bigdelta_entry_gene.get()); bigdelta = np.divide(bigdelta,1000)

        b_gene = (bigdelta - smalldelta/3)*(gamma*G*smalldelta)**2; 
        dif_TEmin = smalldelta + bigdelta; 

        Bval_label_gene.grid_forget()
        Bval_label_gene = Label(frame3, text = np.round(B,2), font=("Helvetica", 12)); Bval_label_gene.grid(row = 7, column = 4)
        TEminval_label_gene.grid_forget()
        TEminval_label_gene = Label(frame3, text = np.round(dif_TEmin*1000,2), font=("Helvetica", 12)); TEminval_label_gene.grid(row = 8, column = 4)
    
    ## change the data to match the FOV specified by the user ##
    # The number or 2D arrays to delete from each axis (at the beginning and end, for each axis)
    to_delete = [int((T1_3D_grad.shape[0] - FOV_gene[0])/2), int((T1_3D_grad.shape[1] - FOV_gene[1])/2), int((T1_3D_grad.shape[2] - FOV_gene[2])/2)]
    # take a sub sample of the data corresponding to the FOV (from the center)
    sh = T1_3D_grad.shape
    T1_3D_grad = T1_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
    T2_3D_grad = T2_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
    M0_3D_grad = M0_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
    B1map_3D_grad = B1map_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
    flipAngleMaprescale_3D_grad = flipAngleMaprescale_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
    t2_star_3D_grad = t2_star_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
    ADC_3D_grad = ADC_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]

    # Calling the generation function
    clean_data, noise_data = create_3D_low_and_high_SNR_data(FOV_gene, Resolution_gene, bd_gene, Pre_def_seq_gene, TR_gene, TE_gene,\
                                                            TI_gene, TI2_gene, Alpha_gene, Generator_noise_factor, T1_3D_grad, T2_3D_grad,\
                                                            M0_3D_grad, B1map_3D_grad, flipAngleMaprescale_3D_grad, t2_star_3D_grad,\
                                                            ADC_3D_grad, b_gene, c_gene, trajectory_gene, ETL_gene, phi_gene, Read, post_traj, post_ETL, post_TE, B0_gene)
    
    clean1 = np.flip(np.rot90(clean_data[int(Data_mat_gene[0]/2),:,:]), axis = 1)
    clean2 = np.rot90(clean_data[:,int(Data_mat_gene[1]/2),:])
    clean3 = np.rot90(clean_data[:,:,int(Data_mat_gene[2]/2)])
    noise1 = np.flip(np.rot90(noise_data[int(Data_mat_gene[0]/2),:,:]), axis = 1)
    noise2 = np.rot90(noise_data[:,int(Data_mat_gene[1]/2),:])
    noise3 = np.rot90(noise_data[:,:,int(Data_mat_gene[2]/2)])
    
    imag_size_width = 5.5
    imag_size_height = 5.5
    fig = plt.figure(figsize=(imag_size_width,imag_size_height))  
    plt.subplot(2, 1, 1)
    plt.imshow(clean1, cmap = 'gray')           
    plt.axis('off')
    plt.subplot(2, 1, 2)
    plt.imshow(noise1, cmap = 'gray')           
    plt.axis('off')
    canvas1 = FigureCanvasTkAgg(fig, root)                          
    canvas1.draw()
    canvas1.get_tk_widget().grid(row = 1, column = 1, rowspan = 4)        
    plt.close()

    fig = plt.figure(figsize=(imag_size_width,imag_size_height))
    plt.subplot(2, 1, 1)
    plt.imshow(clean2, cmap = 'gray')           
    plt.axis('off')
    plt.subplot(2, 1, 2)
    plt.imshow(noise2, cmap = 'gray')           
    plt.axis('off')
    canvas2 = FigureCanvasTkAgg(fig, root)         
    canvas2.draw()
    canvas2.get_tk_widget().grid(row = 1, column = 2, rowspan = 4)        
    plt.close()

    fig = plt.figure(figsize=(imag_size_width,imag_size_height))
    plt.subplot(2, 1, 1)
    plt.imshow(clean3, cmap = 'gray')           
    plt.axis('off')
    plt.subplot(2, 1, 2)
    plt.imshow(noise3, cmap = 'gray')           
    plt.axis('off')
    canvas3 = FigureCanvasTkAgg(fig, root)         
    canvas3.draw()
    canvas3.get_tk_widget().grid(row = 1, column = 3, rowspan = 4)          
    plt.close()
    
    Update_dataset_creation.grid_forget()
    d1 =("High and low SNR " + str(seq_name) + " data sets created \n"
    "Both have shape " + str(Data_mat_gene))
    Update_dataset_creation = Label(root, text = d1, font=("Helvetica", 15))
    Update_dataset_creation.grid(row = 5, column = 3, columnspan = 2, rowspan = 2)
    
    Sagittal_label_gene.grid(row = 0, column = 1)
    Coronal_label_gene.grid(row = 0, column = 2)
    Axial_label_gene.grid(row = 0, column = 3)    
    
def visu_dataset():
    """
    Opens an interactive window to visualise the generated dataset
    Input: None
    Output: None
    """
    global clean_data
    global noise_data
    FOV = [int(FOV1_entry_gene.get()), int(FOV2_entry_gene.get()), int(FOV3_entry_gene.get())]
    Data_mat = clean_data.shape
    fig, ax = plt.subplots(1,2)
    fig3D = pltviewdatasets.interactivePlot(fig, ax, clean_data, noise_data, Data_mat, plotAxis = 2, fov = FOV)
    plt.show()
    
def save_dataset():
    """
    Save the generated dataset as a numpy array
    Input: None
    Output: None
    """
    global clean_data
    global noise_data
    global seq_name
    global FOV_gene
    global Resolution_gene
    global bd_gene
    global Pre_def_seq_gene
    global TR_gene
    global TE_gene 
    global TI_gene
    global TI2_gene
    global Alpha_gene
    global Generator_noise_factor
    global T1_3D_grad
    global T2_3D_grad
    global M0_3D_grad
    global B1map_3D_grad
    global flipAngleMaprescale_3D_grad
    global t2_star_3D_grad 
    global ADC_3D_grad
    global c_gene
    global trajectory_gene
    global ETL_gene
    global phi_gene
    
    post_traj = post_traj_gene.get()
    post_ETL = int(post_TSE_ETL_entry_gene.get())
    post_TE = np.divide(int(post_TSE_TE_entry_gene.get()),1000)
    date_now = datetime.datetime.now().strftime("%Y-%m-%d")
    time_now = datetime.datetime.now().strftime("%H-%M-%S")

    # Will return the directory of the file
    filename_filterdata = filedialog.askdirectory()
    file_path_label = str(filename_filterdata)
    t = 'Clean {} #date {} #time {}'.format(seq_name, date_now, time_now)
    np.save(os.path.join(file_path_label, t), clean_data)
    t = 'Noisy {} #date {} #time {}'.format(seq_name, date_now, time_now)
    np.save(os.path.join(file_path_label, t), noise_data)
    
    # creating (if a file with this name doesn't exists yet) and saving a txt file with a list of all the parameters
    lines = ['Parameters:', 'FOV = {} mm3'.format(FOV_gene), 'Resolution = {} mm3'.format(Resolution_gene), \
             'Bandwidth = {} Hz'.format(bd_gene), 'Sequence = {}'.format('TSE'), 'TR = {} s'.format(TR_gene), \
             'TE = {} s'.format(TE_gene), 'TI = {} s'.format(TI_gene), 'TI2 = {} s'.format(TI2_gene),\
             'Alpha = {}°'.format(Alpha_gene), 'kspace trajectory (for TSE) = {}'.format(trajectory_gene), \
             'ETL = {}'.format(ETL_gene), 'Noise factor = {}'.format(Generator_noise_factor), \
             'Post TSE seq traj = {}'.format(post_traj), 'Post TSE seq ETL = {}'.format(post_ETL), 'Post TSE seq TE = {}'.format(post_TE), \
             'Post TSE seq readout axis = {}'.format(Readout_axis_gene)]
    path = (file_path_label + '/Parameters #{} #{} #{}.txt'.format(seq_name, date_now, time_now))

    with open(path, 'w') as f:
        for line in lines:
            f.write(line)
            f.write('\n')
            
def save_dataset_nifti():
    """
    Save the generated dataset as a nifti file
    Input: None
    Output: None
    """
    post_traj = post_traj_gene.get()
    post_ETL = int(post_TSE_ETL_entry_gene.get())
    post_TE = np.divide(int(post_TSE_TE_entry_gene.get()),1000)
    date_now = datetime.datetime.now().strftime("%Y-%m-%d")
    time_now = datetime.datetime.now().strftime("%H-%M-%S")

    # Will return the directory of the file
    filename_filterdata = filedialog.askdirectory()
    file_path_label = str(filename_filterdata)
    t = 'NIFTI Clean {} #date {} #time {}'.format(seq_name, date_now, time_now)
    new_image = nib.Nifti1Image(clean_data, affine=np.eye(4))
    new_image.to_filename(os.path.join(file_path_label, t))
    t = 'NIFTI Noisy {} #date {} #time {}'.format(seq_name, date_now, time_now)
    new_image = nib.Nifti1Image(noise_data, affine=np.eye(4))
    new_image.to_filename(os.path.join(file_path_label, t))
    
    # creating (is a file with this name doesn't exists yet) and saving a txt file with a list of all the parameters
    lines = ['Parameters:', 'FOV = {} mm3'.format(FOV_gene), 'Resolution = {} mm3'.format(Resolution_gene), \
             'Bandwidth = {} Hz'.format(bd_gene), 'Sequence = {}'.format('TSE'), 'TR = {} s'.format(TR_gene), \
             'TE = {} s'.format(TE_gene), 'TI = {} s'.format(TI_gene), 'TI2 = {} s'.format(TI2_gene),\
             'Alpha = {}°'.format(Alpha_gene), 'kspace trajectory (for TSE) = {}'.format(trajectory_gene), \
             'ETL = {}'.format(ETL_gene), 'Noise factor = {}'.format(Generator_noise_factor), \
             'Post TSE seq traj = {}'.format(post_traj), 'Post TSE seq ETL = {}'.format(post_ETL), 'Post TSE seq TE = {}'.format(post_TE), \
             'Post TSE seq readout axis = {}'.format(Readout_axis_gene)]
    path = (file_path_label + '/NIFTI Parameters #{} #{} #{}.txt'.format(seq_name, date_now, time_now))

    with open(path, 'w') as f:
        for line in lines:
            f.write(line)
            f.write('\n')

#/////// Functions regarding the generation/visualisation/saving of multiple datasets via a dictonnary of sequences defined by the user /////////#      
def add_to_dict():
    """
    Adds a new sequence with all it's parameters to the dictionnary 
    Input: None
    Output: None 
    """
    global dict_multiple_dataset
    global FOV_gene
    global Resolution_gene
    global bd_gene
    global Pre_def_seq_gene
    global TR_gene
    global TE_gene 
    global TI_gene
    global TI2_gene
    global Alpha_gene
    global Generator_noise_factor
    global trajectory_gene
    global ETL_gene
    global phi_gene
    global list_gene
    global dict_drop_gene
    global seq_to_visu
    global drop_to_visu_gene
    
    FOV_gene = [int(FOV1_entry_gene.get()), int(FOV2_entry_gene.get()), int(FOV3_entry_gene.get())]
    Resolution_gene = [float(Res1_entry_gene.get()), float(Res2_entry_gene.get()), float(Res3_entry_gene.get())]
    bd_gene = int(Bandwidth_entry_gene.get())
    Data_mat_gene = [int(Data_mat1_entry_gene.get()), int(Data_mat2_entry_gene.get()), int(Data_mat3_entry_gene.get())]
    TR_gene = np.divide(int(TR_entry_gene.get()),1000) # Divide by 1000 to have the values in milli
    TE_gene = np.divide(int(TE_entry_gene.get()),1000)
    TI_gene = np.divide(int(TI_entry_gene.get()),1000)
    TI2_gene = np.divide(int(TI2_entry_gene.get()),1000)
    Alpha_gene = int(Alpha_entry.get())
    ETL_gene = int(ETL_entry.get())     
    trajectory_gene = traj_gene.get()
    Generator_noise_factor = float(noise_factor_entry_gene.get())
    post_traj = post_traj_gene.get()
    post_ETL = int(post_TSE_ETL_entry_gene.get())
    post_TE = np.divide(int(post_TSE_TE_entry_gene.get()),1000)
    
    l = len(dict_multiple_dataset)
    t = 'Seq{}'.format(l+1)
    
    # adding a new element to the dictonary
    dict_multiple_dataset[t] = [Pre_def_seq_gene, FOV_gene, Resolution_gene, bd_gene, Data_mat_gene, TR_gene, TE_gene, TI_gene, TI2_gene, Alpha_gene, ETL_gene, trajectory_gene, Generator_noise_factor, post_traj, post_ETL, post_TE, Readout_axis_gene]

    # adding new element to list (form visualization in dropdown menu)
    s = '{}, FOV:{}, Resolution:{}, Bandwidth:{}, final data size:{}, TR:{}, TE:{}, TI:{}, TI2:{}, Alpha:{}, ETL:{}, kspace traj:{}, Noise factor:{}, Post TSE seq traj:{}, Post TSE seq ETL:{}, Post TSE seq TE:{}, Post TSE seq readout axis:{}'.format(Pre_def_seq_gene, FOV_gene, Resolution_gene, bd_gene, Data_mat_gene, TR_gene*1000, TE_gene*1000, TI_gene*1000, TI2_gene*1000, Alpha_gene, ETL_gene, trajectory_gene, Generator_noise_factor, post_traj, post_ETL, post_TE*1000, Readout_axis_gene)
    list_gene.append(s)
    
    dict_drop_gene.grid_forget()
    seq_list_gene = StringVar()
    seq_list_gene.set("List of sequences")                         # Default value, or could use; options[0]
    dict_drop_gene = OptionMenu(frame9, seq_list_gene, *list_gene) # Using a list, NEEDS a star in front
    dict_drop_gene.grid(row = 1, column = 1)
    
    seq_to_visu.append('{}'.format(t))
    drop_to_visu_gene.grid_forget()
    drop_to_visu_gene = OptionMenu(frame9, list_to_visu_gene, *seq_to_visu) 
    drop_to_visu_gene.grid(row = 3, column = 1)

def view_dict():    
    """
    Print the sequences that are in the dictionary, thise are the sequences that will be generated when the creat buttons is pressed
    Input: None
    Output: None (print the sequences on terminal)
    """
    for i in range(len(list_gene)): # the + 1 is because the first element of the list is the string "List of sequences"
        if i > 0:
            print(list_gene[i])
    
def create_datasets():
    """
    Generates all the sequences in the dictionnary, each one is comprised of two 3D matrices (low and high SNR)
    Images will not be ploted on the GUI
    Input: None
    Output: None 
    """
    global T1_3D_grad
    global T2_3D_grad
    global M0_3D_grad
    global B1map_3D_grad
    global flipAngleMaprescale_3D_grad
    global t2_star_3D_grad 
    global ADC_3D_grad
    global clean_dict
    global noise_dict
    global Update_datasets_creation
    global num_seq_created_label
    global B0_gene
    
    c_gene = 10
    
    seq = [] # list that will contain all the strings 'seq0', 'seq1',... which are the keys of the dictonary holding the sequences to generate
    num_seq = len(dict_multiple_dataset) # number of sequences to generate
    
    for j in range(num_seq):
        seq.append('Seq{}'.format(j+1))
    
    for i in range(len(seq)):
        l = dict_multiple_dataset[seq[i]]
        phi_gene = 1
        if l[0] == 'SSFP':
            # Phi is specific to SSFP seq
            gamma = 42.58*(10**6)*2*constants.pi # 267538030.37970677 [rad/sT]
            omega = B0_3D * gamma
            center = np.divide(omega.shape,2).astype(int)
            center_freq_value = omega[center[0],center[1],center[2]]
            offset = omega - center_freq_value
            phi_gene = offset * np.divide(int(TR_entry_gene.get()),1000) 
        
        b = 1

        # B and minTE for diffusion sequence
        b_gene = 0
        if Pre_def_seq_gene == 'Dif':
            # Diffusion parameters
            global Bval_label_gene
            global TEminval_label_gene
            gamma =  42.58*(10**6)#*2*constants.pi      # gyromagnetic ratio for hydrogen 42.58 [MHz/T] 

            G = float(G_entry_gene.get()); G = np.divide(G,1000)
            smalldelta = float(smalldelta_entry_gene.get()); smalldelta = np.divide(smalldelta,1000)
            bigdelta = float(bigdelta_entry_gene.get()); bigdelta = np.divide(bigdelta,1000)

            b_gene = (bigdelta - smalldelta/3)*(gamma*G*smalldelta)**2; 
            dif_TEmin = smalldelta + bigdelta; 

            Bval_label_gene.grid_forget()
            Bval_label_gene = Label(frame3, text = np.round(B,2), font=("Helvetica", 12));              Bval_label_gene.grid(row = 7, column = 4)
            TEminval_label_gene.grid_forget()
            TEminval_label_gene = Label(frame3, text = np.round(dif_TEmin*1000,2), font=("Helvetica", 12)); TEminval_label_gene.grid(row = 8, column = 4)
            
        ## change the data to match the FOV specified by the user ##
        # The number or 2D arrays to delete from each axis (at the beginning and end, for each axis)
        FieldOfView = l[1]
        to_delete = [int((T1_3D_grad.shape[0] - FieldOfView[0])/2), int((T1_3D_grad.shape[1] - FieldOfView[1])/2), int((T1_3D_grad.shape[2] - FieldOfView[2])/2)]
        # take a sub sample of the data corresponding to the FOV (from the center)
        sh = T1_3D_grad.shape
        T1_3D_grad = T1_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
        T2_3D_grad = T2_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
        M0_3D_grad = M0_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
        B1map_3D_grad = B1map_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
        flipAngleMaprescale_3D_grad = flipAngleMaprescale_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
        t2_star_3D_grad = t2_star_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
        ADC_3D_grad = ADC_3D_grad[to_delete[0]:(sh[0]-to_delete[0]), to_delete[1]:(sh[1]-to_delete[1]), to_delete[2]:(sh[2]-to_delete[2])]
        
        # Calling the generation function
        clean_data, noise_data = create_3D_low_and_high_SNR_data(l[1], l[2], l[3], l[0], l[5], l[6],\
                                                                l[7], l[8], l[9], l[12], T1_3D_grad, T2_3D_grad,\
                                                                M0_3D_grad, B1map_3D_grad, flipAngleMaprescale_3D_grad, t2_star_3D_grad,\
                                                                ADC_3D_grad, b, c_gene, l[11], l[10], phi_gene, l[16], l[13], l[14], l[15], B0_gene)    
        clean_dict[seq[i]] = clean_data
        noise_dict[seq[i]] = noise_data
        
        y = (seq[i] + ' done on ' + str(len(seq)))
        num_seq_created_label.grid_forget()
        num_seq_created_label = Label(frame9, text = y, font=("Helvetica", 15))
        num_seq_created_label.grid(row = 2, column = 1)
        print(y)
    
    print('All sequences have been generated')
    
    Update_datasets_creation.grid_forget()
    d1 =("All high and low SNR data sets have been created")
    Update_datasets_creation = Label(root, text = d1, font=("Helvetica", 15))
    Update_datasets_creation.grid(row = 7, column = 3, columnspan = 2, rowspan = 2)  
    
def visu_create_datasets():
    """
    Opens an interactive window to visualize one pair of the sequences created
    Input: None
    Output: None 
    """
    global drop_to_visu_gene
    global clean_dict
    global noise_dict
    
    selected = list_to_visu_gene.get()    
    fig, ax = plt.subplots(1,2)
    fig3D = pltviewdatasets.interactivePlot(fig, ax, clean_dict[selected], noise_dict[selected], clean_dict[selected].shape, plotAxis = 2, fov = FOV)
    plt.show()    
      
def save_create_datasets():
    """
    Saves the multiple datasets created as numpy arrays
    Input: None
    Output: None 
    """
    global seq_to_visu
    global clean_dict
    global noise_dict
    global dict_multiple_dataset
    
    date_now = datetime.datetime.now().strftime("%Y-%m-%d")
    time_now = datetime.datetime.now().strftime("%H-%M-%S")

    # Will return the directory of the file
    filename_filterdata = filedialog.askdirectory()
    file_path_label = str(filename_filterdata)
    
    for i in range(1,len(seq_to_visu)):
        key = seq_to_visu[i]
        t = 'Clean {} #date {} #time {}'.format(seq_to_visu[i], date_now, time_now)
        np.save(os.path.join(file_path_label, t), clean_dict[key])
        t = 'Noisy {} #date {} #time {}'.format(seq_to_visu[i], date_now, time_now)
        np.save(os.path.join(file_path_label, t), noise_dict[key])

        l = dict_multiple_dataset[seq_to_visu[i]]

        # creating (is a file with this name doesn't exists yet) and saving a txt file with a list of all the parameters
        lines = ['Parameters:', 'FOV = {} mm3'.format(l[1]), 'Resolution = {} mm3'.format(l[2]), \
                 'Bandwidth = {} Hz'.format(l[3]), 'Sequence = {}'.format('TSE'), 'TR = {} s'.format(l[5]), \
                 'TE = {} s'.format(l[6]), 'TI = {} s'.format(l[7]), 'TI2 = {} s'.format(l[8]),\
                 'Alpha = {}°'.format(l[9]), 'kspace trajectory (for TSE) = {}'.format(l[11]), \
                 'ETL = {}'.format(l[10]), 'Noise factor = {}'.format(l[12]), 'Post TSE seq traj = {}'.format(l[13]),\
                 'Post TSE seq ETL = {}'.format(l[14]), 'Post TSE seq TE = {}'.format(l[15]), 'Post TSE seq readout axis = {}'.format(l[16])]
        path = (file_path_label + '/ Parameters {}, {} #{} #{}.txt'.format(seq_to_visu[i], l[0], date_now, time_now))

        with open(path, 'w') as f:
            for line in lines:
                f.write(line)
                f.write('\n')
                
def save_create_datasets_nifti():
    """
    Saves the multiple datasets created as nifti files
    Input: None
    Output: None 
    """
    global seq_to_visu
    global clean_dict
    global noise_dict
    global dict_multiple_dataset
    
    date_now = datetime.datetime.now().strftime("%Y-%m-%d")
    time_now = datetime.datetime.now().strftime("%H-%M-%S")

    # Will return the directory of the file
    filename_filterdata = filedialog.askdirectory()
    file_path_label = str(filename_filterdata)
    
    for i in range(1,len(seq_to_visu)):
        key = seq_to_visu[i]
        t = 'NIFTI Clean {} #date {} #time {}'.format(seq_to_visu[i], date_now, time_now)
        new_image = nib.Nifti1Image(clean_dict[key], affine=np.eye(4))
        new_image.to_filename(os.path.join(file_path_label, t))
        t = 'NIFTI Noisy {} #date {} #time {}'.format(seq_to_visu[i], date_now, time_now)
        new_image = nib.Nifti1Image(noise_dict[key], affine=np.eye(4))
        new_image.to_filename(os.path.join(file_path_label, t))
        
        l = dict_multiple_dataset[seq_to_visu[i]]

        # creating (is a file with this name doesn't exists yet) and saving a txt file with a list of all the parameters
        lines = ['Parameters:', 'FOV = {} mm3'.format(l[1]), 'Resolution = {} mm3'.format(l[2]), \
                 'Bandwidth = {} Hz'.format(l[3]), 'Sequence = {}'.format('TSE'), 'TR = {} s'.format(l[5]), \
                 'TE = {} s'.format(l[6]), 'TI = {} s'.format(l[7]), 'TI2 = {} s'.format(l[8]),\
                 'Alpha = {}°'.format(l[9]), 'kspace trajectory (for TSE) = {}'.format(l[11]), \
                 'ETL = {}'.format(l[10]), 'Noise factor = {}'.format(l[12]), 'Post TSE seq traj = {}'.format(l[13]),\
                 'Post TSE seq ETL = {}'.format(l[14]), 'Post TSE seq TE = {}'.format(l[15]), 'Post TSE seq readout axis = {}'.format(l[16])]
        path = (file_path_label + '/NIFTI parameters {}, {} #{} #{}.txt'.format(seq_to_visu[i], l[0], date_now, time_now))

        with open(path, 'w') as f:
            for line in lines:
                f.write(line)
                f.write('\n')
                
#/////// Functions reseting generator labels/images/dictionnary/parameters for a generating a new dataset(s) /////////#    
def new_dataset_gene():
    """
    Resets generator images/dictionnary/labels
    Input: None
    Output: None 
    """
    global dict_multiple_dataset
    global Update_dataset_creation
    global Update_datasets_creation
    global list_gene
    global dict_drop_gene
    global num_seq_created_label
    global list_to_visu_gene
    global drop_to_visu_gene
    global seq_to_visu
    global canvas1
    global canvas2
    global canvas3
    
    canvas1.get_tk_widget().destroy()
    canvas2.get_tk_widget().destroy()
    canvas3.get_tk_widget().destroy()
    
    dict_multiple_dataset = dict()
    
    dict_drop_gene.grid_forget()
    list_gene = ['Sequences that will be generated']
    seq_list_gene = StringVar()
    seq_list_gene.set("List of sequences")                         # Default value, or could use; options[0]
    dict_drop_gene = OptionMenu(frame9, seq_list_gene, *list_gene) # Using a list, NEEDS a star in front
    dict_drop_gene.grid(row = 1, column = 1)
    
    seq_to_visu = ['Choose one from list bellow']
    list_to_visu_gene = StringVar()
    list_to_visu_gene.set("Sequence to visualize")                         
    drop_to_visu_gene = OptionMenu(frame9, list_to_visu_gene, *seq_to_visu) 
    drop_to_visu_gene.grid(row = 3, column = 1)
    
    Update_dataset_creation.grid_forget()
    Update_dataset_creation = Label(root, text = "Data set not yet created", font=("Helvetica", 15))
    Update_dataset_creation.grid(row = 5, column = 3, columnspan = 2)
    
    Update_datasets_creation.grid_forget()
    Update_datasets_creation = Label(root, text = "Data sets not yet created", font=("Helvetica", 15))
    Update_datasets_creation.grid(row = 7, column = 3, columnspan = 2, rowspan = 2)
    
    num_seq_created_label.grid_forget()
    num_seq_created_label = Label(frame9, text = 'No sequences generated', font=("Helvetica", 15))
    num_seq_created_label.grid(row = 2, column = 1)
    
    # Read Images
    img1_gene = Image.open(resource_path("Images\Generator1.png"))
    img2_gene = Image.open(resource_path("Images\Generator2.png"))
    img3_gene = Image.open(resource_path("Images\Generator3.png"))
    
    imag_size = 5.5
    fig = plt.figure(figsize=(imag_size,imag_size))                
    plt.imshow(img1_gene)           
    plt.axis('off')
    canvas1 = FigureCanvasTkAgg(fig, root)                          
    canvas1.draw()
    canvas1.get_tk_widget().grid(row = 1, column = 1, rowspan = 4)        
    plt.close()

    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(img2_gene)   
    plt.axis('off')
    canvas2 = FigureCanvasTkAgg(fig, root)         
    canvas2.draw()
    canvas2.get_tk_widget().grid(row = 1, column = 2, rowspan = 4)        
    plt.close()

    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(img3_gene)   
    plt.axis('off')
    canvas3 = FigureCanvasTkAgg(fig, root)         
    canvas3.draw()
    canvas3.get_tk_widget().grid(row = 1, column = 3, rowspan = 4)          
    plt.close()    
    
def reset_param_gene():
    """
    Resets generator parameters
    Input: None
    Output: None 
    """
    global Bval_label_gene
    global TEminval_label_gene
    global FOV_check_gene
    global Res_check_gene 
    global Data_matrix_check_gene
    
    B0_gene = '46'
    Field_strength_gene.set("B0 field strength 46 mT")  
    FOV1_entry_gene.delete(0, END); FOV1_entry_gene.insert(0, '250')
    FOV2_entry_gene.delete(0, END); FOV2_entry_gene.insert(0, '300')
    FOV3_entry_gene.delete(0, END); FOV3_entry_gene.insert(0, '275')
    Res1_entry_gene.delete(0, END); Res1_entry_gene.insert(0, '1.95')
    Res2_entry_gene.delete(0, END); Res2_entry_gene.insert(0, '2.34')
    Res3_entry_gene.delete(0, END); Res3_entry_gene.insert(0, '2.15')
    Data_mat1_entry_gene.delete(0, END); Data_mat1_entry_gene.insert(0, '128')
    Data_mat2_entry_gene.delete(0, END); Data_mat2_entry_gene.insert(0, '128')
    Data_mat3_entry_gene.delete(0, END); Data_mat3_entry_gene.insert(0, '128')
    Bandwidth_entry_gene.delete(0,END); Bandwidth_entry_gene.insert(0, '40000')
    TR_entry_gene.delete(0, END); TR_entry_gene.insert(0, '500')
    TE_entry_gene.delete(0, END); TE_entry_gene.insert(0, '20')
    TI_entry_gene.delete(0, END); TI_entry_gene.insert(0, '250')
    TI2_entry_gene.delete(0, END); TI2_entry_gene.insert(0, '250')
    Alpha_entry_gene.delete(0, END); Alpha_entry_gene.insert(0, '45')
    ETL_entry_gene.delete(0, END); ETL_entry_gene.insert(0, '8')
    G_entry_gene.delete(0, END); G_entry_gene.insert(0, '10')
    smalldelta_entry_gene.delete(0, END); smalldelta_entry_gene.insert(0, '1')
    bigdelta_entry_gene.delete(0, END); bigdelta_entry_gene.insert(0, '2')
    noise_factor_entry_gene.delete(0, END); noise_factor_entry_gene.insert(0, '1')
    traj_gene.set("Linear")   
    Readout_gene.set("Along Sagittal Foot/head")
    post_traj_gene.set("Linear")
    Post_Readout_gene.set("Along Sagittal Foot/head")
    post_TSE_TE_entry_gene.delete(0, END); post_TSE_TE_entry_gene.insert(0, '20')
    post_TSE_ETL_entry_gene.delete(0, END); post_TSE_ETL_entry_gene.insert(0, '8')
    Bval_label_gene.grid_forget()
    Bval_label_gene = Label(frame7, text = "302.18", font=("Helvetica", f))
    TEminval_label_gene.grid_forget()
    TEminval_label_gene = Label(frame7, text = "3.0", font=("Helvetica", f))
    FOV_check_gene = 20625000        # 250*300*275 will be used to see if the FOV was changed by the user
    Res_check_gene = 9.81045         # 1.95*2.34*2.15 will be used to see if the resolution was changed by the user
    Data_matrix_check_gene = 2097152 # 128*128*128 will be used to see if the Data matrix size was changed by the user
    
    if Pre_def_seq_gene == 'Dif':
        Bval_label_gene.grid(row = 7, column = 4)
        TEminval_label_gene.grid(row = 8, column = 4)  
    if Pre_def_seq_gene == 'TSE':  
        showTEeff_gene()
    elif Pre_def_seq_gene == 'IN' or Pre_def_seq_gene == 'Double IN' or Pre_def_seq_gene == 'FLAIR' or Pre_def_seq_gene == 'Dif':
        post_showTEeff_gene()

#/////////////////// Functions updating the readout variable for TSE and post TSE sequences  ///////////////////
def readout_selection_gene(*args):
    """
    Assign the selected readout axis to the TSE axis readout variable
    Input : None
    Output: None
    """
    global Readout_axis_gene
    R = Readout_gene.get()

    if R == 'Along Sagittal Foot/head':
        Readout_axis_gene = 'FH'
        Post_Readout_gene.set("Along Sagittal Foot/head")
    elif R == 'Along Sagittal Anterior/posterior':
        Readout_axis_gene = 'AP'
        Post_Readout_gene.set("Along Sagittal Anterior/posterior")
    elif R == 'Along Coronal Foot/head':
        Readout_axis_gene = 'FH'
        Post_Readout_gene.set("Along Coronal Foot/head")
    elif R == 'Along Coronal Left/right':
        Readout_axis_gene = 'LR'
        Post_Readout_gene.set("Along Coronal Left/right")
    elif R == 'Along Axial Anterior/posterior':
        Readout_axis_gene = 'AP'
        Post_Readout_gene.set("Along Axial Anterior/posterior")
    elif R == 'Along Axial Left/right':
        Readout_axis_gene = 'LR'
        Post_Readout_gene.set("Along Axial Left/right")
        
def post_readout_selection_gene(*args):
    """
    Assign the selected readout axis to the post TSE axis readout variable
    Input : None
    Output: None
    """
    global Readout_axis_gene
    R = Post_Readout_gene.get()

    if R == 'Along Sagittal Foot/head':
        Readout_axis_gene = 'FH'
        Readout_gene.set("Along Sagittal Foot/head")
    elif R == 'Along Sagittal Anterior/posterior':
        Readout_axis_gene = 'AP'
        Readout_gene.set("Along Sagittal Anterior/posterior")
    elif R == 'Along Coronal Foot/head':
        Readout_axis_gene = 'FH'
        Readout_gene.set("Along Coronal Foot/head")
    elif R == 'Along Coronal Left/right':
        Readout_axis_gene = 'LR'
        Readout_gene.set("Along Coronal Left/right")
    elif R == 'Along Axial Anterior/posterior':
        Readout_axis_gene = 'AP'
        Readout_gene.set("Along Axial Anterior/posterior")
    elif R == 'Along Axial Left/right':
        Readout_axis_gene = 'LR'
        Readout_gene.set("Along Axial Left/right")
    
#//// Functions opening and closing the generator ///
def open_gene():
    """
    Opens the generator, display all labels/entries/buttons/frames important for the simulator
    Input: None
    Output: None
    """
    global T1_3D_grad
    global T2_3D_grad
    global M0_3D_grad
    global B1map_3D_grad
    global flipAngleMaprescale_3D_grad
    global t2_star_3D_grad
    global ADC_3D_grad
    global canvas1
    global canvas2
    global canvas3
    global GEN
    
    GEN = 1
    canvas1.get_tk_widget().destroy()
    canvas2.get_tk_widget().destroy()
    
    Exit_gene.grid(row = 9, column = 0) 
    frame6.grid(row = 0, column = 0, rowspan = 2)
    frame7.grid(row = 5, column = 0, columnspan = 2, rowspan = 2)
    frame8.grid(row = 5, column = 2, rowspan = 2)
    frame9.grid(row = 7, column = 2)
    Field_strength_drop_gene.grid(row = 7, column = 0)
    Update_dataset_creation.grid(row = 5, column = 3, columnspan = 2, rowspan = 2)
    Update_datasets_creation.grid(row = 7, column = 3, columnspan = 2, rowspan = 2)
    Reset_gene_button.grid(row = 2, column = 0)
    New_dataset_button.grid(row = 3, column = 0)
    
    # Read Images
    img1_gene = Image.open(resource_path("Images\Generator1.png"))
    img2_gene = Image.open(resource_path("Images\Generator2.png"))
    img3_gene = Image.open(resource_path("Images\Generator3.png"))
    
    imag_size = 5.5
    fig = plt.figure(figsize=(imag_size,imag_size))                
    plt.imshow(img1_gene)           
    plt.axis('off')
    canvas1 = FigureCanvasTkAgg(fig, root)                          
    canvas1.draw()
    canvas1.get_tk_widget().grid(row = 1, column = 1, rowspan = 4)        
    plt.close()

    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(img2_gene)   
    plt.axis('off')
    canvas2 = FigureCanvasTkAgg(fig, root)         
    canvas2.draw()
    canvas2.get_tk_widget().grid(row = 1, column = 2, rowspan = 4)        
    plt.close()

    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(img3_gene)   
    plt.axis('off')
    canvas3 = FigureCanvasTkAgg(fig, root)         
    canvas3.draw()
    canvas3.get_tk_widget().grid(row = 1, column = 3, rowspan = 4)          
    plt.close()
    
    open_simulator.grid_forget()
    open_generator.grid_forget()
    
    # Opening the generator tutorial
    tutorial_gene()
    
def exit_gene():
    """
    Closes the generator, display all labels/entries/buttons/frames important for the generator
    Input: None
    Output: None
    """
    global canvas1
    global canvas2
    global canvas3
    global GEN
    
    canvas1.get_tk_widget().destroy()
    canvas2.get_tk_widget().destroy()
    canvas3.get_tk_widget().destroy()
    GEN = 0
    
    open_simulator.grid(row = 0, column = 1)
    open_generator.grid(row = 0, column = 2)
    
    imag_size = 5.5
    fig = plt.figure(figsize=(imag_size,imag_size))                
    plt.imshow(img1)           
    plt.axis('off')
    canvas1 = FigureCanvasTkAgg(fig, root)                          
    canvas1.draw()
    canvas1.get_tk_widget().grid(row = 1, column = 1, rowspan = 4)        
    plt.close()

    fig = plt.figure(figsize=(imag_size,imag_size))
    plt.imshow(img2)   
    plt.axis('off')
    canvas2 = FigureCanvasTkAgg(fig, root)         
    canvas2.draw()
    canvas2.get_tk_widget().grid(row = 1, column = 2, rowspan = 4)        
    plt.close()
        
    Exit_gene.grid_forget()
    frame6.grid_forget()
    frame7.grid_forget()
    frame8.grid_forget()
    frame9.grid_forget()
    Field_strength_drop_gene.grid_forget()
    Update_dataset_creation.grid_forget()
    Update_datasets_creation.grid_forget()
    Reset_gene_button.grid_forget()
    New_dataset_button.grid_forget()
    Sagittal_label_gene.grid_forget() 
    Coronal_label_gene.grid_forget()
    Axial_label_gene.grid_forget()    
    
#///////////// FRAMES, BUTTONS, ENTRIES, LABELS ///////////////////////#  
    
e = 10
f = 12
# loading the 3D distorted data

T1_3D_grad = np.load(resource_path("Data\T1_3D_cer_lip_grad.npy"))
T2_3D_grad = np.load(resource_path("Data\T2_3D_cer_lip_grad.npy"))
M0_3D_grad = np.load(resource_path("Data\M0_3D_cer_lip_grad.npy"))
B1map_3D_grad = np.load(resource_path("Data\B1_3D_cer_lip_grad.npy"))
flipAngleMaprescale_3D_grad = np.load(resource_path("Data\FlipAngleMaprescale_3D_cer_lip_grad.npy"))
t2_star_3D_grad = np.load(resource_path("Data\T2_star_tensor_3D_cer_lip_grad.npy"))
ADC_3D_grad = np.load(resource_path("Data\ADC_3D_cer_lip_grad.npy"))

# //////////////////////////////////////////////////////////////////// #
#################### Frame of the sequences ############################
# //////////////////////////////////////////////////////////////////// #
frame6 = LabelFrame(root, text = "Sequences to generate", font=("Helvetica", 15))
frame6.grid(row = 0, column = 0, rowspan = 2)

# Label informing user of predefined sequence choosen
lab_gene = Label(frame6, text = "Select one sequence: ", font=("Helvetica", f)).grid(row = 0, column = 0)

# Sequence buttons
SE_button_gene = Button(frame6, text = "Spin echo", font=("Helvetica", f), command = lambda: seq_gene(SE_button_gene.cget('text')))
SE_button_gene.grid(row = 1, column = 0)
GE_button_gene = Button(frame6, text = "Gradient echo", font=("Helvetica", f), command = lambda: seq_gene(GE_button_gene.cget('text')))
GE_button_gene.grid(row = 2, column = 0)
IN_button_gene = Button(frame6, text = "Inversion recovery", font=("Helvetica", f), command = lambda: seq_gene(IN_button_gene.cget('text')))
IN_button_gene.grid(row = 3, column = 0)
double_IN_button_gene = Button(frame6, text = "Double inversion recovery", font=("Helvetica", f), command = lambda: seq_gene(double_IN_button_gene.cget('text')))
double_IN_button_gene.grid(row = 4, column = 0)
FLAIR_button_gene = Button(frame6, text = "Fluid-attenuated inversion recovery", font=("Helvetica", f), command = lambda: seq_gene(FLAIR_button_gene.cget('text')))
FLAIR_button_gene.grid(row = 5, column = 0)
TSE_button_gene = Button(frame6, text = "Turbo spin echo", font=("Helvetica", f), command = lambda: seq_gene(TSE_button_gene.cget('text')))
TSE_button_gene.grid(row = 6, column = 0)
Dif_button_gene = Button(frame6, text = "Diffusion", font=("Helvetica", f), command = lambda: seq_gene(Dif_button_gene.cget('text')))
Dif_button_gene.grid(row = 7, column = 0)
SSFP_button_gene = Button(frame6, text = "Steady-state", font=("Helvetica", f), command = lambda: seq_gene(SSFP_button_gene.cget('text')))
SSFP_button_gene.grid(row = 8, column = 0)

global Pre_def_seq_gene
Pre_def_seq_gene = " "
seq_name = " "
seq_label1_gene = Label(frame6, text = "You have chosen: ", font=("Helvetica", f)).grid(row = 9, column = 0)
seq_label2_gene = Label(frame6, text = Pre_def_seq_gene); seq_label2_gene.grid(row = 10, column = 0); seq_label2_gene.grid_forget()

# //////////////////////////////////////////////////////////////////// #
#################### Frame of the parameters ###########################
# //////////////////////////////////////////////////////////////////// #
frame7 = LabelFrame(root, text = "MRI sequence parameters", font=("Helvetica", 15))
frame7.grid(row = 5, column = 0, columnspan = 2, rowspan = 2)

# Entries of the parameters 
FOV1_entry_gene = Entry(frame7, font=("Helvetica", e));     FOV1_entry_gene.grid(row = 1, column = 1);       FOV1_entry_gene.insert(0, '250')
FOV2_entry_gene = Entry(frame7, font=("Helvetica", e));     FOV2_entry_gene.grid(row = 1, column = 3);       FOV2_entry_gene.insert(0, '300')
FOV3_entry_gene = Entry(frame7, font=("Helvetica", e));     FOV3_entry_gene.grid(row = 1, column = 5);       FOV3_entry_gene.insert(0, '275')
Res1_entry_gene = Entry(frame7, font=("Helvetica", e));     Res1_entry_gene.grid(row = 2, column = 1);       Res1_entry_gene.insert(0, '1.95')
Res2_entry_gene = Entry(frame7, font=("Helvetica", e));     Res2_entry_gene.grid(row = 2, column = 3);       Res2_entry_gene.insert(0, '2.34')
Res3_entry_gene = Entry(frame7, font=("Helvetica", e));     Res3_entry_gene.grid(row = 2, column = 5);       Res3_entry_gene.insert(0, '2.15')
Data_mat1_entry_gene = Entry(frame7, font=("Helvetica", e));Data_mat1_entry_gene.grid(row = 3, column = 1);  Data_mat1_entry_gene.insert(0, '128')
Data_mat2_entry_gene = Entry(frame7, font=("Helvetica", e));Data_mat2_entry_gene.grid(row = 3, column = 3);  Data_mat2_entry_gene.insert(0, '128')
Data_mat3_entry_gene = Entry(frame7, font=("Helvetica", e));Data_mat3_entry_gene.grid(row = 3, column = 5);  Data_mat3_entry_gene.insert(0, '128')
Bandwidth_entry_gene = Entry(frame7, font=("Helvetica", e));Bandwidth_entry_gene.grid(row = 4, column = 1);  Bandwidth_entry_gene.insert(0, '40000')
TR_entry_gene = Entry(frame7, font=("Helvetica", e));       TR_entry_gene.grid(row = 5, column = 1);         TR_entry_gene.insert(0, '500')
# Parameters dependent on the sequence
TE_entry_gene = Entry(frame7, font=("Helvetica", e));       TE_entry_gene.grid(row = 6, column = 1);         TE_entry_gene.insert(0, '20')
TI_entry_gene = Entry(frame7, font=("Helvetica", e));       TI_entry_gene.grid(row = 7, column = 1);         TI_entry_gene.insert(0, '250')
TI2_entry_gene = Entry(frame7, font=("Helvetica", e));      TI2_entry_gene.grid(row = 7, column = 3);        TI2_entry_gene.insert(0, '250')
Alpha_entry_gene = Entry(frame7, font=("Helvetica", e));    Alpha_entry_gene.grid(row = 8, column = 1);      Alpha_entry_gene.insert(0, '45')
ETL_entry_gene = Entry(frame7, font=("Helvetica", e));      ETL_entry_gene.grid(row = 9, column = 1);        ETL_entry_gene.insert(0, '8')
G_entry_gene = Entry(frame7, font=("Helvetica", e));        G_entry_gene.grid(row = 10, column = 1);          G_entry_gene.insert(0, '10')
smalldelta_entry_gene = Entry(frame7, font=("Helvetica",e));smalldelta_entry_gene.grid(row = 11, column = 1);smalldelta_entry_gene.insert(0, '1')
bigdelta_entry_gene = Entry(frame7, font=("Helvetica", e)); bigdelta_entry_gene.grid(row = 12, column = 1);  bigdelta_entry_gene.insert(0, '2')
noise_factor_entry_gene = Entry(frame7, font=("Helvetica",e))
noise_factor_entry_gene.grid(row = 4, column = 5)
noise_factor_entry_gene.insert(0, '1')

FOV_check_gene = 20625000        # 250*300*275 will be used to see if the FOV was changed by the user
Res_check_gene = 9.81045         # 1.95*2.34*2.15 will be used to see if the resolution was changed by the user
Data_matrix_check_gene = 2097152 # 128*128*128 will be used to see if the Data matrix size was changed by the user
  
# Dropdown menu to select kspace trajectory in TSE sequence   
options_gene = [
    "Linear",
    "In-out",
    "Out-in"
]
traj_gene = StringVar()
traj_gene.set(options_gene[0])                                           # Default value, or could use; options[0]
tse_drop_gene = OptionMenu(frame7, traj_gene, *options_gene, command = showTEeff_gene) # Using a list, NEEDS a star in front
tse_drop_gene.grid(row = 9, column = 1)

global Readout_axis_gene
Readout_axis_gene = 'FH'

# Dropdown menu to select readout axis in TSE sequence  
options_gene = [
    "Along Sagittal Foot/head",
    "Along Sagittal Anterior/posterior",
    "Along Coronal Foot/head",
    "Along Coronal Left/right",
    "Along Axial Anterior/posterior",
    "Along Axial Left/right"
]
Readout_gene = StringVar()
Readout_gene.set(options_gene[0])                                          
tse_read_drop_gene = OptionMenu(frame7, Readout_gene, *options_gene, command = readout_selection_gene) 
tse_read_drop_gene.grid(row = 6, column = 5)
tse_read_label_gene = Label(frame7, text = "Readout axe", font=("Helvetica", f)); tse_read_label_gene.grid(row = 6, column = 4)

###### Labels of the parameters ######
# Parameters that will always be shown
readout_label_gene = Label(frame7, text = "Readout gradient ", font=("Helvetica", f));          readout_label_gene.grid(row = 0, column = 1)
phase1_label_gene = Label(frame7, text = "Phase gradient 1", font=("Helvetica", f));            phase1_label_gene.grid(row = 0, column = 3)
phase2_label_gene = Label(frame7, text = "Phase gradient 2", font=("Helvetica", f));            phase2_label_gene.grid(row = 0, column = 5)
FOV1_label_gene = Label(frame7, text = "Field of view (mm) ", font=("Helvetica", f));           FOV1_label_gene.grid(row = 1, column = 0)
FOV2_label_gene = Label(frame7, text = "x", font=("Helvetica", f));                             FOV2_label_gene.grid(row = 1, column = 2)
FOV3_label_gene = Label(frame7, text = "x", font=("Helvetica", f));                             FOV3_label_gene.grid(row = 1, column = 4)
Resolution1_label_gene = Label(frame7, text = "Voxel resolution (mm) ", font=("Helvetica", f)); Resolution1_label_gene.grid(row = 2, column = 0)
Resolution2_label_gene = Label(frame7, text = "x", font=("Helvetica", f));                      Resolution2_label_gene.grid(row = 2, column = 2)
Resolution3_label_gene = Label(frame7, text = "x", font=("Helvetica", f));                      Resolution3_label_gene.grid(row = 2, column = 4)
Data_matrix1_label_gene = Label(frame7, text = "Data matrix ", font=("Helvetica", f));          Data_matrix1_label_gene.grid(row = 3, column = 0)
Data_matrix2_label_gene = Label(frame7, text = "x", font=("Helvetica", f));                     Data_matrix2_label_gene.grid(row = 3, column = 2)
Data_matrix3_label_gene = Label(frame7, text = "x", font=("Helvetica", f));                     Data_matrix3_label_gene.grid(row = 3, column = 4)
Bandwidth_label_gene = Label(frame7, text = "Bandwidth (Hz) ", font=("Helvetica", f));          Bandwidth_label_gene.grid(row = 4, column = 0, pady = 10)
TR_label_gene = Label(frame7, text = "TR (ms) ", font=("Helvetica", f));                        TR_label_gene.grid(row = 5, column = 0)
noise_factor_label_gene = Label(frame7, text = "noise factor ", font=("Helvetica", f));         noise_factor_label_gene.grid(row = 4, column = 4)

# Parameters dependent on the sequence
TE_label_gene = Label(frame7, text = "TE (ms) ", font=("Helvetica", f));                        TE_label_gene.grid(row = 6, column = 0)
TI_label_gene = Label(frame7, text = "TI (ms) ", font=("Helvetica", f));                        TI_label_gene.grid(row = 7, column = 0)
TI2_label_gene = Label(frame7, text = "TI2 (ms) ", font=("Helvetica", f));                      TI2_label_gene.grid(row = 7, column = 2)
Alpha_label_gene = Label(frame7, text = "Alpha ", font=("Helvetica", f));                       Alpha_label_gene.grid(row = 8, column = 0)
ETL_label_gene = Label(frame7, text = "ETL ", font=("Helvetica", f));                           ETL_label_gene.grid(row = 8, column = 0)
Kspace_traj_label_gene = Label(frame7, text = "Kspace trajectory ", font = ("Helvetica", f));   Kspace_traj_label_gene.grid(row = 8, column = 1)    
TEeff_label_gene = Label(frame7, text = "TEeff (ms) ", font = ("Helvetica", f));                TEeff_label_gene.grid(row = 8, column = 2) 
TEeff_val_label_gene = Label(frame7, text = "  ", font = ("Helvetica", f));                     TEeff_val_label_gene.grid(row = 8, column = 3)  
G_label_gene = Label(frame7, text = "G (mT/mm) ", font=("Helvetica", f));                       G_label_gene.grid(row = 10, column = 0)
smalldelta_label_gene = Label(frame7, text = "small delta (ms) ", font=("Helvetica", f));       smalldelta_label_gene.grid(row = 11, column = 0)
bigdelta_label_gene = Label(frame7, text = "big delta (ms) ", font=("Helvetica", f));           bigdelta_label_gene.grid(row = 12, column = 0)
B_label_gene = Label(frame7, text = "b (s/mm2) ", font=("Helvetica", f));                       B_label_gene.grid(row = 10, column = 3)
TEmin_label_gene = Label(frame7, text = "TEmin diffusion (ms) ", font=("Helvetica", f));        TEmin_label_gene.grid(row = 11, column = 3)
Bval_label_gene = Label(frame7, text = "302.18", font=("Helvetica", f));                        Bval_label_gene.grid(row = 10, column = 4)
TEminval_label_gene = Label(frame7, text = "3.0", font=("Helvetica", f));                       TEminval_label_gene.grid(row = 11, column = 4)

# Hide the labels and the entries of specific sequences, they will appear only if the sequence where they are involed is selected
# (Label ---- Entry ---- Dropdown)
TE_label_gene.grid_forget();         TE_entry_gene.grid_forget();           tse_drop_gene.grid_forget()
TI_label_gene.grid_forget();         TI_entry_gene.grid_forget();           tse_read_drop_gene.grid_forget()
TI2_label_gene.grid_forget();        TI2_entry_gene.grid_forget()
Alpha_label_gene.grid_forget();      Alpha_entry_gene.grid_forget()
ETL_label_gene.grid_forget();        ETL_entry_gene.grid_forget()
Kspace_traj_label_gene.grid_forget()
TEeff_label_gene.grid_forget()
TEeff_val_label_gene.grid_forget()
G_label_gene.grid_forget();          G_entry_gene.grid_forget()
smalldelta_label_gene.grid_forget(); smalldelta_entry_gene.grid_forget()
bigdelta_label_gene.grid_forget();   bigdelta_entry_gene.grid_forget()
B_label_gene.grid_forget()
Bval_label_gene.grid_forget()
TEmin_label_gene.grid_forget()       
TEminval_label_gene.grid_forget()
tse_read_label_gene.grid_forget()

# //////////////////////////////////////////////////////////////////// #
#################### Frame of the post sequence TSE ####################
# //////////////////////////////////////////////////////////////////// #
frame15 = LabelFrame(frame7, text = "TSE readout", font=("Helvetica", 15))
frame15.grid(row = 10, column = 7, rowspan = 3 , columnspan = 4)

# Labels
Post_seq_TSE_TE_label_gene = Label(frame15, text = "TE (ms) ", font=("Helvetica", f))          
Post_seq_TSE_TE_label_gene.grid(row = 0, column = 0)
Post_seq_TSE_ETL_label_gene = Label(frame15, text = "ETL ", font=("Helvetica", f))            
Post_seq_TSE_ETL_label_gene.grid(row = 1, column = 0)
post_seq_TSE_traj_label_gene = Label(frame15, text = "Kspace trajectory ", font = ("Helvetica", f))    
post_seq_TSE_traj_label_gene.grid(row = 0, column = 2) 

# Dropdown menu to select kspace trajectory in post sequence TSE   
options = [
    "Linear",
    "In-out",
    "Out-in"
]
post_traj_gene = StringVar()
post_traj_gene.set(options[0])                                          
post_tse_drop_gene = OptionMenu(frame15, post_traj_gene, *options, command = post_showTEeff_gene) 
post_tse_drop_gene.grid(row = 1, column = 2)
# Entries
post_TSE_TE_entry_gene = Entry(frame15,width=entry_width,font=("Helvetica",e));post_TSE_TE_entry_gene.grid(row=0,column=1); post_TSE_TE_entry_gene.insert(0, '20')
post_TSE_ETL_entry_gene = Entry(frame15,width=entry_width,font=("Helvetica",e));post_TSE_ETL_entry_gene.grid(row=1,column=1); post_TSE_ETL_entry_gene.insert(0, '8')

options = [
    "Along Sagittal Foot/head",
    "Along Sagittal Anterior/posterior",
    "Along Coronal Foot/head",
    "Along Coronal Left/right",
    "Along Axial Anterior/posterior",
    "Along Axial Left/right"
]

Post_Readout_gene = StringVar()
Post_Readout_gene.set(options[0])                                            # Default value, or could use; options[0]
post_tse_read_drop_gene = OptionMenu(frame15, Post_Readout_gene, *options, command = post_readout_selection_gene) # Using a list, NEEDS a star in front
post_tse_read_drop_gene.grid(row = 2, column = 2)
post_tse_read_label_gene = Label(frame15, text = "Readout axe", font=("Helvetica", f)); post_tse_read_label_gene.grid(row = 2, column = 0, columnspan = 2)

frame15.grid_forget()

# //////////////////////////////////////////////////////////////////// #
################ Frame for generating a single dataset #################
# //////////////////////////////////////////////////////////////////// #
frame8 = LabelFrame(root, text = "Single data set generation", font=("Helvetica", 15))
frame8.grid(row = 5, column = 2, rowspan = 2)

clean_data = 0
noise_data = 0

# Dataset creation button
Dataset_creation_button = Button(frame8, text = "Create data set", font=("Helvetica", f), command = create_dataset)
Dataset_creation_button.grid(row = 0, column = 0)
# Visualization of created data set button
Dataset_visu_button = Button(frame8, text = "Visualization", font=("Helvetica", f), command = visu_dataset)
Dataset_visu_button.grid(row = 0, column = 1)
# Saving buttons
Dataset_save_button = Button(frame8, text = "Save data set as numpy array", font=("Helvetica", f), command = save_dataset)
Dataset_save_button.grid(row = 0, column = 2)
Dataset_save_button_nifti = Button(frame8, text = "Save data set as nifti file", font=("Helvetica", f), command = save_dataset_nifti)
Dataset_save_button_nifti.grid(row = 1, column = 2)

# //////////////////////////////////////////////////////////////////// #
############## Frame for generating multiple datasets ##################
# //////////////////////////////////////////////////////////////////// #
frame9 = LabelFrame(root, text = "Multiple data sets generation", font=("Helvetica", 15))
frame9.grid(row = 7, column = 2)
num_seq_created_label = Label(frame9, text = 'No sequences generated', font=("Helvetica", 15))
num_seq_created_label.grid(row = 2, column = 1)

# Create the dictonaries
dict_multiple_dataset = dict() # Dictionary containing the list of all sequences (and parameters) to simulate
clean_dict = dict() # dictionary that will contain all the clean data
noise_dict = dict() # dictionary that will contain all the noisy data

# Add new sequence to the dataset to be created (adding new element to dictionary)
Add_seq_to_dict_button = Button(frame9, text = "Add sequence", font=("Helvetica", f), command = add_to_dict)
Add_seq_to_dict_button.grid(row = 0, column = 0)
View_dict_button = Button(frame9, text = "View list", font=("Helvetica", f), command = view_dict)
View_dict_button.grid(row = 0, column = 1)

# Datasets creation button
Datasets_creation_button = Button(frame9, text = "Create", font=("Helvetica", f), command = create_datasets)
Datasets_creation_button.grid(row = 2, column = 0)

# Dropdown menu to select a sequence from the multiple data sets to view
seq_to_visu = ['Choose one from list bellow']
list_to_visu_gene = StringVar()
list_to_visu_gene.set("Sequence to visualize")                         # Default value, or could use; options[0]
drop_to_visu_gene = OptionMenu(frame9, list_to_visu_gene, *seq_to_visu) # Using a list, NEEDS a star in front
drop_to_visu_gene.grid(row = 3, column = 1)

# Datasets visualization buttons
Visu_datasets_button = Button(frame9, text = "Visualize", font=("Helvetica", f), command = visu_create_datasets)
Visu_datasets_button.grid(row = 3, column = 0)

# Save dataset buttons
Save_datasets_button = Button(frame9, text = "Save data sets as numpy arrays", font=("Helvetica", f), command = save_create_datasets)
Save_datasets_button.grid(row = 2, column = 2)
Save_datasets_button_nifti = Button(frame9, text = "Save data sets as nifti files", font=("Helvetica", f), command = save_create_datasets_nifti)
Save_datasets_button_nifti.grid(row = 3, column = 2)

# Dropdown menu to view the sequences that will be generated
list_gene = ['Sequences that will be generated']
seq_list_gene = StringVar()
seq_list_gene.set("List of sequences")                         # Default value, or could use; options[0]
dict_drop_gene = OptionMenu(frame9, seq_list_gene, *list_gene) # Using a list, NEEDS a star in front
dict_drop_gene.grid(row = 1, column = 1)

# //////////////////////////////////////////////////////////////////// #
############# Buttons, entry, labels not in any frame ##################
# //////////////////////////////////////////////////////////////////// #    

Reset_gene_button = Button(root, text = "Reset parameters", font=("Helvetica", f), command = reset_param_gene)
Reset_gene_button.grid(row = 2, column = 0)

New_dataset_button = Button(root, text = "Reset for create new data set", font=("Helvetica", f), command = new_dataset_gene)
New_dataset_button.grid(row = 3, column = 0)

def B0_strength_gene(*arg):
    """
    Drop-down menu that will offer the different B0 strengh to the user
    """
    global B0_gene
    if Field_strength_gene.get() == "B0 field strength 46 mT":
        B0_gene = '46'
    elif Field_strength_gene.get() == "B0 field strength 1.5 T":
        B0_gene = '15'
    elif Field_strength_gene.get() == "B0 field strength 3 T":
        B0_gene = '3' 

# Dropdown menu to select kspace trajectory in post sequence TSE   
options = [
    "B0 field strength 46 mT",
    "B0 field strength 1.5 T",
    "B0 field strength 3 T"
]
Field_strength_gene = StringVar()
Field_strength_gene.set(options[0])                                          
Field_strength_drop_gene = OptionMenu(root, Field_strength_gene, *options, command = B0_strength_gene) 
Field_strength_drop_gene.grid(row = 7, column = 0)
B0_gene = '46'

# Label updating the user if a dataset (datasets) has (have) been created or not
Update_dataset_creation = Label(root, text = "Data set not yet created", font=("Helvetica", 15))
Update_dataset_creation.grid(row = 5, column = 3, columnspan = 2, rowspan = 2)
Update_datasets_creation = Label(root, text = "Data sets not yet created", font=("Helvetica", 15))
Update_datasets_creation.grid(row = 7, column = 3, columnspan = 2, rowspan = 2)

Sagittal_label_gene = Label(root, text = "Sagittal", font=("Helvetica", 15)); Sagittal_label_gene.grid(row = 0, column = 1)
Coronal_label_gene = Label(root, text = "Coronal", font=("Helvetica", 15)); Coronal_label_gene.grid(row = 0, column = 2)
Axial_label_gene = Label(root, text = "Axial", font=("Helvetica", 15)); Axial_label_gene.grid(row = 0, column = 3)

frame6.grid_forget()
frame7.grid_forget()
frame8.grid_forget()
frame9.grid_forget()
Field_strength_drop_gene.grid_forget()
Update_dataset_creation.grid_forget()
Update_datasets_creation.grid_forget()
Reset_gene_button.grid_forget()
New_dataset_button.grid_forget()
Sagittal_label_gene.grid_forget() 
Coronal_label_gene.grid_forget()
Axial_label_gene.grid_forget()
 
# ///////////////////////////////////////////////////////////////////////// #
### Window before simulator in which the user decides which aps to use ######
# ///////////////////////////////////////////////////////////////////////// #    

def import_t1():
    """
    Opens a search window for the user to select a T1 map
    Input: None
    Output: None
    """
    global T1_import
    root.filename = filedialog.askopenfilename() # Will return the directory of the file
    file_path_label = str(root.filename)
    T1_import = np.load(file_path_label)
    s = T1_import.shape
    T1_import = zoom(T1_import, (250/s[0], 300/s[1], 275/s[2]), order=0) # To match (250, 300, 275) shape

def import_t2():
    """
    Opens a search window for the user to select a T2 map
    Input: None
    Output: None
    """
    global T2_import
    root.filename = filedialog.askopenfilename()
    file_path_label = str(root.filename)
    T2_import = np.load(file_path_label)
    s = T2_import.shape
    T2_import = zoom(T2_import, (250/s[0], 300/s[1], 275/s[2]), order=0)

def import_b0():
    """
    Opens a search window for the user to select a B0 map
    Input: None
    Output: None
    """
    global B0_import
    Delta_B0_import = 0
    root.filename = filedialog.askopenfilename()
    file_path_label = str(root.filename)
    B0_import = np.load(file_path_label)

def open_pre_simu():
    """
    Opens a window offering the user to choice of using the predefined T1, T2, B0, etc... maps, or the user can load his/her own maps    
    Input: None
    Output: None
    """
    global canvas1
    global canvas2
    global canvas3
    global continue_to_simu
    
    canvas1.get_tk_widget().destroy()
    canvas2.get_tk_widget().destroy()
    
    open_simulator.grid_forget()
    open_generator.grid_forget()
    
    frame13.grid(row = 1, column = 0, columnspan = 3)
    frame14.grid(row = 1, column = 3, columnspan = 3)
    
    continue_to_simu.grid(row = 2, column = 0, columnspan = 6)

frame13 = LabelFrame(root, text = "Predefine simulator maps", font=("Helvetica", 15))
frame13.grid(row = 1, column = 0, columnspan = 3)

frame14 = LabelFrame(root, text = "Importable maps", font=("Helvetica", 15))
frame14.grid(row = 1, column = 3, columnspan = 3)
    
Pre_simu_label1 = Label(frame13, text = "Pre define simulator maps:", font=("Helvetica", 25)) 
Pre_simu_label2 = Label(frame13, text = "T1, T2, T2* relaxation", font=("Helvetica", 25)) 
Pre_simu_label3 = Label(frame13, text = "B0, delta B0", font=("Helvetica", 25))
Pre_simu_label4 = Label(frame13, text = "Proton density, B1", font=("Helvetica", 25))
Pre_simu_label5 = Label(frame13, text = "Apparent diffusion coefficient", font=("Helvetica", 25))
Pre_simu_label1.grid(row=0,column=0,columnspan=3)
Pre_simu_label2.grid(row=1,column=0,columnspan=3)
Pre_simu_label3.grid(row=2,column=0,columnspan=3)
Pre_simu_label4.grid(row=3,column=0,columnspan=3)
Pre_simu_label5.grid(row=4,column=0,columnspan=3)

Pre_simu_label6 = Label(frame14, text = "Import T1 map [s]", font=("Helvetica", 25))
Pre_simu_label7 = Label(frame14, text = "Import T2 map [s]", font=("Helvetica", 25))
Pre_simu_label8 = Label(frame14, text = "Import B0 map", font=("Helvetica", 25))
Pre_simu_label9 = Label(frame14, text = "Simulator assumes a FOV of 250, 300, 275 [mm]", font=("Helvetica", 25))
Pre_simu_label10 = Label(frame14, text = "Imported maps will be scaled to match those dimensions", font=("Helvetica", 25))
Pre_simu_label6.grid(row=2,column=0,columnspan=3)
Pre_simu_label7.grid(row=3,column=0,columnspan=3)
Pre_simu_label8.grid(row=4,column=0,columnspan=3)
Pre_simu_label9.grid(row=0,column=0,columnspan=3)
Pre_simu_label10.grid(row=1,column=0,columnspan=3)

import_map_T1 = Button(frame14, text = "Import", font=("Helvetica", 15), command = import_t1)
import_map_T1.grid(row = 2, column = 3)
import_map_T2 = Button(frame14, text = "Import", font=("Helvetica", 15), command = import_t2)
import_map_T2.grid(row = 3, column = 3)
import_map_B0 = Button(frame14, text = "Import", font=("Helvetica", 15), command = import_b0)
import_map_B0.grid(row = 4, column = 3)

continue_to_simu = Button(root, text = "Open simulator", font=("Helvetica", 15), command = open_simu)
continue_to_simu.grid(row = 2, column = 0, columnspan = 6)   
continue_to_simu.grid_forget()

frame13.grid_forget()
frame14.grid_forget()

# ///////////////////////////////////////////////////////////////////////// #
############## Opening window, images and buttons ###########################
# ///////////////////////////////////////////////////////////////////////// #

# M0
M0_3D = np.load(resource_path("Data\M0_3D_cerb_lip.npy"))      

# B0
global B0_3D
B0_3D = np.nan_to_num(np.load(resource_path("Data\Tom-B0.npy"))) 
B0_3D = zoom(B0_3D, (4.9, 5.88, 5.39), order=0) # (51,51,51) --> 250/51 = 4.90196078431, 300/51 = 5.88235294118, 275/51 = 5.39215686275
B0_3D = np.divide(B0_3D,1000)
B0_3D_ori = np.copy(B0_3D)

# Figures that will be there when programm is started
img1 = Image.open(resource_path("Images\Slide1.jpg"))
img2 = Image.open(resource_path("Images\Slide2.jpg"))
img1 = img1.resize((900, 800))
img2 = img2.resize((900, 800))

imag_size = 5.5
fig = plt.figure(figsize=(imag_size,imag_size))                 # Create plot
plt.imshow(img1)
plt.axis('off')
canvas1 = FigureCanvasTkAgg(fig, root)                          # Tkinter canvas which contains matplotlib figure
canvas1.draw()
canvas1.get_tk_widget().grid(row = 1, column = 1, rowspan = 4)  # Placing canvas on Tkinter window        
plt.close()

fig = plt.figure(figsize=(imag_size,imag_size))
plt.imshow(img2)   
plt.axis('off')
canvas2 = FigureCanvasTkAgg(fig, root)         
canvas2.draw()
canvas2.get_tk_widget().grid(row = 1, column = 2, rowspan = 4)        
plt.close()

# Simulator
SIM = 0 # variable keeping track if the simulator is open (1) or not (0)
open_simulator = Button(root, text = "MRI simulator", font=("Helvetica", 15), command = open_pre_simu)
open_simulator.grid(row = 0, column = 1)
# Generator
GEN = 0 # variable keeping track if the generator is open (1) or not (0)
open_generator = Button(root, text = "Data generator", font=("Helvetica", 15), command = open_gene)
open_generator.grid(row = 0, column = 2)

# Exit buttons
# Simulator
Exit_simu = Button(root, text = "Exit simulator", font=("Helvetica", f), command = exit_simu)
Exit_simu.grid(row = 9, column = 0) 
Exit_simu.grid_forget()
# Generator
Exit_gene = Button(root, text = "Exit generator", font=("Helvetica", f), command = exit_gene)
Exit_gene.grid(row = 9, column = 0) 
Exit_gene.grid_forget()

# Function called when keyboard buttons are pressed
root.bind('<Key>', UPDATE)

# To make the buttons/labels adjust to window size
frame1.columnconfigure(0, weight=1)
frame1.rowconfigure(0, weight=1)
frame1.rowconfigure(1, weight=1)
frame1.rowconfigure(2, weight=1)
frame1.rowconfigure(3, weight=1)
frame1.rowconfigure(4, weight=1)
frame1.rowconfigure(5, weight=1)
frame1.rowconfigure(6, weight=1)
frame1.rowconfigure(7, weight=1)
frame1.rowconfigure(8, weight=1)
frame2.columnconfigure(0, weight=1)
frame2.rowconfigure(0, weight=1)
frame2.rowconfigure(1, weight=1)
frame2.rowconfigure(2, weight=1)
frame3.columnconfigure(0, weight=1)
frame3.columnconfigure(1, weight=1)
frame3.columnconfigure(2, weight=1)
frame3.columnconfigure(3, weight=1)
frame3.columnconfigure(4, weight=1)
frame3.columnconfigure(5, weight=1)
frame3.columnconfigure(6, weight=1)
frame3.columnconfigure(7, weight=1)
frame3.columnconfigure(8, weight=1)
frame3.rowconfigure(0, weight=1)
frame3.rowconfigure(1, weight=1)
frame3.rowconfigure(2, weight=1)
frame3.rowconfigure(3, weight=1)
frame3.rowconfigure(4, weight=1)
frame3.rowconfigure(5, weight=1)
frame3.rowconfigure(6, weight=1)
frame3.rowconfigure(7, weight=1)
frame4.columnconfigure(0, weight=1)
frame4.columnconfigure(1, weight=1)
frame4.columnconfigure(2, weight=1)
frame4.columnconfigure(3, weight=1)
frame4.columnconfigure(4, weight=1)
frame4.rowconfigure(0, weight=1)
frame4.rowconfigure(1, weight=1)
frame4.rowconfigure(2, weight=1)
frame4.rowconfigure(3, weight=1)
frame4.rowconfigure(4, weight=1)
frame4.rowconfigure(5, weight=1)
frame4.rowconfigure(6, weight=1)
frame4.rowconfigure(7, weight=1)
frame5.columnconfigure(0, weight=1)
frame5.columnconfigure(1, weight=1)
frame5.columnconfigure(2, weight=1)
frame5.rowconfigure(0, weight=1)
frame5.rowconfigure(1, weight=1)
frame6.columnconfigure(0, weight=1)
frame6.rowconfigure(0, weight=1)
frame6.rowconfigure(1, weight=1)
frame6.rowconfigure(2, weight=1)
frame6.rowconfigure(3, weight=1)
frame6.rowconfigure(4, weight=1)
frame6.rowconfigure(5, weight=1)
frame6.rowconfigure(6, weight=1)
frame6.rowconfigure(7, weight=1)
frame6.rowconfigure(8, weight=1)
frame6.rowconfigure(9, weight=1)
frame6.rowconfigure(10, weight=1)
frame7.columnconfigure(0, weight=1)
frame7.columnconfigure(1, weight=1)
frame7.columnconfigure(2, weight=1)
frame7.columnconfigure(3, weight=1)
frame7.columnconfigure(4, weight=1)
frame7.columnconfigure(5, weight=1)
frame7.rowconfigure(0, weight=1)
frame7.rowconfigure(1, weight=1)
frame7.rowconfigure(2, weight=1)
frame7.rowconfigure(3, weight=1)
frame7.rowconfigure(4, weight=1)
frame7.rowconfigure(5, weight=1)
frame7.rowconfigure(6, weight=1)
frame7.rowconfigure(7, weight=1)
frame7.rowconfigure(8, weight=1)
frame7.rowconfigure(9, weight=1)
frame7.rowconfigure(10, weight=1)
frame7.rowconfigure(11, weight=1)
frame8.columnconfigure(0, weight=1)
frame8.columnconfigure(1, weight=1)
frame8.columnconfigure(2, weight=1)
frame8.rowconfigure(0, weight=1)
frame9.columnconfigure(0, weight=1)
frame9.columnconfigure(1, weight=1)
frame9.columnconfigure(2, weight=1)
frame9.rowconfigure(0, weight=1)
frame9.rowconfigure(1, weight=1)
frame9.rowconfigure(2, weight=1)
frame9.rowconfigure(3, weight=1)
frame10.columnconfigure(0, weight=1)
frame10.rowconfigure(0, weight=1)
frame10.rowconfigure(1, weight=1)
frame10.rowconfigure(2, weight=1)
frame10.rowconfigure(3, weight=1)
frame10.rowconfigure(4, weight=1)
frame11.columnconfigure(0, weight=1)
frame11.columnconfigure(1, weight=1)
frame11.rowconfigure(0, weight=1)
frame11.rowconfigure(1, weight=1)
frame11.rowconfigure(2, weight=1)
frame11.rowconfigure(3, weight=1)
frame12.columnconfigure(0, weight=1)
frame12.columnconfigure(1, weight=1)
frame12.columnconfigure(2, weight=1)
frame12.rowconfigure(0, weight=1) 
frame15.columnconfigure(0, weight=1)
frame15.columnconfigure(1, weight=1)
frame15.columnconfigure(2, weight=1)
frame15.rowconfigure(0, weight=1)
frame15.rowconfigure(1, weight=1)
root.columnconfigure(0, weight=1)
root.columnconfigure(1, weight=1)
root.columnconfigure(2, weight=1)
root.columnconfigure(3, weight=1)
root.columnconfigure(4, weight=1)
root.columnconfigure(5, weight=1)
root.rowconfigure(0, weight=1)
root.rowconfigure(1, weight=1)
root.rowconfigure(2, weight=1)
root.rowconfigure(3, weight=1)
root.rowconfigure(4, weight=1)
root.rowconfigure(5, weight=1)
root.rowconfigure(6, weight=1)
root.rowconfigure(7, weight=1)
root.rowconfigure(8, weight=1)
root.rowconfigure(9, weight=1)
root.rowconfigure(10, weight=1)
root.rowconfigure(11, weight=1)
root.rowconfigure(12, weight=1)
root.rowconfigure(13, weight=1)
root.rowconfigure(14, weight=1)
root.rowconfigure(15, weight=1)

root.mainloop()

