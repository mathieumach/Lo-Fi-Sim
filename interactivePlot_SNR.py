import numpy as np
from matplotlib.widgets import Slider, Button, RadioButtons

class interactivePlot(object):
    def __init__(self,fig, ax, X, Data_mat, func, plotAxis = 2, axisLabels = None, fov = None, cmap = 'gray'):
        
#        set various event handlers
        fig.canvas.mpl_connect('button_press_event', self.onClick)
        fig.canvas.mpl_connect('button_release_event', self.onRelease)
        fig.canvas.mpl_connect('motion_notify_event', self.onMotion)
        fig.canvas.mpl_connect('key_press_event', self.keyPress)
        
        self.fig = fig
        self.plotAxis = plotAxis
        self.ax = ax
        self.mouseClicked = None
        self.cmapRange = np.nanmax(X) - np.nanmin(X)
        self.cmapCenter = np.nanmin(X) + self.cmapRange/2
        self.tempCmapRange = self.cmapRange
        self.tempCmapCenter = self.cmapCenter
        self.X = X
        self.Data_mat = Data_mat
        self.axisLabels = axisLabels
        self.colorMap = cmap
        self.box = 'off'
        self.std_box_defined = 'no'
        self.mean_box_defined = 'no'
        self.noise_box_center = [0,0]
        self.mean_box_center = [0,0]
        self.snr = 0
        self.func = func
        self.length = 5
        self.center = [int(Data_mat[0]/2), int(Data_mat[1]/2)]
        self.fig.suptitle("SNR of spin echo axial image ", fontsize=14)
        
        imageData = np.rot90(self.X)
        self.im = ax.imshow(imageData, cmap = self.colorMap,vmin = self.cmapCenter-self.cmapRange/2, vmax= self.cmapCenter+self.cmapRange/2)
        
        # create sliders       
        box_size = np.concatenate([np.linspace(0, 30, 31)])
        self.axSlider_box = fig.add_axes([0.25, 0.01, 0.6, 0.03])
        self.slider_box_size = Slider(ax=self.axSlider_box, label='box size', valmin=1, valstep = box_size, valmax = box_size[-1], valinit=5)  
        
        def update(val):    
            
            if self.box != 'okay':
                print('the std and mean boxes must be defined before computing the SNR!')
            elif self.box == 'okay':
                
                if self.std_box_defined == 'yes' and self.mean_box_defined == 'yes':
                
                    # box size 
                    if self.slider_box_size.val == 0:
                        m = 5
                    else:
                        m = self.slider_box_size.val
                    m = int(m)
                    M = 2*m
                    self.length = m
                    
                    arrstd = np.ones(M)
                    arrm = np.ones(M)
                    arr2 = np.zeros((Data_mat[1], Data_mat[0])) # This array is only zeros with 2 white boxes
                    # std box
                    A = int(self.noise_box_center[1])
                    B = int(self.noise_box_center[0])
                    arr2[A-m:A+m, B-m] = arrstd
                    arr2[A-m:A+m, B+m] = arrstd
                    arr2[A-m, B-m:B+m] = arrstd
                    arr2[A+m, B-m:B+m] = arrstd
                    # mean box
                    a = int(self.mean_box_center[1])
                    b = int(self.mean_box_center[0])
                    arr2[a-m:a+m, b-m] = arrm
                    arr2[a-m:a+m, b+m] = arrm
                    arr2[a-m, b-m:b+m] = arrm
                    arr2[a+m, b-m:b+m] = arrm

                    self.arr2 = arr2
                    self.snr = self.func(np.rot90(self.X),A-m,A+m,B-m,B+m,a-m,a+m,b-m,b+m)
                    self.fig.suptitle("SNR of image " + str(np.round(self.snr,2)), fontsize=14)
                    self.ax.set_title('Right click to start')
                    
                    if self.axisLabels != None:
                        self.ax.set_xlabel(self.axisLabels[1])
                        self.ax.set_ylabel(self.axisLabels[0])
                    self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
                    self.im.set_data(imageData + 10000*self.arr2)
                    self.im.axes.figure.canvas.draw()
                    
                else:
                    print('the std and mean boxes must be defined before computing the SNR!')
                    
            #return self.length, self.noise_box_center, self.mean_box_center
        
        self.slider_box_size.on_changed(update)        
        self.updateSlice()
    
    def keyPress(self, event):
        if event.key == " ": #change orientation on space bar press
            self.plotAxis += 1
            if self.plotAxis > 2:
                self.plotAxis = 0
            self.updateSlice()
    
    def onClick(self, event):
        if event.button == 3:   #Reset image on right click
            if self.box == 'off':
                self.box = 'std'
                self.ax.set_title('Select center of noise box (for standard deviation)')
                self.im.axes.figure.canvas.draw()
                
            elif self.box == 'std':
                self.box = 'mean'
                self.ax.set_title('Select center of mean box')
                self.im.axes.figure.canvas.draw()
                
            elif self.box == 'mean':
                self.box = 'okay'
                self.ax.set_title('SNR will is defined based on the boxes')
                self.updateSlice()
                self.im.axes.figure.canvas.draw()
                
            elif self.box == 'okay':
                self.box = 'off'
                self.std_box_defined = 'no'
                self.mean_box_defined = 'no'
                self.updateSlice()
                
            print('Action is ' + str(self.box))
            
        else:
            self.mouseClicked = event.xdata, event.ydata
            if self.box == 'std':
                self.std_box_defined = 'yes'
                self.noise_box_center[0] = np.round(self.mouseClicked[0])
                self.noise_box_center[1] = np.round(self.mouseClicked[1])
            elif self.box == 'mean':
                self.mean_box_defined = 'yes'
                self.mean_box_center[0] = np.round(self.mouseClicked[0])
                self.mean_box_center[1] = np.round(self.mouseClicked[1])
                
                m = 5
                A = int(self.noise_box_center[1])
                B = int(self.noise_box_center[0])
                a = int(self.mean_box_center[1])
                b = int(self.mean_box_center[0])
                
                self.snr = self.func(np.rot90(self.X),A-m,A+m,B-m,B+m, a-m,a+m,b-m,b+m)
                self.fig.suptitle("SNR of image " + str(np.round(self.snr,2)), fontsize=14)
            self.updateSlice()
        
    def onRelease(self, event):
        self.mouseClicked = None
        self.cmapRange = self.tempCmapRange
        self.cmapCenter = self.tempCmapCenter
        
    def onMotion(self, event):
        if self.mouseClicked == None: return        #if mouse isn't clicked ignore movement

        dx = event.xdata - self.mouseClicked[0]
        dy = event.ydata - self.mouseClicked[1]
        
        normDx = dx/self.mouseClicked[0]
        normDy = dy/self.mouseClicked[1]
        
        self.tempCmapRange = self.cmapRange*(1+normDy)
        self.tempCmapCenter = self.cmapCenter*(1+normDx)
        
        self.im.set_clim(self.tempCmapCenter-self.tempCmapRange/2, self.tempCmapCenter+self.tempCmapRange/2)
        self.im.axes.figure.canvas.draw()
        
    def updateSlice(self):
        imageData = np.rot90(self.X)
        
        if self.box == 'off': # Will not draw the boxes
            if self.axisLabels != None:
                self.ax.set_xlabel(self.axisLabels[1])
                self.ax.set_ylabel(self.axisLabels[0])
            self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
            self.im.set_data(imageData)
            self.ax.set_title('Right click to start')
            self.im.axes.figure.canvas.draw()
        
        elif self.box != 'okay':             
            if self.std_box_defined == 'no' and self.mean_box_defined == 'no':
                if self.axisLabels != None:
                    self.ax.set_xlabel(self.axisLabels[1])
                    self.ax.set_ylabel(self.axisLabels[0])
                self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
                self.im.set_data(imageData)
                self.ax.set_title('Right click to start')
                self.im.axes.figure.canvas.draw()          
            
            elif self.std_box_defined == 'yes' and self.mean_box_defined == 'no':
                
                # Creation of the boxes with center defined by user but start size is 10x10
                M = 10; m = 5
                arrstd = np.ones(M)
                arr2 = np.zeros((self.Data_mat[1], self.Data_mat[0])) # This array is only zeros with 2 white boxes
                # std box
                A = int(self.noise_box_center[1]); B = int(self.noise_box_center[0])
                arr2[A-m:A+m, B-m] = arrstd; arr2[A-m:A+m, B+m] = arrstd
                arr2[A-m, B-m:B+m] = arrstd; arr2[A+m, B-m:B+m] = arrstd
                self.arr2 = arr2
                if self.axisLabels != None:
                    self.ax.set_xlabel(self.axisLabels[1])
                    self.ax.set_ylabel(self.axisLabels[0])
                self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
                self.im.set_data(imageData + 10000*self.arr2)
                self.im.axes.figure.canvas.draw()
                
            elif self.std_box_defined == 'yes' and self.mean_box_defined == 'yes':                
                # Creation of the boxes with center defined by user but start size is 10x10
                M = 10; m = 5
                arrstd = np.ones(M); arrm = np.ones(M)
                arr2 = np.zeros((self.Data_mat[1], self.Data_mat[0])) # This array is only zeros with 2 white boxes
                # std box
                A = int(self.noise_box_center[1]); B = int(self.noise_box_center[0])
                arr2[A-m:A+m, B-m] = arrstd; arr2[A-m:A+m, B+m] = arrstd
                arr2[A-m, B-m:B+m] = arrstd; arr2[A+m, B-m:B+m] = arrstd
                # mean box
                a = int(self.mean_box_center[1]); b = int(self.mean_box_center[0])
                arr2[a-m:a+m, b-m] = arrm; arr2[a-m:a+m, b+m] = arrm
                arr2[a-m, b-m:b+m] = arrm; arr2[a+m, b-m:b+m] = arrm
                self.arr2 = arr2
                if self.axisLabels != None:
                    self.ax.set_xlabel(self.axisLabels[1])
                    self.ax.set_ylabel(self.axisLabels[0])
                self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
                self.im.set_data(imageData + 10000*self.arr2)
                self.im.axes.figure.canvas.draw()
              
        elif self.box == 'okay': # Will draw both boxes
            # Creation of the boxes with center defined by user but start size is 10x10
            M = 10; m = 5
            arrstd = np.ones(M); arrm = np.ones(M)
            arr2 = np.zeros((self.Data_mat[1], self.Data_mat[0])) # This array is only zeros with 2 white boxes
            # std box
            A = int(self.noise_box_center[1]); B = int(self.noise_box_center[0])
            arr2[A-m:A+m, B-m] = arrstd; arr2[A-m:A+m, B+m] = arrstd
            arr2[A-m, B-m:B+m] = arrstd; arr2[A+m, B-m:B+m] = arrstd
            # mean box
            a = int(self.mean_box_center[1]); b = int(self.mean_box_center[0])
            arr2[a-m:a+m, b-m] = arrm; arr2[a-m:a+m, b+m] = arrm
            arr2[a-m, b-m:b+m] = arrm; arr2[a+m, b-m:b+m] = arrm
            self.arr2 = arr2
            if self.axisLabels != None:
                self.ax.set_xlabel(self.axisLabels[1])
                self.ax.set_ylabel(self.axisLabels[0])
            self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
            self.im.set_data(imageData + 10000*self.arr2)
            self.im.axes.figure.canvas.draw()



