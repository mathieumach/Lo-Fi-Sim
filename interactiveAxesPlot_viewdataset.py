import numpy as np
import matplotlib as plt
from skimage.restoration import denoise_nl_means, estimate_sigma # estimate_sigma --> wavelet-based estimator of the (Gaussian) noise standard deviation.
from matplotlib.widgets import Slider, Button

class interactivePlot(object):
    def __init__(self,fig, ax, X, Y, data_mat, plotAxis = 2, axisLabels = None, fov = None, cmap = 'gray'):
        
        #        set various event handlers
        fig.canvas.mpl_connect('scroll_event', self.onScroll)
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
        self.cmapRangey = np.nanmax(Y) - np.nanmin(Y)
        self.cmapCentery = np.nanmin(Y) + self.cmapRangey/2
        self.tempCmapRangey = self.cmapRangey
        self.tempCmapCentery = self.cmapCentery
        self.X = X
        self.Y = Y
        self.axisLabels = axisLabels
        self.colorMap = cmap
        self.original = np.copy(X)
        self.index = 0
        self.data_mat = data_mat
        
        if fov is None:
            self.fov = (1,1,1)
            self.resolution = (1,1,1)
        else:
            self.fov = fov
            self.resolution = np.divide(fov, np.shape(X))

        self.shape = np.shape(X)
        x = self.shape[plotAxis]
        self.ind = x//2
        
        # To change the axis value to mm (the fov)
        ax_pos = [[0, data_mat[0]/2, data_mat[0]-1], [0, data_mat[1]/2, data_mat[1]-1], [0, data_mat[2]/2, data_mat[2]-1]]
        ax_label = [[0, self.fov[0]/2, self.fov[0]], [0, self.fov[1]/2, self.fov[1]], [0, self.fov[2]/2, self.fov[2]]]
        
        if self.plotAxis == 0:
            imageDatax = self.X[self.ind,:,:]
            imageDatay = self.Y[self.ind,:,:]
            self.slices = self.shape[0]
        elif self.plotAxis == 1:
            imageDatax = self.X[:,self.ind,:]
            imageDatay = self.Y[:,self.ind,:]
            self.slices = self.shape[1]
        elif self.plotAxis == 2:
            imageDatax = self.X[:,:,self.ind]
            imageDatay = self.Y[:,:,self.ind]
            self.slices = self.shape[2]
        else:
            print("invalid axis")
            return -1
        ax[0].set_title('Clean data')
        ax[1].set_title('Noisy data')
        self.imx = self.ax[0].imshow(np.rot90(imageDatax), cmap = self.colorMap,vmin = self.cmapCenter-self.cmapRange/2, vmax= self.cmapCenter+self.cmapRange/2)
        self.imy = self.ax[1].imshow(np.rot90(imageDatay), cmap = self.colorMap,vmin = self.cmapCentery-self.cmapRangey/2, vmax= self.cmapCentery+self.cmapRangey/2)
        if self.axisLabels != None:
            self.ax[0].set_xlabel(self.axisLabels[1])
            self.ax[0].set_ylabel(self.axisLabels[0])
            self.ax[1].set_xlabel(self.axisLabels[1])
            self.ax[1].set_ylabel(self.axisLabels[0])

        # Global title
        self.fig.suptitle('Axial axis, Slice ' + str(self.ind) + '/' + str(self.slices), fontsize=14)  
            
    def keyPress(self, event):
        if event.key == " ": #change orientation on space bar press
            self.plotAxis += 1
            if self.plotAxis == 1:
                self.slices = self.shape[1]
                self.fig.suptitle('Coronal axis, Slice ' + str(self.ind) + '/' + str(self.slices), fontsize=14)
            elif self.plotAxis == 2:
                self.slices = self.shape[2]
                self.fig.suptitle('Axial axis, Slice ' + str(self.ind) + '/' + str(self.slices), fontsize=14)
            if self.plotAxis > 2:
                self.plotAxis = 0
                self.slices = self.shape[0]
                self.fig.suptitle('Sagittal axis, Slice ' + str(self.ind) + '/' + str(self.slices), fontsize=14)
            self.updateSlice()
            
    def onScroll(self, event):
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.updateSlice()
    
    def onClick(self, event):
        if event.button == 2: #change orientation on scroll wheel click
            self.plotAxis += 1
            if self.plotAxis > 2:
                self.plotAxis = 0
            self.updateSlice()
        else:
            self.mouseClicked = event.xdata, event.ydata
        
    def onRelease(self, event):
        self.mouseClicked = None
        self.cmapRange = self.tempCmapRange
        self.cmapCenter = self.tempCmapCenter
        self.cmapRangey = self.tempCmapRangey
        self.cmapCentery = self.tempCmapCentery
        
    def onMotion(self, event):
        if self.mouseClicked == None: return        #if mouse isn't clicked ignore movement

        dx = event.xdata - self.mouseClicked[0]
        dy = event.ydata - self.mouseClicked[1]
        
        normDx = dx/self.mouseClicked[0]
        normDy = dy/self.mouseClicked[1]
        
        self.tempCmapRange = self.cmapRange*(1+normDy)
        self.tempCmapCenter = self.cmapCenter*(1+normDx)
        self.tempCmapRangey = self.cmapRangey*(1+normDy)
        self.tempCmapCentery = self.cmapCentery*(1+normDx)
        
        self.imx.set_clim(self.tempCmapCenter-self.tempCmapRange/2, self.tempCmapCenter+self.tempCmapRange/2)
        self.imx.axes.figure.canvas.draw()
        self.imy.set_clim(self.tempCmapCentery-self.tempCmapRangey/2, self.tempCmapCentery+self.tempCmapRangey/2)
        self.imy.axes.figure.canvas.draw()
        
    def updateSlice(self):
        if self.plotAxis == 0:
            self.slices = self.shape[0]
            if self.ind > self.shape[0]:
                self.ind = int(self.shape[0]/2)
            self.fig.suptitle('Sagittal axis, Slice ' + str(self.ind) + '/' + str(self.slices), fontsize=14)
            imageDataX = np.flip(np.rot90(self.original[self.ind,:,:]), axis = 1)
            imageDataY = np.flip(np.rot90(self.Y[self.ind,:,:]), axis = 1)
            if self.axisLabels != None:
                self.ax[0].set_xlabel(self.axisLabels[2])
                self.ax[0].set_ylabel(self.axisLabels[1])            
                self.ax[1].set_xlabel(self.axisLabels[2])
                self.ax[1].set_ylabel(self.axisLabels[1])
        elif self.plotAxis == 1:
            self.slices = self.shape[1]
            if self.ind > self.shape[1]:
                self.ind = int(self.shape[1]/2)
            self.fig.suptitle('Coronal axis, Slice ' + str(self.ind) + '/' + str(self.slices), fontsize=14)
            imageDataX = np.rot90(self.original[:,self.ind,:])
            imageDataY = np.rot90(self.Y[:,self.ind,:])
            if self.axisLabels != None:
                self.ax[0].set_xlabel(self.axisLabels[2])
                self.ax[0].set_ylabel(self.axisLabels[0])
                self.ax[1].set_xlabel(self.axisLabels[2])
                self.ax[1].set_ylabel(self.axisLabels[0])
        else:
            self.slices = self.shape[2]
            if self.ind > self.shape[2]:
                self.ind = int(self.shape[2]/2)
            self.fig.suptitle('Axial axis, Slice ' + str(self.ind) + '/' + str(self.slices), fontsize=14)
            imageDataX = np.rot90(self.original[:,:,self.ind])
            imageDataY = np.rot90(self.Y[:,:,self.ind])
            if self.axisLabels != None:
                self.ax[0].set_xlabel(self.axisLabels[1])
                self.ax[0].set_ylabel(self.axisLabels[0])
                self.ax[1].set_xlabel(self.axisLabels[1])
                self.ax[1].set_ylabel(self.axisLabels[0])
                
        self.imx.set_extent((-0.5,np.size(imageDataX,1) - 0.5,np.size(imageDataX,0) - 0.5,-0.5))
        self.imx.set_data(imageDataX)
        self.imx.axes.figure.canvas.draw()  
        
        self.imy.set_extent((-0.5,np.size(imageDataY,1) - 0.5,np.size(imageDataY,0) - 0.5,-0.5))
        self.imy.set_data(imageDataY)
        self.imy.axes.figure.canvas.draw()




