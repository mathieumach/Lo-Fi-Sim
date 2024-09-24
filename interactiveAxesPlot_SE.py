import numpy as np
import matplotlib.pyplot as plt 
from Sequences import *
from matplotlib.widgets import Slider, Button, RadioButtons

class interactivePlot(object):   
    def __init__(self,fig, ax, X, TR, TE, T1, T2, M0, B1, noise, data_mat, plotAxis = 2, axisLabels = None, fov = None, cmap = 'gray'):
        
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
        self.X = X
        self.axisLabels = ['Left/Right [mm]', 'Anterior/Posterior [mm]', 'Foot/Head [mm]']
        self.colorMap = cmap
        self.original = X
        self.index = 0
        self.TR = TR
        self.TE = TE
        self.T1 = T1
        self.T2 = T2
        self.M0 = M0
        self.noise = noise
        self.B1 = B1
        self.data_mat = data_mat
        
        if fov is None:
            self.fov = (1,1,1)
            self.resolution = (1,1,1)
        else:
            self.fov = fov
            self.resolution = np.divide(fov, np.shape(X))
            
        self.shape = np.shape(X)
        self.slices = np.shape(X)[plotAxis]
        self.ind = self.slices//2

        self.x_ind = np.shape(X)[0]//2
        self.y_ind = np.shape(X)[1]//2
        self.z_ind = np.shape(X)[2]//2
        
        # Global title
        self.fig.suptitle("Action will be used on sagittal axis", fontsize=14)
        
        if TE > TR:
            self.X = np.zeros((data_mat))
            self.fig.suptitle("TE can't be bigger than TR!", fontsize=14)
        
        imageDataX = self.X[self.x_ind,:,:]
        imageDataY = self.X[:,self.y_ind,:]
        imageDataZ = self.X[:,:,self.z_ind]
        
        # To change the axis value to mm (the fov)
        ax_pos = [[0, data_mat[0]/2, data_mat[0]-1], [0, data_mat[1]/2, data_mat[1]-1], [0, data_mat[2]/2, data_mat[2]-1]]
        ax_label = [[0, self.fov[0]/2, self.fov[0]], [0, self.fov[1]/2, self.fov[1]], [0, self.fov[2]/2, self.fov[2]]]
        
        self.imx = ax[0].imshow(np.flip(np.rot90(imageDataX), axis = 1), cmap = self.colorMap,vmin = self.cmapCenter-self.cmapRange/2, vmax= self.cmapCenter+self.cmapRange/2)
        self.ax[0].set_xlabel(self.axisLabels[1]); self.ax[0].set_xticks(ax_pos[1]); self.ax[0].set_xticklabels(ax_label[1])
        self.ax[0].set_ylabel(self.axisLabels[2]); self.ax[0].set_yticks(ax_pos[2]); self.ax[0].set_yticklabels(ax_label[2])
        self.ax[0].set_title('Slice sag %s/%s' % (self.x_ind,self.shape[0]))   
        
        self.imy = ax[1].imshow(np.rot90(imageDataY), cmap = self.colorMap,vmin = self.cmapCenter-self.cmapRange/2, vmax= self.cmapCenter+self.cmapRange/2)
        self.ax[1].set_xlabel(self.axisLabels[0]); self.ax[1].set_xticks(ax_pos[0]); self.ax[1].set_xticklabels(ax_label[0])
        self.ax[1].set_ylabel(self.axisLabels[2]); self.ax[1].set_yticks(ax_pos[2]); self.ax[1].set_yticklabels(ax_label[2])
        self.ax[1].set_title('Slice cor %s/%s' % (self.y_ind,self.shape[1]))
        
        self.imz = ax[2].imshow(np.rot90(imageDataZ), cmap = self.colorMap,vmin = self.cmapCenter-self.cmapRange/2, vmax= self.cmapCenter+self.cmapRange/2)
        self.ax[2].set_xlabel(self.axisLabels[0]); self.ax[2].set_xticks(ax_pos[0]); self.ax[2].set_xticklabels(ax_label[0])
        self.ax[2].set_ylabel(self.axisLabels[1]); self.ax[2].set_yticks(ax_pos[1]); self.ax[2].set_yticklabels(ax_label[1])
        self.ax[2].set_title('Slice ax %s/%s' % (self.z_ind,self.shape[2]))
        
        #create sliders
        # TR
        self.axSlider_TR = fig.add_axes([0.15, 0.08, 0.6, 0.03])
        self.slider_TR = Slider(ax=self.axSlider_TR, label='TR', valmin=0, valmax=5000, valinit=self.TR*1000)
        # TE
        self.axSlider_TE = fig.add_axes([0.15, 0.04, 0.6, 0.03])
        self.slider_TE = Slider(ax=self.axSlider_TE, label='TE', valmin=0, valmax=1000, valinit=self.TE*1000)        
        
        def update(val):

            # TR
            if self.slider_TR.val == 0:
                tr = self.TR
            else:
                tr = np.divide(self.slider_TR.val,1000)
            # TE
            if self.slider_TE.val == 0:
                te = self.TE
            else:
                te = np.divide(self.slider_TE.val,1000)
            
            if te < tr:
                self.X = np.multiply(spin_echo_seq(tr, te, self.T1, self.T2, self.M0),self.B1)
                self.X = self.X + self.noise
                if self.index == 0:
                    self.fig.suptitle("Action will be used on sagittal axis", fontsize=14)
                elif self.index == 1:
                    self.fig.suptitle("Action will be used on coronal axis", fontsize=14)
                elif self.index == 2:
                    self.fig.suptitle("Action will be used on axial axis", fontsize=14)
            elif te > tr:
                self.X = np.zeros((data_mat))
                self.fig.suptitle("TE can't be bigger than TR!", fontsize=14)
                
            self.updateAllSlice()
            
        self.slider_TR.on_changed(update)
        self.slider_TE.on_changed(update)

        #create radiobutton
        self.axRadio = fig.add_axes([0.85, 0.8, 0.15, 0.1])
        self.radio = RadioButtons(self.axRadio, ('no cross line', 'cross line'))
        self.plane = 'no cross line'
        
        def planefunc(label):
            self.plane = label
            if self.plane == 'cross line':
                
                imageDataX = np.flip(np.rot90(self.X[self.x_ind,:,:]), axis = 1)
                s = np.flip(np.rot90(self.X[self.x_ind,:,:]), axis = 1).shape
                arr1 = np.ones(s[0])
                arr2 = np.ones(s[1])
                arr3 = np.zeros((s)) # This array will filled with zeros with a vertical and horizontal white lines
                arr3[0, :] = arr2
                arr3[:, 0] = arr1
                if self.x_ind >= np.size(self.X,0):
                    self.x_ind = np.size(self.X,0)-1
                if self.axisLabels != None:
                    self.ax[0].set_xlabel(self.axisLabels[1])
                    self.ax[0].set_ylabel(self.axisLabels[2])
                self.imx.set_extent((-0.5,np.size(imageDataX,1) - 0.5,np.size(imageDataX,0) - 0.5,-0.5))
                #self.ax[0].set_aspect(self.resolution[1]/self.resolution[2])
                self.slices = np.size(self.X,0)
                self.imx.set_data(imageDataX + 10000*arr3)
                self.ax[0].set_title('Slice sag %s/%s' % (self.x_ind,self.slices))
                self.imx.axes.figure.canvas.draw()
                
                imageDataY = np.rot90(self.X[:,self.y_ind,:])
                s = np.rot90(self.X[:,self.y_ind,:]).shape
                arr1 = np.ones(s[0])
                arr2 = np.ones(s[1])
                arr3 = np.zeros((s)) # This array is only zeros with a vertical and horizontal white lines
                arr3[0, :] = arr2
                arr3[:, 0] = arr1
                if self.y_ind >= np.size(self.X,1):
                    self.y_ind = np.size(self.X,1)-1
                if self.axisLabels != None:
                    self.ax[1].set_xlabel(self.axisLabels[0])
                    self.ax[1].set_ylabel(self.axisLabels[2])
                self.imy.set_extent((-0.5,np.size(imageDataY,1) - 0.5,np.size(imageDataY,0) - 0.5,-0.5))
                #self.ax[1].set_aspect(self.resolution[0]/self.resolution[2])
                self.slices = np.size(self.X,1)
                self.imy.set_data(imageDataY + 10000*arr3)
                self.ax[1].set_title('Slice cor %s/%s' % (self.y_ind,self.slices))
                self.imy.axes.figure.canvas.draw()
                
                imageDataZ = np.rot90(self.X[:,:,self.z_ind])
                s = np.rot90(self.X[:,:,self.z_ind]).shape
                arr1 = np.ones(s[0])
                arr2 = np.ones(s[1])
                arr3 = np.zeros((s)) # This array is only zeros with a vertical and horizontal white lines
                arr3[0, :] = arr2
                arr3[:, 0] = arr1
                if self.z_ind >= np.size(self.X,2):
                    self.z_ind = np.size(self.X,2)-1
                if self.axisLabels != None:
                    self.ax[2].set_xlabel(self.axisLabels[0])
                    self.ax[2].set_ylabel(self.axisLabels[1])
                self.imz.set_extent((-0.5,np.size(imageDataZ,1) - 0.5,np.size(imageDataZ,0) - 0.5,-0.5))
                #self.ax[2].set_aspect(self.resolution[0]/self.resolution[1])
                self.slices = np.size(self.X,2)
                self.imz.set_data(imageDataZ + 10000*arr3)
                self.ax[2].set_title('Slice ax %s/%s' % (self.z_ind,self.slices))
                self.imz.axes.figure.canvas.draw()    
                
            elif self.plane == 'no cross line':
                imageDataX = np.flip(np.rot90(self.X[self.x_ind,:,:]), axis = 1)
                if self.x_ind >= np.size(self.X,0):
                    self.x_ind = np.size(self.X,0)-1
                if self.axisLabels != None:
                    self.ax[0].set_xlabel(self.axisLabels[1])
                    self.ax[0].set_ylabel(self.axisLabels[2])
                self.imx.set_extent((-0.5,np.size(imageDataX,1) - 0.5,np.size(imageDataX,0) - 0.5,-0.5))
                #self.ax[0].set_aspect(self.resolution[1]/self.resolution[2])
                self.slices = np.size(self.X,0)
                self.imx.set_data(imageDataX)
                self.ax[0].set_title('Slice sag %s/%s' % (self.x_ind,self.slices))
                self.imx.axes.figure.canvas.draw()
                
                imageDataY = np.rot90(self.X[:,self.y_ind,:])
                if self.y_ind >= np.size(self.X,1):
                    self.y_ind = np.size(self.X,1)-1
                if self.axisLabels != None:
                    self.ax[1].set_xlabel(self.axisLabels[0])
                    self.ax[1].set_ylabel(self.axisLabels[2])
                self.imy.set_extent((-0.5,np.size(imageDataY,1) - 0.5,np.size(imageDataY,0) - 0.5,-0.5))
                #self.ax[1].set_aspect(self.resolution[0]/self.resolution[2])
                self.slices = np.size(self.X,1)
                self.imy.set_data(imageDataY)
                self.ax[1].set_title('Slice cor %s/%s' % (self.y_ind,self.slices))
                self.imy.axes.figure.canvas.draw()
                
                imageDataZ = np.rot90(self.X[:,:,self.z_ind])
                if self.z_ind >= np.size(self.X,2):
                    self.z_ind = np.size(self.X,2)-1
                if self.axisLabels != None:
                    self.ax[2].set_xlabel(self.axisLabels[0])
                    self.ax[2].set_ylabel(self.axisLabels[1])
                self.imz.set_extent((-0.5,np.size(imageDataZ,1) - 0.5,np.size(imageDataZ,0) - 0.5,-0.5))
                #self.ax[2].set_aspect(self.resolution[0]/self.resolution[1])
                self.slices = np.size(self.X,2)
                self.imz.set_data(imageDataZ)
                self.ax[2].set_title('Slice ax %s/%s' % (self.z_ind,self.slices))
                self.imz.axes.figure.canvas.draw()   
            
        self.radio.on_clicked(planefunc)
        
        self.updateSlice()
            
    def keyPress(self, event):
        if event.key == "right":
            if self.index < 2 :
                self.index += 1
                if self.index == 1:
                    self.fig.suptitle("Action will be used on coronal axis", fontsize=14)
                elif self.index == 2:
                    self.fig.suptitle("Action will be used on axial axis", fontsize=14)
            elif self.index == 2:
                self.index = 0
                self.fig.suptitle("Action will be used on sagittal axis", fontsize=14)
        self.updateSlice()

    def onScroll(self, event):
        if event.button == 'up':
            if self.index == 0:
                self.x_ind = (self.x_ind + 1) % self.slices
            elif self.index == 1:
                self.y_ind = (self.y_ind + 1) % self.slices
            elif self.index == 2:
                self.z_ind = (self.z_ind + 1) % self.slices
        else:
            if self.index == 0:
                self.x_ind = (self.x_ind - 1) % self.slices
            elif self.index == 1:
                self.y_ind = (self.y_ind - 1) % self.slices
            elif self.index == 2:
                self.z_ind = (self.z_ind - 1) % self.slices        
        self.updateSlice()

    def onClick(self, event):
        if event.button == 3:  
            if self.index == 0:
                self.y_ind = self.data_mat[1] - round(event.xdata)
                self.z_ind = self.data_mat[2] - round(event.ydata)
            elif self.index == 1:
                self.x_ind = round(event.xdata)
                self.z_ind = self.data_mat[2] - round(event.ydata)
            elif self.index == 2:
                self.x_ind = round(event.xdata)
                self.y_ind = self.data_mat[1] - round(event.ydata)
            self.updateAllSlice()
            
            
        elif event.button == 2: #change orientation on scroll wheel click
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
        
    def onMotion(self, event):
        if self.mouseClicked == None: #return        #if mouse isn't clicked ignore movement
            if self.plane == 'cross line':
                if event.xdata != None:
                    if event.ydata != None:

                        if self.index == 0:  

                            imageDataX = np.flip(np.rot90(self.X[self.x_ind,:,:]), axis = 1)
                            s = np.flip(np.rot90(self.X[self.x_ind,:,:]), axis = 1).shape
                            arr1 = np.ones(s[0])
                            arr2 = np.ones(s[1])
                            arr3 = np.zeros((s)) # This array is only zeros with a vertical and horizontal white lines
                            arr3[int(event.ydata),:] = arr2
                            arr3[:,int(event.xdata)] = arr1

                            if self.x_ind >= np.size(self.X,0):
                                self.x_ind = np.size(self.X,0)-1
                            if self.axisLabels != None:
                                self.ax[0].set_xlabel(self.axisLabels[1])
                                self.ax[0].set_ylabel(self.axisLabels[2])
                            self.imx.set_extent((-0.5,np.size(imageDataX,1) - 0.5,np.size(imageDataX,0) - 0.5,-0.5))
                            #self.ax[0].set_aspect(self.resolution[1]/self.resolution[2])
                            self.slices = np.size(self.X,0)
                            self.imx.set_data(imageDataX + 10000*arr3)
                            self.ax[0].set_title('Slice sag %s/%s' % (self.x_ind,self.slices))
                            self.imx.axes.figure.canvas.draw()
                            
                            self.y_ind = self.data_mat[1] - int(event.xdata)
                            imageDataY = np.rot90(self.X[:,self.y_ind,:])
                            if self.y_ind >= np.size(self.X,1):
                                self.y_ind = np.size(self.X,1)-1
                            if self.axisLabels != None:
                                self.ax[1].set_xlabel(self.axisLabels[0])
                                self.ax[1].set_ylabel(self.axisLabels[2])
                            self.imy.set_extent((-0.5,np.size(imageDataY,1) - 0.5,np.size(imageDataY,0) - 0.5,-0.5))
                            #self.ax[1].set_aspect(self.resolution[0]/self.resolution[2])
                            self.slices = np.size(self.X,1)
                            self.imy.set_data(imageDataY)
                            self.ax[1].set_title('Slice cor %s/%s' % (self.y_ind,self.slices))
                            self.imy.axes.figure.canvas.draw()
                            
                            self.z_ind = self.data_mat[2] - 1 - int(event.ydata)
                            imageDataZ = np.rot90(self.X[:,:,self.z_ind])
                            if self.z_ind >= np.size(self.X,2):
                                self.z_ind = np.size(self.X,2)-1
                            if self.axisLabels != None:
                                self.ax[2].set_xlabel(self.axisLabels[0])
                                self.ax[2].set_ylabel(self.axisLabels[1])
                            self.imz.set_extent((-0.5,np.size(imageDataZ,1) - 0.5,np.size(imageDataZ,0) - 0.5,-0.5))
                            #self.ax[2].set_aspect(self.resolution[0]/self.resolution[1])
                            self.slices = np.size(self.X,2)
                            self.imz.set_data(imageDataZ)
                            self.ax[2].set_title('Slice ax %s/%s' % (self.z_ind,self.slices))
                            self.imz.axes.figure.canvas.draw()
                            
                        elif self.index == 1:

                            imageDataY = np.rot90(self.X[:,self.y_ind,:])
                            s = np.rot90(self.X[:,self.y_ind,:]).shape
                            arr1 = np.ones(s[0])
                            arr2 = np.ones(s[1])
                            arr3 = np.zeros((s)) # This array is only zeros with a vertical and horizontal white lines
                            arr3[int(event.ydata), :] = arr2
                            arr3[:, int(event.xdata)] = arr1

                            if self.y_ind >= np.size(self.X,1):
                                self.y_ind = np.size(self.X,1)-1
                            if self.axisLabels != None:
                                self.ax[1].set_xlabel(self.axisLabels[0])
                                self.ax[1].set_ylabel(self.axisLabels[2])
                            self.imy.set_extent((-0.5,np.size(imageDataY,1) - 0.5,np.size(imageDataY,0) - 0.5,-0.5))
                            #self.ax[1].set_aspect(self.resolution[0]/self.resolution[2])
                            self.slices = np.size(self.X,1)
                            self.imy.set_data(imageDataY + 10000*arr3)
                            self.ax[1].set_title('Slice cor %s/%s' % (self.y_ind,self.slices))
                            self.imy.axes.figure.canvas.draw()
                            
                            self.x_ind = int(event.xdata)
                            imageDataX = np.flip(np.rot90(self.X[self.x_ind,:,:]), axis = 1)
                            if self.x_ind >= np.size(self.X,0):
                                self.x_ind = np.size(self.X,0)-1
                            if self.axisLabels != None:
                                self.ax[0].set_xlabel(self.axisLabels[1])
                                self.ax[0].set_ylabel(self.axisLabels[2])
                            self.imx.set_extent((-0.5,np.size(imageDataX,1) - 0.5,np.size(imageDataX,0) - 0.5,-0.5))
                            #self.ax[0].set_aspect(self.resolution[1]/self.resolution[2])
                            self.slices = np.size(self.X,0)
                            self.imx.set_data(imageDataX)
                            self.ax[0].set_title('Slice sag %s/%s' % (self.x_ind,self.slices))
                            self.imx.axes.figure.canvas.draw()
                            
                            self.z_ind = self.data_mat[2] - 1 - int(event.ydata)
                            imageDataZ = np.rot90(self.X[:,:,self.z_ind])
                            if self.z_ind >= np.size(self.X,2):
                                self.z_ind = np.size(self.X,2)-1
                            if self.axisLabels != None:
                                self.ax[2].set_xlabel(self.axisLabels[0])
                                self.ax[2].set_ylabel(self.axisLabels[1])
                            self.imz.set_extent((-0.5,np.size(imageDataZ,1) - 0.5,np.size(imageDataZ,0) - 0.5,-0.5))
                            #self.ax[2].set_aspect(self.resolution[0]/self.resolution[1])
                            self.slices = np.size(self.X,2)
                            self.imz.set_data(imageDataZ)
                            self.ax[2].set_title('Slice ax %s/%s' % (self.z_ind,self.slices))
                            self.imz.axes.figure.canvas.draw()
                            
                        elif self.index == 2:

                            imageDataZ = np.rot90(self.X[:,:,self.z_ind])
                            s = np.rot90(self.X[:,:,self.z_ind]).shape
                            arr1 = np.ones(s[0])
                            arr2 = np.ones(s[1])
                            arr3 = np.zeros((s)) # This array is only zeros with a vertical and horizontal white lines
                            arr3[int(event.ydata), :] = arr2
                            arr3[:, int(event.xdata)] = arr1

                            if self.z_ind >= np.size(self.X,2):
                                self.z_ind = np.size(self.X,2)-1
                            if self.axisLabels != None:
                                self.ax[2].set_xlabel(self.axisLabels[0])
                                self.ax[2].set_ylabel(self.axisLabels[1])
                            self.imz.set_extent((-0.5,np.size(imageDataZ,1) - 0.5,np.size(imageDataZ,0) - 0.5,-0.5))
                            #self.ax[2].set_aspect(self.resolution[0]/self.resolution[1])
                            self.slices = np.size(self.X,2)
                            self.imz.set_data(imageDataZ + 10000*arr3)
                            self.ax[2].set_title('Slice ax %s/%s' % (self.z_ind,self.slices))
                            self.imz.axes.figure.canvas.draw()
                            
                            self.x_ind = int(event.xdata)
                            imageDataX = np.flip(np.rot90(self.X[self.x_ind,:,:]), axis = 1)
                            if self.x_ind >= np.size(self.X,0):
                                self.x_ind = np.size(self.X, 0)-1
                            if self.axisLabels != None:
                                self.ax[0].set_xlabel(self.axisLabels[1])
                                self.ax[0].set_ylabel(self.axisLabels[2])
                            self.imx.set_extent((-0.5,np.size(imageDataX,1) - 0.5,np.size(imageDataX,0) - 0.5,-0.5))
                            #self.ax[0].set_aspect(self.resolution[1]/self.resolution[2])
                            self.slices = np.size(self.X,0)
                            self.imx.set_data(imageDataX)
                            self.ax[0].set_title('Slice sag %s/%s' % (self.x_ind,self.slices))
                            self.imx.axes.figure.canvas.draw()
                            
                            self.y_ind = self.data_mat[1] - 1 - int(event.ydata)
                            imageDataY = np.rot90(self.X[:,self.y_ind,:])
                            if self.y_ind >= np.size(self.X,1):
                                self.y_ind = np.size(self.X,1)-1
                            if self.axisLabels != None:
                                self.ax[1].set_xlabel(self.axisLabels[0])
                                self.ax[1].set_ylabel(self.axisLabels[2])
                            self.imy.set_extent((-0.5,np.size(imageDataY,1) - 0.5,np.size(imageDataY,0) - 0.5,-0.5))
                            #self.ax[1].set_aspect(self.resolution[0]/self.resolution[2])
                            self.slices = np.size(self.X,1)
                            self.imy.set_data(imageDataY)
                            self.ax[1].set_title('Slice cor %s/%s' % (self.y_ind,self.slices))
                            self.imy.axes.figure.canvas.draw()
                
        elif self.mouseClicked != None:

            dx = event.xdata - self.mouseClicked[0]
            dy = event.ydata - self.mouseClicked[1]

            normDx = dx/self.mouseClicked[0]
            normDy = dy/self.mouseClicked[1]

            self.tempCmapRange = self.cmapRange*(1+normDy)
            self.tempCmapCenter = self.cmapCenter*(1+normDx)

            if self.index == 0:
                self.imx.set_clim(self.tempCmapCenter-self.tempCmapRange/2, self.tempCmapCenter+self.tempCmapRange/2)
                self.imx.axes.figure.canvas.draw()
            elif self.index == 1:
                self.imy.set_clim(self.tempCmapCenter-self.tempCmapRange/2, self.tempCmapCenter+self.tempCmapRange/2)
                self.imy.axes.figure.canvas.draw()
            elif self.index == 2:
                self.imz.set_clim(self.tempCmapCenter-self.tempCmapRange/2, self.tempCmapCenter+self.tempCmapRange/2)
                self.imz.axes.figure.canvas.draw()
        
    def updateSlice(self):
        
        if self.index == 0:
            if self.x_ind >= np.size(self.X,0):
                self.x_ind = np.size(self.X,0)-1
            imageDataX = np.flip(np.rot90(self.X[self.x_ind,:,:]), axis = 1)
            if self.axisLabels != None:
                self.ax[0].set_xlabel(self.axisLabels[1])
                self.ax[0].set_ylabel(self.axisLabels[2])
            self.imx.set_extent((-0.5,np.size(imageDataX,1) - 0.5,np.size(imageDataX,0) - 0.5,-0.5))
            #self.ax[0].set_aspect(self.resolution[1]/self.resolution[2])
            self.slices = np.size(self.X,0)
            self.imx.set_data(imageDataX)
            self.ax[0].set_title('Slice sag %s/%s' % (self.x_ind,self.slices))
            self.imx.axes.figure.canvas.draw()
        elif self.index == 1:
            if self.y_ind >= np.size(self.X,1):
                self.y_ind = np.size(self.X,1)-1
            imageDataY = np.rot90(self.X[:,self.y_ind,:])
            if self.axisLabels != None:
                self.ax[1].set_xlabel(self.axisLabels[0])
                self.ax[1].set_ylabel(self.axisLabels[2])
            self.imy.set_extent((-0.5,np.size(imageDataY,1) - 0.5,np.size(imageDataY,0) - 0.5,-0.5))
            #self.ax[1].set_aspect(self.resolution[0]/self.resolution[2])
            self.slices = np.size(self.X,1)
            self.imy.set_data(imageDataY)
            self.ax[1].set_title('Slice cor %s/%s' % (self.y_ind,self.slices))
            self.imy.axes.figure.canvas.draw()
        elif self.index == 2:
            if self.z_ind >= np.size(self.X,2):
                self.z_ind = np.size(self.X,2)-1
            imageDataZ = np.rot90(self.X[:,:,self.z_ind])
            if self.axisLabels != None:
                self.ax[2].set_xlabel(self.axisLabels[0])
                self.ax[2].set_ylabel(self.axisLabels[1])
            self.imz.set_extent((-0.5,np.size(imageDataZ,1) - 0.5,np.size(imageDataZ,0) - 0.5,-0.5))
            #self.ax[2].set_aspect(self.resolution[0]/self.resolution[1])
            self.slices = np.size(self.X,2)
            self.imz.set_data(imageDataZ)        
            self.ax[2].set_title('Slice ax %s/%s' % (self.z_ind,self.slices))
            self.imz.axes.figure.canvas.draw()
            
    def updateAllSlice(self):
        
        if self.x_ind >= np.size(self.X,0):
            self.x_ind = np.size(self.X,0)-1
        imageDataX = np.flip(np.rot90(self.X[self.x_ind,:,:]), axis = 1)
        if self.axisLabels != None:
            self.ax[0].set_xlabel(self.axisLabels[1])
            self.ax[0].set_ylabel(self.axisLabels[2])
        self.imx.set_extent((-0.5,np.size(imageDataX,1) - 0.5,np.size(imageDataX,0) - 0.5,-0.5))
        #self.ax[0].set_aspect(self.resolution[1]/self.resolution[2])
        self.slices = np.size(self.X,0)
        self.imx.set_data(imageDataX)
        self.ax[0].set_title('Slice sag %s/%s' % (self.x_ind,self.slices))
        self.imx.axes.figure.canvas.draw()

        if self.y_ind >= np.size(self.X,1):
            self.y_ind = np.size(self.X,1)-1
        imageDataY = np.rot90(self.X[:,self.y_ind,:])
        if self.axisLabels != None:
            self.ax[1].set_xlabel(self.axisLabels[0])
            self.ax[1].set_ylabel(self.axisLabels[2])
        self.imy.set_extent((-0.5,np.size(imageDataY,1) - 0.5,np.size(imageDataY,0) - 0.5,-0.5))
        #self.ax[1].set_aspect(self.resolution[0]/self.resolution[2])
        self.slices = np.size(self.X,1)
        self.imy.set_data(imageDataY)
        self.ax[1].set_title('Slice cor %s/%s' % (self.y_ind,self.slices))
        self.imy.axes.figure.canvas.draw()

        if self.z_ind >= np.size(self.X,2):
            self.z_ind = np.size(self.X,2)-1
        imageDataZ = np.rot90(self.X[:,:,self.z_ind])
        if self.axisLabels != None:
            self.ax[2].set_xlabel(self.axisLabels[0])
            self.ax[2].set_ylabel(self.axisLabels[1])
        self.imz.set_extent((-0.5,np.size(imageDataZ,1) - 0.5,np.size(imageDataZ,0) - 0.5,-0.5))
        #self.ax[2].set_aspect(self.resolution[0]/self.resolution[1])
        self.slices = np.size(self.X,2)
        self.imz.set_data(imageDataZ)        
        self.ax[2].set_title('Slice ax %s/%s' % (self.z_ind,self.slices))
        self.imz.axes.figure.canvas.draw()
        


