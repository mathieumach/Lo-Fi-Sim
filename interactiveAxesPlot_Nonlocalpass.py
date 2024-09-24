import numpy as np
import matplotlib as plt
from skimage.restoration import denoise_nl_means, estimate_sigma # estimate_sigma --> wavelet-based estimator of the (Gaussian) noise standard deviation.
from matplotlib.widgets import Slider, Button, RadioButtons

class interactivePlot(object):
    def __init__(self,fig, ax, X, h, s, d, data_mat, plotAxis = 2, axisLabels = None, fov = None, cmap = 'gray'):
        
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
        self.axisLabels = axisLabels
        self.colorMap = cmap
        self.original = np.copy(X)
        self.index = 0
        self.h = h
        self.s = s
        self.d = d
        self.data_mat = data_mat
        
        # estimate the noise standard deviation from the noisy image
        sigma_estim = estimate_sigma(self.original, channel_axis=None)
        self.sigma_estim = sigma_estim
        
        if fov is None:
            self.fov = (1,1,1)
            self.resolution = (1,1,1)
        else:
            self.fov = fov
            self.resolution = np.divide(fov, np.shape(X))

        self.slices = np.shape(X)[plotAxis]
        self.ind = self.slices//2
        self.x_ind = 0
        self.y_ind = 1
        self.z_ind = 2
        self.filterimages = np.copy(X)
        
        # To change the axis value to mm (the fov)
        ax_pos = [[0, data_mat[0]/2, data_mat[0]-1], [0, data_mat[1]/2, data_mat[1]-1], [0, data_mat[2]/2, data_mat[2]-1]]
        ax_label = [[0, self.fov[0]/2, self.fov[0]], [0, self.fov[1]/2, self.fov[1]], [0, self.fov[2]/2, self.fov[2]]]
        
        if self.plotAxis == 0:
            imageDatax = self.X[self.x_ind]
            denoise_nl_means(self.filterimages[self.x_ind], h=h * sigma_estim, sigma=sigma_estim, fast_mode=False, patch_size=s, patch_distance=d)
        elif self.plotAxis == 1:
            imageDatax = self.X[self.y_ind]
            denoise_nl_means(self.filterimages[self.y_ind], h=h * sigma_estim, sigma=sigma_estim, fast_mode=False, patch_size=s, patch_distance=d)
        elif self.plotAxis == 2:
            imageDatax = self.X[self.z_ind]
            imageDatay = denoise_nl_means(self.filterimages[self.z_ind], h=h * sigma_estim, sigma=sigma_estim, fast_mode=False, patch_size=s, patch_distance=d)
        else:
            print("invalid axis")
            return -1
        self.imx = ax[0].imshow(np.rot90(imageDatax), cmap = self.colorMap,vmin = self.cmapCenter-self.cmapRange/2, vmax= self.cmapCenter+self.cmapRange/2)
        self.imy = ax[1].imshow(np.rot90(imageDatay), cmap = self.colorMap,vmin = self.cmapCenter-self.cmapRange/2, vmax= self.cmapCenter+self.cmapRange/2)

        # Global title
        self.fig.suptitle("Axial axis", fontsize=14)
        
        #create sliders
        # H
        self.axSlider_H = fig.add_axes([0.15, 0.1, 0.6, 0.03])
        self.slider_H = Slider(ax=self.axSlider_H, label='h', valmin=0, valmax=3, valinit=self.h)
        # s
        allowed_s = np.concatenate([np.linspace(0, 20, 21)])
        self.axSlider_s = fig.add_axes([0.15, 0.06, 0.6, 0.03])
        self.slider_s = Slider(ax=self.axSlider_s, label='s', valmin=0, valstep = allowed_s, valmax = allowed_s[-1], valinit=self.s) 
        # d
        allowed_d = np.concatenate([np.linspace(0, 15, 16)])
        self.axSlider_d = fig.add_axes([0.15, 0.02, 0.6, 0.03])
        self.slider_d = Slider(ax=self.axSlider_d, label='d', valmin=0, valstep = allowed_d, valmax = allowed_d[-1], valinit=self.d) 
        
        def update(val):
            # H
            if self.slider_H.val == 0:
                hh = self.H
            else:
                hh = self.slider_H.val
            # s
            if self.slider_s.val == 0:
                ss = self.s
            else:
                ss = self.slider_s.val
            # d
            if self.slider_d.val == 0:
                dd = self.alpha_d
            else:
                dd = self.slider_d.val

            ss = int(ss)
            dd = int(dd)
            
            if self.plotAxis == 0:
                self.filterimages[self.x_ind] = denoise_nl_means(self.original[self.x_ind], h=hh * sigma_estim, sigma=sigma_estim, fast_mode=False, patch_size=ss, patch_distance=dd)
            if self.plotAxis == 1:
                self.filterimages[self.y_ind] = denoise_nl_means(self.original[self.y_ind], h=hh * sigma_estim, sigma=sigma_estim, fast_mode=False, patch_size=ss, patch_distance=dd)
            elif self.plotAxis == 2:
                self.filterimages[self.z_ind] = denoise_nl_means(self.original[self.z_ind], h=hh * sigma_estim, sigma=sigma_estim, fast_mode=False, patch_size=ss, patch_distance=dd)
            
            self.updateAllSlice()
        
        self.slider_H.on_changed(update)
        self.slider_s.on_changed(update)
        self.slider_d.on_changed(update)
            
    def keyPress(self, event):
        hh = self.slider_H.val
        ss = int(self.slider_s.val)
        dd = int(self.slider_d.val)
        if event.key == " ": #change orientation on space bar press
            self.plotAxis += 1
            if self.plotAxis == 1:
                self.fig.suptitle("Coronal axis", fontsize=14)
                self.filterimages[self.y_ind] = denoise_nl_means(self.original[self.y_ind], h=hh * self.sigma_estim, sigma=self.sigma_estim, fast_mode=False, patch_size=ss, patch_distance=dd)
            elif self.plotAxis == 2:
                self.fig.suptitle("Axial axis", fontsize=14)
                self.filterimages[self.z_ind] = denoise_nl_means(self.original[self.z_ind], h=hh * self.sigma_estim, sigma=self.sigma_estim, fast_mode=False, patch_size=ss, patch_distance=dd)
            if self.plotAxis > 2:
                self.plotAxis = 0
                self.fig.suptitle("Sagittal axis", fontsize=14)
                self.filterimages[self.x_ind] = denoise_nl_means(self.original[self.x_ind], h=hh * self.sigma_estim, sigma=self.sigma_estim, fast_mode=False, patch_size=ss, patch_distance=dd)
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
        
    def onMotion(self, event):
        if self.mouseClicked == None: return        #if mouse isn't clicked ignore movement

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
        
    def updateSlice(self):
        if self.plotAxis == 0:
            imageDataX = np.flip(np.rot90(self.original[self.x_ind]), axis = 1)
            if self.axisLabels != None:
                self.ax.set_xlabel(self.axisLabels[2])
                self.ax.set_ylabel(self.axisLabels[1])            
            imageDataY = np.flip(np.rot90(self.filterimages[self.x_ind]), axis = 1)
            if self.axisLabels != None:
                self.ax.set_xlabel(self.axisLabels[2])
                self.ax.set_ylabel(self.axisLabels[1])
            
        elif self.plotAxis == 1:
            imageDataX = np.rot90(self.original[self.y_ind])
            if self.axisLabels != None:
                self.ax.set_xlabel(self.axisLabels[2])
                self.ax.set_ylabel(self.axisLabels[0])
            imageDataY = np.rot90(self.filterimages[self.y_ind])
            if self.axisLabels != None:
                self.ax.set_xlabel(self.axisLabels[2])
                self.ax.set_ylabel(self.axisLabels[0])
            
        else:
            imageDataX = np.rot90(self.original[self.z_ind])
            if self.axisLabels != None:
                self.ax.set_xlabel(self.axisLabels[1])
                self.ax.set_ylabel(self.axisLabels[0])
            imageDataY = np.rot90(self.filterimages[self.z_ind])
            if self.axisLabels != None:
                self.ax.set_xlabel(self.axisLabels[1])
                self.ax.set_ylabel(self.axisLabels[0])

        self.imx.set_extent((-0.5,np.size(imageDataX,1) - 0.5,np.size(imageDataX,0) - 0.5,-0.5))
        self.imx.set_data(imageDataX)
        self.imx.axes.figure.canvas.draw()  
        
        self.imy.set_extent((-0.5,np.size(imageDataY,1) - 0.5,np.size(imageDataY,0) - 0.5,-0.5))
        self.imy.set_data(imageDataY)
        self.imy.axes.figure.canvas.draw()
        
    def updateAllSlice(self):
        if self.plotAxis == 0:
            
            imageDataX = np.flip(np.rot90(self.original[self.x_ind]), axis = 1)
            if self.axisLabels != None:
                self.ax[0].set_xlabel(self.axisLabels[1])
                self.ax[0].set_ylabel(self.axisLabels[2])

            imageDataY = np.flip(np.rot90(self.filterimages[self.x_ind]), axis = 1)
            if self.axisLabels != None:
                self.ax[1].set_xlabel(self.axisLabels[1])
                self.ax[1].set_ylabel(self.axisLabels[2])
            
        elif self.plotAxis == 1:
            
            imageDataX = np.rot90(self.original[self.y_ind])
            if self.axisLabels != None:
                self.ax[0].set_xlabel(self.axisLabels[0])
                self.ax[0].set_ylabel(self.axisLabels[2])

            imageDataY = np.rot90(self.filterimages[self.y_ind])
            if self.axisLabels != None:
                self.ax[1].set_xlabel(self.axisLabels[0])
                self.ax[1].set_ylabel(self.axisLabels[2])
            
        elif self.plotAxis == 2:
    
            imageDataX = np.rot90(self.original[self.z_ind])
            if self.axisLabels != None:
                self.ax[0].set_xlabel(self.axisLabels[0])
                self.ax[0].set_ylabel(self.axisLabels[1])

            imageDataY = np.rot90(self.filterimages[self.z_ind])
            if self.axisLabels != None:
                self.ax[1].set_xlabel(self.axisLabels[0])
                self.ax[1].set_ylabel(self.axisLabels[1])
                
        self.imx.set_data(imageDataX)
        self.imx.axes.figure.canvas.draw()
        self.imy.set_data(imageDataY)
        self.imy.axes.figure.canvas.draw()




