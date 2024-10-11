import numpy as np
from scipy import constants
import matplotlib.pyplot as plt

# Functions of the sequences
def spin_echo_seq(TR, TE, t1, t2, m0):
    # Spin echo formula
    TR_div = np.divide(TR, t1, out=np.zeros_like(t1), where=t1!=0)
    TE_div = np.divide(TE, t2, out=np.zeros_like(t2), where=t2!=0)
    return np.abs(m0 * (1 - np.exp(-TR_div)) * np.exp(-TE_div))   

def Gradient_seq(TR, TE, t1, t2_star, m0, alp):
    # Gradient echo formula
    TR_div = np.divide(-TR, t1, out=np.zeros_like(t1), where=t1!=0)
    TE_div = np.divide(-TE, t2_star, out=np.zeros_like(t2_star), where=t2_star!=0)
    a = np.sin(alp) * (1 - np.exp(TR_div)) * np.exp(TE_div)
    b = 1 - np.cos(alp) * np.exp(TR_div)
    im = np.divide(a, b, out=np.zeros_like(a), where=a!=0)
    return np.abs( m0 * im )

def IN_seq(TR, TI, t1, t2, m0):
    # Inversion recovery formula
    a = 2 * np.exp(-np.divide(TI, t1, out=np.zeros_like(t1), where=t1!=0))
    b = np.exp(-np.divide(TR, t1, out=np.zeros_like(t1), where=t1!=0))
    c = 1-a+b
    return np.abs(m0 * c)

def DoubleInversion_seq(TR, TE, TI1, TI2, t1, t2, m0):
    # Double inversion recovery formula 
    E1 = np.exp(-np.divide(TI1, t1, out=np.zeros_like(t1), where=t1!=0))
    E2 = np.exp(-np.divide(TI2, t1, out=np.zeros_like(t1), where=t1!=0))
    Ec = np.exp(-np.divide(TR, t1, out=np.zeros_like(t1), where=t1!=0))
    Etau = np.exp(np.divide((TE/2), t1, out=np.zeros_like(t1), where=t1!=0))
    Mz = 1 - 2*E2 + 2*E1*E2 - Ec*(2*Etau - 1)
    return np.abs(m0 * Mz)

def Diffusion_seq(TR, t1, t2, m0, b, D):
    # Diffusion formula with predefine b 
    TR_div = np.divide(TR, t1, out=np.zeros_like(t1), where=t1!=0)
    spin = np.abs(m0 * (1 - np.exp(-TR_div)))
    return np.abs(spin * np.exp(-b*D))    

def SSFP_Echo_seq(t1, t2, M0, alpha, phi):
    # SSFP formula
    a = np.abs(np.sin(alpha))*np.sqrt(2 - 2*np.cos(phi))
    b = (1+np.cos(alpha))*(1-np.cos(alpha)) + 2*(1-np.cos(alpha))*np.divide(t1, t2, out=np.zeros_like(t2), where=t2!=0)
    return M0 * (a/b)

def SSFP_new(TR, t1, t2, M0, alpha, phi):
    # SSFP formula
    E1 = np.divide(-TR, t1, out=np.zeros_like(t1), where=t1!=0)
    E2 = np.divide(-TR, t2, out=np.zeros_like(t2), where=t2!=0)
    #E1 = np.divide(-t1, TR, out=np.zeros_like(t1))
    #E2 = np.divide(-t2, TR, out=np.zeros_like(t2))
    #a = (1-E1) * (1-E2*(np.exp(-1j*phi).real)) * np.sin(alpha)
    a = (1-E1) * (1-E2*np.exp(-phi)) * np.sin(alpha)
    C = E2 * (E1-1) * (1 + np.cos(alpha))
    D = (1-E1*np.cos(alpha)) - (E1 - np.cos(alpha))*(E2**2)
    b = C * np.cos(phi) + D
    return M0 * (a/b)

def TSE_seq_R(TE, ETL, t2, met, read, data):
    
    Data_mat_ETL = np.divide(data.shape,ETL).astype(int)
    Data_mat_ETL[Data_mat_ETL == 0] = 1

    T1 = 0
    T2 = 0
    
    t2_flat = t2.flatten()[None]
    ma = t2_flat[t2_flat != 0]
    
    if read == 'LR':               
                
        # number of points in each dimensions
        T1 = data.shape[1] 
        T2 = data.shape[2]
    
        # Time array is just a linear line x=y
        time1 = np.linspace(0,TE*ETL,T1)[:, None] 
        time2 = np.linspace(0,TE*ETL,T2)[:, None]
        sig1 = np.zeros(T1)
        sig2 = np.zeros(T2)
        
        T2decay1 = np.exp((-time1)/ma)
        sig1 = np.sum(T2decay1, 1)  # Normal t2 decay
        del T2decay1
        T2decay2 = np.exp((-time2)/ma)
        sig2 = np.sum(T2decay2, 1)        
        del T2decay2
        
        # kspace trajectories
        linear1 = np.zeros((sig1.shape))
        inout1 = np.zeros((sig1.shape))
        outin1 = np.zeros((sig1.shape))
        linear2 = np.zeros((sig2.shape))
        inout2 = np.zeros((sig2.shape))
        outin2 = np.zeros((sig2.shape))

        even1 = np.arange(0,sig1.shape[0],2)
        odd1 = np.arange(1,sig1.shape[0],2)
        flip1 = np.flip(even1)
        even2 = np.arange(0,sig2.shape[0],2)
        odd2 = np.arange(1,sig2.shape[0],2)
        flip2 = np.flip(even2)

        arr1 = np.concatenate((flip1,odd1))
        arr21 = np.concatenate((odd1,flip1))
        arr2 = np.concatenate((flip2,odd2))
        arr22 = np.concatenate((odd2,flip2))

        for i in range(len(sig1)):
            linear1[i] = sig1[i]
            inout1[i] = sig1[arr1[i]]
            outin1[i] = sig1[arr21[i]]
        for i in range(len(sig2)):
            linear2[i] = sig2[i]
            inout2[i] = sig2[arr2[i]]
            outin2[i] = sig2[arr22[i]]
              
        # Making the steps
        outin_step1 = np.zeros((outin1.shape))
        inout_step1 = np.zeros((inout1.shape))
        linear_step1 = np.zeros((linear1.shape))
        outin_step2 = np.zeros((outin2.shape))
        inout_step2 = np.zeros((inout2.shape))
        linear_step2 = np.zeros((linear2.shape))
        
        for i in range(ETL):
            s = Data_mat_ETL[1]
            outin_step1[i*s:(i+1)*s] = outin1[i*s]
            inout_step1[i*s:(i+1)*s] = inout1[i*s]
            linear_step1[i*s:(i+1)*s] = linear1[i*s]
        for i in range(ETL):
            s = Data_mat_ETL[2]
            outin_step2[i*s:(i+1)*s] = outin2[i*s]
            inout_step2[i*s:(i+1)*s] = inout2[i*s]
            linear_step2[i*s:(i+1)*s] = linear2[i*s]
             
        # Fourier transform of trajectory
        FT1 = 0
        ft1_coef = 0
        FT2 = 0
        ft2_coef = 0
        if met == 'Out-in':    
            FT1 = np.fft.fftshift(np.fft.fft(outin_step1))  
            FT2 = np.fft.fftshift(np.fft.fft(outin_step2))  
        elif met == 'In-out':
            FT1 = np.fft.fftshift(np.fft.fft(inout_step1))
            FT2 = np.fft.fftshift(np.fft.fft(inout_step2))
        elif met == 'Linear':
            FT1 = np.fft.fftshift(np.fft.fft(linear_step1))
            FT2 = np.fft.fftshift(np.fft.fft(linear_step2))
        
        FT1 = np.abs(FT1/np.max(FT1))
        FT2 = np.abs(FT2/np.max(FT2))
        
        ft1_coef = np.sum(FT1)
        ft2_coef = np.sum(FT2)

        im = np.zeros((data.shape))
            
        for j in range(data.shape[0]):
            for i in range(data.shape[1]): # Convolution of each row with the lorentzian
                im[j,i,:] = np.convolve(FT2.real, data[j,i,:], mode='same')

        for j in range(data.shape[0]):
            for i in range(data.shape[2]): # Convolution of each row with the lorentzian
                im[j,:,i] = np.convolve(FT1.real, im[j,:,i], mode='same')
                
        im = np.divide(im,ft1_coef)
        im = np.divide(im,ft2_coef)        
        
    elif read == 'AP':
     
        # number of points in each dimensions
        T1 = data.shape[2] 
        T2 = data.shape[0]
    
        time1 = np.linspace(0,TE*ETL,T1)[:, None] # Time array is just a linear line x=y
        time2 = np.linspace(0,TE*ETL,T2)[:, None]
        sig1 = np.zeros(T1)
        sig2 = np.zeros(T2)
        
        T2decay1 = np.exp((-time1)/ma)
        sig1 = np.sum(T2decay1, 1)  # Normal t2 decay
        del T2decay1
        T2decay2 = np.exp((-time2)/ma)
        sig2 = np.sum(T2decay2, 1)        
        del T2decay2       
    
        # kspace trajectories
        linear1 = np.zeros((sig1.shape))
        inout1 = np.zeros((sig1.shape))
        outin1 = np.zeros((sig1.shape))
        linear2 = np.zeros((sig2.shape))
        inout2 = np.zeros((sig2.shape))
        outin2 = np.zeros((sig2.shape))

        even1 = np.arange(0,sig1.shape[0],2)
        odd1 = np.arange(1,sig1.shape[0],2)
        flip1 = np.flip(even1)
        even2 = np.arange(0,sig2.shape[0],2)
        odd2 = np.arange(1,sig2.shape[0],2)
        flip2 = np.flip(even2)

        arr1 = np.concatenate((flip1,odd1))
        arr21 = np.concatenate((odd1,flip1))
        arr2 = np.concatenate((flip2,odd2))
        arr22 = np.concatenate((odd2,flip2))

        for i in range(len(sig1)):
            linear1[i] = sig1[i]
            inout1[i] = sig1[arr1[i]]
            outin1[i] = sig1[arr21[i]]
        for i in range(len(sig2)):
            linear2[i] = sig2[i]
            inout2[i] = sig2[arr2[i]]
            outin2[i] = sig2[arr22[i]]
              
        # Making the steps
        outin_step1 = np.zeros((outin1.shape))
        inout_step1 = np.zeros((inout1.shape))
        linear_step1 = np.zeros((linear1.shape))
        outin_step2 = np.zeros((outin2.shape))
        inout_step2 = np.zeros((inout2.shape))
        linear_step2 = np.zeros((linear2.shape))
        
        for i in range(ETL):
            s = Data_mat_ETL[2]
            outin_step1[i*s:(i+1)*s] = outin1[i*s]
            inout_step1[i*s:(i+1)*s] = inout1[i*s]
            linear_step1[i*s:(i+1)*s] = linear1[i*s]
        for i in range(ETL):
            s = Data_mat_ETL[0]
            outin_step2[i*s:(i+1)*s] = outin2[i*s]
            inout_step2[i*s:(i+1)*s] = inout2[i*s]
            linear_step2[i*s:(i+1)*s] = linear2[i*s]
             
        # Fourier transform of trajectory
        FT1 = 0
        ft1_coef = 0
        FT2 = 0
        ft2_coef = 0
        if met == 'Out-in':    
            FT1 = np.fft.fftshift(np.fft.fft(outin_step1))  
            FT2 = np.fft.fftshift(np.fft.fft(outin_step2))  
        elif met == 'In-out':
            FT1 = np.fft.fftshift(np.fft.fft(inout_step1))
            FT2 = np.fft.fftshift(np.fft.fft(inout_step2))
        elif met == 'Linear':
            FT1 = np.fft.fftshift(np.fft.fft(linear_step1))
            FT2 = np.fft.fftshift(np.fft.fft(linear_step2))
        
        FT1 = np.abs(FT1/np.max(FT1))
        FT2 = np.abs(FT2/np.max(FT2))
        
        ft1_coef = np.sum(FT1)
        ft2_coef = np.sum(FT2)

        im = np.zeros((data.shape))
            
        for j in range(data.shape[1]):
            for i in range(data.shape[2]): # Convolution of each row with the lorentzian
                im[:,j,i] = np.convolve(FT2.real, data[:,j,i], mode='same')

        for j in range(data.shape[1]):
            for i in range(data.shape[0]): # Convolution of each row with the lorentzian
                im[i,j,:] = np.convolve(FT1.real, im[i,j,:], mode='same')
                
        im = np.divide(im,ft1_coef)
        im = np.divide(im,ft2_coef)
         
    elif read == 'FH':
        
        # number of points in each dimensions
        T1 = data.shape[1] 
        T2 = data.shape[0]
    
        time1 = np.linspace(0,TE*ETL,T1)[:, None] # Time array is just a linear line x=y
        time2 = np.linspace(0,TE*ETL,T2)[:, None]
        sig1 = np.zeros(T1)
        sig2 = np.zeros(T2)
        
        T2decay1 = np.exp((-time1)/ma)
        sig1 = np.sum(T2decay1, 1)  # Normal t2 decay
        del T2decay1
        T2decay2 = np.exp((-time2)/ma)
        sig2 = np.sum(T2decay2, 1)        
        del T2decay2    
    
        # kspace trajectories
        linear1 = np.zeros((sig1.shape))
        inout1 = np.zeros((sig1.shape))
        outin1 = np.zeros((sig1.shape))
        linear2 = np.zeros((sig2.shape))
        inout2 = np.zeros((sig2.shape))
        outin2 = np.zeros((sig2.shape))

        even1 = np.arange(0,sig1.shape[0],2)
        odd1 = np.arange(1,sig1.shape[0],2)
        flip1 = np.flip(even1)
        even2 = np.arange(0,sig2.shape[0],2)
        odd2 = np.arange(1,sig2.shape[0],2)
        flip2 = np.flip(even2)

        arr1 = np.concatenate((flip1,odd1))
        arr21 = np.concatenate((odd1,flip1))
        arr2 = np.concatenate((flip2,odd2))
        arr22 = np.concatenate((odd2,flip2))

        for i in range(len(sig1)):
            linear1[i] = sig1[i]
            inout1[i] = sig1[arr1[i]]
            outin1[i] = sig1[arr21[i]]
        for i in range(len(sig2)):
            linear2[i] = sig2[i]
            inout2[i] = sig2[arr2[i]]
            outin2[i] = sig2[arr22[i]]
              
        # Making the steps
        outin_step1 = np.zeros((outin1.shape))
        inout_step1 = np.zeros((inout1.shape))
        linear_step1 = np.zeros((linear1.shape))
        outin_step2 = np.zeros((outin2.shape))
        inout_step2 = np.zeros((inout2.shape))
        linear_step2 = np.zeros((linear2.shape))
        
        for i in range(ETL):
            s = Data_mat_ETL[1]
            outin_step1[i*s:(i+1)*s] = outin1[i*s]
            inout_step1[i*s:(i+1)*s] = inout1[i*s]
            linear_step1[i*s:(i+1)*s] = linear1[i*s]
        for i in range(ETL):
            s = Data_mat_ETL[0]
            outin_step2[i*s:(i+1)*s] = outin2[i*s]
            inout_step2[i*s:(i+1)*s] = inout2[i*s]
            linear_step2[i*s:(i+1)*s] = linear2[i*s]
             
        # Fourier transform of trajectory
        FT1 = 0
        ft1_coef = 0
        FT2 = 0
        ft2_coef = 0
        if met == 'Out-in':    
            FT1 = np.fft.fftshift(np.fft.fft(outin_step1))  
            FT2 = np.fft.fftshift(np.fft.fft(outin_step2))  
        elif met == 'In-out':
            FT1 = np.fft.fftshift(np.fft.fft(inout_step1))
            FT2 = np.fft.fftshift(np.fft.fft(inout_step2))
        elif met == 'Linear':
            FT1 = np.fft.fftshift(np.fft.fft(linear_step1))
            FT2 = np.fft.fftshift(np.fft.fft(linear_step2))
        
        FT1 = np.abs(FT1/np.max(FT1))
        FT2 = np.abs(FT2/np.max(FT2))
        
        ft1_coef = np.sum(FT1)
        ft2_coef = np.sum(FT2)

        im = np.zeros((data.shape))
            
        for j in range(data.shape[2]):
            for i in range(data.shape[0]): # Convolution of each row with the lorentzian
                im[i,:,j] = np.convolve(FT1.real, data[i,:,j], mode='same')

        for j in range(data.shape[2]):
            for i in range(data.shape[1]): # Convolution of each row with the lorentzian
                im[:,i,j] = np.convolve(FT2.real, im[:,i,j], mode='same')
                
        im = np.divide(im,ft1_coef)
        im = np.divide(im,ft2_coef)
    
    return np.abs(im)
    