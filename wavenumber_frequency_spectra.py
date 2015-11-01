def spec_est3(A,d1,d2,d3):

    l1,l2,l3 = A.shape


    df1 = 1./(l1*d1)
    df2 = 1./(l2*d2)
    df3 = 1./(l3*d3)
    f1Ny = 1./(2*d1)
    f2Ny = 1./(2*d2)
    f3Ny = 1./(2*d3)

    f1 = np.arange(-f1Ny,f1Ny,df1)
    f2 = np.arange(-f2Ny,f2Ny,df2)
    f3 = np.arange(0,l3/2+1)*df3

    # spectral window

    # first, the spatial window
    wx = np.matrix(np.hanning(l1))
    wy =  np.matrix(np.hanning(l2))
    window_s = np.repeat(np.array(wx.T*wy),l3).reshape(l1,l2,l3)

    # now, the time window
    wt = np.hanning(l3)
    window_t = np.repeat(wt,l1*l2).reshape(l3,l2,l1).T

    Ahat = np.fft.rfftn(window_s*window_t*A)
    
    Aabs = 2 * (Ahat*Ahat.conjugate()) / (df1*df2*df3) / ((l1*l2*l3)**2)

    return np.fft.fftshift(Aabs.real,axes=(0,1)),f1,f2,f3

if __name__=='__main__':
    
    import matplotlib.pyplot as plt
    import numpy as np
    import scipy.signal
    import scipy as sp
    import glob, os 
    import seawater as sw
    from netCDF4 import Dataset
    import sys

    plt.close('all')
    plt.rcParams.update({'font.size': 24})

    iz = 100     # vertical level [m]
    #data_path = 'uv_100m_single_file_small.nc' 
    #data_path = "uv_100m_single_file_small_daily_averaged.nc"
    data_path = "uv_0m_single_file_small_anomalies.nc"
    grid_path = './grid.npz' 

    UV = Dataset(data_path)
    grid = np.load(grid_path)

    # wavenumber
    ix,iy,it = UV.variables['u'][:,:,:].shape

    d1, d2, d3 = 1.19,1.17,1.

    iaux = 24*7
    nt = it/(24*7)
    for i in range(nt):
        uaux = UV.variables['u'][:,:,i*iaux:i*iaux+iaux]
        vaux = UV.variables['v'][:,:,i*iaux:i*iaux+iaux]

        if i == 0:
            Eu,k,l,omg = spec_est3(uaux,d1,d2,d3)
            Ev,k,l,omg = spec_est3(vaux,d1,d2,d3)
        else:
            Eua,_,_,_ = spec_est3(uaux,d1,d2,d3)
            Eva,_,_,_ = spec_est3(vaux,d1,d2,d3)

            Eu = Eu + Eua
            Ev = Ev + Eva

    Eu = Eu/nt
    Ev = Ev/nt

    np.savez("wavenumber_frequency_spec.npz", E=(Eu+Ev)/2.,k=k,l=l,omg=omg)

