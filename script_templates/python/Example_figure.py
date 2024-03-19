import numpy as np
import matplotlib.pyplot as plt
from sem2d_read_fault import sem2d_read_fault
from sem2d_read_specgrid import sem2d_read_specgrid

# Modify the path of the output files
datadir = '../simulations/'

grid = sem2d_read_specgrid(datadir)
data = sem2d_read_fault('Flt01', datadir)


plt.figure(1)

plt.subplot(2, 2, 1)
# Slip rate image
plt.imshow(data['v'], extent=[0, (data['nt'] - 1) * data['dt'], data['x'].min() / 1e3, data['x'].max() / 1e3],aspect='auto')
plt.xlabel('Time (s)')
plt.ylabel('Along strike distance (km)')
plt.colorbar().set_label('Slip rate (m/s)')

plt.subplot(2, 2, 2)
# Shear stress image
plt.imshow((data['st0'] + data['st']) / 1e6, extent=[0, (data['nt'] - 1) * data['dt'], data['x'].min() / 1e3, data['x'].max() / 1e3],aspect='auto')
plt.xlabel('Time (s)')
plt.ylabel('Along strike distance (km)')
plt.colorbar().set_label('Shear stress (MPa)')

plt.subplot(2, 2, 3)
p1 = int(data['nx'] / 4.0)
p2 = int(data['nx'] / 2.0)
p3 = int(3 * data['nx'] / 4.0)
plt.plot(np.arange(data['nt']) * data['dt'], data['v'][:,p1])
plt.plot(np.arange(data['nt']) * data['dt'], data['v'][:,p2])
plt.plot(np.arange(data['nt']) * data['dt'], data['v'][:,p3])
plt.xlabel('Time (s)')
plt.ylabel('Slip rate (m/s)')

plt.subplot(2, 2, 4)
plt.plot(np.arange(data['nt']) * data['dt'], (data['st'][:,p1] + data['st0'][p1]) / 1e6)
plt.plot(np.arange(data['nt']) * data['dt'], (data['st'][:,p2] + data['st0'][p2]) / 1e6)
plt.plot(np.arange(data['nt']) * data['dt'], (data['st'][:,p3] + data['st0'][p3]) / 1e6)
plt.xlabel('Time (s)')
plt.ylabel('Shear stress (MPa)')

plt.savefig('./figs/rupture_info.png')

# Check if TP is used
try:
    _ = data['P'].shape
    print('TP is used.')

    plt.figure(2)

    plt.subplot(2, 2, 1)
    # Pore pressure image
    plt.imshow(data['P']/1e6, extent=[0, (data['nt'] - 1) * data['dt'], data['x'].min() / 1e3, data['x'].max() / 1e3],aspect='auto')
    plt.xlabel('Time (s)')
    plt.ylabel('Along strike distance (km)')
    plt.colorbar().set_label('Pore pressure (MPa)')

    plt.subplot(2, 2, 2)
    # Temperature
    plt.imshow(data['T'], extent=[0, (data['nt'] - 1) * data['dt'], data['x'].min() / 1e3, data['x'].max() / 1e3],aspect='auto')
    plt.xlabel('Time (s)')
    plt.ylabel('Along strike distance (km)')
    plt.colorbar().set_label('Temperature (K)')

    plt.subplot(2, 2, 3)
    p1 = int(data['nx'] / 4.0)
    p2 = int(data['nx'] / 2.0)
    p3 = int(3 * data['nx'] / 4.0)
    plt.plot(np.arange(data['nt']) * data['dt'], data['P'][:,p1]/1e6)
    plt.plot(np.arange(data['nt']) * data['dt'], data['P'][:,p2]/1e6)
    plt.plot(np.arange(data['nt']) * data['dt'], data['P'][:,p3]/1e6)
    plt.xlabel('Time (s)')
    plt.ylabel('Pore pressure (MPa)')

    plt.subplot(2, 2, 4)
    plt.plot(np.arange(data['nt']) * data['dt'], data['T'][:,p1])
    plt.plot(np.arange(data['nt']) * data['dt'], data['T'][:,p2])
    plt.plot(np.arange(data['nt']) * data['dt'], data['T'][:,p3])
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
except:
    print('TP is not used.')

plt.show()
