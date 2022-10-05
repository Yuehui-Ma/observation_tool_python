# %%
from tkinter import Grid
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style, quantity_support
plt.style.use(astropy_mpl_style)
quantity_support()
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.table import Table
import random
# %%

def obs_lst(off_coord, source_coord, date_time, lstrg = None):
    obs = EarthLocation(lat=37.3906*u.deg, lon=97.7289*u.deg, height = 3200*u.m) # observatory location: PMO DLH obsetvatory
    lon = 97.7289*u.deg
    time = Time(date_time)  
    off = SkyCoord(off_coord[0], off_coord[1], unit = 'deg', frame = 'galactic')
    source = SkyCoord(source_coord[0], source_coord[1], unit = 'deg', frame = 'galactic')

    delta_time = np.linspace(0, 24, 25)*u.hour
    bj_time = time + delta_time
    lst_time = bj_time.sidereal_time('mean', longitude = lon)
    lst = lst_time.hour
    tickname=lst_time.to_string(sep=':', pad = True, precision = 0, fields=2)

    offaltaz = off.transform_to(AltAz(obstime=bj_time,location=obs))
    offalt = offaltaz.alt.degree
    offaz = offaltaz.az.degree

    sourcealtaz = source.transform_to(AltAz(obstime=bj_time,location=obs))
    sourcealt = sourcealtaz.alt.degree
    sourceaz = sourcealtaz.az.degree

    delta_alt = np.abs(sourcealt - offalt)
    delta_az = np.abs(sourceaz - offaz)

    if np.any(lstrg == None):
        good = (delta_alt<=1) & (sourcealt>=30) & (delta_az<=5) & (sourcealt<80)
    else:
        if lstrg[0]<lstrg[1]:
            good = (delta_alt<=1) & (sourcealt>=30) & (delta_az<=5) & (sourcealt<80) & ((lst_time.hour>=lstrg[0]) & (lst_time.hour<=lstrg[1]))
        elif lstrg[0]>lstrg[1]:
            good = (delta_alt<=1) & (sourcealt>=30) & (delta_az<=5) & (sourcealt<80) & ((lst_time.hour>=lstrg[0]) | ((lst_time.hour<=lstrg[1]) & (lst_time.hour>=0)))
    msk = np.hstack((good, [False, False]))
    msk = msk & np.roll(msk, 1) & np.roll(msk, 2)
    msk  = msk | np.roll(msk, -1) | np.roll(msk, -2)
    if np.sum(msk[:-2])>=3:
        return True, (tickname[msk[:-2]])[0], (tickname[msk[:-2]])[-1]
    else:
        return False, None, None 

# %% 设定
lr = [120, 130]   # Sky area for MWISP survey, the Galactic longitude range.
br = [4.5, 10]      # The Galactic latitude range for test the observable sky areas. 
off = 'off.txt'# input calatog of candidate off position
date_time = '2022-03-11 00:00:00' # The start date-time of the observational plan. Scripts will test the observations for the next 24 hours. 
LST_r = [2, 6]          # The LST time inteval, in which we test the observable sky areas using the given off positions. 
# The Default observable condition is that the source altitude should be higher than 30 degree and lower than 80 degree, the altitude difference between the source and off position should less than 1 degree, and their azimuthal difference should be less than 5 degree. 


l = np.linspace(lr[0], lr[1], np.array((lr[1]-lr[0])/0.5).astype(int)+1)
b = np.linspace(br[0], br[1], np.array((br[1]-br[0])/0.5).astype(int)+1)
map = np.zeros([len(b), len(l)], dtype = bool)
r = lambda:random.randint(0,255) 

fig = plt.figure()
plt.imshow(map, origin='lower', cmap='Blues', extent=[l[-1], l[0], b[0], b[-1]])
ax = plt.gca()
ax.set_xlim(l[-1], l[0])
ax.set_ylim(b[0], b[-1])
ax.set_xlabel(r'$\rm Galactic\ Longitude\ (^{\circ})$')
ax.set_ylabel(r'$\rm Galactic\ Latitude\ (^{\circ})$')

if isinstance(off, str):
    offcat = Table.read(off, format = 'ascii')
    for k in range(len(offcat['col1'].data)):
        color = '#%02X%02X%02X' % (r(), r(), r())
        off_coord = [offcat['col1'].data[k], offcat['col2'].data[k]]
        ax.plot(off_coord[0], off_coord[1], marker = 'o', color = color, markeredgecolor = 'black', zorder = 2, alpha = 0.7)
        offname = 'F'+'%07.3f' % off_coord[0] + '%+07.3f' % off_coord[1]
        f = open(offname+'.cat', 'w')
        for i in range(len(l)):
            for j in range(len(b)):
                source_coord = [l[i], b[j]]
                flag = obs_lst(off_coord, source_coord, date_time, lstrg=LST_r)
                if any(flag):
                    ax.plot(source_coord[0], source_coord[1], marker = 's', color = color)
                    cld = '%04i' % (source_coord[0]*10) + '%+04i' % (source_coord[1]*10)
                    f.write(cld +','+flag[1]+'->'+flag[2]+'\n') 
        f.close()                
else: 
    off_coord = off.copy()
    color = '#%02X%02X%02X' % (r(), r(), r())
    ax.plot(off_coord[0], off_coord[1], marker = 'o', color = color, markeredgecolor = 'black', zorder = 2, alpha = 0.7)
    offname = 'F'+'%07.3f' % off_coord[0] + '%+07.3f' % off_coord[1]
    f = open(offname+'.cat', 'w')
    for i in range(len(l)):
        for j in range(len(b)):
            source_coord = [l[i], b[j]]
            flag = obs_lst(off_coord, source_coord, date_time, lstrg=LST_r)
            if any(flag):
                ax.plot(source_coord[0], source_coord[1], marker = 's', color = color)
                cld = '%04i' % (source_coord[0]*10) + '%+04i' % (source_coord[1]*10)
                f.write(cld +','+flag[1]+'->'+flag[2]+'\n') 
    f.close()
fig.savefig('G'+str(lr[0])+'_obs.pdf', bbox_inches = 'tight')
