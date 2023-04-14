#!/usr/bin/env python
"""
Created on Thu Apr 13 21:23:25 2023

@author: fanhao

 KPOINTS, OUTCAR, PROCAR are needed
"""
import numpy as np
import re
import os
import matplotlib.pyplot as plt
import matplotlib as mpl

############################################################
def get_bandInfo(fname = 'OUTCAR'):
    'This part is copied from pyband: https://github.com/QijingZheng/pyband'
    outcar = [line for line in open(fname) if line.strip()]
    for ii, line in enumerate(outcar):
        if 'NKPTS =' in line:
            nkpts = int(line.split()[3])
            nband = int(line.split()[-1])
        if 'ISPIN  =' in line:
            ispin = int(line.split()[2])
        if "k-points in reciprocal lattice and weights" in line:
            Lvkpts = ii + 1
        if 'reciprocal lattice vectors' in line:
            ibasis = ii + 1
        if 'E-fermi' in line:
            Efermi = float(line.split()[2])
            LineEfermi = ii + 3
            # break
    B = np.array([line.split()[-3:] for line in outcar[ibasis:ibasis+3]], dtype=float)
    # k-points vectors and weights
    tmp = np.array([line.split() for line in outcar[Lvkpts:Lvkpts+nkpts]], dtype=float)
    vkpts = tmp[:,:3]
    wkpts = tmp[:,-1]
    # for ispin = 2, there are two extra lines "spin component..."
    N = (nband + 2) * nkpts * ispin + (ispin - 1) * 2
    bands = []
    # vkpts = []
    for line in outcar[LineEfermi:LineEfermi + N-1]:
        if 'spin component' in line or 'band No.' in line:
            continue
        if 'k-point' in line:
            # vkpts += [line.split()[3:]]
            continue
        #print (line)
        bands.append(float(line.split()[1]))

    bands = np.array(bands, dtype=float).reshape((ispin, nkpts, nband))

    if os.path.isfile('KPOINTS'):
        kp = open('KPOINTS').readlines()

    if os.path.isfile('KPOINTS') and kp[2][0].upper() == 'L':
        Nk_in_seg = int(kp[1].split()[0])
        Nseg = nkpts // Nk_in_seg
        vkpt_diff = np.zeros_like(vkpts, dtype=float)
        
        for ii in range(Nseg):
            start = ii * Nk_in_seg
            end = (ii + 1) * Nk_in_seg
            vkpt_diff[start:end, :] = vkpts[start:end,:] - vkpts[start,:]

        kpt_path = np.linalg.norm(np.dot(vkpt_diff, B), axis=1)
        # kpt_path = np.sqrt(np.sum(np.dot(vkpt_diff, B)**2, axis=1))
        for ii in range(1, Nseg):
            start = ii * Nk_in_seg
            end = (ii + 1) * Nk_in_seg
            kpt_path[start:end] += kpt_path[start-1]

        # kpt_path /= kpt_path[-1]
        kpt_bounds =  np.concatenate((kpt_path[0::Nk_in_seg], [kpt_path[-1],]))
    else:
        # get band path
        vkpt_diff = np.diff(vkpts, axis=0)
        kpt_path = np.zeros(nkpts, dtype=float)
        kpt_path[1:] = np.cumsum(np.linalg.norm(np.dot(vkpt_diff, B), axis=1))
        # kpt_path /= kpt_path[-1]

        # get boundaries of band path
        xx = np.diff(kpt_path)
        kpt_bounds = np.concatenate(([0.0,], kpt_path[1:][np.isclose(xx, 0.0)], [kpt_path[-1],]))

    return kpt_path, bands, Efermi, kpt_bounds

##########################################################################################################
def read_procar_band(fname='PROCAR'):
    fc=open(fname,'r'); line=fc.readline(); line=fc.readline(); 
    list_k=line.split(' '); a=[]
    for word in list_k:
        if len(word)>0 and word != ' ':
            if word[-1]== '\n':
                a.append(word[:-1])
            else:
                a.append(word)
    nkpoints=int(a[3]);nbands=int(a[7]);natoms=int(a[-1])
    k_point=np.zeros((nkpoints , 3))
    band_energy=np.zeros((nkpoints,nbands))
    S=np.zeros((4,nkpoints,nbands))
    fc.close()
    #print (nkpoints, nbands, natoms)
    #print (S[0])
    
    fp=open(fname,'r'); line=fp.readline(); i_kpoint=0; k=1 #klines=[];
    while line:
        if " k-point " in line:
            print('reading kpoints: ', i_kpoint+1)
            kp=[]
            pattern=re.compile(r'(.{8})(.{6}) :(.{14})(.{11})(.{11})     weight = (\d+\.\d+)')
            matchs=pattern.match(line)
            kline=matchs.groups()
            kp.append(float(kline[2]));kp.append(float(kline[3]));kp.append(float(kline[4]))
            #print (kp)
            #klines.append(kline)
            k_point[i_kpoint]=kp
            for i_band in range(nbands):
            #for i_band in range(3):
                fp.readline();k=k+1  #tmp line
                band=fp.readline().rstrip(); k=k+1
                #print (k, band)
                band=band.split(' ')
                j=0
                for word in band:
                    if len(word)>0 and word != ' ':
                        j=j+1
                        if j==5:
                            band_energy[i_kpoint,i_band]=word;
                fp.readline();k=k+1
                fp.readline();k=k+1
                for i in range(4):
                    for i_ion in range(natoms):
                        fp.readline() ; k=k+1#tmp line
                        i_ion=i_ion+1
                    line1=fp.readline().rstrip() ; k=k+1
                    list_s=line1.split(' '); a=[]
                    for word in list_s:
                        if len(word)>0 and word != ' ':
                            a.append(word)
                    S[i,i_kpoint,i_band]=a[-1]

                    
            i_kpoint=i_kpoint+1
        line=fp.readline()
        k=k+1
    fp.close()
    
    fa=open('./k_points','w')
    fa=open('./k_points','a+')
    for i_kpoint in range(nkpoints):
        line='%15.9f%15.9f%15.9f\n' % (k_point[i_kpoint][0], k_point[i_kpoint][1],k_point[i_kpoint][2])
        fa.writelines(line)
    fa.close()
    
    return nkpoints, nbands, band_energy, S


def savedata(kpath, nkpoints, nbands, band_energy, spin_vec):
    #write k_points
    #write band_energy
    fb=open('./band.dat','w')
    fb=open('./band.dat','a+')
    fb.write("#BNAD  1\n")
    for i_band in range(nbands):
        for i_kpoint in range(nkpoints):
            line='%15.9f%15.9f\n' % (kpath[i_kpoint], band_energy[i_kpoint,i_band])
            fb.writelines(line)
        lie='#BNAD  '+str(i_band+2)+'\n'
        fb.write(lie)
    fb.close()

    #write band_energy with spin_vector
    fb=open('./band_spin_vector.dat','w')
    fb=open('./band_spin_vector.dat','a+')
    fb.write("#BNAD  1\n")
    for i_band in range(nbands):
        for i_kpoint in range(nkpoints):
            line='%15.9f%15.9f%15.9f%15.9f%15.9f%15.9f\n' % (kpath[i_kpoint],band_energy[i_kpoint,i_band], spin_vec[1,i_kpoint,i_band], spin_vec[2,i_kpoint,i_band],\
            spin_vec[3,i_kpoint,i_band], spin_vec[0,i_kpoint,i_band])
            fb.writelines(line)
        lie='#BNAD  '+str(i_band+2)+'\n'
        fb.write(lie)
    fb.close()

def bandplot(kpath, bands,  kpt_bounds, vbm_index, metal=False, ylim=None, kpts=None, EnergyWeight=None):
    from matplotlib.collections import LineCollection
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib.ticker import AutoMinorLocator
    '''
    'This part is copied and modified from pyband: https://github.com/QijingZheng/pyband'
    '''
    fig = plt.figure()
    fig.set_size_inches(8, 10)
    ax = plt.subplot(111)
    
    nspin, nkpts, nbands = bands.shape
    if EnergyWeight is not None:
        #norm = mpl.colors.Normalize(vmin=EnergyWeight.min(), vmax=EnergyWeight.max())
        norm = mpl.colors.Normalize(vmin=-0.5, vmax=0.5)
        
        s_m = mpl.cm.ScalarMappable(cmap='bwr', norm=norm)
        s_m.set_array([EnergyWeight])
        
        if nspin==1:
            if metal==True:
                shift=0.0
            else:
                shift=np.amax(bands[0, :, vbm_index-1])
                
            for Ik in range(nkpts):
                for Iband in range(nbands):
                    bands[0, Ik, Iband]-=shift

            for jj in range(nbands):
                x = kpath
                y = bands[0,:,jj]
                z = EnergyWeight[:,jj]

                ax.plot(x, y, lw=2, color='black', zorder=1)

                points = np.array([x, y]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                lc = LineCollection(segments,
                                    # cmap=opts.occLC_cmap, # alpha=0.7,
                                    colors=[s_m.to_rgba(ww) for ww in (z[1:] + z[:-1])/2.]
                                    # norm=plt.Normalize(0, 1)
                                    )
                # lc.set_array((z[1:] + z[:-1]) / 2)
                lc.set_linewidth(2)
                ax.add_collection(lc)

        divider = make_axes_locatable(ax)
        ax_cbar = divider.append_axes('right', size=0.1, pad=0.1)  
        cbar = plt.colorbar(s_m, cax=ax_cbar,orientation='vertical') #, ticks=[0.0, 1.0])
        
    else:
        if nspin==1:
            if metal==True:
                shift=0.0
            else:
                shift=np.amax(bands[0, :, vbm_index-1])
            
            for Ik in range(nkpts):
                for Iband in range(nbands):
                    bands[0, Ik, Iband]-=shift            
        for Ispin in range(nspin):
            for Iband in range(nbands):
                lc = 'black' 
                #lc = None if Iband == 0 else line.get_color()
                line, = ax.plot(kpath, bands[Ispin, :, Iband], lw=2, zorder=0, alpha=0.8, color=lc)
                '''
                if whts is not None:
                    for ii in range(len(opts.occ)):
                        ax.scatter(kpath, bands[Ispin, :, Iband],
                                color='red',
                                s=whts[ii][Ispin,:,Iband] * opts.occMarkerSize[ii],
                                marker=opts.occMarker[ii], zorder=1, lw=0.0,
                                alpha=0.5)
                '''
    for bd in kpt_bounds:
        ax.axvline(x=bd, ls='-', color='k', lw=0.5, alpha=0.5)
        
    ax.set_ylabel('Energy (eV)',  fontsize= 16, labelpad=5)
    
    if ylim is not None:
        ax.set_ylim(ylim[0],  ylim[1])
    ax.set_xlim(kpath.min(), kpath.max())

    ax.set_xticks(kpt_bounds)
    if kpts is not None:
        kname = [x.upper() for x in kpts]
        for ii in range(len(kname)):
            if kname[ii] == 'G':
                kname[ii] = r'$\mathrm{\mathsf{\Gamma}}$'
            else:
                kname[ii] = r'$\mathrm{\mathsf{%s}}$' % kname[ii]
        ax.set_xticklabels(kname, fontsize= 16)
    else:
        ax.set_xticklabels([])

    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    plt.tick_params('y',which='both',direction='in')
    plt.tick_params(bottom=False)
    plt.tight_layout(pad=0.20)
    plt.savefig('band_soc_vec.png', dpi=600)
    

if __name__ == "__main__":
    kpath, bands, efermi, kpt_bounds = get_bandInfo('./OUTCAR')
    nkpoints, nbands, band_energy, spin_vec=read_procar_band(fname='PROCAR')
    savedata(kpath, nkpoints, nbands, band_energy, spin_vec)
    
    #abs(total spin): spin_vec[0]
    #spin x :spin_vec[1]
    #spin y :spin_vec[2]
    #spin z :spin_vec[3]    
    bandplot( kpath, bands, kpt_bounds,vbm_index=40, metal=False, ylim=[-4,4], kpts='GMKG', EnergyWeight=spin_vec[3])













