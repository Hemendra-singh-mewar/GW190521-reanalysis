"""Figures with shared noise conditioning and explicit statistic labels."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from .conditioning import gwpy_asd
COLOURS={'H1':'#bb3a4a','L1':'#167a91','V1':'#cc972d'}

def style():
    plt.rcParams.update({'font.size':10,'figure.dpi':120,'savefig.dpi':180,
        'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.18,
        'pdf.fonttype':42})

def save(fig,name,c):
    fig.tight_layout()
    for extension in ['png','pdf']:
        fig.savefig(c.figures_dir/f'{name}.{extension}',bbox_inches='tight')
    plt.close(fig)

def save_all(table,matches,psds,event_data,conditioned,curves,scan,c,make_q=True):
    style()
    fig,ax=plt.subplots(figsize=(7.5,4.2))
    for ifo,p in psds.items():
        f=np.asarray(p.sample_frequencies); good=(f>=10)&(f<=512)
        ax.loglog(f[good],np.sqrt(np.asarray(p)[good]),color=COLOURS[ifo],label=ifo,lw=.9)
    ax.axvspan(c.f_low,c.f_high,alpha=.06,color='tab:blue')
    ax.set(xlabel='Frequency [Hz]',ylabel=r'ASD [strain / $\sqrt{\mathrm{Hz}}$]')
    ax.legend(frameon=False);save(fig,'detector_psd',c)
    fig,axs=plt.subplots(3,1,figsize=(7.5,5.5),sharex=True)
    for ax,(ifo,d) in zip(axs,conditioned.items()):
        ax.plot(np.asarray(d.times)-c.gps,d.value,color=COLOURS[ifo],lw=.65)
        ax.set_ylabel(ifo);ax.set_xlim(-.3,.3)
    axs[-1].set_xlabel('Detector time relative to catalogue GPS [s]')
    fig.supylabel('Whitened, filtered strain');save(fig,'conditioned_strain',c)
    if make_q:
        fig,axs=plt.subplots(3,1,figsize=(7.5,6),sharex=True,sharey=True)
        for ax,ifo in zip(axs,c.ifos):
            q=event_data[ifo].q_transform(qrange=(4,64),frange=(20,256),
                outseg=(c.gps-.3,c.gps+.3),gps=c.gps,search=.6,
                whiten=gwpy_asd(psds[ifo]),fduration=4,highpass=20)
            im=ax.pcolormesh(np.asarray(q.xindex)-c.gps,np.asarray(q.yindex),q.value.T,
                            cmap='viridis',shading='auto',vmin=0,vmax=25,rasterized=True)
            ax.set_ylabel(ifo+'\nFrequency [Hz]')
            fig.colorbar(im,ax=ax,label='Normalised energy',pad=.015)
        axs[-1].set_xlabel('Detector time relative to catalogue GPS [s]');save(fig,'qtransform',c)
    fig,ax=plt.subplots(figsize=(8,4.7))
    cols=['Independent_peak_bound','Time_aligned','Fixed_response_coherent']
    labels=['Independent peak bound','Common time, free detector phases','Fixed response, common amplitude']
    x=np.arange(len(table));w=.24
    for i,(col,label) in enumerate(zip(cols,labels)):
        ax.bar(x+(i-1)*w,table[col],w,label=label)
    ax.set_xticks(x,table.index,rotation=10,ha='right');ax.set_ylabel('Network statistic')
    ax.legend(frameon=False,fontsize=8);save(fig,'snr_comparison',c)
    fig,ax=plt.subplots(figsize=(6.7,5.1))
    im=ax.imshow(matches,vmin=0,vmax=1,cmap='viridis');ax.grid(False)
    ax.set_xticks(range(len(matches)),matches.columns,rotation=20,ha='right')
    ax.set_yticks(range(len(matches)),matches.index)
    for i in range(len(matches)):
        for j in range(len(matches)):
            v=matches.iloc[i,j];ax.text(j,i,f'{v:.3f}',ha='center',va='center',color='black' if v>.6 else 'white')
    fig.colorbar(im,ax=ax,label='H1 PSD weighted match');save(fig,'match_matrix',c)
    fig,axs=plt.subplots(2,1,figsize=(7.5,5.6),sharex=True)
    for name,curve in curves.items():
        axs[0].plot(curve['offset_s'],curve['time_aligned'],label=name,lw=1.2)
        axs[1].plot(curve['offset_s'],curve['coherent'],label=name,lw=1.2)
    axs[0].set_ylabel('Time aligned statistic');axs[1].set_ylabel('Fixed response statistic')
    axs[1].set_xlabel('Geocentric reference time relative to catalogue GPS [s]')
    axs[0].legend(frameon=False,fontsize=8,ncol=2);save(fig,'network_timeseries',c)
    if scan is not None:
        matrix=scan.pivot(index='q',columns='frequency_hz',values='Time_aligned')
        fig,ax=plt.subplots(figsize=(7.2,4.4))
        im=ax.pcolormesh(matrix.columns,matrix.index,matrix,shading='nearest',cmap='viridis')
        idx=scan['Time_aligned'].idxmax();b=scan.loc[idx]
        ax.plot(b.frequency_hz,b.q,'w+',ms=12,mew=2)
        ax.set(xlabel='Sine Gaussian central frequency [Hz]',ylabel='Quality factor Q')
        fig.colorbar(im,ax=ax,label='Time aligned statistic');save(fig,'sine_gaussian_scan',c)
