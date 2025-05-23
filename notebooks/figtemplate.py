import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl

mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['font.family'] = 'DejaVu Sans'


def cm_to_inch(cm: float) -> float:
    return cm/2.54


##################################################
# Figure Size 
# 8.6 for single column 
# 17.6 for double column 

fig_width = 8.6
fig_height = 5

##################################################
# general Figure Setup
axlinewidth = 0.8
axcolor     ='grey'
axalpha     = 0.7
axtick_major_width  = 0.8
axtick_major_length = 2.4
axtick_minor_width  = 0.4
axtick_minor_length = 1.6

tick_pad        = 2
tick_labelsize  = 5
label_fontsize  = 6
legend_fontsize = 6

panel_label_fontsize = 8
label_fontweight= 'bold'
panel_label_fontweight= 'bold'

##################################################
# Plot Specs
marker_size     = 10
marker_linewidth = 0.7
plot_linewidth  = 1.2


###########################################################################################
### FIGURE SETUP
###########################################################################################

fig = plt.figure(figsize=(cm_to_inch(fig_width), cm_to_inch(fig_height)), dpi=300,facecolor='w',edgecolor='k') 
axes=[]
axes.append(plt.subplot2grid(shape=(1, 2), loc=(0, 0), colspan=1,rowspan=1))
plt.minorticks_on() 
axes.append(plt.subplot2grid(shape=(1, 2), loc=(0, 1), colspan=1,rowspan=1))
plt.minorticks_on() 


xlabel_pos = [0.5,-0.08]
ylabel_pos = [-0.13,0.5]


###########################################################################################
# Panel (a)

ax = axes[0]
# have a look at https://www.color-hex.com for nice color selection
ax.plot(np.linspace(0.2,0.8,10),np.random.normal(size=10),lw=plot_linewidth,zorder=1,color='#FF7901',label='random')


ax.set_xlabel('Translational Positioning',fontsize=label_fontsize,fontweight=label_fontweight)
ax.xaxis.set_label_coords(*xlabel_pos)

ax.set_ylabel(r'Free Energy ($\mathbf{k_{\mathbf{\mathrm{B}}} T}$)',fontsize=label_fontsize,fontweight=label_fontweight)
ax.yaxis.set_label_coords(*ylabel_pos)


###########################################################################################
# Panel (b)
ax = axes[1]

ax.set_xlabel('Translational Positioning',fontsize=label_fontsize,fontweight=label_fontweight)
ax.xaxis.set_label_coords(*xlabel_pos)

ax.set_ylabel(r'Template',fontsize=label_fontsize,fontweight=label_fontweight)
ax.yaxis.set_label_coords(*ylabel_pos)


axes[0].legend(fontsize=legend_fontsize,borderpad=0.2,framealpha=0.8,fancybox=True,handlelength=0.8,handletextpad=0.5,loc='lower left', bbox_to_anchor=(-0.03,0.99),ncol=5,columnspacing=0.8)


######################################################################################################################################################
######################################################################################################################################################
# Panel Labels

ax = axes[0]
ax.text(-0.16, 1.06, '(a)',
    fontsize=panel_label_fontsize,
    horizontalalignment='center',
    verticalalignment='center',
    transform = ax.transAxes,
    fontweight=panel_label_fontweight)

ax = axes[1]
ax.text(-0.16, 1.06, '(b)',
    fontsize=panel_label_fontsize,
    horizontalalignment='center',
    verticalalignment='center',
    transform = ax.transAxes,
    fontweight=panel_label_fontweight)

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################

##############################################
# Axes configs

for ax in axes:
    
    ###############################
    # set major and minor ticks
    # ax.tick_params(axis="both",which='major',direction="in",width=axtick_major_width,length=axtick_major_length,labelsize=tick_labelsize,pad=tick_pad,)
    # ax.tick_params(axis='both',which='minor',direction="in",width=axtick_minor_width,length=axtick_minor_length)
    ax.tick_params(axis="both",which='major',direction="in",width=axtick_major_width,length=axtick_major_length,labelsize=tick_labelsize,pad=tick_pad,color='#cccccc')
    ax.tick_params(axis='both',which='minor',direction="in",width=axtick_minor_width,length=axtick_minor_length,color='#cccccc')

    ###############################
    ax.xaxis.set_ticks_position('both')
    # set ticks right and top
    ax.yaxis.set_ticks_position('both')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(axlinewidth)
        ax.spines[axis].set_color(axcolor)
        ax.spines[axis].set_alpha(axalpha)

##############################################
# Setup subpanels

plt.subplots_adjust(
    left=0.10,
    right=0.97,
    bottom=0.12,
    top=0.90,
    wspace=0.25,
    hspace=0.25
    )

##################################################
# Save Figure
savefn = 'template'

fig.savefig(savefn+'.pdf',dpi=300,transparent=True)
fig.savefig(savefn+'.svg',dpi=300,transparent=True)
fig.savefig(savefn+'.png',dpi=300,transparent=False)
