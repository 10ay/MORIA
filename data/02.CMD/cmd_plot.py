import numpy as np
import matplotlib.pyplot as plt

#Load the MATCHUP files


#xv, yv, mv are the x, y and magnitudes for the V band. Same logic for the I band

xv, yv, mv = np.loadtxt("MATCHUP.F606W.XYM", unpack=True, usecols=(0,1,2))
xi, yi, mi = np.loadtxt("MATCHUP.F814W.XYM.02", unpack=True, usecols=(0,1,2))

#Establish target parameters

xtarg, ytarg = 535.2320, 623.8950
Vtarg, Itarg = -10.2651, -10.3052
VmItarg = Vtarg - Itarg

#Function to find the CMD of the target and get our list of Sim+Ref stars.

def show_cmd_targ():

    #Plotting parameters 

    rad_max, box_max = 300, 300
    mag_range, col_range = 0.50, 0.30
    ref_st_Imx, ref_st_Imn = -8.8, -12.75
    ref_st_Vmx, ref_st_Vmn = -8.7, -12.75

    d = np.sqrt((xi-xtarg)**2 + (yi-ytarg)**2)
    n = np.arange(1, len(mi)+1)

    #Plotting parameter to decide how many stars around the target should be selected 

    vprox = np.abs(mi - Itarg) < mag_range
    cprox = np.abs(mv - mi - VmItarg) < col_range

    #vprox and cprox are boolean arrays. Stars close to our target within a given window 

    u = vprox & cprox #These are the stars close to the target in both magnitude and color space.


    uref = (d < rad_max) & cprox & (mv < ref_st_Vmx) & (mv > ref_st_Vmn) & (mi < ref_st_Imx) & (mi > ref_st_Imn) & (n > 1) #More selective than u

    fig, ax = plt.subplots(2, 2, figsize=(10, 10))
    ax_cmd = ax[0,0]
    ax_xy_I = ax[1,1]
    ax_xy_zoom = ax[1,0]
    ax_xy_V = ax[0, 1]

    # CMD 
    ax_cmd.scatter(mv-mi, mi, s=5, c='k', label='All Stars')
    ax_cmd.scatter((mv-mi)[u], mi[u], s=15, c='purple', label='Selected Stars')
    ax_cmd.scatter([VmItarg], [Itarg], marker='x', lw = 5, s=100, c='grey', label='Target')
    ax_cmd.set_xlim(-0.75, 1.75)
    ax_cmd.set_ylim(-15, -7)   
    ax_cmd.set_xlabel('F606W - F814W')
    ax_cmd.set_ylabel('F814W')
    ax_cmd.legend(loc = 'lower right')
    ax_cmd.set_title('CMD')

    # XY for I filter
    ax_xy_I.scatter(xi, yi, s=5, c='k', label = 'All Stars') # All stars
    #ax_xy_I.scatter(xi[u], yi[u], s=15, c='r') # I avoid this from SM and prefer uref instead
    ax_xy_I.scatter(xi[uref], yi[uref], s=30, c='purple', label = 'Actual selected stars')
    ax_xy_I.scatter([xtarg], [ytarg], marker='x', lw = 5, s=100, c='grey', label = 'Target') #Target
    circle = plt.Circle((xtarg, ytarg), rad_max, color='hotpink', fill=False, lw=2) 
    ax_xy_I.add_patch(circle)
    ax_xy_I.set_xlim(xtarg-box_max, xtarg+box_max)
    ax_xy_I.set_ylim(ytarg-box_max, ytarg+box_max)
    ax_xy_I.set_xlabel('x coord')
    ax_xy_I.set_ylabel('y coord')
    ax_xy_I.legend(loc = 'lower right')
    ax_xy_I.set_title('XY Coord I-Filter')

    # XY for V filter
    ax_xy_V.scatter(xv, yv, s=5, c='k', label = 'All Stars') # All stars
    ax_xy_V.scatter(xv[uref], yv[uref], s=30, c='purple', label = 'Actual selected stars')
    ax_xy_V.scatter([xtarg], [ytarg], marker='x',lw = 5, s=100, c='grey', label = 'Target') #Target
    circle = plt.Circle((xtarg, ytarg), rad_max, color='hotpink', fill=False, lw=2) 
    ax_xy_V.add_patch(circle)
    ax_xy_V.set_xlim(xtarg-box_max, xtarg+box_max)
    ax_xy_V.set_ylim(ytarg-box_max, ytarg+box_max)
    ax_xy_V.set_xlabel('x coord')
    ax_xy_V.set_ylabel('y coord')
    ax_xy_V.legend(loc = 'lower right')
    ax_xy_V.set_title('XY Coord V-Filter')


    # XY zoom for I-filter
    mask_zoom = u & (d < rad_max)
    ax_xy_zoom.scatter(xi, yi, s=5, c='k', label = 'All Stars')
    #ax_xy_zoom.scatter(xi[mask_zoom], yi[mask_zoom], s=15, c='r')
    ax_xy_zoom.scatter(xi[uref], yi[uref], s=30, c='purple', label = 'Actual selected stars')
    ax_xy_zoom.scatter([xtarg], [ytarg], marker='x', lw = 5, s=100, c='grey', label = 'Target')
    circle2 = plt.Circle((xtarg, ytarg), rad_max, color='hotpink', fill=False, lw=2)
    ax_xy_zoom.add_patch(circle2)
    ax_xy_zoom.set_xlim(xtarg-box_max, xtarg+box_max)
    ax_xy_zoom.set_ylim(ytarg-box_max, ytarg+box_max)
    ax_xy_zoom.set_xlabel('x coord')
    ax_xy_zoom.set_ylabel('y coord')
    ax_xy_zoom.legend(loc='lower right')
    ax_xy_zoom.set_title('XY Coord I-Filter Zoom')

    plt.tight_layout()
    plt.savefig("show_cmd_targ.pdf")
    plt.close(fig)

    xu = xi[u & (d < rad_max)]
    yu =  yi[u & (d < rad_max)]
    miu = mi[u & (d < rad_max)]
    mvu = mv[u & (d < rad_max)]

    upsf = np.ones_like(xu, dtype=int)

    assert len(upsf) > 0
#    if len(upsf) > 0 could be another condition

    upsf[0] = 0  # first star (the target) can't be used for PSF 
 

    np.savetxt('NEARBY_SIM_STARS.XYIVB_targ', np.column_stack([xu, yu, miu, mvu, upsf]), fmt='%10.3f %10.3f %8.4f %8.4f %1d', header="xu         yu      miu      mvu   upsf \n")

    np.savetxt('NEARBY_REF_STARS.XYIVB_targ', np.column_stack([xi[uref], yi[uref], mi[uref], mv[uref]]), fmt='%10.3f %10.3f %8.4f %8.4f', header = "xuref      yuref   miuref      mvuref \n")

#Function to give calibration stars
def show_cmd_Cal():
    #Plotting parameters 

    rad_max, box_max = 450, 480
    Vcalc, Icalc = -13.0, -13.0
    mag_range, col_range = 1.9, 1.2
    ref_st_Imx, ref_st_Imn = Icalc + mag_range, Icalc - mag_range
    ref_st_Vmx, ref_st_Vmn = Vcalc + mag_range, Vcalc - mag_range

    d = np.sqrt((xi-xtarg)**2 + (yi-ytarg)**2)
    vprox = np.abs(mi - Icalc) < mag_range
    cprox = np.abs(mv - mi - VmItarg) < col_range

    u = vprox & cprox

    uref = (d < rad_max) & cprox & (mv < ref_st_Vmx) & (mv > ref_st_Vmn) & (mi < ref_st_Imx) & (mi > ref_st_Imn)

    fig, ax = plt.subplots(1, 2, figsize=(10, 5))

    # CMD
    ax_cmd = ax[0]
    ax_cmd.scatter(mv-mi, mi, s=5, c='k', label='All Stars')
#    ax_cmd.scatter((mv-mi)[u], mi[u], s=15, c='r', label='Selected Stars')
    ax_cmd.scatter((mv-mi)[u], mi[u], s=15, c='purple', label='Selected Stars')
    ax_cmd.scatter([VmItarg], [Itarg], marker='x', lw =5, s=100, c='grey', label='Target')
    ax_cmd.set_xlim(-0.75, 1.25)
    ax_cmd.set_ylim(-15, -7)
    ax_cmd.set_xlabel('F606W - F814W')
    ax_cmd.set_ylabel('F814W')
    ax_cmd.set_title('CMD')
    ax_cmd.legend(loc = 'lower right')

    # XY
    ax_xy = ax[1]
    ax_xy.scatter(xi, yi, s=5, c='k', label = 'All Stars')
    #ax_xy.scatter(xi[u], yi[u], s=15, c='purple')
    ax_xy.scatter(xi[uref], yi[uref], s=30, c='purple', label = 'Actual Selected Stars')
    ax_xy.scatter([xtarg], [ytarg], marker='x', lw = 5, s=100, c='grey', label = 'Target')
    circle = plt.Circle((xtarg, ytarg), rad_max, color='hotpink', fill=False, lw=2)
    ax_xy.add_patch(circle)
    ax_xy.set_xlim(xtarg-box_max, xtarg+box_max)
    ax_xy.set_ylim(ytarg-box_max, ytarg+box_max)
    ax_xy.set_xlabel('x coord')
    ax_xy.set_ylabel('y coord')
    ax_xy.set_title('Calibration Stars')
    ax_xy.legend(loc = 'lower right')
    plt.tight_layout()
    plt.savefig('show_cmd_Cal.pdf')
    plt.close(fig)

    xuref = xi[uref]
    yuref =  yi[uref]
    miuref = mi[uref]
    mvuref = mv[uref]

    iMuref = np.arange(1, len(xi) + 1)[uref]
    upsf = np.ones_like(xuref, dtype=int)

    np.savetxt('NOTFAR_CAL_STARS.XYIVB_targ',  np.column_stack([xuref, yuref, miuref, mvuref, upsf, iMuref]), fmt='%10.3f %10.3f %8.4f %8.4f %1d %6d', header=" #  xuref      yuref    miuref   mvuref  upsf iMuref  \n")


show_cmd_targ()
show_cmd_Cal()
