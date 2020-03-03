from modules import *



#"""

hx = h5py.File('selected_infaller_xs.hdf5', 'w')
hy = h5py.File('selected_infaller_ys.hdf5', 'w')
hz = h5py.File('selected_infaller_zs.hdf5', 'w')
#hvx = h5py.File('selected_infaller_vx.hdf5', 'w')
#hvy = h5py.File('selected_infaller_vy.hdf5', 'w')
#hvz = h5py.File('selected_infaller_vz.hdf5', 'w')
hr = h5py.File('selected_infaller_rs.hdf5', 'w')

crange = np.array(np.loadtxt('/home/ppxrh2/Documents/test_pollux/TheThreeHundred/playground/rhaggar/G3X_data/G3X_300_selected_sample_257.txt'), dtype='int')
all_relative_clus_data = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
        + 'Infalling_Groups/MergerTreeAHF_Infalling_Re-written/'
        + 'all_members_halo_data/rel_to_cluster/')
counter = 0

for c_no in crange:
    print(c_no)

    data_xs = h5py.File(all_relative_clus_data + 'xs/CLUSTER_%04d_xs_reltoCLUS.hdf5' % c_no, 'r')
    data_ys = h5py.File(all_relative_clus_data + 'ys/CLUSTER_%04d_ys_reltoCLUS.hdf5' % c_no, 'r')
    data_zs = h5py.File(all_relative_clus_data + 'zs/CLUSTER_%04d_zs_reltoCLUS.hdf5' % c_no, 'r')
    data_vx = h5py.File(all_relative_clus_data + 'vx/CLUSTER_%04d_vx_reltoCLUS.hdf5' % c_no, 'r')
    data_vy = h5py.File(all_relative_clus_data + 'vy/CLUSTER_%04d_vy_reltoCLUS.hdf5' % c_no, 'r')
    data_vz = h5py.File(all_relative_clus_data + 'vz/CLUSTER_%04d_vz_reltoCLUS.hdf5' % c_no, 'r')
    data_ms = h5py.File(all_relative_clus_data + 'ms/CLUSTER_%04d_ms_reltoCLUS.hdf5' % c_no, 'r')
    data_mstar = h5py.File(all_relative_clus_data + 'mstar/CLUSTER_%04d_mstar_reltoCLUS.hdf5' % c_no, 'r')

    keys = np.array(list(data_xs.keys()))
    keys_l = []
    lowsizelim = 1 #non-inclusive
    highsizelim = 10 #non-inclusive
    for i in range(len(keys)):
        if len(list(data_xs[keys[i]].keys())) > lowsizelim:
            keys_l += [keys[i]]
    keys_l = np.array(keys_l)

    for i in range(len(keys_l)):
        keys_l_i = keys_l[i]

        data_xs_n = data_xs[keys_l_i]
        data_ys_n = data_ys[keys_l_i]
        data_zs_n = data_zs[keys_l_i]
        data_vx_n = data_vx[keys_l_i]
        data_vy_n = data_vy[keys_l_i]
        data_vz_n = data_vz[keys_l_i]
        data_ms_n = data_ms[keys_l_i]
        data_mstar_n = data_mstar[keys_l_i]

        keys = np.array(list(data_xs_n.keys()))
        X = np.array(data_xs_n[keys[0]])[0]
        Y = np.array(data_ys_n[keys[0]])[0]
        Z = np.array(data_zs_n[keys[0]])[0]
        R = (X**2. + Y**2. + Z**2.)**0.5
        theta = np.arctan2(Y, X)
        phi = np.arccos(Z/R)

        #############testing to see if enough satellites satisfy mass criteria
        M, MS = [], []
        for k in range(len(keys)):
            M += [int(np.array(data_ms_n[keys[k]])[0])]
            MS += [int(np.array(data_mstar_n[keys[k]])[0])]
        M, MS = np.array(M), np.array(MS)
        RAT = MS/M
        mass_bool = (M > 10.**10.5) * (MS > 10.**9.5) * (RAT < 0.3)
        if np.sum(mass_bool) <= lowsizelim or np.sum(mass_bool) >= highsizelim:
            keys = np.zeros(0)
        #############
        
        if len(keys)>0:
            M = np.array(data_ms_n[keys[0]])[0]
            MS = np.array(data_mstar_n[keys[0]])[0]
            RAT = MS/M
            mass_bool = (M > 10.**10.5) * (MS > 10.**9.5) * (RAT < 0.3)
            if mass_bool==False:
                keys = np.zeros(0)
        

        for j in range(len(keys)):
            xs = np.array(data_xs_n[keys[j]])
            ys = np.array(data_ys_n[keys[j]])
            zs = np.array(data_zs_n[keys[j]])
            vx = np.array(data_vx_n[keys[j]])
            vy = np.array(data_vy_n[keys[j]])
            vz = np.array(data_vz_n[keys[j]])
            ms = np.array(data_ms_n[keys[j]])
            mstar = np.array(data_mstar_n[keys[j]])
            rat = mstar / ms
            rs = (xs**2. + ys**2. + zs**2.)**0.5
            
            #ts = np.arctan2(ys, xs) - theta
            #ps = np.arccos(zs/rs) + (np.pi/2) - phi
            #xs_n = (1./R) * (xs*X + ys*Y + zs*Z*np.cos(ts-theta))
            #ys_n = (1./R) * (ys*X - xs*Y + zs*Z*np.sin(ts-theta))

            #xs_n = rs * np.cos(ts)*np.sin(ps)
            #ys_n = rs * np.sin(ts)*np.sin(ps)
            #zs_n = rs * np.cos(ps)

            s_theta, c_theta = np.sin(theta), np.cos(theta)
            s_phi, c_phi = np.sin(phi), np.cos(phi)

            xs_n = (zs*c_phi) + s_phi*((xs*c_theta)+(ys*s_theta))
            ys_n = (ys*c_theta) - (xs*s_theta)
            zs_n = (zs*s_phi) - c_phi*((xs*c_theta)+(ys*s_theta))

            vx_n = (vz*c_phi) + s_phi*((vx*c_theta)+(vy*s_theta))
            vy_n = (vy*c_theta) - (vx*s_theta)
            vz_n = (vz*s_phi) - c_phi*((vx*c_theta)+(vy*s_theta))

            ms_n = copy.deepcopy(ms)
            mstar_n = copy.deepcopy(mstar)
            rat_n = copy.deepcopy(rat)


            rs_cut_lim = 0
            rs_cut = np.where(rs[1:] > rs[:-1])[0]
            if len(rs_cut)>0:
                rs_cut = rs_cut[0]
                rs_n = rs[rs_cut:]
                rs_cut_n = np.where(rs_n[1:] < rs_n[:-1])[0]
                if len(rs_cut_n)>0:
                    rs_cut_lim = rs_cut + rs_cut_n[0] + 1
                    xs_n = xs_n[:rs_cut_lim]
                    ys_n = ys_n[:rs_cut_lim]
                    zs_n = zs_n[:rs_cut_lim]
                    vx_n = vx[:rs_cut_lim]
                    vy_n = vy[:rs_cut_lim]
                    vz_n = vz[:rs_cut_lim]
                    ms_n = ms[:rs_cut_lim]
                    mstar_n = mstar[:rs_cut_lim]
                    rat_n = rat[:rs_cut_lim]
            mass_bool = (ms_n > 10.**10.5) * (mstar_n > 10.**9.5) * (rat_n < 0.3)
            xs_n, ys_n, zs_n = xs_n[mass_bool], ys_n[mass_bool], zs_n[mass_bool]
            rs_n = (ys_n**2. + zs_n**2.)**0.5
            rs_n[ys_n<0.] = -1.*rs_n[ys_n<0.]
            if j==0:
                flip_dir = 1. #1 for keep, -1 for reverse direction
                if len(ys_n)>1:
                    if ys_n[0] > ys_n[1]:
                    #if vy_n[0] < 0.:
                        flip_dir = -1.
            rs_n *= flip_dir


            if len(xs_n) > 1:
                hx.create_dataset('%07d' % counter, data=xs_n)
                hy.create_dataset('%07d' % counter, data=ys_n)
                hz.create_dataset('%07d' % counter, data=zs_n)
                hr.create_dataset('%07d' % counter, data=rs_n)
                counter += 1
print(counter)
hx.close()
hy.close()
hz.close()
hr.close()
#"""




plt.figure(figsize=(7, 6))
hx = h5py.File('selected_infaller_xs.hdf5', 'r')
hy = h5py.File('selected_infaller_rs.hdf5', 'r')
#hy = h5py.File('selected_infaller_zs.hdf5', 'r')
keys = np.array(list(hx.keys()))
xs_all, ys_all = np.zeros(0), np.zeros(0)
main=0
for i in range(len(keys)):
    xs_all = np.append(xs_all, np.array(hy[keys[i]]))
    ys_all = np.append(ys_all, -1.*np.array(hx[keys[i]]))
    #plt.plot(-1.*np.array(hx[keys[i]]), hy[keys[i]], alpha=0.5, c='b')
    #if main==0:
    #    plt.plot(-1.*np.array(hx[keys[i]]), hy[keys[i]], alpha=1, linewidth=4., c='b')
    #main=1
plt.xlim(-1.5, 2.)
plt.ylim(-1.5, 1.5)
plt.xlabel(r'$x/R_{200}$')
plt.ylabel(r'$y/R_{200}$')
#plt.tight_layout()
#plt.savefig('example_002.png', dpi=400)
#plt.show()

#xbins = np.linspace(-1.5, 1.5, 30)
#ybins = (np.arange(36)-15)/10.
res = 10
xbins = np.linspace(-1.5, 1.5, 6*res)
ybins = (np.arange((7*res)+1)-(3*res))/(2*res)
xvals = np.zeros((len(ybins)-1, 3))
xvals2 = np.zeros((len(ybins)-1, 3))
xvals3 = np.zeros((len(ybins)-1, 3))
xvals4 = np.zeros((len(ybins)-1, 3))
#xhist = np.zeros((len(xbins)-1, len(ybins)-1))
for i in range(len(ybins)-1):
    where = np.where((ys_all > ybins[i]) * (ys_all < ybins[i+1]))[0]
    xvals[i] = find_stdev(xs_all[where], 0.5)
    xvals2[i] = find_stdev(xs_all[where], 0.2)
    xvals3[i] = find_stdev(xs_all[where], 0.05)
    xvals4[i] = find_stdev(xs_all[where], 0.8)
    #xhist[:, i] = np.histogram(xs_all[where], xbins)[0]
    #xhist[:, i] = np.log10(xhist[:, i] + 1.*(xhist[:, i]<1))

xhist = plt.hist2d(ys_all, xs_all, bins=(len(ybins)-1, len(xbins)-1), range=[[-1.5, 2.], [-1.5, 1.5]])[0]

xhist[:, (len(xbins)//2)-1] *= 2.
##########NNOTE: NEEDS NORMALISING TO VOLUME OF ANNULI: BIAS AGAINST R=0 CURRENTLY#############
plt.clf()
xhist = np.transpose(np.log10(xhist + 1.*(xhist<1)))

xbins += 0.5*(xbins[1]-xbins[0])
ybins += 0.5*(ybins[1]-ybins[0])
#plt.plot(ybins[5:-1], xvals[5:,0], c='k', linewidth=3.)
plt.plot(ybins[res:-1], xvals[res:,0]+xvals[res:,1], c='k', linewidth=2.)
plt.plot(ybins[res:-1], xvals[res:,0]-xvals[res:,2], c='k', linewidth=2.)
plt.plot(ybins[res:-1], xvals2[res:,0]+xvals2[res:,1], c='k', linewidth=1.)
plt.plot(ybins[res:-1], xvals2[res:,0]-xvals2[res:,2], c='k', linewidth=1.)
plt.plot(ybins[res:-1], xvals3[res:,0]+xvals3[res:,1], c='k', linewidth=1., linestyle=':')
plt.plot(ybins[res:-1], xvals3[res:,0]-xvals3[res:,2], c='k', linewidth=1., linestyle=':')
#plt.plot(ybins[5:-1], xvals4[5:,0]+xvals4[5:,1], c='k', linewidth=1.)
#plt.plot(ybins[5:-1], xvals4[5:,0]-xvals4[5:,2], c='k', linewidth=1.)

def find_angles(ys, xs):
    ys[:,1] = ys[:,0]+ys[:,1]
    ys[:,2] = ys[:,0]-ys[:,2]
    rs1 = (xs**2. + ys[:,1]**2.)**0.5
    rs2 = (xs**2. + ys[:,2]**2.)**0.5
    rs1, rs2 = rs1[1:], rs2[1:]
    int1 = np.where(rs1<1.)[0][-1]+1
    int2 = np.where(rs2<1.)[0][-1]+1
    xs1 = xs[int1:int1+2]
    xs2 = xs[int2:int2+2]
    ys1 = ys[int1:int1+2, 1]
    ys2 = ys[int2:int2+2, 2]
    ang1 = find_where_r1(xs1, ys1)
    ang2 = find_where_r1(xs2, ys2)
    return ang1+ang2

def find_where_r1(xs, ys):
    """ Find point on line [x0, x1], [y0, y1] where x^2+y^2=1 """
    m = (ys[1]-ys[0]) / (xs[1]-xs[0])
    c = ys[0] - m*xs[0]
    A = m**2. + 1.
    B = 2.*m*c
    C = c**2. - 1.
    x = ((-1.*B) + (B**2. - 4.*A*C)**0.5) / (2.*A)
    y = (1. - x**2.)**0.5
    angle = np.arctan(y/x)*180./np.pi
    return angle

#print(find_angles(xvals4[5:], ybins[5:-1]))
#print(find_angles(xvals[5:], ybins[5:-1]))
#print(find_angles(xvals2[5:], ybins[5:-1]))
#print(find_angles(xvals3[5:], ybins[5:-1]))




rng = np.arange(101)*2.*np.pi/100.
plt.plot(np.sin(rng), np.cos(rng), linewidth=1., c='k', linestyle='--', zorder=999999)
plt.scatter(0., 0., s=120., marker='x', c='r', zorder=1000000)
plt.scatter(-1.15, 0., s=250., marker='>', c='r', zorder=1000000)
plt.plot([-1.4, -1.2], [0., 0.], linewidth=5., c='r', zorder=1000000)

plt.imshow(xhist, extent=(-1.5, 2., -1.5, 1.5), origin='lower', cmap='viridis')
Y, X = np.meshgrid(xbins[:-1], ybins[:-1])
plt.xlabel(r'$x/R_{200}$')
plt.xlabel('hiiiiiiiiiiii roan')
plt.ylabel(r'$y/R_{200}$')
plt.xticks((np.arange(8)-3)*0.5)
plt.yticks((np.arange(7)-3)*0.5)
plt.tight_layout()
#plt.savefig('hist2d_rs_poster.png', dpi=600, facecolor='#333333')
plt.savefig('hist2d_2-9.png', dpi=400)

#plt.show()


hx.close()
hy.close()






"""
#3d plotting stuff
#plt.figure(figsize=(7, 6))
fig = plt.figure(figsize=(10,10))
ax = Axes3D(fig)
ax.set_aspect('equal')
ax.xaxis.pane.fill =False
ax.yaxis.pane.fill =False
ax.zaxis.pane.fill =False
ax.grid(False)
ax.set_axis_off()
ax.view_init(50,90)

u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = np.cos(u)*np.sin(v)
y = np.sin(u)*np.sin(v)
z = np.cos(v)
ax.plot_wireframe(x, y, z, color='k', alpha=0.3, zorder=1000) 

hx = h5py.File('selected_infaller_xs.hdf5', 'r')
hy = h5py.File('selected_infaller_ys.hdf5', 'r')
hz = h5py.File('selected_infaller_zs.hdf5', 'r')
keys = np.array(list(hx.keys()))
xs_all, ys_all = np.zeros(0), np.zeros(0)
main=0


for i in range(len(keys))[-48:-26]:
    xs_all = np.append(xs_all, np.array(hy[keys[i]]))
    ys_all = np.append(ys_all, -1.*np.array(hx[keys[i]]))
    plt.plot(np.array(hx[keys[i]]), hy[keys[i]], hz[keys[i]], alpha=0.8)#, c='b')
    ax.scatter(np.array(hx[keys[i]])[-1], hy[keys[i]][-1], hz[keys[i]][-1], alpha=0.8, s=30., linewidth=0.)#, c='b')
    if main==0:
        plt.plot(np.array(hx[keys[i]]), hy[keys[i]], hz[keys[i]], alpha=1, linewidth=3., c='b', zorder=999)
        ax.scatter(np.array(hx[keys[i]])[-1], hy[keys[i]][-1], hz[keys[i]][-1], alpha=1, c='b', s=50., zorder=999)
    main=1
#plt.xlim(-1.5, 2.)
#plt.ylim(-1.5, 1.5)
plt.scatter(0., 0., s=100., marker='x', c='r', zorder=999)
plt.scatter(1.1, 0., s=100., marker='>', c='r', zorder=999)
ax.set_xlim(-2.8, 0.7)
ax.set_ylim(-1.5, 2.)
ax.set_zlim(-1.5, 2.)
plt.xlabel(r'$x/R_{200}$')
plt.ylabel(r'$y/R_{200}$')
#plt.tight_layout()
plt.savefig('example_3d.png', dpi=400)
plt.show()
print(alphabet)
"""
