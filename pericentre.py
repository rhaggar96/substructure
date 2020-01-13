from modules import * 

group_halo_survivors = True



out_dir = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/Infalling_Groups'
        '/MergerTreeAHF_Infalling_Groups/')
track_dir = ('/home/ppxrh2/Documents/test_pollux/TheThreeHundred/playground/'
        'rhaggar/substructure_data/tracking_data/')

crange = np.array(np.loadtxt('/home/ppxrh2/Documents/test_pollux/TheThree'
        'Hundred/playground/rhaggar/G3X_data/G3X_300_selected_sample_257'
        '.txt'), dtype='int')
#crange = np.array(np.arange(1, 325), dtype='int')

def bound_crit(x):
    """ y = v/vcrit, x = r/R200, find boundness condition """
    x[x>2.5]=2.5
    x[x<0.000001] = 0.000001
    return ((5./(3.*x)) - (2./3.))**0.5


total = ld_arr(track_dir+'forward_tracking.txt')
final_s_in = ld_arr(track_dir+'forward_tracking_finalsnap.txt')

clus_id = np.array(total[:, 4], dtype='int')
if len(crange)<324:
    for i in range(len(clus_id)):
        if not clus_id[i] in crange:
            clus_id[i] = 0
    clus_bl = clus_id>0
else:
    clus_bl = np.array(total[:, 4], dtype='int') > 0

if group_halo_survivors==True: #consider only z=0 groups
    grp_surv_bl = final_s_in[:, 1] > 128*mod
else:
    grp_surv_bl = final_s_in[:, 1] > 0

total = total[clus_bl*grp_surv_bl] 
final_s_in = final_s_in[clus_bl*grp_surv_bl]
clus_id = clus_id[clus_bl*grp_surv_bl]






rs_total, vs_total, id_total = total[:, 0], total[:, 1], total[:, 2]//mod

xmin, xmax, ymin, ymax, steps = 0., 2.6, 0., 2.5, 100


final_s = np.array(final_s_in[:, :2], dtype='int')//mod
final_s_bd = np.array(final_s_in[:, 2], dtype='int')

grp_total = np.char.mod('%15d', total[:, 3])
grp_peri_bool = np.zeros(len(grp_total))
for i in range(len(grp_total)):
    if clus_id[i] != clus_id[i-1]:
        ld_str = ('Group_member_tracking/CLUSTER_%04d/position_wrt_cluster/'
                'CLUSTER_%04d_member_') % (clus_id[i], clus_id[i])
        peri_data_xs_i = h5py.File(out_dir + ld_str + 'xs_cluster.hdf5')
        peri_data_ys_i = h5py.File(out_dir + ld_str + 'ys_cluster.hdf5')
        peri_data_zs_i = h5py.File(out_dir + ld_str + 'zs_cluster.hdf5')
        peri_data_r2_i = h5py.File(out_dir + ld_str + 'r200host_cluster.hdf5')
    peri_data_xs = np.array(peri_data_xs_i[grp_total[i]]['              0'])
    peri_data_ys = np.array(peri_data_ys_i[grp_total[i]]['              0'])
    peri_data_zs = np.array(peri_data_zs_i[grp_total[i]]['              0'])
    peri_data_r2 = np.array(peri_data_r2_i[grp_total[i]]['              0'])
    peri_rs = (peri_data_xs**2. + peri_data_ys**2. + peri_data_zs**2.
            )**0.5 / (0.000001+peri_data_r2)
    if peri_rs[-1] > np.min(peri_rs):
        grp_peri_bool[i] = 1.
grp_peri_bool = grp_peri_bool > 0.5

no_rep = grp_total != np.append(np.zeros(1), grp_total[:-1])
print(len(grp_peri_bool[no_rep]))
print(np.sum(grp_peri_bool[no_rep]))
print(np.sum(grp_peri_bool[no_rep] / len(grp_peri_bool[no_rep])))
print('')



bln_bnd = (final_s[:, 1]==128) * (final_s[:, 0]==128) * (final_s_in[:, 2]==1) * grp_peri_bool
bln_ubnd = (final_s[:, 1]==128) * (final_s[:, 0]==128) * (final_s_in[:, 2]==0) * grp_peri_bool
bln_dest = (final_s[:, 1]==128) * (final_s[:, 0]<128) * grp_peri_bool
print(np.sum(bln_ubnd))
print(np.sum(bln_bnd))
print(np.sum(bln_dest))

rs_t_b, vs_t_b = rs_total[clus_bl][bln_bnd], vs_total[clus_bl][bln_bnd]
rs_t_u, vs_t_u = rs_total[clus_bl][bln_ubnd], vs_total[clus_bl][bln_ubnd]
rs_t_d, vs_t_d = rs_total[clus_bl][bln_dest], vs_total[clus_bl][bln_dest]

bins = np.linspace(0., 3., 31)
bins2d = np.meshgrid(0.5*(bins[:-1]+bins[1:]), 0.5*(bins[:-1]+bins[1:]))
bins2d_x = bins2d[1]
bins2d_y = bins2d[0]
plt.figure()
grid_b = plt.hist2d(rs_t_b, vs_t_b, bins=bins)[0]
grid_u = plt.hist2d(rs_t_u, vs_t_u, bins=bins)[0]
grid_d = plt.hist2d(rs_t_d, vs_t_d, bins=bins)[0]

where = np.where(bins2d_y < bound_crit(bins2d_x))
where = np.where((grid_d>grid_u) * (grid_d>grid_b))
bound = np.zeros((len(grid_b), len(grid_b[0])))
bound[where] = 1. 
where = np.where((grid_b>grid_u) * (grid_b>grid_d))
bound[where] = 2.
where = np.where((grid_u>grid_b) * (grid_u>grid_d))
bound[where] = 3. #1:dest, 2:bound, 3:unbound

plt.close()
fig = plt.figure(figsize=(8, 6))
ax = plt.subplot()
cax = ax.imshow(np.transpose(bound), origin='lower', extent=[0., 3., 0., 3.], 
        cmap='rainbow')
cbar = fig.colorbar(cax, ticks=[0., 1., 2., 3.], boundaries=np.arange(5)-0.5)
cbar.ax.set_yticklabels(['None', 'Destroyed', 'Bound', 'Unbound'])
plt.ylim(0., 2.5)
plt.xticks([0., 0.5, 1., 1.5, 2., 2.5, 3.])
plt.yticks([0., 0.5, 1., 1.5, 2., 2.5])
plt.xlabel(r'$r/R_{200,\rm{group}}$')
plt.ylabel(r'$v/v_{\rm{crit}}$')
#plt.savefig('heatmap.png', dpi=500)
plt.show()


plt.close()
plt.figure(figsize=(7, 6))
plt.scatter(rs_t_d, vs_t_d, s=2., c='b')
plt.scatter(rs_t_b, vs_t_b, s=2., c='y')
plt.scatter(rs_t_u, vs_t_u, s=2., c='r')
plt.xlim(0., 3.)
plt.ylim(0., 2.5)
plt.xticks([0., 0.5, 1., 1.5, 2., 2.5, 3.])
plt.yticks([0., 0.5, 1., 1.5, 2., 2.5])
plt.xlabel(r'$r/R_{200,\rm{group}}$')
plt.ylabel(r'$v/v_{\rm{crit}}$')
plt.show()
#plt.savefig('scatter.png', dpi=500)

print(a)










boolean = (final_s[:, 1]==128) * (final_s[:, 0]<128)
boolean = (final_s[:, 1]==128)# * (final_s[:, 0]==128) * (final_s_in[:, 2]==1)
rs_t, vs_t = rs_total[clus_bl][boolean], vs_total[clus_bl][boolean]
total = total[clus_bl]
rs_f_t, vs_f_t = np.zeros(0), np.zeros(0)
rs_min_t, vs_min_t = np.zeros(0), np.zeros(0)
total_infall_data = total[:, 2:][boolean]
memb_id = np.char.mod('%15d', total_infall_data[:, 0])
grou_id = np.char.mod('%15d', total_infall_data[:, 1])
clus_id = np.array(total_infall_data[:, 2], dtype='int')


for i in range(len(total_infall_data)):
    if clus_id[i] != clus_id[i-1]:
        hf_x = h5py.File(out_dir + 'Group_member_tracking/CLUSTER_%04d/pos'
                'ition_wrt_group/CLUSTER_%04d_member_xs_group.hdf5' % (
                clus_id[i], clus_id[i]), 'r')
        hf_y = h5py.File(out_dir + 'Group_member_tracking/CLUSTER_%04d/pos'
                'ition_wrt_group/CLUSTER_%04d_member_ys_group.hdf5' % (
                clus_id[i], clus_id[i]), 'r')
        hf_z = h5py.File(out_dir + 'Group_member_tracking/CLUSTER_%04d/pos'
                'ition_wrt_group/CLUSTER_%04d_member_zs_group.hdf5' % (
                clus_id[i], clus_id[i]), 'r')
        hf_vx = h5py.File(out_dir + 'Group_member_tracking/CLUSTER_%04d/pos'
                'ition_wrt_group/CLUSTER_%04d_member_vx_group.hdf5' % (
                clus_id[i], clus_id[i]), 'r')
        hf_vy = h5py.File(out_dir + 'Group_member_tracking/CLUSTER_%04d/pos'
                'ition_wrt_group/CLUSTER_%04d_member_vy_group.hdf5' % (
                clus_id[i], clus_id[i]), 'r')
        hf_vz = h5py.File(out_dir + 'Group_member_tracking/CLUSTER_%04d/pos'
                'ition_wrt_group/CLUSTER_%04d_member_vz_group.hdf5' % (
                clus_id[i], clus_id[i]), 'r')
        hf_r2 = h5py.File(out_dir + 'Group_member_tracking/CLUSTER_%04d/pos'
                'ition_wrt_group/CLUSTER_%04d_member_r200host_group.hdf5' % (
                clus_id[i], clus_id[i]), 'r')
        hf_m = h5py.File(out_dir + 'Group_member_tracking/CLUSTER_%04d/pos'
                'ition_wrt_group/CLUSTER_%04d_member_m200host_group.hdf5' % (
                clus_id[i], clus_id[i]), 'r')

    if clus_id[i] in crange:
        xs_i = np.array(hf_x[grou_id[i]+'/'+memb_id[i]])
        ys_i = np.array(hf_y[grou_id[i]+'/'+memb_id[i]])
        zs_i = np.array(hf_z[grou_id[i]+'/'+memb_id[i]])
        vx_i = np.array(hf_vx[grou_id[i]+'/'+memb_id[i]])
        vy_i = np.array(hf_vy[grou_id[i]+'/'+memb_id[i]])
        vz_i = np.array(hf_vz[grou_id[i]+'/'+memb_id[i]])
        r2s_i = np.array(hf_r2[grou_id[i]+'/'+memb_id[i]])
        ms_i = np.array(hf_m[grou_id[i]+'/'+memb_id[i]])

        rs_i = (xs_i**2. + ys_i**2. + zs_i**2.)**0.5 / (0.01+r2s_i)
        #rs_f_t = np.append(rs_f_t, rs_i[rs_i>0.000001][-1])
        min_ind = np.where(rs_i[rs_i>0.000001]==np.min(
                rs_i[rs_i>0.000001]))[0][0]
        rs_min_t = np.append(rs_min_t, rs_i[rs_i>0.000001][min_ind])

        ge = ((1.3271244*10.**20.) * (1.+ms_i)) / (3.0857*10.**19.)
        vs_i = (vx_i**2. + vy_i**2. + vz_i**2.)**0.5 / (
                0.001*(1.2*ge/(0.01+r2s_i))**0.5)
        #vs_f_t = np.append(vs_f_t, vs_i[rs_i>0.000001][-1])
        vs_min_t = np.append(vs_min_t, vs_i[rs_i>0.000001][min_ind])

"""
plt.figure()
bins_arr = np.linspace(0., 2.6, 1000)
#plt.hist(rs_t, bins=bins, alpha=0.5)

from scipy.stats import gaussian_kde
def plot_kde(xs, bw, bins):
    kde = gaussian_kde(xs, bw_method=bw / xs.std(ddof=1))
    return kde.evaluate(bins)

kde_infal = plot_kde(rs_t, 0.1, bins_arr)
kde_destr = plot_kde(rs_f_t, 0.1, bins_arr)
kde_peric = plot_kde(rs_min_t, 0.1, bins_arr)
plt.plot(bins_arr, kde_infal, c='k')
plt.plot(bins_arr, kde_destr, c='r')
plt.plot(bins_arr, kde_peric, c='b')
print(bins_arr[np.where(kde_infal==max(kde_infal))[0][0]])
print(bins_arr[np.where(kde_destr==max(kde_destr))[0][0]])
print(np.median(rs_t))
print(np.median(rs_f_t))
print(np.sum(rs_f_t < 0.5)/len(rs_f_t))

plt.xlim(0., 2.)
plt.ylim(0., 2.5)

#plt.hist(rs_f_t, bins=bins, alpha=0.5)
plt.show()
"""





def density_estimation(m1, m2, crop=False):
    """ Smoothed density for contour plot """
    x, y = np.linspace(xmin, xmax, steps), np.linspace(ymin, ymax, steps)
    X, Y = np.meshgrid(x, y)
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([m1, m2])

    kernel = stats.gaussian_kde(values) 
    Z = np.reshape(kernel(positions).T, X.shape)

    if crop==True:
        bnd = Y < bound_crit(X)
        Z[bnd==False] = 0.
    if crop=='upper':
        X_t = np.copy(X)
        X_t[X_t>2.5]=2.5
        bnd = Y < bound_crit(X_t)
        Z[bnd==True] = 0.

    return X, Y, Z


plt.figure(figsize=(6, 5))
#plt.scatter(rs_t, vs_t, s=1., c='k')
#plt.scatter(rs_f_t, vs_f_t, s=1., c='r')
#plt.scatter(rs_min_t, vs_min_t, s=1., c='b')

X, Y, Z = density_estimation(rs_t, vs_t, True)
plt.contour(X, Y, Z, colors='b', levels=[0.3,0.6,0.9], 
        linewidths=[1., 1.8, 2.5])
#X, Y, Z = density_estimation(rs_f_t[(rs_f_t<xmax)*(vs_f_t<ymax)], 
#        vs_f_t[(rs_f_t<xmax)*(vs_f_t<ymax)])
#plt.contour(X, Y, Z, colors='r', levels=[0.3,0.6,0.9], 
#        linewidths=[1., 1.8, 2.5])
X, Y, Z = density_estimation(rs_min_t[(rs_min_t<xmax)*(vs_min_t<ymax)], 
        vs_min_t[(rs_min_t<xmax)*(vs_min_t<ymax)])
plt.contour(X, Y, Z, colors='r', levels=[0.3,0.6,0.9], 
        linewidths=[1., 1.8, 2.5])
plt.plot((np.arange(1,2501)/1000.),bound_crit(np.arange(1,2501)/1000.), 
        c='k', linewidth=2.)
plt.plot([-1., -1.], [-1., -1.], c='b', linewidth=2., 
        label=r'Unbound at $z=0$, at $z_{\rm{infall}}$')
plt.plot([-1., -1.], [-1., -1.], c='r', linewidth=2., 
        label=r'Unbound at $z=0$, at $z_{\rm{pericentre}}$')
plt.legend(fontsize=14)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.xlabel(r'$r/R_{200,\rm{group}}$')
plt.ylabel(r'$v/v_{\rm{crit}}$')
plt.tight_layout()
#plt.savefig('unbound_pericentre.png', dpi=500)
plt.show()