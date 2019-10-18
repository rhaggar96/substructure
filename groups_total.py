from modules import *


run_infall_finder = False
run_find_bound_groups = False
savefig = False
showfig = True
plot_contours = True
plot_contours_z0 = False
plot_scatter = False
plot_scatter_z0 = False
run_final_state_finder = False

out_dir = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/Infalling_Groups'
        '/MergerTreeAHF_Infalling_Groups/')
fig_dir = ('/home/ppxrh2/Documents/test_pollux/TheThreeHundred/playground/'
        'rhaggar/substructure_data/figures/')
track_dir = ('/home/ppxrh2/Documents/test_pollux/TheThreeHundred/playground/'
        'rhaggar/substructure_data/tracking_data/')
crange = np.array(np.loadtxt('/home/ppxrh2/Documents/test_pollux/TheThree'
        'Hundred/playground/rhaggar/G3X_data/G3X_300_selected_sample_257'
        '.txt'), dtype='int')

p_lim = 20
m_rat = 0.3
mstarlim = 10.**9.5
mlim = 10.**10.5


def find_groups(c, outdir, m_lim=0., ms_lim=0., rat_lim=1.):
    """ Find all the objects which pass within R200 of the cluster, having 
    been outside of R200 at a previous snapshot """

    loaddir = ('/run/media/ppxrh2/166AA4B87A2DD3B7/NewMDCLUSTER_data/reduced_'
            'cluster_info/')
    xs_in  = h5py.File(loaddir+'xs/CLUSTER_%04d_xs' % c)
    xs_in  = h5py.File(loaddir+'xs/CLUSTER_%04d_xs' % c)
    ys_in  = h5py.File(loaddir+'ys/CLUSTER_%04d_ys' % c)
    zs_in  = h5py.File(loaddir+'zs/CLUSTER_%04d_zs' % c)
    ms_in  = h5py.File(loaddir+'ms/CLUSTER_%04d_ms' % c)
    mstar_in  = h5py.File(loaddir+'mstars/CLUSTER_%04d_mstars' % c)
    r200_in = h5py.File(loaddir+'rvirs/CLUSTER_%04d_rvirs' % c)

    halo_ids = ld_arr('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
            'MergerTreeAHF_General_Tree_Comp/NewMDCLUSTER_%04d/snap_128/' % c
            + 'CLUSTER_%04d.txt' % c, dtype='int')[c_ids[c-1]]

    halo_ids = halo_ids[halo_ids > 0]
    keys = np.char.mod(u'%03d', halo_ids // mod)


    j = keys[0]
    xs = np.array(xs_in[j])
    ys = np.array(ys_in[j])
    zs = np.array(zs_in[j])
    ms = np.array(ms_in[j])+1.
    ids = np.array(np.arange(len(xs)), dtype='int')+int(j)*mod+1
    mstars = np.array(mstar_in[j])
    r200s = np.array(r200_in[j])
    
    h_id = halo_ids[0] - (int(j)*mod+1)
    clust = [xs[h_id], ys[h_id], zs[h_id], r200s[h_id]]
    ratio = mstars / ms

    m_bool = (ms >= m_lim) * (mstars >= ms_lim) * (ratio <= rat_lim)
    m_bool[h_id] = False
    xs, ys, zs, ids = xs[m_bool], ys[m_bool], zs[m_bool], ids[m_bool]
        
    xs, ys, zs = xs-clust[0], ys-clust[1], zs-clust[2]
    rs = (xs**2. + ys**2. + zs**2.)**0.5 / clust[3]
    ids_out = ids[rs > 1.]

    ids_infall_t = np.array(np.zeros(0), dtype='int')
    
    for i in range(1, len(keys)):
        j = keys[i]
        #print(j)
        rs_p = rs

        xs = np.array(xs_in[j])
        ys = np.array(ys_in[j])
        zs = np.array(zs_in[j])
        ms = np.array(ms_in[j])+1.
        mstars = np.array(mstar_in[j])
        r200s = np.array(r200_in[j])
        ids = np.array(np.arange(len(xs)), dtype='int')+int(j)*mod+1
        ids_main_p_in = ld_arr('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTree'
                'AHF/MergerTreeAHF_HDF5_MainProg/CLUSTER_%04d_CRMratio2_Main'
                'Prog/CLUSTER_%04d_prog_output.snapshot_' % (c, c) + j
                + '.AHF.MainProg.txt', dtype='int')
        j = int(j)
        ids_main_p = np.array(np.zeros(ids_main_p_in[-1, 0]-j*mod), 
                dtype='int')
        ids_main_p[ids_main_p_in[:, 0]-(j*mod+1)] = ids_main_p_in[:, 1]

        h_id = halo_ids[i] - (int(j)*mod+1)
        clust = [xs[h_id], ys[h_id], zs[h_id], r200s[h_id]]
        ratio = mstars / ms

        m_bool = (ms >= m_lim) * (mstars >= ms_lim) * (ratio <= rat_lim)
        m_bool[h_id] = False
        xs, ys, zs, ids = xs[m_bool], ys[m_bool], zs[m_bool], ids[m_bool]
        
        xs, ys, zs = xs-clust[0], ys-clust[1], zs-clust[2]
        rs = (xs**2. + ys**2. + zs**2.)**0.5 / clust[3]
        

        ids_infall = ids[rs < 1.]
        ids_prev = ids_main_p[ids_infall-(j*mod+1)]
        #ids_prev: ids of where haloes in r200 previously were located 
        #ids_out: haloes in ALL previous snapshots which were outside r200 
        
        for k in range(len(ids_prev)):
            if len(np.where(ids_out==ids_prev[k])[0])==0:
                ids_prev[k] = 0
        
        #ids_infall: ids of every halo which has JUST fallen into the r200
        ids_infall = ids_infall[ids_prev>0]
        ids_infall_t = np.append(ids_infall_t, ids_infall)

        ids_out = np.append(ids_out, ids[rs > 1.])
    
    total_tree = np.array(pd.read_csv('/run/media/ppxrh2/166AA4B87A2DD3B7/'
            'MergerTreeAHF/MergerTreeAHF_ASCII/MergerTree_GadgetX-NewMD'
            'CLUSTER_%04d.txt-CRMratio2' % c, sep='\s+', skiprows=2, 
            usecols=[0], dtype='str')[:-1], dtype='int')[:, 0]

    hf = h5py.File(outdir + 'cluster_%04d_groups.hdf5' % c, 'w')
    i_nonrep = 0
    for i in range(len(ids_infall_t)):
        #print(i)
        id_temp = ids_infall_t[i]
        ids_infall_t[i] = 0
        
        branch = find_branch(np.array([id_temp]), total_tree, keep_sing=False)
        n_infall = 1
        for j in branch[1:]:
            repeat = np.where(ids_infall_t == j)[0]
            if len(repeat) == 1:
                n_infall += 1
                ids_infall_t[repeat[0]] = 0
        if len(branch) > 0:
            #label dataset with index, and number of infalls
            hf.create_dataset('%04d_%02d' % (i_nonrep, n_infall), data=branch)
            i_nonrep += 1
    
    hf.close()

    return None


def find_bound_groups(c, loaddir, m_lim=0., ms_lim=0., rat_lim=1.):
    infalls = h5py.File(loaddir 
            + 'Group_halo_IDs/cluster_%04d_groups.hdf5' % c, 'r')
    keys = np.array(list(infalls.keys()), dtype='str')
    n_infalls = np.array(np.core.defchararray.partition(keys, '_')[:, 2], 
            dtype='int')
    infall_ids_t = np.array([np.array(infalls[keys[i]], dtype='int')[0] 
            for i in np.arange(len(keys))])
    
    infall_snaps = infall_ids_t // mod
    infall_ids = infall_ids_t - (infall_snaps*mod+1)
    infall_snaps_ct = np.array(list(Counter(infall_snaps).values()))
    infall_snaps = np.array(list(Counter(infall_snaps).keys()))
    
    h_info_dir = ('/run/media/ppxrh2/166AA4B87A2DD3B7/NewMDCLUSTER_data/'
            'reduced_cluster_info/')

    xs = h5py.File(h_info_dir + 'xs/CLUSTER_%04d_xs' % c, 'r')
    ys = h5py.File(h_info_dir + 'ys/CLUSTER_%04d_ys' % c, 'r')
    zs = h5py.File(h_info_dir + 'zs/CLUSTER_%04d_zs' % c, 'r')
    vx = h5py.File(h_info_dir + 'vx/CLUSTER_%04d_vx' % c, 'r')
    vy = h5py.File(h_info_dir + 'vy/CLUSTER_%04d_vy' % c, 'r')
    vz = h5py.File(h_info_dir + 'vz/CLUSTER_%04d_vz' % c, 'r')
    ms = h5py.File(h_info_dir + 'ms/CLUSTER_%04d_ms' % c, 'r')
    mstars = h5py.File(h_info_dir + 'mstars/CLUSTER_%04d_mstars' % c, 'r')
    r200s = h5py.File(h_info_dir + 'rvirs/CLUSTER_%04d_rvirs' % c, 'r')

    rs_out, vs_out = np.zeros(0), np.zeros(0) 
    ids_out = np.array(np.zeros(0), dtype='int')
    hosts_out = np.array(np.zeros(0), dtype='int')
    counter = 0
    bound_out = h5py.File(loaddir
            + 'Group_member_IDs/cluster_%04d_members.hdf5' % c, 'w')
    for s in range(len(infall_snaps)):
        snap = '%03d' % infall_snaps[s]
        xs_s, ys_s, zs_s, r200s_s = xs[snap], ys[snap], zs[snap], r200s[snap]
        vx_s, vy_s, vz_s = vx[snap], vy[snap], vz[snap]
        ms_s = np.array(ms[snap], dtype='int') 
        mstars_s = np.array(mstars[snap], dtype='int')
        ratio = (mstars_s+1) / ms_s
        ids_s = np.array(np.arange(len(xs_s)), dtype='int')

        boolean = (ms_s >= m_lim) * (mstars_s >= ms_lim) * (ratio <= rat_lim)
        xs_s, ys_s, zs_s = xs_s[boolean], ys_s[boolean], zs_s[boolean]
        vx_s, vy_s, vz_s = vx_s[boolean], vy_s[boolean], vz_s[boolean]
        ms_s, mstars_s = ms_s[boolean], mstars_s[boolean]
        ids_s, r200s_s = ids_s[boolean], r200s_s[boolean]
        
        for i in range(infall_snaps_ct[s]):
            id_n = infall_ids[counter]
            id_n = np.where(ids_s == id_n)[0][0]
            rs_si = ((xs_s-xs_s[id_n])**2. + (ys_s-ys_s[id_n])**2.
                    + (zs_s-zs_s[id_n])**2.)**0.5
            vs_si = ((vx_s-vx_s[id_n])**2. + (vy_s-vy_s[id_n])**2.
                    + (vz_s-vz_s[id_n])**2.)**0.5
            ms_si = (ms_s < ms_s[id_n])
            ids_si, rs_si, vs_si = ids_s[ms_si], rs_si[ms_si], vs_si[ms_si]
            boundres = bound(vs_si, rs_si, ms_s[id_n], r200s_s[id_n])
            rs_out = np.append(rs_out, boundres[1][boundres[0]])
            vs_out = np.append(vs_out, boundres[2][boundres[0]])
            ids_bound = infall_snaps[s]*mod+1+ids_si[boundres[0]]
            ids_out = np.append(ids_out, ids_bound)
            hosts_out = np.append(hosts_out, np.array(
                    [infall_ids_t[counter]]*len(ids_bound)))
            bound_out.create_dataset('%16d_%03d' % (infall_ids_t[counter], 
                    len(ids_bound)), data=ids_bound)
            counter += 1
    bound_out.close()
    
    return rs_out, vs_out, ids_out, hosts_out


def find_member_final_states(c, grp_memb_dat):
    """ Find the final ID of infalling groups and each of their 
    members, and for those which both survive, find the relative
    r and v. No limits on mass, etc. after z_infall """
    total_tree = np.array(pd.read_csv('/run/media/ppxrh2/166AA4B87A2DD3B7/'
            'MergerTreeAHF/MergerTreeAHF_ASCII/MergerTree_GadgetX-NewMD'
            'CLUSTER_%04d.txt-CRMratio2' % c, sep='\s+', skiprows=2, 
            usecols=[0], dtype='str')[:-1], dtype='int')[:, 0]

    grp_memb_dat = grp_memb_dat[np.array(grp_memb_dat[:, 4], dtype='int')==c]
    memb_init_id = np.array(grp_memb_dat[:, 2], dtype='int')
    host_init_id = np.array(grp_memb_dat[:, 3], dtype='int')

    memb_final_id = np.array(np.zeros(0), dtype='int')
    host_final_id = np.array(np.zeros(0), dtype='int')

    for i in range(len(memb_init_id)):
        memb_final_id = np.append(memb_final_id, find_branch(np.array([
                memb_init_id[i]]), total_tree, keep_sing=True)[-1])
        if host_init_id[i] != host_init_id[i-1]:
            host_final_id = np.append(host_final_id, find_branch(np.array([
                host_init_id[i]]), total_tree, keep_sing=True)[-1])
        else:
            host_final_id = np.append(host_final_id, host_final_id[-1])

    memb_snaps, host_snaps = memb_final_id//mod, host_final_id//mod
    memb_ids = memb_final_id-(mod*memb_snaps+1)
    host_ids = host_final_id-(mod*host_snaps+1)
    survive = (memb_snaps==host_snaps)
    survive = np.where(survive==True)[0]
    
    h_info_dir = ('/run/media/ppxrh2/166AA4B87A2DD3B7/NewMDCLUSTER_data/'
            'reduced_cluster_info/')
    xs = h5py.File(h_info_dir + 'xs/CLUSTER_%04d_xs' % c, 'r')
    ys = h5py.File(h_info_dir + 'ys/CLUSTER_%04d_ys' % c, 'r')
    zs = h5py.File(h_info_dir + 'zs/CLUSTER_%04d_zs' % c, 'r')
    vx = h5py.File(h_info_dir + 'vx/CLUSTER_%04d_vx' % c, 'r')
    vy = h5py.File(h_info_dir + 'vy/CLUSTER_%04d_vy' % c, 'r')
    vz = h5py.File(h_info_dir + 'vz/CLUSTER_%04d_vz' % c, 'r')
    ms = h5py.File(h_info_dir + 'ms/CLUSTER_%04d_ms' % c, 'r')
    r200s = h5py.File(h_info_dir + 'rvirs/CLUSTER_%04d_rvirs' % c, 'r')

    rs_f, vs_f = np.zeros(len(memb_final_id)), np.zeros(len(memb_final_id))
    bound_f = np.array(np.zeros(len(memb_final_id)), dtype='int')
    for i in survive:
        s = memb_snaps[i]
        xs_s = np.array(xs['%03d'%s])
        ys_s = np.array(ys['%03d'%s])
        zs_s = np.array(zs['%03d'%s])
        r_si = ((xs_s[memb_ids[i]]-xs_s[host_ids[i]])**2.
                + (ys_s[memb_ids[i]]-ys_s[host_ids[i]])**2.
                + (zs_s[memb_ids[i]]-zs_s[host_ids[i]])**2.)**0.5
        vx_s = np.array(vx['%03d'%s])
        vy_s = np.array(vy['%03d'%s])
        vz_s = np.array(vz['%03d'%s])
        v_si = ((vx_s[memb_ids[i]]-vx_s[host_ids[i]])**2.
                + (vy_s[memb_ids[i]]-vy_s[host_ids[i]])**2.
                + (vz_s[memb_ids[i]]-vz_s[host_ids[i]])**2.)**0.5
        ms_s, r200s_s = np.array(ms['%03d'%s]), np.array(r200s['%03d'%s])
        bndres = bound(np.array([v_si]), np.array([r_si]), ms_s[host_ids[i]], 
                r200s_s[host_ids[i]])
        rs_f[i], vs_f[i] = bndres[1], bndres[2]
        bound_f[i] += 1*bndres[0]

    return memb_final_id, host_final_id, bound_f>0, rs_f, vs_f


def find_branch(id_a, id_list, keep_sing=False):
    id = id_a[-1]
    tree_search = np.where(id_list==id)[0]
    if len(tree_search)==0:
        result = np.array([id]*keep_sing, dtype='int')
    elif id//mod == id_list[tree_search[0]-1]//mod:
        result = id_a
    elif len(tree_search)==1:
        result = id_a
    else:
        id_a = np.append(id_a, id_list[tree_search[0]-1])
        result = find_branch(id_a, id_list)
    
    return result


def bound(vrels, rrel, m, r200, incgroup=False):
    """ Determine whether a halo is bound to its host """
    rrels = rrel
    host = rrels<0.001
    rrels[host] = 0.001
    ke = 0.5 * (1000.*vrels)**2.
    ge = ((1.3271244*10.**20.) * m) / (3.0857*10.**19.)
    vir = (ke - (ge / rrels)) < (ge / (-2.5*r200))
    v_cr = 0.001*(1.2*ge/r200)**0.5
    if incgroup==False:
        vir[host]=False
    return vir, rrels/r200, vrels/v_cr


def bound_crit(x):
    """ y = v/vcrit, x = r/R200, find boundness condition """
    x[x>2.5]=2.5
    x[x<0.000001] = 0.000001
    return ((5./(3.*x)) - (2./3.))**0.5


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



if run_infall_finder==True:
    for c_id in range(1, 325):
        print(c_id)
        find_groups(c_id, out_dir, mlim, mstarlim, m_rat)

rs_total, vs_total = np.zeros(0), np.zeros(0)
id_total, hs_total = np.zeros(0), np.zeros(0)
clus_nos = np.array(np.zeros(0), dtype='int')
if run_find_bound_groups==True:
    for c_id in crange:
        print(c_id)
        res = find_bound_groups(c_id, out_dir, mlim, mstarlim, m_rat)
        rs_total = np.append(rs_total, res[0])
        vs_total = np.append(vs_total, res[1])
        id_total = np.append(id_total, res[2])
        hs_total = np.append(hs_total, res[3])
        clus_nos = np.append(clus_nos, c_id*np.ones(len(res[0])))
    pd.DataFrame(np.transpose(np.array([np.char.mod('%12.9f', rs_total), 
            np.char.mod('%12.9f', vs_total), np.char.mod('%15d', id_total), 
            np.char.mod('%15d', hs_total), np.char.mod('%04d', clus_nos)]))
            ).to_csv(track_dir+'forward_tracking.txt', index=None, sep='\t', 
            header=['    rs_total', '    vs_total', 'ids_bound_membs', 
            '  id_group_halo', 'clus'])
total = ld_arr(track_dir+'forward_tracking.txt')
rs_total, vs_total, id_total = total[:, 0], total[:, 1], total[:, 2]//mod


steps=100
if plot_contours==True:
    xmin, xmax, ymin, ymax = 0., 2.6, 0., 2.5
    final_s_in = ld_arr(track_dir+'forward_tracking_finalsnap.txt')
    final_s = np.array(final_s_in[:, :2], dtype='int')//mod
    final_s_bd = np.array(final_s_in[:, 2], dtype='int')

    boolean = (final_s[:, 1]==128) * (final_s[:, 0]==128)
    rs_t, vs_t = rs_total[boolean], vs_total[boolean] 
    fnl_s_bd = final_s_bd[boolean]
    plt.figure(figsize=(14, 6))

    plt.subplot(121)
    X, Y, Z = density_estimation(rs_t, vs_t, True)
    plt.contour(X, Y, Z, colors='g', levels=[0.3,0.5,0.7,0.9], 
            linewidths=[1., 1.5, 2., 2.5])
    plt.plot([-1,-1], [-1,-1], c='g', label='All surviving galaxies')
    plt.legend()
    plt.plot((np.arange(1,2501)/1000.),bound_crit(np.arange(1,2501)/1000.), 
            c='k', linewidth=2.)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel(r'$r_{\rm{infall}}/R_{\rm{200,group}}$')
    plt.ylabel(r'$v_{\rm{infall}}/v_{\rm{crit}}$')

    plt.subplot(122)
    X, Y, Z = density_estimation(rs_t[fnl_s_bd>0], 
            vs_t[fnl_s_bd>0], True)
    plt.contour(X, Y, Z, colors='r', levels=[0.3,0.5,0.7,0.9], 
            linewidths=[1., 1.5, 2., 2.5], linestyles='-.')
    X, Y, Z = density_estimation(rs_t[fnl_s_bd==0], 
            vs_t[fnl_s_bd==0], True)
    plt.contour(X, Y, Z, colors='b', levels=[0.3,0.5,0.7,0.9], 
            linewidths=[1., 1.5, 2., 2.5])

    plt.plot([-1,-1], [-1,-1], c='r', label='Bound at $z=0$', linestyle='-.')
    plt.plot([-1,-1], [-1,-1], c='b', label='Unbound at $z=0$')
    plt.legend()
    plt.plot((np.arange(1,2501)/1000.),bound_crit(np.arange(1,2501)/1000.), 
            c='k', linewidth=2.)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel(r'$r_{\rm{infall}}/R_{\rm{200,group}}$')
    plt.ylabel(r'$v_{\rm{infall}}/v_{\rm{crit}}$')

    plt.tight_layout()
    if savefig==True:   
        plt.savefig(figdir+'survivors_infall_params_con.png', dpi=500)
        plt.savefig(figdir+'survivors_infall_params_con.pdf')
    if showfig==True:
        plt.show()
    plt.close()
    

    plt.figure(figsize=(7, 6))
    boolean = (final_s[:, 1]==128) * (final_s[:, 0]<128)
    rs_t, vs_t = rs_total[boolean], vs_total[boolean] 
    X, Y, Z = density_estimation(rs_t, vs_t, True)
    plt.contour(X, Y, Z, colors='k', levels=[0.3,0.5,0.7,0.9], 
            linewidths=[1., 1.5, 2., 2.5])
    plt.plot([-1,-1], [-1,-1], c='k', label='Destroyed at $z=0$')
    plt.legend()
    plt.plot((np.arange(1,2501)/1000.),bound_crit(np.arange(1,2501)/1000.), 
            c='k', linewidth=2.)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel(r'$r_{\rm{infall}}/R_{\rm{200,group}}$')
    plt.ylabel(r'$v_{\rm{infall}}/v_{\rm{crit}}$')

    if savefig==True:   
        plt.savefig(figdir+'destroyed_infall_params_con.png', dpi=500)
        plt.savefig(figdir+'destroyed_infall_params_con.pdf')
    if showfig==True:
        plt.show()
    plt.close()

if plot_scatter==True:
    xmin, xmax, ymin, ymax = 0., 2.6, 0., 2.5
    final_s_in = ld_arr(track_dir+'forward_tracking_finalsnap.txt')
    final_s = np.array(final_s_in[:, :2], dtype='int')//mod
    final_s_bd = np.array(final_s_in[:, 2], dtype='int')

    boolean = (final_s[:, 1]==128) * (final_s[:, 0]==128)
    rs_t, vs_t = rs_total[boolean], vs_total[boolean] 
    fnl_s_bd = final_s_bd[boolean]
    plt.figure(figsize=(14, 6))

    plt.subplot(121)
    plt.scatter(rs_t, vs_t, s=2., c='g', label='All surviving galaxies')
    plt.legend()
    plt.plot((np.arange(1,2501)/1000.),bound_crit(np.arange(1,2501)/1000.), 
            c='k', linewidth=2.)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel(r'$r_{\rm{infall}}/R_{\rm{200,group}}$')
    plt.ylabel(r'$v_{\rm{infall}}/v_{\rm{crit}}$')

    plt.subplot(122)
    plt.scatter(rs_t[fnl_s_bd>0], vs_t[fnl_s_bd>0], s=2., c='r', 
            label='Bound at $z=0$')
    plt.scatter(rs_t[fnl_s_bd==0], vs_t[fnl_s_bd==0], s=8., c='b', 
            marker='x', linewidth=0.8, label='Unbound at $z=0$')
    plt.legend()
    plt.plot((np.arange(1,2501)/1000.),bound_crit(np.arange(1,2501)/1000.), 
            c='k', linewidth=2.)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel(r'$r_{\rm{infall}}/R_{\rm{200,group}}$')
    plt.ylabel(r'$v_{\rm{infall}}/v_{\rm{crit}}$')

    plt.tight_layout()
    if savefig==True:   
        plt.savefig(figdir+'survivors_infall_params_sca.png', dpi=500)
        plt.savefig(figdir+'survivors_infall_params_sca.pdf')
    if showfig==True:
        plt.show()
    plt.close()
    

    boolean = (final_s[:, 1]==128) * (final_s[:, 0]<128)
    rs_t, vs_t = rs_total[boolean], vs_total[boolean] 
    plt.figure(figsize=(7, 6))

    plt.scatter(rs_t, vs_t, s=2., c='k', label='Destroyed at $z=0$')
    plt.legend()
    plt.plot((np.arange(1,2501)/1000.),bound_crit(np.arange(1,2501)/1000.), 
            c='k', linewidth=2.)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel(r'$r_{\rm{infall}}/R_{\rm{200,group}}$')
    plt.ylabel(r'$v_{\rm{infall}}/v_{\rm{crit}}$')

    plt.tight_layout()
    if savefig==True:   
        plt.savefig(figdir+'destroyed_infall_params_sca.png', dpi=500)
        plt.savefig(figdir+'destroyed_infall_params_sca.pdf')
    if showfig==True:
        plt.show()
    plt.close()

if plot_contours_z0==True:
    xmin, xmax, ymin, ymax = 0., 6., 0., 6.
    final_s_in = ld_arr(track_dir+'forward_tracking_finalsnap.txt')
    rs_total, vs_total = final_s_in[:, 3], final_s_in[:, 4]
    final_s = np.array(final_s_in[:, :2], dtype='int')//mod
    final_s_bd = np.array(final_s_in[:, 2], dtype='int')

    boolean = (final_s[:, 1]==128) * (final_s[:, 0]==128)
    rs_t, vs_t = rs_total[boolean], vs_total[boolean] 
    fnl_s_bd = final_s_bd[boolean]
    plt.figure(figsize=(14, 6))

    plt.subplot(121)
    boolean = (rs_t < xmax) * (vs_t < ymax)
    rs_t_all, vs_t_all = rs_t[boolean], vs_t[boolean]
    X, Y, Z = density_estimation(rs_t_all, vs_t_all, False)
    plt.contour(X, Y, Z, colors='g', levels=[0.01, 0.04, 0.1, 0.25], 
            linewidths=[1., 1.5, 2., 2.5])
    plt.plot([-1,-1], [-1,-1], c='g', label='All surviving galaxies')
    plt.legend()
    plt.plot((np.arange(1,2501)/1000.),bound_crit(np.arange(1,2501)/1000.), 
            c='k', linewidth=2.)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel(r'$r_{z=0}/R_{\rm{200,group}}$')
    plt.ylabel(r'$v_{z=0}/v_{\rm{crit}}$')

    plt.subplot(122)
    X, Y, Z = density_estimation(rs_t[fnl_s_bd>0], 
            vs_t[fnl_s_bd>0], True)
    plt.contour(X, Y, Z, colors='r', levels=[0.3,0.5,0.7,0.9], 
            linewidths=[1., 1.5, 2., 2.5], linestyles='-.')
    X, Y, Z = density_estimation(rs_t[boolean*(fnl_s_bd==0)], 
            vs_t[boolean*(fnl_s_bd==0)], 'upper')
    plt.contour(X, Y, Z, colors='b', levels=[0.015, 0.05, 0.1, 0.18], 
            linewidths=[1., 1.5, 2., 2.5])

    plt.plot([-1,-1], [-1,-1], c='r', label='Bound at $z=0$', linestyle='-.')
    plt.plot([-1,-1], [-1,-1], c='b', label='Unbound at $z=0$')
    plt.legend()
    plt.plot((np.arange(1,2501)/1000.),bound_crit(np.arange(1,2501)/1000.), 
            c='k', linewidth=2.)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel(r'$r_{z=0}/R_{\rm{200,group}}$')
    plt.ylabel(r'$v_{z=0}/v_{\rm{crit}}$')

    plt.tight_layout()
    if savefig==True:   
        plt.savefig(figdir+'survivors_z=0_params_con.png', dpi=500)
        plt.savefig(figdir+'survivors_z=0_params_con.pdf')
    if showfig==True:
        plt.show()
    plt.close()

if plot_scatter_z0==True:
    xmin, xmax, ymin, ymax = 0., 6., 0., 6.
    final_s_in = ld_arr(track_dir+'forward_tracking_finalsnap.txt')
    rs_total, vs_total = final_s_in[:, 3], final_s_in[:, 4]
    final_s = np.array(final_s_in[:, :2], dtype='int')//mod
    final_s_bd = np.array(final_s_in[:, 2], dtype='int')

    boolean = (final_s[:, 1]==128) * (final_s[:, 0]==128)
    rs_t, vs_t = rs_total[boolean], vs_total[boolean] 
    fnl_s_bd = final_s_bd[boolean]
    plt.figure(figsize=(14, 6))

    plt.subplot(121)
    boolean = (rs_t < xmax) * (vs_t < ymax)
    rs_t_all, vs_t_all = rs_t[boolean], vs_t[boolean]
    plt.scatter(rs_t_all, vs_t_all, s=2., c='g', 
            label='All surviving galaxies')
    plt.legend(frameon=True, edgecolor='k', framealpha=1.)
    plt.plot((np.arange(1,2501)/1000.),bound_crit(np.arange(1,2501)/1000.), 
            c='k', linewidth=2.)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel(r'$r_{z=0}/R_{\rm{200,group}}$')
    plt.ylabel(r'$v_{z=0}/v_{\rm{crit}}$')

    plt.subplot(122)
    plt.scatter(rs_t[fnl_s_bd>0], vs_t[fnl_s_bd>0], s=2., c='r', 
            label='Bound at $z=0$')
    plt.scatter(rs_t[fnl_s_bd==0], vs_t[fnl_s_bd==0], s=8., c='b', 
            marker='x', linewidth=0.8, label='Unbound at $z=0$')
    plt.legend(frameon=True, edgecolor='k', framealpha=1.)
    plt.plot((np.arange(1,2501)/1000.),bound_crit(np.arange(1,2501)/1000.), 
            c='k', linewidth=2.)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel(r'$r_{z=0}/R_{\rm{200,group}}$')
    plt.ylabel(r'$v_{z=0}/v_{\rm{crit}}$')

    plt.tight_layout()
    if savefig==True:   
        plt.savefig(figdir+'survivors_z=0_params_sca.png', dpi=500)
        plt.savefig(figdir+'survivors_z=0_params_sca.pdf')
    if showfig==True:
        plt.show()
    plt.close()
    

if run_final_state_finder==True:
    final_id = np.array(np.zeros(0), dtype='int')
    final_hs = np.array(np.zeros(0), dtype='int')
    final_bd = np.array(np.zeros(0), dtype='int')
    final_rs = np.array(np.zeros(0), dtype='int')
    final_vs = np.array(np.zeros(0), dtype='int')
    for c_id in crange:
        print(c_id)
        final = find_member_final_states(c_id, total)
        final_id = np.append(final_id, final[0])
        final_hs = np.append(final_hs, final[1])
        final_bd = np.append(final_bd, final[2])
        final_rs = np.append(final_rs, final[3])
        final_vs = np.append(final_vs, final[4])
    pd.DataFrame(np.transpose(np.array([np.char.mod('%15d', final_id), 
            np.char.mod('%15d', final_hs), final_bd, np.char.mod('%12.9f', 
            final_rs), np.char.mod('%12.9f', final_vs)]))
            ).to_csv(track_dir+'forward_tracking_finalsnap.txt', index=None, 
            sep='\t', header=['final_ids_membs', 'final_ids_group', 'bound', 
            '    rs_final', '    vs_final'])


