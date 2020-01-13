from modules import *


run_find_member_trees = False
run_find_member_data = False
run_find_boundness = False
run_tracking_wrt_group = False
run_tracking_wrt_cluster = False


out_dir = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/Infalling_Groups'
        '/MergerTreeAHF_Infalling_Groups/')
data_dir = ('/run/media/ppxrh2/166AA4B87A2DD3B7/NewMDCLUSTER_data/reduced_'
        'cluster_info/')
track_dir = ('/home/ppxrh2/Documents/test_pollux/TheThreeHundred/playground/'
        'rhaggar/substructure_data/tracking_data/')
cluster_ids = ld_arr('/home/ppxrh2/Documents/test_pollux/TheThreeHundred/'
        'playground/rhaggar/G3X_data/G3X_300_host_ids.txt', dtype='int')

crange = np.array(np.loadtxt('/home/ppxrh2/Documents/test_pollux/TheThree'
        'Hundred/playground/rhaggar/G3X_data/G3X_300_selected_sample_257'
        '.txt'), dtype='int')
crange = np.array(np.arange(1, 325), dtype='int')


total = ld_arr(track_dir+'forward_tracking.txt')
rs_total, vs_total, id_total = total[:, 0], total[:, 1], total[:, 2]//mod


def find_member_trees(c, outdir, grp_memb_dat):
    """ Find total branch for group members """
    total_tree = np.array(pd.read_csv('/run/media/ppxrh2/166AA4B87A2DD3B7/'
            'MergerTreeAHF/MergerTreeAHF_ASCII/MergerTree_GadgetX-NewMD'
            'CLUSTER_%04d.txt-CRMratio2' % c, sep='\s+', skiprows=2, 
            usecols=[0], dtype='str')[:-1], dtype='int')[:, 0]

    grp_memb_dat = grp_memb_dat[np.array(grp_memb_dat[:, 4], dtype='int')==c]
    memb_init_id = np.array(grp_memb_dat[:, 2], dtype='int')
    host_init_id = np.array(grp_memb_dat[:, 3], dtype='int')

    memb_final_id = np.array(np.zeros(0), dtype='int')
    host_final_id = np.array(np.zeros(0), dtype='int')

    if not os.path.exists(outdir):
        os.mkdir(outdir)
    hf = h5py.File(outdir + 'CLUSTER_%04d_member_branches.hdf5' % c, 'w')

    for i in range(len(memb_init_id)):
        memb_final_id = np.append(memb_final_id, find_branch(np.array([
                memb_init_id[i]]), total_tree, keep_sing=True)[-1])
        
        if host_init_id[i] != host_init_id[i-1] or len(host_init_id)==1:
            hf.create_dataset('%15d/%15d' % (host_init_id[i], 0), 
                    data=find_branch(np.array([host_init_id[i]]), total_tree,
                    keep_sing=True))
        hf.create_dataset('%15d/%15d' % (host_init_id[i],memb_init_id[i]),
                data=find_branch(np.array([memb_init_id[i]]), total_tree, 
                keep_sing=True))
    hf.close()
        
    return None



def find_member_data(c, outdir, datadir):
    """ Find position, velocity, etc. of haloes in groups """

    ids_g = h5py.File(outdir+'CLUSTER_%04d_member_branches.hdf5' % c, 'r')
    keys = list(ids_g.keys())

    xs = h5py.File(datadir+'xs/CLUSTER_%04d_xs' % c, 'r')
    ys = h5py.File(datadir+'ys/CLUSTER_%04d_ys' % c, 'r')
    zs = h5py.File(datadir+'zs/CLUSTER_%04d_zs' % c, 'r')
    vx = h5py.File(datadir+'vx/CLUSTER_%04d_vx' % c, 'r')
    vy = h5py.File(datadir+'vy/CLUSTER_%04d_vy' % c, 'r')
    vz = h5py.File(datadir+'vz/CLUSTER_%04d_vz' % c, 'r')
    npars = h5py.File(datadir+'npars/CLUSTER_%04d_npars' % c, 'r')
    ms = h5py.File(datadir+'ms/CLUSTER_%04d_ms' % c, 'r')
    mstars = h5py.File(datadir+'mstars/CLUSTER_%04d_mstars' % c, 'r')
    rvirs = h5py.File(datadir+'rvirs/CLUSTER_%04d_rvirs' % c, 'r')

    xs_g = h5py.File(outdir+'CLUSTER_%04d_member_xs.hdf5' % c, 'w')
    ys_g = h5py.File(outdir+'CLUSTER_%04d_member_ys.hdf5' % c, 'w')
    zs_g = h5py.File(outdir+'CLUSTER_%04d_member_zs.hdf5' % c, 'w')
    vx_g = h5py.File(outdir+'CLUSTER_%04d_member_vx.hdf5' % c, 'w')
    vy_g = h5py.File(outdir+'CLUSTER_%04d_member_vy.hdf5' % c, 'w')
    vz_g = h5py.File(outdir+'CLUSTER_%04d_member_vz.hdf5' % c, 'w')
    npars_g = h5py.File(outdir+'CLUSTER_%04d_member_npars.hdf5' % c, 'w')
    ms_g = h5py.File(outdir+'CLUSTER_%04d_member_ms.hdf5' % c, 'w')
    mstars_g = h5py.File(outdir+'CLUSTER_%04d_member_mstars.hdf5' % c, 'w')
    r200s_g = h5py.File(outdir+'CLUSTER_%04d_member_r200s.hdf5' % c, 'w')


    for i in range(len(keys)):
        keys_i = list(ids_g[keys[i]].keys())
        for j in range(len(keys_i)):
            ids_ij = np.array(ids_g[keys[i]][keys_i[j]], dtype='int')
            snaps_ij = ids_ij // mod
            hid_ij = ids_ij - (snaps_ij*mod + 1)
            snaps_ij = np.char.mod('%03d', snaps_ij)

            n_sna = len(ids_ij)
            xs_ij, ys_ij = np.zeros(n_sna), np.zeros(n_sna) 
            zs_ij, vx_ij = np.zeros(n_sna), np.zeros(n_sna)
            vy_ij, vz_ij = np.zeros(n_sna), np.zeros(n_sna)
            npars_ij = np.zeros(n_sna)
            ms_ij = np.zeros(n_sna)
            mstars_ij, r200s_ij = np.zeros(n_sna), np.zeros(n_sna)
            for k in range(n_sna):
                xs_ij[k] = xs[snaps_ij[k]][hid_ij[k]]
                ys_ij[k] = ys[snaps_ij[k]][hid_ij[k]]
                zs_ij[k] = zs[snaps_ij[k]][hid_ij[k]]
                vx_ij[k] = vx[snaps_ij[k]][hid_ij[k]]
                vy_ij[k] = vy[snaps_ij[k]][hid_ij[k]]
                vz_ij[k] = vz[snaps_ij[k]][hid_ij[k]]
                npars_ij[k] = npars[snaps_ij[k]][hid_ij[k]]
                ms_ij[k] = ms[snaps_ij[k]][hid_ij[k]]
                mstars_ij[k] = mstars[snaps_ij[k]][hid_ij[k]]
                r200s_ij[k] = rvirs[snaps_ij[k]][hid_ij[k]]
            xs_g.create_dataset(keys[i]+'/'+keys_i[j], data=xs_ij)
            ys_g.create_dataset(keys[i]+'/'+keys_i[j], data=ys_ij)
            zs_g.create_dataset(keys[i]+'/'+keys_i[j], data=zs_ij)
            vx_g.create_dataset(keys[i]+'/'+keys_i[j], data=vx_ij)
            vy_g.create_dataset(keys[i]+'/'+keys_i[j], data=vy_ij)
            vz_g.create_dataset(keys[i]+'/'+keys_i[j], data=vz_ij)
            npars_g.create_dataset(keys[i]+'/'+keys_i[j], data=npars_ij)
            ms_g.create_dataset(keys[i]+'/'+keys_i[j], data=ms_ij)
            mstars_g.create_dataset(keys[i]+'/'+keys_i[j], data=mstars_ij)
            r200s_g.create_dataset(keys[i]+'/'+keys_i[j], data=r200s_ij)
    
    return None



def find_boundness(c, outdir):
    """ Produce boolean of what objects are bound """

    ids = h5py.File(outdir+'CLUSTER_%04d_member_branches.hdf5' % c, 'r')
    xs = h5py.File(outdir+'CLUSTER_%04d_member_xs.hdf5' % c, 'r')
    ys = h5py.File(outdir+'CLUSTER_%04d_member_ys.hdf5' % c, 'r')
    zs = h5py.File(outdir+'CLUSTER_%04d_member_zs.hdf5' % c, 'r')
    vx = h5py.File(outdir+'CLUSTER_%04d_member_vx.hdf5' % c, 'r')
    vy = h5py.File(outdir+'CLUSTER_%04d_member_vy.hdf5' % c, 'r')
    vz = h5py.File(outdir+'CLUSTER_%04d_member_vz.hdf5' % c, 'r')
    ms = h5py.File(outdir+'CLUSTER_%04d_member_ms.hdf5' % c, 'r')
    mstars = h5py.File(outdir+'CLUSTER_%04d_member_mstars.hdf5' % c, 'r')
    r200s = h5py.File(outdir+'CLUSTER_%04d_member_r200s.hdf5' % c, 'r')

    keys = list(xs.keys())
    bounded = h5py.File(outdir+'CLUSTER_%04d_member_bounded.hdf5' % c, 'w')
    
    for i in range(len(keys)):
        keys_i = list(ids[keys[i]].keys())
        ids_h = np.array(ids[keys[i]][keys_i[0]])//mod
        bounded.create_dataset(keys[i]+'/'+keys_i[0], data=np.array(
                np.ones(len(ids_h)), dtype='int'))
        ids_h -= ids_h[0] #snapshots in which host is found, adjusted to zero
        xs_h, ys_h = np.zeros(ids_h[-1]+1), np.zeros(ids_h[-1]+1) 
        zs_h, r200s_h = np.zeros(ids_h[-1]+1), 0.01+np.zeros(ids_h[-1]+1)
        vx_h, vy_h = np.zeros(ids_h[-1]+1), np.zeros(ids_h[-1]+1) 
        vz_h = np.zeros(ids_h[-1]+1)
        ms_h = np.array(np.zeros(ids_h[-1]+1), dtype='int')

        #data, with spaces left for missing snapshots
        xs_h[ids_h] = np.array(xs[keys[i]][keys_i[0]])
        ys_h[ids_h] = np.array(ys[keys[i]][keys_i[0]])
        zs_h[ids_h] = np.array(zs[keys[i]][keys_i[0]])
        r200s_h[ids_h] = np.array(r200s[keys[i]][keys_i[0]])
        vx_h[ids_h] = np.array(vx[keys[i]][keys_i[0]])
        vy_h[ids_h] = np.array(vy[keys[i]][keys_i[0]])
        vz_h[ids_h] = np.array(vz[keys[i]][keys_i[0]])
        ms_h[ids_h] = np.array(ms[keys[i]][keys_i[0]])
        
        for j in range(1, len(keys_i)):
            ids_ij = np.array(ids[keys[i]][keys_i[j]])//mod
            ids_ij -= ids_ij[0] #snapshots in which host is found, adjusted to zero
            xs_ij, ys_ij = np.zeros(ids_ij[-1]+1), np.zeros(ids_ij[-1]+1) 
            zs_ij, r200s_ij = np.zeros(ids_ij[-1]+1), np.zeros(ids_ij[-1]+1)
            vx_ij, vy_ij = np.zeros(ids_ij[-1]+1), np.zeros(ids_ij[-1]+1) 
            vz_ij = np.zeros(ids_ij[-1]+1)

            #data, with spaces left for missing snapshots
            xs_ij[ids_ij] = np.array(xs[keys[i]][keys_i[j]])
            ys_ij[ids_ij] = np.array(ys[keys[i]][keys_i[j]])
            zs_ij[ids_ij] = np.array(zs[keys[i]][keys_i[j]])
            vx_ij[ids_ij] = np.array(vx[keys[i]][keys_i[j]])
            vy_ij[ids_ij] = np.array(vy[keys[i]][keys_i[j]])
            vz_ij[ids_ij] = np.array(vz[keys[i]][keys_i[j]])
            r200s_ij[ids_ij] = np.array(r200s[keys[i]][keys_i[j]])

            ldiff = len(xs_ij) - len(xs_h)
            if ldiff >= 0:
                xs_h = np.append(xs_h, np.zeros(ldiff))
                ys_h = np.append(ys_h, np.zeros(ldiff))
                zs_h = np.append(zs_h, np.zeros(ldiff))
                vx_h = np.append(vx_h, np.zeros(ldiff))
                vy_h = np.append(vy_h, np.zeros(ldiff))
                vz_h = np.append(vz_h, np.zeros(ldiff))
                r200s_h = np.append(r200s_h, 0.01+np.zeros(ldiff))
                ms_h = np.append(ms_h, 0.01+np.zeros(ldiff))

            r200s_bool = np.array([True]*len(xs_ij))
            r200s_bool[ids_ij]=False
            xs_ij = xs_ij - xs_h[:len(xs_ij)]
            ys_ij = ys_ij - ys_h[:len(xs_ij)]
            zs_ij = zs_ij - zs_h[:len(xs_ij)]
            rs_ij = (xs_ij**2. + ys_ij**2. + zs_ij**2.)**0.5
            rs_ij[rs_ij<0.000001] = 1000000.
            r200s_h_temp = r200s_h[:len(xs_ij)] - r200s_bool * (
                    r200s_h[:len(xs_ij)]-0.01)
            rs_ij = rs_ij / r200s_h_temp

            vx_ij = vx_ij - vx_h[:len(xs_ij)]
            vy_ij = vy_ij - vy_h[:len(xs_ij)]
            vz_ij = vz_ij - vz_h[:len(xs_ij)]
            vs_ij = (vx_ij**2. + vy_ij**2. + vz_ij**2.)**0.5
            bnd = bound(vs_ij, rs_ij, ms_h[:len(xs_ij)], 
                    r200s_h[:len(xs_ij)]) * 1
            bnd[rs_ij > 1000.] = -1
            #bnd:   1 == bound
            #       0 == unbound
            #      -1 == missing snapshot (host or satellite)
            bounded.create_dataset(keys[i]+'/'+keys_i[j], data=bnd)

    return None


def find_pos_relative_group(c, outdir):
    """ Produce positions of all objects in groups, wrt group 
    centre """

    if not os.path.exists(outdir+'position_wrt_group'):
        os.mkdir(outdir+'position_wrt_group')

    tracking_inf = ld_arr(track_dir+'forward_tracking.txt')
    inf_clus = np.array(tracking_inf[:, 4], dtype='int')==c
    inf_grou = np.char.mod('%15d', tracking_inf[:, 3])[inf_clus]
    inf_memb = np.char.mod('%15d', tracking_inf[:, 2])[inf_clus]

    ids = h5py.File(outdir + 'CLUSTER_%04d_member_branches.hdf5' % c, 'r')
    xs = h5py.File(outdir + 'CLUSTER_%04d_member_xs.hdf5' % c, 'r')
    ys = h5py.File(outdir + 'CLUSTER_%04d_member_ys.hdf5' % c, 'r')
    zs = h5py.File(outdir + 'CLUSTER_%04d_member_zs.hdf5' % c, 'r')
    vx = h5py.File(outdir + 'CLUSTER_%04d_member_vx.hdf5' % c, 'r')
    vy = h5py.File(outdir + 'CLUSTER_%04d_member_vy.hdf5' % c, 'r')
    vz = h5py.File(outdir + 'CLUSTER_%04d_member_vz.hdf5' % c, 'r')
    r200s = h5py.File(outdir + 'CLUSTER_%04d_member_r200s.hdf5' % c, 'r')
    ms = h5py.File(outdir + 'CLUSTER_%04d_member_ms.hdf5' % c, 'r')

    hf_x = h5py.File(outdir+'position_wrt_group/CLUSTER_%04d_member_xs_group'
            '.hdf5' % c, 'w')
    hf_y = h5py.File(outdir+'position_wrt_group/CLUSTER_%04d_member_ys_group'
        '.hdf5' % c, 'w')
    hf_z = h5py.File(outdir+'position_wrt_group/CLUSTER_%04d_member_zs_group'
            '.hdf5' % c, 'w')
    hf_vx = h5py.File(outdir+'position_wrt_group/CLUSTER_%04d_member_vx_group'
            '.hdf5' % c, 'w')
    hf_vy = h5py.File(outdir+'position_wrt_group/CLUSTER_%04d_member_vy_group'
            '.hdf5' % c, 'w')
    hf_vz = h5py.File(outdir+'position_wrt_group/CLUSTER_%04d_member_vz_group'
            '.hdf5' % c, 'w')
    hf_r2 = h5py.File(outdir+'position_wrt_group/CLUSTER_%04d_member_r200host'
            '_group.hdf5' % c, 'w')
    hf_m = h5py.File(outdir+'position_wrt_group/CLUSTER_%04d_member_m200host'
            '_group.hdf5' % c, 'w')

    for i in range(len(inf_grou)):
        ids_m = np.array(ids[inf_grou[i]][inf_memb[i]])//mod
        ids_h = np.array(ids[inf_grou[i]]['              0'])//mod
        
        xyzs_memb_ld = np.transpose(np.array([xs[inf_grou[i]][inf_memb[i]], 
                ys[inf_grou[i]][inf_memb[i]], zs[inf_grou[i]][inf_memb[i]],
                vx[inf_grou[i]][inf_memb[i]], vy[inf_grou[i]][inf_memb[i]],
                vz[inf_grou[i]][inf_memb[i]]]))
        xyzs_host_ld = np.transpose(np.array([
                xs[inf_grou[i]]['              0'], 
                ys[inf_grou[i]]['              0'], 
                zs[inf_grou[i]]['              0'],
                vx[inf_grou[i]]['              0'], 
                vy[inf_grou[i]]['              0'], 
                vz[inf_grou[i]]['              0'], 
                r200s[inf_grou[i]]['              0'],
                ms[inf_grou[i]]['              0']]))
        
        xyzs_memb = np.zeros((1+ids_h[-1]-ids_h[0], 6))
        xyzs_host = np.zeros((1+ids_h[-1]-ids_h[0], 8))
        
        ids_m, ids_h = ids_m-ids_m[0], ids_h-ids_m[0]
        ids_m_crop = ids_m < (ids_h[-1]+1)
        xyzs_memb_ld = xyzs_memb_ld[ids_m_crop]

        ids_m = ids_m[ids_m_crop]
        ids_h_bool = np.zeros((len(xyzs_memb), 6)) > 1.
        ids_h_bool[ids_h] = True
        ids_m_bool = np.zeros((len(xyzs_memb), 8)) > 1.
        ids_m_bool[ids_m] = True

        xyzs_memb[ids_m] = xyzs_memb_ld
        xyzs_host[ids_h] = xyzs_host_ld
        
        xyzs_memb *= ids_h_bool
        xyzs_host *= ids_m_bool
        
        #position relative to host group centre
        xyzs_memb_h = (xyzs_memb[:, :3] - xyzs_host[:, :3])
        vels_memb = (xyzs_memb[:, 3:6] - xyzs_host[:, 3:6])
        
        hf_x.create_dataset(inf_grou[i]+'/'+inf_memb[i],data=xyzs_memb_h[:,0])
        hf_y.create_dataset(inf_grou[i]+'/'+inf_memb[i],data=xyzs_memb_h[:,1])
        hf_z.create_dataset(inf_grou[i]+'/'+inf_memb[i],data=xyzs_memb_h[:,2])
        hf_vx.create_dataset(inf_grou[i]+'/'+inf_memb[i],data=vels_memb[:,0])
        hf_vy.create_dataset(inf_grou[i]+'/'+inf_memb[i],data=vels_memb[:,1])
        hf_vz.create_dataset(inf_grou[i]+'/'+inf_memb[i],data=vels_memb[:,2])
        hf_r2.create_dataset(inf_grou[i]+'/'+inf_memb[i],data=xyzs_host[:,6])
        hf_m.create_dataset(inf_grou[i]+'/'+inf_memb[i],data=xyzs_host[:,7])

        #print(xyzs_memb_h[:, 0])
        #print(xyzs_memb_h[:, 1])
        #print(xyzs_memb_h[:, 2])

    hf_x.close()
    hf_y.close()
    hf_z.close()
    hf_vx.close()
    hf_vy.close()
    hf_vz.close()
    hf_r2.close()
    hf_m.close()


    return None

###############################UNFINISHED#########################
def find_pos_relative_cluster(c, outdir, datadir):
    """ Produce positions of all objects in groups, wrt cluster 
    centre """

    if not os.path.exists(outdir+'position_wrt_cluster'):
        os.mkdir(outdir+'position_wrt_cluster')

    tracking_inf = ld_arr(track_dir+'forward_tracking.txt')
    inf_clus = np.array(tracking_inf[:, 4], dtype='int')==c
    inf_grou = np.char.mod('%15d', tracking_inf[:, 3])[inf_clus]
    inf_memb = np.char.mod('%15d', tracking_inf[:, 2])[inf_clus]

    ids = h5py.File(outdir + 'CLUSTER_%04d_member_branches.hdf5' % c, 'r')
    xs = h5py.File(outdir + 'CLUSTER_%04d_member_xs.hdf5' % c, 'r')
    ys = h5py.File(outdir + 'CLUSTER_%04d_member_ys.hdf5' % c, 'r')
    zs = h5py.File(outdir + 'CLUSTER_%04d_member_zs.hdf5' % c, 'r')
    vx = h5py.File(outdir + 'CLUSTER_%04d_member_vx.hdf5' % c, 'r')
    vy = h5py.File(outdir + 'CLUSTER_%04d_member_vy.hdf5' % c, 'r')
    vz = h5py.File(outdir + 'CLUSTER_%04d_member_vz.hdf5' % c, 'r')

    hc_xs = h5py.File(datadir + 'xs/CLUSTER_%04d_xs' % c, 'r')
    hc_ys = h5py.File(datadir + 'ys/CLUSTER_%04d_ys' % c, 'r')
    hc_zs = h5py.File(datadir + 'zs/CLUSTER_%04d_zs' % c, 'r')
    hc_vx = h5py.File(datadir + 'vx/CLUSTER_%04d_vx' % c, 'r')
    hc_vy = h5py.File(datadir + 'vy/CLUSTER_%04d_vy' % c, 'r')
    hc_vz = h5py.File(datadir + 'vz/CLUSTER_%04d_vz' % c, 'r')
    hc_r2 = h5py.File(datadir + 'rvirs/CLUSTER_%04d_rvirs' % c, 'r')
    c_ids = np.append(np.array(np.zeros(20), dtype='int'), 
            host_ids[c-1])
    c_sna = c_ids//mod
    c_ids -= (mod*c_sna)
    c_xs, c_ys, c_zs = np.zeros(129), np.zeros(129), np.zeros(129)
    c_vx, c_vy, c_vz = np.zeros(129), np.zeros(129), np.zeros(129)
    c_r2 = np.zeros(129)
    for i in range(129):
        if c_ids[i] > 0:
            c_xs[i] = hc_xs['%03d'%i][c_ids[i]-1]
            c_ys[i] = hc_ys['%03d'%i][c_ids[i]-1]
            c_zs[i] = hc_zs['%03d'%i][c_ids[i]-1]
            c_vx[i] = hc_vx['%03d'%i][c_ids[i]-1]
            c_vy[i] = hc_vy['%03d'%i][c_ids[i]-1]
            c_vz[i] = hc_vz['%03d'%i][c_ids[i]-1]
            c_r2[i] = hc_r2['%03d'%i][c_ids[i]-1]

    hf_x = h5py.File(outdir+'position_wrt_cluster/CLUSTER_%04d_member_xs_'
            'cluster.hdf5' % c, 'w')
    hf_y = h5py.File(outdir+'position_wrt_cluster/CLUSTER_%04d_member_ys_'
            'cluster.hdf5' % c, 'w')
    hf_z = h5py.File(outdir+'position_wrt_cluster/CLUSTER_%04d_member_zs_'
            'cluster.hdf5' % c, 'w')
    hf_vx = h5py.File(outdir+'position_wrt_cluster/CLUSTER_%04d_member_vx_'
            'cluster.hdf5' % c, 'w')
    hf_vy = h5py.File(outdir+'position_wrt_cluster/CLUSTER_%04d_member_vy_'
            'cluster.hdf5' % c, 'w')
    hf_vz = h5py.File(outdir+'position_wrt_cluster/CLUSTER_%04d_member_vz_'
            'cluster.hdf5' % c, 'w')
    hf_r2 = h5py.File(outdir+'position_wrt_cluster/CLUSTER_%04d_member_'
            'r200host_cluster.hdf5' % c, 'w')

    norepeat = False
    for i in range(len(inf_grou)):
        if inf_grou[i] != inf_grou[i-1] or (max(inf_grou)==min(inf_grou)
                and i==0):
            obj_id = '              0'

            ids_m = np.array(ids[inf_grou[i]][obj_id])//mod
            ids_h = c_sna[ids_m[0]:ids_m[-1]+1]
            ids_h = ids_h[ids_h>0]
        
            xyzs_memb_ld = np.transpose(np.array([xs[inf_grou[i]][obj_id], 
                    ys[inf_grou[i]][obj_id], zs[inf_grou[i]][obj_id],
                    vx[inf_grou[i]][obj_id], vy[inf_grou[i]][obj_id],
                    vz[inf_grou[i]][obj_id]]))
            xyzs_host_ld = np.transpose(np.array([c_xs[ids_h], c_ys[ids_h],
                    c_zs[ids_h], c_vx[ids_h], c_vy[ids_h], c_vz[ids_h], 
                    c_r2[ids_h]]))
                    
                    
            xyzs_memb = np.zeros((1+ids_m[-1]-ids_m[0], 6))
            xyzs_host = np.zeros((1+ids_m[-1]-ids_m[0], 7))
            
            ids_m, ids_h = ids_m-ids_m[0], ids_h-ids_m[0]
            ids_m_crop = ids_m < (ids_h[-1]+1)
            xyzs_memb_ld = xyzs_memb_ld[ids_m_crop]

            ids_m = ids_m[ids_m_crop]
            ids_h_bool = np.zeros((len(xyzs_memb), 6)) > 1.
            ids_h_bool[ids_h] = True
            ids_m_bool = np.zeros((len(xyzs_memb), 7)) > 1.
            ids_m_bool[ids_m] = True

            xyzs_memb[ids_m] = xyzs_memb_ld
            xyzs_host[ids_h] = xyzs_host_ld
            
            xyzs_memb *= ids_h_bool
            xyzs_host *= ids_m_bool
            
            #position relative to cluster centre
            xyzs_memb_h = (xyzs_memb[:, :3] - xyzs_host[:, :3])
            vels_memb = (xyzs_memb[:, 3:6] - xyzs_host[:, 3:6])
            
            hf_x.create_dataset(inf_grou[i]+'/'+obj_id,data=xyzs_memb_h[:,0])
            hf_y.create_dataset(inf_grou[i]+'/'+obj_id,data=xyzs_memb_h[:,1])
            hf_z.create_dataset(inf_grou[i]+'/'+obj_id,data=xyzs_memb_h[:,2])
            hf_vx.create_dataset(inf_grou[i]+'/'+obj_id,data=vels_memb[:,0])
            hf_vy.create_dataset(inf_grou[i]+'/'+obj_id,data=vels_memb[:,1])
            hf_vz.create_dataset(inf_grou[i]+'/'+obj_id,data=vels_memb[:,2])
            hf_r2.create_dataset(inf_grou[i]+'/'+obj_id,data=xyzs_host[:,6])

            
        obj_id = inf_memb[i]

        ids_m = np.array(ids[inf_grou[i]][obj_id])//mod
        ids_h = c_sna[ids_m[0]:ids_m[-1]+1]
        ids_h = ids_h[ids_h>0]
        
        xyzs_memb_ld = np.transpose(np.array([xs[inf_grou[i]][obj_id], 
                ys[inf_grou[i]][obj_id], zs[inf_grou[i]][obj_id],
                vx[inf_grou[i]][obj_id], vy[inf_grou[i]][obj_id],
                vz[inf_grou[i]][obj_id]]))
        xyzs_host_ld = np.transpose(np.array([c_xs[ids_h], c_ys[ids_h], 
                c_zs[ids_h], c_vx[ids_h], c_vy[ids_h], c_vz[ids_h], 
                c_r2[ids_h]]))
                
                
        xyzs_memb = np.zeros((1+ids_m[-1]-ids_m[0], 6))
        xyzs_host = np.zeros((1+ids_m[-1]-ids_m[0], 7))
        
        ids_m, ids_h = ids_m-ids_m[0], ids_h-ids_m[0]
        ids_m_crop = ids_m < (ids_h[-1]+1)
        xyzs_memb_ld = xyzs_memb_ld[ids_m_crop]

        ids_m = ids_m[ids_m_crop]
        ids_h_bool = np.zeros((len(xyzs_memb), 6)) > 1.
        ids_h_bool[ids_h] = True
        ids_m_bool = np.zeros((len(xyzs_memb), 7)) > 1.
        ids_m_bool[ids_m] = True

        xyzs_memb[ids_m] = xyzs_memb_ld
        xyzs_host[ids_h] = xyzs_host_ld
        
        xyzs_memb *= ids_h_bool
        xyzs_host *= ids_m_bool
        
        #position relative to cluster centre
        xyzs_memb_h = (xyzs_memb[:, :3] - xyzs_host[:, :3])
        vels_memb = (xyzs_memb[:, 3:6] - xyzs_host[:, 3:6])
        
        hf_x.create_dataset(inf_grou[i]+'/'+obj_id,data=xyzs_memb_h[:,0])
        hf_y.create_dataset(inf_grou[i]+'/'+obj_id,data=xyzs_memb_h[:,1])
        hf_z.create_dataset(inf_grou[i]+'/'+obj_id,data=xyzs_memb_h[:,2])
        hf_vx.create_dataset(inf_grou[i]+'/'+obj_id,data=vels_memb[:,0])
        hf_vy.create_dataset(inf_grou[i]+'/'+obj_id,data=vels_memb[:,1])
        hf_vz.create_dataset(inf_grou[i]+'/'+obj_id,data=vels_memb[:,2])
        hf_r2.create_dataset(inf_grou[i]+'/'+obj_id,data=xyzs_host[:,6])


    hf_x.close()
    hf_y.close()
    hf_z.close()
    hf_vx.close()
    hf_vy.close()
    hf_vz.close()
    hf_r2.close()


    return None




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


def bound(vrels, rrels, ms, r200s):
    """ Determine whether a halo is bound to its host """
    ke = 0.5 * (1000.*vrels)**2.
    ge = ((1.3271244*10.**20.) * ms) / (3.0857*10.**19.)
    vir = (ke - (ge / rrels)) < (ge / (-2.5*r200s))
    #v_cr = 0.001*(1.2*ge/r200s)**0.5

    return vir#, rrels/r200s, vrels/v_cr




if run_find_member_trees==True:
    for c_id in crange:
        print(c_id)
        dir_name = out_dir+'Group_member_tracking/CLUSTER_%04d/' % c_id
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        find_member_trees(c_id, dir_name, total)

if run_find_member_data==True:
    for c_id in crange:
        print(c_id)
        find_member_data(c_id, out_dir+'Group_member_tracking/CLUSTER_'
                +'%04d/' % c_id, data_dir)


if run_find_boundness==True:
    for c_id in crange:
        print(c_id)
        find_boundness(c_id, out_dir+'Group_member_tracking/CLUSTER_'
                +'%04d/' % c_id)

if run_tracking_wrt_group==True:
    for c_id in crange:
        print(c_id)
        find_pos_relative_group(c_id, out_dir+'Group_member_tracking/'
                'CLUSTER_%04d/' % c_id)

if run_tracking_wrt_cluster==True:
    for c_id in crange:
        print(c_id)
        find_pos_relative_cluster(c_id, out_dir+'Group_member_tracking/'
                'CLUSTER_%04d/' % c_id, data_dir)
