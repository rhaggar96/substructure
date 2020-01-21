from modules import *



run_infall_finder = True
run_infaller_branches = True
run_bounded_at_infall = True
run_all_grp_member_branches = True


all_infalling_objects = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
        + 'Infalling_Groups/MergerTreeAHF_Infalling_Re-written/all_infalling'
        + '_objects/')
all_infalling_branches = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
        + 'Infalling_Groups/MergerTreeAHF_Infalling_Re-written/all_infaller'
        + '_branches/')
all_infalling_bounded = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
        + 'Infalling_Groups/MergerTreeAHF_Infalling_Re-written/all_bounded_at'
        + '_infall_ids/')
all_grp_member_branches = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
        + 'Infalling_Groups/MergerTreeAHF_Infalling_Re-written/all_group_'
        + 'member_branches/')
halo_data = ('/run/media/ppxrh2/166AA4B87A2DD3B7/NewMDCLUSTER_data/reduced_'
        + 'cluster_info/')



def find_infalling(c, outdir):
    """ Find all the objects which pass within R200 of the cluster, 
    having been outside of R200 at a previous snapshot. Does not 
    account for double infalls, sub-subhaloes, etc. """

    loaddir = ('/run/media/ppxrh2/166AA4B87A2DD3B7/NewMDCLUSTER_data/reduced_'
            'cluster_info/')
    xs_in  = h5py.File(loaddir+'xs/CLUSTER_%04d_xs' % c)
    ys_in  = h5py.File(loaddir+'ys/CLUSTER_%04d_ys' % c)
    zs_in  = h5py.File(loaddir+'zs/CLUSTER_%04d_zs' % c)
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
    ids = np.array(np.arange(len(xs)), dtype='int')+int(j)*mod+1
    r200s = np.array(r200_in[j])
    
    h_id = halo_ids[0] - (int(j)*mod+1)
    clust = [xs[h_id], ys[h_id], zs[h_id], r200s[h_id]]
        
    #ids_out: id of every object outside the cluster in the first snapshot
    xs, ys, zs = xs-clust[0], ys-clust[1], zs-clust[2]
    rs = (xs**2. + ys**2. + zs**2.)**0.5 / clust[3]
    ids_out = ids[rs > 1.]


    ids_infall_t = np.array(np.zeros(0), dtype='int')
    
    for i in range(1, len(keys)):
        j = keys[i]

        xs = np.array(xs_in[j])
        ys = np.array(ys_in[j])
        zs = np.array(zs_in[j])
        r200s = np.array(r200_in[j])
        ids = np.array(np.arange(len(xs)), dtype='int')+int(j)*mod+1
        ids_main_p_in = ld_arr('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTree'
                'AHF/MergerTreeAHF_HDF5_MainProg/CLUSTER_%04d_CRMratio2_Main'
                'Prog/CLUSTER_%04d_prog_output.snapshot_' % (c, c) + j
                + '.AHF.MainProg.txt', dtype='int')
        j = int(j)
        ids_main_p = np.array(np.zeros(len(xs)), dtype='int')
        ids_main_p[ids_main_p_in[:, 0]-(j*mod+1)] = ids_main_p_in[:, 1]
        
        h_id = halo_ids[i] - (int(j)*mod+1)
        clust = [xs[h_id], ys[h_id], zs[h_id], r200s[h_id]]

        xs, ys, zs = xs-clust[0], ys-clust[1], zs-clust[2]
        rs = (xs**2. + ys**2. + zs**2.)**0.5 / clust[3]


        #ids_infall: ids of haloes that are now in r200
        #ids_prev: ids of where haloes in r200 previously were located 
        #ids_out: haloes in ALL previous snapshots which were outside r200 
        ids_infall = ids[rs < 1.]
        ids_prev = ids_main_p[ids_infall-(j*mod+1)]

        ids_infall = ids_infall[ids_prev > 0]
        ids_prev = ids_prev[ids_prev > 0]
        for k in range(len(ids_prev)):
            if len(np.where(ids_out==ids_prev[k])[0])==0:
                ids_prev[k] = 0


        #ids_infall: ids of every halo which has JUST fallen into the r200
        ids_infall = ids_infall[ids_prev>0]
        ids_infall_t = np.append(ids_infall_t, ids_infall)

        ids_out = np.append(ids_out, ids[rs > 1.])

    pd.DataFrame(np.char.mod('%15d', ids_infall_t)).to_csv(outdir + 
            '/CLUSTER_%04d_all_infallers.txt' % c, 
            header=['  all_infallers'], index=None)
    

    return None


def find_infaller_evolution(c, indir, outdir):
    """ Finds the subsequent evolution of every object that has fallen 
    into the cluster. Objects with multiple infalls are eliminated, 
    such that only the ID at first infall is recorded. The total number 
    of infalls is stored """

    total_tree = np.array(pd.read_csv('/run/media/ppxrh2/166AA4B87A2DD3B7/'
            'MergerTreeAHF/MergerTreeAHF_ASCII/MergerTree_GadgetX-NewMD'
            'CLUSTER_%04d.txt-CRMratio2' % c, sep='\s+', skiprows=2, 
            usecols=[0], dtype='str')[:-1], dtype='int')[:, 0]
    
    ids_infall_t = ld_arr(indir + '/CLUSTER_%04d_all_infallers.txt' % c, 
            dtype='int')[:, 0]

    hf = h5py.File(outdir + 'CLUSTER_%04d_groups.hdf5' % c, 'w')
    infall_number = np.array(np.ones(len(ids_infall_t)), dtype='int')
    i_nonrep = 0
    for i in range(len(ids_infall_t)):
        branch = find_branch(np.array([ids_infall_t[i]]),total_tree)
        fate = branch[1]
        branch = branch[0]
        for j in branch[1:]:
            repeat = np.where(ids_infall_t == j)[0]
            if len(repeat) == 1:
                infall_number[repeat] = infall_number[i] + 1
        if len(branch) > 0:
            #label dataset with index, final state amd number of infalls
            hf.create_dataset('%04d_%01d_%02d' % (i_nonrep, fate, 
                    infall_number[i]), data=branch)
            i_nonrep += 1
    
    hf.close()

    return None



def find_bound_groups(c, loaddir, datadir, outdir):
    """ Take the full list of all infalling objects (regardless of 
    mass, etc.) and find which objects are bound to them at infall. 
    Outputs lists labelled with an index, number of infalls of the 
    host, and the number of members in the group (including the host). 
    List begins with host ID, followed by IDs of any additional bound 
    objects, infalling or not. No restictions on repeated objects. """
    
    infalls = h5py.File(loaddir + 'CLUSTER_%04d_groups.hdf5' % c, 'r')
    keys = np.array(list(infalls.keys()), dtype='str')

    #ids of everything infalling
    infall_ids_t = np.array([np.array(infalls[keys[i]], dtype='int')[0] 
            for i in np.arange(len(keys))])
    
    infall_snaps = infall_ids_t // mod
    infall_ids = infall_ids_t - (infall_snaps*mod+1) #zero index
    #infall_snaps: snapshots at which objects infall
    #infall_snaps_ct: number of objects at each snapshot
    infall_snaps_ct = np.array(list(Counter(infall_snaps).values()))
    infall_snaps = np.array(list(Counter(infall_snaps).keys()))

    #loading relevant data for finding boundness
    xs = h5py.File(datadir + 'xs/CLUSTER_%04d_xs' % c, 'r')
    ys = h5py.File(datadir + 'ys/CLUSTER_%04d_ys' % c, 'r')
    zs = h5py.File(datadir + 'zs/CLUSTER_%04d_zs' % c, 'r')
    vx = h5py.File(datadir + 'vx/CLUSTER_%04d_vx' % c, 'r')
    vy = h5py.File(datadir + 'vy/CLUSTER_%04d_vy' % c, 'r')
    vz = h5py.File(datadir + 'vz/CLUSTER_%04d_vz' % c, 'r')
    ms = h5py.File(datadir + 'ms/CLUSTER_%04d_ms' % c, 'r')
    r200s = h5py.File(datadir + 'rvirs/CLUSTER_%04d_rvirs' % c, 'r')

    counter = 0
    #file to write lists of bound objects to
    bound_out = h5py.File(outdir+'CLUSTER_%04d_bound_members.hdf5' % c, 'w')
    

    for s in range(len(infall_snaps)):
        snap = '%03d' % infall_snaps[s]
        xs_s, ys_s, zs_s, r200s_s = xs[snap], ys[snap], zs[snap], r200s[snap]
        vx_s, vy_s, vz_s = vx[snap], vy[snap], vz[snap]
        ms_s = np.array(ms[snap], dtype='int') 
        ids_s = np.array(np.arange(len(xs_s)), dtype='int')
        
        for i in range(infall_snaps_ct[s]):
            id_n = infall_ids[counter]
            #relative position and velocity of each halo wrt h
            rs_si = ((xs_s-xs_s[id_n])**2. + (ys_s-ys_s[id_n])**2.
                    + (zs_s-zs_s[id_n])**2.)**0.5
            vs_si = ((vx_s-vx_s[id_n])**2. + (vy_s-vy_s[id_n])**2.
                    + (vz_s-vz_s[id_n])**2.)**0.5

            #everything bounded to main infaller
            boundres = bound(vs_si, rs_si, ms_s[id_n], r200s_s[id_n])[0]

            ids_bound = infall_snaps[s]*mod+1+ids_s[boundres]
            #adding bound objects to infaller
            ids_bound = np.append(np.array([infall_ids_t[counter]]), 
                    ids_bound)

            bound_out.create_dataset(keys[counter]+'_%03d' % len(ids_bound), 
                    data=ids_bound)
            counter += 1
    bound_out.close()
    
    return None



def find_member_trees(c, loaddir, mainbranchdir, outdir):
    """ Finds the full branches of each bound member of an infalling 
    group """

    group_objs = h5py.File(loaddir+'CLUSTER_%04d_bound_members.hdf5' % c, 'r')
    infaller_branch = h5py.File(mainbranchdir+'CLUSTER_%04d_groups.hdf5' % c,
            'r')
    
    total_tree = np.array(pd.read_csv('/run/media/ppxrh2/166AA4B87A2DD3B7/'
            'MergerTreeAHF/MergerTreeAHF_ASCII/MergerTree_GadgetX-NewMD'
            'CLUSTER_%04d.txt-CRMratio2' % c, sep='\s+', skiprows=2, 
            usecols=[0], dtype='str')[:-1], dtype='int')[:, 0]

    #keys including id, fate, infall no, no of members
    keys_init = np.array(list(group_objs.keys()))
    if len(keys_init) > 0:
        keys = np.array(np.zeros((len(keys_init), 4)), dtype='str')
        keys[:, [0,1,2]] = np.core.defchararray.partition(keys_init, '_')
        keys[:, [1,2,3]] = np.core.defchararray.partition(
                keys[:, 2], '_')
        keys[:, [2,3]] = np.core.defchararray.partition(
                keys[:, 3], '_')[:, [0, 2]]

    hf = h5py.File(outdir + 'CLUSTER_%04d_grp_memb_branches.hdf5' % c, 'w')
    for i in range(len(keys_init)):
        #first, save branch of host object
        hf.create_dataset(keys_init[i] + '/000_' + keys[i,1], 
                data=infaller_branch[keys[i,0]+'_'+keys[i,1]+'_'+keys[i,2]])
        if keys[i, 3] != '001':
            for j in range(1, int(keys[i, 3])):
                #for each member, find branch and save
                grp_memb_id = np.array(list(group_objs[keys_init[i]]))[j]
                memb_branch = find_branch(np.array([grp_memb_id]), 
                        total_tree, True)
                memb_fate = memb_branch[1]
                memb_branch = memb_branch[0]
                #save in sub-directory, with id number for each member
                hf.create_dataset(keys_init[i] + '/%03d_%01d' % (j, memb_fate),
                        data=memb_branch)

    hf.close()

    return None









def find_branch(id_a, id_list, keep_sing=False):
    """ Finds the full evolution of a halo, starting with the earliest 
    id """
    id = id_a[-1]
    tree_search = np.where(id_list==id)[0]
    final_state = 0 #survives to z=0
    if len(tree_search)==0:
        result = np.array([id]*keep_sing, dtype='int')
        final_state = 1 #no progenitor - should not be the case for these
        #only the case when id_a is '0', for second infall
    elif id//mod == id_list[tree_search[0]-1]//mod:
        result = id_a
        final_state = 2 #absorbed by a larger halo - not main progenitor
    elif len(tree_search)==1 and id//mod > id_list[tree_search[0]+1]//mod:
        result = id_a
        final_state = 1 #lost by halo finder
    else:
        id_a = np.append(id_a, id_list[tree_search[0]-1])
        result, final_state = find_branch(id_a, id_list)
    
    return result, final_state


def bound(vrels, rrel, m, r200, incgroup=False):
    """ Determine whether a list of haloes are bound to their host """
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



crange = range(1, 325)

if run_infall_finder==True:
    print('run_infall_finder')
    for clus_no in crange:
        print(clus_no)
        find_infalling(clus_no, all_infalling_objects)


if run_infaller_branches==True:
    print('run_infaller_branches')
    for clus_no in crange:
        print(clus_no)
        find_infaller_evolution(clus_no, all_infalling_objects, 
                all_infalling_branches)

if run_bounded_at_infall==True:
    print('run_bounded_at_infall')
    for clus_no in crange:
        print(clus_no)
        find_bound_groups(clus_no, all_infalling_branches, halo_data, 
                all_infalling_bounded)        

if run_all_grp_member_branches==True:
    print('run_all_grp_member_branches')
    for clus_no in crange:
        print(clus_no)
        find_member_trees(clus_no, all_infalling_bounded,
                all_infalling_branches, all_grp_member_branches)

