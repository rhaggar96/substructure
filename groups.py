from modules import *


def bound(vrels, rrel, m, r200):
    rrels = rrel
    rrels[rrels==0.] = 0.001
    ke = 0.5 * (1000.*vrels)**2.
    ge = ((1.3271244*10.**20.) * m) / (3.0857*10.**19.)
    vir = (ke - (ge / rrels)) < (ge / (-2.5*r200))
    v_cr = 0.001*(1.2*ge/r200)**0.5
    return vir, rrels/r200, vrels/v_cr

def find_group_members(c):
    c_id = clus_ids[c-1]
    loaddir = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/MergerTree'
            'AHF_General_Tree_Comp/NewMDCLUSTER_'
            '%04d/snap_128/CLUSTER_%04d_' % (c, c))
    z_list = redshifts[c-1]

    npar_i = ld_arr(loaddir+'npars.txt', dtype='int')[:,-1]
    npar_i = np.where(npar_i < p_lim)[0][0]
    xs_i = ld_arr(loaddir+'xs.txt')[:npar_i]
    ys_i = ld_arr(loaddir+'ys.txt')[:npar_i]
    zs_i = ld_arr(loaddir+'zs.txt')[:npar_i]
    vx_i = ld_arr(loaddir+'vx.txt')[:npar_i]
    vy_i = ld_arr(loaddir+'vy.txt')[:npar_i]
    vz_i = ld_arr(loaddir+'vz.txt')[:npar_i]
    ms_i = ld_arr(loaddir+'ms.txt')[:npar_i]
    starlim = ld_arr(loaddir+'mstars.txt', dtype='int')[:npar_i] <= mstarlim
    rvirs_i = ld_arr(loaddir+'rvirs.txt')[:npar_i]
    subs_i = ld_arr(loaddir+'subs.txt', dtype='int')[:npar_i]
    
    xs_i[starlim], ys_i[starlim], zs_i[starlim] = 0., 0., 0.
    vx_i[starlim], vy_i[starlim], vz_i[starlim] = 0., 0., 0.
    ms_i[starlim], rvirs_i[starlim], subs_i[starlim] = 0., 0., 0.

    bool_t = xs_i[:, -1]>0
    xs_i, ys_i, zs_i = xs_i[bool_t], ys_i[bool_t], zs_i[bool_t]
    vx_i, vy_i, vz_i = vx_i[bool_t], vy_i[bool_t], vz_i[bool_t]
    ms_i, rvirs_i, subs_i = ms_i[bool_t], rvirs_i[bool_t], subs_i[bool_t]
    
    diff = ((xs_i-xs_i[c_id])**2. + (ys_i-ys_i[c_id])**2.
            + (zs_i-zs_i[c_id])**2.)**0.5

    offset = len(z_list) - len(diff[0])

    rs_tot = np.zeros(0)
    vs_tot = np.zeros(0)
    #bool_tot = np.zeros(0)
    #subs_tot = np.zeros(0)
    n_memb = np.array(np.zeros(0), dtype='int')
    r_virs = np.zeros(0)
    z_infa = np.zeros(0)
    for h_id in range(len(xs_i)):
        snap = np.where(diff[h_id]<rvirs_i[c_id])[0]
        if len(snap) > 0:
            snap = snap[0]
        
            xs = (xs_i[:,snap] - xs_i[h_id,snap])
            ys = (ys_i[:,snap] - ys_i[h_id,snap])
            zs = (zs_i[:,snap] - zs_i[h_id,snap])
            rs = (xs**2. + ys**2. + zs**2.)**0.5

            vx = (vx_i[:,snap] - vx_i[h_id,snap])
            vy = (vy_i[:,snap] - vy_i[h_id,snap])
            vz = (vz_i[:,snap] - vz_i[h_id,snap])
            vs = (vx**2. + vy**2. + vz**2.)**0.5

            bs = list(bound(vs, rs, ms_i[h_id,snap], rvirs_i[h_id,snap]))
            #bs[0][h_id]=False
            bs[0][:h_id+1]=False
            bs[0] = bs[0] * (rvirs_i[h_id,snap]>=6.5)
            rs_tot = np.append(rs_tot, bs[1][bs[0]])
            vs_tot = np.append(vs_tot, bs[2][bs[0]])
            #bool_tot = np.append(bool_tot, bs[0])
            #subs_tot = np.append(subs_tot, subs_i[:,snap][bs[0]]==128*mod+1+h_id)
            n_memb = np.append(n_memb, int(np.sum(bs[0])))
            r_virs = np.append(r_virs, rvirs_i[h_id,snap])
            z_infa = np.append(z_infa, z_list[snap+offset])
    
    print 'group members: ' + str(np.sum(n_memb))
    return rs_tot, vs_tot, n_memb, r_virs, z_infa



clus_ids = ld_arr('G3X_data/ds_infor_G3X_progen/DS_G3X_snap_128_center-'
        'cluster_progenitors.txt', dtype='int')[:, 1]-(128*mod+1)  

crange = np.array(np.loadtxt('G3X_data/G3X_300_selected_sample_257.txt'), 
        dtype='int')

redshifts = ld_arr('G3X_data/G3X_300_redshifts.txt')

p_lim = 100
mstarlim = 10.**9.5