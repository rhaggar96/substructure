from modules import *

run_full = False
remove_small_hs = True
save = False
group_sizes = True

p_lim = 100
mstarlim = 10.**9.5

clus_ids = ld_arr(G3X_data+'ds_infor_G3X_progen/DS_G3X_snap_128_center-'
        'cluster_progenitors.txt', dtype='int')[:, 1]-(128*mod+1)  

crange = np.array(np.loadtxt(G3X_data+'G3X_300_selected_sample_257.txt'), 
        dtype='int')

redshifts = ld_arr(G3X_data+'G3X_300_redshifts.txt')

def bound(vrels, rrel, m, r200):
    """ Determine whether a halo is bound to its host """
    rrels = rrel
    rrels[rrels==0.] = 0.001
    ke = 0.5 * (1000.*vrels)**2.
    ge = ((1.3271244*10.**20.) * m) / (3.0857*10.**19.)
    vir = (ke - (ge / rrels)) < (ge / (-2.5*r200))
    v_cr = 0.001*(1.2*ge/r200)**0.5
    return vir, rrels/r200, vrels/v_cr

def bound_crit(x):
    """ y = v/vcrit, x = r/R200, find boundness condition """
    return ((5./(3.*x)) - (2./3.))**0.5

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
            if remove_small_hs == True:
                bs[0][:h_id+1]=False
                bs[0] = bs[0] * (rvirs_i[h_id,snap]>=6.5)
            else:
                bs[0][h_id]=False
            rs_tot = np.append(rs_tot, bs[1][bs[0]])
            vs_tot = np.append(vs_tot, bs[2][bs[0]])
            n_memb = np.append(n_memb, int(np.sum(bs[0])))
            r_virs = np.append(r_virs, rvirs_i[h_id,snap])
            z_infa = np.append(z_infa, z_list[snap+offset])
    
    print 'group members: ' + str(np.sum(n_memb))
    return rs_tot, vs_tot, n_memb, r_virs, z_infa

def find_group_num_dist(nm, lims):
    nm = nm[nm>0]
    total = float(len(nm))
    result = np.zeros(0)
    result = np.append(result, float(np.sum((nm <= lims[0]))) / total)
    for i in range(len(lims)-1):
        result = np.append(result, float(np.sum((nm > lims[i])*(nm <= 
                lims[i+1]))) / total)
    result = np.append(result, float(np.sum((nm > lims[-1]))) / total)
    return result



if run_full == True:
    rs_total = np.zeros(0)
    vs_total = np.zeros(0)
    nm_total = np.zeros(0)
    r2_total = np.zeros(0)
    zi_total = np.zeros(0)
    member_no = np.zeros(0)
    total_no = np.zeros(0)

    for c_val in crange:
        print c_val
        result = find_group_members(c_val)
        rs_total = np.append(rs_total, result[0])
        vs_total = np.append(vs_total, result[1])
        nm_total = np.append(nm_total, result[2])
        r2_total = np.append(r2_total, result[3])
        zi_total = np.append(zi_total, result[4])
        member_no = np.append(member_no, np.sum(result[2]))
        total_no = np.append(total_no, len(result[2]))

    if not os.path.exists('data_out'):
        os.mkdir('data_out')
    pd.DataFrame(np.char.mod('%12.9f', np.transpose(np.array(
            [rs_total, vs_total])))).to_csv('data_out/rs_vs_257.txt', 
            sep='\t', index=None, header=[' r/R200_host', '    v/v_crit'])
    pd.DataFrame(np.char.mod('%6d', nm_total)).to_csv('data_out/n_memb_257'
            '.txt', index=None, header=['n_memb'])
    pd.DataFrame(np.char.mod('%12.7f', r2_total)).to_csv('data_out/r_200s_257'
            '.txt', index=None, header=['       r_200'])
    pd.DataFrame(np.char.mod('%6.3f', zi_total)).to_csv('data_out/z_infall_'
            '257.txt', index=None, header=['z_infa'])
    pd.DataFrame(np.char.mod('%6d', np.transpose(np.array(
            [crange, member_no, total_no])))).to_csv('data_out/groups_257.'
            'txt', sep='\t', index=None, header=['  c_id','groups','infall'])

else:
    rs_total, vs_total = np.transpose(ld_arr('data_out/rs_vs_257.txt'))
    zi_total = ld_arr('data_out/z_infall_257.txt')
    nm_total = ld_arr('data_out/n_memb_257.txt')
    r2_total = ld_arr('data_out/r_200s_257.txt')


plt.figure()
plt.scatter(rs_total*1., vs_total, c='b', s=5.)
plt.plot(np.arange(1, 251)/100., bound_crit(np.arange(1, 251)/100.), c='r', 
        linewidth=2.)
plt.xlim(0., 3.)
plt.ylim(0., 2.5)
plt.xlabel(r'$r/R_{\rm{200,group}}$')
plt.ylabel(r'$v/v_{\rm{crit}}$')
plt.tight_layout()
if save==True:
    plt.savefig('data_out/fig1_masslim_rlim.pdf')
#plt.show()


if group_sizes==True:
    grp_bins = [1, 17]
    print 'Upper limits on group size bins: ' + str(grp_bins)
    print find_group_num_dist(np.array(np.transpose(nm_total), 
            dtype='int')[0], grp_bins)
