from modules import *


relax = ld_arr(('/home/ppxrh2/Documents/test_pollux/TheThreeHundred/playground'
        + '/rhaggar/G3X_data/ds_infor_G3X_progen/'
        + 'DS_G3X_snap_128_center-cluster_progenitors.txt'))[:, [3, 4, 5]]

relax_tot = np.loadtxt(('/home/ppxrh2/Documents/test_pollux/TheThreeHundred/'
        + 'playground/rhaggar/G3X_data/'
        + 'G3X_snap_128_relaxation_param.txt'))
relax_tot = np.transpose(relax_tot)

fs_new = ld_arr(('/home/ppxrh2/Documents/test_pollux/TheThreeHundred/'
        + 'playground/rhaggar/G3X_data/ds_infor_G3X_progen/'
        + 'DS_G3X_snap_128_fs-lim_center-cluster_.txt'))[:, 2]/0.1

relax_new = (3. / (
        (abs(1.-relax[:,0])/0.15)**2. + 
        (relax[:,1]/0.04)**2. +
        (fs_new/0.1)**2. +
        #(relax[:,2]/0.1)**2. +
        0))**0.5


masses = np.loadtxt(('/home/ppxrh2/Documents/test_pollux/TheThreeHundred/'
        + 'playground/rhaggar/Backsplash/G3X_data/'
        + 'G3X_300_cluster_masses.txt'))
masses = np.transpose(np.log10(masses))


def rolling_av(xs, ys, bs=100, spread=20):
    bins = np.transpose(np.linspace(min(xs), max(xs), bs))
    b_wid = (max(xs)-min(xs))*spread/bs
    res = np.zeros((bs, 3))
    for i in range(bs):
        res[i] = bootstrap(ys[(xs>bins[i]-b_wid)*(xs<bins[i]+b_wid)], 100)
        res[i, 0] = np.average(ys[(xs>bins[i]-b_wid)*(xs<bins[i]+b_wid)])
    return bins, res


meds_n = rolling_av(masses, fs_new)

#relax_par = relax_tot
#relax_par = abs(relax[:, 0]-1.)/0.15
#relax_par = relax[:, 1]/0.04
relax_par = relax[:, 2]/0.1


meds = rolling_av(masses, relax_par)

plt.figure(figsize=(6, 5))
#plt.scatter(masses, relax_par, c='b', s=5.)
plt.scatter(masses, fs_new, c='r', s=5.)
boo = (meds[0]>14.8)*(meds[0]<15.2)
plt.plot(meds[0][boo], meds[1][:, 0][boo], c='k', linewidth=2.)
plt.fill_between(meds[0][boo], (meds[1][:, 0][boo]+meds[1][:, 1][boo]), (meds[1][:, 0][boo]-meds[1][:, 2][boo]), color='b', alpha=0.5, linewidth=0.)

boo = (meds_n[0]>14.8)*(meds_n[0]<15.2)
plt.plot(meds_n[0][boo], meds_n[1][:, 0][boo], c='k', linewidth=2.)
plt.fill_between(meds_n[0][boo], (meds_n[1][:, 0][boo]+meds_n[1][:, 1][boo]), (meds_n[1][:, 0][boo]-meds_n[1][:, 2][boo]), color='r', alpha=0.5, linewidth=0.)

#plt.plot(meds_n[0][boo], (meds_n[1][:, 0]/meds[1][:, 0])[boo], c='k', linewidth=3.)

plt.plot([14., 16.], [1., 1.], c='k', linestyle='--')
plt.xlim(14.6, 15.5)
plt.ylim(-0.1, 3.)
plt.xlabel(r'$\textrm{log}_{10}M_{200}$')
plt.ylabel(r'$f_s$')
plt.tight_layout()
#plt.savefig('fs_mean_new.png', dpi=200)
plt.show()