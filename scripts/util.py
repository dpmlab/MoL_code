import numpy as np
from scipy.stats import norm
import deepdish as dd

def findcols(X):
    # X is a dictionary with elements ir and jc
    # returns a list of arrays of locations of nonzero entries in each column
    cols = [None] * (len(X['jc'])-1)
    for i in range(len(cols)):
        cols[i] = X['ir'][X['jc'][i]:X['jc'][i+1]]
    return cols

def SLtoVox(D, SLlist, nv, zeronan=True):
    # D is dict of L, R, with N x arbitrary dims
    # SLlist is dict of L, R list of length N, with vertices for each SL

    Dvox = dict()
    Dcount = dict()
    for hem in ['L', 'R']:
        Dvox[hem] = np.zeros((nv,)+ D[hem].shape[1:])
        Dcount[hem] = np.zeros((nv,)+(1,)*len(D[hem].shape[1:]))
        for i in range(len(SLlist[hem])):
            Dvox[hem][SLlist[hem][i]] += D[hem][i]
            Dcount[hem][SLlist[hem][i]] += 1

        Dcount[hem][Dcount[hem] == 0] = np.nan
        Dvox[hem] = Dvox[hem] / Dcount[hem]

        if zeronan:
            Dvox[hem][np.isnan(Dvox[hem])] = 0

    return Dvox


def nullZ(X):
    # Last dimension of X is nPerm+1, with real data at 0 element
    X_roll = np.rollaxis(X, len(X.shape)-1)
    means = X_roll[1:].mean(axis=0)
    std = X_roll[1:].std(axis=0)
    if len(X.shape) > 1:
        std[std==0] = np.nan
    Z = (X_roll[0] - means) / std
    return Z

def FDR_p(pvals, adjustment = True):
    # Port of AFNI mri_fdrize.c
    assert np.all(pvals>=0) and np.all(pvals<=1)
    pvals[pvals < np.finfo(np.float_).eps] = np.finfo(np.float_).eps
    pvals[pvals == 1] = 1-np.finfo(np.float_).eps
    n = pvals.shape[0]

    qvals = np.zeros((n))
    sorted_ind = np.argsort(pvals)
    sorted_pvals = pvals[sorted_ind]
    qmin = 1.0
    for i in range(n-1,-1,-1):
        qval = (n * sorted_pvals[i])/(i+1)
        if qval > qmin:
            qval = qmin
        else:
            qmin = qval
        qvals[sorted_ind[i]] = qval
    # Estimate number of true positives m1 and adjust q
    if n >= 233 and adjustment==True:
        phist = np.histogram(pvals, bins=20, range=(0, 1))[0]
        sorted_phist = np.sort(phist[3:19])
        if np.sum(sorted_phist) >= 160:
            median4 = n - 20*np.dot(np.array([1, 2, 2, 1]), sorted_phist[6:10])/6
            median6 = n - 20*np.dot(np.array([1, 2, 2, 2, 2, 1]), sorted_phist[5:11])/10
            m1 = min(median4, median6)

            qfac = (n - m1)/n
            if qfac < 0.5:
                qfac = 0.25 + qfac**2
            qvals *= qfac
    return qvals

# Z is dict with L,R, nan indicates invalid verts
# Returns q values for each hem
def FDR_z_hem(Z, adjustment = True):
    z_cat = np.concatenate((Z['L'], Z['R']))
    valid_inds = np.logical_not(np.isnan(z_cat))
    q_cat = np.ones(z_cat.shape[0])
    p_cat = np.ones(z_cat.shape[0])
    p_cat[z_cat>0] = norm.sf(z_cat[z_cat>0])
    p_cat[z_cat<0] = norm.cdf(z_cat[z_cat<0])
    q_cat[valid_inds] = FDR_p(p_cat[valid_inds], adjustment)
    print(p_cat)
    q = dict()
    q['L'] = q_cat[:Z['L'].shape[0]]
    q['R'] = q_cat[Z['R'].shape[0]:]

    return q

# D is SLs x nPerm+1 (first nSL_L are left hemisphere)
def SL_array_to_maps(D, savename, nSL_L, SLlist, nv, adjustment=True):
    hem = dict()
    hem['L'] = D[:nSL_L, :]
    hem['R'] = D[nSL_L:, :] 
    vox = SLtoVox(hem, SLlist, nv, zeronan=False)

    z_vox = dict()
    for hem in ['L', 'R']:
        z_vox[hem] = nullZ(vox[hem])

    q = FDR_z_hem(z_vox, adjustment)

    np.savetxt(savename + '.lh.dset', np.column_stack((np.nan_to_num(z_vox['L']), -np.log10(q['L']))))
    np.savetxt(savename + '.rh.dset', np.column_stack((np.nan_to_num(z_vox['R']), -np.log10(q['R']))))


# downsample event matrix to TR 
def regressor_to_TR(E, TR, nTR, E_dt=1):
    # E is event matrix with shape (time in sec, num of Events), TR is the 1.5, nTR is the number of TRs in a run
    nEvents = E.shape[1]

    # HRF (from AFNI)
    dt = np.arange(0, 15, E_dt)
    p = 8.6
    q = 0.547
    hrf = np.power(dt / (p * q), p) * np.exp(p - dt / q)

    # Convolve event matrix to get design matrix
    design_dt = np.zeros(E.shape)
    for e in range(nEvents):
        design_dt[:, e] = np.convolve(E[:, e], hrf)[:E.shape[0]]

    # Downsample event matrix to TRs
    timepoints = np.linspace(0, (nTR - 1) * TR, nTR)
    design = np.zeros((len(timepoints), nEvents))
    for e in range(nEvents):
        design[:, e] = np.interp(timepoints, np.arange(0, E.shape[0]*E_dt, E_dt),
                                 design_dt[:, e])
        design[:, e] = design[:, e] / np.max(design[:, e])

    return design


sub2subj = {"sub-01":"subj001", "sub-02":"subj002","sub-03":"subj003","sub-04":"subj005",
              "sub-05":"subj006", "sub-06":"subj007", "sub-07":"subj008", "sub-08":"subj009", 
           "sub-09":"subj010","sub-10":"subj011","sub-11":"subj013", "sub-12":"subj014", 
            "sub-15":"subj017", "sub-16":"subj018", "sub-17":"subj019", "sub-18":"subj020",
           "sub-19":"subj021", "sub-20":"subj022", "sub-21":"subj023", "sub-22":"subj024",'sub-23':'subj025',
           'sub-24':'subj026','sub-25':'subj027','sub-26':'subj029','sub-27':'subj031'}
ses2w = {"ses-01":"W2", "ses-02":"W4D1", "ses-03":"W4D2"}

nv = 40962

subjects = ['sub-%.2d'%s for s in range(1,28) if s!=13 and s!=14]
sessions = ['ses-%.2d'%s for s in range(1,3)]
runs = ['run-%s'%str(s) for s in range(1,3)]
TR = 1.5
nTRs = {'Item':302, 'Loci':302, 'Encode':355}

nTRs_w4d2 = {'Item': 156, 'Loci': 156, 'Encode': 182}
SL_lh = list(dd.io.load('SLlist_verydense.lh.h5').values())
SL_rh = list(dd.io.load('SLlist_verydense.rh.h5').values())

SLlist = {'L':SL_lh, "R": SL_rh}
nSL_L = len(SLlist['L'])