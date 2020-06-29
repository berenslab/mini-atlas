import numpy as np
import time
import warnings
import pylab as plt

import glmnet_python
from glmnet import glmnet


###################################################
# Elastic net reduced-rank regression

def elastic_rrr(X, Y, rank=2, lambdau=1, alpha=0.5, max_iter = 100, verbose=0,
                sparsity='row-wise'):

    # in the pure ridge case, analytic solution is available:
    if alpha == 0:
         U,s,V = np.linalg.svd(X, full_matrices=False)
         B = V.T @ np.diag(s/(s**2 + lambdau*X.shape[0])) @ U.T @ Y
         U,s,V = np.linalg.svd(X@B, full_matrices=False)
         w = B @ V.T[:,:rank]
         v = V.T[:,:rank]

         pos = np.argmax(np.abs(v), axis=0)
         flips = np.sign(v[pos, range(v.shape[1])])
         v = v * flips
         w = w * flips

         return (w,v)

    # initialize with PLS direction
    _,_,v = np.linalg.svd(X.T @ Y, full_matrices=False)
    v = v[:rank,:].T
    
    loss = np.zeros(max_iter)
    
    for iter in range(max_iter):
        if rank == 1:
            w = glmnet(x = X.copy(), y = (Y @ v).copy(), alpha = alpha, lambdau = np.array([lambdau]), 
                       standardize = False, intr = False)['beta']
        else: 
            if sparsity=='row-wise':
                w = glmnet(x = X.copy(), y = (Y @ v).copy(), alpha = alpha, lambdau = np.array([lambdau]), 
                           family = "mgaussian", standardize = False, intr = False,
                           standardize_resp = False)['beta']
            else:
                w = []
                for i in range(rank):
                    w.append(glmnet(x = X.copy(), y = (Y @ v[:,i]).copy(), alpha = alpha, lambdau = np.array([lambdau]), 
                             standardize = False, intr = False, standardize_resp = False)['beta'])
            w = np.concatenate(w, axis=1)
                
        if np.all(w==0):
            v = v * 0
            return (w, v)
            
        A = Y.T @ X @ w
        a,c,b = np.linalg.svd(A, full_matrices = False)
        v = a @ b
        pos = np.argmax(np.abs(v), axis=0)
        flips = np.sign(v[pos, range(v.shape[1])])
        v = v * flips
        w = w * flips
        
        loss[iter] = np.sum((Y - X @ w @ v.T)**2)/np.sum(Y**2);        
        
        if iter > 0 and np.abs(loss[iter]-loss[iter-1]) < 1e-6:
            if verbose > 0:
                print('Converged in {} iteration(s)'.format(iter))
            break
        if (iter == max_iter-1) and (verbose > 0):
            print('Did not converge. Losses: ', loss)
    
    return (w, v)

def relaxed_elastic_rrr(X, Y, rank=2, lambdau=1, alpha=0.5, max_iter = 100,
                        sparsity='row-wise'):
    w,v = elastic_rrr(X, Y, rank=rank, lambdau=lambdau, alpha=alpha, 
                      sparsity=sparsity, max_iter=max_iter)

    if alpha==0:   # pure ridge: no need to re-fit
        return (w,v)

    nz = np.sum(np.abs(w), axis=1) != 0    
    wr,vr = elastic_rrr(X[:,nz], Y, rank=rank, lambdau=lambdau, alpha=0,
                        sparsity=sparsity, max_iter=max_iter)
    if np.sum(nz)>=np.shape(w)[1]:
        w[nz,:] = wr
        v = vr
    else:
        w[nz,:][:,:np.sum(nz)] = wr
        w[nz,:][:,np.sum(nz):] = 0
        v[:,:np.sum(nz)] = vr
        v[:,np.sum(nz):] = 0
    
    return (w,v)


###################################################
# Double biplot function
def bibiplot(X, Y, w, v, 
             YdimsNames=np.array([]), YdimsToShow=None,
             XdimsNames=np.array([]), XdimsToShow=None, 
             titles=[], xylim = 3,
             cellTypes=np.array([]), cellTypeColors={}, cellTypeLabels={},
             figsize=(9,4), axes=None):

    if XdimsToShow is None:
        nz = np.sum(np.abs(w), axis=1) != 0
        XdimsToShow = np.where(nz)[0]
    if YdimsToShow is None:
        nz = np.sum(np.abs(v), axis=1) != 0
        YdimsToShow = np.where(nz)[0]
    
    # Project and standardize
    Zx = X @ w[:,:2]
    Zy = Y @ v[:,:2]
    Zx = Zx / np.std(Zx, axis=0)
    Zy = Zy / np.std(Zy, axis=0)
    
    if not axes:
        plt.figure(figsize=figsize)
        plt.subplot(121, aspect='equal')
    else:
        plt.sca(axes[0])
    
    if cellTypes.size == 0:
        plt.scatter(Zx[:,0], Zx[:,1])
    else:
        for u in np.unique(cellTypes):
            if not cellTypeLabels:
                plt.scatter(Zx[cellTypes==u,0], Zx[cellTypes==u,1], color=cellTypeColors[u])
            else:
                plt.scatter(Zx[cellTypes==u,0], Zx[cellTypes==u,1], color=cellTypeColors[u], label=cellTypeLabels[u])
    plt.xlim([-xylim,xylim])
    plt.ylim([-xylim,xylim])
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])
    plt.xlabel('Component 1')
    plt.ylabel('Component 2')
    if titles:
        plt.title(titles[0])
    if cellTypeLabels:
        plt.legend(bbox_to_anchor=(1.35, 1.0))
        
    if XdimsToShow.size > 0:
        scaleFactor = 2
        L = np.corrcoef(np.concatenate((Zx[:,:2], X), axis=1), rowvar=False)[2:,:2]
        for i in XdimsToShow:
            plt.plot([0, scaleFactor*L[i,0]], [0, scaleFactor*L[i,1]], linewidth=1, color=[.4, .4, .4])
            plt.text(scaleFactor*L[i,0]*1.2, scaleFactor*L[i,1]*1.2, XdimsNames[i], 
                     ha='center', va='center', color=[.4, .4, .4], fontsize=10)
        circ = plt.Circle((0,0), radius=scaleFactor, color=[.4, .4, .4], fill=False, linewidth=1)
        plt.gca().add_patch(circ)

    if not axes:
        plt.subplot(122, aspect='equal')
    else:
        if not axes[1]:
            return
        plt.sca(axes[1])
        
    if cellTypes.size == 0:
        plt.scatter(Zy[:,0], Zy[:,1])
    else:
        for u in np.unique(cellTypes):
            plt.scatter(Zy[cellTypes==u,0], Zy[cellTypes==u,1], color=cellTypeColors[u])
    plt.xlim([-xylim,xylim])
    plt.ylim([-xylim,xylim])
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])
    plt.xlabel('Component 1')
    plt.ylabel('Component 2')
    if titles:
        plt.title(titles[1])
    plt.tight_layout()

    if YdimsToShow.size > 0:
        scaleFactor = 2
        L = np.corrcoef(np.concatenate((Zy[:,:2], Y), axis=1), rowvar=False)[2:,:2]
        for i in YdimsToShow:
            plt.plot([0, scaleFactor*L[i,0]], [0, scaleFactor*L[i,1]], linewidth=1, color=[.4, .4, .4])
            plt.text(scaleFactor*L[i,0]*1.2, scaleFactor*L[i,1]*1.2, YdimsNames[i], 
                     ha='center', va='center', color=[.4, .4, .4], fontsize=10)
        circ = plt.Circle((0,0), radius=scaleFactor, color=[.4, .4, .4], fill=False, linewidth=1)
        plt.gca().add_patch(circ)
        
        
###################################################
# Permutation procedures to estimate dimensionality
def dimensionality(X, Y, nrep = 100, seed = 42, axes=None, figsize=(7,2)):

    np.random.seed(seed)

    _,spectrum,_ = np.linalg.svd(X, full_matrices=False)
    spectra = np.zeros((nrep, spectrum.size))
    for rep in range(nrep):
        Xperm = X.copy()
        for i in range(Xperm.shape[1]):
            Xperm[:,i] = Xperm[:,i][np.random.permutation(Xperm.shape[0])]
        _, spectra[rep,:], _ = np.linalg.svd(Xperm, full_matrices=False)
    
    if not axes:
        plt.figure(figsize=figsize)
        plt.subplot(131)
    else:
        plt.sca(axes[0])
    plt.plot(np.arange(1, spectrum.size), spectra[:,:-1].T**2/np.sum(spectrum**2), 'k', linewidth=1)
    plt.plot(np.arange(1, spectrum.size), spectrum[:-1]**2/np.sum(spectrum**2), '.-')
    dimX = np.where(spectrum < np.percentile(spectra, 95, axis=0))[0][0]
    plt.text(plt.xlim()[1]*.2, plt.ylim()[1]*.8, 'X dimensionality: ' + str(dimX), fontsize=8)

    _,spectrum,_ = np.linalg.svd(Y, full_matrices=False)
    spectra = np.zeros((nrep, spectrum.size))
    for rep in range(nrep):
        Xperm = Y.copy()
        for i in range(Xperm.shape[1]):
            Xperm[:,i] = Xperm[:,i][np.random.permutation(Xperm.shape[0])]
        _, spectra[rep,:], _ = np.linalg.svd(Xperm, full_matrices=False)

    showy = True
    if not axes:
        plt.subplot(132)
    else:
        if axes[1]:
            plt.sca(axes[1])
        else:
            showy = False
    if showy:
        plt.plot(np.arange(1, spectrum.size), spectra[:,:-1].T**2/np.sum(spectrum**2), 'k', linewidth=1)
        plt.plot(np.arange(1, spectrum.size), spectrum[:-1]**2/np.sum(spectrum**2), '.-')
        dimY = np.where(spectrum < np.percentile(spectra, 95, axis=0))[0][0]
        plt.text(plt.xlim()[1]*.2, plt.ylim()[1]*.8, 'Y dimensionality: ' + str(dimY), fontsize=8)

    Xz,_,_ = np.linalg.svd(X, full_matrices=False)
    Xz = Xz[:,:dimX]
    yhat = Xz @ Xz.T @ Y
    _,spectrum,_ = np.linalg.svd(yhat, full_matrices=False)
    spectra = np.zeros((nrep, spectrum.size))
    for rep in range(nrep):
        Xz = Xz[np.random.permutation(Xz.shape[0]),:]
        yhat = Xz @ Xz.T @ Y
        _, spectra[rep,:], _ = np.linalg.svd(yhat, full_matrices=False)
    
    if not axes:
        plt.subplot(133)
    else:
        plt.sca(axes[2])
    p = min(dimX, Y.shape[1])
    plt.plot(np.arange(1, p+1), spectra[:,:p].T**2/np.sum(spectrum**2), 'k', linewidth=1)
    plt.plot(np.arange(1, p+1), spectrum[:p]**2/np.sum(spectrum**2), '.-')
    dimRRR = np.where(spectrum > np.percentile(spectra, 95, axis=0))[0][-1]
    plt.text(plt.xlim()[1]*.2, plt.ylim()[1]*.8, 'RRR dimensionality: ' + str(dimRRR), fontsize=8)

    plt.tight_layout()

    
###################################################
# Cross-validation for elastic net reduced-rank regression
def elastic_rrr_cv(X, Y, alphas = np.array([.2, .5, .9]), lambdas = np.array([.01, .1, 1]), 
                   reps=10, folds=10, rank=1, seed=42, sparsity='row-wise'):
    n = X.shape[0]
    r2 = np.zeros((folds, reps, len(lambdas), len(alphas))) * np.nan
    r2_relaxed = np.zeros((folds, reps, len(lambdas), len(alphas))) * np.nan
    corrs = np.zeros((folds, reps, len(lambdas), len(alphas), rank)) * np.nan
    corrs_relaxed = np.zeros((folds, reps, len(lambdas), len(alphas), rank)) * np.nan
    nonzero = np.zeros((folds, reps, len(lambdas), len(alphas))) * np.nan

    # CV repetitions
    np.random.seed(seed)
    t = time.time()
    for rep in range(reps):
        print(rep+1, end='')
        ind = np.random.permutation(n)
        X = X[ind,:]
        Y = Y[ind,:]
        
        # CV folds
        for cvfold in range(folds):
            print('.', end='')

            indtest  = np.arange(cvfold*int(n/folds), (cvfold+1)*int(n/folds))
            indtrain = np.setdiff1d(np.arange(n), indtest)
            Xtrain = np.copy(X[indtrain,:])
            Ytrain = np.copy(Y[indtrain,:])
            Xtest  = np.copy(X[indtest,:])
            Ytest  = np.copy(Y[indtest,:])
            
            # mean centering
            X_mean = np.mean(Xtrain, axis=0)
            Xtrain -= X_mean
            Xtest  -= X_mean
            Y_mean = np.mean(Ytrain, axis=0)
            Ytrain -= Y_mean
            Ytest  -= Y_mean
            
            # loop over regularization parameters
            for i,a in enumerate(lambdas):    
                for j,b in enumerate(alphas):
                    vx,vy = elastic_rrr(Xtrain, Ytrain, lambdau=a, alpha=b, rank=rank, sparsity=sparsity)
                    
                    nz = np.sum(np.abs(vx), axis=1) != 0
                    if np.sum(nz) < rank:
                        continue

                    if np.allclose(np.std(Xtest @ vx, axis=0), 0):
                        continue
                    
                    nonzero[cvfold, rep, i, j] = np.sum(nz)
                    r2[cvfold, rep, i, j] = 1 - np.sum((Ytest - Xtest @ vx @ vy.T)**2) / np.sum(Ytest**2)
                    for r in range(rank):
                        corrs[cvfold, rep, i, j, r] = np.corrcoef(Xtest @ vx[:,r], Ytest @ vy[:,r], rowvar=False)[0,1]
                        
                    # Relaxation
                    vxr,vyr = elastic_rrr(Xtrain[:,nz], Ytrain, lambdau=a, alpha=0, rank=rank, sparsity=sparsity)
                    if np.sum(nz)>=np.shape(vy)[1]:
                        vx[nz,:] = vxr
                        vy = vyr
                    else:
                        vx[nz,:][:,:np.sum(nz)] = vxr
                        vx[nz,:][:,np.sum(nz):] = 0
                        vy[:,:np.sum(nz)] = vyr
                        vy[:,np.sum(nz):] = 0

                    if np.allclose(np.std(Xtest @ vx, axis=0), 0):
                        continue

                    r2_relaxed[cvfold, rep, i, j] = 1 - np.sum((Ytest - Xtest @ vx @ vy.T)**2) / np.sum(Ytest**2)
                    for r in range(rank):
                        corrs_relaxed[cvfold, rep, i, j, r] = np.corrcoef(Xtest @ vx[:,r], Ytest @ vy[:,r], rowvar=False)[0,1]
                    
        print(' ', end='')
    
    t = time.time() - t
    m,s = divmod(t, 60)
    h,m = divmod(m, 60)
    print('Time: {}h {:2.0f}m {:2.0f}s'.format(h,m,s))
    
    return r2, r2_relaxed, nonzero, corrs, corrs_relaxed


###################################################
# Bootstrap selection for elastic net reduced-rank regression
def elastic_rrr_bootstrap(X, Y, rank=1, lambdau = 1.5, alpha = .5, nrep = 100, seed=42):
    np.random.seed(seed)
    ww = np.zeros((X.shape[1], nrep))
    for rep in range(nrep):
        print('.', end='')
        n = np.random.choice(X.shape[0], size = X.shape[0])
        w,v = elastic_rrr(X[n,:], Y[n,:], rank = rank, lambdau = lambdau, alpha = alpha)
        ww[:,rep] = w[:,0]
    print(' ')
    bootCounts = np.sum(ww!=0, axis=1)/nrep
    return bootCounts

####################################################
# Plot CV results
def plot_cv_results(r2=None, r2_relaxed=None, nonzeros=None, corrs=None, corrs_relaxed=None, alphas=None):
    
    # suppressing "mean of empty slice" warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        n = np.nanmean(nonzeros, axis=(0,1))
        cr = np.nanmean(r2_relaxed, axis=(0,1))
        c = np.nanmean(r2, axis=(0,1))
        c1 = np.nanmean(corrs_relaxed, axis=(0,1))[:,:,0]
        if corrs_relaxed.shape[4]>1:
            c2 = np.nanmean(corrs_relaxed, axis=(0,1))

    plt.figure(figsize=(9,4))
    plt.subplot(121)
    plt.plot(n, cr, '.-', linewidth=1)
    plt.gca().set_prop_cycle(None)
    plt.plot(n, c, '.--', linewidth=1, alpha=.5)
    plt.xscale('log')
    plt.xlabel('Number of non-zero genes')
    plt.ylabel('Test R2')
    plt.legend(['$\\alpha='+str(a)+'$' for a in alphas])

    plt.subplot(122)
    plt.plot(n, c1, '.-', linewidth=1)
    if corrs_relaxed.shape[4]>1:
        for k in range(1, corrs_relaxed.shape[4]):
            plt.gca().set_prop_cycle(None)            
            plt.plot(n, c2[:,:,k], '.--', linewidth=1)
    plt.xscale('log')
    plt.xlabel('Number of non-zero genes')
    plt.ylabel('Correlations')
    plt.legend(alphas)
    plt.legend(['$\\alpha='+str(a)+'$' for a in alphas])
    plt.tight_layout()

