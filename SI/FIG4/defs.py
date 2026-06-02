import numpy as np
from numba import njit

@njit(cache=True)
def numbified_version(dummy_data,r):
    C = np.corrcoef(dummy_data.T)
    C = 1/2*(C+C.T) #force symmetrization to ensure real eigenvalues
    eigvals,_ = np.linalg.eigh(C)
    eigvals = np.real(eigvals)
    if r>1:
        eigvals = 1/r*eigvals[eigvals>1e-12] #keep nonzero eigenvalues, properly rescaled.
    
    
    eigvals = np.sort(eigvals) #from smallest to largest
    lam_plus = eigvals[-1]
    lam_minus = eigvals[0]
    return eigvals,lam_minus,lam_plus

print('Importing functions...')
def check_significant_eigenvalues(data):
    A_pair,V_pair,others = data
    r,L,v,M,lam_null_minus,lam_null_plus = others
    ii,A = A_pair
    jj,V = V_pair
    print(ii,jj,'A={}'.format(A),'V={}'.format(V))
    Sigma = np.zeros((L+1,L+1)) #the covariance matrix
    den = V+A+L*v
    Sigma[0,0] = 1 -A/den
    Sigma[0,1:] = -np.sqrt(v*A)/den
    Sigma[1:,0] = -np.sqrt(v*A)/den
    Sigma[1:,1:] = - v/den  #also fills diagonal but it's ok bc we overwrite it in the next line
    np.fill_diagonal(Sigma[1:,1:],
                     1-v/den)
    dummy_data = np.random.multivariate_normal([0]*(L+1),Sigma,size=M)

    # eigvals,_ =np.linalg.eig(np.corrcoef(dummy_data.T))
    # eigvals = np.real(eigvals)
    # if r>1:
    #     eigvals = 1/r*eigvals[eigvals>1e-10] #keep nonzero eigenvalues, properly rescaled.
    
    
    # eigvals = np.sort(eigvals) #from smallest to largest
    # lam_plus = eigvals[-1]
    # lam_minus = eigvals[0]
    
    eigvals,lam_minus,lam_plus = numbified_version(dummy_data,r)
    print(lam_minus,lam_plus,lam_null_minus,lam_null_plus)
    if lam_minus<lam_null_minus:
        if lam_plus>lam_null_plus:
            ''' Both eigenvalues significant '''
            modes =  +2
        else:
            ''' Only low-var eigenvalue from total comm. abu significant '''
            modes = +1
    else:
        if lam_plus>lam_null_plus:
            ''' Only high-var eigenval from A significant '''
            modes = 0 #will not happen anyways
        else:
            ''' None are significant '''
            modes = 0
    print('ii,jj,modes',ii,jj,modes)
    return ii,jj,modes,eigvals

def get_null_eigenvalues(data):
    (L,M,r,v) = data
    dummy_data = np.random.multivariate_normal([0]*(L+1),v*np.eye(L+1),size=M)
    eigvals,_ =np.linalg.eigh(np.corrcoef(dummy_data.T))
    eigvals = np.real(eigvals)
    if r<1:
        pass
    else:
        ''' More species than samples. Append only non-zero eigenvalues, scaled by the factor 1/r such that one can compare with the results from the later matrix '''
        eigvals  = 1/r*eigvals[eigvals>1e-12]
    return eigvals


print('Functions Imported')
