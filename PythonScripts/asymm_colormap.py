import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def normalize(x,a,b): #normalization map: [a,b] --> [0,1]
    if a>b:
        raise ValueError('(a,b) is not an interval')
    return (float(x)-a)/(b-a)

def display_cmap(cmap):
    plt.imshow(np.linspace(0, 100, 256)[None, :],  aspect=25, interpolation='nearest', cmap=cmap) 
    plt.axis('off')

def asymmetric_cmap(data, div_cmap, ref_point=0.0, name= 'asym_cmap'):
    '''
    Input
    -----
      data: data to be visualized (a numpy aray of shape (m,n), a data frame, a list of lists of equal len)
      div_cmap :  a diverging matplotlib or seaborn colormap  (a matplotlib.colors.LinearSegmentedColormap object)
      ref_point is the reference point for data, the threshold of interest  
    '''
    if isinstance(data, pd.DataFrame):
        D = data.values
    elif isinstance(data, np.ma.core.MaskedArray):
        D=np.ma.copy(data)
    else:    
        D=np.asarray(data, dtype=float) 
    D=np.ma.masked_invalid(D)
    
    dmin=np.min(D)
    dmax=np.max(D)
    
    if not (dmin < ref_point < dmax):
        raise ValueError('data are not appropriate for visualization with a diverging colormap')
        
    if dmax-ref_point > ref_point-dmin:
        left=2*ref_point-dmax
        right=dmax
        
        tp=normalize(dmin, left,right)#normalized value of dmin
        refp_norm=normalize(ref_point, left, right) # normalized value of the ref_point in the symmetric interval
                                                   #[left, right]. It is 0.5 
       
        A=tp
        B=1.0
    else:
        left=dmin
        right=2*ref_point-dmin
        
        tp=normalize(dmax, left, right)#normalized value of dmax
        refp_norm=normalize(ref_point, left, right)
        
        A=0.0
        B=tp
    max_lumin_idx=normalize(refp_norm, A, B) # index for the max luminance position in the asymm div cmap
    
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # T is a (256,)  array 
    T = np.hstack([
        np.linspace(A, refp_norm, 128, endpoint=False), 
        np.linspace(refp_norm, B, 128)
    ])
   
    # T_assym is (256,) array 
    T_asymm = np.hstack([
        np.linspace(0, max_lumin_idx, 128, endpoint=False), 
        np.linspace(max_lumin_idx, 1.0, 128)
    ])

    for t, s in zip(T, T_asymm):
        r, g, b, a = div_cmap(t)

        cdict['red'].append((s, r, r))
        cdict['green'].append((s, g, g))
        cdict['blue'].append((s, b, b))
        cdict['alpha'].append((s, a, a))

    asym_cmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=asym_cmap)
    display_cmap(asym_cmap)
    return D, asym_cmap, dmin, dmax

