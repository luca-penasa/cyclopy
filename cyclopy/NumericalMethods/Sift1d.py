from scipy.interpolate import interp1d
import numpy as np

def gaussianBlur(signal, sigma):
    from astropy.convolution import Gaussian1DKernel, convolve
    g = Gaussian1DKernel(stddev=sigma)
    z = convolve(signal, g, boundary='extend')
    return z
    
    
    
def siresize(signal, factor):
    '''
    like imresize in matlab
    '''
    N = len(signal)
    finalN = np.floor(N * factor)
    
    x = np.arange(N)
    newx = np.arange(finalN)/factor
    
    i = interp1d(x,signal, fill_value="extrapolate")
    newy = i(newx)
    
    return newy


def builSigmaTable( Sigma0, OctaveNum, ScalesNum, k):
    Sigma = [[np.nan for j in np.arange(ScalesNum+3)] for i in np.arange(OctaveNum)]
    
    for o in np.arange(OctaveNum):
        print(o)
        if o==0:
            Sigma[o][0] =Sigma0 *0.5
        else:
            Sigma[o][0]= Sigma[o-1][0] * 2 

        for s in np.arange(ScalesNum+2) +1 :
             Sigma[o][s]=Sigma[o][s-1]*k**s;
    
    return Sigma
    
def buildBlurredSignalPyramid(y,  SigmaTable):
    
    OctavesNum = len(SigmaTable)
    ScalesNum = len(SigmaTable[0])
    
    Gaussian=[[np.nan for j in np.arange(ScalesNum)] for i in np.arange(OctavesNum)]

    for o in np.arange(OctavesNum):
        for s in np.arange(ScalesNum): 
            if o==0:              
                resized = siresize(y,2)
                blurred = gaussianBlur(resized,SigmaTable[o][s]);
                Gaussian[o][s]= blurred
            else:
                Gaussian[o][s]=siresize(gaussianBlur(y,SigmaTable[o][s]),0.5**(o-1));
                
    return Gaussian
    
def builDifferenceOfGaussians(Gaussians):
    OctavesNum = len(Gaussians)
    ScalesNum = len(Gaussians[0])
    
    DOG=[[np.nan for j in np.arange(ScalesNum-1)] for i in np.arange(OctavesNum)]

    for o in np.arange(OctavesNum):
        for s in np.arange(ScalesNum - 1):
            DOG[o][s]= Gaussians[o][s+1] - Gaussians[o][s]

    return DOG
    
def findPeaks(DOG):
    OctavesNum = len(DOG)
    ScalesNum = len(DOG[0])
    
    KeyPoint=[[np.nan for j in np.arange(ScalesNum)] for i in np.arange(OctavesNum)]

    for i,o in enumerate(DOG):
        from scipy.ndimage.filters import maximum_filter, minimum_filter
        from scipy.ndimage.morphology import generate_binary_structure, binary_erosion

        neighborhood = generate_binary_structure(2,2)

            #apply the local maximum filter; all pixel of maximal value 
            #in their neighborhood are set to 1
        local_max = maximum_filter(o, footprint=neighborhood)==o

        local_min = minimum_filter(o, footprint=neighborhood)==o
        all_peaks = local_min | local_max

        # all_peaks = np.delete(all_peaks, (0), axis=0)
        # all_peaks = np.delete(all_peaks, (-1), axis=0)

        #no border points
        all_peaks[:,0] = False
        all_peaks[0,:] = False
        all_peaks[:,-1] = False
        all_peaks[-1,:] = False

        KeyPoint[i] = all_peaks
        
    return KeyPoint

def doAccurateLocalizationAndFiltering(KeyPoint, DOG):
    print ("Init kpoints:")
    print ([np.sum(kp) for kp in KeyPoint])


    
    OctavesNum  = len(KeyPoint)
    ScalesNum = len(KeyPoint[0])


    #KeyPoint=[[np.nan for j in np.arange(ScalesNum)] for i in np.arange(OctavesNum)]

    
    from copy import deepcopy
    newKeyPoint = deepcopy(KeyPoint)

    shifts = []

    for o in np.arange(OctavesNum):
        current_octave = []
        for s in np.arange(ScalesNum):
            current_scale = []
            mykeys = KeyPoint[o][s]

            ids = np.where(mykeys == True)[0]

            for xi in ids:
                Dx=(DOG[o][s][xi+1] - DOG[o][s][xi-1]) / 2.0
                Ds=(DOG[o][s+1][xi] - DOG[o][s-1][xi]) / 2.0

                Dxx=DOG[o][s][xi+1] + DOG[o][s][xi-1] - 2*DOG[o][s][xi]           
                Dss=DOG[o][s+1][xi] + DOG[o][s-1][xi] - 2*DOG[o][s][xi]

                Dxs = (DOG[o][s+1][xi+1] + DOG[o][s-1][xi-1]- DOG[o][s-1][xi+1] - DOG[o][s+1][xi-1])/4.0
                Dsx = (DOG[o][s+1][xi+1] + DOG[o][s-1][xi-1]- DOG[o][s-1][xi+1] - DOG[o][s+1][xi-1])/4.0

                vec1 = np.array([[Dxx,Dxs],[Dss,Dsx]])
                vec2 = np.array([Dx,Ds])

                offset2 =  -np.dot( np.linalg.inv(vec1), vec2)
                prod2 = np.dot(vec2, offset2)
                D2 = DOG[o][s][xi] + 1/2* prod2


                # print (np.abs(offset2[0]))



                #filtering
                if abs(D2)<0.03:
                    newKeyPoint[o][s][xi]=False
                else:
                    current_scale.append(offset2[0])

            current_octave.append(current_scale)

        shifts.append(current_octave)

    print ("Filtered kpoints:")
    print ([np.sum(kp) for kp in newKeyPoint])
    
    return newKeyPoint, shifts

def concatenateShifts(shifts):
    all_shifts = []
    for o in shifts:
        for s in o:
            for i in s:
                all_shifts.append(i)

    all_shifts = np.array(all_shifts)
    return all_shifts


def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx


def getNearestValueSignal(xsignal, ysignal, xpos):
    id = find_nearest(xsignal, xpos)
    return ysignal[id]


def octaveScaleToSigmas(octaves, scales, sigma_table):
    oo = np.asarray(octaves, dtype=int)
    ss = np.asarray(scales, dtype=int)
    sigmas = []

    for o, s in zip(oo, ss):
        sigmas.append(sigma_table[o][s])

    return np.array(sigmas)


def keypointsToCoordinates(keypoints, original_signal):
    OctavesNum = len(keypoints)
    ScalesNum = len(keypoints[0]) - 2

    fakex = np.arange(len(keypoints[0][1]))
    out = []
    for o in np.arange(OctavesNum):
        if o == 0:
            factor = 2
        else:
            factor = 0.5 ** (o - 1)

        for s in np.arange(ScalesNum):
            myk = keypoints[o][s + 1]
            x_kps = np.arange(len(myk))[myk] / factor

            for x_kp in x_kps:
                val = getNearestValueSignal(fakex, original_signal, x_kp)
                out.append((x_kp, val, o, s))

    out = np.array(out)

    return out