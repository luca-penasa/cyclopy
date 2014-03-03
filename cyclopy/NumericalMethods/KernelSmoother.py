from __future__ import division
import numpy as np


class KernelSmoother:
    """
    Kernel Smoothing as in "The elements of statistical learning: data mining, inference, and prediction" by Hastie et al. (2009)
    by Springer editions. KdTree can be used for big-to-huge datasets. This option is enabled by default.
    """

    def __init__(self, x, y, use_kdtree=True):
        """
        Parameters
        ----------
        x: array-like with shape (1,N)
            x coordinates
        y: array-like with shape (1, N)
            variable to be smoothed
        use_kdtree: boolean
            specify if the kdtree structure need to be build (if False no kdtree can be used)
        """

        #Add data to self
        sort_ids = np.argsort(x)
        self.x = x[sort_ids]
        self.y = y[sort_ids]

        #Some infos
        self.N = len(self.x)
        self.x_max = self.x[-1]
        self.x_min = self.y[0]
        self.x_range = self.x_max - self.x_min
        self.mean_density = self.N / self.x_range

        #if kdtree is requested or not
        self.use_kdtree = use_kdtree

        if self.use_kdtree == True:
            from scipy.spatial import cKDTree

            print("-> Building KDTree indexing\n")
            self.locator = cKDTree(np.array([self.x]).T)
            print("   ... Done\n")

    def extract_neigh(self, x_pos, distance):
        """
        Given a x coordinate "x_pos" extract points within the distance "distance" using kdtree
        Note that the maximum number of points to retrieve is limited to "distance * mean_density * 20".
        20 is an arbitrary chosen constant, if yor dataset shows big variation in point density maybe you
        should increase it. See initialization function for see how density is computed.
        """

        ndists, nids = self.locator.query([x_pos], k=10000, p=2, distance_upper_bound=distance)
        return ndists[ndists != np.inf], nids[ndists != np.inf]

    def _gaussian(self, x):
        """
        gaussian kernel
        """
        result = 1 / np.sqrt(2 * np.pi) * np.exp(-0.5 * x ** 2)
        return result

    def _epanechnikov(self, x):
        """
        Epanechnikov kernel. "x" MUST be a numpy array
        """
        epanech = lambda x: 3 / 4 * (1 - x ** 2)
        result = np.zeros(len(x))
        mask = np.abs(x) <= 1
        result[mask] = epanech(x[mask])

        return result

    def __call__(self, x, h=1, method='epa', disable_kdtree=False, just_n_count=False, just_variance=False):
        """
        Evaluate the Kernel Smoother

        Parameters
        ----------
        x: array [1, n]
            positions at which evaluate the smoother

        h: positive real
            scaling factor (support radius or bandwidth)

        method: string
            can be: - "gau" for gaussian kernel
                    - "epa" for epanechnikov kernel
            Note! The effective support region for a gaussian is 2/3 times the support region of epanechnikov, so
            computing two series with the same h but different kernels shows different smoothness (around the double).

        disable_kdtree: bool
            permits to temporary disable the kdtrre
            Have no effects if the kdtree was not build during class initialization
        """
        import progressbar as pbar

        pb = pbar.ProgressBar()

        old_kdtree_status = self.use_kdtree  #was kdtree build?

        #temporary kdtree disabling
        if (disable_kdtree == True) & (self.use_kdtree == False):
            print('You are already NOT using the kdtree structure for nn search!')
        elif (disable_kdtree == True) & (self.use_kdtree == True):
            print(
                'Kdtree searching disabled for this evaluation! By default it would be used. Next run will use it if you dont specifically disable it.')
            self.use_kdtree = False

        #kernel selection
        if method == 'gau':
            self.kernel = self._gaussian
            self.k_support = 3  #we will need to catch all values inside a 3*sigma distance (approximation of support region)
        elif method == 'epa':
            self.kernel = self._epanechnikov
            self.k_support = 1
        else:
            raise TypeError ("You must give an existent kernel name: gaussian -> 'gau' or epanechnikov -> 'epa' are supported")

        #evaluator selcetion, depending on kdtree use or not
        if self.use_kdtree == True:
            self.evaluator = self.evaluate_single_kdtree
        else:
            self.evaluator = self.evaluate_single

        if (just_n_count == True) & (self.use_kdtree == True):
            self.evaluator = self.neighbor_counter_kdtree

        elif (just_n_count == True) & (self.use_kdtree == False):
            print("you must enable kdtree for neighbors counting. N-count without kdtree is not supported")

            if (just_variance == True) & (self.use_kdtree == True):
                self.evaluator = self.variance_evaluator

        elif (just_variance == True) & (self.use_kdtree == False):
            print("you must enable kdtree for neighbors counting. N-count without kdtree is not supported")

            #Ensure "x" is an array
        if type(x) != np.ndarray:
            x = np.array([x])

        #Compute result
        result = np.zeros(len(x))
        for i in pb(np.arange(len(x))):
            result[i] = self.evaluator(x[i], h)

        #reset kdtree usage
        self.use_kdtree = old_kdtree_status
        return result

    def neighbor_counter_kdtree(self, x, h):
        #old version: was simple counting - not smoothing
        #r_support = h / 2
        #d, nids = self.extract_neigh(x, r_support)
        #result = len(d)

        r_support = h * self.k_support
        d, nids = self.extract_neigh(x, r_support)
        d = d / h
        nx = self.x[nids]
        ny = np.ones(lenght(nx))  #each point have value 1!

        w = self.kernel(d)
        result = np.sum(w * ny) / np.sum(w)

        return result

    def variance_evaluator(self, x, h):
        r_support = h / 2
        d, nids = self.extract_neigh(x, r_support)
        result = np.var(self.y[nids])
        return result

    def evaluate_single(self, x, h):
        """
        evaluate Kernel Smoother at a given unique "x"
        """

        d = np.abs(self.x - x) / h
        w = self.kernel(d)
        result = np.sum(w * self.y) / np.sum(w)
        return result

    def evaluate_single_kdtree(self, x, h):
        """
        evaluate Kernel Smoother at a given unique "x", using kdtree search structure
        """
        r_support = h * self.k_support
        d, nids = self.extract_neigh(x, r_support)
        d = d / h
        nx = self.x[nids]
        ny = self.y[nids]
        w = self.kernel(d)
        result = np.sum(w * ny) / np.sum(w)
        return result



def smooth_signal(x, y, bandwidth, kdtree=False):
    ks = KernelSmoother(x, y, kdtree)
    return ks(x, bandwidth)
