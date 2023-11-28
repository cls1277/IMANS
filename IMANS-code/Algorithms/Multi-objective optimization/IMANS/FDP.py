import numpy as np
import matplotlib.pyplot as plt

class FastDensityPeaks():


    def __init__(self, verbose=False, autoassign=True, sample_cutoff=False):
        """

        Parameters:
        ----------
        - verbose: if true, the decision graph is shown aftern the fitting
        - autoassign: if true, the centroids are authomatically retrieved from
        the decision graph after the fitting
        - sample_cutoff: if true, the cutoff is retrieved using a sample of the 
        distance matrix

        """
        self.verbose    = verbose
        self.autoassign = autoassign
        self.sample_cutoff = sample_cutoff
        self.centroids  = []


    def fit(self, X, fraction=0.04):
        """

        Retrieve rho, delta and nearest neighbours from the distance matrix of X

        Parameters:
        ----------
        - X (2-dimensional numpy array): the dataset to cluster
        - fraction (float): the fraction of neighbours for retrieving the cutoff dc


        """
        self.N = X.shape[0] # Number of points
        if self.N > 2000:
            self.sample_cutoff = True # For large datasets, get the cutoff from a sample
        self.fraction = fraction 
        self.e_matrix = self._get_euclidean_matrix(X.astype(np.float32)) # Retrieve euclidean matrix
        self.dc  = self._get_cutoff() # Retrieve the cutoff from the fraction of nn
        self.rho = self._get_rho() # Retrieve rho
        self.delta, self.nns = self._get_delta() # Retrieve delta

        if self.autoassign:  # If autoassign, retrieve the centroids from rho/delta _plot
            self.centroids = self._get_centroids_()
            self._assign()

        if self.verbose:  # If verbose, _plot delta as a function of rho
            self._plot()
         
        
    def fit_predict(self, X, fraction=0.04):
        """
        
        Function that calls the fit method and returns the labels (it uses autoassign)
        
        """
        self.fit(X, fraction)
        self.centroids = self._get_centroids_()
        self._assign()
        return self.labels_


    def decision_graph(self, min_rho, min_delta):
        """

        Function to retrieve the centroids manually from the decision graph

        Parameters
        ----------
        - min_rho (float): defines minumum rho value for a point to be considered a centroid
        - min_delta (float): defines minumum delta value for a point to be considered a centroid

        """    
        # Find the centrois from the graph decision
        self.centroids = np.intersect1d(np.where(self.rho > min_rho)[0], np.where(self.delta > min_delta)[0],
                                        assume_unique=True)

        if self.verbose: # If verbose, _plot the decision graph with the grid defined by the user
            self._plot(min_rho, min_delta)
            print('Cluster indexes: ',self.centroids)

        self._assign() # _assign the points to the clusters of the centroids


    def get_halo(self):
        """

        Function that returns the labels with halo (-1)

        """
        self.halo  = self.labels_.copy()
        self.rho_b = np.zeros(self.n_clusters) # the density at the border

        for c in range(self.n_clusters):

            # Find the maximum rho of the boarder points
            cluster_points    = np.where(self.labels_==c)[0]
            no_cluster_points = np.where(self.labels_!=c)[0]

            for i in cluster_points:
                less_dc_ncp = no_cluster_points[np.where(self.e_matrix[i, no_cluster_points]<self.dc)[0]]
                rho_ncp     = self.rho[less_dc_ncp]
                if len(rho_ncp)>0:
                    rho_max = np.max((rho_ncp+self.rho[i])/2.)
                    if rho_max > self.rho_b[c]:
                        self.rho_b[c] = rho_max

            # _assign the points with density inferior than the highest to the halo
            halo_idx = cluster_points[np.where(self.rho[cluster_points] < self.rho_b[c])]
            self.halo[halo_idx] = -1
            
        return self.halo

            

    def _assign(self):
        """

        Function that assigns to each point the label of its nearest neighbour

        """
        if len(self.centroids)==0:
            raise ValueError("You first need to define the centroids! Use the decision_graph\
                              function or set autoassign=True when fit the data")

        self.n_clusters = self.centroids.shape[0]
        self.rho_sort_idxs = np.argsort(self.rho)[::-1] # The order of rho defines the order of assignment
        self.labels_ = np.full(self.N, -1)
        self.labels_[self.centroids] = [n for n in range(self.n_clusters)] # _assign the initial labels_ to the centroids
        for r in self.rho_sort_idxs: # slide along the sorted indexes to _assign to each label the one of its nns
            if self.labels_[r] == -1:
                self.labels_[r] = self.labels_[self.nns[r]] 



    def _plot(self, min_rho=None, min_delta=None):
        """

        Function that plots delta as a function of rho. Automatically called when verbose=True

        """
        plt.scatter(self.rho, self.delta)
        if min_rho and min_delta:
            plt.plot([min_rho, self.rho.max()], [min_delta, min_delta], linewidth=1, color="green")
            plt.plot([min_rho, min_rho], [min_delta, self.delta.max()], linewidth=1, color="green")
        if len(self.centroids) > 0:
            plt.scatter(self.rho[self.centroids], self.delta[self.centroids])
        plt.title('Decision Graph')
        plt.xlabel(r'$\rho$')
        plt.ylabel(r'$\delta$')
        plt.show()


    def _get_euclidean_matrix(self, X):
        """

        Function that computes the (euclidean) distance matrix

        Parameters:
        ----------
        - X (2-dimensional numpy array): the dataset to cluster

        """
        squarred = np.einsum('ij,ij->i', X, X)[:, np.newaxis] # equivalent to np.sum(X**2) but more optimized
        dist = -2 * X @ X.T # dot product optimized for 2-d dimension
        dist += squarred 
        dist += squarred.T
        np.maximum(dist, 0, out=dist)
        np.fill_diagonal(dist, 0) # to avoid small errors and sqrt warnings
        return np.sqrt(dist, out=dist)


    def _get_cutoff(self):
        """

        Function that retrieves the radious of the gaussian kernel (cutoff dc)

        """
        # in the case in which we have a lot of instnaces, the algorithm is performed
        # in a sumbsample of the distances
        if self.sample_cutoff: 
            sorted_distances = np.random.choice(self.e_matrix.ravel(),self.N) # sample N distances
            sorted_distances.sort()
            
        # In the other cases, the cutoff is retrieved with more precision extracting
        # the triangular upper matrix (without diagonal) of the distnaces (1-d array) and sorting it
        else:
            sorted_distances = self.e_matrix[np.triu_indices_from(self.e_matrix, k=1)]
            sorted_distances.sort()
            
        # Now that we have all the distances sorted, we can take one that is around the fraction required
        index = round(0.5+sorted_distances.size*self.fraction) # 0.5 to avoid 0
        return sorted_distances[index]


    def _get_rho(self):   
        """

        Funtion that retrieves the local density using a gaussian kernel with radious dc

        """  
        exp_matrix = np.exp(-(self.e_matrix/self.dc)**2) # Vectorized normalized distances
        return np.sum(exp_matrix, axis=1) - 1 # -1 because the 0 in the diagonal becomes 1 with exp

    
    def _get_delta(self):
        """

        Function that retrieves deltas as the minimum  distance between a point i and any other point
        with higher density (and the nearest neighbour as the indexes of those points in rho)

        """        
        delta = np.zeros(self.N)
        nearest_neighbours = np.zeros(self.N, dtype=np.int64)
        higher_index = 0

        for i in range(self.N):
            
            indexes = np.where(self.rho > self.rho[i])[0] # Indexes of points with higher density than i
            distances = self.e_matrix[i] # Distances of the points with higher density than i
            
            # In general, take the mimimum distance; if no indexes are found, we have 
            # the point with higher density: conventionally take delta = max(d_ij)
            if indexes.size != 0:
                idx  = np.argmin(distances[indexes])
                delta[i] = distances[indexes][idx]
                nearest_neighbours[i] = indexes[idx]
            else:
                higher_index = i
                
        #self.higher_index = higher_index # Save the higher index   
        delta[higher_index] = np.max(delta) # assign to the higher_index the maximum value of delta        
        return delta, nearest_neighbours

    def _get_centroids_(self):
        """

        Function that automatically retrieve the centroids as outlier of decision graph

        """    
        mean = np.mean(self.delta)
        std  = np.std(self.delta)
        z_scores = (self.delta-mean)/std # get the z_scores values in respect to delta
        idx = np.argsort(z_scores)[::-1] # get the indexes of the sorted z_scores
        
        # in the case of a tie, take the as first centroid the element with highest rho
        if self.delta[idx[0]]==self.delta[idx[1]]: 
            if self.rho[idx[1]] > self.rho[idx[0]]:
                t = idx[0]
                idx[0] = idx[1]
                idx[1] = t
    
        centroids = [idx[0]] # The point with highest delta is considered for sure a centroid

        # until when there is a point with z_score at half of the last centroids,
        # consider that point a centroid
        for i in range(1, idx.size): 
            if 2*z_scores[idx[i]] > z_scores[centroids[-1]]: # Discard all the points that are too far from the last centroid
                if (self.delta[idx[i]] == self.delta[idx[0]] or self.delta[idx[i]]>np.pi*(mean+std))\
                    and np.pi*self.rho[idx[i]]>self.rho[idx[0]]: # consider centroids only the points around the first centroid
                    centroids.append(idx[i])
            else:
                break
                
        return np.array(centroids)




