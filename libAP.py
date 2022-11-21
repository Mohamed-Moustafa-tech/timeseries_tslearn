import numpy as np
import operator
#from sklearn.cluster import KMeans
import graph_tool as gt

from utils import *
from tslearn.clustering import TimeSeriesKMeans
import tslearn.metrics
from tslearn.barycenters import euclidean_barycenter, dtw_barycenter_averaging, softdtw_barycenter
from tslearn.clustering.kmeans import _k_init_metric as centroid_init
from tslearn.clustring.utils import _compute_inertia

flatten = lambda l: [item for sublist in l for item in sublist]


class LSOprimizer:
    def __init__(self, GE, G, pat,cluster_centers= None, metric="dtw", L_min, L_max, T = 20, max_iter=100, plot=True, opt_pat=None, k=2,
                 init_size=None, seed=None, verbose=True):
        """
        Given a graph G and gene expression array GE finds the optimal subnetwork in G of size at least L_min and
        at most L_max that can provide the optimal patients clustering in k clusters.
        :param GE: pandas DatFrame with gene expression
        :param G: graphml graph with PPI network
        :param cluster_centers: centroids with default as None
        :param L_min: minimal desired solution subnetwork size
        :param L_max: maximal desired solution subnetwork size
        :param T: temprature parameter for SA
        :param max_iter: maximal allowed number of iterations
        :param plot: convergence plot (True/False)
        :param opt_pat: patients labels (if provided, patients clustering won't be performed
        :param k: number of clusters
        :param init_size: initial subnetwork size (default L_max *2)
        :param seed: seed
        :param verbose: True/False
        """
        self.G = gt.load_graph(G)# changed fom nxgt(G)
        self.T = T
        self.L_min = L_min
        self.L_max = L_max
        self.max_iter = max_iter
        self.plot = plot
        if opt_pat is None:
            self.opt_pat = []
        else:
            self.opt_pat = opt_pat
        self.k = k
        if init_size is None:
            self.init_size = L_max * 2
        else:
            self.init_size = init_size
        self.seed = seed
        self.verbose = verbose
        self.ge = GE
        self.genes = list(self.G.get_vertices())
        self.patients = np.array(list(GE.columns))
        self.metric=metric

    def APUtil(self, u, visited, ap, parent, low, disc, nodes, Time=0):
        """
        A recursive function that find articulation points
        using DFS traversal
        :param u: the vertex to be visited next
        :param visited: keeps track of visited vertices
        :param ap: stores articulation points
        :param parent: stores parent vertices in DFS tree
        :param low: low value
        :param disc: stores discovery times of visited vertices
        :param nodes: current node set

        for more details: https://www.geeksforgeeks.org/articulation-points-or-cut-vertices-in-a-graph/
        """

        # Count of children in current node
        children = 0

        # Mark the current node as visited and print it
        visited[u] = True

        # Initialize discovery time and low value
        disc[u] = Time
        low[u] = Time
        Time += 1

        # for all the vertices adjacent to this vertex
        for v in self.G.vertex(u).out_neighbors():
            # If v is not visited yet, then make it a child of u
            # in DFS tree and recur for it
            if int(v) in nodes:
                if not visited[int(v)]:
                    parent[int(v)] = int(u)
                    children += 1
                    self.APUtil(int(v), visited, ap, parent, low, disc, nodes, Time)

                    # Check if the subtree rooted with v has a connection to
                    # one of the ancestors of u
                    low[u] = min(low[u], low[int(v)])

                    # u is an articulation point in following cases
                    # (1) u is root of DFS tree and has two or more children.
                    if parent[u] == -1 and children > 1:
                        ap[u] = True

                    # (2) If u is not root and low value of one of its child is more
                    # than discovery value of u.
                    if parent[u] != -1 and low[int(v)] >= disc[u]:
                        ap[u] = True

                        # Update low value of u for parent function calls
                elif int(v) != parent[u]:
                    low[u] = min(low[u], disc[int(v)])

                    # The function to do DFS traversal. It uses recursive APUtil()

    def is_AP(self, nodes):
        """
        Checks which nodes in the given set of nodes can NOT be removed without breaking
        disconnecting the induced subgraph
        :param nodes: set of nodes that make an induced subgraph of G
        :return: dictionary where each key is a node and each value indicates if a node is
        removable (articulation point)
        """
        visited = dict()
        disc = dict()
        low = dict()
        parent = dict()
        ap = dict()
        for node in nodes:
            visited[node] = False
            disc[node] = float("Inf")
            low[node] = float("Inf")
            parent[node] = -1
            ap[node] = False

        # Call the recursive helper function
        # to find articulation points
        # in DFS tree rooted with vertex 'i'
        for node in nodes:
            if not visited[node]:
                self.APUtil(node, visited, ap, parent, low, disc, set(nodes))

        return ap

    def score(self,nodes,labels,cluster_centers):
        """
        scores  given solution which is defined as a subnetwork and patient clusters
        :param nodes: list of nodes used in the solution
        :param labels: patient cluster labels
        :return: objective function value
        """
         #Compute distances to center candidates
        GE=self.ge.transpose()
        GE= GE.to_numpy().reshape(71,2,58051)
        #
        initialize centroids
        # function to update the centroids
        # define function called cdist metric which includes the 3 different metrices dtw,softdtw and euclidean
        if metric == "dtw":
                    def metric_fun(x, y):
                        return tslearn.metrics.cdist_dtw(x, y, n_jobs=self.n_jobs,
                                         verbose=self.verbose, **metric_params)
        elif self.metric == "softdtw":
                    def metric_fun(x, y):
                        return cdist_soft_dtw(x, y, **metric_params)
                    
        # intialization of cluster centers
        
        if cluster_centers is None:
            cluster_centers = centroid_init(GE,k, cdist_metric=metric_fun, random_state=rs)
        
        else:
        cluster_centers1 = dtw_barycenter_averaging( X=GE[labels== k],
                                                   barycenter_size=None,
                                                   init_barycenter=cluster_centers,
                                                   metric_params=None,
                                                   verbose=False)

        
      
         if nodes is None:
            Distance = tslearn.metrics.cdist_dtw(GE,cluster_centers1)
            inertia = Distance.min(axis=1).sum()
        else:
            Distance = tslearn.metrics.cdist_dtw(GE[:,:,nodes],cluster_centers1)
            #inertia = Distance.min(axis=1).sum()
            
           inertia= _compute_inertia(Distance, labels, squared = True)
        return inertia

       # vs = []
       # centroids = []
       # for i in range(self.k):
       #     idx = np.asarray(labels == i).nonzero()[0]
       #     vals = np.mean(self.ge[np.ix_(nodes, idx)], axis=1)
       #     centroids.append(np.mean(self.ge[np.ix_(nodes, idx)]))
       #     vs.append(vals)
       # objective = []
       #  for i in range(self.k):
       #     dif = np.mean(np.power((vs[i] - centroids[i]), 2))
        #    objective.append(dif)

       # return np.mean(objective)

    def dfs(self, node, d, visited=None):
        """
        Recursive DFS
        :param node: starting node
        :param d: length of s required subnetwork
        :param visited: should be left empty
        :return: a list of connected nodes of length d
        """

        if visited is None:
            visited = []
        if int(node) not in visited and len(visited) < d:
            visited.append(int(node))
            for neighbour in node.out_neighbors():
                self.dfs(neighbour, d, visited)
        if len(visited) == d:
            return visited

    def get_candidates(self, nodes):
        """
        Outputs first-degree neighbours of given nodes in graph G
        :param nodes: list of nodes that form a subnetwork/solution
        :return: list of first neighbour nodes
        """
        subst_candidates = flatten([[int(n) for n in self.G.get_all_neighbours(x)] for x in nodes])
        subst_candidates = set(subst_candidates).difference(set(nodes))
        return subst_candidates

    def insertion(self, nodes, labels):
        """
        Scores all possible insertions
        :param nodes: current solution
        :param labels: patient clusters labels
        :return: dictionary where key are possible insertions and values are scores
        """
        results = dict()
        size = len(nodes)
        if size < self.L_max:
            candidates = self.get_candidates(nodes)
            for c in candidates:
                nodes_new = nodes + [c]
                sc = self.score(nodes_new, labels)
                results["i_" + str(c)] = sc
        return results

    def deletion(self, nodes, labels, AP):
        """
        Scores all possible deletions
        :param nodes: current solution
        :param labels: patient clusters labels
        :param AP: articulation points (can't be removed since they separate the subnetwork)
        :return: dictionary where key are possible deletions and values are scores
        """
        results = dict()
        size = len(nodes)

        if size > self.L_min:
            for node in nodes:
                if not AP[node]:
                    nodes_new = list(set(nodes).difference({node}))
                    sc = self.score(nodes_new)
                    results["d_" + str(node)] = sc
        return results

    def subst(self, nodes, labels, AP):
        """
        Scores all possible substitutions
        :param nodes: current solution
        :param labels: patient clusters labels
        :param AP: articulation points (can't be removed since they separate the subnetwork)
        :return: dictionary where key are possible substitutions and values are scores
        """
        results = dict()
        size = len(nodes)
        if (size < self.L_max) and (size > self.L_min):
            for node in nodes:
                without_node = set(nodes) - {node}
                candidates = self.get_candidates(list(without_node))
                candidates = candidates - {node}
                for c in candidates:
                    if AP[node]:
                        nodes_new = list(without_node.union({c}))
                        if self.is_connected(nodes_new):
                            sc = self.score(nodes_new, labels)
                            results["s_" + str(node) + "_" + str(c)] = sc

                    else:
                        nodes_new = list(without_node.union({c}))
                        sc = self.score(nodes_new, labels)
                        results["s_" + str(node) + "_" + str(c)] = sc

        return results

    def is_connected(self, nodes):
        """
        Checks if a subgraph of G that consists of the given nodes is connected
        :param nodes: list of nodes
        :return: bool
        """
        sg = self.G.new_vertex_property("bool")
        for node in nodes:
            sg[node] = True
        g = gt.GraphView(self.G, vfilt=sg)

        comp, _ = gt.label_components(g, vprop=sg)
        if len(set(comp.a[nodes])) > 1:
            return False
        else:
            return True


    @staticmethod
    def do_action_nodes(action, nodes):
        """
        Updates the set of nodes given the action
        :param action: a key from the results dictionary that has a description of an action
        :param nodes: previous solution
        :return: new set of nodes
        """
        if len(action.split("_")) == 2:  # inserion or deletion
            act, node = action.split("_")
            node = int(node)
            if act == "i":
                nodes = nodes + [node]
            else:
                nodes = list(set(nodes).difference({node}))
        else:  # substitution
            act, node, cand = action.split("_")
            node = int(node)
            cand = int(cand)
            nodes = nodes + [cand]
            nodes = list(set(nodes).difference({node}))
        return nodes

    @staticmethod
    def to_key(nodes):
        """
        Creates a string representation of nodes
        :param nodes: node list
        :return: string of nodes
        """
        nodes = sorted(nodes)
        nodes = [str(node) for node in nodes]
        nodes = "|".join(nodes)
        return nodes

    @staticmethod
    def do_action_patients(action, labels):
        """
        Modifies patient cluster labels according to the given action
        :param action: a key from results dictionary
        :param labels: old cluster labels
        :return: updated cluster labels
        """
        if len(action.split("_")) == 2:  # add patient to a group
            _, group, idx = action.split("_")
            idx = int(idx)
            group = int(group)
            labels[idx] = group
        else:  # substitution
            _, idx1, idx2 = action.split("_")
            idx1 = int(idx1)
            idx2 = int(idx2)
            old = labels[idx1]
            labels[idx1] = labels[idx2]
            labels[idx2] = old
        return labels

    def ls_on_genes(self, nodes, labels, solutions, score0, T):
        """
        Runs local search on a gene set
        :param nodes: current node set
        :param labels: current patient clusters lables
        :param solutions: dictionary wth previously used solutions
        :param score0: last objective function score
        :param T: temperature for SA

        :return:
        nodes - new set of nodes
        score1 - new score
        move - True if further optimization was possible
        """
        # SUBNETWORK OPTIMIZATION
        move = False  # indicates if the solution feasible
        AP = self.is_AP(nodes)
        results = {**self.insertion(nodes, labels), **self.deletion(nodes, labels, AP),
                   **self.subst(nodes, labels, AP)}
        # first select the highest scoring solution which doesn't lead to the same set of nodes
        while not move:
            action = min(results.items(), key=operator.itemgetter(1))[0]
            score1 = results[action]
            # check if the solution is feasible
            nodes_new = self.do_action_nodes(action, nodes)
            nodes_new = self.to_key(nodes_new)
            if solutions.get(nodes_new) == None:  # solution wasn't used before
                move = True
            else:
                del results[action]
                if len(results) == 0:
                    print("no more feasible solutions")
                    return nodes, score0, move

        delta = score1 - score0
        if delta < 0:  # move on
            print(action)
            print("Score after genes LS {0}".format(score1))
            nodes = self.do_action_nodes(action, nodes)

        else:  # SA
            try:
                val = np.exp(-delta / T)
            except RuntimeError:
                val = 0
            p = np.random.uniform()
            if val > p:  # move on
                print("SA on genes at {0} degrees".format(T))
                print(action)
                print("Score after genes LS".format(score1))
                nodes = self.do_action_nodes(action, nodes)
            else:  # terminate if no improvement in two rounds
                print("too cold for genes SA, no actions taken")
                move = False
                score1 = score0
        return nodes, score1, move

    def ls_on_patients(self, nodes):
        """

        :param nodes: current node set

        :return:
        labels - new patient clusters labels
        score1 - updated objective function score
        """
        # PARTITION OPTIMIZATION
        GE=self.ge.transpose()
        GE= GE.to_numpy().reshape(71,2,58051)
        seed = 0
        labels = TimeSeriesKMeans(n_clusters=2,
                          n_init=2,
                          metric="dtw",
                          max_iter_barycenter=10,
                          random_state=seed).fit_predict(GE[:,:,nodes])
                          #random_state=seed).fit_predict(GE/self.ge[nodes, :].T)
        
        #kmeans = KMeans(n_clusters=self.k, random_state=0).fit(self.ge[nodes, :].T)
        #     else:
        #         centroids = []
        #         for i in range(k):
        #             idx = np.asarray(labels0 == i).nonzero()[0]
        #             vals = np.mean(ge[np.ix_(nodes, idx)], axis = 1)
        #             centroids.append(vals)
        # #            print(vals)
        #         kmeans = KMeans(n_clusters=k, random_state=0, init = np.array(centroids)).fit(ge[nodes, :].T)
        score = self.score(nodes)
        print("Reclustered score {}".format(score))
        return labels, score

    def run_ls(self):
        """
        Runs LS on patients and nodes
        :return:
        best_nodes - optimized node set
        best_labels - optimized label set
        score_max -maximal score
        """

        T0 = self.T
        T = T0
        score_min = 0
        best_nodes = []
        best_labels = []
        #n, m = self.ge.shape
        pats = len(self.patients)//2
        #pats= self.patients -n
        # initialization
        if self.seed is None:
            nodes = []
            # for  whatever reason dfs sometimes returns nothing
            no_type = True
            while no_type:
                nodes = self.dfs(self.G.vertex(np.random.choice(self.genes, 1)[0]), self.init_size) 
                if nodes is not None:
                    no_type = False
        else:
            nodes = self.seed
        
        labels = np.random.choice([0, 1], pats)
        score0 = self.score(nodes)
        start_score = score0
        labels = self.ls_on_patients(nodes)
        #if self.verbose:
            #print('start score is:{} and first set of labels is{}'.format(start_score,labels))

        scores = [start_score]
        solutions = dict()
        nodes_keys = self.to_key(nodes)
        solutions[nodes_keys] = ""
        count = 0
        for it in range(self.max_iter):
            if count == 0:
                labels,score1 = self.ls_on_patients(nodes)
            else:
                labels,score1= self.ls_on_patients(nodes)
            nodes_backup = nodes
            labels_backup = labels
            nodes, score2, move_genes = self.ls_on_genes(nodes, labels, solutions, score1, T)
            if not self.is_connected(nodes):
                print("something is wrong, network is disconnected")
                return nodes_backup, labels_backup, 0
            T = T0 * (0.9 ** it)  # cool down
            if self.verbose:
                print(it)
            solutions[self.to_key(nodes)] = ""
            scores.append(score2)
            #if self.plot:
             #   convergence_plot(scores)

            score0 = score2
            if score2 > score_min:
                score_min = score2
                best_nodes = nodes
                best_labels = labels
                print('iteration number{}, best nodes:{}, best_labels:{}, score_min:{}'.format(it,best_nodes,best_labels,score_min))
            else:
                print('iteration number{}, nodes:{}, labels:{}, score:{}'.format(it,nodes,labels,score2))
            count = count + 1
            if not move_genes:
                break
           
        return best_nodes, best_labels, score_min

