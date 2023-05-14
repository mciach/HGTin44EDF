from ete3 import Tree, NCBITaxa
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from collections import deque, Counter
from warnings import warn
ncbi = NCBITaxa()

class UnrootedForest():
    def __init__(self, T=None, V=None, E=None, edge_lengths=None, edge_supports=None, binary=True):
        """
        Creates a single-tree unrooted forest either from a rooted binary ete3 tree T
        or from sets of vertices and edges corresponding to an unrooted binary tree.
        T has to have unique node names, including internal nodes.
        V contains node names from T.
        Other node attributes are discarded unless supplied to edge_labels or node_labels.     
        An edge is of the form frozenset({i, j}), where i, j are nodes from V.
        edge_supports and edge_lengths are dicts frozenset({i, j}) : float.
        The flag binary (default True) ensures that the tree is binary.
        Setting it to false disables binary checks.  
        """
        self.binary = binary
        if T is not None:
            if binary:
                assert all(n.is_leaf() or len(n.children) == 2 for n in T.traverse()), 'T is non-binary'
            self.nodes = set()  # required for initialization
            self.edges = set()
            self.edge_supports = dict()
            self.edge_lengths = dict()
            self._add_ete3_tree(T)
            if len(T.children) == 2:  # removing root node
                self.reduce_node(T.name)
            if binary:
                for v in self.nodes:
                    assert self.get_degree(v) in {1, 3}, 'Improper tree structure: a node has degree other than 1 or 3'
            else:
                for v in self.nodes:
                    assert self.get_degree(v) >= 1, 'Disjoint tree encountered'
        elif V is not None and E is not None:
            self.nodes, self.edges = V, E
            if binary:
                for v in self.nodes:
                    assert self.get_degree(v) in {0, 1, 3}, 'Improper tree structure: a node has degree other than 1 or 3'
                assert len(E) == len(V) - 1, 'A binary unrooted tree has to have |E| = |V|-1'
            self.edge_supports = edge_supports
            self.edge_lengths = edge_lengths

    def _convert_to_graph(self, T):
        """
        Returns a tuple of nodes and edges (V, E) representing an ete3 tree T.
        """
        E = set()
        V = set()
        for n in T.traverse():
            assert n.name and n.name not in V, "T doesn't have unique node names"
            V.add(n.name)
            if not n.is_root():
                E.add(frozenset([n.name, n.up.name]))
        return(V, E)

    def _add_ete3_tree(self, T):
        """
        Adds a new unrooted tree to the forest, generated from an ete3 tree T. 
        """
        for n in T.traverse():
            if not n.name:
                raise ValueError('Unnamed node')
            if n.name in self.nodes:
                raise ValueError("Node %s already in the forest" % n.name)
            self.nodes.add(n.name)
            if not n.is_root():
                edge = frozenset([n.name, n.up.name])
                self.edges.add(edge)
                self.edge_lengths[edge] = n.dist
                self.edge_supports[edge] = n.support
                if n.is_leaf():
                    assert n.support in {1, 100}, 'Unsupported leaf - impossible, probably wrong tree structure or support parsing'

    def get_connected_components(self):
        visited = set()
        cc = []   # each entry = (V, E) of a component
        def get_component(node, cmp):
            cmp[0].add(node)
            nbve = self.get_adjacent_edges(node)
            for e in nbve:
                if e not in cmp[1]:
                    cmp[1].add(e)
                    n1, n2 = e
                    if n1 == node:
                        get_component(n2, cmp)
                    elif n2 == node:
                        get_component(n1, cmp)
        
        for v in self.nodes:
            if v not in visited:
                cmp = (set(), set())             
                get_component(v, cmp)
                cc.append(cmp)
                visited.update(cmp[0])
        return cc
        
    def _convert_component_to_ete3(self, cmp):
        """
        Returns an ete3 Tree corresponding to an unrooted tree 
        with vertices V == cmp[0] and edges E == cmp[1].  
        The unrooted tree needs to be connected, with every internal
        node of degree 3.  
        """
        if self.binary:
            assert len(cmp[1]) == len(cmp[0]) - 1, 'A binary unrooted tree has to have |E| = |V|-1'
        if len(cmp[0]) == 1:
            return Tree(name=next(iter(cmp[0])), dist=0)
        else: 
            UG = UnrootedForest(V = cmp[0], E = cmp[1])
            ref_v = next(v for v in UG.nodes if UG.get_degree(v) == 1)
            ref_edge = UG.get_adjacent_edges(ref_v)
            assert len(ref_edge) == 1
            ref_edge = ref_edge.pop()
            ref_dist = self.edge_lengths[ref_edge]
            a, b = ref_edge
            assert ref_v in {a, b}
            v_neighb = a if ref_v == b else b  
            T = Tree(dist=0)
            ref_v_in_T = T.add_child(name = ref_v, dist=0.5*ref_dist)
            v_neighb_in_T = T.add_child(name = v_neighb, dist=0.5*ref_dist)
        def generate_tree(T, v, vp):
            """
            Generates tree through DFS from v directed from already visited vp
            """
            nbv = UG.get_neighbours(v)
            for neighbour in nbv:
                if neighbour != vp:
                    distance = self.edge_lengths[frozenset({neighbour, v})]
                    new_node = T.add_child(name = neighbour, dist=distance)
                    generate_tree(new_node, neighbour, v)
        generate_tree(v_neighb_in_T, v_neighb, ref_v) 
        return(T) 
    
    def get_ete3(self):
        """
        Returns the forest as a list of arbitrarily rooted trees.  
        """
        forest = []
        cmps = self.get_connected_components()
        for cmp in cmps:
            T = self._convert_component_to_ete3(cmp)
            forest.append(T)
        return forest
        
    def get_neighbours(self, v):
        assert v in self.nodes, 'Node not found'
        nbv = set()
        for e in self.edges:
            x, y = e
            if v == x:
                nbv.add(y)
            elif v == y:
                nbv.add(x)
        return nbv
    
    def get_adjacent_edges(self, v):
        """
        Returns edges adjacent to the node v. 
        """
        nbve = []
        for e in self.edges:
            if v in e:
                nbve.append(e)
        nbve = set(nbve)
        return nbve
        
    def get_degree(self, v):
        return len(self.get_adjacent_edges(v))
    
    def get_leaves(self):
        return [v for v in self.nodes if self.get_degree(v) == 1]
    
    def reduce_node(self, v):
        """
        Removes a node v and joins its neighbours.  
        The node v needs to have degree 2 or be separated (degree 0).   
        Returns the new edge connecting v's neighbours or None.  
        """
        nbve = self.get_adjacent_edges(v)        
        if len(nbve) == 0:
            return None
        elif len(nbve) == 2:
            new_edge_dist = sum(self.edge_lengths[e] for e in nbve)
            new_edge_support = max(self.edge_supports[e] for e in nbve)
            self.nodes.remove(v)  # This is costly... 
            self.edges -= nbve 
            new_edge = frozenset(nv for e in nbve for nv in e if nv != v)
            assert len(new_edge) == 2
            self.edges.add(new_edge)
            # we won't remove labels of non-existent edges for now
            # because those dicts should be rather small anyways
            self.edge_lengths[new_edge] = new_edge_dist  
            self.edge_supports[new_edge] = new_edge_support
            return new_edge
        else:
            raise ValueError('Can only reduce degree-2 or degree-0 nodes')

        
    def disintegrate(self, threshold):
        """
        Cut edges of the forest until all edges have length
        shorter than the threshold.  
        """
        long_edge = next((e for e in self.edges if self.edge_lengths[e] >= threshold), None)
        while long_edge is not None:
            #print(long_edge)
            v, w = long_edge
            #print('Removing', long_edge)
            self.edges.remove(long_edge)  # This is costly... 
            #print('Reducing', v)
            ne1 = self.reduce_node(v)
            #print('Reducing', w)
            ne2 = self.reduce_node(w)  
            long_edge = next((e for e in self.edges if self.edge_lengths[e] >= threshold), None)
            
            
    def prune_long_leaves(self, threshold):
        """
        Cut leaves of the forest until all edges have length
        shorter than the threshold.  
        """
        leaves = set(self.get_leaves())
        long_edge = next((e for e in self.edges if self.edge_lengths[e] >= threshold and e & leaves), None)
        while long_edge is not None:
            #print(long_edge)
            v, w = long_edge
            #print('Removing', long_edge)
            self.edges.remove(long_edge)  # This is costly... 
            #print('Reducing', v)
            ne1 = self.reduce_node(v)
            #print('Reducing', w)
            ne2 = self.reduce_node(w)  
            leaves.difference_update(long_edge)
            long_edge = next((e for e in self.edges if self.edge_lengths[e] >= threshold and e & leaves), None)


    def get_clade_attaching_edges(self, nodenames):
        """
        Return all edges attaching clades induced by the given set of nodes.
        A clade is a maximum (w.r.t. containment) subtree that contains only
        nodes from the supplied set.
        """
        assert nodenames.issubset(self.nodes), 'Nodes not found'
        all_in_nodenames = {} # dict (node1, node2) : bool whether a directed edge node1 -> node2 points to an all-target subtree
        leaves = self.get_leaves()
        if set(nodenames) == set(leaves):
            return []
        for leaf in leaves:
            nbh = self.get_neighbours(leaf).pop()
            if leaf in nodenames:
                all_in_nodenames[(nbh, leaf)] = True
            else:
                all_in_nodenames[(nbh, leaf)] = False
        def dfs(n1, n2):
            """
            DFS from node n2 ignoring node n1.
            Checks if edge (n1, n2) points to an all-target clade.
            """
            neighbors = [p for p in self.get_neighbours(n2) if p != n1]
            for v in neighbors:
                if (n2, v) not in all_in_nodenames:
                    dfs(n2, v)
            all_in_nodenames[(n1, n2)] = all(all_in_nodenames[(n2, v)] for v in neighbors)
        for leaf in leaves:
            nbh = self.get_neighbours(leaf).pop()
            dfs(leaf, nbh)
        clade_attaching_edges = []  # List of directed edges containing the roots of clades
        for v in self.nodes:
            neighbors = list(self.get_neighbours(v))
            separating_directions = [all_in_nodenames[(v, nb)] or all_in_nodenames[(nb, v)]  for nb in neighbors]
            # Simple version, works only for binary trees:
##            if sum(separating_directions) == 1:
##                clade_attaching_edges.append((v, neighbors[separating_directions.index(True)]))
            # Alternative version, better suited for non-binary trees,
            # but requires symmetrization of separating_directions:
            if any(separating_directions) and not all(separating_directions):
                for nbh_v, sep in zip(neighbors, separating_directions):
                    if sep:
                        clade_attaching_edges.append((v, nbh_v))
        return clade_attaching_edges


    def get_clade_leaves(self, attaching_edge):
        """
        Given a directed edge (w, v), returns leaves reachable from w through v, i.e. leaves of a clade induced by this edge.
        Note that attaching_edge is a tuple of nodes, not an edge as implemented internally in UnrootedForest class.  
        """
        clade_leaves = []
        def dfs_leaf_search(leaf_list, w, v):
            """
            DFS search from v ignoring the direction of w.
            The edge {w, v} is assumed to be in the tree.   
            """
            nbh_of_v = self.get_neighbours(v)
            if len(nbh_of_v) == 1:
                leaf_list.append(v)
            else:
                for neighbour in nbh_of_v:
                    if neighbour != w:
                        dfs_leaf_search(leaf_list, v, neighbour)
        dfs_leaf_search(clade_leaves, *attaching_edge)
        return set(clade_leaves)
    

    @staticmethod
    def get_cdi(analyzed_split, reference_multisplit, verbose=False):
        """
        Returns the value of the clade displacement index based on
        a split induced in the analyzed tree (partitioning of leaf labels into two sets)
        and a multisplit induced in the reference tree (partitioning of leaf
        labels into multiple sets).
        """
        assert len(analyzed_split) == 2, 'The analyzed_split needs to have lenght 2, corresponding to a binary analyzed tree.'
        analyzed_labels = set(l for subset in analyzed_split for l in subset)
        reference_labels = set(l for subset in reference_multisplit for l in subset)
        assert analyzed_labels == reference_labels, 'Inconsistent sets of labels in analyzed and reference splits'
        reference_subtree_sizes = map(len, reference_multisplit)
        left_weights = Counter(analyzed_split[0])
        label_counts = left_weights + Counter(analyzed_split[1])
        if verbose:
            print('label counts:', label_counts)
        left_weights = {label: left_weights[label]/label_counts[label] for label in label_counts}
        weight_sum = sum(left_weights[label] for label in left_weights)
        if verbose:
            print('label left-weights:', left_weights)
        reference_subtree_weights = list(sum(left_weights[label] for label in subset) for subset in reference_multisplit)
        if verbose:
            print('reference subtree left-weights:')
            print(list(reference_subtree_weights))
        cdi_partial_costs = [clade_size - 2*clade_weight for clade_size, clade_weight in zip(reference_subtree_sizes, reference_subtree_weights)]
        if verbose:
            print('subtree costs:')
            print(cdi_partial_costs)
        mincost = min(cdi_partial_costs)
        maxcost = max(cdi_partial_costs)
        mincost_index = cdi_partial_costs.index(mincost)
        maxcost_index = cdi_partial_costs.index(maxcost)
        cdi = weight_sum + mincost + sum(cost for i, cost in enumerate(cdi_partial_costs) if cost < 0 and i != mincost_index and i != maxcost_index)
        return cdi

    
    def get_clade_displacements(self, target_leaves, reference_tree, leaf_name_map, verbose=False):
        """
        Measures displacement of clades induced by target_nodes with respect to
        reference_unrooted_tree.
        The node_name_map gives the mapping between leaves of self and reference.
        The mapping may not be injective (self can have multiple leaves corresponding to the same leaf of the reference).  
        Returns a list of tuples (clade, CDI), where CDI is the clade displacement index
        of the clade in self with respect to the reference tree.
        Each clade is evaluated independently (other target_nodes are masked in both trees).
        The target_leaves need to induce a single clade in the reference tree after mapping with node_name_map.
        However, the induced clade may have multiple neighbours (the reference tree may not be binary).  
        """
        assert self.binary, 'The analyzed tree needs to be binary.'
        clade_cdis = []
        # Note: when looking for the corresponding clade in the reference tree, we want to mask
        # the target_leaves which are not in clade_leaves.
        # This can potentially change the neighbouring clades if some of the masked leaves are a local outgroup of clade_leaves
        # in the reference tree if we analyzed the neighbouring clades of clade induced by clade_leaves.
        # This is why we analyze the composition of the reference_split computed only once using all target_leaves.
        # This is equivalent to getting a split after masking nodes if target_leaves are monophyletic in the reference.
        mapped_target_leaves = set(leaf_name_map[name] for name in target_leaves)
        reference_clade_edge = reference_tree.get_clade_attaching_edges(mapped_target_leaves)
        assert len(reference_clade_edge) == 1, 'Target leaves do not form a single clade in the reference tree.'
        ref_anode, ref_rnode = reference_clade_edge[0]  # attaching and root nodes
        ref_neighbours = reference_tree.get_neighbours(ref_anode)
        reference_multisplit = [reference_tree.get_clade_leaves((ref_anode, neighbour)) for neighbour in ref_neighbours if neighbour != ref_rnode]
        if verbose:
            print('Split in the reference tree induced by clade', target_leaves, ':')
            print(reference_multisplit)
        clade_attaching_edges = self.get_clade_attaching_edges(target_leaves)
        for clade_edge in clade_attaching_edges:
            if verbose:
                print('Clade', clade_edge)
            clade_anode, clade_rnode = clade_edge  # attaching and root nodes 
            clade_leaves = self.get_clade_leaves(clade_edge)
            clade_neighbour_nodes = self.get_neighbours(clade_anode)
            clade_split = [self.get_clade_leaves((clade_anode, neighbour)) for neighbour in clade_neighbour_nodes if neighbour != clade_rnode]
            if verbose:
                print('Split in the analyzed tree induced by edge', clade_edge, ':')
                print(clade_split)
            clade_split = [[label for label in subset if label not in target_leaves] for subset in clade_split]
            if verbose:
                print('Split in the analyzed tree after masking target leaves:')
                print(clade_split)
            mapped_clade_split = [[leaf_name_map[name] for name in clade] for clade in clade_split]
            if verbose:
                print('Split mapped to names of the reference:')
                print(mapped_clade_split)
            cdi = UnrootedForest.get_cdi(mapped_clade_split, reference_multisplit)
            if verbose:
                print('CDI:', cdi)
            clade_cdis.append((clade_edge, cdi))
        return clade_cdis


    def get_optimal_clade_location(self, target_leaves, reference_tree, leaf_name_map, verbose=False):
        """
        Looks for an optimal edge to which a clade composed of target_leaves
        can be attached so that its location is as similar to its location in reference_tree
        as possible.
        target_leaves need to induce a single monophyletic clade in the reference_tree.
        Returns a pair (optimal edge, optimal cdi).
        """
        assert self.binary, 'The analyzed tree needs to be binary.'
        mapped_target_leaves = set(leaf_name_map[name] for name in target_leaves)
        reference_clade_edge = reference_tree.get_clade_attaching_edges(mapped_target_leaves)
        assert len(reference_clade_edge) == 1, 'Target leaves do not form a single clade in the reference tree.'
        ref_anode, ref_rnode = reference_clade_edge[0]  # attaching and root nodes
        ref_neighbours = reference_tree.get_neighbours(ref_anode)
        reference_multisplit = [reference_tree.get_clade_leaves((ref_anode, neighbour)) for neighbour in ref_neighbours if neighbour != ref_rnode]
        if verbose:
            print('Split in the reference tree:')
            print(reference_multisplit)
        cdi_values = []
        for i, edge in enumerate(self.edges):
            if verbose:
                print('Checking edge', edge)
            edge = list(edge)
            subset1 = self.get_clade_leaves((edge[0], edge[1]))
            subset1 -= target_leaves
            if not subset1:
                if verbose:
                    print('CDI = N/A (within target clade); skipping')
                continue
            subset1 = [leaf_name_map[name] for name in subset1]
            subset2 = self.get_clade_leaves((edge[1], edge[0]))
            subset2 -= target_leaves
            if not subset2:
                if verbose:
                    print('CDI = N/A (within target clade); skipping')
                continue
            subset2 = [leaf_name_map[name] for name in subset2]
            mapped_split = [subset1, subset2]
            cdi = UnrootedForest.get_cdi(mapped_split, reference_multisplit)
            if verbose:
                print('CDI = ', cdi)
            if cdi == 0:
                return (edge, 0)  # minimum value of cdi
            else:
                cdi_values.append((frozenset(edge), cdi))
        return min(cdi_values, key = lambda x: x[1])
            
                    
                
        

            
###
###  CLADE DISPLACEMENT FUNCTIONS
###


def update_taxids(taxids):
    """
    Given a list of (possibly obsolete or merged) taxIDs, 
    returns a dictionary that maps them to their current versions.
    """
    taxids = list(map(int, taxids))
    names = ncbi.get_taxid_translator(taxids)
    new_taxids = ncbi.get_name_translator([names[tx] for tx in taxids])
    assert all(len(new_taxids[n]) == 1 for n in new_taxids)
    return {tx: new_taxids[names[tx]][0] for tx in taxids}


def get_nodes_at_depth(T, depth=0):
    """
    Returns a list of nodes with a given distance from the root.
    """
    def go_down(node, curr_depth):
        if curr_depth == depth or node.is_leaf():
            return [node]
        elif curr_depth > depth:
            raise RuntimeError
        else:
            return [n for child in node.children for n in go_down(child, curr_depth + 1)]
    return go_down(T, 0)


def relabel_nodes(G, name_mapping):
    """
    Labels nodes of G using a dictionary name_mapping.   
    Returns a dictionary with inverse mapping such that G == relabel_nodes(G, relabel_nodes(G, name_mapping)).  
    Only the nodes with names in name_mapping get relabeled.  
    """
    reverse_map = dict()
    for n in G.traverse():
        if n.name in name_mapping:
            new_name = name_mapping[n.name]
            reverse_map[new_name] = n.name
            n.name = new_name
    return reverse_map


def add_scientific_names(G):
    """
    Adds a sci_name attribute to leaves of G.
    Leaves of G need to be labelled with taxIDs.
    """
    taxids = {l.name for l in G}
    translator = ncbi.get_taxid_translator(taxids)
    for l in G:
        l.sci_name = translator[int(l.name)]
        
def remove_intermediate_nodes(S):
    """
    Removes all internal nodes with a signle child.
    Works in situ.
    """
    to_remove = [s for s in S.traverse() if len(s.children) == 1 and not s.is_leaf()]
    for s in to_remove:
        s.delete()
    while len(S.children)==1:
        S.children[0].delete()
        

def get_children_taxids(S, parent_taxid):
    """
    Returns taxids of species from S that are descendants of parent_taxid in NCBI taxonomy.
    Leaves of S need to be named with taxIDs.
    S needs to contain intermediate nodes.
    """
    try:
        LCA = S.search_nodes(name=parent_taxid)[0]
    except IndexError:
        warn('Parent taxid %s not found in the tree!' % str(parent_taxid))
        return set()
    return set(l.name for l in LCA)


def unify_trees(G, S, quiet=False):
    """
    Trims both trees so that only the leaves with labels shared between the trees are left.
    Works in situ.
    """
    g_labels = set(l.name for l in G)
    s_labels = set(l.name for l in S)
    common_labels = g_labels & s_labels  # set union
    g_to_remove = [l for l in G if l.name not in common_labels]
    s_to_remove = [l for l in S if l.name not in common_labels]
#     print(len(g_to_remove), 'to remove from G')
#     print(g_to_remove)
#     print(len(s_to_remove), 'to remove from S')
#     print(s_to_remove)
    if not quiet:
        if g_to_remove:
            print(len(g_to_remove), 'to remove from G')
            print(g_to_remove)
        if s_to_remove:
            print(len(s_to_remove), 'to remove from S')
            print(s_to_remove)
            raise RuntimeError('X!')
    if len(g_to_remove) == len(G):
        warn('Tree G is deleted - no leaf matches S!')
    if len(s_to_remove) == len(S):
        warn('Tree S is deleted - no leaf matches G!')
    for n in g_to_remove:
        while not n.is_root():
            if len(n.up.children) > 1:
                n.detach()
                # print('detaching in G', n)
                if len(n) > 1:
                    raise RuntimeError('large in G')
                break
            n = n.up
    for n in s_to_remove:
        while not n.is_root():
            if len(n.up.children) > 1:
                n.detach()
                # print('detaching in S', n)
                if len(n) > 1:
                    raise RuntimeError('large in S')
                break
            n = n.up
            

def LCA(T, leaf_set):
    """
    Implemented myself because the ete3 version sux cox
    """
    # Since we use leaf pointers not leaf labels in this implementation, 
    # we can use the fact that leaves are unique
    # and traverse the tree in postorder, storing how many leaves we're lacking to 
    # be the common ancestor.
    n = len(leaf_set)
    if n == 1:
        return list(leaf_set)[0]
    assert leaf_set.issubset(set(T)), 'Target leaf set is not a subset of leaves of T!'
    number_below = {g: 0 for g in T.traverse()}
    prev_node = T
    for g in T.traverse('postorder'):
        if g.is_leaf():
            if g in leaf_set:
                number_below[g] = 1
        else:
            nb = sum(number_below[c] for c in g.children) 
            if nb == n:
                return g
            else:
                number_below[g] = nb   
    if nb == n:
        return T
    else:
        raise RuntimeError('Unable to find LCA!')
        
        
        
def check_strict_monophyly(T, target_names):
    """
    Returns True if target_names correspond to a strictly monophyletic clade in T.
    Strict monophyly is defined as the existence of a branch such that on one of its sides
    all the leaves belong to target_names, and no leaf from the other side belongs to it.
    If the tree contains only target or only non-target leaves, False is returned as well.
    The tree T may be non-binary and unrooted.
    It is assumed that T does not contain intermediate nodes. 
    
    Parameters
    ----------
    T: ete3 Tree object
    target_names: a set of species names. 
        It is assumed that leaves of T are labelled with species names, without gene identifiers. 
    """
    target_names = set(target_names)
    target_leaves = {l for l in T if l.name in target_names}
    other_leaves = {l for l in T if l.name not in target_names}
    if not target_leaves or not other_leaves:
        return False
    other_names = {l.name for l in other_leaves}
    
    target_lca = LCA(T, target_leaves)
    other_lca = LCA(T, other_leaves)
    names_under_target_lca = {l.name for l in target_lca}
    names_under_other_lca = {l.name for l in other_lca}
    # print(target_names, other_names, names_under_target_lca, names_under_other_lca)
    return target_names == names_under_target_lca or other_names == names_under_other_lca

def check_conditional_monophyly(T, target_names):
    """
    Returns True if there exists a binarization of T containing a strictly monophyletic clade 
    induced by the set of leaf labels target_names. 
    If there is no leaf of T with label from target_names, or if all leaves are labeled with
    target_names, return False.
    The tree may be unrooted.  
    """
    target_names = set(target_names)
    target_leaves = {l for l in T if l.name in target_names}
    other_leaves = {l for l in T if l.name not in target_names}
    if not target_leaves or not other_leaves:
        return False
    other_names = {l.name for l in other_leaves}
    
    target_lca = LCA(T, target_leaves)
    other_lca = LCA(T, other_leaves)
    if target_lca == other_lca:  # mixing of leaves
        return False
    
    # We'll check which of the target_lca and other_lca is below the root of T.
    # The one below the root will be selected as the node to binarize.
    # We iterate over its children and check if they can be partitioned into those
    # containing only target_leaves or only other_leaves.
    
    if target_lca == T:  # other_lca is below the root
        for c in other_lca.children:
            if not (set(c).isdisjoint(target_leaves) or set(c).isdisjoint(other_leaves)):
                return False
    else:
        for c in target_lca.children:
            if not (set(c).isdisjoint(target_leaves) or set(c).isdisjoint(other_leaves)):
                return False
    return True



def get_separating_branch(T, target_names):
    """
    If target_names are strictly monophyletic in T, return the node adjacent to the branch 
    that separates them.  
    Otherwise, return None.
    The separating branch is always given by result.up, regardless of the configuration of
    the tree (i.e. whether target_names' or the outgroup are internal in T, i.e. their LCA is an internal node)
    
    Parameters
    ----------
    T: ete3 Tree object
    target_names: a set of species names. 
        It is assumed that leaves of T are labelled with species names, without gene identifiers. 
    """
    if len(T.children) == 3:
        warn('The tree appears to be unrooted!')
    target_leaves = {l for l in T if l.name in target_names}
    other_leaves = {l for l in T if l.name not in target_names}
    if not target_leaves or not other_leaves:
        return None
    target_names = set(target_names)
    other_names = {l.name for l in other_leaves}
    
    target_lca = LCA(T, target_leaves)
    other_lca = LCA(T, other_leaves)
    names_under_target_lca = {l.name for l in target_lca}
    names_under_other_lca = {l.name for l in other_lca}
    # print(target_names, other_names, names_under_target_lca, names_under_other_lca)
    if target_names == names_under_target_lca:
        return target_lca
    elif other_names == names_under_other_lca:
        return other_lca
    else:
        return None
    
    
    
    
def get_multisplit_from_target(T, target_names, verbose=False):
    """
    Returns a list of clades going outwards from the LCA of target_names. 
    That is, a forest obtained by removing the clade of target_names together with its parent node.
    Effectively, this is an implementation of get_split_from_target for non-binary trees.
    Assumes that T is rooted. Applying this function to an unrooted tree will give improper results.
    The position of the root does not matter, so the tree may be rooted arbitrarily.
    
    Parameters
    ---------
    T: ete3 Tree object
        The tree to process.
        It is assumed that leaves of T are named using species names.
    target_names: iterable
        The set of species names defining the reference clade.
        The target names need to form a strictly monophyletic clade in T. 
        Otherwise, an error is raised.
    
    Returns: list
        A list of lists. Each component contains leaf labels of T corresponding
        to a clade going outwards from the target clade.
    
    """
    assert check_strict_monophyly(T, target_names), 'The target leaves are not monophyletic!'
    if len(T.children) == 3:
        warn('The tree appears to be unrooted!')
    target_leaves = {l for l in T if l.name in target_names}
    non_target_leaves = {l for l in T if l.name not in target_names}
    assert target_leaves, 'No target found in tree!'
    assert len(non_target_leaves)>1, 'Not enough non-target species in tree!'
    # print(target_leaves)
    lca = LCA(T, target_leaves)
    if verbose:
        print(lca)
    ### We need to handle three cases here: target LCA being the root, one of its children, or another node.
    # Case 1:
    if lca == T:
        # This is the most complex case - need to find LCA of non-target species and return its children
        if verbose:
            print('Case 1')
        lca2 = LCA(T, non_target_leaves)
        split = [[l.name for l in c] for c in lca2.children]
    elif lca in T.children and len(T.children) > 2:
        if verbose:
            print('case2')
        split = [[l.name for l in c] for c in T.children if c != lca]
    elif lca in T.children and len(T.children) == 2:
        if verbose:
            print('case3')
        if lca == T.children[0]:
            split = [[l.name for l in c] for c in T.children[1].children]
        else: 
            split = [[l.name for l in c] for c in T.children[0].children]
    else:
        if verbose:
            print('case4')
        split = [[l.name for l in c] for c in lca.get_sisters()]  # sisters of lca
        split.append([l.name for l in T if l not in lca.up])  # remaining leaves
    return split


def multisplit_distance(gene_split, species_split, verbose=False):
    """
    Computes the multisplit distance, defined as the minimal absolute difference 
    of the proportion of each label in the left out-going clade in both trees 
    over all refinements of species_split.
    
    Parameters
    ----------
    split1, split2: iterables
        Compared splits. Each component contains a list of leaf labels.
        The labels in split2 need to be unique, but in split1 they do not. 
        Labels are matched between splits by their name.
        All labels from gene_split and species_split need to correspond to each other exactly.
    """
    assert len(gene_split) == 2, 'gene_split needs to have length 2'
    gene_split_labels = set(l for s in gene_split for l in s)
    species_split_labels = set(l for s in species_split for l in s)
    if not gene_split_labels == species_split_labels:
        raise ValueError('the splits contain different labels:\n%s\n%s' % (str(sorted(gene_split_labels)), 
                                                                           str(sorted(species_split_labels))))
    
    # check if labels in split2 are unique:
    s2 = Counter()
    for cmp in species_split:
        for spc in cmp:
            s2[spc] += 1
    assert set(s2.values()) == {1}, 'non-unique labels in split2!'
    
    # compute label scores for gene_split
    # i.e. the proportion of occurence in the reference split
    left_count = Counter(gene_split[0])  
    total_count = left_count + Counter(gene_split[1])
    gene_label_score = {label: left_count[label]/total_count[label] for label in gene_split_labels}
    if verbose:
        print(gene_label_score)

    distance = 0.
    for cmp in species_split:
        if verbose:
            print(cmp)
        s1 = sum(gene_label_score[label] for label in cmp)
        s2 = len(cmp) - s1
        distance += min(s1, s2)
        if verbose:
            print('current choices:', s1, s2)
            print('current distance:', distance)
    assert distance >= 0, 'negative distance - bug!'
    return distance


## not needed anymore
# def split_leaves(G):
#     """
#     Splits leaves with multiple identifiers, corresponding to identical sequences, into polytomies.
#     """
#     # pattern = '\.[0-9]_[A-Z]'
#     for l in G:
#         splits = re.findall('([A-Z].*?\.[0-9])', l.name)
#         if len(splits) > 1:
#             l.name = '' 
#             for s in splits:
#                 l.add_child(name=s, dist=0)

    
### 
###  BRANCH LENGTH FUNCTIONS
###
            
def identity_to_sps(idv, niter=100):
    """
    Return the branch length (number of substitutions per site)
    corresponding to a given sequence identity between two sequences. 
    Based on Nick V. Grishin (1995). 
    """
    d = -np.log(idv)
    for i in range(niter):
        d = (1-np.exp(-2*d))/(2*idv)
    return d

def sps_to_identity(sps):
    """
    Returns the identity fraction corresponding to a given branch length,
    measured in average substitution per site.
    """
    return (1-np.exp(-2*sps))/(2*sps) 

def delta_variance(sps, idv, n=300):
    return (sps**2)*idv*(1-idv)/(n*(1-idv*(2*sps+1))**2)
    
def idv_std(sps, n=300):
    idv = sps_to_identity(sps)
    return delta_variance(sps, idv, n)**0.5


###
### ALIGNMENT PROCESSING
### 

def process_alignment(aln,  
                      remove_constant_columns=True,
                      max_column_gap_proportion=0.8, 
                      min_row_nb_of_aa=10, 
                      max_row_gap_proportion=0.9):
    """
    Removes columns with gap proportion above threshold and (optionally) only a single amino acid (constant columns).
    Next, discards sequences with less than a given number of non-gap and non-X residues, 
    or with gap or X character proportion above threshold. This is calculated AFTER discarding gappy and constant columns!
    """
    ncol_raw = aln.get_alignment_length()
    nrow_raw = len(aln) 
    is_constant = [False]*ncol_raw  # will stay False if remove_constant_columns == False
    is_gappy = [False]*ncol_raw
    for col_id in range(ncol_raw):
        gap_nb = sum(char == '-' for char in aln[:,col_id])
        if remove_constant_columns:
            unique_aas = set(char.upper() for char in aln[:,col_id] if char != '-' and char != 'x' and char != 'X')
            if len(unique_aas) <= 1:
                is_constant[col_id] = True
        if gap_nb/nrow_raw > max_column_gap_proportion:
            is_gappy[col_id] = True
    
    ncol_processed = ncol_raw - sum(c or g for c,g in zip(is_constant, is_gappy))
    all_sequences = []
    all_IDs = []
    for seq in aln:
        trimmed_sequence = ''.join(char for char, constant, gappy in zip(seq.seq, is_constant, is_gappy) if not constant and not gappy)
        # assert len(trimmed_sequence) == ncol_processed
        # nb_of_gaps = sum(char == '-' for char in trimmed_sequence)
        nb_of_aa = sum(char != '-' and char.lower() != 'x' for char in trimmed_sequence)
        # if nb_of_gaps < gap_proportion_threshold*len(trimmed_sequence):
        if nb_of_aa >= min_row_nb_of_aa and nb_of_aa/ncol_processed >= 1 - max_row_gap_proportion:
            all_sequences.append(trimmed_sequence)
            all_IDs.append(seq.id)
# # In order to remove duplicates: 
#             try:
#                 seq_reference = all_sequences.index(trimmed_sequence)
#             except ValueError:
#                 all_sequences.append(trimmed_sequence)
#                 all_IDs.append(seq.id)  
#             else:
#                 all_IDs[seq_reference] = all_IDs[seq_reference] + ',' + seq.id
    new_aln = MultipleSeqAlignment([SeqRecord(Seq(seq), seqid, '', '') for seq, seqid in zip(all_sequences, all_IDs)])         
    return new_aln


def get_aa_counts(aln):
    cts = []
    for seq in aln:
        cts.append(sum(aa!='-' for aa in seq))
    return cts

def discard_constant_columns(aln):
    is_constant = [False]*aln.get_alignment_length()
    for col_id in range(aln.get_alignment_length()):
        unique_aa = set(aln[:,col_id])
        unique_aa.discard('-')
        if len(unique_aa) == 1:
            is_constant[col_id] = True
    for seq in aln:
        seq.seq = ''.join(char for char, constant in zip(seq.seq, is_constant) if not constant)


###
### HGT ANALYSIS
### 


        
if __name__ == '__main__':
    ### Unrooted forest testing
    G = Tree("(((A, B)AB, (C, D)CD:2)ABCD:2, ((E:2, F)EF, G)GEF)R;", format=1)
    print(G.get_ascii(attributes=['name', 'dist']))
    # Mapping node names
    namemap = {n.name: i for i,n in enumerate(G.traverse('preorder')) if not n.is_root()}
    invmap = relabel_nodes(G, namemap)
    print(invmap)
    print(G.get_ascii(attributes=['name', 'dist']))
    relabel_nodes(G, invmap)
    print(G.get_ascii(attributes=['name', 'dist']))
    print('Pruning long leaves:')
    UG = UnrootedForest(G)
    # UG.nodes
    UG.prune_long_leaves(2)
    forest = UG.get_ete3()  
    len(forest)
    for T in forest:
        print(T.get_ascii(attributes=['name', 'dist']))

    print('Cutting long branches:')
    UG = UnrootedForest(G)
    # UG.nodes
    UG.disintegrate(2)
    forest = UG.get_ete3()  
    len(forest)
    for T in forest:
        print(T.get_ascii(attributes=['name', 'dist']))

    
    print('Detecting clades in a single tree:')
    G = Tree('(F_3, (((B11_1, B12_1)B1_1, (B21_1, B22_1)B2_1)1, (((B11_2, F_2)B1F, B12_2)B1_2, (B21_2, (B22_2, (F_1, F_4)FinB)B2F)B2_2)2)BF)R;', format=1)
    UG = UnrootedForest(G)
    print(G)
    clade = set(['F_1', 'F_2', 'F_3', 'B11_1', 'B12_1', 'B21_1'])
    test_map = UG.get_clade_attaching_edges(clade)
    print('Clades induced by', clade)
    print(test_map)
    print('Leaves of first clade:')
    print(UG.get_clade_leaves(test_map[0]))
    # print('Leaves of second clade:')
    # print(UG.get_clade_leaves(test_map[1]))

    print('Detecting displacements, test 1:')
    S = Tree("(F, ((B11, B12)B1, (B21, B22)B2)B)R;", format=1)
    name_map = {c.name : c.name.split('_')[0] for c in G}
    print('G:')
    print(G.get_ascii(attributes=['name']))
    print('S:')
    print(S.get_ascii(attributes=['name']))
    clade = set(['F_1', 'F_2', 'F_3', 'F_4'])
    US = UnrootedForest(S)
    print('Clade displacements:')
    print(UG.get_clade_displacements(clade, US, name_map))
    print('Optimal location:')
    print(UG.get_optimal_clade_location(clade, US, name_map, True))

    print('Detecting displacements, test 2:')
    G = Tree('((B1, (F1_1, B2)F1B2)FB12, ((B31, B32)B3, ((B41_1, ((F1_2, F2_2)F12, B42_1)B42_1F)B4_1F, (B41_2, B42_2)B4_2)B4F)B34F)R;', format=1)
    S = Tree('((F1, F2)F, (B1, B2, (B31, B32)B3, (B41, B42)B4)B)R;', format=1)
    print(G.get_ascii(attributes=['name']))
    print(S.get_ascii(attributes=['name']))
    name_map = {c.name : c.name.split('_')[0] for c in G}
    UG = UnrootedForest(G)
    US = UnrootedForest(S, binary=False)
    clade=set(['F1_1', 'F1_2', 'F2_2'])
    print('Clade displacements:')
    print(UG.get_clade_displacements(clade, US, name_map))
    
    ### Alignment processing testing
    test_seqs = ['AAA-M-KKK', 'ALL-X-LKL', '---KX---L', '---Kx---L']
    test_aln = MultipleSeqAlignment([SeqRecord(Seq(seq), str(seqid), '', '') for seqid, seq in enumerate(test_seqs)])         
    print(test_aln)
    proc_test_aln = process_alignment(test_aln, min_row_nb_of_aa=1, 
                                  max_column_gap_proportion=.8, max_row_gap_proportion=.9,
                                  remove_constant_columns=True)
    print(proc_test_aln)


