
import graph_tool as gt
import graph_tool.draw
import numpy as np
import pandas as pd
import scipy.sparse
import scipy as sp
from collections import defaultdict
from tqdm import tqdm
from functools import reduce
import operator
from itertools import product
from collections import namedtuple
from sklearn.decomposition import non_negative_factorization
from functools import partial
import itertools
from warnings import warn
from collections import defaultdict
import difflib
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
from lib.pandas_util import idxwhere
from itertools import starmap
import contextlib
import warnings
from lib.util import info

# Graph generation

def path_to_edgelist(path):
    u = path[0]
    edges = []
    for v in path[1:]:
        edges.append((u, v))
        u = v
    return edges


def is_subseq(x, y):
    # Borrowed from https://stackoverflow.com/a/24017747/1951857
    return all(any(c == ch for c in y) for ch in x)


def single_stranded_graph_from_merged_paths(paths, lengths, depths):
    g = gt.Graph()
    all_edges = []
    for p in paths:
        all_edges.extend(path_to_edgelist(p))
    g.add_edge_list(set(all_edges))
    g.vp['depth'] = g.new_vp('vector<float>')
    g.vp.depth.set_2d_array(depths)
    g.vp['length'] = g.new_vp('int', lengths)  
    g.gp['nsample'] = g.new_gp('int', len(depths))
    g.vp['sequence'] = g.new_vp('object', vals=[[k] for k in range(g.num_vertices())])
    g.ep['flow'] = g.new_ep('vector<float>', val=[1] * g.gp.nsample)
    return g


def single_stranded_graph_with_simulated_depth(paths, depths, length=None, scale_depth_by=1):
    if length is None:
        length = defaultdict(lambda: 1)

    nvertices = max(itertools.chain(*paths.values())) + 1
    nsamples = max(len(x) for x in depths.values())
    
    vertex_length = np.array([length[i] for i in range(nvertices)])

    # FIXME: Dimensions seem to break in weird ways when nsamples==1.
    # NOTE: Hint: That weirdness DOESN'T happen when I manually set g.gp.nsample=1 after
    # graph building. So maybe it's something about nsamples above being set wrong?
    expected_depths = np.zeros((nsamples, nvertices))
    for p in paths:
        expected_depths[:, paths[p]] += np.outer(np.array(depths[p]) * scale_depth_by, np.ones(len(paths[p])))
    
    # TODO: Consider using nbinom
    # # See docs for sp.stats.nbinom
    # sigma_sq = expected_depths + dispersion * expected_depths**2
    # p = expected_depths / sigma_sq
    # n = expected_depths**2 / (sigma_sq - expected_depths)
    _depths = sp.stats.poisson(mu=expected_depths).rvs()

    g = single_stranded_graph_from_merged_paths(
        paths.values(),
        depths=_depths,
        lengths=vertex_length,
    )
    return g

# Graph statistics

def depth_matrix(g, vs=None, samples=None):
    if vs is None:
        vs = g.get_vertices()
    if samples is None:
        samples = np.arange(g.gp.nsample)
    depth = g.vp.depth.get_2d_array(samples)
    return depth[:, vs]


def ep_as_adjaceny(ep):
    return sp.sparse.csc_array(gt.spectral.adjacency(ep.get_graph(), weight=ep)).toarray()


def total_length_x_depth(g):
    return (depth_matrix(g) * g.vp.length.a).sum()


def edit_ratio(ref, query):
    diff = difflib.SequenceMatcher(a=ref, b=query)
    return diff.ratio() * (len(diff.a) + len(diff.b)) / (2 * len(diff.b))


def vertex_description(g, refs=None):
    if refs is None:
        refs = []
    vertex = pd.DataFrame(dict(
        in_degree=g.degree_property_map('in').a,
        out_degree=g.degree_property_map('out').a,
        length=g.vp.length.a,
    ))
    depth = pd.DataFrame(depth_matrix(g).T)
    depth.rename(columns=lambda i: f"d{i}")
    return vertex.join(depth)

def scale_ep(ep, maximum=2):
    g = ep.get_graph()
    return g.new_edge_property('float', vals=ep.a * maximum / ep.a.max())

def scale_vp(vp, maximum=10):
    g = vp.get_graph()
    return g.new_vertex_property('float', vals=vp.a * maximum / vp.a.max())

# Visualization

def draw_graph(g, output_size=(300, 300), ink_scale=0.8, **kwargs):
    kwargs = dict(
        vertex_fill_color=g.new_vertex_property('float', vals=np.linspace(0, 1, num=max(g.get_vertices()) + 1)),
        vertex_text=g.vertex_index,
    ) | kwargs
    return gt.draw.graph_draw(g, output_size=output_size, ink_scale=ink_scale, **kwargs)


def dotplot(pathA, pathB, ax=None, **scatter_kws):
    if ax is None:
        ax = plt.gca()
    pathA = np.asanyarray(pathA)
    pathB = np.asanyarray(pathB)
    length = max(len(pathA), len(pathB)) + 1
    pathA = np.pad(pathA, (0, length - len(pathA)), constant_values=-1)
    pathA = np.pad(pathA, (0, length - len(pathA)), constant_values=-1)
    match = sp.spatial.distance.cdist(pathA.reshape((-1, 1)), pathB.reshape((-1, 1)), metric=lambda x, y: x == y)
    ax.scatter(*np.where(match.T), **(dict(marker='o', s=1) | scatter_kws))
    ax.set_aspect('equal')

# Unitigs

def edge_has_no_siblings(g):
    "Check whether upstream or downstream sibling edges exist for every edge."
    vs = g.get_vertices()
    v_in_degree = g.degree_property_map('in')
    v_out_degree = g.degree_property_map('out')
    e_num_in_siblings = gt.edge_endpoint_property(g, v_in_degree, 'target')
    e_num_out_siblings = gt.edge_endpoint_property(g, v_out_degree, 'source')
    e_has_no_sibling_edges = g.new_edge_property('bool', (e_num_in_siblings.a <= 1) & (e_num_out_siblings.a <= 1))
    return e_has_no_sibling_edges


def label_maximal_unitigs(g):
    "Assign unitig indices to vertices in maximal unitigs."
    no_sibling_edges = edge_has_no_siblings(g)
    g_filt = gt.GraphView(
        g,
        efilt=no_sibling_edges,
        directed=True
    )
    labels, counts = gt.topology.label_components(g_filt, directed=False)
    return labels, counts, g_filt


def maximal_unitigs(g, include_singletons=False):
    no_sibling_edges = edge_has_no_siblings(g)
    g_no_sibling_edges = gt.GraphView(
        g,
        efilt=no_sibling_edges,
        directed=True
    )
    
    # NOTE: The below is equivalent to
    # > isolated_loops = list(gt.topology.all_circuits(g_no_sibling_edges))
    labels0, counts0 = gt.topology.label_components(g_no_sibling_edges, directed=False)
    in_degree = g_no_sibling_edges.degree_property_map('in').a
    out_degree = g_no_sibling_edges.degree_property_map('out').a
    all_vs = np.arange(len(labels0.a))
    isolated_loops = []
    for i, _ in enumerate(counts0):
        unitig_mask = labels0.a == i
        if (in_degree[unitig_mask] == 1).all() and (out_degree[unitig_mask] == 1).all():
            isolated_loops.append(all_vs[unitig_mask])
    
    out = []
    _mark_isolated_loops = []
    for vs in isolated_loops:
        vs = list(vs)
        if include_singletons or (len(vs) > 1):
            out.append((vs, True))
        _mark_isolated_loops.extend(vs)
    isolated_loop_mask = g_no_sibling_edges.new_vertex_property('bool', val=1)
    isolated_loop_mask.a[_mark_isolated_loops] = 0
    
    g_no_siblings_no_loops = gt.GraphView(
        g_no_sibling_edges,
        vfilt=isolated_loop_mask,
        directed=True
    )
    labels1, counts1 = gt.topology.label_components(g_no_siblings_no_loops, directed=False)
    tsort_idx = gt.topology.topological_sort(g_no_siblings_no_loops)
    tsort_labels = labels1.a[tsort_idx]
    for i, c in enumerate(counts1):
        if (c == 1) and not include_singletons:
            continue
        vs = list(tsort_idx[tsort_labels == i])
        if g.edge(vs[-1], vs[0]):
            out.append((vs, True))
        else:
            out.append((vs, False))
    return out


def list_unitig_neighbors(g, vs):
    "The in and out neighbors of a unitig path."
    all_ins = set(g.get_in_neighbors(vs[0]))
    all_outs = set(g.get_out_neighbors(vs[-1]))
    return list(set(all_ins) - set(vs)), list(set(all_outs) - set(vs))


def mutate_add_compressed_unitig_vertex(g, vs, is_cycle, old_depths=None, drop_vs=False):
    if old_depths is None:
        old_depths = depth_matrix(g, vs)
    v = int(g.add_vertex())
    nsample = g.gp.nsample
    in_neighbors, out_neighbors = list_unitig_neighbors(g, vs)
    g.add_edge_list((neighbor, v) for neighbor in in_neighbors)
    g.add_edge_list((v, neighbor) for neighbor in out_neighbors)
    g.vp.length.a[v] = g.vp.length.a[vs].sum()
    g.vp.sequence[v] = reduce(operator.add, (g.vp.sequence[u] for u in vs), [])
    new_depth = (
        (
            old_depths
            * g.vp.length.a[vs]
        ).sum(1) / g.vp.length.a[v]
    )
    g.vp.depth[v] = new_depth
    # FIXME: Consider dropping this assert.
    assert np.allclose(
        (old_depths * g.vp.length.a[vs]).sum(1),
        # FIXME: Does new_depth need to be reshaped before multiplication?
        (new_depth.reshape((-1, 1)) * g.vp.length.a[v]).sum(1)
    )
    if is_cycle:
        g.add_edge(v, v)
    for old_v in vs:
        g.clear_vertex(old_v)
    return g


def mutate_compress_all_unitigs(g):
    unitig_list = maximal_unitigs(g, include_singletons=False)
    all_vs = []
    old_depths = depth_matrix(g)
    for i, (vs, is_cycle) in enumerate(unitig_list):
        assert len(vs) > 1  # All len(vs) == 1 should be removed by include_singletons=False.
        mutate_add_compressed_unitig_vertex(g, vs, is_cycle, old_depths=old_depths[:, vs])
        all_vs.extend(vs)
    g.remove_vertex(set(all_vs), fast=True)
    # I think, but am not sure, that the number of nodes removed will always equal the number of edges removed.
    return g


def label_self_looping_vertices(g):
    return (
        gt.incident_edges_op(
            g,
            direction='in',
            op='max',
            eprop=gt.topology.label_self_loops(g, mark_only=True),
        )
    )

def describe_nodes(g):
    description = pd.DataFrame(dict(
        length=g.vp.length.a,
        in_degree=g.degree_property_map('in').a,
        out_degree=g.degree_property_map('out').a,
        circular=label_self_looping_vertices(g).a,
        sequence=[g.vp.sequence[v] for v in g.iter_vertices()],
    ))
    depth = pd.DataFrame(depth_matrix(g).T)
    return description.join(depth).assign(magnitude=depth.sum(1) * description.length)


def mutate_extract_singletons(g):
    nodes = describe_nodes(g)
    singletons = idxwhere((nodes.in_degree + nodes.out_degree - 2 * nodes.circular) == 0)
    g.remove_vertex(singletons, fast=True)
    return g, nodes.loc[singletons].reset_index(drop=True)

# Flows

def estimate_flow_old(f0, d, weight=None, eps=1e-2, maxiter=100):
    if weight is None:
        weight = np.ones_like(d)

    loss_hist = [np.finfo('float').max]
    f = f0
    for step_i in range(maxiter):
        f_out = f
        f_total_out = f_out.sum(1)
        d_error_out = f_total_out - d
        
        f_in = f_out.T
        f_total_in = f_in.sum(1)
        d_error_in = f_total_in - d
        
        loss_hist.append(np.square(d_error_out).sum() + np.square(d_error_in).sum())
        if loss_hist[-1] == 0:
            break  # This should only happen if d is all 0's.
        loss_ratio = (loss_hist[-2] - loss_hist[-1]) / loss_hist[-2]
        if loss_ratio < eps:
            break
            
        # NOTE: Because of errstate, this function is NOT threadsafe.
        # TODO: Determine if it's safe across multiple processes.
        with np.errstate(divide='ignore', over='ignore'):
            allocation_out = f_out.T * np.nan_to_num(1 / f_total_out, posinf=1, nan=0)
        allocated_d_error_out = (allocation_out * d_error_out).T
        with np.errstate(divide='ignore', over='ignore'):
            allocation_in = f_in.T * np.nan_to_num(1 / f_total_in, posinf=1, nan=0)
        allocated_d_error_in = (allocation_in * d_error_in).T

        # The final step is calculated as a average of the in and out error, weighted
        # by the node weight.
        inv_weight = 1 / weight
        mean_allocated_d_error = (
            ((allocated_d_error_in * inv_weight).T + (allocated_d_error_out * inv_weight))
            * (1 / (inv_weight.reshape((-1, 1)) + (inv_weight.reshape((1, -1)))))
        )
        
        f = (f_out - mean_allocated_d_error)
        # Very rarely floating point precision results in very small, negative values for f.
        assert f.min() >= -1e-20
        # Replace these with 0.
        f[f < 0] = 0
    else:
        warn(f"loss_ratio < eps ({eps}) not achieved in maxiter ({maxiter}) steps. Final loss_ratio={loss_ratio}. Final loss={loss_hist[-1]}.")
    return f


# This is used in the parallel implementation of estimate_all_flows.
def _estimate_flow_old(kwargs):
    return estimate_flow_old(**kwargs)


def estimate_all_flows_old(g, eps=1e-3, maxiter=1000, use_weights=True, jobs=1):
    if jobs != 1:
        pool = Pool(jobs)
        map_f = pool.map
    else:
        map_f = lambda *args, **kwargs: list(map(*args, **kwargs))
    flows = []
    if use_weights:
        weight = g.vp.length.a
    else:
        weight = None
    f0 = sp.sparse.csr_array(gt.spectral.adjacency(g))
    dd = depth_matrix(g)
    flows = map_f(
        _estimate_flow_old,
        [
            dict(
                f0=f0, d=dd[sample_idx], weight=weight, eps=eps, maxiter=maxiter
            )
            for sample_idx in range(g.gp.nsample)
        ],
    )
    return flows


def mutate_add_flows_old(g, flows):
    props = []
    for sample_idx, f in enumerate(flows):
        p = g.new_edge_property('float', val=0)
        for i, j in g.get_edges():
            p[g.edge(i, j)] = f[j, i]
        props.append(p)
    props = gt.group_vector_property(props)
    g.ep['flow'] = props
    return g


def estimate_all_flows(g, eps=0.001, maxiter=1000, use_weights=True):
    if use_weights:
        weight = g.vp.length
    else:
        weight = g.new_vertex_property('float', val=1)
    target_vertex_weight = gt.edge_endpoint_property(g, weight, 'target')
    source_vertex_weight = gt.edge_endpoint_property(g, weight, 'source')
    all_flows = []
    for i in range(g.gp.nsample):
        depth = gt.ungroup_vector_property(g.vp.depth, [i])[0]
        flow = g.new_edge_property('float', val=1)
        flow.a[:] = 1
        loss_hist = [np.finfo('float').max]
        for _ in range(maxiter):
            total_in_flow = gt.incident_edges_op(g, 'in', 'sum', flow)
            in_flow_error = g.new_vertex_property('float', vals=depth.a - total_in_flow.a)
            target_vertex_total_inflow = gt.edge_endpoint_property(g, total_in_flow, 'target')
            target_vertex_error = gt.edge_endpoint_property(g, in_flow_error, 'target')
            with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
                target_vertex_alloc = np.nan_to_num(flow.a / target_vertex_total_inflow.a, posinf=1, nan=0)
            target_vertex_alloc_error = target_vertex_alloc * target_vertex_error.a

            total_out_flow = gt.incident_edges_op(g, 'out', 'sum', flow)
            out_flow_error = g.new_vertex_property('float', vals=depth.a - total_out_flow.a)
            source_vertex_total_outflow = gt.edge_endpoint_property(g, total_out_flow, 'source')
            source_vertex_error = gt.edge_endpoint_property(g, out_flow_error, 'source')
            with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
                source_vertex_alloc = np.nan_to_num(flow.a / source_vertex_total_outflow.a, posinf=1, nan=0)
            source_vertex_alloc_error = source_vertex_alloc * source_vertex_error.a

            loss_hist.append(np.square(in_flow_error.a).sum() + np.square(out_flow_error.a).sum())
            if loss_hist[-1] == 0:
                break  # This should only happen if d is all 0's.
            loss_ratio = (loss_hist[-2] - loss_hist[-1]) / loss_hist[-2]
            if loss_ratio < eps:
                break
            
            with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
                # NOTE: Some values of (source_vertex_weight.a + target_vertex_weight.a)
                # are 0 because these two edge_properties include edge indices
                # for non-existent edges.
                # TODO: Consider running gt.reindex_edges to get rid of these.
                mean_flow_error = g.new_edge_property(
                    'float',
                    vals=(
                        (source_vertex_alloc_error * source_vertex_weight.a)
                        +
                        (target_vertex_alloc_error * target_vertex_weight.a)
                    )
                    / (source_vertex_weight.a + target_vertex_weight.a)
                )
            flow = g.new_edge_property('float', vals=flow.a + mean_flow_error.a)
        all_flows.append(flow)
    return all_flows


def mutate_add_flows(g, flows):
    g.ep['flow'] = gt.group_vector_property(flows)

# Node splitting

Split = namedtuple('Split', ['u', 'v', 'w'])


def splits_from_sparse_encoding(g, v, threshold=1.):
    # Compile tables
    in_neighbors = list(sorted(g.get_in_neighbors(v)))
    num_in_neighbors = len(in_neighbors)
    in_neighbors_label = {k: v for k, v in enumerate(in_neighbors)}
    in_neighbors_onehot = {k: v for k, v in zip(in_neighbors, np.eye(num_in_neighbors))}
    in_neighbors_onehot[None] = np.zeros(num_in_neighbors)

    out_neighbors = list(sorted(g.get_out_neighbors(v)))
    num_out_neighbors = len(out_neighbors)
    out_neighbors_label = {k: v for k, v in enumerate(out_neighbors)}
    out_neighbors_onehot = {k: v for k, v in zip(out_neighbors, np.eye(num_out_neighbors))}
    out_neighbors_onehot[None] = np.zeros(num_out_neighbors)

    # Build observation matrix.
    in_neighbor_flow = []
    for u in in_neighbors:
        in_neighbor_flow.append(g.ep.flow[g.edge(u, v)])
    out_neighbor_flow = []
    for w in out_neighbors:
        out_neighbor_flow.append(g.ep.flow[g.edge(v, w)])
    depth_row = g.vp.depth[v].a
    obs = np.stack(in_neighbor_flow + [depth_row] + out_neighbor_flow).T
    assert obs.min() > -1e-20
    obs[obs < 0] = 0

    # Build code matrix.
    in_neighbor_code = []
    out_neighbor_code = []
    split_idx = {}
    for i, (u, w) in enumerate(product(in_neighbors + [None], out_neighbors + [None])):
        if (u, w) == (None, None):
            # NOTE: I treat the naked vertex specially.
            naked_vertex_idx = i
            pass
        in_neighbor_code.append(in_neighbors_onehot[u])
        out_neighbor_code.append(out_neighbors_onehot[w])
        split_idx[i] = (u, w)
    in_neighbor_code = np.stack(in_neighbor_code)
    out_neighbor_code = np.stack(out_neighbor_code)
    unnormalized_code = np.concatenate([
        in_neighbor_code,
        np.ones((in_neighbor_code.shape[0], 1)),
        out_neighbor_code
    ], axis=1)
    code_magnitude = np.sqrt(np.square(unnormalized_code).sum(1, keepdims=True))
    code = unnormalized_code / code_magnitude

    # Group Matching Pursuit (GMP)
    # Inspired by https://arxiv.org/pdf/1812.10538.pdf
    resid = obs
    atoms = []
    dictionary = np.zeros_like(code)
    normalized_encoding = np.zeros((obs.shape[0], dictionary.shape[0]))
    for _ in range(code.shape[0]):
        loss = np.abs(resid).sum()
        dot = resid @ code.T
        # TODO: Decide how to decide atoms
        # next_atom = np.square(dot).sum(0).argmax()
        next_atom = dot.argmax() % code.shape[0]
        naked_vertex_dot = dot[:, naked_vertex_idx]
        if dot[:, next_atom].sum() <= threshold:
            break
        if next_atom in atoms:
            break
        if next_atom == naked_vertex_idx:
            break  # TODO: Decide if this is a good stopping criterion.
        atoms.append(next_atom)
        dictionary[atoms[-1]] = code[atoms[-1]]
        normalized_encoding, _, _ = non_negative_factorization(obs, n_components=dictionary.shape[0], H=dictionary, update_H=False, alpha_W=0)
        resid = obs - normalized_encoding @ code

    # Iterate through selected atoms as splits.
    encoding = normalized_encoding / code_magnitude.T
    out = []
    for i in atoms:
        if i == naked_vertex_idx:
            # Don't return the naked vertex here. It'll be returned
            # later.
            continue
        u, w = split_idx[i]
        out.append((Split(u, v, w), g.vp.length[v], encoding[:, i]))
    # Remaining depth must also include any depth assigned to the naked encoding.
    remaining_depth = g.vp.depth[v] - encoding.sum(1) + encoding[:, naked_vertex_idx]
    if not np.allclose(remaining_depth, 0):
        out.append((Split(None, v, None), g.vp.length[v], remaining_depth))
    return out
        
        
def build_tables_from_splits(split_list, start_idx):
    """Generate edges to and from new, split vertices.
    
    Note that if splits from adjacent parents are not
    reciprocated, no new edge is produced.
    
    """
    split_idx = {}
    upstream = defaultdict(list)
    downstream = defaultdict(list)
    length = []
    depth = []
    for idx, (split, l, d) in enumerate(split_list, start=start_idx):
        u, v, w = split
        split_idx[split] = idx
        upstream[(v, w)].append(split)
        downstream[(u, v)].append(split)
        depth.append(d)
        length.append(l)
    assert len(split_list) == len(split_idx)
    return split_idx, upstream, downstream, np.array(length), np.array(depth)
        
        
def new_edges_from_splits(split_list, split_idx, upstream, downstream, start_idx):
    out = []
    for v, (split, _, _) in enumerate(split_list, start=start_idx):
        u_old, v_old, w_old = split
        v = split_idx[split]
        
        # Upstream edges
        if u_old is not None:
            out.append((u_old, v))
        for upstream_split in upstream[(u_old, v_old)]:
            u = split_idx[upstream_split]
            if u is not None:
                out.append((u, v))
            
        # Downstream edges
        if w_old is not None:
            out.append((v, w_old))
        for downstream_split in downstream[(v_old, w_old)]:
            w = split_idx[downstream_split]
            if w is not None:
                out.append((v, w))
    return out
            
            
def mutate_apply_splits(g, split_list):
    """Add edges and drop any parent vertices that were split.
    
    """
    start_idx = g.num_vertices(ignore_filter=True)
    split_idx, upstream, downstream, lengths, depths = (
        build_tables_from_splits(split_list, start_idx=start_idx)
    )
    edges_to_add = list(set(new_edges_from_splits(
        split_list, split_idx, upstream, downstream, start_idx
    )))
    
    # NOTE: Without adding vertices before edges I can get an IndexError
    # running `g.vp.length.a[np.arange(len(lengths)) + start_idx]`.
    # I believe this is because one or more split nodes
    # are without any new edges, and therefore these nodes don't
    # get added implicitly by `g.add_edge_list`.
    # When these accidentally hidden nodes are the highest valued
    # ones, they also don't get implicitly added due to their index.
    # The result is that I'm missing nodes that should actually exist.
    # NOTE: This line returns an unassigned generator. I _think_ all the
    # nodes are still added, but it's not entirely clear.
    # UPDATE: I'm sure the nodes are still added because of the
    # assert afterwards.
    g.add_vertex(n=max(split_idx.values()) - max(g.get_vertices()))
    g.add_edge_list(set(edges_to_add))
    assert max(g.get_vertices()) == max(split_idx.values())
    
    g.vp.length.a[np.arange(len(lengths)) + start_idx] = lengths
    new_depth = depth_matrix(g)
    new_depth[:, np.arange(len(depths)) + start_idx] = depths.T
    g.vp.depth.set_2d_array(new_depth)
    for split, _, _ in split_list:
        g.vp.sequence[split_idx[split]] = g.vp.sequence[split.v]
    # for v in g.iter_vertices():
    #     assert g.vp.sequence[v]
    vertices_to_drop = set(split.v for (split, _, _) in split_list)
    g.remove_vertex(vertices_to_drop, fast=True)
    return g


def splits_for_all_vertices_old(g, split_func):
    out = []
    for v in g.vertices():
        if (v.in_degree() < 2) and (v.out_degree() < 2):
            continue
        else:
            out.extend(split_func(g, v))
    return out


def splits_for_all_vertices(g, split_func, jobs=1):
    if jobs != 1:
        pool = Pool(jobs)
        map_f = pool.starmap
    else:
        map_f = lambda *args, **kwargs: list(starmap(*args, **kwargs))
        
    # NOTE: I _think_ this is why my splits sometimes have integer u and w,
    # but vertex(v).
    # TODO: Consider using int(v) instead of v, for consistency with other
    # code.
    non_linear_vs = [int(v) for v in g.vertices() if (v.in_degree() >= 2) or (v.out_degree() >= 2)]
    splits = map_f(split_func, [(g, v) for v in non_linear_vs])
    return list(itertools.chain.from_iterable(splits))
            

def mutate_split_all_nodes(g, split_func):
    # TODO: Vertices with <= 1 local path will be split into just
    # themselves pointing trivially at their neighbors.
    # These can be dropped as they are effectively a no-op.
    split_list = list(splits_for_all_vertices(g, split_func))
    if len(split_list) > 0:
        g = mutate_apply_splits(g, split_list=split_list)
    return g


def mutate_run_workflow(g, thresh, maxiter=20):
    num_edges, num_vertices = g.num_edges(), g.num_vertices()
    info(f"START (vertices: {num_vertices} edges: {num_edges}).")
    mutate_compress_all_unitigs(g)
    num_edges, num_vertices = g.num_edges(), g.num_vertices()
    info(f"Finished compressing unitigs (vertices: {num_vertices} edges: {num_edges}).")
    extracted_seqs = []
    for i in range(maxiter):
        previous_num_edges, previous_num_vertices = num_edges, num_vertices
        mutate_add_flows(g, estimate_all_flows(g, use_weights=True))
        info("Finished estimating flow.")
        mutate_split_all_nodes(g, partial(splits_from_sparse_encoding, threshold=thresh))
        num_edges, num_vertices = g.num_edges(), g.num_vertices()
        info(f"Finished splitting nodes (vertices: {num_vertices} edges: {num_edges}).")
        mutate_compress_all_unitigs(g)
        num_edges, num_vertices = g.num_edges(), g.num_vertices()
        info(f"Finished compressing unitigs (vertices: {num_vertices} edges: {num_edges}).")
        _, singletons = mutate_extract_singletons(g)
        singletons['extraction_round'] = i
        extracted_seqs.append(singletons)
        num_edges, num_vertices = g.num_edges(), g.num_vertices()
        info(f"Finished extracting singletons (vertices: {num_vertices} edges: {num_edges}).")
        if (
            (num_vertices == 0)
            or (
                (num_edges == previous_num_edges)
                and (num_vertices == previous_num_vertices)
            )
        ):
            info(f"DONE: Converged at round {i}.")
            break
    else:
        info("DONE: Maximum iteration reached.")
    last_extraction = describe_nodes(g)
    last_extraction['extraction_round'] = -1
    num_last_extraction = last_extraction.shape[0]
    info(f"Final graph has {num_last_extraction} vertices remaining.")
    extracted_seqs = pd.concat(extracted_seqs + [last_extraction]).reset_index(drop=True)
    return g, extracted_seqs.sort_values('magnitude', ascending=False)
