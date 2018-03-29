"""Classes for handling non-convex polygons
"""

from __future__ import division, print_function
import numpy as np
from .helpers import *
from .convex import *
import copy

__all__ = ['approximate_convex_decomposition']

# Decompose a non-convex shape into approximately convex hulls
# From Lien & Amato (2006): Approximate convex decomposition of polygons


class _viz_tree_node(object):
    def __init__(self, point, index, parent, total_dist):
        self._point = point
        self._parent = parent
        self._total_dist = total_dist
        self._children = []

    def dist_to_pt(self, point):
        return np.linalg.norm(self._point - point)

    def add_child(self, point, index):
        d = self.dist_to_pt(point)
        self._children.append(_viz_tree_node(point, index,
                                             self, self._total_dist + d))

    def get_greatest_dist(self):
        if len(self._children) == 0:
            return self._point, self._total_dist
        else:
            maxdist = self._total_dist
            bestpt = self._point
            for c in self._children:
                npt, md = c.get_greatest_dist()
                if md > maxdist:
                    maxdist = md
                    bestpt = npt
            return bestpt, maxdist

    def get_leaves(self):
        if len(self._children) == 0:
            return [self._point]
        else:
            leaves = []
            for c in self._children:
                leaves += c.get_leaves
            return leaves


class _viz_tree(object):
    def __init__(self, poly, head_node_idx=0):
        self._head = _viz_tree_node(poly[head_node_idx], head_node_idx, None, 0)
        self._poly = poly
        rem_vertices = range(len(self._poly))
        rem_vertices.remove(head_node_idx)
        while len(rem_vertices) > 0:
            leaves = self._head.get_leaves()
            # Find whether each remaining vertex can get added

    def _segment_intersects_poly(self, idx1, idx2):
        # If they're right next to each other, there's a common path
        if abs(idx1 - idx2) == 1:
            return False
        else:
            return False  # TO BE FIXED _ NOT RIGHT!!!

# Calculate straight-line concavity for a notch in a pocket
# Input:
#  pocket: a list of vertices forming the pocket


def _sl_concavity(pocket, bridge=None):
    if bridge is None:
        bridge = [pocket[0], pocket[-1]]
        dists = [0] + [distance_point_2_seg(p, bridge) for p in pocket[1:-1]] + [0]
    else:
        dists = [distance_point_2_seg(p, bridge) for p in pocket]
    maxdist = 0
    furthestpt = None
    for i, d in enumerate(dists):
        if d > maxdist:
            maxdist = d
            furthestpt = pocket[i]
    return maxdist, furthestpt


# As above, but for shortest-path concavity
def _sp_concavity(pocket, bridge):
    if bridge is None:
        bridge = [pocket[0], pocket[-1]]
    # First, split the pocket into three parts, based on the parallel
    # projections from the ends of the pocket
    pocket = copy.deepcopy(pocket)
    if check_clockwise(pocket):
        pocket.reverse()
    b_minus = pocket[0]
    b_plus = pocket[-1]
    tan_dir = [-(b_plus[1] - b_minus[1]), b_plus[0] - b_minus[0]]
    scaling_factor = max(map(max, pocket))
    tan_ray_minus = b_minus + tan_dir * scaling_factor
    tan_ray_plus = b_plus + tan_dir * scaling_factor
    # Go through the points until we hit an intersection with the tangent line
    i = 1
    b_m_spl = None
    b_p_spl = None
    while b_m_spl is None:
        m_spl_pt = find_intersection_point(b_minus, tan_ray_minus,
                                           pocket[i - 1], pocket[i])
        if m_spl_pt is not None:
            b_m_spl = i
            P_minus = pocket[:b_m_spl] + [m_spl_pt]
            rem_pocket = [m_spl_pt] + pocket[b_m_spl:]
        else:
            i += 1
    i = 1
    while b_p_spl is None:
        p_spl_pt = find_intersection_point(b_plus, tan_ray_plus,
                                           rem_pocket[i - 1], rem_pocket[i])
        if p_spl_pt is not None:
            b_p_spl = b_m_spl + i
            P_mid = rem_pocket[:i] + [p_spl_pt]
            P_plus = [p_spl_pt] + rem_pocket[i:]
    return 0  # TO BE DONE LATER


def _acd_inner(hole, shell, tau, witness_fnc):
    # Find the antinoal pair for this hole, and point to use for cut point
    node, antinode = _find_antinodal(hole, witness_fnc)
    ndist, nvert = _find_min_vertex_dist(node, shell)
    andist, anvert = _find_min_vertex_dist(antinode, shell)
    if ndist <= andist:
        cut_vert = node
        cut_shell = nvert
        other_vert = antinode
        other_shell = anvert
    else:
        cut_vert = antinode
        cut_shell = anvert
        other_vert = node
        other_shell = nvert

    # Try to form a cut with the vertices
    for vtx, v_cut in [(cut_shell, cut_vert), (other_shell, other_vert)]:
        vidx = _find_point_in_vlist(vtx, shell)
        v_pre = shell[(vidx - 1) % len(shell)]
        v_post = shell[(vidx + 1) % len(shell)]
        if _is_resolved(vidx, shell, v_cut, [v_pre, vidx, v_post]):
            # This is good = cut open the shell and add the hole
            if not check_clockwise(hole):
                hole.reverse()
            hidx = _find_point_in_vlist(v_cut, hole)
            comb = shell[:(vidx + 1)] + hole[hidx:] + hole[:(hidx + 1)] + shell[vidx:]
            return comb
    raise Exception("Uh oh -- can't match either cut point (thought I'd never get here)")


def _acd_outer(vertices, tau, witness_fnc):
    d, witness, pocket = _find_witness(vertices, witness_fnc)
    if d < tau:
        return [vertices]  # We're below threshold or it's already convex
    cut_v = _find_cut_heur(witness, vertices, pocket)
    poly_1 = []
    poly_2 = []
    vidx = 0
    on_poly_2 = False
    for vidx in range(len(vertices)):
        this_v = vertices[vidx]
        # Are we at a cut point? If so, add to both polygons
        if all(this_v == witness) or all(this_v == cut_v):
            poly_1.append(this_v)
            poly_2.append(this_v)
            on_poly_2 = not on_poly_2
        else:
            if on_poly_2:
                poly_2.append(this_v)
            else:
                poly_1.append(this_v)
    return _acd_outer(poly_1, tau, witness_fnc) + _acd_outer(poly_2, tau, witness_fnc)


# Find the point with the highest measure of concavity
# Returns the concavity measure, the point
def _find_witness(vertices, witness_fnc):
    wrapped = gift_wrap(vertices)
    if len(wrapped) == len(vertices):
        # We have a convex hull already
        return 0, None, []
    vertices = copy.copy(vertices)
    # Wrap around to align wrapped and vertices
    while not all(wrapped[0] == vertices[0]):
        vertices = vertices[1:] + [vertices[0]]
    # Find the areas with most convexity
    # Make pockets and find the greatest witness distance
    gidx = 0
    pockets = []
    next_pocket = []
    for vi in range(len(vertices)):
        vtx = vertices[vi]
        if all(wrapped[gidx] == vtx):
            # End of a pocket
            if len(next_pocket) > 0:
                next_pocket.append(vtx)
                pockets.append(next_pocket)
                next_pocket = []
            # Special case: start of a pocket
            if gidx < (len(wrapped) - 1):
                if not all(wrapped[gidx + 1] == vertices[vi + 1]):
                    next_pocket = [vtx]
            # Otherwise part of the outer hull
            gidx = (gidx + 1) % len(wrapped)
        else:
            # Otherwise, start carving out a pocket
            next_pocket.append(vtx)
    if len(next_pocket) > 0:
        pockets.append(next_pocket)
    # Find the witness point -- that with the greatest distance
    maxdist = 0
    witness = None
    pock = None
    for p in pockets:
        d, wp = witness_fnc(p)
        if d > maxdist:
            maxdist = d
            witness = wp
            pock = p
    return maxdist, witness, pock

# Find a point to cut


def _find_cut_heur(cut_vertex, vertices, pocket):
    cut_idx = _find_point_in_vlist(cut_vertex, vertices)
    nvs = len(vertices)
    v_pre = vertices[(cut_idx - 1) % nvs]
    v_post = vertices[(cut_idx + 1) % nvs]
    ptr = v_post
    ptr_idx = cut_idx + 1
    best_v = None

    max_score = 0
    while not all(ptr == cut_vertex):
        vec = ptr - cut_vertex
        score = 1. / np.linalg.norm(vec)  # Find the shortest length that works
        if score > max_score:
            if _is_resolved(cut_idx, vertices, ptr, pocket):
                max_score = score
                best_v = ptr
        ptr_idx = (ptr_idx + 1) % len(vertices)
        ptr = vertices[ptr_idx]
    return best_v

# Does this make two hulls when split along a point?
# Can't be one of the pocket points


def _is_resolved(vidx, vertices, v_cut, pocket):
    vtx = vertices[vidx]
    nvs = len(vertices)
    v_pre = vertices[(vidx - 1) % nvs]
    v_post = vertices[(vidx + 1) % nvs]
    # Can't be one of the adjacent vertices
    if any([all(v_cut == v) for v in pocket]):
        return False
    # Otherwise must be on the left of either the pre->v or v->post segments
    if not (point_on_left(v_cut, [v_pre, vtx]) or point_on_left(v_cut, [vtx, v_post])):
        return False
    # And finally, must not intersect with any other line segments
    for i in range(1, nvs):
        v_before = vertices[i - 1]
        v_after = vertices[i]
        if not (all(vtx == v_before) or all(vtx == v_after) or all(v_cut == v_before) or all(v_cut == v_after)):
            if lines_intersect(vtx, v_cut, v_before, v_after):
                return False
    return True

# Finds the antinodal pairs in a hole
# Note: naive implementation -- not the most efficient


def _find_antinodal(hole, witness_fnc):
    node = None
    antinode = None
    maxdist = 0
    for i in range(len(hole)):
        n = hole[i]
        d, an = witness_fnc(hole, [n, n])
        if d > maxdist:
            maxdist = d
            node = n
            antinode = an
    return [node, antinode]

# Finds the closest point on a polygon to a give point


def _find_min_vertex_dist(p, poly):
    mindist = None
    minvert = None
    for v in poly:
        d = np.linalg.norm(p - v)
        if mindist is None:
            mindist = d
            minvert = v
        elif d < mindist:
            mindist = d
            minvert = v
    return mindist, minvert

# Finds the minimum distance between two polygons


def _find_min_poly_dist(poly1, poly2):
    dists = [_find_min_vertex_dist(v, poly2)[0] for v in poly1]
    return min(dists)


def _find_point_in_vlist(pt, vlist):
    for i, v in enumerate(vlist):
        if all(pt == v):
            return i


def approximate_convex_decomposition(vertex_lists, tau,
                                     witness_type="shortest_length"):
    """Decomposes a non-convex shape into approimate convex hulls

    Uses ACD method of From Lien & Amato (2006): Approximate convex
    decomposition of polygons

    Args:
        vertex_lists (list): A list of list of vertices. The first position
        must always be the outer hull. Further are holes
        tau (float): minimum allowable convexity angle for decomposition
        witness_type (str): currently not implemented

    Returns:
        A list of vertex lists, each of which is a component convex hull
    """
    if witness_type == 'shortest_length':
        witness_fnc = _sl_concavity
    else:
        raise NotImplementedError("Only shortest_length is an acceptable witness_type")
    outer = vertex_lists[0]
    if len(vertex_lists) > 1:
        holes = vertex_lists[1:]
        # First, reorder holes to get at the ones closest to the boundary first
        hole_dists = [_find_min_poly_dist(h, outer) for h in holes]
        ordered_holes = [h for (d, h) in sorted(zip(hole_dists, holes))]
        for h in ordered_holes:
            outer = _acd_inner(h, outer, tau, witness_fnc)
    return _acd_outer(outer, tau, witness_fnc)
