from __future__ import division
import numpy as np
import copy

__all__ = ['check_clockwise', 'check_counterclockwise','point_in_poly', 'euclid_dist', 'lines_intersect',
           'find_intersection_point', 'angle_between', 'find_mutually_visible', 'point_on_line',
           'distance_point_2_line', 'distance_point_2_seg', 'point_on_left']


# Returns the euclidean distance between two points
def euclid_dist(p1, p2):
    dx = p1[0]-p2[0]
    dy = p1[1]-p2[1]
    return np.sqrt(dx*dx + dy*dy)

# Gets the angle between two vectors
def angle_between(v1, v2):
    v1_n = v1 / np.linalg.norm(v1)
    v2_n = v2 / np.linalg.norm(v2)
    return np.arccos(np.clip(np.dot(v1_n,v2_n), -1., 1.))

# Returns a boolean to determine whether a set of points winds in a clockwise order
# Adapted from http://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
def check_clockwise(vertices):
    tot = 0
    nv = len(vertices)
    i = nv-1
    j = 0
    while j < nv:
        tot += (vertices[j][0]-vertices[i][0]) * (vertices[j][1]+vertices[i][1])
        i = j
        j += 1
    return tot > 0

# As above, but reversed
def check_counterclockwise(vertices):
    tot = 0
    nv = len(vertices)
    i = nv-1
    j = 0
    while j < nv:
        tot += (vertices[j][0]-vertices[i][0]) * (vertices[j][1]+vertices[i][1])
        i = j
        j += 1
    return tot < 0

# Returns boolean determining if a point is within a given convex polygon
# Adapted from http://stackoverflow.com/questions/1119627/how-to-test-if-a-point-is-inside-of-a-convex-polygon-in-2d-integer-coordinates
def point_in_poly(point, vertices):
    pside = None
    for i in range(len(vertices)):
        a = vertices[i]
        b = vertices[(i+1) % len(vertices)]
        affine_seg = b - a
        affine_pt = point - a
        xprod = np.cross(affine_seg, affine_pt)
        if xprod < 0:
            cside = -1
        elif xprod > 0:
            cside = 1
        else:
            return False # On an edge (which counts as inside)
        if pside is None:
            pside = cside
        else:
            if pside != cside:
                return False
    return True

# Returns whether a point lies on a line segment (within tiny tolerance)
def point_on_line(point, seg_p1, seg_p2, tol = 1e-6):
    # Special case: vertical lines:
    if seg_p1[0] == seg_p2[0]:
        if point[0] != seg_p1[0]:
            return False
        else:
            miny = min(seg_p1[1], seg_p2[1])
            maxy = max(seg_p1[1], seg_p2[1])
            return miny <= point[1] <= maxy
    prop_along = (point[0]-seg_p1[0]) / (seg_p2[0]-seg_p1[0])
    if prop_along < 0 or prop_along > 1:
        return False
    guess_y = seg_p1[1] + prop_along*(seg_p2[1] - seg_p1[1])
    return abs(guess_y - point[1]) < tol

# Checks whether the two line segments ab and cd intersect (just whether -- not where)
# From http://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
def _ccw(a, b, c):
    return (c[1]-a[1])*(b[0]-a[0]) > (b[1]-a[1])*(c[0]-a[0])

def lines_intersect(a, b, c, d):
    return _ccw(a, c, d) != _ccw(b, c, d) and _ccw(a, b, c) != _ccw(a, b, d)

# As above, but also calculates the point of intersection
# Returns None if no intersection point exists
def find_intersection_point(a, b, c, d):
    # From http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
    p = np.array(a)
    r = np.array(b) - p
    q = np.array(c)
    s = np.array(d) - q

    qp_cross_r = np.cross(q-p, r)
    r_cross_s = np.cross(r, s)
    if qp_cross_r == 0 and r_cross_s == 0:
        # Parallel and collinear
        if p[0] < q[0]:
            return p
        else:
            return q
    if r_cross_s == 0 and qp_cross_r != 0:
        # Parallel and not collinear
        return None
    u = qp_cross_r / r_cross_s
    t = np.cross(q-p, s) / r_cross_s
    if 0 <= t <= 1 and 0 <= u <= 1:
        return p + t*r
    else:
        return None

# Finds whether a point is on the left of a given line defined by a segment
def point_on_left(pt, seg):
    sa, sb = seg
    cprod = np.cross(sb - sa, pt - sa)
    return cprod > 0


# Finds mututally visible points for doing clipping with interior holes
# Returns: [index of outer shell point, index of inner point, outer point, inner point]
# Logic based on https://www.geometrictools.com/Documentation/TriangulationByEarClipping.pdf
def find_mutually_visible(exterior_hull, interior_hull):
    # Find the interior point with the greatest x
    xs_in = [p[0] for p in interior_hull]
    gx_in_idx = [i for i, x in enumerate(xs_in) if x == max(xs_in)][0]
    gx_in = interior_hull[gx_in_idx]
    # Check to see whether the ray from gx_in to the right intersects with any edges of the outer poly -- find nearest
    i = len(exterior_hull) - 1
    closest_seg = None
    closest_x = None
    closest_pt = None
    for j in range(len(exterior_hull)):
        hp1 = exterior_hull[i]
        hp2 = exterior_hull[j]
        ray_x = max(abs(hp1[0]), abs(hp2[0])) + gx_in[0]
        isects = find_intersection_point(gx_in, [ray_x, gx_in[1]], hp1, hp2)
        if isects is not None:
            xdiff = isects[0] - gx_in[0]
            if closest_seg is None:
                closest_x = xdiff
                closest_seg = [i,j]
                closest_pt = isects
            elif xdiff < closest_x:
                closest_x = xdiff
                closest_seg = [i, j]
                closest_pt = isects
        i = j
    # Does this point happen to be on the polygon? Then it's the mutually visible one
    for i, v in enumerate(exterior_hull):
        if all(closest_pt == v):
            return i, gx_in_idx, closest_pt, gx_in
    # If not, pick the point on the containing edge with the greatest x-value
    hp1 = exterior_hull[closest_seg[0]]
    hp2 = exterior_hull[closest_seg[1]]
    if hp1[0] > hp2[0]:
        hp = hp1
        hp_idx = closest_seg[0]
    else:
        hp = hp2
        hp_idx = closest_seg[1]
    # Search to make sure there are no vertices from the hull in the interior of the triangle made of the ray and
    # hull point found
    tri = [gx_in, closest_pt, hp]
    ps_in = []
    for i, pt in enumerate(exterior_hull):
        if i not in closest_seg:
            if point_in_poly(pt, tri):
                ps_in.append((i, pt))
    # If there aren't any in, use that endpoint
    if len(ps_in) == 0:
        return hp_idx, gx_in_idx, hp, gx_in
    # Otherwise, find the vertex in there that minimizes the angle between the angle with the inner point as center
    min_ang = 99999
    min_idx = None
    min_pt = None
    ray_vec = closest_pt - gx_in
    for i, p in ps_in:
        pdiff = p - gx_in
        ang = angle_between(ray_vec, pdiff)
        if abs(ang) < min_ang:
            min_ang = abs(ang)
            min_idx = i
            min_pt = p
    return min_idx, hp_idx, min_pt, gx_in


# Finds the minimum distance and closest point between a point and a line
def distance_point_2_line(point, seg):
    dseg = seg[1] - seg[0]
    dpt = point - seg[0]
    proj = (np.dot(dpt, dseg) / np.dot(dseg, dseg)) * dseg
    dist = np.linalg.norm(dpt, proj)
    return dist, seg[0]+proj

# Finds the minimum distance between a point and a line segment
def distance_point_2_seg(point, seg):
    dseg = seg[1] - seg[0]
    # Segment is actually a point
    if all(dseg == np.array([0,0])):
        return np.linalg.norm(point - seg[0])
    dpt = point - seg[0]
    prop_proj = (np.dot(dpt, dseg) / np.dot(dseg, dseg))
    if prop_proj < 0:
        return np.linalg.norm(point - seg[0])
    elif prop_proj > 1:
        return np.linalg.norm(point - seg[1])
    else:
        proj = dseg*prop_proj
        return np.linalg.norm(dpt - proj)
