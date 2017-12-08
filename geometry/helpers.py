"""Helper functions that are used to support geometric functions
"""

from __future__ import division, print_function
import numpy as np
import copy

__all__ = ['check_clockwise', 'check_counterclockwise', 'point_in_poly',
           'euclid_dist', 'lines_intersect', 'find_intersection_point',
           'angle_between', 'find_mutually_visible', 'point_on_line',
           'distance_point_2_line', 'distance_point_2_seg', 'point_on_left',
           'point_in_concave_poly', 'point_on_infinite_line']


def euclid_dist(p1, p2):
    """Returns the euclidean distance between two points

    Args:
        p1 ([float, float]): First point
        p2 ([float, float]): Second point

    Returns:
        The euclidean distance between the points
    """
    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]
    return np.sqrt(dx * dx + dy * dy)


def angle_between(v1, v2):
    """Returns the angle from v1 to v2

    Args:
        v1 ([float, float]): First vector
        v2 ([float, float]): Second vector

    Returns:
        The angle (in radians) between v1 and v2
    """
    v1_n = v1 / np.linalg.norm(v1)
    v2_n = v2 / np.linalg.norm(v2)
    return np.arccos(np.clip(np.dot(v1_n, v2_n), -1., 1.))


def check_clockwise(vertices):
    """Checks whether a set of vertices wind in clockwise order

    Adapted from http://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order

    Args:
        vertices (list): A list of [x,y] vertices

    Returns:
        bool indicating if the points wind in a clockwise order
    """
    tot = 0
    nv = len(vertices)
    i = nv - 1
    j = 0
    while j < nv:
        tot += (vertices[j][0] - vertices[i][0]) * (vertices[j][1] + vertices[i][1])
        i = j
        j += 1
    return tot > 0


def check_counterclockwise(vertices):
    """Checks whether a set of vertices wind in counter-clockwise order

    Adapted from http://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order

    Args:
        vertices (list): A list of [x,y] vertices

    Returns:
        bool indicating if the points wind in a counter-clockwise order
    """
    tot = 0
    nv = len(vertices)
    i = nv - 1
    j = 0
    while j < nv:
        tot += (vertices[j][0] - vertices[i][0]) * (vertices[j][1] + vertices[i][1])
        i = j
        j += 1
    return tot < 0


def point_in_poly(point, vertices, count_edges=False):
    """Determines whether a point is in a convex polygon

    Adapted from http://stackoverflow.com/questions/1119627/how-to-test-if-a-point-is-inside-of-a-convex-polygon-in-2d-integer-coordinates

    Args:
        point ([float, float]): The (x,y) point to test
        vertices (list): A list of (x,y) vertices forming a convex hull
        count_edges (bool): Should points on the edge be part of the poly?

    Returns:
        bool indicating whether the point is on the inside of the polygon
    """
    pside = None
    for i in range(len(vertices)):
        a = vertices[i]
        b = vertices[(i + 1) % len(vertices)]
        affine_seg = b - a
        affine_pt = point - a
        xprod = np.cross(affine_seg, affine_pt)
        if xprod < 0:
            cside = -1
        elif xprod > 0:
            cside = 1
        else:
            return count_edges  # On an edge
        if pside is None:
            pside = cside
        else:
            if pside != cside:
                return False
    return True

def point_in_concave_poly(point, vertices, count_edges=False):
    """Determines whether a point is in a non-convex polygon

    Adapted from http://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/

    Args:
        point ([float, float]): The (x,y) point to test
        vertices (list): A list of (x,y) vertices forming a non-convex hull
        count_edges (bool): Should points on the edge be part of the poly?

    Returns:
        bool indicating whether the point is on the inside of the polygon
    """
    max_vx = max([x for x, y in vertices])
    p_ray = point + np.array([max_vx, 0])
    n_intersects = 0
    for i in range(len(vertices)):
        v1 = vertices[i]
        v2 = vertices[(i + 1) % len(vertices)]
        if count_edges:
            if point_on_line(point, v1, v2):
                return True
        if lines_intersect(point, p_ray, v1, v2):
            n_intersects += 1
    return (n_intersects % 2) == 1

def point_on_infinite_line(point, p1, p2, tol=1e-6):
    """Determines whether a point lies on an infinite line

    Allows for a small tolerance

    Args:
        point ([float, float]): the point to test
        p1 ([float, float]): point 1 defining the line
        p2 ([float, float]): point 2 defining the line
        tol (float): the maximum distance of the point from the 1d line

    Returns:
        bool indicating whether the point lies within tol of the infinite line
    """
    # Special case: vertical lines:
    if p1[0] == p2[0]:
        if point[0] != p1[0]:
            return False
        else:
            miny = min(p1[1], p2[1])
            maxy = max(p1[1], p2[1])
            return miny <= point[1] <= maxy
    prop_along = (point[0] - p1[0]) / (p2[0] - p1[0])
    guess_y = p1[1] + prop_along * (p2[1] - p1[1])
    return abs(guess_y - point[1]) < tol

def point_on_line(point, seg_p1, seg_p2, tol=1e-6):
    """Determines whether a point lies on a line segment

    Allows for a small tolerance

    Args:
        point ([float, float]): the point to test
        seg_p1 ([float, float]): the beginning of the line segment
        seg_p2 ([float, float]): the end of the line segment
        tol (float): the maximum distance of the point from the 1d line

    Returns:
        bool indicating whether the point lies within tol of the line
    """
    # Special case: vertical lines:
    if seg_p1[0] == seg_p2[0]:
        if point[0] != seg_p1[0]:
            return False
        else:
            miny = min(seg_p1[1], seg_p2[1])
            maxy = max(seg_p1[1], seg_p2[1])
            return miny <= point[1] <= maxy
    prop_along = (point[0] - seg_p1[0]) / (seg_p2[0] - seg_p1[0])
    if prop_along < 0 or prop_along > 1:
        return False
    guess_y = seg_p1[1] + prop_along * (seg_p2[1] - seg_p1[1])
    return abs(guess_y - point[1]) < tol


def _ccw(a, b, c):
    return (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])


def lines_intersect(a, b, c, d):
    """Checks whether the two line segments ab and cd intersect

    From http://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/

    Args:
        a, b, c, d ([float,float]): endpoints of the two segments

    Returns:
        bool if the lines intersect
    """
    return _ccw(a, c, d) != _ccw(b, c, d) and _ccw(a, b, c) != _ccw(a, b, d)


def find_intersection_point(a, b, c, d):
    """Finds the point of intersections between segments ab and cd

    From http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect

    Args:
        a, b, c, d ([float,float]): endpoints of the two segments

    Returns:
        None if no intersection, [x,y] point if intersection exists
    """
    p = np.array(a)
    r = np.array(b) - p
    q = np.array(c)
    s = np.array(d) - q

    qp_cross_r = np.cross(q - p, r)
    r_cross_s = np.cross(r, s)
    if qp_cross_r == 0 and r_cross_s == 0:
        # Make sure they overlap
        if r[0] > 0: # AB goes to the right
            if a[0] < c[0] < b[0]:
                return q
            elif a[0] < d[0] < b[0]:
                return q+s
            else:
                return None
        elif r[0] < 0: # AB goes to the left
            if b[0] < c[0] < a[0]:
                return q
            elif b[0] < d[0] < a[0]:
                return q+s
            else:
                return None
        else: # Vertical
            if r[1] > 0:
                if a[1] < c[1] < b[1]:
                    return q
                elif a[1] < d[1] < b[1]:
                    return q+s
                else:
                    return None
            elif r[1] < 0:
                if b[1] < c[1] < a[1]:
                    return q
                elif b[1] < d[1] < a[1]:
                    return q+s
                else:
                    return None
    if r_cross_s == 0 and qp_cross_r != 0:
        # Parallel and not collinear
        return None
    u = qp_cross_r / r_cross_s
    t = np.cross(q - p, s) / r_cross_s
    if 0 <= t <= 1 and 0 <= u <= 1:
        return p + t * r
    else:
        return None


def point_on_left(pt, seg):
    """Finds whether a point is on the left of a given line defined by a segment

    Args:
        pt ([float, float]): (x,y) point to test
        seg ([[float, float], [float, float]]): line segment

    Returns:
        bool value if the point is to the left of the line segment
    """
    sa, sb = seg
    cprod = np.cross(sb - sa, pt - sa)
    return cprod > 0


def find_mutually_visible(exterior_hull, interior_hull):
    """Finds mututally visible points for doing clipping with interior holes

    Logic based on https://www.geometrictools.com/Documentation/TriangulationByEarClipping.pdf

    Args:
        exterior_hull (list): set of vertices in exterior hull
        interior_hull (list): set of vertices in interior hull

    Returns:
        A list of four items:
        * Index of the point on the exterior hull
        * Index of the point on the interior hull
        * (x,y) point on the exterior hull
        * (x,y) point on the interior hull
    """
    # Find the interior point with the greatest x
    xs_in = [p[0] for p in interior_hull]
    gx_in_idx = [i for i, x in enumerate(xs_in) if x == max(xs_in)][0]
    gx_in = interior_hull[gx_in_idx]
    # Check to see whether the ray from gx_in to the right intersects with any
    # edges of the outer poly -- find nearest
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
                closest_seg = [i, j]
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
    # Otherwise, find the vertex in there that minimizes the angle between the
    # angle with the inner point as center
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


def distance_point_2_line(point, seg):
    """Finds the minimum distance and closest point between a point and a line

    Args:
        point ([float, float]): (x,y) point to test
        seg ([[float, float], [float, float]]): two points defining the line

    Returns:
        A list of two items:
        * Distance between the point and line
        * The (x,y) value on the line that is the closest point
    """
    dseg = seg[1] - seg[0]
    dpt = point - seg[0]
    proj = (np.dot(dpt, dseg) / np.dot(dseg, dseg)) * dseg
    dist = np.linalg.norm(dpt, proj)
    return dist, seg[0] + proj


def distance_point_2_seg(point, seg):
    """Finds the minimum distance and closest point between a point and a segment

    Args:
        point ([float, float]): (x,y) point to test
        seg ([[float, float], [float, float]]): endpoints of the segment

    Returns:
        A list of two items:
        * Distance between the point and segment
        * The (x,y) value on the segment that is the closest point
    """
    dseg = seg[1] - seg[0]
    # Segment is actually a point
    if all(dseg == np.array([0, 0])):
        return np.linalg.norm(point - seg[0])
    dpt = point - seg[0]
    prop_proj = (np.dot(dpt, dseg) / np.dot(dseg, dseg))
    if prop_proj < 0:
        return np.linalg.norm(point - seg[0])
    elif prop_proj > 1:
        return np.linalg.norm(point - seg[1])
    else:
        proj = dseg * prop_proj
        return np.linalg.norm(dpt - proj)
