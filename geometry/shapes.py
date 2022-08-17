"""Functions that work on collections of shapes
"""

from typing import Tuple, Annotated, Dict
import numpy as np
from .convex import convex_area, convex_centroid

__all__ = ['recenter_polygon', 'centroid_for_shapes',
           'centroid_for_uncomputed_shapes', 'recenter_system',
           'rescale_and_recenter_system', 'rotate_polygon',
           'rotate_system', 'mirror_polygon', 'mirror_system',
           'find_concave_outline']

def recenter_polygon(vertices: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    """Returns a new convex polygon with centroid at (0,0)

    Args:
        vertices (list): list of (x,y) vertices of convex polygon

    Returns:
        A list just like the input with the recentered vertices (but possibly
        transformed into numpy arrays)
    """
    centroid = convex_centroid(vertices)
    new_verts = []
    for v in vertices:
        v = np.array(v)
        new_verts.append(v - centroid)
    return new_verts

def centroid_for_shapes(centroids: List[Tuple[float, float]],
                        areas: List[float] = None) -> Tuple[float, float]:
    """Calculates the centroid for a set of shapes

    Requires pre-computed centroids and areas

    Args:
        centroids (list): list of (x,y) centroids for each shape
        areas (list): list of areas (floats) for each shape (if not given,
        assumes they are all equal)

    Returns:
        The (x,y) position of the weighted centroid (as np.array)
    """
    gc = np.zeros(2)
    area = 0
    if areas is None:
        areas = np.ones(len(centroids))
    for pc, a in zip(centroids, areas):
        gc += np.array(pc)*a
        area += a
    gc /= area
    return np.array(gc)


def centroid_for_uncomputed_shapes(shape_list: List[List[Tuple[float, float]]]) -> Tuple[float, float]:
    """Like centroid_for_shapes but calculates centroids & areas

    Args:
        shape_list (list): a list of list of vertices (one for each shape)

    Returns:
        The (x,y) position of the weighted centroid (as np.array)
    """
    centroids = []
    areas = []
    for s in shape_list:
        centroids.append(convex_centroid(s))
        areas.append(convex_area(s))
    return centroid_for_shapes(centroids, areas)


def recenter_system(shape_list: List[List[Tuple[float, float]]]):
    """Recenters a set of shapes around the centroid of all of them

    Args:
        shape_list (list): a list of list of vertices (one for each shape)

    Returns:
        List of two items:
        * Similar format as input, but transformed so that calculating the
        centroid_for_uncomputed_shapes() on that list returns (0,0)
        * The grand centroid for the system in original coordinates
    """
    centroids = []
    areas = []
    new_shapes = []
    # Decompose each of the individual shapes
    for s in shape_list:
        c = convex_centroid(s)
        a = convex_area(s)
        new_s = []
        for v in s:
            new_s.append(np.array(v) - c)
        centroids.append(c)
        areas.append(a)
        new_shapes.append(new_s)
    # Find the grand centroid & new centers of each shape
    center = centroid_for_shapes(centroids, areas)
    re_centroids = [c - center for c in centroids]
    # Go back and change the vertices of each shape
    final_shapes = []
    for ns,c in zip(new_shapes, re_centroids):
        final_shapes.append([s+c for s in ns])
    return final_shapes, center


def rescale_and_recenter_system(shape_list: List[List[Tuple[float, float]]],
                                total_area: float):
    """Recenters a set of shapes and resizes them to have a total fixed area

    Args:
        shape_list (list): a list of list of vertices (one for each shape)
        total_area (float): the area to fix the shapes to

    Returns:
    List of two items:
        * Similar format as input, but transformed so that calculating the
        `centroid_for_uncomputed_shapes()` on that list returns (0,0) and summing
        the areas gets to `total_area`
        * The grand centroid for the system in original coordinates
    """
    centroids = []
    areas = []
    new_shapes = []
    # Decompose each of the individual shapes
    for s in shape_list:
        c = convex_centroid(s)
        a = convex_area(s)
        new_s = []
        for v in s:
            new_s.append(np.array(v) - c)
        centroids.append(c)
        areas.append(a)
        new_shapes.append(new_s)
    # Find the grand centroid & new centers of each shape
    center = centroid_for_shapes(centroids, areas)
    re_centroids = [c - center for c in centroids]
    # Find rescaling factor
    tot_a = sum(areas)
    dim_scale = np.sqrt(total_area / tot_a)
    # Go back and change the vertices of each shape
    final_shapes = []
    for ns,c in zip(new_shapes, re_centroids):
        final_shapes.append([(s+c)*dim_scale for s in ns])
    return final_shapes, center

def rotate_polygon(vertices, angle, center_point = [0., 0.]):
    """Rotates a shape around a given point (the origin)

    Args:
        vertices (list): A list of (x,y) vertices
        angle (float): Angle in radians to rotate counterclockwise
        center_point ([float, float]): (x,y) point to rotate around

    Returns:
        A list of vertices rotated around the center point
    """
    np_o = np.array(center_point)
    np_vs = [np.array(v) - np_o for v in vertices]
    rot_mat = np.array([[np.cos(angle), -np.sin(angle)],
                        [np.sin(angle), np.cos(angle)]])
    return [np.dot(rot_mat, v)+np_o for v in np_vs]

def rotate_system(shape_list, angle, center_point = None):
    """Rotates a set of shapes around a given point

    If no center point is given, assume the center of mass of the shape

    Args:
        shape_list (list): A list of list of (x,y) vertices
        angle (float): Angle in radians to rotate counterclockwise
        center_point ([float, float]): (x,y) point to rotate around

    Returns:
        A new shape list with rotated vertices
    """
    if center_point is None:
        center_point = centroid_for_uncomputed_shapes(shape_list)
    return [rotate_polygon(s, angle, center_point) for s in shape_list]

def mirror_polygon(vertices: List[Tuple[float, float]],
                   axes: Tuple[bool, bool]=(False, True),
                   center_point: Tuple[float, float]=None):
    """Mirrors a polygon around an x or y line

    If center_point is None, mirror around the center of the shape

    Args:
        vertices (list): A list of (x,y) vertices
        axes ([bool, bool]): Whether to mirror around the (x,y) axes
        center_point ([float, float]): (x,y) point to mirror around

    Returns:
        A new polygon with rotated vertices
    """
    if center_point is None:
        center_point = convex_centroid(vertices)
    xm = -1 if axes[0] else 1
    ym = -1 if axes[1] else 1
    return [np.array([xm*(v[0]-center_point[0])+center_point[0],
                      ym*(v[1]-center_point[1])+center_point[1]]) for v
            in vertices]

def mirror_system(shape_list: List[List[Tuple[float, float]]],
                  axes: Tuple[bool, bool]=(False, True),
                  center_point: Tuple[float, float]=None):
    """Mirrors a polygon around an x or y line

    Mirrors around the center of the system if center_point is None

    Args:
        shape_list (list): A list of list of (x,y) vertices
        axes ([bool, bool]): Whether to mirror around the (x,y) axes
        center_point ([float, float]): (x,y) point to mirror around

    Returns:
        A new shape list with rotated vertices
    """
    if center_point is None:
        center_point = centroid_for_uncomputed_shapes(shape_list)
    return [mirror_polygon(s, axes, center_point) for s in shape_list]


def _point_equal(p1, p2):
    return p1[0]==p2[0] and p1[1] == p2[1]

def _arr_eq(a1, a2):
    return all(_point_equal(p1,p2) for p1, p2 in zip(a1, a2))

def find_concave_outline(shape_list: List[List[Tuple[float, float]]]):
    """Find the outline of a set of shapes

    Assuming all shapes have edges in common with other shapes where they touch,
    provides a set of vertices for drawing the outline

    Args:
        shape_list (list): A list of list of (x,y) vertices

    Returns:
        A list of ordered (x,y) vertices for drawing an outline
    """
    # Find the most lower-right point
    current_shape = shape_list[0]
    current_pt = current_shape[0]
    test_idx = 1
    next_test_dir = 1
    for s in shape_list:
        for i in range(len(s)):
            p = s[i]
            if ((p[0] < current_pt[0]) or
                    (p[0] == current_pt[0] and p[1] < current_pt[1])):
                # Replace
                current_pt = p
                current_shape = s
                test_idx = (i+1) % len(s)
                next_test_dir = 1
    vertex_list = [current_pt]
    # Keep going until you reach back to the first point
    while not _point_equal(current_shape[test_idx], vertex_list[0]):
        # Iterate through all the shapes to try to find a matching edge
        checking = True
        for s in (s for s in shape_list if not _arr_eq(s, current_shape)):
            if checking:  # Way to break out if match found
                for i in range(len(s)):
                    spt = s[i]
                    if _point_equal(current_pt, spt):
                        spt_after = s[(i+1) % len(s)]
                        spt_before = s[(i-1) % len(s)]
                        test_pt = current_shape[test_idx]
                        if _point_equal(test_pt, spt_after):
                            test_idx = (i-1) % len(s)
                            next_test_dir = -1
                            current_shape = s
                            checking = False
                        elif _point_equal(test_pt, spt_before):
                            test_idx = (i+1) % len(s)
                            next_test_dir = 1
                            current_shape = s
                            checking = False
        # Have you exhausted all shapes?
        if checking:
            current_pt = current_shape[test_idx]
            vertex_list.append(current_pt)
            test_idx += next_test_dir
            test_idx %= len(current_shape)
    return vertex_list
