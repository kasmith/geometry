"""Functions that work on collections of shapes
"""

from __future__ import division, print_function
import numpy as np
from convex import convex_area, convex_centroid

__all__ = ['recenter_polygon', 'centroid_for_shapes',
           'centroid_for_uncomputed_shapes', 'recenter_system',
           'rescale_and_recenter_system']

def recenter_polygon(vertices):
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

def centroid_for_shapes(centroids, areas = None):
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


def centroid_for_uncomputed_shapes(shape_list):
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


def recenter_system(shape_list):
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


def rescale_and_recenter_system(shape_list, total_area):
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
