"""Functions that help with convex hulls
"""

from typing import Tuple, Annotated, Dict
import numpy as np
from .helpers import *

__all__ = ['gift_wrap', 'convex_area', 'convex_centroid']


def gift_wrap(vertices: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    """Simple function to gift wrap a set of vertices to get the convex hull

    Uses the gift-wrapping algorithm

    Args:
        vertices (list): A list of (x,y) points to gift wrap

    Returns:
        The CCW set of vertices forming the convex hull of the input
    """
    # Get the left-most vertex
    pt = vertices[0]
    for v in vertices:
        if v[0] < pt[0]:
            pt = v
    running = True
    hull = []
    i = 0
    while running:
        hull.append(pt)
        ep = vertices[0]
        for v in vertices:
            isleft = point_on_left(v, [pt, ep])
            if all(ep == pt) or isleft:
                ep = v
        i += 1
        pt = ep
        running = any(hull[0] != pt)
    hull.reverse()
    return hull

def convex_area(vertices: List[Tuple[float, float]]) -> float:
    """Returns the area of a convex polygon

    Args:
        vertices (list): list of (x,y) vertices of convex polygon

    Returns:
        The area of the polygon
    """
    off = vertices[0]
    twicearea = 0
    nverts = len(vertices)
    for i in range(nverts):
        j = (i - 1) % nverts
        v1 = vertices[i]
        v2 = vertices[j]
        twicearea += ((v1[0]-off[0]) * (v2[1]-off[1]) -
                      (v2[0]-off[0]) * (v1[1]-off[1]))
    return twicearea / 2.

def convex_centroid(vertices: List[Tuple[float, float]]) -> Tuple[float, float]:
    """Returns the centroid of a convex polygon

    Args:
        vertices (list): list of (x,y) vertices of convex polygon

    Returns:
        The (x,y) position of the center of the polygon
    """
    tsum = 0
    vsum = np.zeros(2)
    arrv = [np.array(v) for v in vertices]
    for i in range(len(arrv)):
        v1 = arrv[i]
        v2 = arrv[(i+1) % len(arrv)]
        cross = np.cross(v1, v2)
        tsum += cross
        vsum += (v1+v2) * cross
    return vsum * (1/(3*tsum))
