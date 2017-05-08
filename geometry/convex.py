import numpy as np
from helpers import *

__all__ = ['gift_wrap']

# Simple function to gift wrap a set of vertices to get the convex hull
def gift_wrap(vertices):
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

