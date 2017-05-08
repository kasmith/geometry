import numpy as np
from helpers import *

__all__ = ['ear_clip', 'ear_clip_with_holes']


# Uses ear clipping to perform triangulation of a non-convex vertex set
def ear_clip(vertices):
    assert len(vertices) >= 3, "Requires at least 3 points to triangulate"
    is_clockwise = check_clockwise(vertices)
    trilist = []
    vl = map(np.array, vertices) # Copy and make into numpy arrays (for easy calculation)
    # Go through and keep clipping ears off
    while len(vl) > 3:
        i = 0
        running = True
        lnvl = len(vl)
        while running:
            p1 = vl[i]
            p2 = vl[(i+1) % lnvl]
            p3 = vl[(i+2) % lnvl]
            # Do these points make an ear?
            # First, the p1-p3 line must be inside the polygon
            is_left = np.cross(p3 - p1, p2 - p1) > 0
            if (is_left and is_clockwise) or (not is_left and not is_clockwise):
                # Next, are no other points inside the triangle
                is_ear = True
                j = 0
                while j < len(vl) and is_ear:
                    if j != i and j != ((i+1) % lnvl) and j != ((i+2) % lnvl):
                        # Check if there are points inside -- if so, it's not an ear
                        if point_in_poly(vl[j], [p1,p2,p3]):
                            is_ear = False
                        # Check if there are points on the vertices -- if so, it's not an ear (unless it's a copy)
                        if point_on_line(vl[j], p1, p2) or point_on_line(vl[j], p2, p3) or point_on_line(vl[j], p1, p3):
                            if not (all(vl[j] == p1) or all(vl[j]==p2) or all(vl[j]==p3)):
                                is_ear = False
                    j += 1
                if(is_ear):
                    # Pop that triangle onto the list, remove the ear point
                    trilist.append([p1,p2,p3])
                    del vl[(i+1) % lnvl]
                    running = False
            i += 1
    trilist.append(vl)
    return trilist


# Does ear clipping with holes -- returns the triangle list from clipping and the new vertex list
def ear_clip_with_holes(outer_shell, hole_list):
    # Copy each and make into numpy arrays
    outer_shell = map(np.array, outer_shell)
    hole_list = [map(np.array, hole) for hole in hole_list]
    # Make sure that holes are reverse wound from outer_shell
    # So outer_shell goes ccw, holes go cw
    if check_clockwise(outer_shell):
        outer_shell.reverse()
    for hole in hole_list:
        if check_counterclockwise(hole):
            hole.reverse()
        # Find mutually visible points
        oidx, hidx, opt, hpt = find_mutually_visible(outer_shell, hole)
        # Edge case - the cut points on each hull and next are on the same line -- skip that point
        if opt[1] == hpt[1] == outer_shell[(oidx + 1) % len(outer_shell)][1]:
            outer_shell = outer_shell[:(oidx+1)] + hole[hidx:] + hole[:(hidx+1)] + outer_shell[(oidx+1):]
        # Cut the outer shell list and insert the hole vertices
        else:
            outer_shell = outer_shell[:(oidx+1)] + hole[hidx:] + hole[:(hidx+1)] + outer_shell[oidx:]
    # Now do ear clipping on the full structure
    return ear_clip(outer_shell), outer_shell
