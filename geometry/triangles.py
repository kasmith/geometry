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


# Does ear clipping with holes
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
        oidx, hidx, _, _ = find_mutually_visible(outer_shell, hole)
        # Cut the outer shell list and insert the hole vertices
        outer_shell = outer_shell[:(oidx+1)] + hole[hidx:] + hole[:(hidx+1)] + outer_shell[oidx:]
    # Now do ear clipping on the full structure
    return ear_clip(outer_shell)

if __name__ == '__main__':
    vlist = map(np.array, [(79, 78), (493, 78), (493, 114), (678, 114), (678, 163), (928, 163), (928, 269), (934, 269), (934, 491),
               (974, 491), (974, 596), (909, 596), (909, 560), (663, 560), (663, 491), (699, 491), (699, 401),
               (656, 401), (656, 469), (588, 469), (588, 401), (551, 401), (551, 365), (438, 365), (438, 431),
               (153, 431), (153, 258), (79, 258)])

    #vlist = map(np.array, [(0,0), (4,0), (4,4), (3,4), (3,2), (1,2), (1,3), (2,3), (2,4), (0,4)])
    print ear_clip(vlist)