VERBOSE = True

progname = 'dev_intersection'
f = open( '../' + progname + '.f90', 'w' )
f.write('program ' + progname + '\n')

f.writelines( [l for l in open('main.f90').readlines()] )


#f.write('\n\ncontains\n')

listfiles = ['add_intersection_curve.f90',
             'append_vector.f90',
             'characterize_tangential_intersection_point.f90',
             'check_unicity.f90',
             'classify_border_surface_intersection_point.f90',
             'classify_endpoint.f90',
             'diffgeom_intersection_curve.f90',
             'find_collineal_corners.f90',
             'find_collineal_points.f90',
             'inherit_points.f90',
             'intersect_all_surfaces.f90',
             'intersect_border_surface.f90',
             'intersect_curve_surface.f90',
             'newton_curve_surface.f90',
             'intersect_elsewhere.f90',
             'intersect_gaussmaps_elsewhere.f90',
             'intersect_simple_surfaces.f90',
             'intersect_surface_pair.f90',
             'loop_detection_criterion.f90',
             'merge_intersection_data.f90',
             'rearrange_for_separability_test.f90',
             'transfer_intersection_curves.f90']


for filename in listfiles:
    if VERBOSE:
        print '\t * <---', filename
    f.write('\n\n\n\n\n')
    f.writelines( [l for l in open(filename).readlines()] )

f.write('end program ' + progname + '\n')

f.close()
