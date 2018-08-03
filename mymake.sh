python src/intersection/makemodule.py

# modules
for f in src/modules/*.f90
do
    gfortran -g -Wall -Wextra -fbounds-check -Iobj/ -Jobj/ -c $f -o obj/$(basename $f .f90).o
    echo "Compiling $(basename $f)"
done

echo ""

## TESTS

#echo "...union_intersection_arrays.out"
#gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/union_intersection_arrays.f90 -o test/obj/union_intersection_arrays.o
#gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/union_intersection_arrays.out test/obj/union_intersection_arrays.o obj/mod_util.o


#echo "...tree.out"
#gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/tree.f90 -o test/obj/tree.o
#gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/tree.out test/obj/tree.o obj/mod_util.o obj/mod_math.o obj/mod_obb.o obj/mod_chebyshev2.o obj/mod_bernstein2.o obj/mod_polynomial.o obj/mod_regiontree.o -Lsrc/dfftpack -ldfftpack


#echo "...chebyshev2.out"
#gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/chebyshev2.f90 -o test/obj/chebyshev2.o
#gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/chebyshev2.out test/obj/chebyshev2.o obj/mod_util.o obj/mod_math.o obj/mod_chebyshev.o obj/mod_bernstein.o obj/mod_polynomial.o -Lsrc/dfftpack -ldfftpack


#echo "...bernstein2.out"
#gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/bernstein2.f90 -o test/obj/bernstein2.o
#gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/bernstein2.out test/obj/bernstein2.o obj/mod_util.o obj/mod_math.o obj/mod_chebyshev.o obj/mod_bernstein.o obj/mod_polynomial.o -Lsrc/dfftpack -ldfftpack


#echo "...newton_curve_surface_singular.out"
#gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/newton_curve_surface_singular.f90 -o test/obj/newton_curve_surface_singular.o
#gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/newton_curve_surface_singular.out test/obj/newton_curve_surface_singular.o obj/mod_util.o obj/mod_math.o obj/mod_linalg.o obj/mod_chebyshev.o obj/mod_bernstein.o obj/mod_polynomial.o obj/mod_diffgeom.o obj/mod_tolerances.o -Lsrc/dfftpack -ldfftpack


#echo "...separation.out"
#gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/separation.f90 -o test/obj/separation.o
#gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/separation.out test/obj/separation.o obj/mod_util.o obj/mod_math.o obj/mod_linprog.o obj/mod_geometry.o obj/mod_separation.o 


#echo "...linearprogramming.out"
#gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/linearprogramming.f90 -o test/obj/linearprogramming.o
#gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/linearprogramming.out test/obj/linearprogramming.o obj/mod_util.o obj/mod_math.o obj/mod_linprog.o


#echo "...convex_hull.out"
#gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/convex_hull.f90 -o test/obj/convex_hull.o
#gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/convex_hull.out test/obj/convex_hull.o obj/mod_util.o obj/mod_math.o obj/mod_tolerances.o obj/mod_geometry.o 


#echo "...linearalgebra.out"
#gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/linearalgebra.f90 -o test/obj/linearalgebra.o
#gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/linearalgebra.out test/obj/linearalgebra.o obj/mod_util.o obj/mod_math.o obj/mod_linalg.o


echo "...dev_intersection.out"
gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/dev_intersection.f90 -o test/obj/dev_intersection.o
gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/dev_intersection.out test/obj/dev_intersection.o obj/mod_util.o obj/mod_math.o obj/mod_linalg.o obj/mod_chebyshev.o obj/mod_bernstein.o obj/mod_polynomial.o obj/mod_diffgeom.o obj/mod_linprog.o obj/mod_geometry.o obj/mod_separation.o obj/mod_obb.o obj/mod_regiontree.o obj/mod_tolerances.o obj/mod_types_intersection.o obj/mod_intersection.o -Lsrc/dfftpack -ldfftpack


#echo "...polylines.out"
#gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/polylines.f90 -o test/obj/polylines.o
#gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/polylines.out test/obj/polylines.o obj/mod_util.o obj/mod_math.o

echo "...jouke.out"
gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/jouke.f90 -o test/obj/jouke.o
gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/jouke.out test/obj/jouke.o obj/mod_util.o obj/mod_math.o obj/mod_linalg.o obj/mod_chebyshev.o obj/mod_bernstein.o obj/mod_polynomial.o obj/mod_diffgeom.o obj/mod_linprog.o obj/mod_geometry.o obj/mod_separation.o obj/mod_obb.o obj/mod_regiontree.o obj/mod_tolerances.o obj/mod_types_intersection.o obj/mod_intersection.o -Lsrc/dfftpack -ldfftpack


echo "...diffgeom.out"
gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/diffgeom.f90 -o test/obj/diffgeom.o
gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/diffgeom.out test/obj/diffgeom.o obj/mod_util.o obj/mod_math.o obj/mod_linalg.o obj/mod_chebyshev.o obj/mod_bernstein.o obj/mod_polynomial.o obj/mod_diffgeom.o obj/mod_linprog.o obj/mod_geometry.o obj/mod_separation.o obj/mod_obb.o obj/mod_regiontree.o obj/mod_tolerances.o obj/mod_types_intersection.o obj/mod_intersection.o -Lsrc/dfftpack -ldfftpack
