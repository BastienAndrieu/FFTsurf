verbose=true

python src/chebyshev/makemodule.py
python src/bernstein/makemodule.py
#python src/intersection/makemodule.py

# modules
for f in src/modules/*.f90
do
    gfortran -Wall -Wextra -fbounds-check -g -Iobj/ -Jobj/ -c $f -o obj/$(basename $f .f90).o
    if $verbose; then
	echo "Compiling $(basename $f)"
    fi
done


if $verbose; then
    echo ""
fi
# intersections
for f in src/intersection/*.f90
do
    gfortran -Wall -Wextra -fbounds-check -g -Iobj/ -c $f -o obj/intersection/$(basename $f .f90).o
    if $verbose; then
	echo "Compiling $(basename $f)"
    fi
done
if $verbose; then
    echo "Library intersection"
fi
ar cr obj/intersection/libintersection.a obj/intersection/*.o
chmod 777 obj/intersection/libintersection.a

if $verbose; then
    echo ""
fi
# tests
gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/pre_intersection.f90 -o test/obj/pre_intersection.o
gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/pre_intersection.out test/obj/pre_intersection.o obj/mod_util.o  obj/mod_math.o obj/mod_chebyshev.o obj/mod_bernstein.o obj/mod_diffgeom.o obj/mod_obb.o obj/mod_errors_intersection.o obj/mod_linearprogramming.o obj/mod_geometry.o obj/mod_separation.o obj/mod_tolerances.o -Lsrc/dfftpack -ldfftpack -Lobj/intersection -lintersection
#-L/home/bandrieu/FFTsurf/src/dfftpack -ldfftpack -L/home/bandrieu/FFTsurf/obj/intersection -lintersection


#obj/mod_intersection.o 
#-L/stck/bandrieu/FFTsurf/obj/intersection/libintersection.a
#-L/stck/bandrieu/FFTsurf/src/dfftpack/libdfftpack.a


gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/complexite_fcht.f90 -o test/obj/complexite_fcht.o
gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/complexite_fcht.out test/obj/complexite_fcht.o obj/mod_util.o obj/mod_math.o obj/mod_chebyshev.o -Lsrc/dfftpack -ldfftpack

gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/union_intersection_arrays.f90 -o test/obj/union_intersection_arrays.o
gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/union_intersection_arrays.out test/obj/union_intersection_arrays.o obj/mod_util.o

gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/decasteljau.f90 -o test/obj/decasteljau.o
gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/decasteljau.out test/obj/decasteljau.o obj/mod_util.o obj/mod_math.o obj/mod_bernstein.o obj/mod_obb.o


gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/pre_intersection_bezier.f90 -o test/obj/pre_intersection_bezier.o
gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/pre_intersection_bezier.out test/obj/pre_intersection_bezier.o obj/mod_util.o obj/mod_math.o obj/mod_chebyshev.o obj/mod_bernstein.o obj/mod_diffgeom.o obj/mod_obb.o obj/mod_errors_intersection.o obj/mod_linearprogramming.o obj/mod_geometry.o obj/mod_separation.o -Lsrc/dfftpack -ldfftpack -Lobj/intersection -lintersection


gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/subdiv2_bezier.f90 -o test/obj/subdiv2_bezier.o
gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/subdiv2_bezier.out test/obj/subdiv2_bezier.o obj/mod_util.o obj/mod_math.o obj/mod_chebyshev.o obj/mod_bernstein.o obj/mod_obb.o -Lsrc/dfftpack -ldfftpack



gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/dev_intersection_simple_surface.f90 -o test/obj/dev_intersection_simple_surface.o
gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/dev_intersection_simple_surface.out test/obj/dev_intersection_simple_surface.o obj/mod_util.o  obj/mod_math.o obj/mod_chebyshev2.o obj/mod_bernstein2.o obj/mod_polynomial.o obj/mod_diffgeom2.o obj/mod_linearprogramming.o obj/mod_geometry.o obj/mod_separation.o obj/mod_obb.o obj/mod_regiontree.o obj/mod_tolerances.o -Lsrc/dfftpack -ldfftpack

gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/perfo_subdivision.f90 -o test/obj/perfo_subdivision.o
gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/perfo_subdivision.out test/obj/perfo_subdivision.o obj/mod_util.o obj/mod_math.o obj/mod_chebyshev.o obj/mod_bernstein.o -Lsrc/dfftpack -ldfftpack



gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/tree.f90 -o test/obj/tree.o
gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/tree.out test/obj/tree.o obj/mod_util.o obj/mod_math.o obj/mod_obb.o obj/mod_chebyshev2.o obj/mod_bernstein2.o obj/mod_polynomial.o obj/mod_regiontree.o -Lsrc/dfftpack -ldfftpack



gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/chebyshev2.f90 -o test/obj/chebyshev2.o
gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/chebyshev2.out test/obj/chebyshev2.o obj/mod_util.o obj/mod_math.o obj/mod_chebyshev2.o obj/mod_bernstein2.o obj/mod_polynomial.o -Lsrc/dfftpack -ldfftpack


gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/bernstein2.f90 -o test/obj/bernstein2.o
gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/bernstein2.out test/obj/bernstein2.o obj/mod_util.o obj/mod_math.o obj/mod_chebyshev2.o obj/mod_bernstein2.o obj/mod_polynomial.o -Lsrc/dfftpack -ldfftpack


gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/debugobb.f90 -o test/obj/debugobb.o
gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/debugobb.out test/obj/debugobb.o obj/mod_util.o obj/mod_math.o obj/mod_obb.o


gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/dev_intersection_surface_surface.f90 -o test/obj/dev_intersection_surface_surface.o
gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/dev_intersection_surface_surface.out test/obj/dev_intersection_surface_surface.o obj/mod_util.o  obj/mod_math.o obj/mod_chebyshev2.o obj/mod_bernstein2.o obj/mod_polynomial.o obj/mod_diffgeom2.o obj/mod_linearprogramming.o obj/mod_geometry.o obj/mod_separation.o obj/mod_obb.o obj/mod_regiontree.o obj/mod_tolerances.o -Lsrc/dfftpack -ldfftpack


gfortran -Wall -Wextra -fbounds-check -g -c -Iobj/ test/src/newton_curve_surface_singular.f90 -o test/obj/newton_curve_surface_singular.o
gfortran -Wall -Wextra -fbacktrace -fbounds-check -g -o test/newton_curve_surface_singular.out test/obj/newton_curve_surface_singular.o obj/mod_util.o  obj/mod_math.o obj/mod_chebyshev2.o obj/mod_bernstein2.o obj/mod_polynomial.o obj/mod_diffgeom2.o obj/mod_tolerances.o -Lsrc/dfftpack -ldfftpack
