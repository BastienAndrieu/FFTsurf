GFORTRAN module version '10' created from src/modules/mod_bernstein.f90
MD5:f35fcdb5ccdc48297b286b850f459ea7 -- If you edit this, you'll get what you deserve.

(() () () () () () () () () () () () () () () () () () () () () () ()
() () () ())

()

(('cross' 'mod_math' 2) ('print_mat' 'mod_math' 3 4) (
'ptr_bernstein_series1' 'mod_bernstein' 5) ('ptr_bernstein_series2'
'mod_bernstein' 6) ('type_bernstein_series1' 'mod_bernstein' 7) (
'type_bernstein_series2' 'mod_bernstein' 8) ('type_matrix' 'mod_math' 9))

()

()

()

(5 'Ptr_bernstein_series1' 'mod_bernstein' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((10 'ptr' (DERIVED 7 0 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0
90496578)
6 'Ptr_bernstein_series2' 'mod_bernstein' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((11 'ptr' (DERIVED 8 0 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0
90496579)
7 'Type_bernstein_series1' 'mod_bernstein' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOC_COMP) (UNKNOWN 0 0 0 0 UNKNOWN ())
0 0 () () 0 ((12 'degr' (INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0')) (13 'dim' (INTEGER 4 0 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0')) (14
'coef' (REAL 8 0 0 0 REAL ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 90841458)
8 'Type_bernstein_series2' 'mod_bernstein' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOC_COMP) (UNKNOWN 0 0 0 0 UNKNOWN ())
0 0 () () 0 ((15 'degr' (INTEGER 4 0 0 0 INTEGER ()) (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 0
INTEGER ()) 0 '2')) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 0 INTEGER
()) 0 '0')) (16 'dim' (INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0')) (17 'coef' (REAL 8 0 0 0
REAL ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()))
PUBLIC (() () () ()) () 0 0 90841459)
9 'Type_matrix' 'mod_math' '' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 ALLOC_COMP) (UNKNOWN 0 0 0 0 UNKNOWN ()) 0 0 () () 0
((18 'mat' (REAL 8 0 0 0 REAL ()) (2 0 DEFERRED () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE
DIMENSION) UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 42503820)
19 'ab2n1p1' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 FUNCTION IMPLICIT_PURE) (REAL 8 0 0 0 REAL ()) 20 0 (21 22
23) () 19 () () () 0 0)
24 'berndiff' 'mod_bernstein' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE) (UNKNOWN 0 0 0 0
UNKNOWN ()) 25 0 (26 27 28 29) () 0 () () () 0 0)
30 'berndiff1' 'mod_bernstein' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ()) 31
0 (32 33) () 0 () () () 0 0)
34 'berndiff2' 'mod_bernstein' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
0 UNKNOWN ()) 35 0 (36 37 38) () 0 () () () 0 0)
2 'cross_1' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 DIMENSION FUNCTION ALWAYS_EXPLICIT) (REAL 8 0 0 0 REAL ())
39 0 (40 41) (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1')
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '3')) 42 () () () 0 0)
43 'de_casteljau' 'mod_bernstein' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE ALWAYS_EXPLICIT) (
UNKNOWN 0 0 0 0 UNKNOWN ()) 44 0 (45 46 47 48 49 50 51) () 0 () () () 0
0)
52 'diff_angle' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE) (REAL 8 0 0 0 REAL ()) 53 0 (
54 55) () 52 () () () 0 0)
56 'dnchoosek' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION) (REAL 8 0 0 0 REAL ()) 57 0 (58 59) () 56 ()
() () 0 0)
60 'epsilon' '(intrinsic)' '' 1 ((PROCEDURE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 FUNCTION) (UNKNOWN 0 0 0 0 UNKNOWN ()) 0 0 () () 60
() () () 0 0)
61 'factln' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 FUNCTION) (REAL 8 0 0 0 REAL ()) 62 0 (63) () 61 () () () 0
0)
64 'gammaln' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 FUNCTION) (REAL 8 0 0 0 REAL ()) 65 0 (66) () 64 () () () 0
0)
67 'gauss_elim' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 0 UNKNOWN ())
68 0 (69 70) () 0 () () () 0 0)
71 'householder_vector' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE ALWAYS_EXPLICIT) (
UNKNOWN 0 0 0 0 UNKNOWN ()) 72 0 (73 74 75 76) () 0 () () () 0 0)
77 'identity_matrix' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 DIMENSION FUNCTION IMPLICIT_PURE
ALWAYS_EXPLICIT) (REAL 8 0 0 0 REAL ()) 78 0 (79) (2 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (VARIABLE (INTEGER 4 0 0 0
INTEGER ()) 0 79 ()) (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (
VARIABLE (INTEGER 4 0 0 0 INTEGER ()) 0 79 ())) 77 () () () 0 0)
80 'is_in_interval' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE) (LOGICAL 4 0 0 0
LOGICAL ()) 81 0 (82 83 84) () 80 () () () 0 0)
85 'linsolve_qr' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ()) 86 0 (87 88 89
90 91 92) () 0 () () () 0 0)
93 'linsolve_qr_rank_deficient' 'mod_math' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0
UNKNOWN ()) 94 0 (95 96 97 98 99 100) () 0 () () () 0 0)
101 'linspace' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 DIMENSION FUNCTION ALWAYS_EXPLICIT) (REAL 8 0 0 0 REAL
()) 102 0 (103 104 105) (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER
()) 0 '1') (VARIABLE (INTEGER 4 0 0 0 INTEGER ()) 0 105 ())) 101 () () ()
0 0)
106 'matheps' 'mod_math' '' 1 ((PARAMETER UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN IMPLICIT-SAVE 0 0) (REAL 8 0 0 0 REAL ()) 0 0 () (CONSTANT (
REAL 8 0 0 0 REAL ()) 0 '0.a0000000000000@-12') () 0 () () () 0 0)
107 'mathpi' 'mod_math' '' 1 ((PARAMETER UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN IMPLICIT-SAVE 0 0) (REAL 8 0 0 0 REAL ()) 0 0 () (CONSTANT (
REAL 8 0 0 0 REAL ()) 0 '0.3243f6a8885a30@1') () 0 () () () 0 0)
108 'mathpr' 'mod_math' '' 1 ((PARAMETER UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN IMPLICIT-SAVE 0 0) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '8') () 0 () () () 0 0)
109 'mean_angle' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE ALWAYS_EXPLICIT) (REAL 8 0 0 0
REAL ()) 110 0 (111) () 109 () () () 0 0)
112 'mod_bernstein' 'mod_bernstein' '' 1 ((MODULE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) (UNKNOWN 0 0 0 0 UNKNOWN ()) 0 0 () ()
0 () () () 0 0)
113 'mod_math' 'mod_math' '' 1 ((MODULE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) (UNKNOWN 0 0 0 0 UNKNOWN ()) 0 0 () () 0 () () () 0
0)
114 'myfact' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 SUBROUTINE RECURSIVE) (UNKNOWN 0 0 0 0 UNKNOWN ()) 115 0 (
116 117) () 0 () () () 0 0)
118 'myfactorial' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION) (INTEGER 8 0 0 0 INTEGER ()) 119 0 (120) ()
121 () () () 0 0)
122 'n1p12ab' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE) (REAL 8 0 0 0 REAL ()) 123 0 (
124 125 126) () 122 () () () 0 0)
127 'nd_box_constraint' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE ALWAYS_EXPLICIT) (
UNKNOWN 0 0 0 0 UNKNOWN ()) 128 0 (129 130 131 132 133) () 0 () () () 0
0)
134 'norm2_1' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION ALWAYS_EXPLICIT) (REAL 8 0 0 0 REAL ()) 135 0
(136) () 137 () () () 0 0)
138 'norm2_n' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 DIMENSION FUNCTION ALWAYS_EXPLICIT) (REAL 8 0 0 0 REAL
()) 139 0 (140) (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1')
(FUNCTION (INTEGER 4 0 0 0 INTEGER ()) 0 141 (('' (VARIABLE (REAL 8 0 0
0 REAL ()) 2 140 ((ARRAY (FULL 2 2 2))))) ('' (CONSTANT (INTEGER 4 0 0 0
INTEGER ()) 0 '2')) ('' ())) '' 0 'size')) 142 () () () 0 0)
143 'null' '(intrinsic)' '' 1 ((PROCEDURE UNKNOWN-INTENT INTRINSIC-PROC
UNKNOWN UNKNOWN 0 0 FUNCTION) (UNKNOWN 0 0 0 0 UNKNOWN ()) 0 0 () () 0 ()
() () 0 0)
144 'outer_product' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 DIMENSION FUNCTION ALWAYS_EXPLICIT) (REAL 8
0 0 0 REAL ()) 145 0 (146 147) (2 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 0
INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0 INTEGER ()) 0 148 (('' (
VARIABLE (REAL 8 0 0 0 REAL ()) 1 146 ((ARRAY (FULL 1 2))))) ('' ()) (''
())) '' 0 'size') (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (
FUNCTION (INTEGER 4 0 0 0 INTEGER ()) 0 148 (('' (VARIABLE (REAL 8 0 0 0
REAL ()) 1 147 ((ARRAY (FULL 1 2))))) ('' ()) ('' ())) '' 0 'size')) 144
() () () 0 0)
3 'print_mat_d' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 0 UNKNOWN ())
149 0 (150) () 0 () () () 0 0)
4 'print_mat_r' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 0 UNKNOWN ())
151 0 (152) () 0 () () () 0 0)
153 'ptr_bernstein_series1' 'mod_bernstein' '' 1 ((PROCEDURE
UNKNOWN-INTENT UNKNOWN-PROC DECL UNKNOWN 0 0 FUNCTION GENERIC) (UNKNOWN
0 0 0 0 UNKNOWN ()) 0 0 () () 0 () () () 0 0)
154 'ptr_bernstein_series2' 'mod_bernstein' '' 1 ((PROCEDURE
UNKNOWN-INTENT UNKNOWN-PROC DECL UNKNOWN 0 0 FUNCTION GENERIC) (UNKNOWN
0 0 0 0 UNKNOWN ()) 0 0 () () 0 () () () 0 0)
155 'qr_colpiv' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 0 UNKNOWN ())
156 0 (157 158 159 160 161) () 0 () () () 0 0)
162 'reset_bernstein_series1' 'mod_bernstein' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE) (
UNKNOWN 0 0 0 0 UNKNOWN ()) 163 0 (164 165 166) () 0 () () () 0 0)
167 'reset_bernstein_series2' 'mod_bernstein' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE) (
UNKNOWN 0 0 0 0 UNKNOWN ()) 168 0 (169 170 171) () 0 () () () 0 0)
172 'selected_real_kind' '(intrinsic)' '' 1 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (UNKNOWN 0 0 0 0 UNKNOWN ())
0 0 () () 172 () () () 0 0)
173 'solve_2x2' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE) (UNKNOWN 0 0 0 0 UNKNOWN ())
174 0 (175 176 177 178 179 180 181) () 0 () () () 0 0)
182 'solve_nxn' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 0 UNKNOWN ())
183 0 (184 185 186 187) () 0 () () () 0 0)
188 'solve_up_tri' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE) (UNKNOWN 0 0 0 0
UNKNOWN ()) 189 0 (190 191 192 193 194) () 0 () () () 0 0)
195 'subdiv2' 'mod_bernstein' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
0 UNKNOWN ()) 196 0 (197 198 199 200 201 202) () 0 () () () 0 0)
203 'subdiv2_along_u' 'mod_bernstein' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ())
204 0 (205 206 207 208) () 0 () () () 0 0)
209 'subdiv2_along_v' 'mod_bernstein' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ())
210 0 (211 212 213 214) () 0 () () () 0 0)
215 'swap_rows' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE ALWAYS_EXPLICIT) (UNKNOWN 0 0
0 0 UNKNOWN ()) 216 0 (217 218 219) () 0 () () () 0 0)
220 'triple_product' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (REAL 8 0 0 0 REAL ()) 221 0 (
222 223 224) () 225 () () () 0 0)
226 'type_bernstein_series1' 'mod_bernstein' '' 1 ((PROCEDURE
UNKNOWN-INTENT UNKNOWN-PROC DECL UNKNOWN 0 0 FUNCTION GENERIC) (UNKNOWN
0 0 0 0 UNKNOWN ()) 0 0 () () 0 () () () 0 0)
227 'type_bernstein_series2' 'mod_bernstein' '' 1 ((PROCEDURE
UNKNOWN-INTENT UNKNOWN-PROC DECL UNKNOWN 0 0 FUNCTION GENERIC) (UNKNOWN
0 0 0 0 UNKNOWN ()) 0 0 () () 0 () () () 0 0)
228 'type_matrix' 'mod_math' '' 1 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC DECL UNKNOWN 0 0 FUNCTION GENERIC) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 () () () 0 0)
229 'write_bernstein_series1' 'mod_bernstein' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0
UNKNOWN ()) 230 0 (231 232) () 0 () () () 0 0)
233 'write_bernstein_series2' 'mod_bernstein' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0
UNKNOWN ()) 234 0 (235 236) () 0 () () () 0 0)
21 'x' '' '' 20 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
22 'a' '' '' 20 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
23 'b' '' '' 20 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
26 'b' '' '' 25 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (2 0 EXPLICIT (CONSTANT (INTEGER 4
0 0 0 INTEGER ()) 0 '1') (OP (INTEGER 4 0 0 0 INTEGER ()) 0 PLUS (
VARIABLE (INTEGER 4 0 0 0 INTEGER ()) 0 28 ()) (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1')) (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (
VARIABLE (INTEGER 4 0 0 0 INTEGER ()) 0 29 ())) 0 () () () 0 0)
27 'd' '' '' 25 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (2 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0 INTEGER ())
0 237 (('' (VARIABLE (INTEGER 4 0 0 0 INTEGER ()) 0 28 ())) ('' (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1'))) '__max_i4' 0 'max') (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (VARIABLE (INTEGER 4 0 0 0
INTEGER ()) 0 29 ())) 0 () () () 0 0)
28 'degr' '' '' 25 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
29 'dim' '' '' 25 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
32 'b' '' '' 31 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
DERIVED 7 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
33 'd' '' '' 31 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
DERIVED 7 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
36 'b' '' '' 35 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
DERIVED 8 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
37 'du' '' '' 35 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(DERIVED 8 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
38 'dv' '' '' 35 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (DERIVED 8 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
40 'u' '' '' 39 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 0
INTEGER ()) 0 '3')) 0 () () () 0 0)
41 'v' '' '' 39 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 0
INTEGER ()) 0 '3')) 0 () () () 0 0)
42 'w' '' '' 39 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION RESULT ALWAYS_EXPLICIT) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '3')) 0 () () () 0 0)
45 'b' '' '' 44 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (2 0 EXPLICIT (CONSTANT (INTEGER 4
0 0 0 INTEGER ()) 0 '1') (OP (INTEGER 4 0 0 0 INTEGER ()) 0 PLUS (
VARIABLE (INTEGER 4 0 0 0 INTEGER ()) 0 47 ()) (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1')) (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (
VARIABLE (INTEGER 4 0 0 0 INTEGER ()) 0 48 ())) 0 () () () 0 0)
46 't' '' '' 44 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
47 'degr' '' '' 44 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
48 'dim' '' '' 44 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
49 'f' '' '' 44 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION OPTIONAL DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (VARIABLE (INTEGER 4 0 0 0
INTEGER ()) 0 48 ())) 0 () () () 0 0)
50 'bl' '' '' 44 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION OPTIONAL DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (2 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (OP (INTEGER 4 0 0 0
INTEGER ()) 0 PLUS (VARIABLE (INTEGER 4 0 0 0 INTEGER ()) 0 47 ()) (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1')) (CONSTANT (INTEGER 4 0 0 0
INTEGER ()) 0 '1') (VARIABLE (INTEGER 4 0 0 0 INTEGER ()) 0 48 ())) 0 ()
() () 0 0)
51 'br' '' '' 44 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION OPTIONAL DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (2 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (OP (INTEGER 4 0 0 0
INTEGER ()) 0 PLUS (VARIABLE (INTEGER 4 0 0 0 INTEGER ()) 0 47 ()) (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1')) (CONSTANT (INTEGER 4 0 0 0
INTEGER ()) 0 '1') (VARIABLE (INTEGER 4 0 0 0 INTEGER ()) 0 48 ())) 0 ()
() () 0 0)
54 'a1' '' '' 53 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
55 'a2' '' '' 53 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
58 'n' '' '' 57 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
59 'k' '' '' 57 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
63 'n' '' '' 62 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
66 'xx' '' '' 65 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
69 'a' '' '' 68 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
70 'singular' '' '' 68 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0 0)
73 'x' '' '' 72 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
74 'i' '' '' 72 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
75 'j' '' '' 72 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
76 'v' '' '' 72 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0 INTEGER ())
0 238 (('' (VARIABLE (REAL 8 0 0 0 REAL ()) 1 73 ((ARRAY (FULL 1 2)))))
('' ()) ('' ())) '' 0 'size')) 0 () () () 0 0)
79 'n' '' '' 78 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
82 'x' '' '' 81 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
83 'a' '' '' 81 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
84 'b' '' '' 81 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
87 'x' '' '' 86 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (VARIABLE (INTEGER 4 0 0 0
INTEGER ()) 0 91 ())) 0 () () () 0 0)
88 'a' '' '' 86 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (2 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (VARIABLE (INTEGER 4 0 0 0
INTEGER ()) 0 90 ()) (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (
VARIABLE (INTEGER 4 0 0 0 INTEGER ()) 0 91 ())) 0 () () () 0 0)
89 'b' '' '' 86 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (VARIABLE (INTEGER 4 0 0 0
INTEGER ()) 0 90 ())) 0 () () () 0 0)
90 'm' '' '' 86 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
91 'n' '' '' 86 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
92 'singular' '' '' 86 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0 0)
95 'y' '' '' 94 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') (VARIABLE (INTEGER 4 0 0 0 INTEGER ())
0 99 ())) 0 () () () 0 0)
96 'r' '' '' 94 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (2 0 EXPLICIT (CONSTANT (INTEGER 4
0 0 0 INTEGER ()) 0 '1') (VARIABLE (INTEGER 4 0 0 0 INTEGER ()) 0 98 ())
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (VARIABLE (INTEGER 4 0 0 0
INTEGER ()) 0 99 ())) 0 () () () 0 0)
97 'c' '' '' 94 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (VARIABLE (INTEGER 4 0 0 0
INTEGER ()) 0 100 ())) 0 () () () 0 0)
98 'm' '' '' 94 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
99 'n' '' '' 94 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
100 'rank' '' '' 94 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
103 'a' '' '' 102 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
104 'b' '' '' 102 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
105 'n' '' '' 102 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
111 'ang' '' '' 110 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
116 'n' '' '' 115 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
117 'f' '' '' 115 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 8 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
120 'n' '' '' 119 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
121 'f' '' '' 119 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 RESULT) (INTEGER 8 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
124 'x' '' '' 123 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
125 'a' '' '' 123 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
126 'b' '' '' 123 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
129 'x' '' '' 128 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
130 'lowerb' '' '' 128 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
131 'upperb' '' '' 128 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
132 'dx' '' '' 128 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
133 'lambda' '' '' 128 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
136 'u' '' '' 135 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
137 'm' '' '' 135 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 RESULT ALWAYS_EXPLICIT) (REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () ()
0 0)
140 'u' '' '' 139 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
141 'size' '(intrinsic)' '' 139 ((PROCEDURE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 FUNCTION) (REAL 4 0 0 0 REAL ()) 0 0 () () 141 () ()
() 0 0)
142 'm' '' '' 139 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (
INTEGER 4 0 0 0 INTEGER ()) 0 141 (('' (VARIABLE (REAL 8 0 0 0 REAL ())
2 140 ((ARRAY (FULL 2 2 2))))) ('' (CONSTANT (INTEGER 4 0 0 0 INTEGER ())
0 '2')) ('' ())) '' 0 'size')) 0 () () () 0 0)
146 'u' '' '' 145 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
147 'v' '' '' 145 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
148 'size' '(intrinsic)' '' 145 ((PROCEDURE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 FUNCTION) (REAL 4 0 0 0 REAL ()) 0 0 () () 148 () ()
() 0 0)
150 'a' '' '' 149 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
152 'a' '' '' 151 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 4 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
157 'a' '' '' 156 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
158 'q' '' '' 156 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (2 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0 INTEGER ())
0 239 (('' (VARIABLE (REAL 8 0 0 0 REAL ()) 2 157 ((ARRAY (FULL 2 2 2)))))
('' (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1')) ('' ())) '' 0 'size')
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0
INTEGER ()) 0 239 (('' (VARIABLE (REAL 8 0 0 0 REAL ()) 2 157 ((ARRAY (
FULL 2 2 2))))) ('' (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1')) ('' ()))
'' 0 'size')) 0 () () () 0 0)
159 'r' '' '' 156 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (2 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0 INTEGER ())
0 239 (('' (VARIABLE (REAL 8 0 0 0 REAL ()) 2 157 ((ARRAY (FULL 2 2 2)))))
('' (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1')) ('' ())) '' 0 'size')
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0
INTEGER ()) 0 239 (('' (VARIABLE (REAL 8 0 0 0 REAL ()) 2 157 ((ARRAY (
FULL 2 2 2))))) ('' (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '2')) ('' ()))
'' 0 'size')) 0 () () () 0 0)
160 'p' '' '' 156 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (2 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0 INTEGER ())
0 239 (('' (VARIABLE (REAL 8 0 0 0 REAL ()) 2 157 ((ARRAY (FULL 2 2 2)))))
('' (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '2')) ('' ())) '' 0 'size')
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0
INTEGER ()) 0 239 (('' (VARIABLE (REAL 8 0 0 0 REAL ()) 2 157 ((ARRAY (
FULL 2 2 2))))) ('' (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '2')) ('' ()))
'' 0 'size')) 0 () () () 0 0)
161 'rank' '' '' 156 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
164 'b' '' '' 163 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(DERIVED 7 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
165 'degr' '' '' 163 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
166 'dim' '' '' 163 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
169 'b' '' '' 168 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(DERIVED 8 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
170 'degr' '' '' 168 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 0
INTEGER ()) 0 '2')) 0 () () () 0 0)
171 'dim' '' '' 168 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
175 'x' '' '' 174 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 0 INTEGER ())
0 '2')) 0 () () () 0 0)
176 'a11' '' '' 174 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
177 'a12' '' '' 174 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
178 'a21' '' '' 174 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
179 'a22' '' '' 174 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
180 'b' '' '' 174 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 0 INTEGER ())
0 '2')) 0 () () () 0 0)
181 'singular' '' '' 174 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0 0)
184 'x' '' '' 183 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0 INTEGER ())
0 240 (('' (VARIABLE (REAL 8 0 0 0 REAL ()) 2 185 ((ARRAY (FULL 2 2 2)))))
('' (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '2')) ('' ())) '' 0 'size'))
0 () () () 0 0)
185 'a' '' '' 183 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
186 'b' '' '' 183 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
187 'singular' '' '' 183 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0 0)
190 'x' '' '' 189 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') (VARIABLE (INTEGER 4 0 0 0 INTEGER ())
0 194 ())) 0 () () () 0 0)
191 'u' '' '' 189 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (2 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') (VARIABLE (INTEGER 4 0 0 0 INTEGER ())
0 193 ()) (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (VARIABLE (
INTEGER 4 0 0 0 INTEGER ()) 0 194 ())) 0 () () () 0 0)
192 'b' '' '' 189 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') (VARIABLE (INTEGER 4 0 0 0 INTEGER ())
0 193 ())) 0 () () () 0 0)
193 'm' '' '' 189 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
194 'n' '' '' 189 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
197 'b' '' '' 196 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(DERIVED 8 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
198 'uv' '' '' 196 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 0 INTEGER ())
0 '2')) 0 () () () 0 0)
199 'bsw' '' '' 196 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (DERIVED 8 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
200 'bse' '' '' 196 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (DERIVED 8 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
201 'bnw' '' '' 196 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (DERIVED 8 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
202 'bne' '' '' 196 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (DERIVED 8 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
205 'b' '' '' 204 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(DERIVED 8 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
206 'u' '' '' 204 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
207 'bw' '' '' 204 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(DERIVED 8 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
208 'be' '' '' 204 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(DERIVED 8 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
211 'b' '' '' 210 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(DERIVED 8 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
212 'v' '' '' 210 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
213 'bs' '' '' 210 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(DERIVED 8 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
214 'bn' '' '' 210 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(DERIVED 8 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
217 'a' '' '' 216 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
218 'i' '' '' 216 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
219 'j' '' '' 216 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
222 'a' '' '' 221 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 0
INTEGER ()) 0 '3')) 0 () () () 0 0)
223 'b' '' '' 221 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 0
INTEGER ()) 0 '3')) 0 () () () 0 0)
224 'c' '' '' 221 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 0
INTEGER ()) 0 '3')) 0 () () () 0 0)
225 'd' '' '' 221 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 RESULT) (REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
231 'b' '' '' 230 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(DERIVED 7 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
232 'filename' '' '' 230 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
235 'b' '' '' 234 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(DERIVED 8 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
236 'filename' '' '' 234 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
237 'max' '(intrinsic)' '' 25 ((PROCEDURE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 FUNCTION) (UNKNOWN 0 0 0 0 UNKNOWN ()) 0 0 () () 237
() () () 0 0)
238 'size' '(intrinsic)' '' 72 ((PROCEDURE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 FUNCTION) (REAL 4 0 0 0 REAL ()) 0 0 () () 238 () ()
() 0 0)
239 'size' '(intrinsic)' '' 156 ((PROCEDURE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 FUNCTION) (REAL 4 0 0 0 REAL ()) 0 0 () () 239 () ()
() 0 0)
240 'size' '(intrinsic)' '' 183 ((PROCEDURE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 FUNCTION) (REAL 4 0 0 0 REAL ()) 0 0 () () 240 () ()
() 0 0)
)

('Ptr_bernstein_series1' 0 5 'Ptr_bernstein_series2' 0 6
'Type_bernstein_series1' 0 7 'Type_bernstein_series2' 0 8 'Type_matrix'
0 9 'ab2n1p1' 0 19 'berndiff' 0 24 'berndiff1' 0 30 'berndiff2' 0 34
'cross_1' 0 2 'de_casteljau' 0 43 'diff_angle' 0 52 'dnchoosek' 0 56
'epsilon' 0 60 'factln' 0 61 'gammaln' 0 64 'gauss_elim' 0 67
'householder_vector' 0 71 'identity_matrix' 0 77 'is_in_interval' 0 80
'linsolve_qr' 0 85 'linsolve_qr_rank_deficient' 0 93 'linspace' 0 101
'matheps' 0 106 'mathpi' 0 107 'mathpr' 0 108 'mean_angle' 0 109
'mod_bernstein' 0 112 'mod_math' 0 113 'myfact' 0 114 'myfactorial' 0
118 'n1p12ab' 0 122 'nd_box_constraint' 0 127 'norm2_1' 0 134 'norm2_n'
0 138 'null' 0 143 'outer_product' 0 144 'print_mat_d' 0 3 'print_mat_r'
0 4 'ptr_bernstein_series1' 0 153 'ptr_bernstein_series2' 0 154
'qr_colpiv' 0 155 'reset_bernstein_series1' 0 162
'reset_bernstein_series2' 0 167 'selected_real_kind' 0 172 'solve_2x2' 0
173 'solve_nxn' 0 182 'solve_up_tri' 0 188 'subdiv2' 0 195
'subdiv2_along_u' 0 203 'subdiv2_along_v' 0 209 'swap_rows' 0 215
'triple_product' 0 220 'type_bernstein_series1' 0 226
'type_bernstein_series2' 0 227 'type_matrix' 0 228
'write_bernstein_series1' 0 229 'write_bernstein_series2' 0 233)
