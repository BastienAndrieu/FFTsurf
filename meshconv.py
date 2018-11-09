import sys, getopt
import numpy as np
import os

##### OPTIONS #####
# -f faces (mandatory)
# -uv uv coordinates (optional)
# -xyz xyz coords. (optional)
# -poly polynomial (optional)
# -o output (mandatory)

# au moins un des deux arguments uv ou xyz doit etre renseigne
# si f,uv seuls alors maillage 2d plan texture
# si f,xyz seuls alors maillage surface non texture
# si f,uv,poly alors maillage surface texture avec xyz=poly(uv)
# si f,uv,xyz alors maillage surface texture (poly ignore si present)
# comprendre extension fichier output:
#    .obj, .mesh (MFEM ou INRIA), .stl, .vtk
##################


float_fmt = "08.6e"

#############################################################
def read_points( p, filename, dim, label="" ):
    p = []
    f = open( filename, "r" )
    for l in f:
        pi = [float(a) for a in l.split() ]
        if len(pi) != dim:
            print "read_points: "+label+" must be "+str(dim)+"-dimensional"
            sys.exit(2)
        p.append( pi )
    f.close()
    return p
#############################################################
def read_polynomial( filename ):
    f = open(filename,"r")
    m,n,p = [int(a) for a in f.readline().split()]
    coef = np.zeros((m,n,p))
    for k in range(p):
        for j in range(n):
            for i in range(m):
                coef[i,j,k] = float(f.readline())
    
    f.close()
    return coef
#############################################################
def write_obj( xyz, uv, faces, filename ):
    f = open( filename, "w" )
    f.write( filename + "\n\n\n" )
    
    # vertices
    for v in xyz:
        f.write( "v ")
        for x in v:
            f.write( str(x) + " " )
        f.write( "\n" )
    f.write( "\n" )
    
    # texture coordinates
    if len(uv) > 0:
        for vt in uv:
            f.write( "vt ")
            for u in vt:
                f.write( str(u) + " " )
            f.write( "\n" )
        f.write( "\n" )

    # faces
    if len(uv) > 0:
        for fa in faces:
            f.write( "f ")
            for i in fa:
                f.write( str(i) + "/" + str(i) + " " )
            f.write( "\n" )
        f.write( "\n" )
    else:
        for fa in faces:
            f.write( "f ")
            for i in fa:
                f.write( str(i) + " " )
            f.write( "\n" )
        f.write( "\n" )

    f.close()
    return
#############################################################
def write_MFEM_mesh( xyz, faces, filename):
    f = open( filename, "w" )
    f.write( "MFEM mesh v1.0\n\ndimension\n2\n\n" )
    
    # faces
    f.write( "elements\n" + str(len(faces)) + "\n" )
    for fa in faces:
        f.write( "1 " + str(len(fa)-1) + " " )
        for i in fa:
            f.write( str(i-1) + " " )
        f.write( "\n" )
    f.write( "\n" )
    
    # boundary edges
    bed = boundary_edges( faces, len(xyz) )
    #bed = []
    f.write( "boundary\n" + str(len(bed)) + "\n" )
    for e in bed:
        f.write( "1 1 " + str(e[0]-1) + " " + str(e[1]-1) + "\n" )
    f.write( "\n" )
    
    # vertices
    f.write( "vertices\n" + str(len(xyz)) + "\n" )
    f.write( "\n" )
    f.write( "nodes\nFiniteElementSpace\nFiniteElementCollection: Linear\n" )
    f.write( "VDim: " + str(len(xyz[0])) + "\n" )
    f.write( "Ordering: 1\n" )
    for v in xyz:
        for x in v:
            f.write( format( x, float_fmt ) + " " )
        f.write( "\n" )
    
    f.close()
    return
#############################################################
def write_INRIA_mesh( xyz, faces, filename ):
    f = open( filename, "w" )
    f.write( "MeshVersionFormatted 2\n" )
    f.write( "Dimension\n3\n\n" )
    
    # vertices
    f.write( "Vertices\n" + str(len(xyz)) + "\n" )
    for v in xyz:
        for x in v:
            f.write( format( x, float_fmt ) + " " )
        f.write( "0\n" )
    f.write( "\n" )
    
    # (boundary) edges
    bed = boundary_edges( faces, len(xyz) )
    #bed = []
    f.write( "Edges\n" + str(len(bed)) + "\n" )
    for e in bed:
        f.write( str(e[0]) + " " + str(e[1]) + " 0\n" )
    f.write( "\n" )
    
    # faces
    tri = []
    qua = []
    for fa in faces:
        if len(fa) == 3:
            tri.append( fa )
        elif len(fa) == 4:
            qua.append( fa )
    
    if len(tri) > 0:
        f.write( "Triangles\n" + str(len(tri)) + "\n" )
        for fa in tri:
            for i in fa:
                f.write( str(i) + " " )
            f.write( "0\n" )
        f.write( "\n" )

    if len(qua) > 0:
        f.write( "Quadrangles\n" + str(len(qua)) + "\n" )
        for fa in qua:
            for i in fa:
                f.write( str(i) + " " )
            f.write( "0\n" )
        f.write( "\n" )

    f.write( "\nEnd" )
    f.close()
    return
#############################################################
def write_stl( xyz, faces, filename ):
    f = open( filename, "w" )
    f.write( "solid " + filename + "\n\n" )

    for fa in faces:
        if len(fa) != 3:
            print( "STL format only supports triangles" )
            sys.exit(2)
        # compute normal
        u = np.array( xyz[fa[1]-1] ) - np.array( xyz[fa[0]-1] )
        v = np.array( xyz[fa[2]-1] ) - np.array( xyz[fa[0]-1] )
        w = np.cross( u, v )
        w = w / np.sqrt( np.sum( np.power(w,2) ) )

        f.write( "facet normal " )
        for x in w:
             f.write( format( x, float_fmt ) + " " )
        f.write( "\n" )
        f.write( "    outer loop\n")
        for i in fa:
            f.write( "        vertex ")
            for x in xyz[i-1]:
                f.write( format( x, float_fmt ) + " " )
            f.write( "\n" )
        f.write( "    endloop\n")
        f.write( "endfacet\n\n")
    f.write( "\nendsolid " + filename )
    f.close()
    return
#############################################################
def write_vtk( xyz, faces, filename ):
    f = open( filename, "w" )
    f.write( "# vtk DataFile Version 2.0\n" )
    f.write( filename[0:-4] + "\n" )
    f.write( "ASCII\n" )
    f.write( "DATASET UNSTRUCTURED_GRID\n" )
    
    f.write( "POINTS " + str(len(xyz)) + " double\n" )
    for v in xyz:
        for x in v:
            f.write( format( x, float_fmt ) + " " )
        f.write( "\n" )
        
    f.write( "CELLS " + str(len(faces)) + " " + str(4*len(faces)) + "\n" )
    for t in faces:
        f.write( "3 " )
        for i in t:
            f.write( str(i-1) + " " )
        f.write( "\n" )

    f.write( "CELL_TYPES " + str(len(faces)) + "\n" )
    for t in faces:
        f.write( "5\n" )
    f.close()
    return
#############################################################
"""
def boundary_edges( faces, nv ):
    nhe = nv*( 1 + (nv-1)/2 ) - 1
    hashtab = -np.ones( (nhe,2), dtype=int )
    isb = []
    for i, f in enumerate(faces):
        m = len(f)
        isbf = []
        for j, v, in enumerate(f):
            e = [ f[j], f[(j+1)%m] ]
            hid = hash_integer_pair( e[0]-1, e[1]-1 )
            if hashtab[hid,0] < 0:
                hashtab[hid,0] = i
                hashtab[hid,1] = j
                isbf.append( True )
            else:
                isbf.append( False )
                isb[hashtab[hid,0]][hashtab[hid,1]] = False
        isb.append( isbf )
    bed = []
    for i, f in enumerate(faces):
        for j in range(len(f)):
            if isb[i][j]:
                bed.append( [ f[j], f[(j+1)%len(f)] ] )
    return bed
"""
def boundary_edges( faces, nv ):
    h = []
    isb = []
    for i, f in enumerate(faces):
        m = len(f)
        isbf = []
        for j in range(m):
            e = [f[j], f[(j+1)%m]]
            ih = hash_integer_pair(e[0]-1, e[1]-1)
            added = False
            for jh in range(len(h)):
                if h[jh][0] == ih:
                    # collision
                    isbf.append(False)
                    isb[h[jh][2]][h[jh][1]] = False
                    h.pop(jh)
                    added = True
                    break
                elif h[jh][0] < ih:
                    h.insert(jh, [ih,j,i])
                    isbf.append(True)
                    added = True
                    break
            if not added:
                h.append([ih,j,i])
                isbf.append(True)
        isb.append(isbf)
    
    bed = []
    for i, f in enumerate(faces):
        m = len(f)
        for j in range(m):
            if isb[i][j]:
                bed.append([f[j], f[(j+1)%m]])
    return bed

#############################################################
def hash_integer_pair( i, j ):
    if i == j:
        return -1
    else:
        p = min(i,j)
        q = max(i,j)
        return p + q*(q-1)/2
#############################################################



helpmsg = "-f <faces> [-u <uv-coords>] [-x <xyz-coords>] [-p <polynomial>] -o <output mesh>"


def main(argv):
    inputfaces = ""
    inputuv = ""
    inputxyz = ""
    inputpoly = ""
    outputmesh = ""
    
    faces = []
    uv = []
    xyz = []
    poly = []
    try :
        opts, args = getopt.getopt( 
            argv, 
            "hf:u:x:p:o:",
            ["help", "faces", "uv", "xyz", "poly", "output"])
    except getopt.GetoptError:
        #usage()
        print( helpmsg )
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print helpmsg
            sys.exit()
        elif opt in ("-f", "--faces"):
            inputfaces = arg
        elif opt in ("-u", "--uv"):
            inputuv = arg
        elif opt in ("-x", "--xyz"):
            inputxyz = arg
        elif opt in ("-p", "--poly"):
            inputpoly = arg
        elif opt in ("-o", "--output"):
            outputmesh = arg
    
    """
    print "input faces =", inputfaces
    print "input uv    =", inputuv
    print "input xyz   =", inputxyz
    print "input poly  =", inputpoly
    print "output mesh =", outputmesh
    """
    
    if not inputfaces:
        print "faces mandatory"
        sys.exit(2)
        
    if not outputmesh:
        print "output mesh mandatory"
        sys.exit(2)
    
    # check output mesh extension
    name, extension = os.path.splitext( outputmesh )#[1]
    if extension not in (".obj", ".stl", ".mesh", ".vtk"):
        print "mesh format "+ extension + " not supported (wavefront OBJ, STL, INRIA mesh, MFEM mesh, VTK)"
        sys.exit(2)
        
    # read faces
    faces = []
    f = open( inputfaces, "r" )
    for l in f:
        faces.append( [int(a) for a in l.split() ] )
    f.close()
    
    # read uv-coords if specified
    if inputuv:
        #uv = read_points( uv, inputuv, 2, "uv coordinates" )
        uv = np.loadtxt( inputuv )
        
    # read xyz-coords if specified
    if inputxyz:
        #xyz = read_points( xyz, inputxyz, 3, "xyz coordinates" )
        xyz = np.loadtxt( inputxyz )
    
    # read polynomial
    if inputpoly:
        poly = read_polynomial( inputpoly )
        
    if inputuv:
        if not inputxyz:
            if inputpoly:
                #uv = np.asarray( uv )
                xyz = np.polynomial.chebyshev.chebval2d( 
                        uv[:,0],
                        uv[:,1], 
                        poly )
                xyz = zip(xyz[0], xyz[1], xyz[2])
            else:
                xyz = np.hstack( ( uv, np.zeros( (len(uv),1) ) ) )
    else:
        if not inputxyz:
            print "either uv or xyz coordinates must be specified"
            sys.exit(2)
    
    if extension == ".obj":        
        write_obj( xyz, uv, faces, outputmesh )
    elif extension == ".mesh":
        write_MFEM_mesh( xyz, faces, name + "_mfem" + extension )
        write_INRIA_mesh( xyz, faces, name + "_inria" + extension )
    elif extension == ".stl":
        write_stl( xyz, faces, outputmesh )
    elif extension == ".vtk":
        write_vtk( xyz, faces, outputmesh )  
    else:
        print "output mesh format not supported"
        sys.exit(2)
        
    
if __name__ == "__main__":
    main(sys.argv[1:])
