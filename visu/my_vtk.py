float_fmt = "08.6e"

##################################################
def write_vtk_mesh(filename, xyz, faces, faceref):
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

    f.write( "CELL_DATA " + str(len(faces)) + "\n" )
    f.write( "SCALARS hyperface int\n" )
    f.write( "LOOKUP_TABLE default\n" )
    for r in faceref:
        f.write( str(r) + "\n" )

    f.close()
    return


##################################################
