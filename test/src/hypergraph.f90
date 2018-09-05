program hypergraph

  use mod_util
  use mod_brep
  use mod_hypergraph

  implicit none

  integer                              :: fid
  type(type_intersection_data), target :: interdata
  type(type_BREP), target              :: brep
  logical, allocatable                 :: feat_edg(:), feat_ver(:)
  integer, allocatable                 :: v2x(:), valence(:)
  type(type_hyperface), allocatable    :: hyperfaces(:)
  integer                              :: nhf, nv2x
  integer                              :: i, j, he(6)
  
  call get_free_unit(fid)
  open(unit=fid, file='hypergraph/data.dat', action='read')
  
  read (fid,*) ! nv
  read (fid,*) brep%dcel%nv
  read (fid,*) ! ne
  read (fid,*) brep%ne
  read (fid,*) ! nhe
  read (fid,*) brep%dcel%ne
  read (fid,*) ! n
  read (fid,*) brep%dcel%nf
  
  allocate(brep%dcel%vertices(brep%dcel%nv))
  allocate(brep%edges(brep%ne))
  interdata%nc = brep%ne
  allocate(interdata%curves(brep%ne))
  allocate(brep%dcel%halfedges(brep%dcel%ne))
  allocate(brep%dcel%faces(brep%dcel%nf))

  read (fid,*) ! xy
  do i = 1,brep%dcel%nv
     read (fid,*) 
  end do

  read (fid,*) ! e2h
  do i = 1,brep%ne
     read (fid,*) j
     brep%edges(i)%halfedge => brep%dcel%halfedges(j)
  end do

  read (fid,*) ! esmooth
  do i = 1,brep%ne
     brep%edges(i)%curve => interdata%curves(i)
     read (fid,*) j
     if ( j == 0 ) then
        brep%edges(i)%curve%smooth = .false.
     else
        brep%edges(i)%curve%smooth = .true.
     end if
  end do

  read (fid,*) ! he
  do i = 1,brep%dcel%ne
     read (fid,*) he
     brep%dcel%halfedges(i)%face = he(1)
     if ( he(2) > 0 ) brep%dcel%halfedges(i)%twin => brep%dcel%halfedges(he(2))
     brep%dcel%halfedges(i)%prev => brep%dcel%halfedges(he(3))
     brep%dcel%halfedges(i)%next => brep%dcel%halfedges(he(4))
     brep%dcel%halfedges(i)%orig = he(5)
     brep%dcel%halfedges(i)%edge = he(6)
  end do

  read (fid,*) ! face_outer
  do i = 1,brep%dcel%nf
     read (fid,*) brep%dcel%faces(i)%outer
  end do

  read (fid,*) ! face_inner
  do i = 1,brep%dcel%nf
     read (fid,*) brep%dcel%faces(i)%ninner
     if ( brep%dcel%faces(i)%ninner > 0 ) then
        allocate(brep%dcel%faces(i)%inner(brep%dcel%faces(i)%ninner))
        read (fid,*) brep%dcel%faces(i)%inner
     end if
  end do

  read (fid,*) ! v2h
  do i = 1,brep%dcel%nv
     read (fid,*) brep%dcel%vertices(i)%halfedge
  end do

  read (fid,*) ! edge tangents
  do i = 1,brep%ne
     read (fid,*) brep%edges(i)%tangents
  end do
  
  close(fid)

  print *,'EDGE VERTICES:'
  do i = 1,brep%ne
     print *, i, ' :', get_BREPedge_vertices(brep%edges(i))
  end do

  print *,'VERTEX FACES :'
  do i = 1,brep%dcel%nv
     call get_v2f(brep%dcel, i, v2x, nv2x)
     print *, i, ' :', v2x(1:nv2x)
  end do

  print *,'VERTEX EDGES :'
  do i = 1,brep%dcel%nv
     call get_v2e(brep%dcel, i, v2x, nv2x)
     print *, i, ' :', v2x(1:nv2x)
  end do
  
  allocate(feat_edg(brep%ne))
  call get_feature_edges(brep, feat_edg)
  print *,'FEATURE EDGES :'
  do i = 1,brep%ne
     print *, i, ' :', feat_edg(i)
  end do

  call get_hyperfaces(brep, feat_edg, hyperfaces, nhf)

  print *,'HYPERFACES :'
  do i = 1,nhf
     print *,hyperfaces(i)%faces(1:hyperfaces(i)%nf)
  end do

  allocate(feat_ver(brep%dcel%nv), valence(brep%dcel%nv))
  call get_feature_vertices(brep, feat_edg, feat_ver, valence)
  print *,'FEATURE VERTICES :'
  do i = 1,brep%dcel%nv
     print *, i, feat_ver(i), valence(i)
  end do

  
end program hypergraph
