! This module contains the routines to partition mesh.
! This library is copied and modified from the original
! SPECFEM3D package (Komatitsch and Tromp 1999, Peter et al. 2011)
! REVISION
!   HNG, Dec 09,2011; HNG, Jul 12,2011; HNG, Apr 09,2010
module partmesh_scotch
use global
use partmesh_library
use string_library,only:str2int,str2real
implicit none

!include './scotchf.h'
include "ptscotchf.h"

! number of partitions
integer :: npart ! e.g. 4 for partitioning for 4 CPUs or 4 processes

! mesh arrays
integer,dimension(:),allocatable :: part

integer,dimension(:),allocatable :: xadj
integer,dimension(:),allocatable :: adjncy
integer,dimension(:),allocatable :: nnodeelmt
integer,dimension(:),allocatable :: nodeselmt
integer,dimension(:),allocatable :: elmt_weight

integer,dimension(:),pointer :: glob2loc_elmt
integer,dimension(:),pointer :: glob2loc_nodes_npart
integer,dimension(:),pointer :: glob2loc_nodes_parts
integer,dimension(:),pointer :: glob2loc_nodes

integer,dimension(:),pointer :: tab_size_interfaces, tab_interfaces
integer,dimension(:),allocatable :: my_interfaces
integer,dimension(:),allocatable :: my_nb_interfaces
integer ::  ninterfaces
integer :: my_ninterface

! max number of elements that contain the same node.
integer :: nsize !integer(long).
integer :: nb_edges

integer :: max_neighbour ! Real maximum number of neighbours per element
! integer(long). Majoration of the maximum number of neighbours per element
integer :: sup_neighbour

integer :: nnode_loc, nelmt_loc
integer :: num_elmnt, num_node

!! boundaries
!integer :: nelmt2D_xmin,nelmt2D_xmax,nelmt2D_ymin,nelmt2D_ymax,nelmt2D_bottom, &
!nelmt2D_top
!integer,dimension(:),allocatable :: ibelm_xmin,ibelm_xmax,ibelm_ymin,          &
!ibelm_ymax,ibelm_bottom,ibelm_top
!integer,dimension(:),allocatable :: nodes_ibelm_xmin,nodes_ibelm_xmax,         &
!nodes_ibelm_ymin
!integer,dimension(:,:),allocatable :: nodes_ibelm_ymax,nodes_ibelm_bottom,     &
!nodes_ibelm_top
!
!! moho surface (optional)
!integer :: nelmt2D_moho
!integer,dimension(:),allocatable :: ibelm_moho
!integer,dimension(:,:),allocatable :: nodes_ibelm_moho

character(len=60) :: lineword(9)
character(len=256) :: prname

!logical,dimension(:),allocatable :: mask_nodeselmt
!integer,dimension(:),allocatable :: used_nodeselmt

real(kind=kreal),dimension(SCOTCH_GRAPHDIM) :: scotchgraph !double precision
real(kind=kreal),dimension(SCOTCH_STRATDIM) :: scotchstrat !double precision

!pll
real(kind=kreal),dimension(:,:),allocatable :: mat_prop ! double precision
integer,dimension(:),allocatable :: undef_mat_domain
integer :: nmatblk_undef,imat
character (len=30),dimension(:,:),allocatable :: undef_mat_prop

! default mesh file directory
character(len=256) :: out_phead ! output path and header

contains
!-------------------------------------------------------------------------------

! reads in mesh files
subroutine read_mesh_files
implicit none

print*,'matfile: ',trim(adjustl(matfile))

! reads material definitions
!
! note: format of nummaterial_velocity_file must be
!
! #(1)material_domain_id #(2)material_id  #(3)gam  #(4)ym  #(5)nu  
! #(6)phi  #(7)anisotropy_flag
!
! where
!     material_domain_id : 1=elastic / 2=acoustic / 3=poroelastic
!     material_id      : number of material/volume
!     rho or gam       : density or unit weight
!     ym               : Young's modulus
!     nu               : Poisson's ratio
!     phi              : Angle of friction
!     coh              : Cohesion sternght
!     psi              : Angle of dilation
nmatblk_undef = 0
allocate(mat_prop(6,nmatblk))
allocate(undef_mat_prop(6,nmatblk_undef),undef_mat_domain(nmatblk))
do imat=1,nmatblk
  if(isdensity)then
    mat_prop(1,imat) = rho_blk(imat)
  else
    mat_prop(1,imat) = gam_blk(imat)
  endif
  mat_prop(2,imat) = ym_blk(imat)
  mat_prop(3,imat) = nu_blk(imat)
  mat_prop(4,imat) = phi_blk(imat)
  mat_prop(5,imat) = coh_blk(imat)
  mat_prop(6,imat) = psi_blk(imat)
enddo

end subroutine read_mesh_files
!===============================================================================

! checks valence of nodes
subroutine check_valence
implicit none
integer :: ispec,inode
logical,dimension(:),allocatable :: mask_nodeselmt
integer,dimension(:),allocatable :: used_nodeselmt

allocate(mask_nodeselmt(nnode))
allocate(used_nodeselmt(nnode))
mask_nodeselmt(:) = .false.
used_nodeselmt(:) = 0
do ispec = 1, nelmt
  do inode = 1, NGNOD
    mask_nodeselmt(g_num(inode,ispec)) = .true.
    used_nodeselmt(g_num(inode,ispec)) = used_nodeselmt(g_num(inode,ispec)) + 1
  enddo
enddo
print *, 'nodes valence: '
print *, '  min = ',minval(used_nodeselmt(:)),'max = ', maxval(used_nodeselmt(:))
do inode = 1, nnode
  if (.not. mask_nodeselmt(inode)) then
    stop 'ERROR : nodes not used.'
  endif
enddo
nsize = maxval(used_nodeselmt(:)) ! max number of element sharing a node
sup_neighbour = NGNOD * nsize - (NGNOD + (NGNOD/2 - 1)*nfaces)
print*, '  nsize = ',nsize, 'sup_neighbour = ', sup_neighbour

deallocate(mask_nodeselmt,used_nodeselmt)
end subroutine check_valence
!===============================================================================

! divides model into partitions using scotch library functions
subroutine scotch_partitioning
implicit none
integer :: istat

g_num = g_num - 1

! determines maximum neighbors based on 1 common node
allocate(xadj(1:nelmt+1))
allocate(adjncy(1:sup_neighbour*nelmt))
allocate(nnodeelmt(1:nnode))
allocate(nodeselmt(1:nsize*nnode))
call mesh2dual_ncommonnode(nelmt,nnode,nsize,sup_neighbour,g_num,xadj,      &
adjncy,nnodeelmt,nodeselmt,max_neighbour,1)
print*, 'mesh2dual: '
print*, '  max_neighbour = ',max_neighbour
print*, '  sup_neighbour = ', sup_neighbour
!print*,xadj
!print*,adjncy
nb_edges = xadj(nelmt+1)

! allocates & initializes partioning of elements
allocate(part(1:nelmt))
part(:) = -1

! initializes
! elements load array
allocate(elmt_weight(1:nelmt))

! uniform load
elmt_weight(:) = 1

print*,'mat id:',minval(mat_id),maxval(mat_id)
! in case of acoustic/elastic simulation, weights elements accordingly
call acoustic_elastic_load(elmt_weight,nelmt,nmatblk,mat_id,mat_domain)

! SCOTCH partitioning
call scotchfstratinit (scotchstrat(1), istat)
  if (istat /= 0) then
    stop 'ERROR : MAIN : Cannot initialize strat'
endif

! no need to use this for default strategy
!call scotchfstratgraphmap (scotchstrat(1), trim(scotch_strategy), istat)
!  if (istat /= 0) then
!    stop 'ERROR : MAIN : Cannot build strat'
!endif

call scotchfgraphinit (scotchgraph(1), istat)
if (istat /= 0) then
    stop 'ERROR : MAIN : Cannot initialize graph'
endif

! fills graph structure : see user manual (scotch_user5.1.pdf, page 72/73)
! arguments: #(1) graph_structure     #(2) baseval(either 0/1) 
!            #(3) number_of_vertices  #(4) adjacency_index_array 
!            #(5) adjacency_end_index_array (optional)
!            #(6) vertex_load_array (optional) #(7) vertex_label_array
!            #(7) number_of_arcs                    #(8) adjacency_array
!            #(9) arc_load_array (optional)      #(10) istator
call scotchfgraphbuild (scotchgraph(1), 0, nelmt, &
                      xadj(1), xadj(1), &
                      elmt_weight (1), xadj (1), &
                      nb_edges, adjncy(1), &
                      adjncy(1), istat)

! w/out element load, but adjacency array
!call scotchfgraphbuild (scotchgraph (1), 0, nelmt, &
!                      xadj (1), xadj (1), &
!                      xadj (1), xadj (1), &
!                      nb_edges, adjncy (1), &
!                      adjncy (1), istat)


if (istat /= 0) then
    stop 'ERROR : MAIN : Cannot build graph'
endif

call scotchfgraphcheck (scotchgraph (1), istat)
if (istat /= 0) then
    stop 'ERROR : MAIN : Invalid check'
endif

call scotchfgraphpart (scotchgraph (1), npart, scotchstrat(1),part(1),istat)
if (istat /= 0) then
    stop 'ERROR : MAIN : Cannot part graph'
endif

call scotchfgraphexit (scotchgraph (1), istat)
if (istat /= 0) then
    stop 'ERROR : MAIN : Cannot destroy graph'
endif

call scotchfstratexit (scotchstrat(1), istat)
if (istat /= 0) then
    stop 'ERROR : MAIN : Cannot destroy strat'
endif


! re-partitioning puts acoustic-elastic coupled elements into same partition
!  integer  :: nfaces_coupled
!  integer, dimension(:,:), allocatable  :: faces_coupled
!    call acoustic_elastic_repartitioning (nelmt, nnode, g_num, &
!                   nmatblk, mat_id , mat_prop, &
!                   sup_neighbour, nsize, &
!                   npart, part, nfaces_coupled, faces_coupled)


! local number of each element for each partition
call construct_glob2loc_elmt(nelmt, part, glob2loc_elmt,npart)

! local number of each node for each partition
call construct_glob2loc_nodes(nelmt,nnode,nsize,nnodeelmt,nodeselmt,   &
part,glob2loc_nodes_npart,glob2loc_nodes_parts,glob2loc_nodes,npart)

! mpi interfaces
! MPI interfaces
! acoustic/elastic/poroelastic boundaries will be split into different 
! MPI partitions
call build_interfaces(nelmt,sup_neighbour,part,g_num, &
                         xadj,adjncy,tab_interfaces, &
                         tab_size_interfaces,ninterfaces, &
                         npart)

!! call detect_ghost(out_phead,nelmt,nnode,part+1,g_num+1,npart,              &
!! glob2loc_elmt(0:nelmt-1)+1) ! I need all indices starting from 1 not 0
!! acoustic/elastic boundaries WILL BE SEPARATED into different MPI partitions
!! TODO:WARNING:call below may not be necessary
!call construct_interfaces(nelmt,sup_neighbour,part,g_num,xadj,adjncy,         &
!tab_interfaces,tab_size_interfaces,ninterfaces,npart)
!!or: uncomment if you want acoustic/elastic boundaries NOT to be separated into
!!different MPI partitions
!!call Construct_interfaces_no_ac_el_sep(nelmt,sup_neighbour,part,g_num,xadj,  &
!!adjncy,tab_interfaces,tab_size_interfaces,ninterfaces,nmatblk,          &
!!mat_prop(3,:),mat_id,npart)
end subroutine scotch_partitioning
!===============================================================================

! writes out new Databases files for each partition
subroutine write_mesh_databases
implicit none
integer :: ipart
!character(len=20) :: format_str
!character(len=80) :: out_fname

allocate(my_interfaces(0:ninterfaces-1))
allocate(my_nb_interfaces(0:ninterfaces-1))

! writes out Database file for each partition
do ipart = 0, npart-1

  ! opens output file
  !write(out_fname, fmt=format_str)trim(out_phead)//'proc',ipart,'_database'
  !print*,trim(out_fname); stop
  !open(unit=15,file=trim(out_fname),&
  !     status='unknown', action='write', form='formatted', iostat = istat)
  !if( istat /= 0 ) then
  ! print*,'error file open:',trim(out_fname)
  ! stop
  !endif

  ! gets number of nodes
  call write_glob2loc_nodes_database(out_path,coordfile,ipart,         &
  nnode_loc,g_coord,glob2loc_nodes_npart,glob2loc_nodes_parts,           &
  glob2loc_nodes,nnode,npart)

  ! gets number of spectral elements
  call write_partition_database(out_path,confile,idfile,ipart,nelmt_loc,nelmt, &
  g_num,glob2loc_elmt,glob2loc_nodes_npart,glob2loc_nodes_parts,            &
  glob2loc_nodes,part,mat_id,npart)

  ! Below portion is needed only if we want to have indvidual material list file
  ! for each partition.
  !call write_material_properties_database(out_path,matfile,nmatblk,      &
  !mat_domain,type_blk,mfile_blk,mat_prop,nmatblk_viscoelas,  &
  !imat_visco,visco_rheology,nmaxwell,muratio_blk,viscosity_blk, &
  !nwmat,waterid,ipart,npart)

  ! gets number of MPI interfaces
  call write_interfaces_database(out_phead,tab_interfaces,tab_size_interfaces, &
                                 ipart,ninterfaces,my_ninterface,my_interfaces,&
                                 my_nb_interfaces,glob2loc_elmt,               &
                                 glob2loc_nodes_npart,glob2loc_nodes_parts,    &
                                 glob2loc_nodes, 1, npart)

  ! writes out MPI interfaces elements
  !if (my_ninterface == 0) then
  ! avoids problem with maxval for empty array my_nb_interfaces
  ! write(IIN_database) my_ninterface, 0       
  !else
  ! write(IIN_database) my_ninterface, maxval(my_nb_interfaces)
  !endif

  call write_interfaces_database(out_phead,tab_interfaces,tab_size_interfaces, &
                                 ipart,ninterfaces,my_ninterface,my_interfaces,&
                                 my_nb_interfaces,glob2loc_elmt,               &
                                 glob2loc_nodes_npart,glob2loc_nodes_parts,    &
                                 glob2loc_nodes,2,npart)
enddo

deallocate(mat_id,g_coord,xadj,adjncy,nnodeelmt,nodeselmt,elmt_weight)
deallocate(tab_size_interfaces,tab_interfaces,my_interfaces,my_nb_interfaces)

deallocate(mat_prop)
deallocate(undef_mat_prop,undef_mat_domain)
if(nwmat>0)deallocate(waterid)

! write BC files
write(*,'(a)',advance='no')'writing bc files...'

if(ISDISP_DOF)then
  call write_ssbc(out_path,inp_path,uxfile,part+1,npart,glob2loc_elmt(0:nelmt-1)+1)
  call write_ssbc(out_path,inp_path,uyfile,part+1,npart,glob2loc_elmt(0:nelmt-1)+1)
  call write_ssbc(out_path,inp_path,uzfile,part+1,npart,glob2loc_elmt(0:nelmt-1)+1)
endif
! Free surface file
call write_ssbcFS(out_path,inp_path,fsfile,part+1,npart,glob2loc_elmt(0:nelmt-1)+1)
if(infbc)then
  call write_ssbcINF(out_path,inp_path,inffile,part+1,npart,glob2loc_elmt(0:nelmt-1)+1)
endif
write(*,'(a)')'complete!'

! write traction files
if(istraction)then
  write(*,'(a)',advance='no')'writing traction files...'
  call write_traction(out_path,inp_path,trfile,part+1,npart,                   &
  glob2loc_elmt(0:nelmt-1)+1)
  write(*,'(a)')'complete!'
endif

! write mtraction files
if(ismtraction)then
  write(*,'(a)',advance='no')'writing mtraction files...'
  call write_traction(out_path,inp_path,mtrfile,part+1,npart,                   &
  glob2loc_elmt(0:nelmt-1)+1,ismagnetpar=.true.)
  write(*,'(a)')'complete!'
endif

! write fault slip files
if(iseqsource.and.eqsource_type==3)then
  write(*,'(a)',advance='no')'writing fault slip files...'
  call write_faultslip(out_path,inp_path,faultslipfile_plus,part+1,npart,      &
  glob2loc_elmt(0:nelmt-1)+1)
  call write_faultslip(out_path,inp_path,faultslipfile_minus,part+1,npart,     &
  glob2loc_elmt(0:nelmt-1)+1)
  write(*,'(a)')'complete!'
endif
deallocate(glob2loc_nodes_npart,glob2loc_nodes_parts,glob2loc_nodes)

!!deallocate(ibelm_xmin,nodes_ibelm_xmin)
!!deallocate(ibelm_xmax,nodes_ibelm_xmax)
!!deallocate(ibelm_ymin,nodes_ibelm_ymin)
!
!write(*,'(a)',advance='yes')'finding interfaces...'
!! I need all indices starting from 1 not 0
!call find_interface(nelmt,nnode,part+1,g_num+1,npart)
!write(*,'(a)')'complete!'
!deallocate(g_num,part)
!write(*,'(a)',advance='yes')'writing interfaces...'
!! I need all indices starting from 1 not 0
!call detect_ghost(out_phead,nnode,npart,max_neighbour,                   &
!glob2loc_elmt(0:nelmt-1)+1)
!write(*,'(a)')'complete!'
!deallocate(glob2loc_elmt)

write(*,*)'number of partitions: ',npart
!print*, 'finished successfully'
!write(*,*)'-----------------------------------'
end subroutine write_mesh_databases
!===============================================================================

end module partmesh_scotch
!===============================================================================
