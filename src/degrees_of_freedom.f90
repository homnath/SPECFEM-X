! DESCRIPTION
!  This module conains the routines to set and analyze the degrees of freedoms.
! DEVELOPER
!  Hom Nath Gharti, Princeton University
! REVISION
!  HNG, Feb 19, 2016; HNG, Jul 12,2011; HNG, Apr 09,2010; HNG, Dec 08,2010
! TODO
!  - implement power law
module dof
implicit none
character(len=250),private :: myfname=' => degrees_of_freedom.f90'
character(len=500),private :: errsrc

contains
!-------------------------------------------------------------------------------

! This subroutine sets number of degrees of freedom and IDs of nodal dof  
subroutine initialize_dof()
use global,only:edofu,edofphi,idofu,nenode,idofphi,nndofu,nndofphi,nndof,      &
                nedofu,nedofphi,nedof,ISDISP_DOF,ISPOT_DOF
implicit none
integer :: i_dof,idof
! total number of degrees of freedom per node
nndof=0
nedofu=0
nedofphi=0
nedof=0

! dof IDs
idof=0
idofu=0
idofphi=0

! displacement
if(ISDISP_DOF)then
  nndof=nndof+nndofu
  nedofu=NNDOFU*nenode
  nedof=nedof+nedofu
  do i_dof=1,nndofu
    idof=idof+1
    idofu(i_dof)=idof
  enddo
  allocate(edofu(nedofu))
endif
! gravity
idof=idofu(nndofu)
if(ISPOT_DOF)then
  nndof=nndof+nndofphi
  nedofphi=NNDOFPHI*nenode
  nedof=nedof+nedofphi
  do i_dof=1,nndofphi
    idof=idof+1
    idofphi(i_dof)=idof
  enddo
  allocate(edofphi(nedofphi))
endif
end subroutine initialize_dof
!===============================================================================

! This subroutine sets IDs for the elemental degrees of freedom for
! u and \phi which may be used to map the elemental matrices
subroutine set_element_dof()
use global,only:ISDISP_DOF,ISPOT_DOF,nedofu,nedofphi,ngll,nndofu, &
                edofu,edofphi
implicit none
integer :: i,iu(NNDOFU),iphi,j,nu
integer :: iu0,iphi0

! order ux,uy,uz,\phi
edofu=-9999
edofphi=-9999

iu0=0
iphi0=0

iphi=0
nu=0
iu=0
do i=1,NGLL
  if(ISDISP_DOF)then
    iu(1)=iu0+1
    nu=nu+1
    edofu(nu)=iu(1)
    do j=2,NNDOFU
      nu=nu+1
      iu(j)=iu(j-1)+1
      edofu(nu)=iu(j)
    enddo
    iu0=iu(NNDOFU) ! this will be overwritten if POT_DOF is present
    iphi0=iu(NNDOFU)
  endif
  if(ISPOT_DOF)then
    iphi=iphi0+1
    edofphi(i)=iphi

    iu0=iphi ! this will be overwritten if DISP_DOF is present
    iphi0=iphi
  endif
enddo
!print*,edofu
!print*,edofphi; stop
return
end subroutine set_element_dof
!===============================================================================

! compute mapping of u  to face elemental matrices
subroutine set_face_vecdof(nfgll,ncomp,nxtra,imapuf)                           
implicit none                                                                    
integer,intent(in) :: nfgll,ncomp,nxtra                                          
integer,intent(out) :: imapuf(ncomp*nfgll)                                       
integer :: i,j,iu(ncomp),nu                                                       
                                                                                 
! u must be in the first order, i.e., order ux,uy,uz,phi                                 
nu=0; iu(1)=1                                                                    
do i=1,nfgll                                                                     
  nu=nu+1                                                                        
  imapuf(nu)=iu(1)                                                               
  do j=2,ncomp                                                                   
    nu=nu+1                                                                      
    iu(j)=iu(j-1)+1                                                              
    imapuf(nu)=iu(j)                                                             
  enddo                                                                          
  iu(1)=iu(ncomp)+nxtra+1
enddo
return
end subroutine set_face_vecdof
!===============================================================================                                     
                                                                                 
! compute mapping of scalar, i.e., phi to face elemental matrices
subroutine set_face_scaldof(nfgll,npos,nxtra,imapsf)                           
implicit none                                                                    
integer,intent(in) :: nfgll,npos,nxtra                                           
integer,intent(out) :: imapsf(nfgll)                                             
integer :: i,is,nu                                                               
                                                                                 
! u must be in first order, i.e., order ux,uy,uz,phi                                 
nu=0; is=npos                                                                    
do i=1,nfgll                                                                     
  nu=nu+1                                                                        
  imapsf(nu)=is                                                                  
  is=is+npos+nxtra                                                               
enddo                                                                            
return                                                                           
end subroutine set_face_scaldof                                                
!===============================================================================

! This subroutine activates degrees of freedoms.
subroutine activate_dof(errcode,errtag)
use global
use math_constants,only:ZERO
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: ios
integer :: i_elmt,imat,mdomain
integer :: inodes(ngll)

errtag="ERROR: unknown!"
errcode=-1

! Initialize all DOFs to OFF
gdof=0

! Prescribed DOFs will be made OFF in apply_bc routine.
! Displacement
if(ISDISP_DOF)then
  do i_elmt=1,nelmt
    inodes=g_num(:,i_elmt)
    imat=mat_id(i_elmt)
    mdomain=mat_domain(imat)
    ! elastic domain
    if(mdomain==ELASTIC_DOMAIN)then
      gdof(idofu,inodes)=1
    ! viselastic domain
    elseif(mdomain==VISCOELASTIC_DOMAIN)then
      gdof(idofu,inodes)=1
    ! acoustic domain
    elseif(mdomain==ACOUSTIC_DOMAIN)then
      write(errtag,'(a)')'ERROR: acoustic domain not supported!'
      return
    ! trasnition infinite/infinite domain: only gravity 
    elseif(mdomain==ELASTIC_TRINFDOMAIN .or. &
           mdomain==VISCOELASTIC_TRINFDOMAIN .or. &
           mdomain==ELASTIC_INFDOMAIN .or. &
           mdomain==VISCOELASTIC_INFDOMAIN)then
      if(infbc)then
        if(.not.isempty_blk(imat))then
        ! This will also have displacement DOF
          gdof(idofu,inodes)=1
        endif
      endif
    else
      write(errtag,'(a)')'ERROR: unsupported material domain!'
      return
    endif
  enddo
  !gdof(idofu,:)=1
endif

! Gravity
if(ISPOT_DOF)then
  ! gravity exists everywhere
  gdof(idofphi,:)=1
endif

errcode=0

end subroutine activate_dof
!===============================================================================

! This subroutine finalizes the global degrees of freedom IDs.
subroutine finalize_gdof(errcode,errtag)
use global,only:gdof,neq,nndof,nnode,g_num,part_path,proc_str,file_head
use global,only:myrank,nedof
implicit none
integer,intent(out) :: errcode
character(len=250) :: ofname
character(len=250),intent(out) :: errtag
integer :: i,istat,j

errtag="ERROR: unknown!"
errcode=-1
! Compute modified gdof
neq=0
do j=1,ubound(gdof,2)
  do i=1,ubound(gdof,1)
    if(gdof(i,j)/=0)then
      neq=neq+1
      gdof(i,j)=neq
    endif
  enddo
enddo

ofname='tmp/'//trim(file_head)//'_gdof'//trim(adjustl(proc_str))
open(unit=22,file=trim(ofname),access='stream',form='unformatted', &
status='replace',action='write',iostat=istat)
if (istat /= 0)then
  write(errtag,'(a)')'ERROR: output file "'//trim(ofname)//'" cannot be opened!'
  stop
endif
write(22)nnode
write(22)neq
write(22)gdof
close(22)
ofname='tmp/'//trim(file_head)//'_gnum'//trim(adjustl(proc_str))
open(unit=22,file=trim(ofname),access='stream',form='unformatted', &
status='replace',action='write',iostat=istat)
if (istat /= 0)then
  write(errtag,'(a)')'ERROR: output file "'//trim(ofname)//'" cannot be opened!'
  stop
endif
write(22)nnode
write(22)g_num
close(22)

! Compute nodal to global
errcode=0
return
end subroutine finalize_gdof
!===============================================================================

end module dof
!===============================================================================
