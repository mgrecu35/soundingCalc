! ############################################################
!            Satellite Data Simulation Unit v2
!                 -  def_file_list_CReSS -
!
! This suboutine defines the list of input CRM files for the
! Cloud Resolving Storm Simulator (CReSS). Please see the
! instructions in "rd_CRM_CReSS.f90" to properly activate this
! subroutine.
!
 Subroutine def_file_list_CReSS
!`
   Use MAINparameters
   Implicit NONE
!
! The list of files containing snapshots from CReSS-simulated 
! Baiu-front systems. The simulated data are NOT included in the
! the SDSU package.
! 
! Redefine this list when the files are replaced by user-provided geophysical
! parameters. 
! 
   nmodel = 4    ; allocate ( CRM_file_list(nmodel) )
   CRM_file_list(1) = 'BAIUa/BAIU2009062812Z.200906290000Z.bin'
   CRM_file_list(2) = 'BAIUa/BAIU2009062812Z.200906290600Z.bin'
   CRM_file_list(3) = 'BAIUa/BAIU2009062812Z.200906291200Z.bin'
   CRM_file_list(4) = 'BAIUa/BAIU2009062812Z.200906291800Z.bin'
!
   Return
 End Subroutine def_file_list_CReSS
