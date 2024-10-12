! ############################################################
!            Satellite Data Simulation Unit v2
!                 -  def_file_list_CReSS -
!
! This suboutine defines the list of input CRM files for the
! Cloud Resolving Storm Simulator (CReSS).
! TO CALL THIS SUBROUTINE, REPLACE "def_file_list" WITH 
! "def_file_list_CReSS" IN "SDSU-v2.f90" AND IN "beamconv/BEAM_SIM.f90". 
! Note that the CReSS output nomenclature changes from a simulation 
! to another, so this subroutine only offers a template and may not be
! applicable to all CReSS outputs.
!
 Subroutine def_file_list_CReSS
!`
   Use MAINparameters
   Implicit NONE
!
! Local variables
!
   Integer :: ifile, iargc
!
! The list of files containing snapshots from CReSS-simulated 
! Baiu-front systems. The simulated data are NOT included in the
! the SDSU package.
! 
! Redefine this list when the files are replaced by user-provided geophysical
! parameters. 
! 
   nmodel = 37   ; allocate ( CRM_file_list(nmodel) )
   CRM_file_list(1)  = 'BAIUa/BAIU2009062812Z.dmp00000000.united.bin'
   CRM_file_list(2)  = 'BAIUa/BAIU2009062812Z.dmp00003600.united.bin'
   CRM_file_list(3)  = 'BAIUa/BAIU2009062812Z.dmp00007200.united.bin'
   CRM_file_list(4)  = 'BAIUa/BAIU2009062812Z.dmp00010800.united.bin'
   CRM_file_list(5)  = 'BAIUa/BAIU2009062812Z.dmp00014400.united.bin'
   CRM_file_list(6)  = 'BAIUa/BAIU2009062812Z.dmp00018000.united.bin'
   CRM_file_list(7)  = 'BAIUa/BAIU2009062812Z.dmp00021600.united.bin'
   CRM_file_list(8)  = 'BAIUa/BAIU2009062812Z.dmp00025200.united.bin'
   CRM_file_list(9)  = 'BAIUa/BAIU2009062812Z.dmp00028800.united.bin'
   CRM_file_list(10) = 'BAIUa/BAIU2009062812Z.dmp00032400.united.bin'
   CRM_file_list(11) = 'BAIUa/BAIU2009062812Z.dmp00036000.united.bin'
   CRM_file_list(12) = 'BAIUa/BAIU2009062812Z.dmp00039600.united.bin'
   CRM_file_list(13) = 'BAIUa/BAIU2009062812Z.dmp00043200.united.bin'
   CRM_file_list(14) = 'BAIUa/BAIU2009062812Z.dmp00046800.united.bin'
   CRM_file_list(15) = 'BAIUa/BAIU2009062812Z.dmp00050400.united.bin'
   CRM_file_list(16) = 'BAIUa/BAIU2009062812Z.dmp00054000.united.bin'
   CRM_file_list(17) = 'BAIUa/BAIU2009062812Z.dmp00057600.united.bin'
   CRM_file_list(18) = 'BAIUa/BAIU2009062812Z.dmp00061200.united.bin'
   CRM_file_list(19) = 'BAIUa/BAIU2009062812Z.dmp00064800.united.bin'
   CRM_file_list(20) = 'BAIUa/BAIU2009062812Z.dmp00068400.united.bin'
   CRM_file_list(21) = 'BAIUa/BAIU2009062812Z.dmp00072000.united.bin'
   CRM_file_list(22) = 'BAIUa/BAIU2009062812Z.dmp00075600.united.bin'
   CRM_file_list(23) = 'BAIUa/BAIU2009062812Z.dmp00079200.united.bin'
   CRM_file_list(24) = 'BAIUa/BAIU2009062812Z.dmp00082800.united.bin'
   CRM_file_list(25) = 'BAIUa/BAIU2009062812Z.dmp00086400.united.bin'
   CRM_file_list(26) = 'BAIUa/BAIU2009062812Z.dmp00090000.united.bin'
   CRM_file_list(27) = 'BAIUa/BAIU2009062812Z.dmp00093600.united.bin'
   CRM_file_list(28) = 'BAIUa/BAIU2009062812Z.dmp00097200.united.bin'
   CRM_file_list(29) = 'BAIUa/BAIU2009062812Z.dmp00100800.united.bin'
   CRM_file_list(30) = 'BAIUa/BAIU2009062812Z.dmp00104400.united.bin'
   CRM_file_list(31) = 'BAIUa/BAIU2009062812Z.dmp00108000.united.bin'
   CRM_file_list(32) = 'BAIUa/BAIU2009062812Z.dmp00111600.united.bin'
   CRM_file_list(33) = 'BAIUa/BAIU2009062812Z.dmp00115200.united.bin'
   CRM_file_list(34) = 'BAIUa/BAIU2009062812Z.dmp00118800.united.bin'
   CRM_file_list(35) = 'BAIUa/BAIU2009062812Z.dmp00122400.united.bin'
   CRM_file_list(36) = 'BAIUa/BAIU2009062812Z.dmp00126000.united.bin'
   CRM_file_list(37) = 'BAIUa/BAIU2009062812Z.dmp00129600.united.bin'
!
! / Added for SDSU v2.1.2 /
!
! Define CRM_file_list to be '2BREADIN' if the user wants the input file list
! to be read in from command line arguments rather than specified internally. 
!
!    nmodel = 1  ; CRM_file_list = '2BREADIN'
!
! When this option is activated, provide your file list as
! % sdsu.x FILE1 FILE2 FILE3 ...
! 'nmodel' will be overwritten by the number of the arguments provided.
!
   If ( CRM_file_list(1) == '2BREADIN' ) then
!
      nmodel = iargc()
      Do ifile = 1, nmodel
         Call getarg ( ifile, CRM_file_list(ifile) )
      End Do
!
   End If
!
   Return
 End Subroutine def_file_list_CReSS
