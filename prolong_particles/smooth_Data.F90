
Module smooth_Data

  use tree, ONLY: maxblocks_tr, nchild

#include "constants.h"  
#include "Flash.h" 

  logical, save :: SmoothParticles = .false.
  
  integer, save :: smooth_mype
  
  logical, save :: dont_skip_smoothed_particles = .false.
  
  integer, save :: lrefine_smooth_max
  integer, save :: lrefine_smooth_min
  real, save :: smooth_radius_factor
  
  logical, save :: smooth_initialized = .false.
  
  integer :: upper_lrefine, lower_lrefine
  
  !Data structures for restriction and prolongation

  integer, allocatable, dimension(:), save       :: send_prolong_req
  real, allocatable, dimension(:,:,:,:,:), save  :: send_prolong_data
  real, allocatable, dimension(:,:,:),save       :: recv_prolong_data

  integer, save                :: hg_restrict_n1, hg_restrict_n2, hg_restrict_n3
  integer, save                ::  nbuf_prolong, nbuf_restrict

  integer, save                :: nxl1, nxl2, nyl1, nyl2, nzl1, nzl2
  integer, save                :: nxr1, nxr2, nyr1, nyr2, nzr1, nzr2

  real, allocatable, save      :: send_restrict_data(:,:,:,:)
  real, allocatable, save      :: recv_restrict_data(:,:,:)
  integer, allocatable, save      :: send_restrict_req(:)  

  real,allocatable,dimension(:,:,:), save :: Px
  real,allocatable,dimension(:,:,:), save :: Py
  real,allocatable,dimension(:,:,:), save :: Pz

  integer, dimension(8),save   :: n1off, n2off, n3off

  integer, save :: nmax


  ! Ranges of interior indices for blocks.
  
  integer,save :: hg_ili, hg_iui, hg_jli, hg_jui, hg_kli, hg_kui
  
  ! Ranges of exterior indices for blocks.
  
  integer,save :: hg_ile, hg_iue, hg_jle, hg_jue, hg_kle, hg_kue

  integer, save :: maxMeshRefineLevel, minMeshRefineLevel

  integer, allocatable,dimension(:), save :: smooth_saveNodeType
  logical, allocatable,dimension(:), save :: smooth_saveNewChild

  ! particle array for purpose of sending particles
  ! located in blocks with lrefine>lsmooth to 
  ! coarse parent block and get acceleration
  real, save, allocatable, dimension(:,:) :: particles_smooth
  integer, save :: num_smooth
  integer, save :: max_num_smooth_part
  integer, save :: num_smooth_part_props
  
  integer, save :: blk_og
  integer, save :: blk_coarse
  integer, save :: proc_og
  integer, save :: proc_coarse
  integer, save :: posx
  integer, save :: posy
  integer, save :: posz
  integer, save :: accx
  integer, save :: accy
  integer, save :: accz
  integer, save :: mass
  integer, save :: pid
  integer, save :: lref_dest
  
  integer,save,dimension(MAXBLOCKS) :: smooth_blkList
  integer,save :: smooth_blkCount

  real, allocatable, dimension(:,:) :: destBuf, sourceBuf

  integer, dimension(MDIM) :: smooth_posAttrib
  
end Module smooth_Data

