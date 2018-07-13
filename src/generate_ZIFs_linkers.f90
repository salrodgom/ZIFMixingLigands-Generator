module datatype
 use, intrinsic :: iso_fortran_env
 implicit none
 integer, parameter :: r32  = REAL32
 integer, parameter :: r64  = REAL64
 integer, parameter :: r128 = REAL128
 integer, parameter :: i32  = INT32
 integer, parameter :: i64  = INT64
end module datatype
!
module mod_random
! module for pseudo random numbers
 use datatype
 implicit none
 private
 public init_random_seed, randint, r4_uniform
 contains
 subroutine init_random_seed( )
  implicit none
  integer(i32), allocatable :: seed(:)
  integer(i32)              :: i, n, un, istat, dt(8), pid
  integer(i64)              :: t
  call random_seed(size = n)
  allocate(seed(n))
  write(6,'(a)')"Random seed:"
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
   read(un) seed
   close(un)
   write(6,'(a)')"OS provides a random number generator"
  else
   ! Fallback to XOR:ing the current time and pid. The PID is
   ! useful in case one launches multiple instances of the same
   ! program in parallel.
   call system_clock(t)
   if (t == 0) then
      call date_and_time(values=dt)
      t = (dt(1) - 1970_i64) * 365_i64 * 24 * 60 * 60 * 1000 &
           + dt(2) * 31_i64 * 24 * 60 * 60 * 1000 &
           + dt(3) * 24_i64 * 60 * 60 * 1000 &
           + dt(5) * 60 * 60 * 1000 &
           + dt(6) * 60 * 1000 + dt(7) * 1000 &
           + dt(8)
   end if
   pid = getpid()
   t = ieor(t, int(pid, kind(t)))
   do i = 1, n
      seed(i) = lcg(t)
   end do
   write(6,'(a)')"Fallback to the current time and pid."
  end if
  call random_seed(put=seed)
  write(6,*) seed
 contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729_i64
    else
       s = mod(s, 4294967296_i64)
    end if
    s = mod(s * 279470273_i64, 4294967291_i64)
    lcg = int(mod(s, int(huge(0), i64)), kind(0))
  end function lcg
 end subroutine init_random_seed
 integer(i64) function randint(i,j)
  integer(i64),intent(in)    :: i,j
  real(i64)                  :: r
  call random_number(r)
  randint = i + floor((j+1_i64-i)*r)
 end function randint
! 
 real(r64) function r4_uniform(x,y)
  implicit none
  real(r64),intent(in)       :: x,y
  real(r64)                  :: r
  call random_number(r)
  r4_uniform=(r*(y-x))+x
  return
 end function r4_uniform
end module mod_random

program Mixing_MOF_generator
 use datatype
 use mod_random
 implicit none
 integer(i32)                                  :: i_32
 integer(i64)                                  :: i,j,k,l,h,z
 integer(i64)                                  :: m,ierr,nn,iiii,n
 integer(i64)                                  :: linker_type_max = 10_i64
 integer(i64)                                  :: n_exchanged = 0_i64,  n_exchanged_try = 0_i64
 integer(i64)                                  :: n_rotated = 0_i64,    n_rotated_try   = 0_i64
 integer(i64)                                  :: n_expanded = 0_i64,   n_expanded_try  = 0_i64
 integer(i64)                                  :: n_accepted = 0_i64
 integer(i64)                                  :: ii,jj,kk
 integer(i64)                                  :: num_args
 integer(i64)                                  :: n_atoms = 0_i64,n_nodes=0_i64,n_linkers=0_i64,n_metals=0_i64
 integer(i64),parameter                        :: max_number_tries = 1000_i64
 integer(i64)                                  :: n_files = 1_i64
 integer(i64)                                  :: mc_steps, mc_max_steps=0_i64
 integer(i64)                                  :: linker_type_number = 1_i64
 integer(i64),dimension(:),  allocatable       :: genome
 integer(i64),dimension(:),  allocatable       :: histogram_molar_fraction
 integer(i64),dimension(:,:),allocatable       :: ensemble
!
 real(r64)                                     :: ExchangeProbability     = 1.0_r64
 real(r64)                                     :: VolumeChangeProbability = 0.0_r64
 real(r64)                                     :: RotationProbability     = 1.0_r64
 real(r64)                                     :: rrr,ppp,qqq,xyz
 real(r64)                                     :: atom(3),ouratom(3)
 real(r64),   parameter                        :: r_min_criteria_connectivity=0.56_r64
 real(r64),   parameter                        :: r_min_criteria_overlap = 1.14_r64
 real(r64)                                     :: cell_0(1:6) = 0.0_r64, rv(3,3),vr(3,3)
 real(r64)                                     :: cell_1(1:6) = 0.0_r64, cell_2(1:6) = 0.0_r64
 real(r64),   allocatable                      :: eee(:)
 real(r64),   dimension(:),allocatable         :: linker_type_molar_fraction
 real(r64),   allocatable                      :: backup_positions_linker(:,:)
 real(r64)                                     :: molar_fraction
!
 logical                                       :: flag = .true., wrong
 logical,     dimension(:,:), allocatable      :: a
!
 character(len=2)                              :: main_metal_element = "Zn"
 character(len=3)                              :: topology = "Xxx", nodes_code = "Xxx"
 character(len=20)                             :: spam
 character(len=100)                            :: FolderDataBase = " "
 character(len=100)                            :: CIFFilename=" "
 character(len=100)                            :: filename=" "
 character(len=80)                             :: string_stop_head= "_atom_site_occupancy"
 character(len=100)                            :: line, string
 character(len=100), dimension(:), allocatable :: args
 character(len=3),   dimension(:), allocatable :: linker_type
 character(len=3)                              :: code
 character(len=1)                              :: ss = "-"
!
 type                                          :: cluster
  character(len=3)                            :: code
  integer(i64)                                :: n_components
  integer(i64)                                :: id
  real(r64)                                   :: com(1:3)
  real(r64)                                   :: component_xcrystal(1:3,1:500)
  real(r64)                                   :: component_dcrystal(1:3,1:500)
  real(r64)                                   :: component_size(500)
  character(len=2)                            :: component_label(500)
  character(len=2)                            :: element(500)
  real(r64)                                   :: comfortably
  logical                                     :: virtual
  real(r64)                                   :: population
 end type                                      
 type(cluster),allocatable                     :: nodes(:), linkers(:)
! Generate Seed
 call init_random_seed( )
! Renormalise probabilities
 ppp = ExchangeProbability + RotationProbability + VolumeChangeProbability
 ExchangeProbability     = ExchangeProbability / ppp
 VolumeChangeProbability = VolumeChangeProbability / ppp + ExchangeProbability
 RotationProbability     = RotationProbability / ppp     + ExchangeProbability + VolumeChangeProbability
! Read inputs from command argument line
 num_args = command_argument_count()
 allocate(args(num_args))
 do i_32 = 1, num_args
  call get_command_argument(i_32,args(i_32))
 end do
 do i=1,num_args
  select case(args(i))
   case ('-h','--help')
    call print_help()
    stop
   case ('-l','--linker')
    read(args(i+1),*) linker_type_number 
    allocate(linker_type(1:linker_type_number))
    allocate(linker_type_molar_fraction(1:linker_type_number))
    k=0
    do j = 1, linker_type_number
     k=k+2
     read(args(i+k),'(a)') linker_type(j)  
     read(args(i+k+1),*) linker_type_molar_fraction(j) 
    end do
   case ('-t','--topol')
    read(args(i+1),'(3a)') topology
  end select
 end do
 if(allocated(linker_type).eqv..false.)then
  linker_type_number=1
  allocate(linker_type(1:linker_type_number))
  allocate(linker_type_molar_fraction(1:linker_type_number))
  linker_type(1)='imi'
  linker_type_molar_fraction(1)=1.0_r64
 end if
 rrr=sum(linker_type_molar_fraction)
 linker_type_molar_fraction=linker_type_molar_fraction/rrr
 write(6,'(80a)')('=',j=1,80)
 write(6,*)( linker_type(j),j=1,linker_type_number)
 write(6,*)( linker_type_molar_fraction(j)/rrr ,j=1,linker_type_number)
 write(6,'(80a)')('=',j=1,80)
 allocate(histogram_molar_fraction(linker_type_number))
 call system("if [ -f tmp  ] ; then rm tmp  ; touch tmp  ; fi")
 call system("if [ -f list ] ; then rm list ; touch list ; fi")
 do j=1,linker_type_number
  i=len(trim(linker_type(j)))
  k=LEN(TRIM(topology))
  string="ls tbp_database/???_"//topology(1:k)//"_cif_gin_all/???_???_"//linker_type(j)(1:i)//"_"//topology(1:k)//"_*.cif > tmp"
  call system(string)
  write(line,*) linker_type_molar_fraction(j)
  k=len(trim(line))
  string="awk -v type="      //linker_type(j)(1:i)//&
   " -v r="      //line(4:k)//" '{print $1,type,r}' tmp >> list"
  call system(string)
 end do
 open(111,file="list",iostat=ierr)
 n_files=0
 do
  read(111,'(a)',iostat=ierr) line 
  if(ierr/=0) exit
  read(line(1 :12),'(a)') FolderDataBase
  read(line( 1:54),'(a)') CIFFilename
  read(line(55:),*)       code,molar_fraction
  FolderDataBase=adjustl(trim( FolderDataBase ))
  CIFFilename=adjustl(trim(CIFFilename))
  write(6,'(a,1x,a3,1x,f14.7)')trim(CIFFilename),trim(code),molar_fraction
  n_files=n_files+1
  if(n_files==1)then
  !Investigate number of Metals -> number of ligands -> atoms/ligands
   open(121,file=CIFFilename,iostat=ierr,status='old')
   if(ierr/=0)stop 'File does not found'
   do 
    read(121,'(a)',iostat=ierr) line
    if(ierr/=0)exit
    if(line(1:2)==main_metal_element.or.line(3:4)==main_metal_element)then
     n_metals=n_metals+1
    end if
   end do
   close(121)
   call check_topology_composition(topology,n_metals,n_nodes,n_linkers,nodes_code)
  end if
 end do
 rewind(111)
 allocate(nodes(n_nodes))
 allocate(linkers(n_files*n_linkers))
 allocate(ensemble(n_linkers,n_files))
 allocate( a(n_nodes, n_linkers) )
 write(6,'(80a)')('=',j=1,80)
 write(6,*)topology
 write(6,*)n_metals, 'metals/uc'
 write(6,*)n_nodes, 'nodes/uc'
 write(6,*)n_linkers,'linkers/uc'
 write(6,*)n_files, 'files detected'
 write(6,*)n_files*n_linkers,' total linkers spected'
 h=0
 do i=1,linker_type_number
  h=h+int(linker_type_molar_fraction(i)*n_linkers)
 end do
 write(6,*)(int(linker_type_molar_fraction(i)*n_linkers),i=1,linker_type_number),n_linkers-h
 write(6,*)"=================#"
 z=0 ! count total number of nodes
 h=0 ! count total number of linkers
 write(6,'(80a)')('=',j=1,80)
 do i=1,n_files
  read(111,'(a)',iostat=ierr) line
  read(line(1 :12),'(a)') FolderDataBase
  read(line( 1:54),'(a)') CIFFilename
  read(line(55:),*)       code,molar_fraction
  write(6,'(a1,1x,a)')'#',trim(CIFfilename)
  open(100,file=trim(CIFfilename),status='old',iostat=ierr)
  if(ierr/=0)stop 'File does not found'
  read_cif: do
   read(100,'(a)',iostat=ierr) line
   if(ierr/=0) exit read_cif
   if(line(1:14)=="_cell_length_a")then
    read(line,*)spam,cell_0(1)
    cycle read_cif
   end if
   if(line(1:14)=="_cell_length_b")then
    read(line,*)spam,cell_0(2)
    cycle read_cif
   end if
   if(line(1:14)=="_cell_length_c")then
    read(line,*)spam,cell_0(3)
    cycle read_cif
   end if
   if(line(1:17)=="_cell_angle_alpha")then
    read(line,*)spam,cell_0(4)
    cycle read_cif
   end if
   if(line(1:16)=="_cell_angle_beta")then
    read(line,*)spam,cell_0(5)
    cycle read_cif
   end if
   if(line(1:17)=="_cell_angle_gamma")then
    read(line,*)spam,cell_0(6)
    cycle read_cif
   end if
   if(line(1:)==string_stop_head) exit read_cif
   end do read_cif
  call cell(rv,vr,cell_0)
  n_atoms=0
  read_natoms: do
   read(100,'(a)',iostat=ierr) line
   if(ierr/=0) exit read_natoms
   n_atoms=n_atoms+1
  end do read_natoms
  rewind(100)
  write(6,'(80a)')('=',l=1,80)
  write(6,*) 'Atoms:', n_atoms,CIFFilename
  write(6,*) 'Atoms/linker', (n_atoms-n_nodes)/n_linkers
  do
   read(100,'(a)') line 
   if(line(1:)==string_stop_head) exit
  end do
  do j=1,n_nodes
   z=z+1
   nodes(z)%id=z 
   nodes(z)%code = nodes_code
   nodes(z)%population=1.0
   call check_cluster_type( nodes(z)%code , nodes(z)%n_components )
   nodes(z)%virtual=.false.
   write(6,'(a)')'=============='
   write(6,'(a,1x,i3,1x,a,1x,a)')"Node: ",z,"Code:", nodes(z)%code
   write(6,'(a,1x,i3)')"# of components:", nodes(z)%n_components
   do k=1,nodes(z)%n_components
    read(100,*) nodes(z)%component_label(k),&
     ( nodes(z)%component_xcrystal(m,k),m=1,3),rrr
    nodes(z)%component_dcrystal(1,k) = rv(1,1)*nodes(z)%component_xcrystal(1,k) + &
                                       rv(1,2)*nodes(z)%component_xcrystal(2,k) + &
                                       rv(1,3)*nodes(z)%component_xcrystal(3,k)
    nodes(z)%component_dcrystal(2,k) = rv(2,1)*nodes(z)%component_xcrystal(1,k) + &
                                       rv(2,2)*nodes(z)%component_xcrystal(2,k) + &
                                       rv(2,3)*nodes(z)%component_xcrystal(3,k)
    nodes(z)%component_dcrystal(3,k) = rv(3,1)*nodes(z)%component_xcrystal(1,k) + &
                                       rv(3,2)*nodes(z)%component_xcrystal(2,k) + &
                                       rv(3,3)*nodes(z)%component_xcrystal(3,k)
    call check_atom_type(nodes(z)%component_label(k),nodes(z)%component_size(k),nodes(z)%element(k))
   end do
   call com_cluster( nodes(z), atom )
   nodes(z)%com(1)=atom(1)
   nodes(z)%com(2)=atom(2)
   nodes(z)%com(3)=atom(3)
   write(6,'(a,1x,3(f14.7,1x))') "Centre of mass:",( nodes(z)%com(k), k=1,3 )
  end do
  z=0
  do j=1,n_linkers
   h=h+1                           ! count total number of linkers
   linkers(h)%id=h
   linkers(h)%code=code
   linkers(h)%population=molar_fraction
   call check_cluster_type( linkers(h)%code , linkers(h)%n_components )
   linkers(h)%virtual=.true.
   ensemble(j,i)=linkers(h)%id     ! registro de linker en file
   write(6,'(a)')'=============='
   write(6,'(a,1x,i3,1x,a,1x,a)')"Ligand:",h,"Code:", linkers(h)%code
   write(6,'(a,1x,i3)')"# of components:", linkers(h)%n_components
   do k=1,linkers(h)%n_components
    read(100,*) linkers(h)%component_label(k),&
     ( linkers(h)%component_xcrystal(m,k),m=1,3),rrr
    linkers(h)%component_dcrystal(1,k) = rv(1,1)*linkers(h)%component_xcrystal(1,k) + &
                                         rv(1,2)*linkers(h)%component_xcrystal(2,k) + &
                                         rv(1,3)*linkers(h)%component_xcrystal(3,k)
    linkers(h)%component_dcrystal(2,k) = rv(2,1)*linkers(h)%component_xcrystal(1,k) + &
                                         rv(2,2)*linkers(h)%component_xcrystal(2,k) + &
                                         rv(2,3)*linkers(h)%component_xcrystal(3,k)
    linkers(h)%component_dcrystal(3,k) = rv(3,1)*linkers(h)%component_xcrystal(1,k) + &
                                         rv(3,2)*linkers(h)%component_xcrystal(2,k) + &
                                         rv(3,3)*linkers(h)%component_xcrystal(3,k)
    call check_atom_type(linkers(h)%component_label(k),linkers(h)%component_size(k),linkers(h)%element(k))
   end do
   call com_cluster( linkers(h), atom )
   linkers(h)%com(1)=atom(1)
   linkers(h)%com(2)=atom(2)
   linkers(h)%com(3)=atom(3)
   write(6,'(a,1x,3(f14.7,1x))') "Centre of mass:",( linkers(h)%com(k), k=1,3 )
  end do
  !h=0
  close(100)
 end do
 close(111)
 call system('rm list tmp')
 ! Detect the smallest ligand type:
 ii=500
 h=1
 do j=1,linker_type_number
  call check_cluster_type( linker_type(j), jj )
  if( jj < ii ) then
   h  = j
   ii = jj
  end if
 end do
 j=0
 allocate(genome(n_linkers))
 genome=0 ! lista de linkers
 nn=0
 write(6,'(a,a,a)') 'Allocating the structure with the smallest available lingand type (',linker_type( h ),'):'
 write(6,'(a)') '% Completed        Configuration'
 linkers%virtual = .true.
 add_linkers: do 
  if(j==n_linkers) exit add_linkers
  j=j+1
  if(nn>=max_number_tries) then
   j=0       ! if the number of tries is bigger than a limit, 
   genome=0  ! we start again with the first ligand.
   nn=0
   linkers%virtual=.true. ! all the linkers are virtuals again!
   write(6,'(a)')'[warning] Relocation this ligands is hard, check the input'
   cycle add_linkers
  end if
  flag = .true.
  do while ( flag )
   k=ensemble(j,randint(1_i64,n_files))
   flag = .false.
   if ( linkers(k)%code /= linker_type( h ) ) flag = .true.
   if ( linkers(k)%virtual .eqv. .false. )    flag = .true.
  end do
  if(j==1)then
   genome(1)=k
   linkers(k)%virtual=.false.
   nn=nn+1
   write(6,*)100*j/real(n_linkers),genome(j),nn, linkers(k)%code, linkers(k)%virtual
   cycle add_linkers
  end if
  genome(j)=k
  linkers(k)%virtual=.false.
  write(6,*)100*j/real(n_linkers),genome(j),nn, linkers(k)%code, linkers(k)%virtual
  nn= 0
 end do add_linkers
 call cluster_connectivity(a)
 write(6,'(a)')"Connectivity table:"
 do ii=1,n_nodes
  h=0
  do jj=1,n_linkers
   if( a(ii,jj) ) h=h+1
  end do
  write(6,*)( a(ii,jj), jj=1, n_linkers ), ":", h
 end do
 write(6,'(a)') "Scanning possible overlaps among atoms of different linkers [...]"
 if ( overlap() ) write(6,'(a)') "[Warnning] Overlap in initial configuration."
 write(6,'(80a)')('=',l=1,80)
 write(6,'(20(a5,1x))')( linker_type(i),i=1,linker_type_number )
 write(6,'(20(f10.8,1x))') ( linker_type_molar_fraction(i), i=1,linker_type_number )
 write(6,'((1000(i6,1x)))') ( genome(i),i=1,n_linkers )
 write(6,'(80a)')('=',l=1,80)
 call update_comfortably()
 mc_steps=0
 call writeCIFFile_from_clusters()
 rrr=0.0
 mc_max_steps=1000*n_linkers
 write(6,'(a,1x,i8)') "Performing MC, max number of steps:", mc_max_steps
 open(678,file="log.txt")
 cell_2(1:6)=cell_0(1:6)
 cell_1(1:6)=cell_0(1:6)
 mc_exchange_linkers: do mc_steps=1, mc_max_steps
  l=randint(1_i64,n_linkers)               ! l:= position of linker in genome
  j=genome( l )                            ! j:= old linker
  if(mc_steps<=100*n_linkers)then
   xyz=0.0_r64
  else
   !xyz=0.1_r64
   xyz=r4_uniform(0.0_r64,1.0_r64)
  end if
  if( xyz < 0.3333333333 ) then
   ! Apply Exchange Move
   n_exchanged_try = n_exchanged_try + 1
   allocate(eee(0:1))
   eee(0)=cost_per_linker(l) + cost_molar()
   k=ensemble(l,randint(1_i64,n_files))     ! k:= new posible linker
   do while (k==j)
    k=ensemble(l,randint(1_i64,n_files))
   end do
   genome(l) = k                     ! test linker-k
   linkers(j)%virtual=.true.         ! j <- virtual
   linkers(k)%virtual=.false.        ! k <- real
   eee(1) = cost_per_linker(l) + cost_molar()
   if( eee(1) >= eee(0)  )then
    genome(l) = j                    ! rechazo el cambio
    linkers(k)%virtual=.true.        ! k <- virtual
    linkers(j)%virtual=.false.       ! j <- real
   else
    n_exchanged=n_exchanged+1
    n_accepted=n_accepted + 1
    ppp=cost()/real_n_atoms()
    if(ppp-rrr<0.0) then
     write(6,'((i5,1x,e20.10,1x,e20.10,1x,e20.10,1x,a,1x,f14.7,1x,a,1000(f14.7,1x)))')&
     mc_steps,ppp,ppp-rrr,cost_molar()/real_n_atoms(),&
     'exchange move:',n_exchanged/real(n_exchanged_try),&
     'molar fractions:',(histogram_molar_fraction(i)/real(n_linkers),i=1,linker_type_number),&
     n_accepted/real(mc_steps)
     write(678,'((i5,1x,e20.10,1x,e20.10,1x,e20.10,1x,a,1000(f14.7,1x)))')&
     mc_steps,ppp,ppp-rrr,cost_molar()/real_n_atoms(),&
    'molar fractions:',(histogram_molar_fraction(i)/real(n_linkers),i=1,linker_type_number)
    end if
    rrr=ppp
    call com_cluster( linkers(j), atom )
    call writeCIFFile_from_clusters()
   end if
   deallocate( eee )
  else if ( xyz < 0.6666666666666 ) then
   n_expanded_try = n_expanded_try + 1
   allocate(eee(0:1))
   cell_1(1:6)=cell_0(1:6)
   eee(0)=cost()/real_n_atoms()
   call cell_expansion()
   call cell(rv,vr,cell_0)
   eee(1)=cost()/real_n_atoms()
   if( eee(1) >= eee(0)  )then
    cell_0(1:6)=cell_1(1:6)
    call cell(rv,vr,cell_0)
   else
    n_expanded=n_expanded+1
    n_accepted=n_accepted + 1
    ppp=cost()/real_n_atoms()
    if(ppp-rrr<0.0) then
     ss = "+"
     if(cell_0(1)<=cell_1(1)) ss = "-"
     write(6,'((i5,1x,e20.10,1x,e20.10,1x,e20.10,1x,a,a,a,1x,f14.7,1x,a,1000(f14.7,1x)))')&
     mc_steps,ppp,ppp-rrr,cost_molar()/real_n_atoms(),&
     'volume ',ss,' move:',n_expanded/real(n_expanded_try),&
     'molar fractions:',(histogram_molar_fraction(i)/real(n_linkers),i=1,linker_type_number),&
     n_accepted/real(mc_steps)
     write(678,'((i5,1x,e20.10,1x,e20.10,1x,e20.10,1x,a,1000(f14.7,1x)))')&
     mc_steps,ppp,ppp-rrr,cost_molar()/real_n_atoms(),&
    'molar fractions:',(histogram_molar_fraction(i)/real(n_linkers),i=1,linker_type_number)
    end if
    rrr=ppp
    call writeCIFFile_from_clusters()
   end if
   deallocate( eee )
  else if ( xyz < 1.0 ) then
   allocate(eee(0:1))
   allocate( backup_positions_linker(1:3, 1: linkers(j)%n_components ) )
   n_rotated_try = n_rotated_try + 1
   eee(0)=cost_per_linker(l) + cost_molar()
   backup_positions_linker(1:3,1:linkers(j)%n_components) = &
    linkers(j)%component_xcrystal(1:3,1:linkers(j)%n_components)
   !call Random_Rotation_Axe_NN_atoms( n_nodes, n_linkers, a, l , wrong)
   call ArbitraryRotateLinker( l , wrong)
   eee(1)=cost_per_linker(l) + cost_molar()
! {{ debugging: refuse always:
!   eee(1)=eee(0)+1.0_r64
! }}
   if( eee(1) >= eee(0) .or. wrong  )then                       ! Rechazo el cambio
    do i=1, linkers(j)%n_components
     do kk=1,3
      linkers(j)%component_xcrystal(kk,i) = backup_positions_linker(kk,i)
     end do
     linkers(j)%component_dcrystal(1,i) = rv(1,1)*linkers(j)%component_xcrystal(1,i) + &
                                          rv(1,2)*linkers(j)%component_xcrystal(2,i) + &
                                          rv(1,3)*linkers(j)%component_xcrystal(3,i)
     linkers(j)%component_dcrystal(2,i) = rv(2,1)*linkers(j)%component_xcrystal(1,i) + &
                                          rv(2,2)*linkers(j)%component_xcrystal(2,i) + &
                                          rv(2,3)*linkers(j)%component_xcrystal(3,i)
     linkers(j)%component_dcrystal(3,i) = rv(3,1)*linkers(j)%component_xcrystal(1,i) + &
                                          rv(3,2)*linkers(j)%component_xcrystal(2,i) + &
                                          rv(3,3)*linkers(j)%component_xcrystal(3,i)
    end do
   else
    n_rotated=n_rotated + 1
    n_accepted=n_accepted + 1
    ppp=cost()/real_n_atoms()
    if(ppp-rrr<0.0) then
     write(6,'((i5,1x,e20.10,1x,e20.10,1x,e20.10,1x,a,1x,f14.7,1x,a,1000(f14.7,1x)))')&
     mc_steps,ppp,ppp-rrr,cost_molar()/real_n_atoms(),&
     'rotation move:',n_rotated/real(n_rotated_try),&
     'molar fractions:',(histogram_molar_fraction(i)/real(n_linkers),i=1,linker_type_number),&
     n_accepted/real(mc_steps)
     write(678,'((i5,1x,e20.10,1x,e20.10,1x,e20.10,1x,a,1000(f14.7,1x)))')&
     mc_steps,ppp,ppp-rrr,cost_molar()/real_n_atoms(),&
    'molar fractions:',(histogram_molar_fraction(i)/real(n_linkers),i=1,linker_type_number)
    end if
    rrr=ppp
    call com_cluster( linkers(j), atom )
    call writeCIFFile_from_clusters()
   end if
   deallocate( backup_positions_linker )
   deallocate( eee )
  else
   write(6,*) xyz
   STOP "Wrong Probabilities"
  end if
 end do mc_exchange_linkers
 cell_0(1:6)=cell_2(1:6)
 write(6,'((i5,1x,e20.10,1x,e20.10,1x,e20.10,1x,a,1000(f14.7,1x)))')&
  mc_max_steps,ppp,ppp-rrr,cost_molar()/real_n_atoms(),&
  'molar fractions:',(histogram_molar_fraction(i)/real(n_linkers),i=1,linker_type_number)
 write(678,'((i5,1x,e20.10,1x,e20.10,1x,e20.10,1x,a,1000(f14.7,1x)))')&
  mc_max_steps,ppp,ppp-rrr,cost_molar()/real_n_atoms(),&
  'molar fractions:',(histogram_molar_fraction(i)/real(n_linkers),i=1,linker_type_number)
 write(6,'(1000(i6,1x))')(genome(i),i=1,n_linkers)
! finish program
 close(678)
 if ( overlap() ) write(6,'(a)') "[Warnning] Overlap in last configuration."
 stop "Finish program"
 contains
!
 logical function overlap()
  implicit none
  integer(i64) :: i,j,ii,jj,kk
  real(r64)    :: rrr,r3(1:3)
  overlap = .false.
  scan_overlap: do i=1,n_linkers
   do j=i+1,n_linkers
    do ii=1,linkers(genome(i))%n_components
     do jj=1,linkers(genome(j))%n_components
      do kk=1,3
       atom(kk)=linkers(genome(i))%component_xcrystal(kk,ii)
       ouratom(kk)=linkers(genome(j))%component_xcrystal(kk,jj)
      end do
      call make_distances(.false.,cell_0,ouratom,atom,rv,r3,rrr)
      if( rrr <= r_min_criteria_overlap ) then
       overlap=.true.
       exit scan_overlap
      end if
     end do
    end do
   end do
  end do scan_overlap
  return
 end function overlap
!
subroutine com_cluster( trozo, r )
 implicit none
 type(cluster),intent(inout) :: trozo
 integer(i64)                :: j,k
 real(r64),intent(out)       :: r(1:3)
 real(r64)                   :: atom(3),ouratom(3),x1
 real(r64),allocatable       :: cluster_x(:,:),cluster_d(:,:)
 allocate(  cluster_x(1:3,1:trozo%n_components ) )
 allocate(  cluster_d(1:3,1:trozo%n_components ) )
 forall ( k=1:3, j=1:trozo%n_components ) 
  cluster_x(k,j)=trozo%component_xcrystal(k,j)
 end forall
 do j=2,trozo%n_components
  forall ( k=1:3 )
   atom(k)   = cluster_x(k,j)
   ouratom(k)= cluster_x(k,j-1) ! { atomo pivote }
  end forall
  call make_distances(.true.,cell_0,atom,ouratom,rv,r,x1)
  forall ( k=1:3 )
   cluster_x(k,j) = r(k)
  end forall
 end do
 forall ( k=1:3 , j=1:trozo%n_components )
  cluster_d(k,j)=rv(k,1)*cluster_x(1,j)+rv(k,2)*cluster_x(2,j)+rv(k,3)*cluster_x(3,j)
 end forall
 do k=1,3
  r(k) = sum( cluster_d(k,1:trozo%n_components) )/real( trozo%n_components )
 end do
 atom(1:3)=r(1:3)
 forall ( k=1:3)
  r(k) = vr(k,1)*atom(1)+vr(k,2)*atom(2)+vr(k,3)*atom(3)
 end forall
 trozo%com(1:3)=r(1:3)
 deallocate( cluster_x, cluster_d )
 return
end subroutine com_cluster
!
real(r64) function min_distance_among_clusters()
 implicit none 
 integer(i64) :: i,j,k
 real(r64)    :: atom(3),ouratom(3),r(3),s
 min_distance_among_clusters=9999999999.999
 do i=1,n_nodes
  do j=i+1,n_nodes
   call com_cluster(nodes(i),atom)
   call com_cluster(nodes(j),ouratom)
   !forall (k=1:3)
   ! atom(k)=nodes(i)%com(k)
   ! ouratom(k)=linkers(j)%com(k)
   !end forall
   call make_distances(.false.,cell_0,atom,ouratom,rv,r,s)
   if( s < min_distance_among_clusters ) then
    min_distance_among_clusters = s
   end if
  end do
 end do
 min_distance_among_clusters = 0.56 + min_distance_among_clusters
 return
end function  min_distance_among_clusters
!
subroutine cluster_connectivity(a)
 implicit none
 logical,intent(out) :: a(n_nodes,n_linkers)
 integer(i64)        :: i,j,k
 real(r64)           :: r(3), s, atom(3), ouratom(3), ss
 ss = min_distance_among_clusters()
 a=.false.
 do i=1,n_nodes
  do j=1,n_linkers
   call com_cluster(nodes(i),atom)
   call com_cluster(linkers(j),ouratom)
   !nodes(i)%com(1:3) = atom(1:3)
   !linkers(j)%com(1:3)=ouratom(1:3)
   call make_distances(.false.,cell_0,atom,ouratom,rv,r,s)
   if( s < ss ) then
    a(i,j)=.true.
   end if
  end do
 end do 
 return
end subroutine
!
 subroutine linker2axe( n_nodes,n_linkers,a, linker, r, v )
  implicit none
  integer(i64),intent(in)  :: n_nodes, n_linkers, linker
  logical,     intent(in)  :: a(n_nodes,n_linkers)
  real(r64), intent(out)   :: v(1:3), r(1:2,1:3)
  real(r64)                :: absvec,atom(3),ouratom(3),r1(3),s
  integer(i64)             :: counter, ii, k,h,jj
  counter = 0
  !call  cluster_connectivity(a)
  found_nodes: do ii=1,n_nodes
   if( a( ii, linker ) ) then
    counter=counter+1
    forall (k=1:3) 
     r(counter, k) = nodes(ii)%com(k)
    end forall
    if (counter==2) exit found_nodes
   end if
  end do found_nodes
  atom(1:3)=r(1,1:3)
  ouratom(1:3)=r(2,1:3)
  call make_distances(.true.,cell_0,atom,ouratom,rv,r(2,1:3),absvec)
  absvec = sqrt( (r(1,1)-r(2,1))**2 + (r(1,2)-r(2,2))**2 + (r(1,3)-r(2,3))**2 )
  do k=1,3
   v(k) = ( r(1, k ) - r(2, k) ) / absvec
  end do
  if(counter/=2) then
   write(6,*) "counter:", counter
   write(6,'(a)')"Connectivity table:"
   do ii=1,n_nodes
    h=0
    do jj=1,n_linkers
     if( a(ii,jj) ) h=h+1
    end do
    write(6,*)( a(ii,jj), jj=1, n_linkers ), ":", h
   end do
   stop "WRONG connectivity"
  end if
  return
 end subroutine linker2axe
!
 subroutine cell_expansion()
  implicit none
  real(r64) :: r
  r = r4_uniform(0.95_r64, 1.05_r64)
  cell_0(1:3)=cell_0(1:3)*r
  return
 end subroutine cell_expansion
!
 subroutine linker2points( linker, r1, r2 )
  implicit none
  integer(i64),intent(in)  :: linker
  real(r64),   intent(out) :: r1(1:3), r2(1:3)
  real(r64)                :: s,t, atom(3), ouratom(3), r(1:2,1:3)
  integer(i64)             :: counter, i, k, h, ii, jj
  counter = 0_i64
  t = min_distance_among_clusters() 
  write(6,*)" "
  call com_cluster(linkers(linker),ouratom)
  write(6,*) linkers(linker)%com(1:3)
  !linkers(j)%com(1:3)=ouratom(1:3)
  found_nodes_: do i=1,n_nodes
   call com_cluster(nodes(i),atom)
   write(6,*) nodes(i)%com(1:3)
   !nodes(i)%com(1:3) = atom(1:3)
   call make_distances(.false.,cell_0,atom,ouratom,rv,r,s)
   write(6,*) i,linker,s,t
   if( s <= t ) then
    counter = counter + 1
    do k=1,3
     r(counter,k) = nodes(i)%com(k)
    end do
    if (counter==2) exit found_nodes_
   end if
  end do found_nodes_
  r1(1:3)=r(1,1:3)
  r2(1:3)=r(2,1:3)
  atom(1:3)=r1(1:3)
  ouratom(1:3)=r2(1:3)
  call make_distances(.true.,cell_0,atom,ouratom,rv,r2,s)
  if(counter/=2) then
   write(6,*) "counter:", counter
   call cluster_connectivity(a)
   write(6,'(a)')"Connectivity table:"
   do ii=1,n_nodes
    h=0
    do jj=1,n_linkers
     if( a(ii,jj) ) h=h+1
    end do
    write(6,*)( a(ii,jj), jj=1, n_linkers ), ":", h
   end do
   stop "WRONG connectivity"
  end if
  return
 end subroutine linker2points
!
 subroutine ArbitraryRotateLinker(linker, wrong)
  implicit none
  integer(i64), intent(in)  :: linker
  logical,      intent(out) :: wrong
  integer(i64)              :: i,j,k
  real(r64), parameter      :: MaximumRotation = 45.0_r64
  real(r64), parameter      :: pi = ACOS(-1.0_r64)
  real(r64)                 :: angle, rot(3,3), r1(3), r2(3), x,y,z,r,s,c,trace
  real(r64)                 :: u,v,w,o,p,q,absvec, DEGTORAD
  real(r64)                 :: uv(3),ori(3)
  !write(6,*) "Enter in rotation subroutine", linker
  DEGTORAD = pi/180.0_r64
  wrong=.false.
  angle = MaximumRotation * r4_uniform(0.0_r64, 1.0_r64)
  !write(6,*) "Angle:", angle
  s = sin( angle * DEGTORAD )
  c = cos( angle * DEGTORAD )
  call linker2points( linker, r1, r2 )
! r1 : reciprocal space, (o,p,q) : direct space
  o = r1(1)*rv(1,1) + r1(2)*rv(1,2) + r1(3)*rv(1,3)
  p = r1(1)*rv(2,1) + r1(2)*rv(2,2) + r1(3)*rv(2,3)
  q = r1(1)*rv(3,1) + r1(2)*rv(3,2) + r1(3)*rv(3,3)
! r2 : reciprocal space, (u,v,w) : direct space
  u = r2(1)*rv(1,1) + r2(2)*rv(1,2) + r2(3)*rv(1,3)
  v = r2(1)*rv(2,1) + r2(2)*rv(2,2) + r2(3)*rv(2,3)
  w = r2(1)*rv(3,1) + r2(2)*rv(3,2) + r2(3)*rv(3,3)
  !write(6,'(3(f14.7,1x))')o,p,q
  !write(6,'(3(f14.7,1x))')u,v,w
  call com_cluster( linkers(linker), ori )
  x = ori(1)
  y = ori(2)
  z = ori(3)
  ori(1) = x*rv(1,1) + y*rv(1,2) + z*rv(1,3)
  ori(2) = x*rv(2,1) + y*rv(2,2) + z*rv(2,3)
  ori(3) = x*rv(3,1) + y*rv(3,2) + z*rv(3,3)
  !write(6,'(3(f14.7,1x))')(ori(k),k=1,3)
! Unit vector in direct space
  x = u - o
  y = v - p
  z = w - q
  absvec = sqrt( x*x + y*y + z*z)
  x = x/absvec
  y = y/absvec
  z = z/absvec
!
  !write(6,'(4(f14.7,1x))')x,y,z,x*x + y*y + z*z
!
  rot(1,1)=c + (1.0_r64-c)*x*x
  rot(1,2)=(1.0_r64-c)*x*y - z*s
  rot(1,3)=(1.0_r64-c)*x*z + y*s
  rot(2,1)=(1.0_r64-c)*x*y + z*s
  rot(2,2)=c + (1.0_r64-c)*y*y
  rot(2,3)=(1.0_r64-c)*y*z - x*s
  rot(3,1)=(1.0_r64-c)*x*z - y*s
  rot(3,2)=(1.0_r64-c)*y*z + x*s
  rot(3,3)=c + (1.0_r64-c)*z*z
  trace = rot(1,1)+rot(2,2)+rot(3,3)
! {{ debuging, write matrix:
   !write(6,'(a)')"Rotation Matrix:"
   !write(6,'(a,1x,3(f14.7,1x),a)') "/",( rot(1,ii), ii=1,3),"\"
   !write(6,'(a,1x,3(f14.7,1x),a)') "|",( rot(2,ii), ii=1,3),"|"
   !write(6,'(a,1x,3(f14.7,1x),a)') "\",( rot(3,ii), ii=1,3),"/"
   !write(6,'(a)') " "
   !write(6,'(a)')"Transformation Matrix:"
   !write(6,'(a,1x,3(f14.7,1x),a)') "/",( rv(1,ii), ii=1,3),"\"
   !write(6,'(a,1x,3(f14.7,1x),a)') "|",( rv(2,ii), ii=1,3),"|"
   !write(6,'(a,1x,3(f14.7,1x),a)') "\",( rv(3,ii), ii=1,3),"/"
!  debugging, check Trace of the transformation
   !write( 6,*) "Traze:", trace
  if( abs( trace - (1.0_r64+2.0_r64*c) ) >= 1.0e-6_r64 ) STOP "Wrong Rotation Matrix"
  do i=1,linkers(linker)%n_components
! (x,y,z): direct space, (o,p,q): direct space
   x=linkers(linker)%component_xcrystal(1,i)*rv(1,1)+&
     linkers(linker)%component_xcrystal(2,i)*rv(1,2)+&
     linkers(linker)%component_xcrystal(3,i)*rv(1,3)
   y=linkers(linker)%component_xcrystal(1,i)*rv(2,1)+&
     linkers(linker)%component_xcrystal(2,i)*rv(2,2)+&
     linkers(linker)%component_xcrystal(3,i)*rv(2,3)
   z=linkers(linker)%component_xcrystal(1,i)*rv(3,1)+&
     linkers(linker)%component_xcrystal(2,i)*rv(3,2)+&
     linkers(linker)%component_xcrystal(3,i)*rv(3,3)
!  Move to the origin
   x=x-ori(1)
   y=y-ori(2)
   z=z-ori(3)
!  Apply rotation in direct space
   o=x*rot(1,1)+y*rot(1,2)+z*rot(1,3) 
   p=x*rot(2,1)+y*rot(2,2)+z*rot(2,3)
   q=x*rot(3,1)+y*rot(3,2)+z*rot(3,3)
!  Translate axis and rotated point back to original location
   linkers(linker)%component_dcrystal(1,i) = o + ori(1)
   linkers(linker)%component_dcrystal(2,i) = p + ori(2)
   linkers(linker)%component_dcrystal(3,i) = q + ori(3)
!  Go back to the reciprocal space
   linkers(linker)%component_xcrystal(1,i) =  vr(1,1)*linkers(linker)%component_dcrystal(1,i)+&
                                              vr(1,2)*linkers(linker)%component_dcrystal(2,i)+&
                                              vr(1,3)*linkers(linker)%component_dcrystal(3,i)
   linkers(linker)%component_xcrystal(1,i) =  vr(2,1)*linkers(linker)%component_dcrystal(1,i)+&
                                              vr(2,2)*linkers(linker)%component_dcrystal(2,i)+&
                                              vr(2,3)*linkers(linker)%component_dcrystal(3,i)
   linkers(linker)%component_xcrystal(1,i) =  vr(3,1)*linkers(linker)%component_dcrystal(1,i)+&
                                              vr(3,2)*linkers(linker)%component_dcrystal(2,i)+&
                                              vr(3,3)*linkers(linker)%component_dcrystal(3,i)
  end do
  call com_cluster( linkers(linker), atom )
  o=min_distance_among_clusters()
  detect_wrong_rotation: do i=1,2
   if ( i==1 ) call make_distances(.false.,cell_0,atom,r1,rv,ouratom,s)
   if ( i==2 ) call make_distances(.false.,cell_0,atom,r2,rv,ouratom,s)
   if ( s > o ) then
    wrong = .true.
    !write(6,*) i, s, o
    !write(6,'(a)')"Terrible Rotation detected !!"
    exit detect_wrong_rotation
   end if
  end do detect_wrong_rotation
  !linkers(linker)%element(1:linkers(linker)%n_components)="S "
  !linkers(linker)%component_label(1:linkers(linker)%n_components)="S "
  return
 end subroutine ArbitraryRotateLinker
!
 subroutine Random_Rotation_Axe_NN_atoms(n_nodes,n_linkers,a,linker,wrong)
  ! https://en.wikipedia.org/wiki/Rotation_matrix#cite_ref-5
  implicit none
  integer(i64),intent(in)   ::   n_nodes,n_linkers,linker
  integer                   ::   w=900
  logical,     intent(in)   ::   a(n_nodes,n_linkers)
  integer(i64)              ::   ii,jj,ll
  real(r64)                 ::   MaximumRotation = 45.0_r64
  real(r64), parameter      ::   pi = ACOS(-1.0_r64)
  real(r64)                 ::   DEGTORAD
  real(r64)                 ::   rot(1:3,1:3), angle_change, rr(1:2,1:3)
  real(r64)                 ::   u(1:3),x,y,z,r,s,c,o,p,q, trace
  logical,intent(out)       ::   wrong
! {{ debugging:
  !open( w,file="rotation.cif")
  !write(w,'(a)')'data_cif'
  !write(w,'(a)')"_audit_creation_method 'generated by iGOR'"
  !write(w,'(a)')"_symmetry_space_group_name_H-M     'P 1'"
  !write(w,'(a)')'_symmetry_Int_Tables_number        1'
  !write(w,'(a)')'_symmetry_cell_setting             triclinic'
  !write(w,'(a,f14.7)')'_cell_length_a                 ',cell_0(1)
  !write(w,'(a,f14.7)')'_cell_length_b                 ',cell_0(2)
  !write(w,'(a,f14.7)')'_cell_length_c                 ',cell_0(3)
  !write(w,'(a,f14.7)')'_cell_angle_alpha              ',cell_0(4)
  !write(w,'(a,f14.7)')'_cell_angle_beta               ',cell_0(5)
  !write(w,'(a,f14.7)')'_cell_angle_gamma              ',cell_0(6)
  !write(w,'(a)')'loop_'
  !write(w,'(a)')'_atom_site_type_symbol'
  !write(w,'(a)')'_atom_site_label'
  !write(w,'(a)')'_atom_site_fract_x'
  !write(w,'(a)')'_atom_site_fract_y'
  !write(w,'(a)')'_atom_site_fract_z'
  !write(w,'(a)')'_atom_site_symmetry_multiplicity'
! }}
  DEGTORAD = pi/180.0_r64
  wrong=.false.
  !write(6,'(a)') "================================="
  !write(6,'(a,1x,i5)')"Selected ligand:", linker
  !write(6,*) ( linkers(linker)%com(ii), ii=1,3)  
  !do while(wrong)
  r = r4_uniform(-1.0_r64, 1.0_r64)
  !r = 1.0_i64
  angle_change=MaximumRotation*r
  if ( abs(angle_change - 90.0_r64) <= 1.0e-2 .or. abs(angle_change + 90.0_r64) <= 1.0e-2 ) then
   s = 1.0
   c = 0.0
  else
   s = sin( angle_change*DEGTORAD )
   c = cos( angle_change*DEGTORAD )
  end if
  call linker2axe( n_nodes, n_linkers, a, linker, rr, u )
! {{ debugging:
  !write(w,'(a,1x,a,1x,3(f14.7,1x),i4)') "Xe","Xe",( rr(1,ii), ii=1,3), 0
  !do ll=1,10
   !write(w,'(a,1x,a,1x,3(f14.7,1x),i4)') "H","H",rr(1,1)+ll*(rr(1,1)-rr(2,1))/10.0,&
   ! rr(1,2)+ll*(rr(1,2)-rr(2,2))/10.0,rr(1,3)+ll*(rr(1,3)-rr(2,3))/10.0,0
  !end do
  !write(w,'(a,1x,a,1x,3(f14.7,1x),i4)') "Xe","Xe",( rr(2,ii), ii=1,3), 0
! }}
  x=u(1)
  y=u(2)
  z=u(3)
! {{ debugging
   !write(6,*)'Rotation:'
   !write(6,*)'Axe:', x,y,z, x**2+y**2+z**2
   !write(6,*)'Angle:',MaximumRotation*r,'deg'
! }}
  rot(1,1)=c + x*x*(1-c)
  rot(1,2)=x*y*(1-c) - z*s
  rot(1,3)=x*z*(1-c) + y*s
  rot(2,1)=y*x*(1-c) + z*s
  rot(2,2)=c + y*y*(1-c)
  rot(2,3)=y*z*(1-c) - x*s
  rot(3,1)=z*x*(1-c) - y*s
  rot(3,2)=z*y*(1-c) + x*s
  rot(3,3)=c + z*z*(1-c)
  trace = rot(1,1)+rot(2,2)+rot(3,3)
! {{ debuging, write matrix:
   !write(6,'(a)')"Rotation Matrix:"
   !write(6,'(a,1x,3(f14.7,1x),a)') "/",( rot(1,ii), ii=1,3),"\"
   !write(6,'(a,1x,3(f14.7,1x),a)') "|",( rot(2,ii), ii=1,3),"|"
   !write(6,'(a,1x,3(f14.7,1x),a)') "\",( rot(3,ii), ii=1,3),"/"
   !write(6,'(a)') " "
   !write(6,'(a)')"Transformation Matrix:"
   !write(6,'(a,1x,3(f14.7,1x),a)') "/",( rv(1,ii), ii=1,3),"\"
   !write(6,'(a,1x,3(f14.7,1x),a)') "|",( rv(2,ii), ii=1,3),"|"
   !write(6,'(a,1x,3(f14.7,1x),a)') "\",( rv(3,ii), ii=1,3),"/"
!  debugging, check Trace of the transformation
   !write( 6,*) "Traze:", trace
   !write(6,*) "Checking Axe:"
   !if( abs(trace+1) <= 1e-4 .or. s <= 1e-6 ) then ! s = 0
   ! write(6,*) "+/-",sqrt( 0.5*(1_i64 + rot(1,1)) )
   ! write(6,*) "+/-",sqrt( 0.5*(1_i64 + rot(2,2)) )
   ! write(6,*) "+/-",sqrt( 0.5*(1_i64 + rot(3,3)) )
   !else 
   ! write(6,*) (rot(3,2) - rot(2,3) )/ sqrt( (3-trace)*(1+trace) )
    !if( abs( (rot(3,2) - rot(2,3) )/ sqrt( (3-trace)*(1+trace) ) - x) >= 1e-5 ) STOP "Wrong Calculated Axe X"
   ! write(6,*) (rot(1,3) - rot(3,1) )/ sqrt( (3-trace)*(1+trace) )
    !if( abs( (rot(1,3) - rot(3,1) )/ sqrt( (3-trace)*(1+trace) ) - y) >= 1e-5 ) STOP "Wrong Calculated Axe Y"
   ! write(6,*) (rot(2,1) - rot(1,2) )/ sqrt( (3-trace)*(1+trace) )
    !if( abs( (rot(2,1) - rot(1,2) )/ sqrt( (3-trace)*(1+trace) ) - z) >= 1e-5 ) STOP "Wrong Calculated Axe Z"
   !end if
! }}
  if( abs( trace - (1+2*c) ) >= 1e-6 ) STOP "Wrong Rotation Matrix"
! {{
  !write(6,'(a)') "Transformation of coordinates:"
  do ii=1,linkers(linker)%n_components
   !write(6,'(i3,1x,a,1x,3(f14.7,1x))') ii,linkers(linker)%element(ii),&
   ! ( linkers(linker)%component_xcrystal(jj,ii), jj=1,3)
   !write(w,'(a,1x,a,1x,3(f14.7,1x),i4)') "O","O",&
   ! ( linkers(linker)%component_xcrystal(jj,ii), jj=1,3),1
   x=linkers(linker)%component_xcrystal(1,ii)*rv(1,1)+&
     linkers(linker)%component_xcrystal(2,ii)*rv(1,2)+&
     linkers(linker)%component_xcrystal(3,ii)*rv(1,3)
   y=linkers(linker)%component_xcrystal(1,ii)*rv(2,1)+&
     linkers(linker)%component_xcrystal(2,ii)*rv(2,2)+&
     linkers(linker)%component_xcrystal(3,ii)*rv(2,3)
   z=linkers(linker)%component_xcrystal(1,ii)*rv(3,1)+&
     linkers(linker)%component_xcrystal(2,ii)*rv(3,2)+&
     linkers(linker)%component_xcrystal(3,ii)*rv(3,3)
   !write(6,'(i3,1x,a,1x,3(f14.7,1x))') ii,linkers(linker)%element(ii), x, y, z
   o=x*rot(1,1)+y*rot(1,2)+z*rot(1,3)
   p=x*rot(2,1)+y*rot(2,2)+z*rot(2,3)
   q=x*rot(3,1)+y*rot(3,2)+z*rot(3,3)
! 
   linkers(linker)%component_dcrystal(1,ii) = o
   linkers(linker)%component_dcrystal(2,ii) = p
   linkers(linker)%component_dcrystal(3,ii) = q
! 
   linkers(linker)%component_xcrystal(1,ii) = vr(1,1)*o+vr(1,2)*p+vr(1,3)*q
   linkers(linker)%component_xcrystal(2,ii) = vr(2,1)*o+vr(2,2)*p+vr(2,3)*q
   linkers(linker)%component_xcrystal(3,ii) = vr(3,1)*o+vr(3,2)*p+vr(3,3)*q
   !write(6,'(i3,1x,a,1x,3(f14.7,1x))') ii,linkers(linker)%element(ii),&
   ! ( linkers(linker)%component_xcrystal(jj,ii), jj=1,3)
   !write(w,'(a,1x,a,1x,3(f14.7,1x),i3)') "S","S",&
   ! ( linkers(linker)%component_xcrystal(jj,ii), jj=1,3),2
   !write(6,'(i3,1x,a,1x,3(f14.7,1x))') ii,linkers(linker)%element(ii), o, p, q
  end do
  call com_cluster( linkers(linker), atom )
  o=min_distance_among_clusters()
  do ii=1,2
   call make_distances(.false.,cell_0,atom,rr(ii,1:3),rv,ouratom,s) 
   if ( s > o ) then
    wrong = .true.
    !write(6,*) ii, s, o
    !write(6,'(a)')"Terrible Rotation detected !!"
   end if
  end do
  !close(w)
  return
 end subroutine Random_Rotation_Axe_NN_atoms
!
 subroutine update_molar_fraction()
  implicit none
  integer(i64)         :: abc
  integer(i64)         :: iabc,ijk
  histogram_molar_fraction=0
  do iabc=1,n_linkers
   ijk=genome(iabc)
   do abc=1,linker_type_number
    if(linkers(ijk)%code==linker_type(abc))then
     histogram_molar_fraction(abc)=histogram_molar_fraction(abc)+1
    end if
   end do
  end do
  return
 end subroutine update_molar_fraction
!
 real(r64) function cost_expansion()
  implicit none
  integer(i64) :: i
  real(r64)    :: kv = 1.0e9_r64
  cost_expansion = 0.0_r64
  do i=1,3
   cost_expansion=cost_expansion+0.5*kv*(cell_0(i)-cell_2(i))**2
  end do
  return
 end function cost_expansion
!
 real(r64) function cost_molar()
  implicit none
  real(r64), parameter  ::  molar_constant=1.0e9_r64
  integer(i64)          ::  abc
  call update_molar_fraction()
  cost_molar=0.0_r64
  do abc=1,linker_type_number
   cost_molar=cost_molar+&
   molar_constant*(histogram_molar_fraction(abc)/real(n_linkers)-linker_type_molar_fraction(abc))**2
  end do
  return
 end function  cost_molar
!
 integer(i64) function real_n_atoms()
  implicit none
  integer(i64)              :: i
  real_n_atoms = 0
  do i=1,n_linkers
   real_n_atoms = real_n_atoms + linkers(genome(i))%n_components
  end do
  return
 end function real_n_atoms
!
 subroutine update_comfortably()
  ! Well-being allocated linkers in the ensamble:
  implicit none
  integer(i64) :: i
  real(r64)    :: base_comfort = 0.0
  real(r64)    :: comfort(n_linkers),total_comfort
  real(r64)    :: infinity = huge(0.0)
  comfort=0.0
  total_comfort=0.0
  do i=1,n_linkers
   linkers(genome(i))%comfortably = cost_per_linker(i)
   if( linkers(genome(i))%comfortably /= linkers(genome(i))%comfortably ) then
   ! By definition, NAN is not equal to anything, even itself.
    linkers(genome(i))%comfortably=infinity
   end if
   comfort(i) = linkers(genome(i))%comfortably
  end do
  base_comfort = minval(comfort)
  comfort(1:n_linkers)=comfort(1:n_linkers)-base_comfort
  total_comfort= sum(comfort) 
  ! check to see if your variable is infinity by:
  if( total_comfort > infinity ) total_comfort = infinity
  linkers(genome(1:n_linkers))%comfortably = comfort(1:n_linkers)/total_comfort
  !write(6,*)'Base:',base_comfort,'Total comfort:',total_comfort ! debugging
  !do i=1,n_linkers
  ! write(6,*)'Comfort:',i,genome(i),linkers(genome(i))%code,linkers(genome(i))%comfortably,cost_per_linker(i)
  !end do
  return
 end subroutine update_comfortably
!
 subroutine choose_by_molar_fraction(pivot,choosen)
  implicit none
  integer(i64)             :: ii,k
  integer(i64),intent(in)  :: pivot
  integer(i64),intent(out) :: choosen
  real(r64)                :: sum_weight, eta
  real(i64)                :: weight(n_files)
  do ii=1,n_files
   weight(ii)=linkers(ensemble(pivot,ii))%population 
  end do
  sum_weight = sum( weight(1:n_files ) )
  eta = r4_uniform(0.0_r64, sum_weight)
  do choosen = 1, n_files
    if ( eta <  linkers(ensemble(pivot,choosen))%population ) exit
    eta = eta - linkers(ensemble(pivot,choosen))%population
  end do 
  return 
 end subroutine choose_by_molar_fraction
!
 subroutine choose_by_comfort(choosen)
  implicit none
  integer(i64),intent(inout) :: choosen
  integer(i64)               :: ii
  real(r64)                  :: comfort(n_linkers)
  real(r64)                  :: sum_weight, eta
  call update_comfortably()
  sum_weight = sum( linkers(genome(1:n_linkers ))%comfortably  )  
  eta = r4_uniform(0.0_r64, sum_weight)
  do choosen = 1, n_linkers
    if ( eta <  linkers(genome(choosen))%comfortably ) exit
    eta = eta - linkers(genome(choosen))%comfortably
  end do
  if(choosen>n_linkers) stop 'choosen > n_linkers'
  return
 end subroutine choose_by_comfort
!
 real(r64) function cost_per_linker(identi)
  ! 11 may 2017
  implicit none
  integer(i64),intent(in) :: identi
  integer(i64)            :: ii,jj,k ,i
  real(r64)               :: r,r3(1:3)
  real(r64)               :: ell=0.5
  real(r64)               :: rll=3.9
  real(r64)               :: rml=0.0
  real(r64)               :: eml=0.0
  real(r64)               :: infinite = 1.0e34
  cost_per_linker=0.0
  do i=1,n_linkers
   if(i/=identi)then
    do ii=1,linkers(genome(i))%n_components
     do jj=1,linkers(genome(identi))%n_components
      forall (k=1:3)
       atom(k)=linkers(genome(identi))%component_xcrystal(k,jj)
       ouratom(k)=linkers(genome(i))%component_xcrystal(k,ii)
      end forall
      call make_distances(.false.,cell_0,ouratom,atom,rv,r3,r)
      cost_per_linker = cost_per_linker + ell*((rll/r)**12-2*(rll/r)**6)
     end do
    end do
   end if
  end do
  if( cost_per_linker > infinite ) cost_per_linker = infinite
  if( cost_per_linker /= cost_per_linker ) cost_per_linker= infinite
  return
 end function cost_per_linker
!
 real(r64) function cost_exchange()
  implicit none 
  integer(i64)               :: i
  cost_exchange=0.0
  do i=1,n_linkers
   cost_exchange= cost_exchange+cost_per_linker(i)
  end do
  return
 end function cost_exchange
!
 real(r64) function cost()
  implicit none
  real(r64) :: infinite = 1.0e34
  cost = cost_exchange() + cost_molar() + cost_expansion()
  if ( cost > infinite ) cost = infinite
  return
 end function  cost
!
 subroutine check_topology_composition(top,mmm,nnn,lll,node_code)
  implicit none
  character(len=3),intent(in)  :: top
  integer(i64),intent(in)           :: mmm
  integer(i64),intent(out)          :: nnn
  integer(i64),intent(out)          :: lll
  character(len=3),intent(out) :: node_code
  select case(top)
   case("LP1","LP2")   ! IRMOF-10 interprenetrate large pore
    ! Zn4O(BDC)3
    node_code="oxo"
    nnn = int( mmm/4.0 )
    lll = 3 * nnn
   case("CP1","CP2")   ! IRMOF-10 interpenetrate  narrow pore
    ! Zn4O(BDC)3
    node_code="oxo"
    nnn = int( mmm/4.0 )
    lll = 3 * nnn
   case("I10")   ! IRMOF-10 non-interpenetrate cubic form
    ! Zn4O(BDC)3
    node_code="oxo"
    nnn = int( mmm/4.0 )
    lll = 3 * nnn
   case default  ! ZIF with IZA code
    node_code="Zn_"
    ! ZnN2( XIm )2
    !write(6,'(a)') "detecting Zeolitic Imidazolate Framework" 
    nnn = mmm
    lll = 2 * nnn
  end select 
 end subroutine check_topology_composition
!
 subroutine check_cluster_type(lll,nnn)
  implicit none
  character(len=3),intent(in)  :: lll
  integer(i64)    ,intent(out) :: nnn
  check_cluster: select case (lll)
   case("oxo")
    nnn=23
   case("Zn_")
    nnn=1
   case("li1")
    nnn=32
   case("li2")
    nnn=20
   case("li3")
    nnn=84
   case("li4")
    nnn=58
   case("li5")
    nnn=60
   case("li6")
    nnn=72
   case("li7")
    nnn=84
   case("li8")
    nnn=31
   case("li9")
    nnn=31
   case("l10")
    nnn=42
   case("l11")
    nnn=42
   case("l12")
    nnn=40
   case("imi")
    nnn=8
    write(6,'(a)') "detecting imidazol"
   case("mim")
    nnn=11
   case("bim")
    nnn=14
   case default
    nnn=0
    write(6,'(a)') "Linker type does not found"
    stop "[Error] Composition"
  end select check_cluster
  write(6,'(i3,1x,a)') nnn, "atoms per linker"
  return
 end subroutine check_cluster_type
!
 subroutine check_atom_type(lll,s,element)
  implicit none
  real(r64),intent(out)            :: s
  character(len=2),intent(in) :: lll
  character(len=2),intent(out):: element
  check_atom: select case (lll)
   case('H ',' H','HO','H1','H2','H3')
    rrr=0.320
    element="H "
   case('C ','C1':'C9',' C')
    s=0.720
    element="C "
   case('N ',' N','N1':'N9')
    s=0.7
    element="N "
   case('O ','O1':'O9','OH',' O','O_')
    s=0.7
    element="O "
   case('Si')
    s=1.14
    element="Si"
   case('Al')
    s=1.14
    element="Al"
   case('Xe')
    s=0.0001
    element="Xe"
   case('Zn')
    s=1.6
    element="Zn"
   case('Cl')
    s=1.00
    element="Cl"
   case default
    write(6,'(a)')"[Error]"
    write(6,*) lll
    STOP 'Atom unknowed'
  end select check_atom
  return
 end subroutine check_atom_type
!
 subroutine writeCIFFile_from_clusters()
  implicit none
  integer(i64),parameter :: u=234
  integer(i64)           :: j
  character(len=80) :: outfilename ='output.cif'
  character(len=80) :: strfff,lineff 
  write(lineff,'(a1,i1,a)')'(',linker_type_number,'(a3))'
  !write(6,'(a)') line
  write(strfff,adjustl(trim(lineff)))(adjustl(trim(linker_type(ii))), ii=1,linker_type_number )
  write(outfilename,'(a,a,a,a)')adjustl(trim(topology)),'_',adjustl(trim(strfff)),'_mixture.cif' 
  open(u,file=outfilename,iostat=ierr)
  write(u,'(a)')'data_cif'
  write(u,'(a)')"_audit_creation_method 'generated by iGOR'"
  write(u,'(a)')"_symmetry_space_group_name_H-M     'P 1'"
  write(u,'(a)')'_symmetry_Int_Tables_number        1'
  write(u,'(a)')'_symmetry_cell_setting             triclinic'
  write(u,'(a,f14.7)')'_cell_length_a                 ',cell_0(1)
  write(u,'(a,f14.7)')'_cell_length_b                 ',cell_0(2)
  write(u,'(a,f14.7)')'_cell_length_c                 ',cell_0(3)
  write(u,'(a,f14.7)')'_cell_angle_alpha              ',cell_0(4)
  write(u,'(a,f14.7)')'_cell_angle_beta               ',cell_0(5)
  write(u,'(a,f14.7)')'_cell_angle_gamma              ',cell_0(6)
  write(u,'(a)')'loop_'
  write(u,'(a)')'_atom_site_type_symbol'
  write(u,'(a)')'_atom_site_label'
  write(u,'(a)')'_atom_site_fract_x'
  write(u,'(a)')'_atom_site_fract_y'
  write(u,'(a)')'_atom_site_fract_z'
  write(u,'(a)')'_atom_site_symmetry_multiplicity'
  do i=1,n_nodes
   if( nodes(i)%virtual.eqv..false.)then
    do j=1,nodes(i)%n_components
     write(u,'(a2,1x,a2,1x,3(f14.7,1x),i5)')nodes(i)%element(j),nodes(i)%component_label(j),&
      (nodes(i)%component_xcrystal(m,j),m=1,3),nodes(i)%id
    end do
   end if
  end do
  do i=1,n_linkers
   if( linkers(genome(i))%virtual.eqv..false.)then
    do j=1,linkers(genome(i))%n_components
     write(u,'(a2,1x,a2,1x,3(f14.7,1x),i5)')linkers(genome(i))%element(j),&
       linkers(genome(i))%component_label(j),&
      (linkers(genome(i))%component_xcrystal(m,j),m=1,3),linkers(genome(i))%id
    end do
   end if
  end do
  close(u)
 end subroutine

 PURE INTEGER(i64) FUNCTION Clen(s)      ! returns same result as LEN unless:
  CHARACTER(*),INTENT(IN) :: s       ! last non-blank char is null
  INTEGER(i64) :: i
  Clen = LEN(s)
  i = LEN_TRIM(s)
  IF (s(i:i) == CHAR(0)) Clen = i-1  ! len of C string
 END FUNCTION Clen

 PURE INTEGER(i64) FUNCTION Clen_trim(s) ! returns same result as LEN_TRIM unless:
  CHARACTER(*),INTENT(IN) :: s       ! last char non-blank is null, if true:
  INTEGER(i64) :: i                       ! then len of C string is returned, note:
                                     ! Ctrim is only user of this function
  i = LEN_TRIM(s) ; Clen_trim = i
  IF (s(i:i) == CHAR(0)) Clen_trim = Clen(s)   ! len of C string
 END FUNCTION Clen_trim

 SUBROUTINE cell(rv,vr,cell_0)
  implicit none
  integer(i64) :: i,j
  real(r64), intent(in)  :: cell_0(6)
  real(r64), intent(out) :: rv(3,3),vr(3,3)
  real(r64), parameter   :: pi = ACOS(-1.0)
  real(r64) :: alp,bet
  real(r64) :: cosa,cosb,cosg
  real(r64) :: gam,sing
  real(r64) :: DEGTORAD
  DEGTORAD=pi/180.0_r64
  IF(cell_0(4) == 90.0_r64) THEN
    cosa = 0.0_r64
  ELSE
    ALP=cell_0(4)*degtorad
    COSA=cos(ALP)
  ENDIF
  IF(cell_0(5) == 90.0_r64) THEN
    cosb = 0.0_r64
  ELSE
    bet = cell_0(5)*degtorad
    cosb = cos(bet)
  ENDIF
  IF(cell_0(6) == 90.0_r64) then
    sing = 1.0_r64
    cosg = 0.0_r64
  ELSE
    gam = cell_0(6)*degtorad
    sing = sin(gam)
    cosg = cos(gam)
  ENDIF
  rv(1,1) = cell_0(1)
  rv(1,2) = cell_0(2)*cosg
  rv(1,3) = cell_0(3)*cosb
  rv(2,1) = 0.0_r64
  rv(2,2) = cell_0(2)*sing
  rv(2,3) = cell_0(3)*(cosa - cosb*cosg)/sing
  rv(3,1) = 0.0_r64
  rv(3,2) = 0.0_r64
  rv(3,3) = sqrt( cell_0(3)*cell_0(3) - rv(1,3)*rv(1,3) - rv(2,3)*rv(2,3))
  call inverse(rv,vr,3_i64)
  !WRITE(6,'(a)') 'Cell:'
  !WRITE(6,'(6F14.7)')( cell_0(j), j=1,6 )
  !WRITE(6,'(a)')'Linear Transformation Operator:'
  !DO i=1,3
  ! WRITE(6,'(F14.7,F14.7,F14.7)')( rv(i,j), j=1,3 )
  !ENDDO
  !WRITE(6,'(a)')'----------------------------------------'
  !WRITE(6,'(a)')'Inverse Linear Transformation Operator:'
  !DO i=1,3
  ! WRITE(6,'(F14.7,F14.7,F14.7)')( vr(i,j), j=1,3 )
  !ENDDO
  !WRITE(6,'(a)')'----------------------------------------'
  RETURN
 END SUBROUTINE cell
!
 SUBROUTINE uncell(rv,cell_0)
  implicit none
  real(r64),intent(out)   :: cell_0(6)
  real(r64),intent(in)    :: rv(3,3)
  integer(i64)            :: i,j
  real(r64)               :: temp(6)
  REAL(r64)               :: radtodeg
  REAL(r64), PARAMETER    :: pi=ACOS(-1.0) 
  radtodeg=180.0/PI
  do i = 1,3
    temp(i) = 0.0
    do j = 1,3
      temp(i) = temp(i) + rv(j,i)**2
    enddo
    temp(i) = sqrt(temp(i))
  enddo
  cell_0(1) = abs(temp(1))
  cell_0(2) = abs(temp(2))
  cell_0(3) = abs(temp(3))
  do i = 1,3
    temp(3+i) = 0.0
  enddo
  do j = 1,3
    temp(4) = temp(4) + rv(j,2)*rv(j,3)
    temp(5) = temp(5) + rv(j,1)*rv(j,3)
    temp(6) = temp(6) + rv(j,1)*rv(j,2)
  enddo
  temp(4) = temp(4)/(temp(2)*temp(3))
  temp(5) = temp(5)/(temp(1)*temp(3))
  temp(6) = temp(6)/(temp(1)*temp(2))
  cell_0(4) = radtodeg*acos(temp(4))
  cell_0(5) = radtodeg*acos(temp(5))
  cell_0(6) = radtodeg*acos(temp(6))
  DO i=4,6
     if (abs(cell_0(i) - 90.0 ).lt.0.00001) cell_0(i) = 90.0
     if (abs(cell_0(i) - 120.0).lt.0.00001) cell_0(i) = 120.0
  ENDDO
  RETURN
 END SUBROUTINE uncell
!
 SUBROUTINE inverse(a,c,n)
  implicit none
  integer(i64) ::  n
  real(r64)    ::  a(n,n), c(n,n)
  real(r64)    ::  L(n,n), U(n,n), b(n), d(n), x(n)
  real(r64)    ::  coeff
  integer(i64) ::  i, j, k
  L=0.0
  U=0.0
  b=0.0
  do k=1, n-1
    do i=k+1,n
       coeff=a(i,k)/a(k,k)
       L(i,k) = coeff
       do j=k+1,n
          a(i,j) = a(i,j)-coeff*a(k,j)
       end do
    end do
  end do
  do i=1,n
   L(i,i) = 1.0
  end do
  do j=1,n
   do i=1,j
     U(i,j) = a(i,j)
   end do
  end do
  do k=1,n
   b(k)=1.0
   d(1) = b(1)
   do i=2,n
     d(i)=b(i)
     do j=1,i-1
       d(i) = d(i) - L(i,j)*d(j)
     end do
   end do
   x(n)=d(n)/U(n,n)
   do i = n-1,1,-1
     x(i) = d(i)
     do j=n,i+1,-1
       x(i)=x(i)-U(i,j)*x(j)
     end do
     x(i) = x(i)/u(i,i)
   end do
   do i=1,n
     c(i,k) = x(i)
   end do
   b(k)=0.0
  end do
  RETURN
 END SUBROUTINE inverse
!
 subroutine make_dist_matrix(n,cell_0,rv,vr,x,dist_matrix)
  implicit none
  integer(i64),intent(in)    :: n
  real(r64)   ,intent(in)    :: cell_0(6),rv(3,3),vr(3,3),x(3,n)
  real(r64)   ,intent(out)   :: dist_matrix(n,n)
  integer(i64)               :: i,j,k
  real(r64)                  :: r1(3),r2(3),r3(3),s
  DO i=1,n
     dist_matrix(i,i)=0.0
     DO j=i+1,n
        forall ( k=1:3 )
         r1(k)=x(k,i)
         r2(k)=x(k,j)
        end forall
        call make_distances(.false.,cell_0,r1,r2,rv,r3,s)
        dist_matrix(i,j)=s
        dist_matrix(j,i)=dist_matrix(i,j)
     END DO
  END DO
  return
 end subroutine make_dist_matrix
!
 SUBROUTINE make_distances(flag,cell_0,r2,r1,rv,r3,dist)
  IMPLICIT NONE
  REAL(r64),    intent(in)  :: r1(3),r2(3),rv(3,3),cell_0(6)    ! coordenadas y matriz de cambio
  REAL(r64),    intent(out) :: dist,r3(1:3)
  REAL(r64)                 :: d_image(1:27),image(3,27)        ! array de distancias
  INTEGER(i64)              :: k,l,m,n,o,i,j                    ! variables mudas
  REAL(r64)                 :: atom(3),ouratom(3)               ! coordenadas preparadas
  logical                   :: flag
  real(r64)                 :: phi
  k=0
  do l=-1,1
   do m=-1,1
      do n=-1,1
         k = k + 1
         ouratom(1) = r1(1)
         ouratom(2) = r1(2)
         ouratom(3) = r1(3)
         atom(1) = r2(1) + l
         atom(2) = r2(2) + m
         atom(3) = r2(3) + n
         d_image(k) = distance(atom,ouratom,rv)
         forall ( i=1:3)
           image(i,k) = atom(i)
         end forall
     enddo
   enddo
  enddo
  dist=MINVAL(d_image)
  if(flag)then
   phi=1000.0
   k=1
   do l=1,27
    if(d_image(l)<=phi)then
!       PRINT*,d_image(l),( image(m,l), m=1,3 )
      phi=d_image(l) ! seleccionamos el parametro menor
      k=l            ! y el contador correspondiente.
    endif
   enddo
   forall ( l=1:3)
     r3(l)=image(l,k)
   end forall
  else
   r3(1:3)=0.0
  end if
  return
 END SUBROUTINE make_distances
!
 REAL(r64) FUNCTION DISTANCE(atom,ouratom,rv)
  IMPLICIT NONE
  INTEGER(i64) :: j
  REAL(r64) :: atom(3),ouratom(3),per_atom(3),dist(3),o_atom(3),o_ouratom(3)
  REAL(r64) :: rv(3,3)
  FORALL ( j=1:3 )
   o_ouratom(j) = rv(j,1)*ouratom(1)  + rv(j,2)*ouratom(2)  + rv(j,3)*ouratom(3)
   o_atom(j)    = rv(j,1)*atom(1) + rv(j,2)*atom(2) + rv(j,3)*atom(3)
   dist(j) = o_ouratom(j) - o_atom(j)
  END FORALL
  DISTANCE = sqrt(dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3))
  return
 END FUNCTION
!
 subroutine print_help()
  print '(a)', '  -h, --help   print usage information and exit'
  print '(a)', '  -t, --topol  Topology'
  print '(a)', '  -l, --linker [imi] imi, mimi (real molar fraction [1.0])'
  print '(a)', 'Example:'
  print '(a)', '$ ./generate_ZIFs_linkers -t SOD -l 2 im 25.0 mim 75.0'
 end subroutine print_help
end program Mixing_MOF_generator
