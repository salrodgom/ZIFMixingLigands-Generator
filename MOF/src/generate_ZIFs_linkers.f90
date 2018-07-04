module mod_random
! module for pseudo random numbers
 implicit none
 private
 public init_random_seed, randint, r4_uniform
 contains
 subroutine init_random_seed(seed)
  implicit none
  integer, intent(out) :: seed
! local
  integer   day,hour,i4_huge,milli,minute,month,second,year
  parameter (i4_huge=2147483647)
  double precision temp
  character*(10) time
  character*(8) date
  call date_and_time (date,time)
  read (date,'(i4,i2,i2)')year,month,day
  read (time,'(i2,i2,i2,1x,i3)')hour,minute,second,milli
  temp=0.0D+00
  temp=temp+dble(month-1)/11.0D+00
  temp=temp+dble(day-1)/30.0D+00
  temp=temp+dble(hour)/23.0D+00
  temp=temp+dble(minute)/59.0D+00
  temp=temp+dble(second)/59.0D+00
  temp=temp+dble(milli)/999.0D+00
  temp=temp/6.0D+00
  doext: do
    if(temp<=0.0D+00 )then
       temp=temp+1.0D+00
       cycle doext
    else
       exit doext
    end if
  enddo doext
  doext2: do
    if (1.0D+00<temp) then
       temp=temp-1.0D+00
       cycle doext2
    else
       exit doext2
    end if
  end do doext2
  seed=int(dble(i4_huge)*temp)
  if(seed == 0)       seed = 1
  if(seed == i4_huge) seed = seed-1
  return
 end subroutine init_random_seed
!
 integer function randint(i,j,seed)
  integer,intent(in) :: i,j,seed
  real               :: r
  CALL RANDOM_NUMBER(r)
  randint=int(r*(j+1-i))+i
 end function randint
! 
 real function r4_uniform(x,y,seed)
  implicit none
  real,intent(in)    :: x,y 
  real               :: r
  integer,intent(in) :: seed
  CALL RANDOM_NUMBER(r)
  r4_uniform=(r*(y-x))+x
  return
 end function r4_uniform
end module

program zif_generator
 use iso_fortran_env
 use mod_random
! use topology_agents,only generate
 implicit none
 integer                                       :: i,j,k,l,h,z,m,ierr,nn,seed,iiii
 integer                                       :: linker_type_max=10
 integer                                       :: ii,jj,kk
 real                                          :: rrr,ppp,qqq
 real                                          :: atom(3),ouratom(3)
 integer                                       :: num_args
 real,parameter                                :: r_min_criteria_connectivity=0.56
 real,parameter                                :: r_min_criteria_overlap = 1.0
 integer,parameter                             :: max_number_tries=1000
 integer                                       :: n_atoms = 0,n_nodes=0,n_linkers=0,n_metals=0
 real                                          :: cell_0(1:6) = 0.0, rv(3,3),vr(3,3)
 integer                                       :: n_files=1
 integer                                       :: mc_steps, mc_max_steps=0
  logical                                       :: flag = .true.
 character(len=3)                              :: topology = "Xxx", nodes_code = "Xxx"
 character(len=20)                             :: spam
 character(len=100)                            :: CIFFilename=" "
 character(len=100)                            :: filename=" "
 character(len=80)                             :: string_stop_head= "_atom_site_occupancy"
 character(len=100)                            :: line, string
 character(len=100), dimension(:), allocatable :: args
 character(len=3),dimension(:),allocatable     :: linker_type
 integer                                       :: linker_type_number = 1
 real, allocatable                             :: eee(:)
 real,dimension(:),allocatable                 :: linker_type_molar_fraction
 integer,dimension(:),allocatable              :: genome
 integer,dimension(:),allocatable              :: histogram_molar_fraction
 character(len=3)                              :: code
 real                                          :: molar_fraction
 type                                          :: cluster
  character(len=3)                            :: code
  integer                                     :: n_components
  integer                                     :: id
  real                                        :: component_xcrystal(1:3,1:500)
  real                                        :: component_size(500)
  character(len=2)                            :: component_label(500)
  character(len=2)                            :: element(500)
  real                                        :: comfortably
  logical                                     :: virtual
  real                                        :: population
 end type                                      
 type(cluster),allocatable                     :: nodes(:), linkers(:)
 integer,allocatable                           :: ensemble(:,:)
 call init_random_seed(seed)
 num_args = command_argument_count()
 allocate(args(num_args))
 do i = 1, num_args
  call get_command_argument(i,args(i))
 end do
 do i=1,num_args
  select case(args(i))
   case ('-h','--help')
    call print_help()
    stop
   case ('-l','--linker')
    read(args(i+1),*) linker_type_number                 ! 1
    allocate(linker_type(1:linker_type_number))
    allocate(linker_type_molar_fraction(1:linker_type_number))
    k=0
    do j = 1, linker_type_number
     k=k+2
     read(args(i+k),'(a)') linker_type(j)                ! 2,4
     read(args(i+k+1),*) linker_type_molar_fraction(j) ! 3,5
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
  linker_type_molar_fraction(1)=1.0
 end if
 rrr=sum(linker_type_molar_fraction)
 linker_type_molar_fraction=linker_type_molar_fraction/rrr
 write(6,'(80a)')('=',j=1,80)
 write(6,*)( linker_type(j),j=1,linker_type_number)
 write(6,*)( linker_type_molar_fraction(j)/rrr ,j=1,linker_type_number)
 write(6,'(80a)')('=',j=1,80)
 allocate(histogram_molar_fraction(linker_type_number))
 call system("if [ -f tmp  ] ; then rm tmp  ; touch tmp  ; fi")
 write(6,'(a)') "if [ -f tmp  ] ; then rm tmp  ; touch tmp  ; fi"
 call system("if [ -f list ] ; then rm list ; touch list ; fi")
 write(6,'(a)') "if [ -f list ] ; then rm list ; touch list ; fi"
 do j=1,linker_type_number
  i=len(trim(linker_type(j)))
  k=LEN(TRIM(topology))
  string="ls ???_"//topology(1:k)//"_cif_gin_all/???_???_"//linker_type(j)(1:i)//"_"//topology(1:k)//"_*.cif > tmp"
  call system(string)
  write(6,'(a)')string
  write(line,*) linker_type_molar_fraction(j)
  k=len(trim(line))
  string="awk -v type="      //linker_type(j)(1:i)//&
   " -v r="      //line(4:k)//" '{print $1,type,r}' tmp >> list"
  write(6,'(a)')string
  call system(string)
 end do
 open(111,file="list",iostat=ierr)
 n_files=0
 do
  read(111,'(a)',iostat=ierr) line 
  if(ierr/=0) exit
  read(line(1:41),'(a)') CIFFilename
  CIFFilename=adjustl(trim(CIFFilename))
  write(6,'(a,1x,a1)')CIFFilename,'#'
  read(line(42:),*) code,molar_fraction
  write(6,'(a,1x,a3,1x,f14.7)')trim(CIFFilename),trim(code),molar_fraction
  n_files=n_files+1
  if(n_files==1)then
  !Investigate number of Zn -> number of ligands -> atoms/ligands
   open(121,file=CIFFilename,iostat=ierr,status='old')
   if(ierr/=0)stop 'File does not found'
   do 
    read(121,'(a)',iostat=ierr) line
    if(ierr/=0)exit
    if(line(1:2)=='Zn'.or.line(3:4)=='Zn')then
     write(6,*) "Detecting Zn cluster"
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
 mc_max_steps=1000*n_linkers
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
  read(line(1:41),'(a)') CIFFilename
  read(line(42:) ,*) code,molar_fraction
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
    call check_atom_type(nodes(z)%component_label(k),nodes(z)%component_size(k),nodes(z)%element(k))
    write(6,*)nodes(z)%component_label(k),( nodes(z)%component_xcrystal(m,k),m=1,3)
   end do
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
    call check_atom_type(linkers(h)%component_label(k),linkers(h)%component_size(k),linkers(h)%element(k))
    write(6,*)linkers(h)%component_label(k),( linkers(h)%component_xcrystal(m,k),m=1,3)
   end do
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
   k=ensemble(j,randint(1,n_files,seed))
   flag = .false.
   if ( linkers(k)%code /= linker_type( h ) ) flag = .true.
   if ( linkers(k)%virtual .eqv. .false. )    flag = .true.
   !write( 6,*)  linkers(k)%code, linker_type( h ), linkers(k)%virtual, flag
  end do
  if(j==1)then
   genome(1)=k
   linkers(k)%virtual=.false.
   nn=nn+1
   write(6,*)100*j/real(n_linkers),genome(j),nn, linkers(k)%code, linkers(k)%virtual
   cycle add_linkers
  end if
  !do l=1,j-1 ! scan previous linkers
  ! if(genome(l)==k) then
  !  j=j-1
  !  nn=nn+1
  !  cycle add_linkers ! is it really new?
  ! end if
  !end do
  ! great!
  genome(j)=k
  linkers(k)%virtual=.false.
  write(6,*)100*j/real(n_linkers),genome(j),nn, linkers(k)%code, linkers(k)%virtual
  nn= 0
 end do add_linkers
 write(6,'(a)') "Scanning possible overlaps among atoms of different linkers [...]"
 if ( overlap() ) write(6,'(a)') "[Warnning] Overlap in initial configuration."
 write(6,'(80a)')('=',l=1,80)
 write(6,'(20(a5,1x))')( linker_type(i),i=1,linker_type_number )
 write(6,'(20(f10.8,1x))') ( linker_type_molar_fraction(i), i=1,linker_type_number )
 write(6,'((1000(i6,1x)))') ( genome(i),i=1,n_linkers )
 write(6,'(80a)')('=',l=1,80)
 call update_comfortably()
 mc_steps=0
 write(6,'(1000(i6,1x))')(genome(i),i=1,n_linkers)
 call writeCIFFile_from_clusters()
 !stop
 rrr=0.0
 write(6,'(a,1x,i5)') "Performing MC, max number of steps:", mc_max_steps
 mc_exchange_linkers: do mc_steps=1,mc_max_steps
  allocate(eee(0:1))
  l=randint(1,n_linkers,seed)               ! l:= position of linker in genome
  eee(0)=cost_per_linker(l) + cost_molar()
  j=genome( l )                             ! j:= old linker
  k=ensemble(l,randint(1,n_files,seed))     ! k:= new posible linker
  do while (k==j)
   k=ensemble(l,randint(1,n_files,seed))
  end do
  genome(l) = k                     ! test linker-k
  linkers(j)%virtual=.true.         ! j <- virtual
  linkers(k)%virtual=.false.        ! k <- real
  eee(1) = cost_per_linker(l) + cost_molar()
  if( eee(1) >= eee(0)  )then
   genome(l) = j                    ! rechazo el cambio
   linkers(k)%virtual=.true.        ! k <- virtual
   linkers(j)%virtual=.false.       ! j <- real
  !end if
  else
   ppp=cost()/real_n_atoms()
   if(ppp-rrr<0.0) then
    write(6,'((i5,1x,e20.10,1x,e20.10,1x,e20.10,1x,a,1000(f14.7,1x)))')&
    mc_steps,ppp,ppp-rrr,cost_molar()/real_n_atoms(),&
   'molar fractions:',(histogram_molar_fraction(i)/real(n_linkers),i=1,linker_type_number)
   end if
   rrr=ppp
  end if
  call writeCIFFile_from_clusters()
  deallocate( eee )
 end do mc_exchange_linkers
 write(6,'((i5,1x,e20.10,1x,e20.10,1x,e20.10,1x,a,1000(f14.7,1x)))')&
  mc_max_steps,ppp,ppp-rrr,cost_molar()/real_n_atoms(),&
  'molar fractions:',(histogram_molar_fraction(i)/real(n_linkers),i=1,linker_type_number)
 write(6,'(1000(i6,1x))')(genome(i),i=1,n_linkers)
! finish program
 if ( overlap() ) write(6,'(a)') "[Warnning] Overlap in last configuration."
 stop "Finish program"
 contains
!
 logical function overlap()
  implicit none
  integer :: i,j,ii,jj,kk
  real    :: rrr
  overlap = .false.
  scan_overlap: do i=1,n_linkers
   do j=i+1,n_linkers
    do ii=1,linkers(genome(i))%n_components
     do jj=1,linkers(genome(j))%n_components
      do kk=1,3
       atom(kk)=linkers(genome(i))%component_xcrystal(kk,ii)
       ouratom(kk)=linkers(genome(j))%component_xcrystal(kk,jj)
      end do
      call make_distances(cell_0,ouratom,atom,rv,rrr)
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
 subroutine regrow_move()
  integer :: l,ll,j,jj,k,kk,iiii,jjjj,h,hh
  allocate(eee(0:n_files))
  call choose_by_comfort( l )
  ll = randint(1,n_linkers,seed)
  do while ( ll==l .or. linkers(genome(l))%code == linkers(genome(ll))%code) 
   ll = randint(1,n_linkers,seed)
   write(6,*) linkers(genome(l))%code
   write(6,*) linkers(genome(ll))%code
  end do
  j = genome( l  )
  jj= genome( ll )
  eee(0)=cost_per_linker(l) + cost_per_linker(ii) + cost_molar()
  do iiii=1,n_files
   !
   k=ensemble(l,iiii)
   genome(l) = k
   linkers(j)%virtual=.true.
   linkers(k)%virtual=.false.
   !
   kk=ensemble( ll, iiii )
   genome( ll) = kk
   linkers(jj)%virtual = .true.
   !
   eee(iiii)= cost_per_linker(l) + cost_per_linker(l) + cost_molar()
   !
   write(6,*) iiii
   write(6,*) linkers(j)%code, linkers(k)%code
   write(6,*) linkers(jj)%code, linkers(kk)%code
   write(6,*) eee(iiii)
   !
   if( eee(iiii) < eee(0) )then
    eee(0)=eee(iiii)
    h=k
    hh=kk
   end if
   genome(l) = j
   genome(ll) = jj
   linkers(k)%virtual=.true.
   linkers(j)%virtual=.false.
   linkers(kk)%virtual=.true.
   linkers(jj)%virtual=.false.
  end do
  genome(l) = h
  genome(ll) = hh
  linkers(j)%virtual=.true.
  linkers(h)%virtual=.false.
  linkers(jj)%virtual=.true.
  linkers(hh)%virtual=.false.
  ppp=cost() 
  spam="Regrow move  " 
  deallocate( eee )
  return
 end subroutine regrow_move
!
 subroutine exchange_rotation_move( )
  allocate(eee(0:n_files))
  eee(0)=cost_per_linker(l) + cost_molar()
  j=genome( l )
  h=genome( l )
  do iiii=1,n_files
   ! recorro todo el colectivo de ligandos:
   k=ensemble(l,iiii)
   !
   genome(l) = k
   linkers(j)%virtual=.true.
   linkers(k)%virtual=.false.
   ! mido:
   eee(iiii)= cost_per_linker(l) + cost_molar()
   !write(6,*)'File:', iiii, 'ligand configuration:',k,'by',j,'energy:',eee(iiii),'ligand:',l,'Total Energy:',cost()
   if( eee(iiii) < eee(0) )then
    !write(6,*)'Detecting candidate',k,' [...]'
    eee(0)=eee(iiii)
    h=k
   end if
   ! 
   genome(l) = j
   linkers(k)%virtual=.true.
   linkers(j)%virtual=.false.
  end do
  spam=" "
  if ( linkers(h)%code(1:3) == linkers(genome(l))%code(1:3) ) then
   spam="Rotation move"
  else
   spam="Exchange move"
  end if
  genome(l) = h
  linkers(j)%virtual=.true.
  linkers(h)%virtual=.false.
  !
  deallocate(eee)
  ppp=cost()
  genome(l) = j
  linkers(k)%virtual=.true.
  linkers(j)%virtual=.false.
  return
 end subroutine exchange_rotation_move
!
 subroutine update_molar_fraction()
  implicit none
  integer         :: abc
  integer         :: iabc,ijk
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
 real function cost_molar()
  implicit none
  real, parameter  ::  molar_constant=1e10
  integer          ::  abc
  call update_molar_fraction()
  cost_molar=0.0
  do abc=1,linker_type_number
   cost_molar=cost_molar+&
   molar_constant*(histogram_molar_fraction(abc)/real(n_linkers)-linker_type_molar_fraction(abc))**2
  end do
  return
 end function  cost_molar
!
 integer function real_n_atoms()
  implicit none
  integer              :: i
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
  integer :: i
  real    :: base_comfort = 0.0
  real    :: comfort(n_linkers),total_comfort
  real    :: infinity = huge(0.0)
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
  integer             :: ii,k
  integer,intent(in)  :: pivot
  integer,intent(out) :: choosen
  real                :: sum_weight, eta
  real                :: weight(n_files)
  do ii=1,n_files
   weight(ii)=linkers(ensemble(pivot,ii))%population 
  end do
  sum_weight = sum( weight(1:n_files ) )
  eta = r4_uniform(0.0, sum_weight, seed)
  do choosen = 1, n_files
    if ( eta <  linkers(ensemble(pivot,choosen))%population ) exit
    eta = eta - linkers(ensemble(pivot,choosen))%population
  end do 
  !do ii=1,n_files  ! debugging
  ! write(6,*)ii,linkers(ensemble(pivot,ii))%population,linkers(ensemble(pivot,ii))%code,&
  ! linkers(ensemble(pivot,ii))%n_components,choosen,eta
  !end do
  return 
 end subroutine choose_by_molar_fraction
 subroutine choose_by_comfort(choosen)
  implicit none
  integer,intent(inout) :: choosen
  integer :: ii
  real :: comfort(n_linkers)
  real :: sum_weight, eta
  call update_comfortably()
  sum_weight = sum( linkers(genome(1:n_linkers ))%comfortably  )  
  eta = r4_uniform(0.0, sum_weight, seed)
  do choosen = 1, n_linkers
    if ( eta <  linkers(genome(choosen))%comfortably ) exit
    eta = eta - linkers(genome(choosen))%comfortably
  end do
  !do ii=1,n_linkers  ! debugging
  ! write(6,*)ii,linkers(genome(ii))%comfortably,choosen
  !end do
  !if(choosen>n_linkers) stop 'choosen > n_linkers'
  return
 end subroutine choose_by_comfort
 real function cost_per_linker(identi)
  ! 11 may 2017
  implicit none
  integer,intent(in) :: identi
  integer            :: ii,jj,k ,i
  real               :: r
  real               :: ell=0.5
  real               :: rll=3.9
  real               :: rml=0.0
  real               :: eml=0.0
  real               :: infinite = 1.0e34
  cost_per_linker=0.0
  do i=1,n_linkers
   if(i/=identi)then
    do ii=1,linkers(genome(i))%n_components
     do jj=1,linkers(genome(identi))%n_components
      forall (k=1:3)
       atom(k)=linkers(genome(identi))%component_xcrystal(k,jj)
       ouratom(k)=linkers(genome(i))%component_xcrystal(k,ii)
      end forall
      call make_distances(cell_0,ouratom,atom,rv,r)
      cost_per_linker = cost_per_linker + ell*((rll/r)**12-2*(rll/r)**6)
     end do
    end do
   end if
  end do
  !do i=1,n_nodes
  ! do ii=1,linkers(genome(identi))%n_components
  !  forall (k=1:3)
  !   atom(k)=linkers(genome(identi))%component_xcrystal(k,ii)
  !   ouratom(k)=nodes(i)%component_xcrystal(k,ii)
  !  end forall
  !  call make_distances(cell_0,ouratom,atom,rv,r)
  !  cost_per_linker = cost_per_linker + eml*((rml/r)**12-2*(rml/r)**6)
  ! end do
  !end do
  if( cost_per_linker > infinite ) cost_per_linker = infinite
  if( cost_per_linker /= cost_per_linker ) cost_per_linker= infinite
  return
 end function cost_per_linker
 real function cost_exchange()
  implicit none 
  integer               :: i
  cost_exchange=0.0
  do i=1,n_linkers
   !if(linkers(genome(i))%virtual.eqv..true.) then
   ! write(6,*)( linkers(genome(i))%virtual,ii=1,n_linkers )
   ! stop 'genome with virtual linker'
   !end if
   !do j=i+1,n_linkers
   ! !if(linkers(genome(j))%virtual.eqv..true.)then
   ! ! write(6,*)( linkers(genome(j))%virtual,ii=1,n_linkers )
   ! ! stop 'genome with virtual linker'
   ! !end if
   ! do ii=1,linkers(genome(i))%n_components
   !  do jj=1,linkers(genome(j))%n_components
   !   forall (k=1:3)
   !    atom(k)=linkers(genome(j))%component_xcrystal(k,jj)
   !    ouratom(k)=linkers(genome(i))%component_xcrystal(k,ii)
   !   end forall
   !   call make_distances(cell_0,ouratom,atom,rv,r)
   !   cost_exchange = cost_exchange + 0.04*((2.5/r)**12-(2.5/r)**6)
   !  end do
   ! end do
   !end do
   !write(6,*) cost_per_linker(i)
   cost_exchange= cost_exchange+cost_per_linker(i)
  end do
  return
 end function cost_exchange
 real function cost()
  implicit none
  real :: infinite = 1.0e34
  cost = cost_exchange() + cost_molar()
  if ( cost > infinite ) cost = infinite
  return
 end function  cost
!
 subroutine check_topology_composition(top,mmm,nnn,lll,node_code)
  implicit none
  character(len=3),intent(in)  :: top
  integer,intent(in)           :: mmm
  integer,intent(out)          :: nnn
  integer,intent(out)          :: lll
  character(len=3),intent(out) :: node_code
  select case(top)
   case("LP1","LP2")   ! IRMOF-10 interprenetrate large pore
    ! Zn4O(BDC)3
    node_code="oxo"
    nnn = int( nnn/4.0 )
    lll = 3 * nnn
   case("CP1","CP2")   ! IRMOF-10 interpenetrate  narrow pore
    ! Zn4O(BDC)3
    node_code="oxo"
    nnn = int( nnn/4.0 )
    lll = 3 * nnn
   case("I10")   ! IRMOF-10 non-interpenetrate cubic form
    ! Zn4O(BDC)3
    node_code="oxo"
    nnn = int( nnn/4.0 )
    lll = 3 * nnn
   case default  ! ZIF with IZA code
    node_code="Zn_"
    ! ZnN2( XIm )2
    write(6,'(a)') "detecting Zeolitic Imidazolate Framework" 
    nnn = mmm
    lll = 2 * nnn
  end select 
 end subroutine check_topology_composition
!
 subroutine check_cluster_type(lll,nnn)
  implicit none
  character(len=3),intent(in)  :: lll
  integer         ,intent(out) :: nnn
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
  real,intent(out)            :: s
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
    write(6,'(a)')"============="
    write(6,*) lll
    STOP 'Atom unknowed'
  end select check_atom
  return
 end subroutine check_atom_type
 subroutine writeCIFFile_from_clusters()
  implicit none
  integer,parameter :: u=234
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

 PURE INTEGER FUNCTION Clen(s)      ! returns same result as LEN unless:
 CHARACTER(*),INTENT(IN) :: s       ! last non-blank char is null
 INTEGER :: i
 Clen = LEN(s)
 i = LEN_TRIM(s)
 IF (s(i:i) == CHAR(0)) Clen = i-1  ! len of C string
 END FUNCTION Clen

 PURE INTEGER FUNCTION Clen_trim(s) ! returns same result as LEN_TRIM unless:
 CHARACTER(*),INTENT(IN) :: s       ! last char non-blank is null, if true:
 INTEGER :: i                       ! then len of C string is returned, note:
                                    ! Ctrim is only user of this function
 i = LEN_TRIM(s) ; Clen_trim = i
 IF (s(i:i) == CHAR(0)) Clen_trim = Clen(s)   ! len of C string
 END FUNCTION Clen_trim

 SUBROUTINE cell(rv,vr,cell_0)
 implicit none
 integer :: i,j
 real, intent(in)  :: cell_0(6)
 real, intent(out) :: rv(3,3),vr(3,3)
 real, parameter   :: pi = ACOS(-1.0)
 real :: alp,bet
 real :: cosa,cosb,cosg
 real :: gam,sing
 real :: DEGTORAD
 DEGTORAD=pi/180.0
 IF(cell_0(4) == 90.0) THEN
   cosa = 0.0
 ELSE
   ALP=cell_0(4)*degtorad
   COSA=cos(ALP)
 ENDIF
 IF(cell_0(5) == 90.0) THEN
   cosb = 0.0
 ELSE
   bet = cell_0(5)*degtorad
   cosb = cos(bet)
 ENDIF
 IF(cell_0(6) == 90.0) then
   sing = 1.0
   cosg = 0.0
 ELSE
   gam = cell_0(6)*degtorad
   sing = sin(gam)
   cosg = cos(gam)
 ENDIF
 rv(1,1) = cell_0(1)
 rv(1,2) = cell_0(2)*cosg
 rv(1,3) = cell_0(3)*cosb
 rv(2,1) = 0.0
 rv(2,2) = cell_0(2)*sing
 rv(2,3) = cell_0(3)*(cosa - cosb*cosg)/sing
 rv(3,1) = 0.0
 rv(3,2) = 0.0
 rv(3,3) = sqrt( cell_0(3)*cell_0(3) - rv(1,3)*rv(1,3) - rv(2,3)*rv(2,3))
 call inverse(rv,vr,3)
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
  real,intent(out)   :: cell_0(6)
  real,intent(in)    :: rv(3,3)
  integer            :: i,j
  real               :: temp(6)
  REAL               :: radtodeg
  REAL, PARAMETER    :: pi=ACOS(-1.0) 
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
 integer n
 real a(n,n), c(n,n)
 real L(n,n), U(n,n), b(n), d(n), x(n)
 real coeff
 integer i, j, k
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
subroutine make_dist_matrix(n,cell_0,rv,vr,x,dist_matrix)
 implicit none
 integer,intent(in) :: n
 real,intent(in)    :: cell_0(6),rv(3,3),vr(3,3),x(3,n)
 real,intent(out)   :: dist_matrix(n,n)
 integer            :: i,j,k
 real               :: r1(3),r2(3),s
 DO i=1,n
    dist_matrix(i,i)=0.0
    DO j=i+1,n
       forall ( k=1:3 )
        r1(k)=x(k,i)
        r2(k)=x(k,j)
       end forall
       call make_distances(cell_0,r1,r2,rv,s)
       dist_matrix(i,j)=s
       dist_matrix(j,i)=dist_matrix(i,j)
    END DO
 END DO
 return
end subroutine make_dist_matrix
!
 SUBROUTINE make_distances(cell_0,r2,r1,rv,dist)
 IMPLICIT NONE
 REAL,    intent(in)  :: r1(3),r2(3),rv(3,3),cell_0(6)    ! coordenadas y matriz de cambio
 REAL,    intent(out) :: dist
 REAL                 :: d_image(1:27),image(3,27)        ! array de distancias
 INTEGER              :: k,l,m,n,o,i,j                    ! variables mudas
 REAL                 :: atom(3),ouratom(3)               ! coordenadas preparadas
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
  RETURN
 END SUBROUTINE
!
 REAL FUNCTION DISTANCE(atom,ouratom,rv)
  IMPLICIT NONE
  INTEGER :: j
  REAL :: atom(3),ouratom(3),per_atom(3),dist(3),o_atom(3),o_ouratom(3)
  REAL :: rv(3,3)
  FORALL ( j=1:3 )
   o_ouratom(j) = rv(j,1)*ouratom(1)  + rv(j,2)*ouratom(2)  + rv(j,3)*ouratom(3)
   o_atom(j)    = rv(j,1)*atom(1) + rv(j,2)*atom(2) + rv(j,3)*atom(3)
   dist(j) = o_ouratom(j) - o_atom(j)
  END FORALL
  DISTANCE = sqrt(dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3))
  DISTANCE = sqrt(dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3))
 END FUNCTION
 subroutine print_help()
    print '(a)', '  -h, --help   print usage information and exit'
    print '(a)', '  -t, --topol  Topology'
    print '(a)', '  -l, --linker [imi] imi, mimi (real molar fraction [1.0])'
    print '(a)', 'Example:'
    print '(a)', '$ ./generate_ZIFs_linkers -t SOD -l 2 im 25.0 mim 75.0'
 end subroutine print_help
end program zif_generator
