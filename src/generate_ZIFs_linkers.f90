module mod_random
! module for pseudo random numbers
 implicit none
 private
 public init_random_seed, randint, r4_uniform
contains
 subroutine init_random_seed(seed)
  implicit none
  integer, intent(out) :: seed
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

 integer function randint(i,j,seed)
  real               ::  a,b
  integer,intent(in) ::  i,j,seed
  a = real(i)
  b = real(j)
  randint=int(r4_uniform(a,b+1.0,seed))
 end function randint

 real function r4_uniform(b1,b2,seed)
  implicit none
  real b1,b2
  integer i4_huge,k,seed
  parameter (i4_huge=2147483647)
  if(seed == 0) then
   write(*,'(b1)')'R4_UNIFORM - Fatal error!'
   write(*,'(b1)')'Input value of SEED = 0.'
   stop '[ERROR] Chiquitan chiquitintatan'
  end if
  k=seed/127773
  seed=16807*(seed-k*17773)-k*2836
  if(seed<0) then
    seed=seed+i4_huge
  endif
  r4_uniform=b1+(b2-b1)*real(dble(seed)* 4.656612875D-10)
  return
 end function r4_uniform
end module

program zif_generator
 use iso_fortran_env
 use mod_random
! use topology_agents,only generate
 implicit none
 integer             :: i,j,k,l,h,m,ierr,nn,seed
 integer             :: linker_type_max=4
 integer             :: ii,jj,kk
 real                :: rrr,ppp,qqq
 real                :: atom(3),ouratom(3)
 integer             :: num_args
 real,parameter      :: r_min_criteria_connectivity=0.56
 integer,parameter   :: max_number_tries=10000
 integer             :: n_atoms = 0,n_nodes=0,n_linkers
 real                :: cell_0(1:6) = 0.0, rv(3,3),vr(3,3)
 integer             :: linker_type_number = 1, n_files=1
 character(len=3)    :: topology = "SOD"
 character(len=20)   :: spam
 character(len=100)  :: CIFFilename=" "
 character(len=100)  :: filename=" "
 character(len=80)   :: string_stop_head= "_atom_site_occupancy"
 character(len=100)  :: line,string
 character(len=100), dimension(:), allocatable :: args
 character(len=3),dimension(:),allocatable     :: linker_type
 real,dimension(:),allocatable                 :: linker_type_molar_fraction
 integer,dimension(:),allocatable              :: genome
 integer,dimension(:),allocatable              :: histogram_molar_fraction
 character(len=3) :: code
 real             :: molar_fraction
 type                          :: cluster
  character(len=3)             :: code
  integer                      :: n_components
  integer                      :: id
  real                         :: component_xcrystal(1:3,1:100)
  character(len=2)             :: component_label(100)
  real                         :: fitness
  logical                      :: virtual
 end type
 type(cluster),allocatable     :: linkers(:) 
 type(cluster),allocatable     :: nodes(:)
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
     read(args(i+k),'(4a)') linker_type(j)                ! 2,4
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
  linker_type(1)='im '
  linker_type_molar_fraction(1)=1.0
 end if
 rrr=sum(linker_type_molar_fraction)
 write(6,'(80a)')('=',j=1,80)
 write(6,*)( linker_type(j),j=1,linker_type_number)
 write(6,*)( linker_type_molar_fraction(j)/rrr ,j=1,linker_type_number)
 write(6,'(80a)')('=',j=1,80)
 allocate(histogram_molar_fraction(linker_type_number)
 call system("if [ -f tmp  ] ; then rm tmp  ; touch tmp  ; fi")
 write(6,'(a)') "if [ -f tmp  ] ; then rm tmp  ; touch tmp  ; fi"
 call system("if [ -f list ] ; then rm list ; touch list ; fi")
 write(6,'(a)') "if [ -f list ] ; then rm list ; touch list ; fi"
 do j=1,linker_type_number
  i=len(trim(linker_type(j)))
  k=LEN(TRIM(topology))
  string="ls zif_"//topology(1:k)//"_cif_gin_all/zif_tbp_"//linker_type(j)(1:i)//"_"//topology(1:k)//"_*.cif > tmp"
  call system(string)
  write(6,'(a)')string
  write(line,*) linker_type_molar_fraction(j)
  k=len(trim(line))
  string="awk -v type="//linker_type(j)(1:i)//&
   " -v r="//line(4:k)//" '{print $1,type,r}' tmp >> list"
  write(6,'(a)')string
  call system(string)
 end do
 open(111,file="list",iostat=ierr)
 n_files=0
 do
  read(111,'(a)',iostat=ierr) line 
  if(ierr/=0) exit
  read(line(1:41),'(a)') CIFFilename
  read(line(42:) ,*) code,molar_fraction
  write(6,'(a,1x,a,1x,f10.5)')trim(CIFFilename),trim(code),molar_fraction
  n_files=n_files+1
  if(n_files==1)then
  !Investigate number of Zn -> number of ligands -> atoms/ligands
   open(121,file=CIFFilename,iostat=ierr,status='old')
   if(ierr/=0)stop 'File does not found'
   do 
    read(121,'(a)',iostat=ierr) line
    if(ierr/=0)exit
    if(line(1:4)=='  Zn') n_nodes=n_nodes+1
   end do
   close(121)
   n_linkers=n_nodes*2
  end if
 end do
 rewind(111)
 allocate(nodes(n_nodes))
 allocate(linkers(n_files*n_linkers))
 write(6,'(80a)')('=',j=1,80)
 write(6,*)n_nodes, 'nodes/uc'
 write(6,*)n_linkers,'linkers/uc'
 write(6,*)n_files, 'files detected'
 write(6,*)n_files*n_linkers,' total linkers spected'
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
   read(100,'(a)') line
   if(i==1) then
    nodes(j)%id=j
    nodes(j)%code='Zn'
    nodes(j)%n_components=1
    nodes(j)%virtual=.false.
    read(line,*) nodes(j)%component_label(1),(nodes(j)%component_xcrystal(m,1),m=1,3),rrr
   end if
  end do
  write(6,*)( nodes(j)%component_label(1),j=1,n_nodes)
  do j=1,n_linkers
   h=h+1 ! count total number of linkers
   write(6,'(80a)')('=',l=1,80)
   linkers(h)%id=h
   linkers(h)%code=code
   linkers(h)%n_components=int((n_atoms-n_nodes)/n_linkers)
   linkers(h)%virtual=.true.
   do k=1,int((n_atoms-n_nodes)/n_linkers)
    read(100,*) linkers(h)%component_label(k),&
     ( linkers(h)%component_xcrystal(m,k),m=1,3),rrr
   end do
   write(6,*)( linkers(h)%component_label(k),k=1,int((n_atoms-n_nodes)/n_linkers) )
  end do
  close(100)
 end do
 close(111)
 call system('rm list tmp')
 j=0
 allocate(genome(n_linkers))
 genome=0
 nn=0
 add_linkers: do 
  if(nn==max_number_tries) then
   j=0       ! if the number of tries is bigger than a limit, 
   genome=0  ! we start again with the first ligand.
   nn=0
   linkers%virtual=.true. ! all the linkers are virtuals again!
   write(6,'(a)')'[warning] Relocation this ligands is hard, check the input'
   cycle add_linkers
  end if
  if(j==n_linkers) exit add_linkers
  j=j+1
  k=randint(1,h,seed)  ! new linker
  if(j==1)then
   genome(1)=k
   linkers(k)%virtual=.false.
   nn=nn+1
   cycle add_linkers
  end if
! overlap?
  do l=1,j-1 ! scan previous linkers
   if(genome(l)==k) then
    j=j-1
    nn=nn+1
    cycle add_linkers ! is it really new?
   end if
   do ii=1,linkers(genome(l))%n_components    
    do jj=1,linkers(k)%n_components
     do kk=1,3
      atom(kk)=linkers(genome(l))%component_xcrystal(kk,ii)
      ouratom(kk)=linkers(k)%component_xcrystal(kk,jj)
     end do
     call make_distances(cell_0,ouratom,atom,rv,rrr)
     if (rrr<=1.3) then
      j=j-1
      nn=nn+1
      cycle add_linkers         ! overlap !!!
     end if
    end do
   end do
  end do
  ! great!
  nn=0
  genome(j)=k
  linkers(k)%virtual=.false.
 end do add_linkers
 write(6,'(80a)')('=',l=1,80)
 write(6,*)( linker_type(i),i=1,linker_type_number )
 write(6,'((24(i3,1x)))') ( genome(i),i=1,n_linkers )
 write(6,'(80a)')('=',l=1,80)
 mc_exchange_linkers: do i=1,2
  l=randint(1,n_linkers,seed)
  j=genome(l)
  k=randint(1,h,seed)
  do while (k==j)
   k=randint(1,h,seed)
  end do
  
 end do mc_exchange_linkers
 call writeCIFFile_from_clusters()
 stop
 contains
 subroutine writeCIFFile_from_clusters()
  implicit none
  integer,parameter :: u=234
  character(len=10) :: outfilename ='output.cif'
  !write(line,'(a1,i1,a)')'(',linker_type_number,'(a3))'
  !write(6,'(a)') line
  !write(string,line)(adjustl(trim(linker_type(ii))), ii=1,linker_type_number )
  !write(filename,'(a,a,a,a)')adjustl(trim(topology)),'_',adjustl(trim(string)),'_mixture.cif' 
  !write(6,'(a)') string
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
  write(u,'(a)')'_atom_site_label'
  write(u,'(a)')'_atom_site_fract_x'
  write(u,'(a)')'_atom_site_fract_y'
  write(u,'(a)')'_atom_site_fract_z'
  do i=1,n_nodes
   if( nodes(i)%virtual.eqv..false.)then
    do j=1,nodes(i)%n_components
     write(u,*) nodes(i)%component_label(j),(nodes(i)%component_xcrystal(m,j),m=1,3)
    end do
   end if
  end do
  do i=1,n_linkers
   if( linkers(genome(i))%virtual.eqv..false.)then
    do j=1,linkers(genome(i))%n_components
     write(u,*) linkers(genome(i))%component_label(j),&
      (linkers(genome(i))%component_xcrystal(m,j),m=1,3)
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
 END FUNCTION
 subroutine print_help()
    print '(a)', '  -h, --help   print usage information and exit'
    print '(a)', '  -t, --topol  Topology'
    print '(a)', '  -l, --linker [imi] imi, mimi (real molar fraction [1.0])'
    print '(a)', 'Example:'
    print '(a)', '$ ./generate_ZIFs_linkers -t SOD -l 2 im 25.0 mim 75.0'
 end subroutine print_help
end program zif_generator
