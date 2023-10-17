program main
use mpi
implicit none
!-----------------------------------------------------------------------
!integer,parameter,dimension(3)::	topo		=[3,3,3]
																	
integer,parameter,dimension(3)::	snapList=[10000,30000,50000]
																	
logical,parameter::								isDebug	=.false.

integer,parameter:: q=18				,&
										nBlocks=3		,&
										time=100000	,&
										freqDisp=1
										
integer,parameter,dimension(2)::	m=[4,2]					!,&
																	!span=[10,14]	

										
double precision,parameter:: 	D=6.0d0				,&
															tauC=0.8d0		,&
															rho0=1.0d0		,&
															uInlet=0.01d0	,&
															
															haf=0.5d0			,&
															one=1.0d0			,&
															two=2.0d0			,&
															zer=0.0d0			,&
															quar=0.25d0		,&
															one36th=1.0d0/36.0d0	,&
															pi=4.0d0*datan(one)
!------------------------MPI Variables----------------------------------
integer,dimension(3)		:: bs,es,bn,en,siz,pP,des,src,myDim,myCoord,period
integer		:: linex, liney, linez, yz1p, xz1p, xy1p, yz1pP, yz1pM, xz1pP, xz1pM, xy1pP, xy1pM
integer,dimension(5)		:: blockLen, xpDisp, xmDisp, ypDisp, ymDisp, zpDisp, zmDisp
integer		:: comm3d,nproc, id, ierr, sizDP, STATUS(mpi_status_size),yz(4),xz(4),xy(4)
logical 	:: myperiod(3), bottomId, topId, eastId, westId, northId, southId, masterId, sphereId,isRed,isGreen,isBlue
integer,allocatable,dimension(:,:) 	::tArr,bb,ee
integer::color,id1,nproc1,commSplit,topo(3),myZone,nProcZone(3)										
!-----------------------------------------------------------------------
integer:: i,j,k,ii,jj,kk,a,a1,ts,ia,ja,ka,iter1,iter2,s,t,solNum,mRel(2)
integer:: ex(0:q),ey(0:q),ez(0:q),kb(0:q),tmpI1, tmpI2
double precision:: tmp1, tmp2, tmp3, tmp4, fx_t, fy_t, fz_t, fx(2), fy(2), fz(2), tStart, tEnd
double precision:: Cd, Cl, Cs, wt(0:q), rhoAvg, rhoAvgAll, xiL(3), yiL(3),xPosObj
!coarse variables
double precision		::ux, uy, uz, rho

type zone
double precision,allocatable,dimension(:,:,:,:)	::f, ft, feq
double precision,allocatable,dimension(:,:,:,:,:)	::ftA
double precision							::tau
integer,dimension(3)					::np,leftFace,rightFace,centr,span
end type zone

type(zone),dimension(nBlocks)::z
integer,allocatable,dimension(:,:,:)						::isn
double precision,allocatable,dimension(:)					::xIN, yIN, zIN, xOUT, yOUT, zOUT
double precision,allocatable,dimension(:,:)				::fINyz,fINxz,fINxy,fOUTyz,fOUTxz,fOUTxy

character(len=50)	:: filename
!=======================================================================
CALL mpi_init(ierr)
CALL mpi_comm_rank(mpi_comm_world,id,ierr)
CALL mpi_comm_size(mpi_comm_world,nproc,ierr)
masterId=(id==0)

z(3)%np=[299,179,179]

z(1)%span=[10,10,10]! -2 to exclude buffer,-1 becoz span is 1 less than ny
z(2)%span=[14,14,14]!span should be even

z(2)%np=m(2)*z(2)%span+3 !includes 2 buffer nodes
z(1)%np=m(1)*z(1)%span+3

xPosObj=75
z(3)%centr(1)=xPosObj
z(3)%centr(2:3)=z(3)%np(2:3)/2+1
z(2)%centr=z(2)%np/2+1
z(1)%centr=z(1)%np/2+1

mRel=[m(1)/m(2),m(2)]

z(3)%tau=tauC
z(1:2)%tau = haf + m(1:2)*(z(3)%tau - haf)

z(1:2)%leftFace (1)	= z(3)%centr(1)-z(1:2)%span(1)/2
z(1:2)%rightFace(1)	= z(3)%centr(1)+z(1:2)%span(1)/2
z(1:2)%leftFace (2)	= z(3)%centr(2)-z(1:2)%span(2)/2
z(1:2)%rightFace(2)	= z(3)%centr(2)+z(1:2)%span(2)/2
z(1:2)%leftFace (3)	= z(3)%centr(3)-z(1:2)%span(3)/2
z(1:2)%rightFace(3)	= z(3)%centr(3)+z(1:2)%span(3)/2

if(isDebug .and. masterId) then
	write(*,*) "Re: ",3.0d0*D*uInlet/(z(3)%tau-haf)
	write(*,*) "----------------------"
	write(*,*) "grid ratio    :",m
	write(*,*) "rel grid ratio:",mRel
	write(*,*) "----------------------"
	write(*,*) "zone3.np: ",id,z(3)%np
	write(*,*) "zone2.np: ",id,z(2)%np
	write(*,*) "zone1.np: ",id,z(1)%np
	write(*,*) "----------------------"
	write(*,*) "zone3.span: ",id,z(3)%span
	write(*,*) "zone2.span: ",id,z(2)%span
	write(*,*) "zone1.span: ",id,z(1)%span
	write(*,*) "----------------------"
	write(*,*) "zone3.centr: ",id,z(3)%centr
	write(*,*) "zone2.centr: ",id,z(2)%centr
	write(*,*) "zone1.centr: ",id,z(1)%centr
	write(*,*) "----------------------"
	write(*,*) "zone3.leftface : ",id,z(3)%leftFace
	write(*,*) "zone3.rightface: ",id,z(3)%rightFace
	write(*,*) "zone2.leftface : ",id,z(2)%leftFace
	write(*,*) "zone2.rightface: ",id,z(2)%rightFace
	write(*,*) "zone1.leftface : ",id,z(1)%leftFace
	write(*,*) "zone1.rightface: ",id,z(1)%rightFace
	write(*,*) "----------------------"
	write(*,*) "zone3.tau: ",id,z(3)%tau
	write(*,*) "zone2.tau: ",id,z(2)%tau
	write(*,*) "zone1.tau: ",id,z(1)%tau
	write(*,*) "----------------------"
endif
!==============================MPI======================================
nProcZone=[27,27,27]
!color=id/nProcZone(1)+max(0,(id-nProcZone(1))/nProcZone(2))+max(0,(id-nProcZone(1)-nProcZone(2))/nProcZone(3))

do i=0,sum(nProcZone)-1
	if(i==id .and. i .lt. nProcZone(1)) color=3
	if(i==id .and. i .ge. nProcZone(1) .and. i .lt. (nProcZone(1)+nProcZone(2)) ) color=2
	if(i==id .and. i .ge. (nProcZone(1)+nProcZone(2)) .and. i .lt. sum(nProcZone) ) color=1
enddo

!color=3-color
!~ do i=0,nProc
!~ 	if(i==id)	write(*,*) id,color
!~ 	call sleep(1)
!~ enddo
!~ stop
!~ color=id/27
CALL mpi_comm_split(mpi_comm_world,color,id,commSplit,ierr)
CALL mpi_comm_rank(commSplit,id1,ierr)
CALL mpi_comm_size(commSplit,nproc1,ierr)
!write(*,*) id,id1,color,nproc,nproc1
!call sleep(2)

if(isRedID) then
	myZone=3
	topo=[3,3,3]
	period	=[.false. , .false. , .false.]
endif
	
if(isGreenID) then
	myZone=2
	topo=[3,3,3]
	period	=[.false. , .false. , .false.]
endif
	
if(isBlueID) then
	myZone=1
	topo=[3,3,3]
	period	=[.false. , .false. , .false.]
endif


call mpi_cart_create(commSplit,3,topo,period,.false.,comm3d,ierr)
call mpi_cart_get(comm3d,3,myDim,myperiod,mycoord,ierr)

call mpi_cart_shift(comm3d,0,1,src(1),des(1),ierr)
call mpi_cart_shift(comm3d,1,1,src(2),des(2),ierr)
call mpi_cart_shift(comm3d,2,1,src(3),des(3),ierr)

allocate(tArr(0:maxval(myDim)-1,3))
allocate(bb(3,0:nproc1-1),ee(3,0:nproc1-1))

eastId	=(mycoord(1)==myDim(1)-1)
westId	=(mycoord(1)==0)
northId	=(mycoord(2)==myDim(2)-1)
southId	=(mycoord(2)==0)
topId		=(mycoord(3)==myDim(3)-1)
bottomId=(mycoord(3)==0)

yz=-1
xz=-1
xy=-1

if(des(2) /= -1 .and. des(3) /= -1) yz(1)	=	des(2)+1
if(src(2) /= -1 .and. des(3) /= -1) yz(2)	=	src(2)+1
if(src(2) /= -1 .and. src(3) /= -1) yz(3)	=	src(2)-1
if(des(2) /= -1 .and. src(3) /= -1) yz(4)	=	des(2)-1

if(des(1) /= -1 .and. des(3) /= -1) xz(1)	=	des(1)+1
if(src(1) /= -1 .and. des(3) /= -1) xz(2)	=	src(1)+1
if(src(1) /= -1 .and. src(3) /= -1) xz(3)	=	src(1)-1
if(des(1) /= -1 .and. src(3) /= -1) xz(4)	=	des(1)-1

if(des(1) /= -1 .and. des(2) /= -1) xy(1)	=	des(2)+topo(2)*topo(3)
if(src(1) /= -1 .and. des(2) /= -1) xy(2)	=	des(2)-topo(2)*topo(3)
if(src(1) /= -1 .and. src(2) /= -1) xy(3)	=	src(2)-topo(2)*topo(3)
if(des(1) /= -1 .and. src(2) /= -1) xy(4)	=	src(2)+topo(2)*topo(3)

do k=1,3

	tmpI1=z(myZone)%np(k)/myDim(k)
	tmpI2=myDim(k)-mod(z(myZone)%np(k),myDim(k))

	tArr(0,1)=1
	do i=0,myDim(k)-1
		if(i==tmpI2) tmpI1=tmpI1+1
		tArr(i,2)=tArr(i,1)+tmpI1-1
		tArr(i,3)=tArr(i,2)-tArr(i,1)+1
		if(i==myDim(k)-1) exit
		tArr(i+1,1)=tArr(i,2)+1
	enddo

	do i=0,myDim(k)-1
		if(i==mycoord(k)) then
			bs(k)=tArr(i,1)
			es(k)=tArr(i,2)
			siz(k)=tArr(i,3)
		endif
	enddo

	bn(k)=bs(k)
	en(k)=es(k)
	
	if((westId .and. k==1) .or. (southId 	.and. k==2) .or. (bottomId .and. k==3)) then
		bn(k)=bs(k)+1
		siz(k)=siz(k)-1
	endif
		
	if((eastId .and. k==1) .or. (northId 	.and. k==2) .or. (topId 	 .and. k==3))	then
		en(k)=es(k)-1
		siz(k)=siz(k)-1
	endif
		
enddo
!-----------------------------------------------------------------------
allocate(z(myZone)%f  (bn(1)-1:en(1)+1, bn(2)-1:en(2)+1, bn(3)-1:en(3)+1, 0:q))
allocate(z(myZone)%ft (bn(1)-1:en(1)+1, bn(2)-1:en(2)+1, bn(3)-1:en(3)+1, 0:q))
allocate(z(myZone)%feq(bn(1)-1:en(1)+1, bn(2)-1:en(2)+1, bn(3)-1:en(3)+1, 0:q))

pP=[(en(1)+1)-(bn(1)-1)+1, (en(2)+1)-(bn(2)-1)+1, (en(3)+1)-(bn(3)-1)+1]
!-----------------------------------------------------------------------
CALL mpi_type_extent (mpi_double_precision,sizDP,ierr)

CALL mpi_type_vector (siz(1), 1, 1					, mpi_double_precision,linex ,ierr)
CALL mpi_type_vector (siz(2), 1, pP(1)			, mpi_double_precision,liney ,ierr)
CALL mpi_type_vector (siz(3), 1, pP(1)*pP(2), mpi_double_precision,linez ,ierr)

CALL mpi_type_commit(linex,ierr)
CALL mpi_type_commit(liney,ierr)
CALL mpi_type_commit(linez,ierr)

CALL mpi_type_hvector(siz(3), 1			, pP(1)*pP(2)*sizDP	, liney								,yz1p,ierr)
CALL mpi_type_vector (siz(3), siz(1), pP(1)*pP(2)				, mpi_double_precision,xz1p,ierr)
CALL mpi_type_vector (siz(2), siz(1), pP(1)							, mpi_double_precision,xy1p,ierr)

CALL mpi_type_commit(yz1p,ierr)
CALL mpi_type_commit(xz1p,ierr)
CALL mpi_type_commit(xy1p,ierr)

blockLen=[1,1,1,1,1]
xpDisp=[1,7 ,9 ,11,13]*product(pP)*sizDP
xmDisp=[2,8 ,10,12,14]*product(pP)*sizDP
ypDisp=[3,7 ,8 ,15,17]*product(pP)*sizDP
ymDisp=[4,9 ,10,16,18]*product(pP)*sizDP
zpDisp=[5,11,12,15,16]*product(pP)*sizDP
zmDisp=[6,13,14,17,18]*product(pP)*sizDP

call mpi_type_hindexed(5,blockLen,xpDisp,yz1p,yz1pP,ierr)
call mpi_type_hindexed(5,blockLen,xmDisp,yz1p,yz1pM,ierr)
call mpi_type_hindexed(5,blockLen,ypDisp,xz1p,xz1pP,ierr)
call mpi_type_hindexed(5,blockLen,ymDisp,xz1p,xz1pM,ierr)
call mpi_type_hindexed(5,blockLen,zpDisp,xy1p,xy1pP,ierr)
call mpi_type_hindexed(5,blockLen,zmDisp,xy1p,xy1pM,ierr)

call mpi_type_commit(yz1pP,ierr)
call mpi_type_commit(yz1pM,ierr)
call mpi_type_commit(xz1pP,ierr)
call mpi_type_commit(xz1pM,ierr)
call mpi_type_commit(xy1pP,ierr)
call mpi_type_commit(xy1pM,ierr)

CALL mpi_barrier(mpi_comm_world,ierr)

if(isDebug) then

	do j=1,3
		if(j==color) then
			do i=0,nproc-1
				if(i==id) then
					write(*,'(A,3(I3),3(A,I5,I5))') "bs-es : ",id,color,id1," | ",bs(1),es(1)," | ",bs(2),es(2)," | ",bs(3),es(3)
					write(*,'(A,3(I3),3(A,I5,I5))') "bn-en : ",id,color,id1," | ",bn(1),en(1)," | ",bn(2),en(2)," | ",bn(3),en(3)
					write(*,'(A,3(I3),3(A,I10  ))') "mySpan: ",id,color,id1," | ",      pP(1)," | ",      pP(2)," | ",      pP(3)
				endif
				call sleep(1)
			enddo
		endif
	enddo
	
	call sleep(1)
	if(masterId) 	write(*,*) "----------------------"
	call sleep(1)
	
	do j=1,3
		if(j==color) then
			if(masterId) 	write(*,*) "Master Proc(s): ",id,color,id1
			if(eastId)		write(*,*) "East   Proc(s): ",id,color,id1
			if(westId) 		write(*,*) "West   Proc(s): ",id,color,id1
			if(northId) 	write(*,*) "North  Proc(s): ",id,color,id1
			if(southId) 	write(*,*) "South  Proc(s): ",id,color,id1
			if(topId) 		write(*,*) "Top    Proc(s): ",id,color,id1
			if(bottomId) 	write(*,*) "Bottom Proc(s): ",id,color,id1
		endif
		call sleep(5)
	enddo
	
	call sleep(1)
	if(masterId) 	write(*,*) "----------------------"

endif
!---------This works if sphere lies within one processor----------------
sphereId = .false.

outermost: do k=bn(3),en(3)
	do j=bn(2),en(2)
		do i=bn(1),en(1)
		
			tmp1 = ((i-z(myZone)%centr(1))**two + (j-z(myZone)%centr(2))**two + (k-z(myZone)%centr(3))**two)**haf 

			if (tmp1 .le. haf*D) then
				sphereId = .true.
				exit outermost
			endif
			
		enddo
	enddo
enddo outermost

if(sphereId) write(*,'(A,3(I4),A)') "Processor ",id,color,id1," has the object"!isDebug .and.
stop
!-----------------------------------------------------------------------
if(sphereId .and. isBlueID) then
!~ 	allocate(z(1)%f  (z(1)%np(1),z(1)%np(2),z(1)%np(3),0:q))
!~ 	allocate(z(1)%ft (z(1)%np(1),z(1)%np(2),z(1)%np(3),0:q))
!~ 	allocate(z(1)%feq(z(1)%np(1),z(1)%np(2),z(1)%np(3),0:q))
!~ 	allocate(z(1)%ftA(z(1)%np(1),z(1)%np(2),z(1)%np(3),0:q,3))
	
!~ 	allocate(z(2)%f  (z(2)%np(1),z(2)%np(2),z(2)%np(3),0:q))
!~ 	allocate(z(2)%ft (z(2)%np(1),z(2)%np(2),z(2)%np(3),0:q))
!~ 	allocate(z(2)%feq(z(2)%np(1),z(2)%np(2),z(2)%np(3),0:q))
!~ 	allocate(z(2)%ftA(z(2)%np(1),z(2)%np(2),z(2)%np(3),0:q,3))	
	
	allocate(isn(z(myZone)%np(1),z(myZone)%np(2),z(myZone)%np(3)))
endif
!----------------------------------------------------------------------
ex(0)=0;       ey(0)=0;       ez(0)=0;
ex(1)=1;       ey(1)=0;       ez(1)=0;
ex(2)=-1;      ey(2)=0;       ez(2)=0;
ex(3)=0;       ey(3)=1;       ez(3)=0;
ex(4)=0;       ey(4)=-1;      ez(4)=0;
ex(5)=0;       ey(5)=0;       ez(5)=1;
ex(6)=0;       ey(6)=0;       ez(6)=-1;
ex(7)=1;       ey(7)=1;       ez(7)=0;
ex(8)=-1;      ey(8)=1;       ez(8)=0;
ex(9)=1;       ey(9)=-1;      ez(9)=0;
ex(10)=-1;     ey(10)=-1;     ez(10)=0;
ex(11)=1;      ey(11)=0;      ez(11)=1;
ex(12)=-1;     ey(12)=0;      ez(12)=1;
ex(13)=1;      ey(13)=0;      ez(13)=-1;
ex(14)=-1;     ey(14)=0;      ez(14)=-1;
ex(15)=0;      ey(15)=1;      ez(15)=1;
ex(16)=0;      ey(16)=-1;     ez(16)=1;
ex(17)=0;      ey(17)=1;      ez(17)=-1;
ex(18)=0;      ey(18)=-1;     ez(18)=-1;
!----------------------------------------------------------------------	
do a=0,q
	if(a .eq. 0) then
	wt(a) = 12*one36th
	elseif(a .ge. 1 .and. a .le. 6)	then
	wt(a) = 2*one36th
	else
	wt(a) = one36th
	endif
enddo
!----------------------------------------------------------------------
do a=0,q 
	do a1=a,q
		if( ex(a)+ex(a1) .eq. 0 .and. ey(a)+ey(a1) .eq. 0 .and. ez(a)+ez(a1) .eq. 0) then
			kb(a)	= a1
			kb(a1)=  a
		endif
	enddo
enddo
!----------------------------------------------------------------------
if(isRedID) then
	do k=bs(3),es(3)
		do j=bs(2),es(2) 
			do i=bs(1),es(1)
				do a=0,q
					tmp1				= uInlet*ex(a)
					tmp2				= uInlet*uInlet
					z(myZone)%f (i,j,k,a)	= wt(a)*rho0*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)!commenting this line causes non-zero Cl,Cs
					z(myZone)%ft(i,j,k,a)	= wt(a)*rho0*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
				enddo
			enddo
		enddo
	enddo
endif
!-----------------------------------
if(isGreenID .or. isBlueID) then
	do k=1,z(myZone)%np(3)
		do j=1,z(myZone)%np(2)
			do i=1,z(myZone)%np(1)
				
				do a=0,q
					tmp1				= uInlet*ex(a)  
					tmp2				= uInlet*uInlet  
					!z(1)%f (i,j,k,a  )	= wt(a)*rho0*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
					z(myZone)%ftA(i,j,k,a,2) = wt(a)*rho0*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
					z(myZone)%ftA(i,j,k,a,1) = wt(a)*rho0*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
				enddo
				
			enddo
		enddo
	enddo
endif
!~ !-----------------------------------
!~ if(sphereId) then
!~ 	do k=1,z(2)%np(3)
!~ 		do j=1,z(2)%np(2)
!~ 			do i=1,z(2)%np(1)
				
!~ 				do a=0,q
!~ 					tmp1				= uInlet*ex(a)  
!~ 					tmp2				= uInlet*uInlet  
!~ 					!z(2)%f (i,j,k,a  )	= wt(a)*rho0*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
!~ 					z(2)%ftA(i,j,k,a,2) = wt(a)*rho0*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
!~ 					z(2)%ftA(i,j,k,a,1) = wt(a)*rho0*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
!~ 				enddo
				
!~ 			enddo
!~ 		enddo
!~ 	enddo
!~ endif
!----------------------------------------------------------------------
if(sphereId .and. isBlueID) then
!---------------------------Detect Sphere-------------------------------
	do k=2,z(myZone)%np(3)-1
		do j=2,z(myZone)%np(2)-1
			do i=2,z(myZone)%np(1)-1
			
				tmp1 = ((i-z(myZone)%centr(1))**two + (j-z(myZone)%centr(2))**two + (k-z(myZone)%centr(3))**two)**haf 

				if (tmp1 .le. m(1)*haf*D) then
					isn(i,j,k) = 1
				else
					isn(i,j,k) = 0
				endif
				
			enddo
		enddo
	enddo
	
	fx(1) = zer
	fy(1) = zer
	fz(1) = zer
	
	open(unit=10,file="tRhoCdCsCl.dat")
	write(10,'(A)') "Variables=time,rho,Cd,Cs,Cl"
	
	write(*,'(A)') "Time loop execution starts"
endif

solNum = 0
!----------------------------------------------------------------------
do ts=1,time

	if(sphereId .and. isBlueID .and. mod(ts-1,freqDisp)==0) call cpu_time(tStart)
	
	if(isRedID)
		do k=bn(3),en(3) !2,nz-1 !Streaming
			do j=bn(2),en(2) !2,ny-1
				do i=bn(1),en(1) !2,nx-1
				
					if(i>(z(myZone-1)%leftFace(1)+1) .and. i<(z(myZone-1)%rightFace(1)-1) .and. &
						 j>(z(myZone-1)%leftFace(2)+1) .and. j<(z(myZone-1)%rightFace(2)-1) .and. &
						 k>(z(myZone-1)%leftFace(3)+1) .and. k<(z(myZone-1)%rightFace(3)-1)) cycle
						 
					do a=0,q
						ia = i+ex(a) 
						ja = j+ey(a)
						ka = k+ez(a)
						
						z(myZone)%f(ia,ja,ka,a) = z(myZone)%ft(i,j,k,a)
					enddo
				enddo
			enddo
		enddo
		!---------------------------------------------------------------------
		call mpi_sendrecv(z(myZone)%f(bn(1),bn(2),en(3)+1,0),1,xy1pP,des(3),50,z(myZone)%f(bn(1),bn(2),bn(3),0),1,xy1pP,src(3),50,commSplit,STATUS,ierr)
		call mpi_sendrecv(z(myZone)%f(bn(1),bn(2),bn(3)-1,0),1,xy1pM,src(3),50,z(myZone)%f(bn(1),bn(2),en(3),0),1,xy1pM,des(3),50,commSplit,STATUS,ierr)

		call mpi_sendrecv(z(myZone)%f(bn(1),en(2)+1,bn(3),0),1,xz1pP,des(2),50,z(myZone)%f(bn(1),bn(2),bn(3),0),1,xz1pP,src(2),50,commSplit,STATUS,ierr)
		call mpi_sendrecv(z(myZone)%f(bn(1),bn(2)-1,bn(3),0),1,xz1pM,src(2),50,z(myZone)%f(bn(1),en(2),bn(3),0),1,xz1pM,des(2),50,commSplit,STATUS,ierr)

		call mpi_sendrecv(z(myZone)%f(en(1)+1,bn(2),bn(3),0),1,yz1pP,des(1),50,z(myZone)%f(bn(1),bn(2),bn(3),0),1,yz1pP,src(1),50,commSplit,STATUS,ierr)
		call mpi_sendrecv(z(myZone)%f(bn(1)-1,bn(2),bn(3),0),1,yz1pM,src(1),50,z(myZone)%f(en(1),bn(2),bn(3),0),1,yz1pM,des(1),50,commSplit,STATUS,ierr)
		!----------------------------------------------------------------------
		call mpi_sendrecv(z(myZone)%f(bn(1),en(2)+1,en(3)+1,15),1,linex,yz(1),50,z(myZone)%f(bn(1),bn(2),bn(3),15),1,linex,yz(myZone),50,commSplit,STATUS,ierr)
		call mpi_sendrecv(z(myZone)%f(bn(1),bn(2)-1,en(3)+1,16),1,linex,yz(2),50,z(myZone)%f(bn(1),en(2),bn(3),16),1,linex,yz(4),50,commSplit,STATUS,ierr)
		call mpi_sendrecv(z(myZone)%f(bn(1),bn(2)-1,bn(3)-1,18),1,linex,yz(myZone),50,z(myZone)%f(bn(1),en(2),en(3),18),1,linex,yz(1),50,commSplit,STATUS,ierr)
		call mpi_sendrecv(z(myZone)%f(bn(1),en(2)+1,bn(3)-1,17),1,linex,yz(4),50,z(myZone)%f(bn(1),bn(2),en(3),17),1,linex,yz(2),50,commSplit,STATUS,ierr)
		
		call mpi_sendrecv(z(myZone)%f(en(1)+1,bn(2),en(3)+1,11),1,liney,xz(1),50,z(myZone)%f(bn(1),bn(2),bn(3),11),1,liney,xz(myZone),50,commSplit,STATUS,ierr)
		call mpi_sendrecv(z(myZone)%f(bn(1)-1,bn(2),en(3)+1,12),1,liney,xz(2),50,z(myZone)%f(en(1),bn(2),bn(3),12),1,liney,xz(4),50,commSplit,STATUS,ierr)
		call mpi_sendrecv(z(myZone)%f(bn(1)-1,bn(2),bn(3)-1,14),1,liney,xz(myZone),50,z(myZone)%f(en(1),bn(2),en(3),14),1,liney,xz(1),50,commSplit,STATUS,ierr)
		call mpi_sendrecv(z(myZone)%f(en(1)+1,bn(2),bn(3)-1,13),1,liney,xz(4),50,z(myZone)%f(bn(1),bn(2),en(3),13),1,liney,xz(2),50,commSplit,STATUS,ierr)

		call mpi_sendrecv(z(myZone)%f(en(1)+1,en(2)+1,bn(3),7 ),1,linez,xy(1),50,z(myZone)%f(bn(1),bn(2),bn(3),7 ),1,linez,xy(3),50,commSplit,STATUS,ierr)
		call mpi_sendrecv(z(myZone)%f(bn(1)-1,en(2)+1,bn(3),8 ),1,linez,xy(2),50,z(myZone)%f(en(1),bn(2),bn(3),8 ),1,linez,xy(4),50,commSplit,STATUS,ierr)
		call mpi_sendrecv(z(myZone)%f(bn(1)-1,bn(2)-1,bn(3),10),1,linez,xy(3),50,z(myZone)%f(en(1),en(2),bn(3),10),1,linez,xy(1),50,commSplit,STATUS,ierr)
		call mpi_sendrecv(z(myZone)%f(en(1)+1,bn(2)-1,bn(3),9 ),1,linez,xy(4),50,z(myZone)%f(bn(1),en(2),bn(3),9 ),1,linez,xy(2),50,commSplit,STATUS,ierr)	
	!-----------------------------------------------------------------------
		if(southId) then
			j=2 
			do k=bn(3),en(3) !2,nz-1
				do i=bn(1),en(1) !2,nx-1
					z(myZone)%f(i,j,k,3 )  = z(myZone)%f(i,j-1,k,4 )  
					z(myZone)%f(i,j,k,7 )  = z(myZone)%f(i,j-1,k,9 )  
					z(myZone)%f(i,j,k,8 )  = z(myZone)%f(i,j-1,k,10)  
					z(myZone)%f(i,j,k,15)  = z(myZone)%f(i,j-1,k,16)  
					z(myZone)%f(i,j,k,17)  = z(myZone)%f(i,j-1,k,18)  
				enddo
			enddo
		endif

		if(northId) then
			j=z(myZone)%np(2)-1 
			do k=bn(3),en(3) !2,nz-1
				do i=bn(1),en(1) !2,nx-1
					z(myZone)%f(i,j,k,4 )  = z(myZone)%f(i,j+1,k,3 )  
					z(myZone)%f(i,j,k,9 )  = z(myZone)%f(i,j+1,k,7 )  
					z(myZone)%f(i,j,k,10)  = z(myZone)%f(i,j+1,k,8 )  
					z(myZone)%f(i,j,k,16)  = z(myZone)%f(i,j+1,k,15)  
					z(myZone)%f(i,j,k,18)  = z(myZone)%f(i,j+1,k,17)  
				enddo
			enddo
		endif

		if(bottomId) then
			k=2
			do j=bn(2),en(2) !2,ny-1
				do i=bn(1),en(1) !2,nx-1
					z(myZone)%f(i,j,k,5 )  = z(myZone)%f(i,j,k-1,6 ) 
					z(myZone)%f(i,j,k,11)  = z(myZone)%f(i,j,k-1,13)  
					z(myZone)%f(i,j,k,12)  = z(myZone)%f(i,j,k-1,14)  
					z(myZone)%f(i,j,k,15)  = z(myZone)%f(i,j,k-1,17)  
					z(myZone)%f(i,j,k,16)  = z(myZone)%f(i,j,k-1,18)  
				enddo
			enddo
		endif

		if(topId) then
			k=z(myZone)%np(3)-1 
			do j=bn(2),en(2) !2,ny-1
				do i=bn(1),en(1) !2,nx-1
					z(myZone)%f(i,j,k,6 )  = z(myZone)%f(i,j,k+1,5 ) 
					z(myZone)%f(i,j,k,13)  = z(myZone)%f(i,j,k+1,11)  
					z(myZone)%f(i,j,k,14)  = z(myZone)%f(i,j,k+1,12)  
					z(myZone)%f(i,j,k,17)  = z(myZone)%f(i,j,k+1,15)  
					z(myZone)%f(i,j,k,18)  = z(myZone)%f(i,j,k+1,16)  
				enddo
			enddo
		endif
		
		if(westId) then
			i=2
			do k=bn(3),en(3) !2,nz-1
				do j=bn(2),en(2) !2,ny-1
					tmp1 = zer
					tmp2 = zer
					do a=0,q
						if (ex(a) == 	0)	tmp2 = tmp2 + z(myZone)%f(i,j,k,a)
						if (ex(a) == -1)	tmp1 = tmp1 + z(myZone)%f(i,j,k,a) 
					enddo
					
					rho = (tmp2 + 2.0*tmp1)/(1.0-uInlet) 
					z(myZone)%f  (i,j,k,1 )  = z(myZone)%f(i,j,k,kb(1 ))  + (1.0/3.0)*(rho*uInlet)
					z(myZone)%f  (i,j,k,7 )  = z(myZone)%f(i,j,k,kb(7 ))  - (quar*(z(myZone)%f(i,j,k,3)  - z(myZone)%f(i,j,k,4) )) + (1.0/6.0)*(rho*uInlet)
					z(myZone)%f  (i,j,k,9 )  = z(myZone)%f(i,j,k,kb(9 ))  + (quar*(z(myZone)%f(i,j,k,3)  - z(myZone)%f(i,j,k,4) )) + (1.0/6.0)*(rho*uInlet)
					z(myZone)%f  (i,j,k,11)  = z(myZone)%f(i,j,k,kb(11))  - (quar*(z(myZone)%f(i,j,k,5)  - z(myZone)%f(i,j,k,6) )) + (1.0/6.0)*(rho*uInlet)  
					z(myZone)%f  (i,j,k,13)  = z(myZone)%f(i,j,k,kb(13))  + (quar*(z(myZone)%f(i,j,k,5)  - z(myZone)%f(i,j,k,6) )) + (1.0/6.0)*(rho*uInlet) 
				enddo
			enddo
		endif
		
		if(eastId) then
			i=z(myZone)%np(1)-1
			do k=bn(3),en(3) !2,nz-1
				do j=bn(2),en(2) !2,ny-1
					tmp1 = zer
					tmp2 = zer
					do a=0,q
						if (ex(a) == 0)	tmp2 = tmp2 + z(myZone)%f(i,j,k,a) 
						if (ex(a) == 1)	tmp1 = tmp1 + z(myZone)%f(i,j,k,a) 
					enddo
					
					ux = (tmp2 + 2.0*tmp1)/rho0 - 1.0 
					z(myZone)%f(i,j,k,kb(1 ))  = z(myZone)%f(i,j,k,1 )  - (1.0/3.0)*(rho0*ux) 
					z(myZone)%f(i,j,k,kb(7 ))  = z(myZone)%f(i,j,k,7 )  + (quar*(z(myZone)%f(i,j,k,3)  - z(myZone)%f(i,j,k,4) )) - (1.0/6.0)*(rho0*ux)  
					z(myZone)%f(i,j,k,kb(9 ))  = z(myZone)%f(i,j,k,9 )  - (quar*(z(myZone)%f(i,j,k,3)  - z(myZone)%f(i,j,k,4) )) - (1.0/6.0)*(rho0*ux)  
					z(myZone)%f(i,j,k,kb(11))  = z(myZone)%f(i,j,k,11)  + (quar*(z(myZone)%f(i,j,k,5)  - z(myZone)%f(i,j,k,6) )) - (1.0/6.0)*(rho0*ux)    
					z(myZone)%f(i,j,k,kb(13))  = z(myZone)%f(i,j,k,13)  - (quar*(z(myZone)%f(i,j,k,5)  - z(myZone)%f(i,j,k,6) )) - (1.0/6.0)*(rho0*ux)   
				enddo
			enddo
		endif
	!-----------------------------------------------------------------------
		rhoAvg = zer
		ii=0
		
		do k=bn(3),en(3) !2,nz-1
			do j=bn(2),en(2) !2,ny-1
				do i=bn(1),en(1) !2,nx-1
				
					if(i>(z(myZone-1)%leftFace(1)) .and. i<(z(myZone-1)%rightFace(1)) .and. &
						 j>(z(myZone-1)%leftFace(2)) .and. j<(z(myZone-1)%rightFace(2)) .and. &
						 k>(z(myZone-1)%leftFace(3)) .and. k<(z(myZone-1)%rightFace(3))) cycle
				
					rho = zer
					ux	= zer 
					uy	= zer
					uz	= zer
					
					do a=0,q
						rho= rho+ z(myZone)%f(i,j,k,a) 
						ux = ux + z(myZone)%f(i,j,k,a) *ex(a)
						uy = uy + z(myZone)%f(i,j,k,a) *ey(a)
						uz = uz + z(myZone)%f(i,j,k,a) *ez(a)
					enddo
					ux  = ux/rho
					uy  = uy/rho
					uz  = uz/rho

					rhoAvg = rhoAvg + rho
					ii=ii+1
					
					do a=0,q
						tmp1				= ux*ex(a)  + uy*ey(a) + uz*ez(a)  
						tmp2				= ux**two + uy**two + uz**two   
						z(myZone)%feq(i,j,k,a)  = wt(a)*rho*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2) 
						z(myZone)%ft (i,j,k,a)  = z(myZone)%f(i,j,k,a)  - (z(myZone)%f(i,j,k,a) -z(myZone)%feq(i,j,k,a) )/z(myZone)%tau 
					enddo				
					
				enddo
			enddo
		enddo
		
		CALL MPI_Reduce(rhoAvg,rhoAvgAll,1,mpi_double_precision,mpi_sum,0,commSplit,ierr)
		CALL MPI_Reduce(ii,jj,1,mpi_integer,mpi_sum,0,commSplit,ierr)
		
		if (sphereId .and. mod(ts,freqDisp) == 0) then
			rhoAvgAll =rhoAvgAll/jj!(product(z(3)%np-2)-(z(2)%rightFace(1)-z(2)%leftFace(1)-1)*(z(2)%rightFace(2)-z(2)%leftFace(2)-1)*(z(2)%rightFace(3)-z(2)%leftFace(3)-1))
			write(*,'(A,I8,F8.4)') "zone3.Rho:",	ts, rhoAvgAll
		endif
	endif !isRed ends
!-------Transfer post-collision dist to fine block boundary------------	
	if(isGreen) then
		allocate(xIN (z(myZone)%span(1)+1)) !hardwired
		allocate(yIN (z(myZone)%span(2)+1))
		allocate(zIN (z(myZone)%span(3)+1))
		allocate(xOUT(z(myZone)%np(1)-2))
		allocate(yOUT(z(myZone)%np(2)-2))
		allocate(zOUT(z(myZone)%np(3)-2))

		allocate(fINyz (z(myZone)%span(2)+1,z(myZone)%span(3)+1))
		allocate(fINxz (z(myZone)%span(1)+1,z(myZone)%span(3)+1))
		allocate(fINxy (z(myZone)%span(1)+1,z(myZone)%span(2)+1))
		allocate(fOUTyz(z(myZone)%np(2)-2,z(myZone)%np(3)-2))
		allocate(fOUTxz(z(myZone)%np(1)-2,z(myZone)%np(3)-2))
		allocate(fOUTxy(z(myZone)%np(1)-2,z(myZone)%np(2)-2))	
	
		s=1
		do i=2,z(myZone)%np(1)-1,mRel(2)
			xIN(s)=i
			s=s+1
		enddo
		
		s=1
		do i=2,z(myZone)%np(1)-1
			xOUT(s)=i
			s=s+1
		enddo
		
		s=1
		do j=2,z(myZone)%np(2)-1,mRel(2)
			yIN(s)=j
			s=s+1
		enddo
		
		s=1
		do j=2,z(myZone)%np(2)-1
			yOUT(s)=j
			s=s+1
		enddo
		
		s=1
		do k=2,z(myZone)%np(3)-1,mRel(2)
			zIN(s)=k
			s=s+1
		enddo
		
		s=1
		do k=2,z(myZone)%np(3)-1
			zOUT(s)=k
			s=s+1
		enddo
!~ 		!-----------------------
!~ 		fINyz=5.0d0
!~ 		call spline2D(yIN,zIN,fINyz,yOUT,zOUT,fOUTyz)
!~ 		print*,fOUTyz
!~ 		stop
!~ 		!-----------------------
!~ 		open(unit=15,file="beforeInt.dat")
!~ 		open(unit=16,file="afterInt.dat")

		do a=0,q
			!-----------------------------west yz-plane-------------------------
			do k=2,z(myZone)%np(3)-1,mRel(2)
				do j=2,z(myZone)%np(2)-1,mRel(2)
					i=2
					tmp1 = z(myZone+1)%feq(z(myZone)%leftFace(1)+(i-2)/mRel(2),z(myZone)%leftFace(2)+(j-2)/mRel(2),z(myZone)%leftFace(3)+(k-2)/mRel(2),a)
					tmp2 = z(myZone+1)%ft (z(myZone)%leftFace(1)+(i-2)/mRel(2),z(myZone)%leftFace(2)+(j-2)/mRel(2),z(myZone)%leftFace(3)+(k-2)/mRel(2),a)
					z(myZone)%ftA(i,j,k,a,3) = tmp1 + (z(myZone)%tau-one)/(mRel(2)*(z(myZone+1)%tau-one))*(tmp2-tmp1)
					fINyz((j-2)/mRel(2)+1,(k-2)/mRel(2)+1) = z(myZone)%ftA(i,j,k,a,3)
				enddo
			enddo
			
			call spline2D(yIN,zIN,fINyz,yOUT,zOUT,fOUTyz)

			do k=2,z(myZone)%np(3)-1
				do j=2,z(myZone)%np(2)-1
					i=2
					s=j-1
					t=k-1
					z(myZone)%ftA(i,j,k,a,3)=fOUTyz(s,t)
				enddo
			enddo
			!-----------------------------east yz-plane-------------------------
			do k=2,z(myZone)%np(3)-1,mRel(2)
				do j=2,z(myZone)%np(2)-1,mRel(2)
					i=z(myZone)%np(1)-1
					tmp1 = z(myZone+1)%feq(z(myZone)%leftFace(1)+(i-2)/mRel(2),z(myZone)%leftFace(2)+(j-2)/mRel(2),z(myZone)%leftFace(3)+(k-2)/mRel(2),a)
					tmp2 = z(myZone+1)%ft (z(myZone)%leftFace(1)+(i-2)/mRel(2),z(myZone)%leftFace(2)+(j-2)/mRel(2),z(myZone)%leftFace(3)+(k-2)/mRel(2),a)
					z(myZone)%ftA(i,j,k,a,3) = tmp1 + (z(myZone)%tau-one)/(mRel(2)*(z(myZone+1)%tau-one))*(tmp2-tmp1)
					fINyz((j-2)/mRel(2)+1,(k-2)/mRel(2)+1) = z(myZone)%ftA(i,j,k,a,3)
				enddo
			enddo	
			call spline2D(yIN,zIN,fINyz,yOUT,zOUT,fOUTyz)
			
			do k=2,z(myZone)%np(3)-1
				do j=2,z(myZone)%np(2)-1
					i=z(myZone)%np(1)-1
					s=j-1
					t=k-1
					z(myZone)%ftA(i,j,k,a,3)=fOUTyz(s,t)
				enddo
			enddo
			!------------------------south xz-plane-----------------------------
			do k=2,z(myZone)%np(3)-1,mRel(2)
				do i=2,z(myZone)%np(1)-1,mRel(2)
					j=2
					tmp1 = z(myZone+1)%feq(z(myZone)%leftFace(1)+(i-2)/mRel(2),z(myZone)%leftFace(2)+(j-2)/mRel(2),z(myZone)%leftFace(3)+(k-2)/mRel(2),a)
					tmp2 = z(myZone+1)%ft (z(myZone)%leftFace(1)+(i-2)/mRel(2),z(myZone)%leftFace(2)+(j-2)/mRel(2),z(myZone)%leftFace(3)+(k-2)/mRel(2),a)
					z(myZone)%ftA(i,j,k,a,3) = tmp1 + (z(myZone)%tau-one)/(mRel(2)*(z(myZone+1)%tau-one))*(tmp2-tmp1)
					fINxz((i-2)/mRel(2)+1,(k-2)/mRel(2)+1) = z(myZone)%ftA(i,j,k,a,3)
				enddo
			enddo	
			call spline2D(xIN,zIN,fINxz,xOUT,zOUT,fOUTxz)
			
			do k=2,z(myZone)%np(3)-1
				do i=2,z(myZone)%np(1)-1
					j=2
					s=i-1
					t=k-1
					z(myZone)%ftA(i,j,k,a,3)=fOUTxz(s,t)
				enddo
			enddo			
			!---------------------north xz-plane--------------------------------
			do k=2,z(myZone)%np(3)-1,mRel(2)
				do i=2,z(myZone)%np(1)-1,mRel(2)
					j=z(myZone)%np(2)-1
					tmp1 = z(myZone+1)%feq(z(myZone)%leftFace(1)+(i-2)/mRel(2),z(myZone)%leftFace(2)+(j-2)/mRel(2),z(myZone)%leftFace(3)+(k-2)/mRel(2),a)
					tmp2 = z(myZone+1)%ft (z(myZone)%leftFace(1)+(i-2)/mRel(2),z(myZone)%leftFace(2)+(j-2)/mRel(2),z(myZone)%leftFace(3)+(k-2)/mRel(2),a)
					z(myZone)%ftA(i,j,k,a,3) = tmp1 + (z(myZone)%tau-one)/(mRel(2)*(z(myZone+1)%tau-one))*(tmp2-tmp1)
					fINxz((i-2)/mRel(2)+1,(k-2)/mRel(2)+1) = z(myZone)%ftA(i,j,k,a,3)
				enddo
			enddo	
			call spline2D(xIN,zIN,fINxz,xOUT,zOUT,fOUTxz)
			
			do k=2,z(myZone)%np(3)-1
				do i=2,z(myZone)%np(1)-1
					j=z(myZone)%np(2)-1
					s=i-1
					t=k-1
					z(myZone)%ftA(i,j,k,a,3)=fOUTxz(s,t)
				enddo
			enddo	
			!-----------------------bottom xy-plane-----------------------------
			do j=2,z(myZone)%np(2)-1,mRel(2)
				do i=2,z(myZone)%np(1)-1,mRel(2)
					k=2
					tmp1 = z(myZone+1)%feq(z(myZone)%leftFace(1)+(i-2)/mRel(2),z(myZone)%leftFace(2)+(j-2)/mRel(2),z(myZone)%leftFace(3)+(k-2)/mRel(2),a)
					tmp2 = z(myZone+1)%ft (z(myZone)%leftFace(1)+(i-2)/mRel(2),z(myZone)%leftFace(2)+(j-2)/mRel(2),z(myZone)%leftFace(3)+(k-2)/mRel(2),a)
					z(myZone)%ftA(i,j,k,a,3) = tmp1 + (z(myZone)%tau-one)/(mRel(2)*(z(myZone+1)%tau-one))*(tmp2-tmp1)
					fINxy((i-2)/mRel(2)+1,(j-2)/mRel(2)+1) = z(myZone)%ftA(i,j,k,a,3)
				enddo
			enddo	
			call spline2D(xIN,yIN,fINxy,xOUT,yOUT,fOUTxy)
				
			do j=2,z(myZone)%np(2)-1
				do i=2,z(myZone)%np(1)-1
					k=2
					s=i-1
					t=j-1
					z(myZone)%ftA(i,j,k,a,3)=fOUTxy(s,t)
				enddo
			enddo				
			!-------------------------top xy-plane------------------------------
			do j=2,z(myZone)%np(2)-1,mRel(2)
				do i=2,z(myZone)%np(1)-1,mRel(2)
					k=z(myZone)%np(3)-1
					tmp1 = z(myZone+1)%feq(z(myZone)%leftFace(1)+(i-2)/mRel(2),z(myZone)%leftFace(2)+(j-2)/mRel(2),z(myZone)%leftFace(3)+(k-2)/mRel(2),a)
					tmp2 = z(myZone+1)%ft (z(myZone)%leftFace(1)+(i-2)/mRel(2),z(myZone)%leftFace(2)+(j-2)/mRel(2),z(myZone)%leftFace(3)+(k-2)/mRel(2),a)
					z(myZone)%ftA(i,j,k,a,3) = tmp1 + (z(myZone)%tau-one)/(mRel(2)*(z(myZone+1)%tau-one))*(tmp2-tmp1)
					fINxy((i-2)/mRel(2)+1,(j-2)/mRel(2)+1) = z(myZone)%ftA(i,j,k,a,3)
				enddo
			enddo	
			call spline2D(xIN,yIN,fINxy,xOUT,yOUT,fOUTxy)
			
			do j=2,z(myZone)%np(2)-1
				do i=2,z(myZone)%np(1)-1
					k=z(myZone)%np(3)-1
					s=i-1
					t=j-1
					z(myZone)%ftA(i,j,k,a,3)=fOUTxy(s,t)
				enddo
			enddo	
			
		enddo
		
		deallocate(xIN )
		deallocate(yIN )
		deallocate(zIN )
		deallocate(xOUT)
		deallocate(yOUT)
		deallocate(zOUT)

		deallocate(fINyz)
		deallocate(fINxz)
		deallocate(fINxy)
		deallocate(fOUTyz)
		deallocate(fOUTxz)
		deallocate(fOUTxy)				
	!-----------------------------------------------------------------------
		do k=bn(3),en(3)!2,z(myZone)%np(3)-1 !copy ft_f(1) to ft_f_
			do j=bn(2),en(2)!2,z(myZone)%np(2)-1
				do i=bn(1),en(1)!2,z(myZone)%np(1)-1
					do a=0,q			
						z(myZone)%ft(i,j,k,a) = z(myZone)%ftA(i,j,k,a,2) 
					enddo
				enddo
			enddo
		enddo
		
		do iter2=1,mRel(2)
		
			do k=bn(3),en(3)!2,z(myZone)%np(3)-1 !streaming
				do j=bn(2),en(2)!2,z(myZone)%np(2)-1
					do i=bn(1),en(1)!2,z(myZone)%np(1)-1
						do a=0,q	
							ia = i+ex(a)  
							ja = j+ey(a)
							ka = k+ez(a)

							z(myZone)%f(ia,ja,ka,a) = z(myZone)%ft(i,j,k,a) 
						enddo
					enddo	
				enddo
			enddo
				
				!---------------------------------------------------------------------
			call mpi_sendrecv(z(myZone)%f(bn(1),bn(2),en(3)+1,0),1,xy1pP,des(3),50,z(myZone)%f(bn(1),bn(2),bn(3),0),1,xy1pP,src(3),50,commSplit,STATUS,ierr)
			call mpi_sendrecv(z(myZone)%f(bn(1),bn(2),bn(3)-1,0),1,xy1pM,src(3),50,z(myZone)%f(bn(1),bn(2),en(3),0),1,xy1pM,des(3),50,commSplit,STATUS,ierr)

			call mpi_sendrecv(z(myZone)%f(bn(1),en(2)+1,bn(3),0),1,xz1pP,des(2),50,z(myZone)%f(bn(1),bn(2),bn(3),0),1,xz1pP,src(2),50,commSplit,STATUS,ierr)
			call mpi_sendrecv(z(myZone)%f(bn(1),bn(2)-1,bn(3),0),1,xz1pM,src(2),50,z(myZone)%f(bn(1),en(2),bn(3),0),1,xz1pM,des(2),50,commSplit,STATUS,ierr)

			call mpi_sendrecv(z(myZone)%f(en(1)+1,bn(2),bn(3),0),1,yz1pP,des(1),50,z(myZone)%f(bn(1),bn(2),bn(3),0),1,yz1pP,src(1),50,commSplit,STATUS,ierr)
			call mpi_sendrecv(z(myZone)%f(bn(1)-1,bn(2),bn(3),0),1,yz1pM,src(1),50,z(myZone)%f(en(1),bn(2),bn(3),0),1,yz1pM,des(1),50,commSplit,STATUS,ierr)
			!----------------------------------------------------------------------
			call mpi_sendrecv(z(myZone)%f(bn(1),en(2)+1,en(3)+1,15),1,linex,yz(1),50,z(myZone)%f(bn(1),bn(2),bn(3),15),1,linex,yz(myZone),50,commSplit,STATUS,ierr)
			call mpi_sendrecv(z(myZone)%f(bn(1),bn(2)-1,en(3)+1,16),1,linex,yz(2),50,z(myZone)%f(bn(1),en(2),bn(3),16),1,linex,yz(4),50,commSplit,STATUS,ierr)
			call mpi_sendrecv(z(myZone)%f(bn(1),bn(2)-1,bn(3)-1,18),1,linex,yz(myZone),50,z(myZone)%f(bn(1),en(2),en(3),18),1,linex,yz(1),50,commSplit,STATUS,ierr)
			call mpi_sendrecv(z(myZone)%f(bn(1),en(2)+1,bn(3)-1,17),1,linex,yz(4),50,z(myZone)%f(bn(1),bn(2),en(3),17),1,linex,yz(2),50,commSplit,STATUS,ierr)
			
			call mpi_sendrecv(z(myZone)%f(en(1)+1,bn(2),en(3)+1,11),1,liney,xz(1),50,z(myZone)%f(bn(1),bn(2),bn(3),11),1,liney,xz(myZone),50,commSplit,STATUS,ierr)
			call mpi_sendrecv(z(myZone)%f(bn(1)-1,bn(2),en(3)+1,12),1,liney,xz(2),50,z(myZone)%f(en(1),bn(2),bn(3),12),1,liney,xz(4),50,commSplit,STATUS,ierr)
			call mpi_sendrecv(z(myZone)%f(bn(1)-1,bn(2),bn(3)-1,14),1,liney,xz(myZone),50,z(myZone)%f(en(1),bn(2),en(3),14),1,liney,xz(1),50,commSplit,STATUS,ierr)
			call mpi_sendrecv(z(myZone)%f(en(1)+1,bn(2),bn(3)-1,13),1,liney,xz(4),50,z(myZone)%f(bn(1),bn(2),en(3),13),1,liney,xz(2),50,commSplit,STATUS,ierr)

			call mpi_sendrecv(z(myZone)%f(en(1)+1,en(2)+1,bn(3),7 ),1,linez,xy(1),50,z(myZone)%f(bn(1),bn(2),bn(3),7 ),1,linez,xy(3),50,commSplit,STATUS,ierr)
			call mpi_sendrecv(z(myZone)%f(bn(1)-1,en(2)+1,bn(3),8 ),1,linez,xy(2),50,z(myZone)%f(en(1),bn(2),bn(3),8 ),1,linez,xy(4),50,commSplit,STATUS,ierr)
			call mpi_sendrecv(z(myZone)%f(bn(1)-1,bn(2)-1,bn(3),10),1,linez,xy(3),50,z(myZone)%f(en(1),en(2),bn(3),10),1,linez,xy(1),50,commSplit,STATUS,ierr)
			call mpi_sendrecv(z(myZone)%f(en(1)+1,bn(2)-1,bn(3),9 ),1,linez,xy(4),50,z(myZone)%f(bn(1),en(2),bn(3),9 ),1,linez,xy(2),50,commSplit,STATUS,ierr)	

			!----------------------------------------------------------------------
			rhoAvg = zer
			ii=0
			
			do k=b2n(3),e2n(3)!3,z(myZone)%np(3)-2
				do j=b2n(2),e2n(2)!3,z(myZone)%np(2)-2
					do i=b2n(1),e2n(1)!3,z(myZone)%np(1)-2
						rho = zer
						ux = zer 
						uy = zer
						uz = zer 
						
						do a=0,q 
							rho = rho + z(myZone)%f(i,j,k,a) 
							ux = ux + z(myZone)%f(i,j,k,a)*ex(a) 
							uy = uy + z(myZone)%f(i,j,k,a)*ey(a) 
							uz = uz + z(myZone)%f(i,j,k,a)*ez(a)  
						enddo
						ux 	= ux/rho  
						uy 	= uy/rho  
						uz 	= uz/rho  

						rhoAvg = rhoAvg + rho 
						ii=ii+1
						
						do a=0,q
							tmp1					 = ux*ex(a) + uy*ey(a) + uz*ez(a)  
							tmp2					 = ux**two + uy**two + uz**two 
							z(myZone)%feq(i,j,k,a) = wt(a)*rho*(1.0d0 + 3.0d0*tmp1 + 4.5d0*tmp1**two - 1.5d0*tmp2) 
							z(myZone)%ft(i,j,k,a) = z(myZone)%f(i,j,k,a) - (z(myZone)%f(i,j,k,a)-z(myZone)%feq(i,j,k,a))/z(myZone)%tau 
						enddo						
						
					enddo	
				enddo
			enddo
			
			CALL MPI_Reduce(rhoAvg,rhoAvgAll,1,mpi_double_precision,mpi_sum,0,commSplit,ierr)
			CALL MPI_Reduce(ii,jj,1,mpi_integer,mpi_sum,0,commSplit,ierr)
			
			if (sphereId .and. mod(ts,freqDisp) == 0) then
				rhoAvgAll =rhoAvgAll/jj
				write(*,'(A,I8,F8.4)') "zone2.Rho:",	ts, rhoAvgAll
			endif 
				
			if(iter2==mRel(2)) goto 101	
	!==============Temporal integration on fine block boundary==============
			xiL(1)=1.0d0 !(1) is (n-1)th time z
			xiL(2)=2.0d0 !(2) is (n  )th time z
			xiL(3)=3.0d0 !(3) is (n+1)th time z
!~ 			!---------debug-----------------
!~ 			yiL(1)=6!ft_f(i,j,k,a,1) 
!~ 			yiL(2)=9!ft_f(i,j,k,a,2) 
!~ 			yiL(3)=14!ft_f(i,j,k,a,3) 
!~ 			print*,lagrng3p(xiL,yiL,2.5d0)
!~ 			stop
	!-----------------------------west yz-plane-----------------------------
			if(westId) then
				do k=bn(3),en(3)!2,z(2)%np(3)-1
					do j=bn(2),en(2)!2,z(2)%np(2)-1
						i=2
						do a=0,q
							yiL(1)=z(myZone)%ftA(i,j,k,a,1) 
							yiL(2)=z(myZone)%ftA(i,j,k,a,2) 
							yiL(3)=z(myZone)%ftA(i,j,k,a,3) 
							z(myZone)%ft(i,j,k,a) = lagrng3p(xiL,yiL,two+dble(iter2)/mRel(2))
						enddo
					enddo
				enddo
			endif	
	!-----------------------------east yz-plane-----------------------------	
			if(eastId) then
				do k=bn(3),en(3)!2,z(myZone)%np(3)-1
					do j=bn(2),en(2)!2,z(myZone)%np(2)-1
						i=z(myZone)%np(1)-1
						do a=0,q
							yiL(1)=z(myZone)%ftA(i,j,k,a,1) 
							yiL(2)=z(myZone)%ftA(i,j,k,a,2) 
							yiL(3)=z(myZone)%ftA(i,j,k,a,3) 
							z(myZone)%ft(i,j,k,a) = lagrng3p(xiL,yiL,two+dble(iter2)/mRel(2))
						enddo
					enddo
				enddo
			endif
	!----------------------------south xz-plane-----------------------------			
			if(southId) then
				do k=bn(3),en(3)!2,z(myZone)%np(3)-1
					do i=bn(1),en(1)!2,z(myZone)%np(1)-1
						j=2
						do a=0,q
							yiL(1)=z(myZone)%ftA(i,j,k,a,1) 
							yiL(2)=z(myZone)%ftA(i,j,k,a,2) 
							yiL(3)=z(myZone)%ftA(i,j,k,a,3) 
							z(myZone)%ft(i,j,k,a) = lagrng3p(xiL,yiL,two+dble(iter2)/mRel(2))
						enddo
					enddo
				enddo
			endif
	!----------------------------north xz-plane-----------------------------	
			if(northId) then
				do k=bn(3),en(3)!2,z(myZone)%np(3)-1
					do i=bn(1),en(1)!2,z(myZone)%np(1)-1
						j=z(myZone)%np(2)-1
						do a=0,q
							yiL(1)=z(myZone)%ftA(i,j,k,a,1) 
							yiL(2)=z(myZone)%ftA(i,j,k,a,2) 
							yiL(3)=z(myZone)%ftA(i,j,k,a,3) 
							z(myZone)%ft(i,j,k,a) = lagrng3p(xiL,yiL,two+dble(iter2)/mRel(2))
						enddo
					enddo
				enddo
			endif	
	!----------------------------bottom xy-plane----------------------------		
			if(bottomId) then
				do j=bn(2),en(2)!2,z(myZone)%np(2)-1
					do i=bn(1),en(1)!2,z(myZone)%np(1)-1
						k=2
						do a=0,q
							yiL(1)=z(myZone)%ftA(i,j,k,a,1) 
							yiL(2)=z(myZone)%ftA(i,j,k,a,2) 
							yiL(3)=z(myZone)%ftA(i,j,k,a,3) 
							z(myZone)%ft(i,j,k,a) = lagrng3p(xiL,yiL,two+dble(iter2)/mRel(2)) 
						enddo
					enddo
				enddo	
			endif
	!------------------------------top xy-plane-----------------------------		
			if(topId) then
				do j=bn(2),en(2)!2,z(myZone)%np(2)-1
					do i=bn(1),en(1)!2,z(myZone)%np(1)-1
						k=z(myZone)%np(3)-1
						do a=0,q 
							yiL(1)=z(myZone)%ftA(i,j,k,a,1) 
							yiL(2)=z(myZone)%ftA(i,j,k,a,2) 
							yiL(3)=z(myZone)%ftA(i,j,k,a,3) 
							z(myZone)%ft(i,j,k,a) = lagrng3p(xiL,yiL,two+dble(iter2)/mRel(2)) 
						enddo
					enddo
				enddo	
			endif
			!=========Start finer m loop here========
			101 continue
			
			allocate(xIN (mRel(1)*z(1)%span(1)+1)) !hardwired
			allocate(yIN (mRel(1)*z(1)%span(2)+1))
			allocate(zIN (mRel(1)*z(1)%span(3)+1))
			allocate(xOUT(z(1)%np(1)-2))
			allocate(yOUT(z(1)%np(2)-2))
			allocate(zOUT(z(1)%np(3)-2))

			allocate(fINyz (mRel(1)*z(1)%span(2)+1,mRel(1)*z(1)%span(3)+1))
			allocate(fINxz (mRel(1)*z(1)%span(1)+1,mRel(1)*z(1)%span(3)+1))
			allocate(fINxy (mRel(1)*z(1)%span(1)+1,mRel(1)*z(1)%span(2)+1))
			allocate(fOUTyz(z(1)%np(2)-2,z(1)%np(3)-2))
			allocate(fOUTxz(z(1)%np(1)-2,z(1)%np(3)-2))
			allocate(fOUTxy(z(1)%np(1)-2,z(1)%np(2)-2))	
		
			s=1
			do i=2,z(1)%np(1)-1,mRel(1)
				xIN(s)=i
				s=s+1
			enddo
			
			s=1
			do i=2,z(1)%np(1)-1
				xOUT(s)=i
				s=s+1
			enddo
			
			s=1
			do j=2,z(1)%np(2)-1,mRel(1)
				yIN(s)=j
				s=s+1
			enddo
			
			s=1
			do j=2,z(1)%np(2)-1
				yOUT(s)=j
				s=s+1
			enddo
			
			s=1
			do k=2,z(1)%np(3)-1,mRel(1)
				zIN(s)=k
				s=s+1
			enddo
			
			s=1
			do k=2,z(1)%np(3)-1
				zOUT(s)=k
				s=s+1
			enddo
	!~ 		!-----------------------
	!~ 		fINyz=5.0d0
	!~ 		call spline2D(yIN,zIN,fINyz,yOUT,zOUT,fOUTyz)
	!~ 		print*,fOUTyz
	!~ 		stop
	!~ 		!-----------------------
	!~ 		open(unit=15,file="beforeInt.dat")
	!~ 		open(unit=16,file="afterInt.dat")

			do a=0,q
				!-----------------------------west yz-plane-------------------------
				do k=2,z(1)%np(3)-1,mRel(1)
					do j=2,z(1)%np(2)-1,mRel(1)
						i=2
						ii=(z(1)%leftFace(1)-z(2)%leftFace(1))*mRel(1)+2+(i-2)/mRel(1)
						jj=(z(1)%leftFace(2)-z(2)%leftFace(2))*mRel(1)+2+(j-2)/mRel(1)
						kk=(z(1)%leftFace(3)-z(2)%leftFace(3))*mRel(1)+2+(k-2)/mRel(1)
						tmp1 = z(2)%feq(ii,jj,kk,a)
						tmp2 = z(2)%ft (ii,jj,kk,a)
						z(1)%ftA(i,j,k,a,3) = tmp1 + (z(1)%tau-one)/(mRel(1)*(z(2)%tau-one))*(tmp2-tmp1)
						fINyz((j-2)/mRel(1)+1,(k-2)/mRel(1)+1) = z(1)%ftA(i,j,k,a,3)
					enddo
				enddo
				
				call spline2D(yIN,zIN,fINyz,yOUT,zOUT,fOUTyz)

				do k=2,z(1)%np(3)-1
					do j=2,z(1)%np(2)-1
						i=2
						s=j-1
						t=k-1
						z(1)%ftA(i,j,k,a,3)=fOUTyz(s,t)
					enddo
				enddo
				!-----------------------------east yz-plane-------------------------
				do k=2,z(1)%np(3)-1,mRel(1)
					do j=2,z(1)%np(2)-1,mRel(1)
						i=z(1)%np(1)-1
						ii=(z(1)%leftFace(1)-z(2)%leftFace(1))*mRel(1)+2+(i-2)/mRel(1)
						jj=(z(1)%leftFace(2)-z(2)%leftFace(2))*mRel(1)+2+(j-2)/mRel(1)
						kk=(z(1)%leftFace(3)-z(2)%leftFace(3))*mRel(1)+2+(k-2)/mRel(1)
						tmp1 = z(2)%feq(ii,jj,kk,a)
						tmp2 = z(2)%ft(ii,jj,kk,a)
						z(1)%ftA(i,j,k,a,3) = tmp1 + (z(1)%tau-one)/(mRel(1)*(z(2)%tau-one))*(tmp2-tmp1)
						fINyz((j-2)/mRel(1)+1,(k-2)/mRel(1)+1) = z(1)%ftA(i,j,k,a,3)
					enddo
				enddo	
				call spline2D(yIN,zIN,fINyz,yOUT,zOUT,fOUTyz)
				
				do k=2,z(1)%np(3)-1
					do j=2,z(1)%np(2)-1
						i=z(1)%np(1)-1
						s=j-1
						t=k-1
						z(1)%ftA(i,j,k,a,3)=fOUTyz(s,t)
					enddo
				enddo
				!------------------------south xz-plane-----------------------------
				do k=2,z(1)%np(3)-1,mRel(1)
					do i=2,z(1)%np(1)-1,mRel(1)
						j=2
						ii=(z(1)%leftFace(1)-z(2)%leftFace(1))*mRel(1)+2+(i-2)/mRel(1)
						jj=(z(1)%leftFace(2)-z(2)%leftFace(2))*mRel(1)+2+(j-2)/mRel(1)
						kk=(z(1)%leftFace(3)-z(2)%leftFace(3))*mRel(1)+2+(k-2)/mRel(1)
						tmp1 = z(2)%feq(ii,jj,kk,a)
						tmp2 = z(2)%ft(ii,jj,kk,a)
						z(1)%ftA(i,j,k,a,3) = tmp1 + (z(1)%tau-one)/(mRel(1)*(z(2)%tau-one))*(tmp2-tmp1)
						fINxz((i-2)/mRel(1)+1,(k-2)/mRel(1)+1) = z(1)%ftA(i,j,k,a,3)
					enddo
				enddo	
				call spline2D(xIN,zIN,fINxz,xOUT,zOUT,fOUTxz)
				
				do k=2,z(1)%np(3)-1
					do i=2,z(1)%np(1)-1
						j=2
						s=i-1
						t=k-1
						z(1)%ftA(i,j,k,a,3)=fOUTxz(s,t)
					enddo
				enddo			
				!---------------------north xz-plane--------------------------------
				do k=2,z(1)%np(3)-1,mRel(1)
					do i=2,z(1)%np(1)-1,mRel(1)
						j=z(1)%np(2)-1
						ii=(z(1)%leftFace(1)-z(2)%leftFace(1))*mRel(1)+2+(i-2)/mRel(1)
						jj=(z(1)%leftFace(2)-z(2)%leftFace(2))*mRel(1)+2+(j-2)/mRel(1)
						kk=(z(1)%leftFace(3)-z(2)%leftFace(3))*mRel(1)+2+(k-2)/mRel(1)
						tmp1 = z(2)%feq(ii,jj,kk,a)
						tmp2 = z(2)%ft(ii,jj,kk,a)
						z(1)%ftA(i,j,k,a,3) = tmp1 + (z(1)%tau-one)/(mRel(1)*(z(2)%tau-one))*(tmp2-tmp1)
						fINxz((i-2)/mRel(1)+1,(k-2)/mRel(1)+1) = z(1)%ftA(i,j,k,a,3)
					enddo
				enddo	
				call spline2D(xIN,zIN,fINxz,xOUT,zOUT,fOUTxz)
				
				do k=2,z(1)%np(3)-1
					do i=2,z(1)%np(1)-1
						j=z(1)%np(2)-1
						s=i-1
						t=k-1
						z(1)%ftA(i,j,k,a,3)=fOUTxz(s,t)
					enddo
				enddo	
				!-----------------------bottom xy-plane-----------------------------
				do j=2,z(1)%np(2)-1,mRel(1)
					do i=2,z(1)%np(1)-1,mRel(1)
						k=2
						ii=(z(1)%leftFace(1)-z(2)%leftFace(1))*mRel(1)+2+(i-2)/mRel(1)
						jj=(z(1)%leftFace(2)-z(2)%leftFace(2))*mRel(1)+2+(j-2)/mRel(1)
						kk=(z(1)%leftFace(3)-z(2)%leftFace(3))*mRel(1)+2+(k-2)/mRel(1)
						tmp1 = z(2)%feq(ii,jj,kk,a)
						tmp2 = z(2)%ft(ii,jj,kk,a)
						z(1)%ftA(i,j,k,a,3) = tmp1 + (z(1)%tau-one)/(mRel(1)*(z(2)%tau-one))*(tmp2-tmp1)
						fINxy((i-2)/mRel(1)+1,(j-2)/mRel(1)+1) = z(1)%ftA(i,j,k,a,3)
					enddo
				enddo	
				call spline2D(xIN,yIN,fINxy,xOUT,yOUT,fOUTxy)
					
				do j=2,z(1)%np(2)-1
					do i=2,z(1)%np(1)-1
						k=2
						s=i-1
						t=j-1
						z(1)%ftA(i,j,k,a,3)=fOUTxy(s,t)
					enddo
				enddo				
				!-------------------------top xy-plane------------------------------
				do j=2,z(1)%np(2)-1,mRel(1)
					do i=2,z(1)%np(1)-1,mRel(1)
						k=z(1)%np(3)-1
						ii=(z(1)%leftFace(1)-z(2)%leftFace(1))*mRel(1)+2+(i-2)/mRel(1)
						jj=(z(1)%leftFace(2)-z(2)%leftFace(2))*mRel(1)+2+(j-2)/mRel(1)
						kk=(z(1)%leftFace(3)-z(2)%leftFace(3))*mRel(1)+2+(k-2)/mRel(1)
						tmp1 = z(2)%feq(ii,jj,kk,a)
						tmp2 = z(2)%ft(ii,jj,kk,a)
						z(1)%ftA(i,j,k,a,3) = tmp1 + (z(1)%tau-one)/(mRel(1)*(z(2)%tau-one))*(tmp2-tmp1)
						fINxy((i-2)/mRel(1)+1,(j-2)/mRel(1)+1) = z(1)%ftA(i,j,k,a,3)
					enddo
				enddo	
				call spline2D(xIN,yIN,fINxy,xOUT,yOUT,fOUTxy)
				
				do j=2,z(1)%np(2)-1
					do i=2,z(1)%np(1)-1
						k=z(1)%np(3)-1
						s=i-1
						t=j-1
						z(1)%ftA(i,j,k,a,3)=fOUTxy(s,t)
					enddo
				enddo	
				
			enddo
			
			deallocate(xIN )
			deallocate(yIN )
			deallocate(zIN )
			deallocate(xOUT)
			deallocate(yOUT)
			deallocate(zOUT)

			deallocate(fINyz)
			deallocate(fINxz)
			deallocate(fINxy)
			deallocate(fOUTyz)
			deallocate(fOUTxz)
			deallocate(fOUTxy)				
		!-----------------------------------------------------------------------
			do k=2,z(1)%np(3)-1 !copy ft_f(1) to ft_f_
				do j=2,z(1)%np(2)-1
					do i=2,z(1)%np(1)-1
						do a=0,q			
							z(1)%ft(i,j,k,a) = z(1)%ftA(i,j,k,a,2) 
						enddo
					enddo
				enddo
			enddo
			
			do iter1=1,mRel(1)
			
				do k=2,z(1)%np(3)-1 !streaming
					do j=2,z(1)%np(2)-1
						do i=2,z(1)%np(1)-1
							do a=0,q	
								ia = i+ex(a)  
								ja = j+ey(a)
								ka = k+ez(a)

								z(1)%f(ia,ja,ka,a) = z(1)%ft(i,j,k,a) 
							enddo
						enddo	
					enddo
				enddo
			!-------------------------------------------
				fx(2) = zer
				fy(2) = zer
				fz(2) = zer
				
				do k=2,z(1)%np(3)-1 !Sphere BC
					do j=2,z(1)%np(2)-1
						do i=2,z(1)%np(1)-1
						
							if (isn(i,j,k) == 0) then
								do a=0,q
									ia = i+ex(a)  
									ja = j+ey(a) 
									ka = k+ez(a)
									
									if (isn(ia,ja,ka)==1) then					
										z(1)%f(i ,j ,k ,kb(a)) = z(1)%ft(i ,j, k ,a    ) 
										z(1)%f(ia,ja,ka,a    ) = z(1)%ft(ia,ja,ka,kb(a)) 
										
										fx(2) = fx(2) + ex(a)*2.0*(-z(1)%ft(ia,ja,ka,kb(a)) + z(1)%ft(i,j,k,a))/(m(1)**two) 						
										fy(2) = fy(2) + ey(a)*2.0*(-z(1)%ft(ia,ja,ka,kb(a)) + z(1)%ft(i,j,k,a))/(m(1)**two) 						
										fz(2) = fz(2) + ez(a)*2.0*(-z(1)%ft(ia,ja,ka,kb(a)) + z(1)%ft(i,j,k,a))/(m(1)**two) 
									endif
								enddo
							endif
							
						enddo	
					enddo
				enddo

				!----------------------------------------------------------------------
				fx_t 	= haf*(fx(1)  +  fx(2)) 
				fy_t 	= haf*(fy(1)  +  fy(2)) 
				fz_t 	= haf*(fz(1)  +  fz(2)) 
				
				Cd = fx_t/(haf*rho0*uInlet*uInlet*pi*D**two/4.0d0) 
				Cs = fy_t/(haf*rho0*uInlet*uInlet*pi*D**two/4.0d0) 
				Cl = fz_t/(haf*rho0*uInlet*uInlet*pi*D**two/4.0d0) 
				!----------------------------------------------------------------------
				fx(1) = fx(2) 
				fy(1) = fy(2)
				fz(1) = fz(2)
				!----------------------------------------------------------------------
				rhoAvg = zer

				do k=3,z(1)%np(3)-2
					do j=3,z(1)%np(2)-2
						do i=3,z(1)%np(1)-2
							rho = zer
							ux = zer 
							uy = zer
							uz = zer 
							
							do a=0,q 
								rho= rho+ z(1)%f(i,j,k,a) 
								ux = ux + z(1)%f(i,j,k,a)*ex(a) 
								uy = uy + z(1)%f(i,j,k,a)*ey(a) 
								uz = uz + z(1)%f(i,j,k,a)*ez(a)  
							enddo 
							ux 	= ux/rho  
							uy 	= uy/rho  
							uz 	= uz/rho  

							rhoAvg = rhoAvg + rho 
							
							do a=0,q
								tmp1					 = ux*ex(a) + uy*ey(a) + uz*ez(a)  
								tmp2					 = ux**two + uy**two + uz**two 
								z(1)%feq(i,j,k,a) = wt(a)*rho*(1.0d0 + 3.0d0*tmp1 + 4.5d0*tmp1**two - 1.5d0*tmp2) 
								z(1)%ft(i,j,k,a) = z(1)%f(i,j,k,a) - (z(1)%f(i,j,k,a)-z(1)%feq(i,j,k,a))/z(1)%tau 
							enddo						
							
						enddo	
					enddo
				enddo
				
				rhoAvg = rhoAvg/((z(1)%np(1)-4)*(z(1)%np(2)-4)*(z(1)%np(3)-4)) 
				
		!--------------------------------------------------------------------		
				if (mod(ts,freqDisp) == 0 .and. iter1==mRel(1) .and. iter2==mRel(2)) then
					write(10,'(I8,1X,4(F10.4,1X))') ts, rhoAvg, Cd, Cs, Cl
				endif	 
					
				if(iter1==mRel(1)) exit	
		!==============Temporal integration on fine block boundary==============
				xiL(1)=1.0d0 !(1) is (n-1)th time z
				xiL(2)=2.0d0 !(2) is (n  )th time z
				xiL(3)=3.0d0 !(3) is (n+1)th time z
	!~ 			!---------debug-----------------
	!~ 			yiL(1)=6!ft_f(i,j,k,a,1) 
	!~ 			yiL(2)=9!ft_f(i,j,k,a,2) 
	!~ 			yiL(3)=14!ft_f(i,j,k,a,3) 
	!~ 			print*,lagrng3p(xiL,yiL,2.5d0)
	!~ 			stop
		!-----------------------------west yz-plane-----------------------------
				do k=2,z(1)%np(3)-1
					do j=2,z(1)%np(2)-1
						i=2
						do a=0,q
							yiL(1)=z(1)%ftA(i,j,k,a,1) 
							yiL(2)=z(1)%ftA(i,j,k,a,2) 
							yiL(3)=z(1)%ftA(i,j,k,a,3) 
							z(1)%ft(i,j,k,a) = lagrng3p(xiL,yiL,two+dble(iter1)/mRel(1))
						enddo
					enddo
				enddo	
		!-----------------------------east yz-plane-----------------------------	
				do k=2,z(1)%np(3)-1
					do j=2,z(1)%np(2)-1
						i=z(1)%np(1)-1
						do a=0,q
							yiL(1)=z(1)%ftA(i,j,k,a,1) 
							yiL(2)=z(1)%ftA(i,j,k,a,2) 
							yiL(3)=z(1)%ftA(i,j,k,a,3) 
							z(1)%ft(i,j,k,a) = lagrng3p(xiL,yiL,two+dble(iter1)/mRel(1))
						enddo
					enddo
				enddo	
		!----------------------------south xz-plane-----------------------------			
				do k=2,z(1)%np(3)-1
					do i=2,z(1)%np(1)-1
						j=2
						do a=0,q
							yiL(1)=z(1)%ftA(i,j,k,a,1) 
							yiL(2)=z(1)%ftA(i,j,k,a,2) 
							yiL(3)=z(1)%ftA(i,j,k,a,3) 
							z(1)%ft(i,j,k,a) = lagrng3p(xiL,yiL,two+dble(iter1)/mRel(1))
						enddo
					enddo
				enddo	
		!----------------------------north xz-plane-----------------------------	
				do k=2,z(1)%np(3)-1
					do i=2,z(1)%np(1)-1
						j=z(1)%np(2)-1
						do a=0,q
							yiL(1)=z(1)%ftA(i,j,k,a,1) 
							yiL(2)=z(1)%ftA(i,j,k,a,2) 
							yiL(3)=z(1)%ftA(i,j,k,a,3) 
							z(1)%ft(i,j,k,a) = lagrng3p(xiL,yiL,two+dble(iter1)/mRel(1))
						enddo
					enddo
				enddo			
		!----------------------------bottom xy-plane----------------------------		
				do j=2,z(1)%np(2)-1
					do i=2,z(1)%np(1)-1
						k=2
						do a=0,q
							yiL(1)=z(1)%ftA(i,j,k,a,1) 
							yiL(2)=z(1)%ftA(i,j,k,a,2) 
							yiL(3)=z(1)%ftA(i,j,k,a,3) 
							z(1)%ft(i,j,k,a) = lagrng3p(xiL,yiL,two+dble(iter1)/mRel(1)) 
						enddo
					enddo
				enddo	
		!------------------------------top xy-plane-----------------------------		
				do j=2,z(1)%np(2)-1
					do i=2,z(1)%np(1)-1
						k=z(1)%np(3)-1
						do a=0,q 
							yiL(1)=z(1)%ftA(i,j,k,a,1) 
							yiL(2)=z(1)%ftA(i,j,k,a,2) 
							yiL(3)=z(1)%ftA(i,j,k,a,3) 
							z(1)%ft(i,j,k,a) = lagrng3p(xiL,yiL,two+dble(iter1)/mRel(1)) 
						enddo
					enddo
				enddo	
				
			enddo!iter1 loop ends
			!---------------------------------------------
			do k=3,z(1)%np(3)-2
				do j=3,z(1)%np(2)-2
					do i=3,z(1)%np(1)-2
						do a=0,q
							z(1)%ftA(i,j,k,a,3)	= z(1)%ft(i,j,k,a) 
						enddo
					enddo	
				enddo
			enddo

		!=========Tranfer post-collision dist on coarse block boundary==========
		!-----------------------------west yz-plane-----------------------------	
			do k=2+mRel(1),z(1)%np(3)-1-mRel(1),mRel(1)
				do j=2+mRel(1),z(1)%np(2)-1-mRel(1),mRel(1)
					i=2+mRel(1)
					do a=0,q
						ii=(z(1)%leftFace(1)  -z(2)%leftFace(1))*mRel(1)+2  +(i-2)/mRel(1)
						jj=(z(1)%leftFace(2) -z(2)%leftFace(2))*mRel(1)+2 +(j-2)/mRel(1)
						kk=(z(1)%leftFace(3)-z(2)%leftFace(3))*mRel(1)+2+(k-2)/mRel(1)
						z(2)%ft(ii,jj,kk,a) = z(1)%feq(i,j,k,a) + mRel(1)*(z(2)%tau-1)/(z(1)%tau-1)*(z(1)%ftA(i,j,k,a,3)-z(1)%feq(i,j,k,a)) 
						enddo
				enddo	
			enddo
		!-----------------------------east yz-plane-----------------------------		
			do k=2+mRel(1),z(1)%np(3)-1-mRel(1),mRel(1)
				do j=2+mRel(1),z(1)%np(2)-1-mRel(1),mRel(1)
					i=z(1)%np(1)-1-mRel(1)
					do a=0,q
						ii=(z(1)%leftFace(1)  -z(2)%leftFace(1))*mRel(1)+2  +(i-2)/mRel(1)
						jj=(z(1)%leftFace(2) -z(2)%leftFace(2))*mRel(1)+2 +(j-2)/mRel(1)
						kk=(z(1)%leftFace(3)-z(2)%leftFace(3))*mRel(1)+2+(k-2)/mRel(1)
						z(2)%ft(ii,jj,kk,a) = z(1)%feq(i,j,k,a) + mRel(1)*(z(2)%tau-1)/(z(1)%tau-1)*(z(1)%ftA(i,j,k,a,3)-z(1)%feq(i,j,k,a)) 
					enddo
				enddo
			enddo	
		!-----------------------------south xz-plane-----------------------------	
			do k=2+mRel(1),z(1)%np(3)-1-mRel(1),mRel(1)
				do i=2+mRel(1),z(1)%np(1)-1-mRel(1),mRel(1)
					j=2+mRel(1)
					do a=0,q
						ii=(z(1)%leftFace(1)  -z(2)%leftFace(1))*mRel(1)+2  +(i-2)/mRel(1)
						jj=(z(1)%leftFace(2) -z(2)%leftFace(2))*mRel(1)+2 +(j-2)/mRel(1)
						kk=(z(1)%leftFace(3)-z(2)%leftFace(3))*mRel(1)+2+(k-2)/mRel(1)
						z(2)%ft(ii,jj,kk,a) = z(1)%feq(i,j,k,a) + mRel(1)*(z(2)%tau-1)/(z(1)%tau-1)*(z(1)%ftA(i,j,k,a,3)-z(1)%feq(i,j,k,a))
 					enddo
				enddo	
			enddo
		!-----------------------------north xz-plane----------------------------		
			do k=2+mRel(1),z(1)%np(3)-1-mRel(1),mRel(1)
				do i=2+mRel(1),z(1)%np(1)-1-mRel(1),mRel(1)
					j=z(1)%np(2)-1-mRel(1)
					do a=0,q
						ii=(z(1)%leftFace(1)-z(2)%leftFace(1))*mRel(1)+2+(i-2)/mRel(1)
						jj=(z(1)%leftFace(2)-z(2)%leftFace(2))*mRel(1)+2+(j-2)/mRel(1)
						kk=(z(1)%leftFace(3)-z(2)%leftFace(3))*mRel(1)+2+(k-2)/mRel(1)
						z(2)%ft(ii,jj,kk,a) = z(1)%feq(i,j,k,a) + mRel(1)*(z(2)%tau-1)/(z(1)%tau-1)*(z(1)%ftA(i,j,k,a,3)-z(1)%feq(i,j,k,a)) 
					enddo
				enddo
			enddo	
		!-----------------------------bottom xy-plane---------------------------		
			do j=2+mRel(1),z(1)%np(2)-1-mRel(1),mRel(1)	
				do i=2+mRel(1),z(1)%np(1)-1-mRel(1),mRel(1)
					k=2+mRel(1)
					do a=0,q
						ii=(z(1)%leftFace(1)  -z(2)%leftFace(1))*mRel(1)+2  +(i-2)/mRel(1)
						jj=(z(1)%leftFace(2) -z(2)%leftFace(2))*mRel(1)+2 +(j-2)/mRel(1)
						kk=(z(1)%leftFace(3)-z(2)%leftFace(3))*mRel(1)+2+(k-2)/mRel(1)
						z(2)%ft(ii,jj,kk,a) = z(1)%feq(i,j,k,a) + mRel(1)*(z(2)%tau-1)/(z(1)%tau-1)*(z(1)%ftA(i,j,k,a,3)-z(1)%feq(i,j,k,a)) 
					enddo
				enddo
			enddo	
		!------------------------------top xy-plane-----------------------------		
			do j=2+mRel(1),z(1)%np(2)-1-mRel(1),mRel(1)	
				do i=2+mRel(1),z(1)%np(1)-1-mRel(1),mRel(1)
					k=z(1)%np(3)-1-mRel(1)
					do a=0,q
						ii=(z(1)%leftFace(1)  -z(2)%leftFace(1))*mRel(1)+2  +(i-2)/mRel(1)
						jj=(z(1)%leftFace(2) -z(2)%leftFace(2))*mRel(1)+2 +(j-2)/mRel(1)
						kk=(z(1)%leftFace(3)-z(2)%leftFace(3))*mRel(1)+2+(k-2)/mRel(1)
						z(2)%ft(ii,jj,kk,a) = z(1)%feq(i,j,k,a) + mRel(1)*(z(2)%tau-1)/(z(1)%tau-1)*(z(1)%ftA(i,j,k,a,3)-z(1)%feq(i,j,k,a)) 
					enddo
				enddo
			enddo	
		!-----------------------------------------------------------------------
			do k=2,z(1)%np(3)-1
				do j=2,z(1)%np(2)-1
					do i=2,z(1)%np(1)-1
						do a=0,q
							z(1)%ftA(i,j,k,a,1) = z(1)%ftA(i,j,k,a,2)
							z(1)%ftA(i,j,k,a,2) = z(1)%ftA(i,j,k,a,3)
						enddo
					enddo	
				enddo
			enddo			
				
			!==========End finer m loop here=========
						
		enddo!iter2 loop ends
		!---------------------------------------------
		do k=3,z(2)%np(3)-2
			do j=3,z(2)%np(2)-2
				do i=3,z(2)%np(1)-2
					do a=0,q
						z(2)%ftA(i,j,k,a,3)	= z(2)%ft(i,j,k,a) 
					enddo
				enddo	
			enddo
		enddo

	!=========Tranfer post-collision dist on coarse block boundary==========
	!-----------------------------west yz-plane-----------------------------	
		do k=2+m(2),z(2)%np(3)-1-m(2),m(2)
			do j=2+m(2),z(2)%np(2)-1-m(2),m(2)
				i=2+m(2)
				do a=0,q
					z(3)%ft(z(2)%leftFace(1)+(i-2)/m(2),z(2)%leftFace(2)+(j-2)/m(2),z(2)%leftFace(3)+(k-2)/m(2),a) = z(2)%feq(i,j,k,a) + m(2)*(z(3)%tau-1)/(z(2)%tau-1)*(z(2)%ftA(i,j,k,a,3)-z(2)%feq(i,j,k,a)) 
				enddo
			enddo	
		enddo
	!-----------------------------east yz-plane-----------------------------		
		do k=2+m(2),z(2)%np(3)-1-m(2),m(2)
			do j=2+m(2),z(2)%np(2)-1-m(2),m(2)
				i=z(2)%np(1)-1-m(2)
				do a=0,q
					z(3)%ft(z(2)%leftFace(1)+(i-2)/m(2),z(2)%leftFace(2)+(j-2)/m(2),z(2)%leftFace(3)+(k-2)/m(2),a) = z(2)%feq(i,j,k,a) + m(2)*(z(3)%tau-1)/(z(2)%tau-1)*(z(2)%ftA(i,j,k,a,3)-z(2)%feq(i,j,k,a)) 
				enddo
			enddo
		enddo	
	!-----------------------------south xz-plane-----------------------------	
		do k=2+m(2),z(2)%np(3)-1-m(2),m(2)
			do i=2+m(2),z(2)%np(1)-1-m(2),m(2)
				j=2+m(2)
				do a=0,q
					z(3)%ft(z(2)%leftFace(1)+(i-2)/m(2),z(2)%leftFace(2)+(j-2)/m(2),z(2)%leftFace(3)+(k-2)/m(2),a) = z(2)%feq(i,j,k,a) + m(2)*(z(3)%tau-1)/(z(2)%tau-1)*(z(2)%ftA(i,j,k,a,3)-z(2)%feq(i,j,k,a)) 
				enddo
			enddo	
		enddo
	!-----------------------------north xz-plane----------------------------		
		do k=2+m(2),z(2)%np(3)-1-m(2),m(2)
			do i=2+m(2),z(2)%np(1)-1-m(2),m(2)
				j=z(2)%np(2)-1-m(2)
				do a=0,q
					z(3)%ft(z(2)%leftFace(1)+(i-2)/m(2),z(2)%leftFace(2)+(j-2)/m(2),z(2)%leftFace(3)+(k-2)/m(2),a) = z(2)%feq(i,j,k,a) + m(2)*(z(3)%tau-1)/(z(2)%tau-1)*(z(2)%ftA(i,j,k,a,3)-z(2)%feq(i,j,k,a)) 
				enddo
			enddo
		enddo	
	!-----------------------------bottom xy-plane---------------------------		
		do j=2+m(2),z(2)%np(2)-1-m(2),m(2)	
			do i=2+m(2),z(2)%np(1)-1-m(2),m(2)
				k=2+m(2)
				do a=0,q
					z(3)%ft(z(2)%leftFace(1)+(i-2)/m(2),z(2)%leftFace(2)+(j-2)/m(2),z(2)%leftFace(3)+(k-2)/m(2),a) = z(2)%feq(i,j,k,a) + m(2)*(z(3)%tau-1)/(z(2)%tau-1)*(z(2)%ftA(i,j,k,a,3)-z(2)%feq(i,j,k,a)) 
				enddo
			enddo
		enddo	
	!------------------------------top xy-plane-----------------------------		
		do j=2+m(2),z(2)%np(2)-1-m(2),m(2)	
			do i=2+m(2),z(2)%np(1)-1-m(2),m(2)
				k=z(2)%np(3)-1-m(2)
				do a=0,q
					z(3)%ft(z(2)%leftFace(1)+(i-2)/m(2),z(2)%leftFace(2)+(j-2)/m(2),z(2)%leftFace(3)+(k-2)/m(2),a) = z(2)%feq(i,j,k,a) + m(2)*(z(3)%tau-1)/(z(2)%tau-1)*(z(2)%ftA(i,j,k,a,3)-z(2)%feq(i,j,k,a)) 
				enddo
			enddo
		enddo	
	!-----------------------------------------------------------------------
		do k=2,z(2)%np(3)-1
			do j=2,z(2)%np(2)-1
				do i=2,z(2)%np(1)-1
					do a=0,q
						z(2)%ftA(i,j,k,a,1) = z(2)%ftA(i,j,k,a,2)
						z(2)%ftA(i,j,k,a,2) = z(2)%ftA(i,j,k,a,3)
					enddo
				enddo	
			enddo
		enddo
	
	endif !if(sphereId) ends	
!-----------------------------------------------------------------------
	if(ts == snapList(solNum+1)) then !mod(ts,time/(noOfSnaps-1))
	
		solNum = solNum + 1
!=======================================================================
		write(filename,'(2(A,I2.2),A)') "S",solNum,"_",id+1,".dat"
		open(unit=12,file=filename)
		
		write(12,'(A,I8,A)') 'Title="',ts,'"'
		write(12,'(A)') "Variables=x,y,z,region,u,v,w,rho"
	
		if(sphereId) then
!-------------------------------------1---------------------------------
			write(12,'(3(A,I6))') "Zone I= ",en(1)-bn(1)+1,	&
															 " ,J= ",en(2)-bn(2)+1,			&
															 " ,K= ",en(3)-bn(3)+1

			do k=bn(3),en(3) !2,nz-1
				do j=bn(2),en(2) !2,ny-1
					do i=bn(1),en(1)!bn(1),westFace
					
						tmp1 = zer
						tmp2 = zer 
						tmp3 = zer
						tmp4 = zer 
						
						do a=0,q 
							tmp1 = tmp1 + z(3)%f(i,j,k,a) 
							tmp2 = tmp2 + z(3)%f(i,j,k,a)*ex(a) 
							tmp3 = tmp3 + z(3)%f(i,j,k,a)*ey(a) 
							tmp4 = tmp4 + z(3)%f(i,j,k,a)*ez(a)  
						enddo

						rho	= tmp1  
						ux 	= tmp2/tmp1  
						uy 	= tmp3/tmp1  
						uz 	= tmp4/tmp1 					
					
						write(12,'(4(I4,1X),3(F8.4,1X),F15.8)') i,j,k,0,ux/uInlet,uy/uInlet,uz/uInlet,(rho-rho0)*100.0d0 
					enddo
					write(12,*)
				enddo
				write(12,*)
			enddo
!-----------------------------------2-----------------------------------	
			write(12,'(3(A,I6))') "Zone I= ",z(2)%np(1)-2,	&
															 " ,J= ",z(2)%np(2)-2,	&
															 " ,K= ",z(2)%np(3)-2

			do k=2,z(2)%np(3)-1
				do j=2,z(2)%np(2)-1
					do i=2,z(2)%np(1)-1
					
						tmp1 = zer
						tmp2 = zer 
						tmp3 = zer
						tmp4 = zer 
						
						do a=0,q 
							tmp1 = tmp1 + z(2)%f(i,j,k,a) 
							tmp2 = tmp2 + z(2)%f(i,j,k,a)*ex(a) 
							tmp3 = tmp3 + z(2)%f(i,j,k,a)*ey(a) 
							tmp4 = tmp4 + z(2)%f(i,j,k,a)*ez(a)  
						enddo

						rho	= tmp1  
						ux 	= tmp2/tmp1  
						uy 	= tmp3/tmp1  
						uz 	= tmp4/tmp1 					
					
						write(12,'(3(F8.2,1X),I4,1X,3(F8.4,1X),F15.8)') (z(2)%leftFace(1)+dble(i-2)/m(2)),(z(2)%leftFace(2)+dble(j-2)/m(2)),(z(2)%leftFace(3)+dble(k-2)/m(2)),&
																															0,ux/uInlet,uy/uInlet,uz/uInlet,(rho-rho0)*100.0d0 
					enddo
					write(12,*)
				enddo
				write(12,*)
			enddo
!-----------------------------------3-----------------------------------	
			write(12,'(3(A,I6))') "Zone I= ",z(1)%np(1)-2,	&
															 " ,J= ",z(1)%np(2)-2,	&
															 " ,K= ",z(1)%np(3)-2

			do k=2,z(1)%np(3)-1
				do j=2,z(1)%np(2)-1
					do i=2,z(1)%np(1)-1
					
						tmp1 = zer
						tmp2 = zer 
						tmp3 = zer
						tmp4 = zer 
						
						do a=0,q 
							tmp1 = tmp1 + z(1)%f(i,j,k,a) 
							tmp2 = tmp2 + z(1)%f(i,j,k,a)*ex(a) 
							tmp3 = tmp3 + z(1)%f(i,j,k,a)*ey(a) 
							tmp4 = tmp4 + z(1)%f(i,j,k,a)*ez(a)  
						enddo

						rho	= tmp1  
						ux 	= tmp2/tmp1  
						uy 	= tmp3/tmp1  
						uz 	= tmp4/tmp1 					
					
						write(12,'(3(F8.2,1X),I4,1X,3(F8.4,1X),F15.8)') (z(1)%leftFace(1)+dble(i-2)/m(1)),&
																														(z(1)%leftFace(2)+dble(j-2)/m(1)),&
																														(z(1)%leftFace(3)+dble(k-2)/m(1)),&
																															isn(i,j,k),ux/uInlet,uy/uInlet,uz/uInlet,(rho-rho0)*100.0d0 
					enddo
					write(12,*)
				enddo
				write(12,*)
			enddo
							
		else
!-------------------------------------0---------------------------------
			write(12,'(3(A,I6))') "Zone I= ",en(1)-bn(1)+1,	&
															 " ,J= ",en(2)-bn(2)+1,	&
															 " ,K= ",en(3)-bn(3)+1

			do k=bn(3),en(3)
				do j=bn(2),en(2)
					do i=bn(1),en(1)
					
						tmp1 = zer
						tmp2 = zer 
						tmp3 = zer
						tmp4 = zer 
						
						do a=0,q 
							tmp1 = tmp1 + z(3)%f(i,j,k,a) 
							tmp2 = tmp2 + z(3)%f(i,j,k,a)*ex(a) 
							tmp3 = tmp3 + z(3)%f(i,j,k,a)*ey(a) 
							tmp4 = tmp4 + z(3)%f(i,j,k,a)*ez(a)  
						enddo

						rho	= tmp1  
						ux 	= tmp2/tmp1  
						uy 	= tmp3/tmp1  
						uz 	= tmp4/tmp1 					
					
						write(12,'(4(I4,1X),3(F8.4,1X),F15.8)') i,j,k,0,ux/uInlet,uy/uInlet,uz/uInlet,(rho-rho0)*100.0d0
					enddo
					write(12,*)
				enddo
				write(12,*)
			enddo
		
		endif						
!-----------------------------------------------------------------------		
		close(12)
!-------------------------------------------------------------------
		call mpi_gather(bs(1),3,mpi_integer,bb(1,0),3,mpi_integer,0,mpi_comm_world,ierr)
		call mpi_gather(es(1),3,mpi_integer,ee(1,0),3,mpi_integer,0,mpi_comm_world,ierr)
!-----------------------Write Domain Info---------------------------
		if(masterId) then
			open(unit=12,file="domInfo.txt")
			write(12,'(3(i4))') z(3)%np(1),z(3)%np(2),z(3)%np(3)
			write(12,'(i3)') nproc
			do i=0,nproc-1
				write(12,'(6(i4,1x))') bb(1,i),ee(1,i),bb(2,i),ee(2,i),bb(3,i),ee(3,i)
			enddo
			close(12)	
		endif
!=======================================================================		
		if(masterId) write(*,'(A,I8)') "Solution recorded at time ",ts
 	endif
!----------------------------------------------------------------------
	if (sphereId .and. mod(ts,freqDisp) == 0) then
		call cpu_time(tEnd)
		write(*,'(I3,I8,2X,A,F8.2,A,I4,A)') id, ts," | CPU Time: ",tEnd-tStart," sec for ",freqDisp," timesteps"
	endif	

enddo!Time loop Ends

if(sphereId) close(10) 

contains
function lagrng3p(x,y,xp) result(yp)
implicit none

double precision,dimension(:),intent(in)::x,y
double precision,intent(in)::xp 
integer:: k,j,p=3
double precision:: yp,tmp
	
yp =0.0d0
do k=1,p
	tmp = 1.0d0
	do j=1,p
		if(j==k) cycle
		tmp = tmp*(xp-x(j))/(x(k)-x(j))
	enddo
	yp = yp + tmp*y(k)
enddo

end function lagrng3p

subroutine spline2D(xi,yi,zi,xo,yo,zo)
!Developer	:Parvez Ahmad<pahmed333@gmail.com>
!Code for 2D spline interpolation
!INPUT
!	xi	:x-vector of known data
!	yi	:y-vector of known data
!	zi	:known data on xi-yi grid
!	xo	:x-vector of uknown data
! yo	:y-vector of uknown data
!OUTPUT
!	zo	:unknown data on xo-yo grid
implicit none

double precision,dimension(:),intent(in):: xi,yi,xo,yo
double precision,dimension(:,:),intent(in):: zi
double precision,dimension(:,:),intent(out):: zo

integer:: i,j,k,p,iu,iv,nSx,nSy
double precision,allocatable,dimension(:)::ax,bx,cx,dx,solX
double precision,allocatable,dimension(:)::ay,by,cy,dy,solY
double precision,allocatable,dimension(:,:,:,:)::C
double precision:: u,v,lhs

nSx=size(xi)-1
nSy=size(yi)-1

allocate(ax(nSx),bx(nSx+1),cx(nSx),dx(nSx+1),solX(nSx+1))
allocate(ay(nSy),by(nSy+1),cy(nSy),dy(nSy+1),solY(nSy+1))
allocate(C(nSx+1,nSy+1,4,4))

!-----------SET 1---------------
do i=1,nSx
	ax(i)=1.0;
	bx(i)=4.0;
	cx(i)=1.0;
enddo

ax(nSx)=0.0;
bx(1)=1.0; bx(nSx+1)=1.0;
cx(1)=0.0;

do j=1,nSy+1
	do i=2,nSx
		dx(i)=3.0*( zi(i+1,j) - 2*zi(i,j) + zi(i-1,j) );
	enddo
	dx(1)		=0.0;
	dx(nSx+1)	=0.0;
	call TDMA(bx,ax,cx,solX,dx)
	
	do i=1,nSx+1 
		C(i,j,3,1)=solX(i);
		!printf("%6.4f ",solX(i));
	enddo
	!printf("\n");
	do i=1,nSx
		C(i,j,1,1)=zi(i,j);
		C(i,j,2,1)=zi(i+1,j)-zi(i,j)-C(i,j,3,1)-(C(i+1,j,3,1)-C(i,j,3,1))/3.0;
		C(i,j,4,1)=(C(i+1,j,3,1)-C(i,j,3,1))/3.0;
	enddo
enddo

!-------SET 5---------
do j=1,nSy
	ay(j)=1.0;
	by(j)=4.0;
	cy(j)=1.0;
enddo
ay(nSy)=0.0;
by(1)=1.0; by(nSy+1)=1.0;
cy(1)=0.0;

do i=1,nSx+1
	do j=2,nSy
		dy(j)=3*( zi(i,j+1) - 2*zi(i,j) + zi(i,j-1) );
	enddo
	dy(1)		=0.0;
	dy(nSy+1)	=0.0;
	call TDMA(by,ay,cy,solY,dy)
	
	do j=1,nSy+1
		C(i,j,1,3)=solY(j)
	enddo
	
	do j=1,nSy
		C(i,j,1,2)=zi(i,j+1)-zi(i,j)-C(i,j,1,3)-(C(i,j+1,1,3)-C(i,j,1,3))/3.0;
		C(i,j,1,4)=(C(i,j+1,1,3)-C(i,j,1,3))/3.0;
	enddo
enddo
!-------SET 3---------
do i=1,nSx
	ax(i)=1.0;
	bx(i)=4.0;
	cx(i)=1.0;
enddo
ax(nSx)=0.0;
bx(1)=1.0; bx(nSx+1)=1.0;
cx(1)=0.0;

do j=1,nSy
	do i=2,nSx
		dx(i)=3*(  C(i+1,j,1,3) - 2*C(i,j,1,3) + C(i-1,j,1,3)  );
	enddo
	dx(1)		=0.0;
	dx(nSx+1)	=0.0;
	call TDMA(bx,ax,cx,solX,dx)
	
	do i=1,nSx+1 
		C(i,j,3,3)=solX(i)
	enddo
	
	do i=1,nSx
		C(i,j,2,3)=C(i+1,j,1,3)-C(i,j,1,3)-C(i,j,3,3)-(C(i+1,j,3,3)-C(i,j,3,3))/3.0;
		C(i,j,4,3)=(C(i+1,j,3,3)-C(i,j,3,3))/3.0;
	enddo
enddo
!-----------------------SET 7; C32,C34--------------------------------
do i=1,nSx
	do j=1,nSy
		C(i,j,3,2)=C(i,j+1,3,1)-C(i,j,3,1)-C(i,j,3,3)-(C(i,j+1,3,3)-C(i,j,3,3))/3.0;
		C(i,j,3,4)=(C(i,j+1,3,3)-C(i,j,3,3))/3.0;
	enddo
enddo
! -----------------------SET 4; C24,C44--------------------------------
do i=1,nSx
	do j=1,nSy
		C(i,j,2,4)=(C(i,j+1,2,3)-C(i,j,2,3))/3.0;
		C(i,j,4,4)=(C(i,j+1,4,3)-C(i,j,4,3))/3.0;
	enddo
enddo
!------------------------SET 8; C42-----------------------------------
do i=1,nSx
	do j=1,nSy
		C(i,j,4,2)=(C(i+1,j,3,2)-C(i,j,3,2))/3.0;
	enddo
enddo
! -----------------------SET 2; C22------------------------------------
do i=1,nSx-1
	do j=1,nSy
		C(i,j,2,2)=C(i+1,j,1,2)-C(i,j,1,2)-C(i,j,3,2)-C(i,j,4,2);
	enddo
enddo

do j=1,nSy-1
	C(nSx,j,2,2)=C(nSx,j+1,2,1)-C(nSx,j,2,1)-C(nSx,j,2,3)-C(nSx,j,2,4);
enddo

lhs=0.0;
do i=1,4
	do j=1,4
		if(i==2 .and. j==2) cycle
		lhs = lhs + C(nSx,nSy,i,j);
	enddo
enddo
C(nSx,nSy,2,2)=zi(nSx+1,nSy+1)-lhs;
!---------------------------------------------------
!~ do j=1,nSy+1
!~ 	do i=1,nSx+1
!~ 		write(*,'(I3,I3)') i,j
!~ 		do k=1,4 
!~ 			write(*,'(4(F8.4))') C(i,j,k,1),C(i,j,k,2),C(i,j,k,3),C(i,j,k,4)
!~ 		enddo
!~ 	write(*,*) "------------------------------------------"
!~ 	enddo
!~ enddo
!-----------------Interpolant Evaluation on Fine Mesh-----------------
!zo=0.0d0
do i=1,nSx
	do j=1,nSy
		do k=1,size(xo)
			do p=1,size(yo)
				if(xo(k) .ge. xi(i) .and. xo(k) .le. xi(i+1) .and. yo(p) .ge. yi(j) .and. yo(p) .le. yi(j+1)) then
					zo(k,p)=0.0d0
					u=(xo(k)-xi(i))/(xi(i+1)-xi(i))
					v=(yo(p)-yi(j))/(yi(j+1)-yi(j))
					
					do iu=0,3
						do iv=0,3
							zo(k,p) = zo(k,p) + (u**iu)*C(i,j,iu+1,iv+1)*(v**iv)  !  (1 u u^2 u^3)*C(:,:,i,j)*(1 v v^2 v^3)';
						enddo
					enddo
				endif
			enddo
		enddo
	enddo
enddo
!~ write(*,*) "=================================================="
!~ do i=1,nSx
!~ 	do j=1,nSy
!~ 		write(*,'(I3,I3,F8.4)') i,j,zo(i,j)
!~ 	enddo
!~ 	write(*,*)
!~ enddo

end subroutine spline2D


SUBROUTINE TDMA(d,d1,d2,x,b)
!Developer	:Parvez Ahmad <pahmed333@gmail.com>
!Solves a tridiagonal matrix using TDMA
!INPUT
!	d	:array of diagonal elements
!	d1:array d1 of sub-diagonal elements
!	d2:array of super-diagonal elements
!	b	:array of rhs elements 
!OUTPUT
!	x	:solution vector

IMPLICIT NONE
DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::D,D1,D2,B
DOUBLE PRECISION,DIMENSION(:),INTENT(OUT)	::X
INTEGER							::I,N
DOUBLE PRECISION				::M(SIZE(B),4),A(SIZE(B)),B1(SIZE(B))
N=SIZE(B)

M(2:N,1)=D1
M(:,2)=D
M(1:N-1,3)=D2
M(:,4)=B

!THOMAS ALGORITHM STARTS
A(1)=M(1,2)
B1(1)=M(1,4)/M(1,2);
DO I=2,N
A(I)=M(I,2)-(M(I,1)*M(I-1,3)/A(I-1));
B1(I)=(M(I,4)-M(I,1)*B1(I-1))/A(I);
END DO

X(N)=B1(N);
DO I=N-1,1,-1
X(I)=B1(I)-(M(I,3)*X(I+1)/A(I));
END DO

END SUBROUTINE TDMA

end program main
