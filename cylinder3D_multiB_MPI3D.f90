program main
use mpi
implicit none
!-----------------------------------------------------------------------
integer,parameter,dimension(3)::	nP			=[300,80,180],	&
																	topo		=[2,1,3]
																	
integer,parameter,dimension(4)::	snapList=[1,10000,50000,100000]
																	
logical,parameter,dimension(3)::	period	=[.false. , .true. , .false.]
logical,parameter::								isDebug	=.true.

integer,parameter:: q=18				,&
										time=100000	,&
										freqDisp=1	,&
										cx_c=110		,&
										!cy_c=nP(2)/2,&
										cz_c=nP(3)/2,&
										
										m=2					,&
										nx_f=79		,&
										ny_f=nP(2)-2+(nP(2)-1)*(m-1)	,&
										nz_f=79		,&
					
										cx=nx_f/2+1	,&
										!cy=ny_f/2+1	,&
										cz=nz_f/2+1
										
double precision,parameter:: 	D=15.0d0			,&
															tau_c=0.56d0	,&
															rho0=1.0d0		,&
															uInlet=0.08d0	,&
															
															haf=0.5d0			,&
															one=1.0d0			,&
															two=2.0d0			,&
															zer=0.0d0			,&
															quar=0.25d0		,&
															one36th=1.0d0/36.0d0	,&
															pi=4.0d0*datan(one)
															
!------------------------MPI Variables----------------------------------
integer,dimension(3)		:: bs,es,bn,en,siz,pP,des,src,myDim,myCoord
integer		:: linex, liney, linez, yz1p, xz1p, xy1p, yz1pP, yz1pM, xz1pP, xz1pM, xy1pP, xy1pM
!integer		:: xz1pFine, xz1pFineP, xz1pFineM
integer,dimension(5)		:: blockLen, xpDisp, xmDisp, ypDisp, ymDisp, zpDisp, zmDisp
integer		:: comm3d,nproc, id, ierr, sizDP, STATUS(mpi_status_size),yz(4),xz(4),xy(4)
logical 	:: myperiod(3), bottomId, topId, eastId, westId, northId, southId, masterId, cylinderId
integer,allocatable,dimension(:,:) 	::tArr,bb,ee															
!-----------------------------------------------------------------------
integer:: i,j,k,a,a1,ts,ia,ja,ka,iter,s,t,nxx,nyy,nzz,nxx_f,nyy_f,nzz_f,solNum
integer:: ex(0:q),ey(0:q),ez(0:q),kb(0:q),tmpI1, tmpI2, counter, counterAll
double precision:: tmp1, tmp2, tmp3, tmp4, fx_t, fy_t, fz_t, fx(2), fy(2), fz(2), tStart, tEnd
double precision:: Cd, Cl, Cs, wt(0:q), rhoAvg, rhoAvgAll, xiL(3), yiL(3)
!coarse variables
double precision,allocatable,dimension(:,:,:)		::ux_c, uy_c, uz_c, rho_c
double precision,allocatable,dimension(:,:,:,:)	::f_c, ft_c, feq_c
!fine variables
integer																					::westFace, eastFace, southFace, northFace, bottomFace, topFace 
integer,allocatable,dimension(:,:)							::isn
double precision							::tau_f
double precision,allocatable,dimension(:)					::xIN, yIN, zIN, xOUT, yOUT, zOUT
double precision,allocatable,dimension(:,:)				::fINyz,fINxz,fINxy,fOUTyz,fOUTxz,fOUTxy
double precision,allocatable,dimension(:,:,:)			::ux_f, uy_f, uz_f, rho_f
double precision,allocatable,dimension(:,:,:,:)		::f_f, ft_f_, feq_f
double precision,allocatable,dimension(:,:,:,:,:)	::ft_f

character(len=50)	:: filename
!==============================MPI======================================
CALL mpi_init(ierr)
CALL mpi_comm_rank(mpi_comm_world,id,ierr)
CALL mpi_comm_size(mpi_comm_world,nproc,ierr)

call mpi_cart_create(mpi_comm_world,3,topo,period,.false.,comm3d,ierr)
call mpi_cart_get(comm3d,3,myDim,myperiod,mycoord,ierr)

call mpi_cart_shift(comm3d,0,1,src(1),des(1),ierr)
call mpi_cart_shift(comm3d,1,1,src(2),des(2),ierr)
call mpi_cart_shift(comm3d,2,1,src(3),des(3),ierr)

allocate(tArr(0:maxval(myDim)-1,3))
allocate(bb(3,0:nproc-1),ee(3,0:nproc-1))

masterId=(id==0)
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

	tmpI1=nP(k)/myDim(k)
	tmpI2=myDim(k)-mod(nP(k),myDim(k))

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
allocate(f_c  (bn(1)-1:en(1)+1, bn(2)-1:en(2)+1, bn(3)-1:en(3)+1, 0:q))
allocate(ft_c (bn(1)-1:en(1)+1, bn(2)-1:en(2)+1, bn(3)-1:en(3)+1, 0:q))
allocate(feq_c(bn(1)-1:en(1)+1, bn(2)-1:en(2)+1, bn(3)-1:en(3)+1, 0:q))
allocate(ux_c (bs(1):es(1), bs(2):es(2), bs(3):es(3)))
allocate(uy_c (bs(1):es(1), bs(2):es(2), bs(3):es(3)))
allocate(uz_c (bs(1):es(1), bs(2):es(2), bs(3):es(3)))
allocate(rho_c(bs(1):es(1), bs(2):es(2), bs(3):es(3)))

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
!---------This works if cylinder lies within one processor--------------
cylinderId = .false.

outermost: do k=bn(3),en(3)
	!do j=bn(2),en(2)
		do i=bn(1),en(1)
		
			tmp1 = ((i-cx_c)**two + (k-cz_c)**two)**haf 

			if (tmp1 .le. haf*D) then
				cylinderId = .true.
				exit outermost
			endif
			
		enddo
	!enddo
enddo outermost

if(isDebug .and. cylinderId) write(*,'(A,I3)') "Cylinder processor ",id
!-----------------------------------------------------------------------
if(cylinderId) then
	allocate(xIN (nx_f/2)) !hardwired
	allocate(yIN (ny_f/2))
	allocate(zIN (nz_f/2))
	allocate(xOUT(nx_f-2))
	allocate(yOUT(ny_f-2))
	allocate(zOUT(nz_f-2))

	allocate(fINyz (ny_f/2,nz_f/2))
	!allocate(fINxz (nx_f/2,nz_f/2))
	allocate(fINxy (nx_f/2,ny_f/2))
	allocate(fOUTyz(ny_f-2,nz_f-2))
	!allocate(fOUTxz(nx_f-2,nz_f-2))
	allocate(fOUTxy(nx_f-2,ny_f-2))

	allocate(ux_f (nx_f,ny_f,nz_f))
	allocate(uy_f (nx_f,ny_f,nz_f))
	allocate(uz_f (nx_f,ny_f,nz_f))
	allocate(rho_f(nx_f,ny_f,nz_f))

	allocate(f_f  (nx_f,ny_f,nz_f,0:q))
	allocate(ft_f_(nx_f,ny_f,nz_f,0:q))
	allocate(feq_f(nx_f,ny_f,nz_f,0:q))

	allocate(ft_f(nx_f,ny_f,nz_f,0:q,3))
	allocate(isn(nx_f,nz_f))
endif
!-----------------------------------------------------------------------
tau_f = haf + m*(tau_c - haf)

westFace		=cx_c-(nx_f+1)/(2*m)
eastFace		=cx_c+(nx_f+1)/(2*m)
!southFace		=cy_c-(ny_f+1)/(2*m)
!northFace		=cy_c+(ny_f+1)/(2*m)
bottomFace	=cz_c-(nz_f+1)/(2*m)
topFace			=cz_c+(nz_f+1)/(2*m)

if(isDebug .and. cylinderId) write(*,'(A,4(1X,I4))') "Face: west,east,bottom,top: ",westFace,eastFace,bottomFace,topFace
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
do k=bs(3),es(3) !1,nz
	do j=bs(2),es(2) !1,ny 
		do i=bs(1),es(1) !1,nx 
			!if(i>westFace+1 .and. i<eastFace-1 .and. j>bottomFace+1 .and. j<topFace-1) continue 
			do a=0,q
				tmp1				= uInlet*ex(a)
				tmp2				= uInlet*uInlet
				f_c (i,j,k,a)	= wt(a)*rho0*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
				ft_c(i,j,k,a)	= wt(a)*rho0*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
			enddo
		enddo
	enddo
enddo
!-----------------------------------
if(cylinderId) then
	do k=1,nz_f
		do j=1,ny_f
			do i=1,nx_f
				
				do a=0,q
					tmp1				= uInlet*ex(a)  
					tmp2				= uInlet*uInlet  
					f_f (i,j,k,a  )	= wt(a)*rho0*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
					ft_f(i,j,k,a,2) = wt(a)*rho0*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
					ft_f(i,j,k,a,1) = wt(a)*rho0*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
				enddo
				
			enddo
		enddo
	enddo
endif
!----------------------------------------------------------------------
if(cylinderId) then
!---------------------------Detect cylinder-----------------------------
	do k=1,nz_f
		!do j=2,ny_f-1
			do i=1,nx_f
			
				tmp1 = ((i-cx)**two + (k-cz)**two)**haf 

				if (tmp1 .le. m*haf*D) then
					isn(i,k) = 1
				else
					isn(i,k) = 0
				endif
				
			enddo
		!enddo
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

	if(cylinderId .and. mod(ts-1,freqDisp)==0) call cpu_time(tStart)
	
	do k=bn(3),en(3) !2,nz-1 !Streaming
		do j=bn(2),en(2) !2,ny-1
			do i=bn(1),en(1) !2,nx-1
			
				if(i>(westFace  +2) .and. i<(eastFace -2) .and. &
					 !j>(southFace +2) .and. j<(northFace-2) .and. &
					 k>(bottomFace+2) .and. k<(topFace  -2)) cycle
					 
				do a=0,q
					ia = i+ex(a) 
					ja = j+ey(a)
					ka = k+ez(a)
					
					f_c(ia,ja,ka,a) = ft_c(i,j,k,a)
				enddo
			enddo
		enddo
	enddo
	!---------------------------------------------------------------------
	call mpi_sendrecv(f_c(bn(1),bn(2),en(3)+1,0),1,xy1pP,des(3),50,f_c(bn(1),bn(2),bn(3),0),1,xy1pP,src(3),50,mpi_comm_world,STATUS,ierr)
	call mpi_sendrecv(f_c(bn(1),bn(2),bn(3)-1,0),1,xy1pM,src(3),50,f_c(bn(1),bn(2),en(3),0),1,xy1pM,des(3),50,mpi_comm_world,STATUS,ierr)

	call mpi_sendrecv(f_c(bn(1),en(2)+1,bn(3),0),1,xz1pP,des(2),50,f_c(bn(1),bn(2),bn(3),0),1,xz1pP,src(2),50,mpi_comm_world,STATUS,ierr)
	call mpi_sendrecv(f_c(bn(1),bn(2)-1,bn(3),0),1,xz1pM,src(2),50,f_c(bn(1),en(2),bn(3),0),1,xz1pM,des(2),50,mpi_comm_world,STATUS,ierr)

	call mpi_sendrecv(f_c(en(1)+1,bn(2),bn(3),0),1,yz1pP,des(1),50,f_c(bn(1),bn(2),bn(3),0),1,yz1pP,src(1),50,mpi_comm_world,STATUS,ierr)
	call mpi_sendrecv(f_c(bn(1)-1,bn(2),bn(3),0),1,yz1pM,src(1),50,f_c(en(1),bn(2),bn(3),0),1,yz1pM,des(1),50,mpi_comm_world,STATUS,ierr)
	!----------------------------------------------------------------------
	call mpi_sendrecv(f_c(bn(1),en(2)+1,en(3)+1,15),1,linex,yz(1),50,f_c(bn(1),bn(2),bn(3),15),1,linex,yz(3),50,mpi_comm_world,STATUS,ierr)
	call mpi_sendrecv(f_c(bn(1),bn(2)-1,en(3)+1,16),1,linex,yz(2),50,f_c(bn(1),en(2),bn(3),16),1,linex,yz(4),50,mpi_comm_world,STATUS,ierr)
	call mpi_sendrecv(f_c(bn(1),bn(2)-1,bn(3)-1,18),1,linex,yz(3),50,f_c(bn(1),en(2),en(3),18),1,linex,yz(1),50,mpi_comm_world,STATUS,ierr)
	call mpi_sendrecv(f_c(bn(1),en(2)+1,bn(3)-1,17),1,linex,yz(4),50,f_c(bn(1),bn(2),en(3),17),1,linex,yz(2),50,mpi_comm_world,STATUS,ierr)
	
	call mpi_sendrecv(f_c(en(1)+1,bn(2),en(3)+1,11),1,liney,xz(1),50,f_c(bn(1),bn(2),bn(3),11),1,liney,xz(3),50,mpi_comm_world,STATUS,ierr)
	call mpi_sendrecv(f_c(bn(1)-1,bn(2),en(3)+1,12),1,liney,xz(2),50,f_c(en(1),bn(2),bn(3),12),1,liney,xz(4),50,mpi_comm_world,STATUS,ierr)
	call mpi_sendrecv(f_c(bn(1)-1,bn(2),bn(3)-1,14),1,liney,xz(3),50,f_c(en(1),bn(2),en(3),14),1,liney,xz(1),50,mpi_comm_world,STATUS,ierr)
	call mpi_sendrecv(f_c(en(1)+1,bn(2),bn(3)-1,13),1,liney,xz(4),50,f_c(bn(1),bn(2),en(3),13),1,liney,xz(2),50,mpi_comm_world,STATUS,ierr)

	call mpi_sendrecv(f_c(en(1)+1,en(2)+1,bn(3),7 ),1,linez,xy(1),50,f_c(bn(1),bn(2),bn(3),7 ),1,linez,xy(3),50,mpi_comm_world,STATUS,ierr)
	call mpi_sendrecv(f_c(bn(1)-1,en(2)+1,bn(3),8 ),1,linez,xy(2),50,f_c(en(1),bn(2),bn(3),8 ),1,linez,xy(4),50,mpi_comm_world,STATUS,ierr)
	call mpi_sendrecv(f_c(bn(1)-1,bn(2)-1,bn(3),10),1,linez,xy(3),50,f_c(en(1),en(2),bn(3),10),1,linez,xy(1),50,mpi_comm_world,STATUS,ierr)
	call mpi_sendrecv(f_c(en(1)+1,bn(2)-1,bn(3),9 ),1,linez,xy(4),50,f_c(bn(1),en(2),bn(3),9 ),1,linez,xy(2),50,mpi_comm_world,STATUS,ierr)	
!-----------------------------------------------------------------------
!~ 	if(southId) then
!~ 		j=2 
!~ 		do k=bn(3),en(3) !2,nz-1
!~ 			do i=bn(1),en(1) !2,nx-1
!~ 				f_c(i,j,k,3 ) = f_c(i,j-1,k,4 ) 
!~ 				f_c(i,j,k,7 ) = f_c(i,j-1,k,9 ) 
!~ 				f_c(i,j,k,8 ) = f_c(i,j-1,k,10) 
!~ 				f_c(i,j,k,15) = f_c(i,j-1,k,16) 
!~ 				f_c(i,j,k,17) = f_c(i,j-1,k,18) 
!~ 			enddo
!~ 		enddo
!~ 	endif

!~ 	if(northId) then
!~ 		j=nP(2)-1 
!~ 		do k=bn(3),en(3) !2,nz-1
!~ 			do i=bn(1),en(1) !2,nx-1
!~ 				f_c(i,j,k,4 ) = f_c(i,j+1,k,3 ) 
!~ 				f_c(i,j,k,9 ) = f_c(i,j+1,k,7 ) 
!~ 				f_c(i,j,k,10) = f_c(i,j+1,k,8 ) 
!~ 				f_c(i,j,k,16) = f_c(i,j+1,k,15) 
!~ 				f_c(i,j,k,18) = f_c(i,j+1,k,17) 
!~ 			enddo
!~ 		enddo
!~ 	endif

	if(bottomId) then
		k=2
		do j=bn(2),en(2) !2,ny-1
			do i=bn(1),en(1) !2,nx-1
				f_c(i,j,k,5 ) = f_c(i,j,k-1,6 )
				f_c(i,j,k,11) = f_c(i,j,k-1,13) 
				f_c(i,j,k,12) = f_c(i,j,k-1,14) 
				f_c(i,j,k,15) = f_c(i,j,k-1,17) 
				f_c(i,j,k,16) = f_c(i,j,k-1,18) 
			enddo
		enddo
	endif

	if(topId) then
		k=nP(3)-1 
		do j=bn(2),en(2) !2,ny-1
			do i=bn(1),en(1) !2,nx-1
				f_c(i,j,k,6 ) = f_c(i,j,k+1,5 )
				f_c(i,j,k,13) = f_c(i,j,k+1,11) 
				f_c(i,j,k,14) = f_c(i,j,k+1,12) 
				f_c(i,j,k,17) = f_c(i,j,k+1,15) 
				f_c(i,j,k,18) = f_c(i,j,k+1,16) 
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
					if (ex(a) == 	0)	tmp2 = tmp2 + f_c(i,j,k,a)
					if (ex(a) == -1)	tmp1 = tmp1 + f_c(i,j,k,a) 
				enddo
				
				rho_c(i,j,k   ) = (tmp2 + 2.0*tmp1)/(1.0-uInlet) 
				f_c  (i,j,k,1 ) = f_c(i,j,k,kb(1 )) + (1.0/3.0)*(rho_c(i,j,k)*uInlet)
				f_c  (i,j,k,7 ) = f_c(i,j,k,kb(7 )) - (quar*(f_c(i,j,k,3) - f_c(i,j,k,4))) + (1.0/6.0)*(rho_c(i,j,k)*uInlet)
				f_c  (i,j,k,9 ) = f_c(i,j,k,kb(9 )) + (quar*(f_c(i,j,k,3) - f_c(i,j,k,4))) + (1.0/6.0)*(rho_c(i,j,k)*uInlet)
				f_c  (i,j,k,11) = f_c(i,j,k,kb(11)) - (quar*(f_c(i,j,k,5) - f_c(i,j,k,6))) + (1.0/6.0)*(rho_c(i,j,k)*uInlet)  
				f_c  (i,j,k,13) = f_c(i,j,k,kb(13)) + (quar*(f_c(i,j,k,5) - f_c(i,j,k,6))) + (1.0/6.0)*(rho_c(i,j,k)*uInlet) 
			enddo
		enddo
	endif
	
	if(eastId) then
		i=nP(1)-1
		do k=bn(3),en(3) !2,nz-1
			do j=bn(2),en(2) !2,ny-1
				tmp1 = zer
				tmp2 = zer
				do a=0,q
					if (ex(a) == 0)	tmp2 = tmp2 + f_c(i,j,k,a) 
					if (ex(a) == 1)	tmp1 = tmp1 + f_c(i,j,k,a) 
				enddo
				
				ux_c(i,j,k       ) = (tmp2 + 2.0*tmp1)/rho0 - 1.0 
				f_c (i,j,k,kb(1 )) = f_c(i,j,k,1 ) - (1.0/3.0)*(rho0*ux_c(i,j,k)) 
				f_c (i,j,k,kb(7 )) = f_c(i,j,k,7 ) + (quar*(f_c(i,j,k,3) - f_c(i,j,k,4))) - (1.0/6.0)*(rho0*ux_c(i,j,k))  
				f_c (i,j,k,kb(9 )) = f_c(i,j,k,9 ) - (quar*(f_c(i,j,k,3) - f_c(i,j,k,4))) - (1.0/6.0)*(rho0*ux_c(i,j,k))  
				f_c (i,j,k,kb(11)) = f_c(i,j,k,11) + (quar*(f_c(i,j,k,5) - f_c(i,j,k,6))) - (1.0/6.0)*(rho0*ux_c(i,j,k))    
				f_c (i,j,k,kb(13)) = f_c(i,j,k,13) - (quar*(f_c(i,j,k,5) - f_c(i,j,k,6))) - (1.0/6.0)*(rho0*ux_c(i,j,k))   
			enddo
		enddo
	endif
!-----------------------------------------------------------------------
	rhoAvg	= zer
	counter = 0
	
	do k=bn(3),en(3) !2,nz-1
		do j=bn(2),en(2) !2,ny-1
			do i=bn(1),en(1) !2,nx-1
			
				if(i>(westFace  +1) .and. i<(eastFace -1) .and. &
					 !j>(southFace +1) .and. j<(northFace-1) .and. &
					 k>(bottomFace+1) .and. k<(topFace  -1)) cycle
			
				tmp1 = zer
				tmp2 = zer 
				tmp3 = zer
				tmp4 = zer
				
				do a=0,q
					tmp1 = tmp1 + f_c(i,j,k,a)
					tmp2 = tmp2 + f_c(i,j,k,a)*ex(a)
					tmp3 = tmp3 + f_c(i,j,k,a)*ey(a)
					tmp4 = tmp4 + f_c(i,j,k,a)*ez(a)
				enddo

				rho_c(i,j,k) = tmp1
				ux_c (i,j,k) = tmp2/tmp1
				uy_c (i,j,k) = tmp3/tmp1
				uz_c (i,j,k) = tmp4/tmp1

				rhoAvg	= rhoAvg + tmp1
				counter = counter + 1
			enddo
		enddo
	enddo
	
	CALL MPI_Reduce(rhoAvg,rhoAvgAll,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
	CALL MPI_Reduce(counter,counterAll,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
	
	if (masterId .and. mod(ts,freqDisp) == 0) then
		rhoAvgAll =rhoAvgAll/counterAll!(product(nP-2)-(eastFace-westFace-3)*(northFace-southFace-3)*(topFace-bottomFace-3))
		write(*,'(A,I8,F8.4)') "Rho(Coarse) ",	ts, rhoAvgAll
	endif
!----------------------------------------------------------------------
	do k=bn(3),en(3) !2,nz-1 !coarse block Collision
		do j=bn(2),en(2) !2,ny-1
			do i=bn(1),en(1) !2,nx-1
			
				if(i>(westFace  +1) .and. i<(eastFace -1) .and. &
					 !j>(southFace +1) .and. j<(northFace-1) .and. &
					 k>(bottomFace+1) .and. k<(topFace  -1)) cycle
					 
				do a=0,q
					tmp1				= ux_c(i,j,k)*ex(a)  + uy_c(i,j,k)*ey(a) + uz_c(i,j,k)*ez(a)  
					tmp2				= ux_c(i,j,k)**two + uy_c(i,j,k)**two + uz_c(i,j,k)**two   
					feq_c(i,j,k,a)= wt(a)*rho_c(i,j,k)*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2) 
					ft_c (i,j,k,a)= f_c(i,j,k,a) - (f_c(i,j,k,a)-feq_c(i,j,k,a))/tau_c 
				enddo
			enddo	
		enddo
	enddo
!-------Transfer post-collision dist to fine block boundary------------	
	if(cylinderId) then
		s=1
		do i=2,nx_f-1,m
			xIN(s)=i
			s=s+1
		enddo
		
		s=1
		do i=2,nx_f-1
			xOUT(s)=i
			s=s+1
		enddo
		
		s=1
		do j=2,ny_f-1,m
			yIN(s)=j
			s=s+1
		enddo
		
		s=1
		do j=2,ny_f-1
			yOUT(s)=j
			s=s+1
		enddo
		
		s=1
		do k=2,nz_f-1,m
			zIN(s)=k
			s=s+1
		enddo
		
		s=1
		do k=2,nz_f-1
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
			do k=2,nz_f-1,m
				do j=2,ny_f-1,m
					i=2
					tmp1 = feq_c(westFace+i/m,1+j/m,bottomFace+k/m,a)
					tmp2 = ft_c (westFace+i/m,1+j/m,bottomFace+k/m,a)
					ft_f(i,j,k,a,3) = tmp1 + (tau_f-one)/(m*(tau_c-one))*(tmp2-tmp1)
					fINyz(j/m,k/m) = ft_f(i,j,k,a,3)
				enddo
			enddo
			
!~ 			if(a==8) then
!~ 				open(unit=15,file="interp1.dat")
!~ 				do k=2,nz_f-1
!~ 					do j=2,ny_f-1
!~ 						i=2
!~ 						write(15,'(4(I3,1X),F10.6)') i,j,k,a,ft_f(i,j,k,a,3)
!~ 					enddo
!~ 				enddo
!~ 				close(15)
!~ 			endif
			
			call spline2D(yIN,zIN,fINyz,yOUT,zOUT,fOUTyz)
			
			do k=3,nz_f-2,m
				do j=3,ny_f-2,m
					i=2
					s=j-1
					t=k-1
					ft_f(i,j  ,k  ,a,3)=fOUTyz(s  ,t  )
					ft_f(i,j-1,k  ,a,3)=fOUTyz(s-1,t  )
					ft_f(i,j+1,k  ,a,3)=fOUTyz(s+1,t  )
					ft_f(i,j  ,k-1,a,3)=fOUTyz(s  ,t-1)
					ft_f(i,j  ,k+1,a,3)=fOUTyz(s  ,t+1)
				enddo
			enddo
			
!~ 			if(a==8) then
!~ 				open(unit=15,file="interp2.dat")
!~ 				do k=2,nz_f-1
!~ 					do j=2,ny_f-1
!~ 						i=2
!~ 						write(15,'(4(I3,1X),F10.6)') i,j,k,a,ft_f(i,j,k,a,3)
!~ 					enddo
!~ 				enddo
!~ 				close(15)
!~ 				stop
!~ 			endif
			!-----------------------------east yz-plane-------------------------
			do k=2,nz_f-1,m
				do j=2,ny_f-1,m
					i=nx_f-1
					tmp1 = feq_c(westFace+i/m,1+j/m,bottomFace+k/m,a)
					tmp2 = ft_c (westFace+i/m,1+j/m,bottomFace+k/m,a)
					ft_f(i,j,k,a,3) = tmp1 + (tau_f-one)/(m*(tau_c-one))*(tmp2-tmp1)
					fINyz(j/m,k/m) = ft_f(i,j,k,a,3)
				enddo
			enddo	
			call spline2D(yIN,zIN,fINyz,yOUT,zOUT,fOUTyz)
			
			do k=3,nz_f-2,m
				do j=3,ny_f-2,m
					i=nx_f-1
					s=j-1
					t=k-1
					ft_f(i,j  ,k  ,a,3)=fOUTyz(s  ,t  )
					ft_f(i,j-1,k  ,a,3)=fOUTyz(s-1,t  )
					ft_f(i,j+1,k  ,a,3)=fOUTyz(s+1,t  )
					ft_f(i,j  ,k-1,a,3)=fOUTyz(s  ,t-1)
					ft_f(i,j  ,k+1,a,3)=fOUTyz(s  ,t+1)
				enddo
			enddo
!~ 			!------------------------south xz-plane-----------------------------
!~ 			do k=2,nz_f-1,m
!~ 				do i=2,nx_f-1,m
!~ 					j=2
!~ 					tmp1 = feq_c(westFace+i/m, southFace+j/m, bottomFace+k/m, a)
!~ 					tmp2 = ft_c (westFace+i/m, southFace+j/m, bottomFace+k/m, a)
!~ 					ft_f(i,j,k,a,3) = tmp1 + (tau_f-one)/(m*(tau_c-one))*(tmp2-tmp1)
!~ 					fINxz(i/m,k/m) = ft_f(i,j,k,a,3)
!~ 				enddo
!~ 			enddo	
!~ 			call spline2D(xIN,zIN,fINxz,xOUT,zOUT,fOUTxz)
			
!~ 			do k=3,nz_f-2,m
!~ 				do i=3,nx_f-2,m
!~ 					j=2
!~ 					s=i-1
!~ 					t=k-1
!~ 					ft_f(i  ,j,k  ,a,3)=fOUTxz(s  ,t  )
!~ 					ft_f(i-1,j,k  ,a,3)=fOUTxz(s-1,t  )
!~ 					ft_f(i+1,j,k  ,a,3)=fOUTxz(s+1,t  )
!~ 					ft_f(i  ,j,k-1,a,3)=fOUTxz(s  ,t-1)
!~ 					ft_f(i  ,j,k+1,a,3)=fOUTxz(s  ,t+1)
!~ 				enddo
!~ 			enddo			
!~ 			!---------------------north xz-plane--------------------------------
!~ 			do k=2,nz_f-1,m
!~ 				do i=2,nx_f-1,m
!~ 					j=ny_f-1
!~ 					tmp1 = feq_c(westFace+i/m, southFace+j/m, bottomFace+k/m, a)
!~ 					tmp2 = ft_c (westFace+i/m, southFace+j/m, bottomFace+k/m, a)
!~ 					ft_f(i,j,k,a,3) = tmp1 + (tau_f-one)/(m*(tau_c-one))*(tmp2-tmp1)
!~ 					fINxz(i/m,k/m) = ft_f(i,j,k,a,3)
!~ 				enddo
!~ 			enddo	
!~ 			call spline2D(xIN,zIN,fINxz,xOUT,zOUT,fOUTxz)
			
!~ 			do k=3,nz_f-2,m
!~ 				do i=3,nx_f-2,m
!~ 					j=ny_f-1
!~ 					s=i-1
!~ 					t=k-1
!~ 					ft_f(i  ,j,k  ,a,3)=fOUTxz(s  ,t  )
!~ 					ft_f(i-1,j,k  ,a,3)=fOUTxz(s-1,t  )
!~ 					ft_f(i+1,j,k  ,a,3)=fOUTxz(s+1,t  )
!~ 					ft_f(i  ,j,k-1,a,3)=fOUTxz(s  ,t-1)
!~ 					ft_f(i  ,j,k+1,a,3)=fOUTxz(s  ,t+1)
!~ 				enddo
!~ 			enddo	
			!-----------------------bottom xy-plane-----------------------------
			do j=2,ny_f-1,m
				do i=2,nx_f-1,m
					k=2
					tmp1 = feq_c(westFace+i/m,1+j/m,bottomFace+k/m,a)
					tmp2 = ft_c (westFace+i/m,1+j/m,bottomFace+k/m,a)
					ft_f(i,j,k,a,3) = tmp1 + (tau_f-one)/(m*(tau_c-one))*(tmp2-tmp1)
					fINxy(i/m,j/m) = ft_f(i,j,k,a,3)
				enddo
			enddo	
			call spline2D(xIN,yIN,fINxy,xOUT,yOUT,fOUTxy)
				
			do j=3,ny_f-2,m
				do i=3,nx_f-2,m
					k=2
					s=i-1
					t=j-1
					ft_f(i  ,j  ,k,a,3)=fOUTxy(s  ,t  )
					ft_f(i-1,j  ,k,a,3)=fOUTxy(s-1,t  )
					ft_f(i+1,j  ,k,a,3)=fOUTxy(s+1,t  )
					ft_f(i  ,j-1,k,a,3)=fOUTxy(s  ,t-1)
					ft_f(i  ,j+1,k,a,3)=fOUTxy(s  ,t+1)
				enddo
			enddo				
			!-------------------------top xy-plane------------------------------
			do j=2,ny_f-1,m
				do i=2,nx_f-1,m
					k=nz_f-1
					tmp1 = feq_c(westFace+i/m,1+j/m,bottomFace+k/m,a)
					tmp2 = ft_c (westFace+i/m,1+j/m,bottomFace+k/m,a)
					ft_f(i,j,k,a,3) = tmp1 + (tau_f-one)/(m*(tau_c-one))*(tmp2-tmp1)
					fINxy(i/m,j/m) = ft_f(i,j,k,a,3)
				enddo
			enddo	
			call spline2D(xIN,yIN,fINxy,xOUT,yOUT,fOUTxy)
			
			do j=3,ny_f-2,m
				do i=3,nx_f-2,m
					k=nz_f-1
					s=i-1
					t=j-1
					ft_f(i  ,j  ,k,a,3)=fOUTxy(s  ,t  )
					ft_f(i-1,j  ,k,a,3)=fOUTxy(s-1,t  )
					ft_f(i+1,j  ,k,a,3)=fOUTxy(s+1,t  )
					ft_f(i  ,j-1,k,a,3)=fOUTxy(s  ,t-1)
					ft_f(i  ,j+1,k,a,3)=fOUTxy(s  ,t+1)
				enddo
			enddo	
			
		enddo
	!-----------------------------------------------------------------------
		do k=2,nz_f-1 !copy ft_f(1) to ft_f_
			do j=2,ny_f-1
				do i=2,nx_f-1
					do a=0,q			
						ft_f_(i,j,k,a) = ft_f(i,j,k,a,2) 
					enddo
				enddo
			enddo
		enddo
		
		do iter=1,m
		
			do k=2,nz_f-1 !streaming
				do j=2,ny_f-1
					do i=2,nx_f-1
						do a=0,q	
							ia = i+ex(a)  
							ja = j+ey(a)
							ka = k+ez(a)
							
							if(ja<2) 			ja=ny_f-1
							if(ja>ny_f-1)	ja=2

							f_f(ia,ja,ka,a) = ft_f_(i,j,k,a) 
						enddo
					enddo	
				enddo
			enddo
			!-------------------------------------------
			fx(2) = zer
			fy(2) = zer
			fz(2) = zer
			
			do k=2,nz_f-1 !cylinder BC
				do j=2,ny_f-1
					do i=2,nx_f-1
					
						if (isn(i,k) == 0) then
							do a=0,q
								ia = i+ex(a)  
								ja = j+ey(a) 
								ka = k+ez(a)
								
								if (isn(ia,ka)==1) then					
									f_f(i ,j ,k ,kb(a)) = ft_f_(i ,j, k ,a    ) 
									f_f(ia,ja,ka,a    ) = ft_f_(ia,ja,ka,kb(a)) 
									
									fx(2) = fx(2) + haf*ex(a)*2.0*(-ft_f_(ia,ja,ka,kb(a)) + ft_f_(i,j,k,a)) 						
									fy(2) = fy(2) + haf*ey(a)*2.0*(-ft_f_(ia,ja,ka,kb(a)) + ft_f_(i,j,k,a)) 						
									fz(2) = fz(2) + haf*ez(a)*2.0*(-ft_f_(ia,ja,ka,kb(a)) + ft_f_(i,j,k,a)) 
								endif
							enddo
						endif
						
					enddo	
				enddo
			enddo
			!----------------------------------------------------------------------
			rhoAvg 	= zer
			counter = 0
			
			do k=2,nz_f-1
				do j=2,ny_f-1
					do i=2,nx_f-1
						tmp1 = zer
						tmp2 = zer 
						tmp3 = zer
						tmp4 = zer 
						
						do a=0,q 
							tmp1 = tmp1 + f_f(i,j,k,a) 
							tmp2 = tmp2 + f_f(i,j,k,a)*ex(a) 
							tmp3 = tmp3 + f_f(i,j,k,a)*ey(a) 
							tmp4 = tmp4 + f_f(i,j,k,a)*ez(a)  
						enddo

						rho_f(i,j,k)	= tmp1  
						ux_f (i,j,k)	= tmp2/tmp1  
						uy_f (i,j,k)	= tmp3/tmp1  
						uz_f (i,j,k)	= tmp4/tmp1  

						rhoAvg 	= rhoAvg + tmp1
						counter = counter + 1
					enddo	
				enddo
			enddo
			
			rhoAvg = rhoAvg/counter!((nx_f-2)*(ny_f-2)*(nz_f-2)) 
	!----------------------------------------------------------------------
			do i=3,nx_f-2 !Collision
				do j=2,ny_f-1
					do k=3,nz_f-2
					
						do a=0,q
							tmp1					 = ux_f(i,j,k)*ex(a) + uy_f(i,j,k)*ey(a) + uz_f(i,j,k)*ez(a)  
							tmp2					 = ux_f(i,j,k)**two + uy_f(i,j,k)**two + uz_f(i,j,k)**two 
							feq_f(i,j,k,a) = wt(a)*rho_f(i,j,k)*(1.0d0 + 3.0d0*tmp1 + 4.5d0*tmp1**two - 1.5d0*tmp2) 
							ft_f_(i,j,k,a) = f_f(i,j,k,a) - (f_f(i,j,k,a)-feq_f(i,j,k,a))/tau_f 
						enddo
						
					enddo	
				enddo
			enddo
	!----------------------------------------------------------------------
			fx_t 	= haf*(fx(1)  +  fx(2)) 
			fy_t 	= haf*(fy(1)  +  fy(2)) 
			fz_t 	= haf*(fz(1)  +  fz(2)) 
			
			Cd = fx_t/(haf*rho0*uInlet*uInlet*D*(ny_f-2)) 
			Cs = fy_t/(haf*rho0*uInlet*uInlet*D*(ny_f-2)) 
			Cl = fz_t/(haf*rho0*uInlet*uInlet*D*(ny_f-2)) 
	!--------------------------------------------------------------------		
			if (mod(ts,freqDisp) == 0 .and. iter==m) then
				write(10,'(I8,1X,4(F12.8,1X))') ts, rhoAvg, Cd, Cs, Cl
			endif
	!----------------------------------------------------------------------
			fx(1) = fx(2) 
			fy(1) = fy(2)
			fz(1) = fz(2) 
				
			if(iter==m) exit	
	!==============Temporal integration on fine block boundary==============
			xiL(1)=1.0d0 !(1) is (n-1)th time level
			xiL(2)=2.0d0 !(2) is (n  )th time level
			xiL(3)=3.0d0 !(3) is (n+1)th time level
!~ 			!---------debug-----------------
!~ 			yiL(1)=6!ft_f(i,j,k,a,1) 
!~ 			yiL(2)=9!ft_f(i,j,k,a,2) 
!~ 			yiL(3)=14!ft_f(i,j,k,a,3) 
!~ 			print*,lagrng3p(xiL,yiL,2.5d0)
!~ 			stop
	!-----------------------------west yz-plane-----------------------------
			do k=2,nz_f-1
				do j=2,ny_f-1
					i=2
					do a=0,q
						yiL(1)=ft_f(i,j,k,a,1) 
						yiL(2)=ft_f(i,j,k,a,2) 
						yiL(3)=ft_f(i,j,k,a,3) 
						ft_f_(i,j,k,a) = lagrng3p(xiL,yiL,2.5d0)
					enddo
				enddo
			enddo	
	!-----------------------------east yz-plane-----------------------------	
			do k=2,nz_f-1
				do j=2,ny_f-1
					i=nx_f-1
					do a=0,q
						yiL(1)=ft_f(i,j,k,a,1) 
						yiL(2)=ft_f(i,j,k,a,2) 
						yiL(3)=ft_f(i,j,k,a,3) 
						ft_f_(i,j,k,a) = lagrng3p(xiL,yiL,2.5d0) !two+(double)iter/m 
					enddo
				enddo
			enddo	
!~ 	!----------------------------south xz-plane-----------------------------			
!~ 			do k=2,nz_f-1
!~ 				do i=2,nx_f-1
!~ 					j=2
!~ 					do a=0,q
!~ 						yiL(1)=ft_f(i,j,k,a,1) 
!~ 						yiL(2)=ft_f(i,j,k,a,2) 
!~ 						yiL(3)=ft_f(i,j,k,a,3) 
!~ 						ft_f_(i,j,k,a) = lagrng3p(xiL,yiL,2.5d0)
!~ 					enddo
!~ 				enddo
!~ 			enddo	
!~ 	!----------------------------north xz-plane-----------------------------	
!~ 			do k=2,nz_f-1
!~ 				do i=2,nx_f-1
!~ 					j=ny_f-1
!~ 					do a=0,q
!~ 						yiL(1)=ft_f(i,j,k,a,1) 
!~ 						yiL(2)=ft_f(i,j,k,a,2) 
!~ 						yiL(3)=ft_f(i,j,k,a,3) 
!~ 						ft_f_(i,j,k,a) = lagrng3p(xiL,yiL,2.5d0)
!~ 					enddo
!~ 				enddo
!~ 			enddo			
	!----------------------------bottom xy-plane----------------------------		
			do j=2,ny_f-1
				do i=2,nx_f-1
					k=2
					do a=0,q
						yiL(1)=ft_f(i,j,k,a,1) 
						yiL(2)=ft_f(i,j,k,a,2) 
						yiL(3)=ft_f(i,j,k,a,3) 
						ft_f_(i,j,k,a) = lagrng3p(xiL,yiL,2.5d0) 
					enddo
				enddo
			enddo	
	!------------------------------top xy-plane-----------------------------		
			do j=2,ny_f-1
				do i=2,nx_f-1
					k=nz_f-1
					do a=0,q 
						yiL(1)=ft_f(i,j,k,a,1) 
						yiL(2)=ft_f(i,j,k,a,2) 
						yiL(3)=ft_f(i,j,k,a,3) 
						ft_f_(i,j,k,a) = lagrng3p(xiL,yiL,2.5d0) 
					enddo
				enddo
			enddo	
			
		enddo!iter loop ends
		!---------------------------------------------
		do k=3,nz_f-2
			do j=2,ny_f-1
				do i=3,nx_f-2
					do a=0,q
						ft_f(i,j,k,a,3)	= ft_f_(i,j,k,a) 
					enddo
				enddo	
			enddo
		enddo

	!=========Tranfer post-collision dist on coarse block boundary==========
	!-----------------------------west yz-plane-----------------------------	
		do k=4,nz_f-3,m
			do j=2,ny_f-1,m
				i=4
				do a=0,q
					ft_c(westFace+i/m,1+j/m,bottomFace+k/m,a) = feq_f(i,j,k,a) + m*(tau_c-1)/(tau_f-1)*(ft_f(i,j,k,a,3)-feq_f(i,j,k,a)) 
				enddo
			enddo	
		enddo
	!-----------------------------east yz-plane-----------------------------		
		do k=4,nz_f-3,m
			do j=2,ny_f-1,m
				i=nx_f-3
				do a=0,q
					ft_c(westFace+i/m,1+j/m,bottomFace+k/m,a) = feq_f(i,j,k,a) + m*(tau_c-1)/(tau_f-1)*(ft_f(i,j,k,a,3)-feq_f(i,j,k,a)) 
				enddo
			enddo
		enddo	
!~ 	!-----------------------------south xz-plane-----------------------------	
!~ 		do k=4,nz_f-3,m
!~ 			do i=4,nx_f-3,m
!~ 				j=4
!~ 				do a=0,q
!~ 					ft_c(westFace+i/m,southFace+j/m,bottomFace+k/m,a) = feq_f(i,j,k,a) + m*(tau_c-1)/(tau_f-1)*(ft_f(i,j,k,a,3)-feq_f(i,j,k,a)) 
!~ 				enddo
!~ 			enddo	
!~ 		enddo
!~ 	!-----------------------------north xz-plane----------------------------		
!~ 		do k=4,nz_f-3,m
!~ 			do i=4,nx_f-3,m
!~ 				j=ny_f-3
!~ 				do a=0,q
!~ 					ft_c(westFace+i/m,southFace+j/m,bottomFace+k/m,a) = feq_f(i,j,k,a) + m*(tau_c-1)/(tau_f-1)*(ft_f(i,j,k,a,3)-feq_f(i,j,k,a)) 
!~ 				enddo
!~ 			enddo
!~ 		enddo	
	!-----------------------------bottom xy-plane---------------------------		
		do j=2,ny_f-1,m	
			do i=4,nx_f-3,m
				k=4
				do a=0,q
					ft_c(westFace+i/m,1+j/m,bottomFace+k/m,a) = feq_f(i,j,k,a) + m*(tau_c-1)/(tau_f-1)*(ft_f(i,j,k,a,3)-feq_f(i,j,k,a)) 
				enddo
			enddo
		enddo	
	!------------------------------top xy-plane-----------------------------		
		do j=2,ny_f-1,m	
			do i=4,nx_f-3,m
				k=nz_f-3
				do a=0,q
					ft_c(westFace+i/m,1+j/m,bottomFace+k/m,a) = feq_f(i,j,k,a) + m*(tau_c-1)/(tau_f-1)*(ft_f(i,j,k,a,3)-feq_f(i,j,k,a)) 
				enddo
			enddo
		enddo	
	!-----------------------------------------------------------------------
		do k=2,nz_f-1
			do j=2,ny_f-1
				do i=2,nx_f-1
					do a=0,q
						ft_f(i,j,k,a,1) = ft_f(i,j,k,a,2)
						ft_f(i,j,k,a,2) = ft_f(i,j,k,a,3)
					enddo
				enddo	
			enddo
		enddo
	
	endif !if(cylinderId) ends	
!-----------------------------------------------------------------------
	if(ts == snapList(solNum+1)) then !mod(ts,time/(noOfSnaps-1))
	
		solNum = solNum + 1
!=======================================================================
		write(filename,'(2(A,I2.2),A)') "S",solNum,"_",id+1,".dat"
		open(unit=12,file=filename)
		
		write(12,'(A,I8,A)') 'Title="',ts,'"'
		write(12,'(A)') "Variables=x,y,z,region,u,v,w,rho"
	
		if(cylinderId) then
	!-------------------------------------1---------------------------------
			write(12,'(3(A,I6))') "Zone I= ",en(1)-bn(1)+1,	&
															 " ,J= ",en(2)-bn(2)+1,			&
															 " ,K= ",en(3)-bn(3)+1

			do k=bn(3),en(3) !2,nz-1
				do j=bn(2),en(2) !2,ny-1
					do i=bn(1),en(1)!bn(1),westFace+1
						write(12,'(4(I4,1X),4(F8.4,1X))') i,j,k,0,ux_c(i,j,k)/uInlet,uy_c(i,j,k)/uInlet,uz_c(i,j,k)/uInlet,&
																							(rho_c(i,j,k)-rho0)/(haf*rho0*uInlet**two) 
					enddo
					write(12,*)
				enddo
				write(12,*)
			enddo
!~ 	!-------------------------------------2---------------------------------		
!~ 			write(12,'(3(A,I6))') "Zone I= ",(eastFace-1)-(westFace+1)+1,	&
!~ 															 " ,J= ",en(2)-bn(2)+1,								&
!~ 															 " ,K= ",en(3)-(topFace-1)+1

!~ 			do k=topFace-1,en(3)
!~ 				do j=bn(2),en(2) !2,ny-1
!~ 					do i=westFace+1,eastFace-1
!~ 						write(12,'(4(I4,1X),4(F8.4,1X))')  i,j,k,0,ux_c(i,j,k)/uInlet,uy_c(i,j,k)/uInlet,uz_c(i,j,k)/uInlet,&
!~ 																							(rho_c(i,j,k)-rho0)/(haf*rho0*uInlet**two)  
!~ 					enddo
!~ 					write(12,*)
!~ 				enddo	
!~ 				write(12,*)
!~ 			enddo
!~ 	!-------------------------------------3---------------------------------
!~ 			write(12,'(3(A,I6))')"Zone I= ",en(1)-(eastFace-1)+1,	&
!~ 															" ,J= ",en(2)-bn(2)+1,				&
!~ 															" ,K= ",en(3)-bn(3)+1

!~ 			do k=bn(3),en(3) !2,nz-1
!~ 				do j=bn(2),en(2) !2,ny-1
!~ 					do i=eastFace-1,en(1) !nx-1
!~ 						write(12,'(4(I4,1X),4(F8.4,1X))') i,j,k,0,ux_c(i,j,k)/uInlet,uy_c(i,j,k)/uInlet,uz_c(i,j,k)/uInlet,&
!~ 																							(rho_c(i,j,k)-rho0)/(haf*rho0*uInlet**two) 
!~ 					enddo
!~ 					write(12,*)
!~ 				enddo	
!~ 				write(12,*)
!~ 			enddo
!~ 	!------------------------------------4----------------------------------		
!~ 			write(12,'(3(A,I6))') "Zone I= ",(eastFace-1)-(westFace+1)+1,	&
!~ 															 " ,J= ",en(2)-bn(2)+1,								&
!~ 															 " ,K= ",(bottomFace+1)-bn(3)+1

!~ 			do k=bn(3),bottomFace+1
!~ 				do j=bn(2),en(2) !2,ny-1
!~ 					do i=westFace+1,eastFace-1
!~ 						write(12,'(4(I4,1X),4(F8.4,1X))') i,j,k,0,ux_c(i,j,k)/uInlet,uy_c(i,j,k)/uInlet,uz_c(i,j,k)/uInlet,&
!~ 																							(rho_c(i,j,k)-rho0)/(haf*rho0*uInlet**two) 
!~ 					enddo
!~ 					write(12,*)
!~ 				enddo	
!~ 				write(12,*)
!~ 			enddo
!~ 	!-------------------------------------5----------------------------------		
!~ 			write(12,'(3(A,I6))') "Zone I= ",(eastFace-1)-(westFace+1)+1,	&
!~ 															 " ,J= ",(southFace+1)-bn(2)+1,						&
!~ 															 " ,K= ",(topFace-1)-(bottomFace+1)+1

!~ 			do k=bottomFace+1,topface-1
!~ 				do j=bn(2),southFace+1
!~ 					do i=westFace+1,eastFace-1
!~ 						write(12,'(4(I4,1X),4(F8.4,1X))') i,j,k,0,ux_c(i,j,k)/uInlet,uy_c(i,j,k)/uInlet,uz_c(i,j,k)/uInlet,&
!~ 																							(rho_c(i,j,k)-rho0)/(haf*rho0*uInlet**two) 
!~ 					enddo
!~ 					write(12,*)
!~ 				enddo	
!~ 				write(12,*)
!~ 			enddo
!~ 	!-------------------------------------6---------------------------------		
!~ 			write(12,'(3(A,I6))') "Zone I= ",(eastFace-1)-(westFace+1)+1,	&
!~ 															 " ,J= ",en(2)-(northFace-1)+1,				&
!~ 															 " ,K= ",(topFace-1)-(bottomFace+1)+1

!~ 			do k=bottomFace+1,topFace-1
!~ 				do j=northFace-1,en(2)
!~ 					do i=westFace+1,eastFace-1
!~ 						write(12,'(4(I4,1X),4(F8.4,1X))') i,j,k,0,ux_c(i,j,k)/uInlet,uy_c(i,j,k)/uInlet,uz_c(i,j,k)/uInlet,&
!~ 																							(rho_c(i,j,k)-rho0)/(haf*rho0*uInlet**two) 
!~ 					enddo
!~ 					write(12,*)
!~ 				enddo	
!~ 				write(12,*)
!~ 			enddo				
!-----------------------------------2-----------------------------------	
			write(12,'(3(A,I6))') "Zone I= ",nx_f-2,	&
															 " ,J= ",ny_f-2,	&
															 " ,K= ",nz_f-2

			do k=2,nz_f-1
				do j=2,ny_f-1
					do i=2,nx_f-1
						write(12,'(3(F8.2,1X),I4,1X,4(F8.4,1X))') (westFace+dble(i)/m),(1+dble(j)/m),(bottomFace+dble(k)/m),&
																											isn(i,k),ux_f(i,j,k)/uInlet,uy_f(i,j,k)/uInlet,uz_f(i,j,k)/uInlet,&
																											(rho_f(i,j,k)-rho0)/(haf*rho0*uInlet**two) 
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
						write(12,'(4(I4,1X),4(F8.4,1X))') i,j,k,0,ux_c(i,j,k)/uInlet,uy_c(i,j,k)/uInlet,uz_c(i,j,k)/uInlet,&
																							(rho_c(i,j,k)-rho0)/(haf*rho0*uInlet**two) 
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
			write(12,'(3(i4))') nP(1),nP(2),nP(3)
			write(12,'(i3)') nproc
			do i=0,nproc-1
				write(12,'(6(i4,1x))') bb(1,i),ee(1,i),bb(2,i),ee(2,i),bb(3,i),ee(3,i)
			enddo
			close(12)	
		endif
!=======================================================================		
		if(masterId) write(*,'(A,I8)') "Solution recorded at time ",ts
 	endif
!~ 		!-------Patch fine block boundary for display------------
!~ 		!left fine block boundary
!~ 		for(i=1,j=1; j<=ny_f; j+=m)
!~ 		{
!~ 			ii=(i-1)/m;
!~ 			jj=(j-1)/m;
!~ 			ux_f(i,j) = ux_c(westFace+ii,bottomFace+jj);
!~ 			xi(jj)=j;
!~ 			yi(jj)=ux_f(i,j);
!~ 			if(j==ny_f) break;
!~ 			xo(jj)=j+1;!hardwired
!~ 		}
		
!~ 		spline2d(xi,yi,xo,yo,SIZEOF(xi),SIZEOF(xo));
!~ 		for(i=1,j=1; j<ny_f; j+=m)
!~ 		{
!~ 			ux_f(i,j+1)=yo((j-1)/2);!hardwired
!~ 		}	
!~ 		!right fine block boundary
!~ 		for(i=nx_f,j=1; j<=ny_f; j+=m)
!~ 		{
!~ 			ii=(i-1)/m;
!~ 			jj=(j-1)/m;
!~ 			ux_f(i,j) = ux_c(westFace+ii,bottomFace+jj);
!~ 			xi(jj)=j;
!~ 			yi(jj)=ux_f(i,j);
!~ 			if(j==ny_f) break;
!~ 			xo(jj)=j+1;!hardwired
!~ 		}
!~ 		spline2d(xi,yi,xo,yo,SIZEOF(xi),SIZEOF(xo));
		
!~ 		for(i=nx_f,j=1; j<ny_f; j+=m)
!~ 		{
!~ 			ux_f(i,j+1)=yo((j-1)/2);
!~ 		}
!~ 		!bottomFace fine block boundary
!~ 		for(j=1,i=1; i<=nx_f; i+=m)
!~ 		{
!~ 			ii=(i-1)/m;
!~ 			jj=(j-1)/m;
!~ 			ux_f(i,j) = ux_c(westFace+ii,bottomFace+jj);
!~ 			xi(ii)=i;
!~ 			yi(ii)=ux_f(i,j);
!~ 			if(i==nx_f) break;
!~ 			xo(ii)=i+1;!hardwired
!~ 		}
!~ 		spline2d(xi,yi,xo,yo,SIZEOF(xi),SIZEOF(xo));
		
!~ 		for(j=1,i=1; i<nx_f; i+=m)
!~ 		{
!~ 			ux_f(i+1,j)=yo((i-1)/2);
!~ 		}		
!~ 		!topFace fine block boundary
!~ 		for(j=ny_f,i=1; i<=nx_f; i+=m)
!~ 		{
!~ 			ii=(i-1)/m;
!~ 			jj=(j-1)/m;
!~ 			ux_f(i,j) = ux_c(westFace+ii,bottomFace+jj);
!~ 			xi(ii)=i;
!~ 			yi(ii)=ux_f(i,j);
!~ 			if(i==nx_f) break;
!~ 			xo(ii)=i+1;!hardwired
!~ 		}
!~ 		spline2d(xi,yi,xo,yo,SIZEOF(xi),SIZEOF(xo));
		
!~ 		for(j=ny_f,i=1; i<nx_f; i+=m)
!~ 		{
!~ 			ux_f(i+1,j)=yo((i-1)/2);
!~ 		}
!~ 		!-----------------	
!=======================================================================
!~ 		write(filename,'(A,I2,A)') "fine",solNum,".dat"
!~ 		open(unit=12,file=filename)
		
!~ 		write(12,*) "Variables=x,y,region,u,v,rho"
!~ 		write(12,*) "Zone I=",nx_f,"J=",ny_f 

!~ 		do j=2,ny_f+1
!~ 			do i=2,nx_f+1
!~ 				write(12,'(2(F8.2,1X),I4,1X,3(F8.4,1X))')  (westFace+(double)(i)/m),(bottomFace+(double)(j)/m),isn(i,j),ux_f(i,j),uy_f(i,j),rho_f(i,j) 
!~ 			enddo
!~ 			write(*,*)
!~ 		enddo
		
!~ 		close(12)
		
!~ 		write(*,*) "Fine recorded at time ",ts
!----------------------------------------------------------------------
	if (cylinderId .and. mod(ts,freqDisp) == 0) then
		call cpu_time(tEnd)
		write(*,'(I3,I8,2X,A,F8.2,A,I4,A)') id, ts," | CPU Time: ",tEnd-tStart," sec for ",freqDisp," timesteps"
	endif	

enddo!Time loop Ends

if(cylinderId) close(10) 

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
