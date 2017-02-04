!*************************************空间刚架程序设计***********************************************************************************
!******************************主程序部分************************************************************************************************
program main
!je为储存单元两端节点编号数组，je(3,ne)为单元虚节点编号,nod为节点自由度编号数组
integer,dimension(:,:),allocatable::je,nod                
!cfd为约束自由度数组,efnl为非节点载荷单元储存数组，p为天空线一维定位数组
integer,dimension(:),allocatable::cfd,efnl,p            
!cx,cy,cz为节点坐标数组，a,e为单元横截面积和弹性模量，fq为外载荷数组
!cfnl为非节点载荷数组,eiy,eiz为单元极惯性矩,eg为剪切模量，ej为截面极惯性矩,sky为天空线刚度矩阵一维储存数组
real,dimension(:),allocatable::cx,cy,cz,a,f,fq,e,cfnl,eiy,eiz,eg,ej,sky                
!ke为单刚，lg为单元定位数组,fne为单元内力数组,t为转换矩阵，t1为转换矩阵的转置,at为矢量变换矩阵
dimension ke(12,12),lg(12),fne(12,1),t(12,12),t1(12,12),at(3,3)  
!cos,sin为单元定位向量，l为单元长度      
real ke,cos,sin,l,t,t1,at 
!nn为节点数，ne为单元数，nd为总自由度数（原始自由度数，不考虑约束情况），nfd为被约束的自由度数,nf为非节点载荷数 
!nfq为受节点力的节点数，ii为与总体自由度i的所以相关自由度的最小值，nwk为天空线数组长度                              
integer nn,ne,nd,nfd,lg,nf,nfq,ii,nwk
                    
open(5,file='shuju.txt')

read(5,*) nn,ne,nd,nfd,nf,nfq               !************            
open(6,file='result.txt')                               !*****输入输出节点数，单元数，原始总自由度数，被约束的自由度数
write(6,7) nn,ne,nd,nfd,nf,nfq                          !*****有非节点载的荷单元数，受节点力的节点数
7 format(2x,'nn  ne  nd  nfd  nf  nfq'/6i4) !************

!对动态数组赋予确定范围宽度     
allocate(cx(nn+ne),cy(nn+ne),cz(nn+ne),je(3,ne),nod(6,nn),a(ne),e(ne),eiy(ne),eiz(ne),eg(ne),&
ej(ne),cfd(nfd),p(nd-nfd),f(nd-nfd),fq(nd),cfnl(12*nf),efnl(nf))

write(6,*)           !输出一行空格

read(5,*) (cx(i),cy(i),cz(i),i=1,nn)       !***************
write(6,8)                                                !
8 format(3x,'NODE   X-COORD   Y-COORD   Z-COORD')         !*****输入输出节点坐标
write(6,9) (i,cx(i),cy(i),cz(i),i=1,nn)                   !
9 format(2x,i4,1x,3f10.5)                  !*************** 
 
write(6,*)           !输出一行空格

!******************输入输出单元节点编号，弹性模量，横截面积，两个极惯性矩，剪切模量，横截面极惯性矩**************** 
read(5,*) ((je(i,j),i=1,3),e(j),a(j),eiy(j),eiz(j),eg(j),ej(j),j=1,ne)      
write(6,10)                                                             !
10 format(3x,'ELEMENT   NODE1   NODE2   NODEV        E          AREA          EIX         EIZ         EG          EJ')  
write(6,11) (j,(je(i,j),i=1,3),e(j),a(j),eiy(j),eiz(j),eg(j),ej(j),j=1,ne)                  
11 format(3x,i4,3x,3i7,2x,6e13.5) 

write(6,*)           !输出一行空格

read(5,*) (cx(i),cy(i),cz(i),i=nn+1,nn+ne)!************
write(6,25)                                           !
25 format(3x,'V-NODE   X-COORD   Y-COORD   Z-COORD')  !*****输入输出各单元虚节点坐标
write(6,26) (i,cx(i),cy(i),cz(i),i=nn+1,nn+ne)        !
26 format(4x,i4,1x,3f10.5)                !************  
          
write(6,*)           !输出一行空格

read(5,*) (cfd(i),i=1,nfd)   !********
write(6,12)                          !
12 format(5x,'约束自由度号:')        !*****输入输出被约束的原始自由度编号
write(6,13) (cfd(i),i=1,nfd)         !
13 format(2x,10i3)           !********

write(6,*)           !输出一行空格

fq=0.0       !节点外载先赋零
write(6,*) '                       单元节点载荷'               !*********************************
read(5,*) (j,(fq(6*j+k-6),k=1,6),i=1,nfq)      !j为节点编号                                     !
write(6,14)                                                                                     !   
14 format(3x,'NODE      X-LOAD      Y-LOAD       Z-LOAD         MX           MY          MZ')   !******输入输出节点载荷
write(6,15) (i,(fq(6*i+j-6),j=1,6),i=1,nn)                                                      !
15 format(3x,i3,1x,6e13.5)                                     !*********************************

write(6,*)           !输出一行空格

!**********************************输入输出单元非节点载荷**************************************************
write(6,*) '                       单元非节点载荷'
read(5,*) (efnl(i),(cfnl(12*i-12+j),j=1,12),i=1,nf)   
write(6,16)                                                                                            
16 format(3x,'ELEMENT             RX           RY           RZ           MX          MY           MZ')            
write(6,17) (efnl(i),(cfnl(12*i-12+j),j=1,12),i=1,nf)     
17 format(3x,i4,'(LEFT)',3x,6e13.5/6x,'(RIGHT)',3x,6e13.5) 

!******计算各节点的新自由度编号（被约束住的自由度赋零）**********
call cnod(nod,cfd,nfd,nn)


!**********计算总刚矩阵处于上三角的每一列最小非零元素的行号和每一列在天空线矩阵中的元素个数
do i=1,nd-nfd
 ii=nd-nfd
 do j=1,ne
 call clg(lg,nod,je,j) 
  do m=1,12
  i1=lg(m)
   if(i1.eq.i) then
   call cbh(lg,i,ii)
   endif
  enddo
 enddo
p(i)=i-ii+1
enddo


!***********计算天空线一维定位数组并确定天空线数组长度************
p(1)=1
do i=2,nd-nfd
p(i)=p(i-1)+p(i)
enddo
nwk=p(nd-nfd)
allocate(sky(nwk))
   
                                                                      
!************计算总刚的天空线存储数组**************
sky=0.0

do i=1,ne 
ke=0.0                                                       

call clg(lg,nod,je,i)     !调用子程序计算lg                                                                           
call cat(at,l,i,nn,je,cx,cy,cz)      !调用子程序cat计算单元矢量变换矩阵at和长度l                             
call cdg(ke,l,e(i),a(i),eiy(i),eiz(i),eg(i),ej(i))     !调用子程序cdg计算局部坐标系下的单刚   
call caut(at,t,t1)             !调用子程序caut计算转换矩阵t及其转置t1                         
ke=matmul(t1,ke)                                                                              
ke=matmul(ke,t)                          !计算整体坐标系下的单刚                              
call caq(efnl,cfnl,t1,fq,i,je,nf,nd,ne) !调用子程序caq计算修正后的节点外载荷数组                                                     
 !***********确定sky*************                                               
  do m=1,12
   if(lg(m).gt.0) then
   j=p(lg(m))
   sky(j)=sky(j)+ke(m,m)                                                                                  
    do n=1,12
	 if(lg(n).gt.0) then                                                                                 
      if(lg(n).lt.lg(m)) then
	  i1=p(lg(m))+lg(n)-lg(m) 
	  sky(i1)=sky(i1)+ke(m,n)
	  endif 
	 endif                                                     
    enddo
   endif	                                                                                       
  enddo  
                                                                                  
enddo                                                           


!***********调用子程序计算修正的外载f（去除零自由度所对应的外载）*******
f=0.0
call cf(f,fq,cfd,nfd,nd)

!*********调用子程序解线性方程组求解位移*********
call skylineLDL(sky,p,f,nwk,nd-nfd)

!**********调用子程序计算修正后的各节点位移*******
fq=0.0
call cfq(f,fq,cfd,nfd,nd)


write(6,*)           !输出一行空格
 
write(6,*) '                   节点位移与转角'        !*************************************************
write(6,18)                                                                                      !
18 format(3x,'NODE      X-DISP       Y-DISP      Z-DISP        ROT-X        ROT-Y       ROT-Z')  !*****输出节点位移和转角
write(6,19) (i,(fq(6*i-6+j),j=1,6),i=1,nn)                                                       !
19 format(3x,i3,1x,6e13.5)                      !*************************************************


  
!********************************************计算并输出单元内力******************************************************  
write(6,*) '                       单元内力'                    
write(6,22)                                                
22 format(3x,'ELEMENT              FX           FY           FZ           MX           MY          MZ')                      
do i=1,ne 
ke=0.0                                                                                           
call cat(at,l,i,nn,je,cx,cy,cz)      !调用子程序cat计算单元矢量变换矩阵at和长度l                                                                
call caut(at,t,t1)                          !调用子程序caut计算转换矩阵t及其转置t1                                                    
call cdg(ke,l,e(i),a(i),eiy(i),eiz(i),eg(i),ej(i))       !调用子程序cdg计算局部坐标系下的单刚                                                                                  
call caf(fne,ke,t,fq,je,i,cfnl,efnl,nf,nd) !调用caf子程序计算单元内力                                
write(6,23) i,(fne(j,1),j=1,12)                                                                       
23 format(3x,i4,'(LEFT)',3x,6e13.5/6x,'(RIGHT)',3x,6e13.5)                                                                           
enddo                                                                                                              
end  
                                                     
!*********************************子程序部分*************************************************************************************
!子程序1：计算矢量变换矩阵
subroutine  cat(at,l,ie,nn,je,x,y,z)
dimension at(3,3),je(3,1),x(1),y(1),z(1)
real at,x,y,z,a1,a2,a3,l,b1,b2,b3,s1,s2,s3,u1,u2,u3,uu
integer je,nn
i=je(1,ie)
j=je(2,ie)
k=je(3,ie)
a1=x(j)-x(i)
a2=y(j)-y(i)
a3=z(j)-z(i)
l=sqrt(a1*a1+a2*a2+a3*a3)
at(1,1)=a1/l
at(1,2)=a2/l
at(1,3)=a3/l
b1=x(k)-x(i)
b2=y(k)-y(i)
b3=z(k)-z(i)
s1=a2*b3-a3*b2
s2=a3*b1-a1*b3
s3=a1*b2-a2*b1
u1=s2*a3-s3*a2
u2=s3*a1-s1*a3
u3=s1*a2-s2*a1
uu=sqrt(u1*u1+u2*u2+u3*u3)
at(2,1)=u1/uu
at(2,2)=u2/uu
at(2,3)=u3/uu
at(3,1)=at(1,2)*at(2,3)-at(1,3)*at(2,2)
at(3,2)=at(1,3)*at(2,1)-at(1,1)*at(2,3)
at(3,3)=at(1,1)*at(2,2)-at(1,2)*at(2,1)
end
!子程序2：计算单元转换矩阵和转换矩阵的转置
subroutine caut(at,t,t1)
dimension t(12,12),t1(12,12),at(3,3)
real t,t1,at
t=0.0
t1=0.0
do i=1,3
 do j=1,3
 t(i,j)=at(i,j)
 t(i+3,j+3)=at(i,j)
 t(i+6,j+6)=at(i,j)
 t(i+9,j+9)=at(i,j)
 enddo
enddo
do i=1,12
 do j=1,12
 t1(i,j)=t(j,i)
 enddo
enddo
end
!子程序3：局部坐标系下单刚的计算
subroutine cdg(ke,l,e,a,eiy,eiz,eg,ej)
dimension ke(12,12)
real l,e,a,eiy,eiz,eg,ej,ke
ke=0.0
ke(1,1)=e*a/l
ke(1,7)=-ke(1,1)
ke(2,2)=12*e*eiz/(l**3)
ke(2,6)=6*e*eiz/(l**2)
ke(2,8)=-ke(2,2)
ke(2,12)=ke(2,6)
ke(3,3)=12*e*eiy/(l**3)
ke(3,5)=-6*e*eiy/(l**2)
ke(3,9)=-ke(3,3)
ke(3,11)=ke(3,5)
ke(4,4)=eg*ej/l
ke(4,10)=-ke(4,4)
ke(5,5)=4*e*eiy/l
ke(5,9)=6*e*eiy/(l**2)
ke(5,11)=2*e*eiy/l
ke(6,6)=4*e*eiz/l
ke(6,8)=-6*e*eiz/(l**2)
ke(6,12)=2*e*eiz/l
ke(7,7)=e*a/l
ke(8,8)=12*e*eiz/(l**3)
ke(8,12)=-6*e*eiz/(l**2)
ke(9,9)=12*e*eiy/(l**3)
ke(9,11)=6*e*eiy/(l**2)
ke(10,10)=eg*ej/l
ke(11,11)=4*e*eiy/l
ke(12,12)=4*e*eiz/l
do i=1,12
 do j=1,12
 ke(j,i)=ke(i,j)
 enddo
enddo
end
!子程序4：将非节点载荷等效为节点载荷并求出修正后的外载荷
subroutine caq(efnl,cfnl,t1,fq,ie,je,nf,nd,ne)
dimension efnl(nf),cfnl(12*nf),t1(12,12),q(12,1),fq(nd),je(3,ne)
integer nf,nd,efnl,je
real q,cfnl,t1,fq
do i=1,nf
k=efnl(i)
 if(ie.eq.k) then
  do j=1,12
  q(j,1)=-cfnl(12*i-12+j)
  enddo
q=matmul(t1,q)
i2=je(1,k)
i3=je(2,k)
  do m=1,6
  fq(6*i2+m-6)=fq(6*i2+m-6)+q(m,1)
  fq(6*i3+m-6)=fq(6*i3+m-6)+q(m+6,1)
  enddo
 endif
enddo
end


!子程序5：计算节点自由度编号
subroutine cnod(nod,cfd,nfd,nn)
dimension nod(6,nn),cfd(nfd),a(6)
integer nod,cfd,nfd,nn,a,nuf
nuf=0
do i=1,nn
a=1
 do j=1,nfd
 j1=cfd(j)
  do k=1,6
   if((6*i+k-6).eq.j1) then
   nod(k,i)=0
   a(k)=0
   endif
  enddo
 enddo
 do m=1,6
  if(a(m)) then
  nuf=nuf+1
  nod(m,i)=nuf
  endif
 enddo
enddo  
end


!子程序6：计算单元节点自由度与总体自由度的关联数组
subroutine clg(lg,nod,je,ie)
dimension lg(12),nod(6,1),je(3,1)
integer lg,nod,je,ie
do i=1,2
 i1=je(i,ie)
 do j=1,6
 lg(6*i+j-6)=nod(j,i1)
 enddo
enddo
end

!子程序7：计算总刚每一列的最小非零元素的行编号
subroutine cbh(lg,ie,ii)
dimension lg(12)
integer lg,ie,ii
do i=1,12
 if(lg(i).gt.0) then
  if(lg(i).lt.ii) then
  ii=lg(i)
  endif
 endif
enddo
end

!子程序8：在已知节点力数组的情况下求解解线性方程组所需的外载数组
subroutine cf(f,fq,cfd,nfd,nd)
dimension cfd(nfd),fq(nd),f(nd-nfd)
integer cfd,nfd,nd,a,p,d
real f,fq
a=0
p=0
if(cfd(1)) then
d=2
a=1
else
d=1
endif
do i=d,nd
    p=0
	do i1=1,nfd
	j1=cfd(i1)
		if(i.eq.j1) then
		p=1
		endif
    enddo
	if(p) then
	a=a+1
	else
	f(i-a)=fq(i)
	endif
enddo
end


!子程序9：在已知线性方程组解的情况下求节点位移
subroutine cfq(f,fq,cfd,nfd,nd)
dimension cfd(nfd),fq(nd),f(nd-nfd)
integer cfd,nfd,nd,a,p,d
real f,fq
a=0
p=0
if(cfd(1)) then
d=2
a=1
else
d=1
endif
do i=d,nd
    p=0
	do i1=1,nfd
	j1=cfd(i1)
		if(i.eq.j1) then
		p=1
		endif
    enddo
	if(p) then
	a=a+1
	else
	fq(i)=f(i-a)
	endif
enddo
end

!子程序10：计算单元节点内力
subroutine caf(fne,ke,t,fq,je,ie,cfnl,efnl,nf,nd)
dimension fne(12,1),ke(12,12),t(12,12),fq(nd),q(12,1),je(3,1),efnl(nf),cfnl(12*nf),a(12,12)
real fne,ke,t,fq,q,cfnl,a
integer ie,je,nf,nd,efnl
i1=je(1,ie)
i2=je(2,ie)
do i=1,6
q(i,1)=fq(i1*6-6+i)
q(i+6,1)=fq(i2*6-6+i)
enddo
a=matmul(ke,t)
fne=matmul(a,q)
do i=1,nf
j=efnl(i)
 if(j.eq.ie) then
  do k=1,12
  fne(k,1)=fne(k,1)+cfnl(12*i-12+k)
  enddo
 endif
enddo
end

!子程序11：天空线的LDL分解法（解线性方程组）
subroutine skylineLDL(sky,P,R,Nsky,NP)
real sky(Nsky),R(NP)
integer P(NP),nr(NP)
real sigma,term

nr(1)=1
do,i=2,NP
nr(i)=i-(P(i)-P(i-1))+1
enddo
!************三角分解，求L(ij)和d(ij)*******
do,j=2,NP
!==================求g(ij)================== 
 do,i=nr(j),j-1
 sigma=0
  nrmax=max(nr(i),nr(j))
   do,k=nrmax,i-1
   term=sky(P(i)-(i-k))*sky(P(j)-(j-k))
   sigma=sigma+term
   enddo
 sky(P(j)-(j-i))=sky(P(j)-(j-i))-sigma
 enddo
!=================求d(jj)===================
 sigma=0
  do,k=nr(j),j-1
  term=(sky(P(j)-(j-k)))**2/sky(P(k))
  sigma=sigma+term
  enddo
 sky(P(j))=sky(P(j))-sigma
!==================求L(ij)================== 
 do,i=nr(j),j-1
 sky(P(j)-(j-i))=sky(P(j)-(j-i))/sky(P(i))
 enddo
enddo
!************前消，求{V}**********
do,j=2,NP
 sigma=0
  do,i=nr(j),j-1
  term=sky(P(j)-(j-i))*R(i)
  sigma=sigma+term
  enddo
 R(j)=R(j)-sigma
enddo
!************回代，求{U}***********
!求{V-bar}的值
do,i=1,NP
 R(i)=R(i)/sky(P(i))
enddo
!求{U}的值
do,i=NP,2,-1
 do,j=nr(i),i-1
 R(j)=R(j)-sky(P(i)-(i-j))*R(i)
 enddo
enddo
end