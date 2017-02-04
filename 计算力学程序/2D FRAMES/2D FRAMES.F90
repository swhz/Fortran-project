!********************************************平面刚架程序设计***************************************************************
!******************************主程序部分************************************************************************************************
program main
!定义动态数组
integer,dimension(:,:),allocatable::je,nod                 !je为储存单元两端节点编号数组，nod为节点自由度编号数组

integer,dimension(:),allocatable::cfd,efnl,p               !cfd为约束自由度数组,efnl为非节点载荷单元储存数组，p为天空线一维定位数组

real,dimension(:),allocatable::cx,cy,a,f,fq,e,cfnl,ei,sky  !cx,cy为节点坐标数组，a,e为单元横截面积和弹性模量，f为外载荷数组,
                                                           !fr为计算所得节点力数组，cfnl为非节点载荷数组,ei为单元极惯性矩
														   !sky为天空线刚度矩阵一维储存数组

dimension ke(6,6),lg(6),fne(6,1),t(6,6),t1(6,6)            !ke为单刚，lg为单元定位数组,fne为单元内力数组,t为转换矩阵，t1为转换矩阵的转置

real ke,cos,sin,l,t,t1,fne                                    !cos,sin为单元定位向量，l为单元长度

integer nn,ne,nd,nfd,lg,nf,nfq,ii,nwk           !nn为节点数，ne为单元数，nd为原始总自由度数（不考虑约束），nfd为约束自由度数,nf为非节点载荷数
                                                !nfq为受节点力的节点数，ii为与总体自由度i的所以相关自由度的最小值，nwk为天空线数组长度

open(5,file='shuju.txt')

read(5,*) nn,ne,nd,nfd,nf,nfq                !**************             
open(6,file='result.txt')                                  !*****输入输出节点数，单元数，原始总自由度数
write(6,7) nn,ne,nd,nfd,nf,nfq                             !*****被约束的自由度数，非节点载荷数，受节点力的节点数
7 format(2x,'nn  ne  nd  nfd  nf  nfq'/6i4)  !************** 

!给动态数组赋予确定范围边界   
allocate(cx(nn),cy(nn),je(2,ne),nod(3,nn),a(ne),e(ne),ei(ne),cfd(nfd),fq(nd),f(nd-nfd),p(nd-nfd),cfnl(6*nf),efnl(nf))

write(6,*)     !输出一行空格

read(5,*) (cx(i),cy(i),i=1,nn)       !***************
write(6,8)                                          !
8 format(3x,'NODE   X-COORD   Y-COORD')             !*****输入输出节点坐标
write(6,9) (i,cx(i),cy(i),i=1,nn)                   !
9 format(2x,i4,1x,2f10.5)            !***************
  
write(6,*)      !输出一行空格
    
read(5,*) ((je(i,j),i=1,2),e(j),a(j),ei(j),j=1,ne)         !****************
write(6,10)                                                                 !
10 format(3x,'ELEMENT  NODE1  NODE2         E          AREA         ei')    !*****输入输出单元节点编号，弹性模量，横截面积，极惯性矩
write(6,11) (j,(je(i,j),i=1,2),e(j),a(j),ei(j),j=1,ne)                      !
11 format(3x,i4,2i8,1x,3e13.5)                              !****************

write(6,*)       !输出一行空格

read(5,*) (cfd(i),i=1,nfd)   !********
write(6,12)                          !
12 format(5x,'约束自由度号:')        !*****输入输出被约束的原始自由度编号
write(6,13) (cfd(i),i=1,nfd)         !
13 format(2x,10i3)           !******** 

write(6,*)       !输出一行空格

fq=0.0                                        !************************
write(6,*) '                       单元节点载荷'                      !
read(5,*) (j,fq(3*j-2),fq(3*j-1),fq(3*j),i=1,nfq) !j为节点编号        !
write(6,14)                                                           !
14 format(3x,'NODE      X-LOAD      Y-LOAD       MOMENT')             !******输入输出节点载荷
write(6,15) (i,fq(3*i-2),fq(3*i-1),fq(3*i),i=1,nn)                    !
15 format(3x,i3,1x,3e13.5)                    !************************

write(6,*)       !输出一行空格    

write(6,*) '                       单元非节点载荷'
read(5,*) (efnl(i),cfnl(6*i-5),cfnl(6*i-4),cfnl(6*i-3),cfnl(6*i-2),cfnl(6*i-1),cfnl(6*i),i=1,nf)    !*********
write(6,16)                                                                                                  !
16 format(3x,'ELEMENT        RXI          RYI          MI          RXJ          RYJ          MJ')            !*****输入输出单元非节点载荷
write(6,17) (efnl(i),cfnl(6*i-5),cfnl(6*i-4),cfnl(6*i-3),cfnl(6*i-2),cfnl(6*i-1),cfnl(6*i),i=1,nf)           !
17 format(3x,i4,3x,6e13.5)                                                                          !*********

!******计算各节点的新自由度编号（被约束住的自由度赋零）**********
call cnod(nod,cfd,nfd,nn)


!**********计算总刚矩阵处于上三角的每一列最小非零元素的行号和每一列在天空线矩阵中的元素个数
do i=1,nd-nfd
 ii=nd-nfd
 do j=1,ne
 call clg(lg,nod,je,j) 
  do m=1,6
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
call elcs(cos,sin,l,i,je,cx,cy)      !调用子程序elcs计算单元cos,sin,l                           
call cdg(ke,l,e(i),a(i),ei(i))     !调用子程序cdg计算局部坐标系下的单刚          
call caut(cos,sin,t,t1)             !调用子程序caut计算转换矩阵t及其转置t1      
ke=matmul(t1,ke)                         !计算整体坐                            
ke=matmul(ke,t)                          !标系下的单刚                          
call caq(efnl,cfnl,t1,fq,i,je,nf,nd,ne)!调用子程序caq计算修正后的节点外载荷数组                                                              
 !***********确定sky*************                                               
  do m=1,6
   if(lg(m).gt.0) then
   j=p(lg(m))
   sky(j)=sky(j)+ke(m,m)                                                                                  
    do n=1,6
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

write(6,*)       !输出一行空格

write(6,*) '                   节点位移'    !**************
write(6,18)                                               !
18 format(3x,'NODE      X-DISP      Y-DISP        ROT')   !*****输出节点位移和转角
write(6,19) (i,fq(3*i-2),fq(3*i-1),fq(3*i),i=1,nn)        !
19 format(3x,i3,1x,3e13.5)                  !**************

                 
write(6,*)      !输出一行空格
  
write(6,*) '                       单元内力'                    
write(6,22)                                                !******************************************
22 format(3x,'ELEMENT        FXI          FYI          MZI          FXJ          FYJ          MZJ')  !                      
do i=1,ne                                                                                            !
call elcs(cos,sin,l,i,je,cx,cy)                                                                      !
call caut(cos,sin,t,t1)                                                                              !
call cdg(ke,l,e(i),a(i),ei(i))                                                                       !*****计算并输出单元内力                     
call caf(fne,ke,t,fq,je,i,cfnl,efnl,nf,nd) !调用caf子程序计算单元内力                                !
write(6,23) i,(fne(j,1),j=1,6)                                                                       !
23 format(3x,i4,3x,6e13.5)                                                                           !
enddo                                                                                                !                 
end                                                        !****************************************** 
!*********************************子程序部分*************************************************************************************
!子程序1：计算单元转换矩阵和转换矩阵的转置
subroutine caut(cos,sin,t,t1)
dimension t(6,6),t1(6,6)
real cos,sin,t,t1
t=0.0
t1=0.0
t(1,1)=cos
t(1,2)=sin
t(2,1)=-sin
t(2,2)=cos
t(3,3)=1.0
do i=4,6
 do j=4,6
 t(i,j)=t(i-3,j-3)
 enddo
enddo
do i=1,6
 do j=1,6
 t1(i,j)=t(j,i)
 enddo
enddo
end
!子程序2：计算单元长度l，单元cos,sin
subroutine elcs(cos,sin,l,ie,je,x,y)
dimension je(2,1),x(1),y(1)
real cos,sin,l,s1,s2,x,y
integer je
i1=je(1,ie)
j1=je(2,ie)
s1=x(j1)-x(i1)
s2=y(j1)-y(i1)
l=sqrt(s1*s1+s2*s2)
cos=s1/l
sin=s2/l 
end
!子程序3：局部坐标系下单刚的计算
subroutine cdg(ke,l,e,a,ei)
dimension ke(6,6)
real l,e,a,ei,ke
ke=0.0
ke(1,1)=e*a/l
ke(1,4)=-e*a/l
ke(2,2)=12*e*ei/(l**3)
ke(2,3)=6*e*ei/(l**2)
ke(2,5)=-ke(2,2)
ke(2,6)=ke(2,3)
ke(3,3)=4*e*ei/l
ke(3,5)=-ke(2,3)
ke(3,6)=2*e*ei/l
ke(4,4)=e*a/l
ke(5,5)=ke(2,2)
ke(5,6)=-ke(2,3)
ke(6,6)=ke(3,3)
do i=1,6
 do j=1,6
 ke(j,i)=ke(i,j)
 enddo
enddo
end
!子程序4：将非节点载荷等效为节点载荷并求出修正后的外载荷
subroutine caq(efnl,cfnl,t1,fq,ie,je,nf,nd,ne)
dimension efnl(nf),cfnl(6*nf),t1(6,6),q(6,1),fq(nd),je(2,ne)
integer nf,nd,efnl,je
real q,cfnl,t1,fq
do i=1,nf
k=efnl(i)
 if(ie.eq.k) then
  do j=1,6
  q(j,1)=-cfnl(6*i-6+j)
  enddo
q=matmul(t1,q)
i1=efnl(i)
i2=je(1,i1)
i3=je(2,i1)
  do k=1,3
  fq(3*i2+k-3)=fq(3*i2+k-3)+q(k,1)
  fq(3*i3+k-3)=fq(3*i3+k-3)+q(k+3,1)
  enddo
 endif
enddo
end
!子程序5：计算单元节点内力
subroutine caf(fne,ke,t,fq,je,ie,cfnl,efnl,nf,nd)
dimension fne(6,1),ke(6,6),t(6,6),fq(nd),q(6,1),je(2,1),efnl(nf),cfnl(6*nf),a(6,6)
real fne,ke,t,fq,q,cfnl,a
integer ie,je,nf,nd,efnl
i1=je(1,ie)
i2=je(2,ie)
do i=1,3
q(i,1)=fq(i1*3-3+i)
q(i+3,1)=fq(i2*3-3+i)
enddo
a=matmul(ke,t)
fne=matmul(a,q)
do i=1,nf
j=efnl(i)
 if(j.eq.ie) then
  do k=1,6
  fne(k,1)=fne(k,1)+cfnl(6*i-6+k)
  enddo
 endif
enddo
end

!子程序6：计算节点自由度编号
subroutine cnod(nod,cfd,nfd,nn)
dimension nod(3,nn),cfd(nfd),a(3)
integer nod,cfd,nfd,nn,a,nuf
nuf=0
do i=1,nn
a=1
 do j=1,nfd
 j1=cfd(j)
  do k=1,3
   if((3*i+k-3).eq.j1) then
   nod(k,i)=0
   a(k)=0
   endif
  enddo
 enddo
 do m=1,3
  if(a(m)) then
  nuf=nuf+1
  nod(m,i)=nuf
  endif
 enddo
enddo  
end


!子程序7：计算单元节点自由度与总体自由度的关联数组
subroutine clg(lg,nod,je,ie)
dimension lg(6),nod(3,1),je(2,1)
integer lg,nod,je,ie
do i=1,2
 i1=je(i,ie)
 do j=1,3
 lg(3*i+j-3)=nod(j,i1)
 enddo
enddo
end

!子程序8：计算总刚每一列的最小非零元素的行编号
subroutine cbh(lg,ie,ii)
dimension lg(6)
integer lg,ie,ii
do i=1,6
 if(lg(i).gt.0) then
  if(lg(i).lt.ii) then
  ii=lg(i)
  endif
 endif
enddo
end

!子程序9：在已知节点力数组的情况下求解解线性方程组所需的外载数组
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


!子程序10：在已知线性方程组解的情况下求节点位移
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
