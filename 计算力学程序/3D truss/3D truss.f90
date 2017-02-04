!*****************************************空间桁架程序设计***************************************
!******************************主程序部分************************************************************************************************
program main
integer,dimension(:,:),allocatable::je              !je(2,ne)储存单元两端节点编号
integer,dimension(:),allocatable::cfd               !cfd为约束自由度数组
real,dimension(:),allocatable::cx,cy,cz,a,f,e,fne   !cx,cy,cz为节点坐标数组，a,e为单元横截面积和弹性模量，f为外载荷数组，fne为单元内力数组
real,dimension(:,:),allocatable::k                  !k为总刚
dimension ke(6,6),lg(6)                             !ke为单刚，lg为单元定位数组 
real ke,lx,ly,lz,l                                  !lx,ly,lz为单元定位向量，l为单元长度
integer nn,ne,nd,nfd,lg,c,i,j,p                     !nn为节点数，ne为单元数，nd为原始总自由度数（不考虑约束），nfd为被约束的自由度数目 
 
open(5,file='shuju.txt')

read(5,*) nn,ne,nd,nfd            !**************
open(6,file='result.txt')                       !
write(6,7) nn,ne,nd,nfd                         !*****输入输出节点数，单元数，原始总自由度数，被约束的自由度数目
7 format(2x,'nn  ne  nd  nfd'/4i4)!*************!  

!给动态数组赋予确定边界范围   
allocate(cx(nn),cy(nn),cz(nn),je(2,ne),fne(ne),a(ne),e(ne),k(nd,nd),cfd(nfd),f(nd))

write(6,*)  !输出一行空格

read(5,*) (cx(i),cy(i),cz(i),i=1,nn) !***************
write(6,8)                                          !
8 format(3x,'NODE   X-COORD   Y-COORD   Z-COORD')   !*****输入输出节点坐标
write(6,9) (i,cx(i),cy(i),cz(i),i=1,nn)             !
9 format(2x,i4,1x,3f10.5)            !***************

write(6,*)    !输出一行空格

read(5,*) ((je(i,j),i=1,2),e(j),a(j),j=1,ne)     !****************
write(6,10)                                                      !
10 format(3x,'ELEMENT   NODE1   NODE2         E         AREA')   !*****输入输出单元节点编号，弹性模量，横截面积
write(6,11) (j,(je(i,j),i=1,2),e(j),a(j),j=1,ne)                 !
11 format(3x,i4,3x,2i7,2x,2e13.5)                !****************

write(6,*)    !输出一行空格

read(5,*) (cfd(i),i=1,nfd)   !********
write(6,12)                          !
12 format(5x,'约束自由度号:')        !*****输入输出被约束的原始自由度编号
write(6,13) (cfd(i),i=1,nfd)         !
13 format(2x,10i3)           !********

write(6,*)    !输出一行空格

read(5,*) (f(i),i=1,nd)            !************************
write(6,14)                                                !
14 format(3x,'NODE      X-LOAD       Y-LOAD       Z-LOAD') !******输入输出节点载荷
write(6,15) (i,f(3*i-2),f(3*i-1),f(3*i),i=1,nn)            !
15 format(3x,i3,1x,3e13.5)         !************************

k=0.0        !k置零

do i=1,ne         !************************************************************                                                                     !
lg=0                                                                          !                                                                    !
lg(1)=(je(1,i)-1)*3+1                                                         !
lg(2)=lg(1)+1                                                                 !
lg(3)=lg(1)+2                                                                 !
lg(4)=(je(2,i)-1)*3+1                                                         !
lg(5)=lg(4)+1                                                                 !
lg(6)=lg(4)+2                                                                 !*****组装总刚 
call elxyz(lx,ly,lz,l,i,je,cx,cy,cz)  !调用子程序elxyz计算单元lx，ly，lz，l   !
call cdg(ke,lx,ly,lz,l,e(i)*a(i))     !调用子程序cdg计算单刚                  !
do m=1,6                                                                      !
do n=1,6                                                                      !
k(lg(m),lg(n))=k(lg(m),lg(n))+ke(m,n)                                         !
enddo                                                                         !
enddo                                                                         !
enddo             !************************************************************  
                                                         
do i=1,nfd         !*************
c=cfd(i)                        !
k(c,c)=(k(c,c)+1)*1.0e15        !*****乘大数法引入约束条件
enddo              !*************         

call gauss(k,f,nd) !调用gauss子程序解线性方程组

do i=1,nfd !*****
c=cfd(i)        !
f(c)=0          !*****约束位移置零
enddo      !*****  

write(6,*)     !输出一行空格
 
write(6,16)                        !***********************
16 format(3x,'NODE      X-DISP      Y-DISP      Z-DISP')  !
write(6,17) (i,f(3*i-2),f(3*i-1),f(3*i),i=1,nn)           !*****输出节点位移
17 format(3x,i3,1x,3e13.5)         !***********************
                        
do p=1,ne                                              !********************************                      
call elxyz(lx,ly,lz,l,p,je,cx,cy,cz)                                                   !
i=je(1,p)                                                                              !
j=je(2,p)                                                                              !*****计算单元内力
fne(p)=-e(p)*a(p)*(lx*(f(3*i-2)-f(3*j-2))+ly*(f(3*i-1)-f(3*j-1))+lz*(f(3*i)-f(3*j)))/l !
enddo                                                  !********************************

write(6,*)      !输出一行空格

write(6,18)           !*****************
18 format(3x,'ELEMENT      FORCE')     !
write(6,19) (i,fne(i),i=1,ne)          !*****输出单元内力
19 format(3x,i4,3x,e13.5)              !
end                   !***************** 

!*********************************子程序部分*************************************************************************************
!子程序1：计算单元长度l，单元lx,ly,lz
subroutine elxyz(lx,ly,lz,l,ie,je,x,y,z)
dimension je(2,1),x(1),y(1),z(1)
real lx,ly,lz,l,s1,s2,s3,x,y,z
integer i,j,je
i=je(1,ie)
j=je(2,ie)
s1=x(j)-x(i)
s2=y(j)-y(i)
s3=z(j)-z(i)
l=sqrt(s1*s1+s2*s2+s3*s3)
lx=s1/l
ly=s2/l
lz=s3/l
return 
end
!子程序2：3D truss单刚计算子程序
subroutine cdg(ke,lx,ly,lz,l,ea)
dimension ke(6,6)
real ke,lx,ly,lz,l,ea
ke(1,1)=ea*lx*lx/l
ke(2,2)=ea*ly*ly/l
ke(3,3)=ea*lz*lz/l
do i=4,6
ke(i,i)=ke(i-3,i-3)
enddo
ke(1,2)=ea*lx*ly/l
ke(1,3)=ea*lx*lz/l
do i=4,6
ke(1,i)=-ke(1,i-3)
enddo
ke(2,3)=ea*ly*lz/l
ke(2,4)=-ke(1,2)
ke(2,5)=-ke(2,2)
ke(2,6)=-ke(2,3)
ke(3,4)=-ke(1,3)
ke(3,5)=ke(2,6)
ke(3,6)=-ke(3,3)
ke(4,5)=ke(1,2)
ke(4,6)=ke(1,3)
ke(5,6)=ea*ly*lz/l
do i=1,6
do j=1,6
ke(j,i)=ke(i,j)
enddo
enddo
return
end
!子程序3：高斯消去法
subroutine gauss(a,b,n)
dimension a(n,n),b(n)
do i=1,n
i1=i+1
do j=i1,n
a(i,j)=a(i,j)/a(i,i)
enddo
b(i)=b(i)/a(i,i)
a(i,i)=1.0
do j=i1,n
do m=i1,n
a(j,m)=a(j,m)-a(j,i)*a(i,m)
enddo
b(j)=b(j)-a(j,i)*b(i)
enddo
enddo
do i=n-1,1,-1
do j=i+1,n
b(i)=b(i)-a(i,j)*b(j)
enddo
enddo
return
end
