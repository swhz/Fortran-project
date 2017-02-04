!*****************************************�ռ���ܳ������***************************************
!******************************�����򲿷�************************************************************************************************
program main
integer,dimension(:,:),allocatable::je              !je(2,ne)���浥Ԫ���˽ڵ���
integer,dimension(:),allocatable::cfd               !cfdΪԼ�����ɶ�����
real,dimension(:),allocatable::cx,cy,cz,a,f,e,fne   !cx,cy,czΪ�ڵ��������飬a,eΪ��Ԫ�������͵���ģ����fΪ���غ����飬fneΪ��Ԫ��������
real,dimension(:,:),allocatable::k                  !kΪ�ܸ�
dimension ke(6,6),lg(6)                             !keΪ���գ�lgΪ��Ԫ��λ���� 
real ke,lx,ly,lz,l                                  !lx,ly,lzΪ��Ԫ��λ������lΪ��Ԫ����
integer nn,ne,nd,nfd,lg,c,i,j,p                     !nnΪ�ڵ�����neΪ��Ԫ����ndΪԭʼ�����ɶ�����������Լ������nfdΪ��Լ�������ɶ���Ŀ 
 
open(5,file='shuju.txt')

read(5,*) nn,ne,nd,nfd            !**************
open(6,file='result.txt')                       !
write(6,7) nn,ne,nd,nfd                         !*****��������ڵ�������Ԫ����ԭʼ�����ɶ�������Լ�������ɶ���Ŀ
7 format(2x,'nn  ne  nd  nfd'/4i4)!*************!  

!����̬���鸳��ȷ���߽緶Χ   
allocate(cx(nn),cy(nn),cz(nn),je(2,ne),fne(ne),a(ne),e(ne),k(nd,nd),cfd(nfd),f(nd))

write(6,*)  !���һ�пո�

read(5,*) (cx(i),cy(i),cz(i),i=1,nn) !***************
write(6,8)                                          !
8 format(3x,'NODE   X-COORD   Y-COORD   Z-COORD')   !*****��������ڵ�����
write(6,9) (i,cx(i),cy(i),cz(i),i=1,nn)             !
9 format(2x,i4,1x,3f10.5)            !***************

write(6,*)    !���һ�пո�

read(5,*) ((je(i,j),i=1,2),e(j),a(j),j=1,ne)     !****************
write(6,10)                                                      !
10 format(3x,'ELEMENT   NODE1   NODE2         E         AREA')   !*****���������Ԫ�ڵ��ţ�����ģ����������
write(6,11) (j,(je(i,j),i=1,2),e(j),a(j),j=1,ne)                 !
11 format(3x,i4,3x,2i7,2x,2e13.5)                !****************

write(6,*)    !���һ�пո�

read(5,*) (cfd(i),i=1,nfd)   !********
write(6,12)                          !
12 format(5x,'Լ�����ɶȺ�:')        !*****���������Լ����ԭʼ���ɶȱ��
write(6,13) (cfd(i),i=1,nfd)         !
13 format(2x,10i3)           !********

write(6,*)    !���һ�пո�

read(5,*) (f(i),i=1,nd)            !************************
write(6,14)                                                !
14 format(3x,'NODE      X-LOAD       Y-LOAD       Z-LOAD') !******��������ڵ��غ�
write(6,15) (i,f(3*i-2),f(3*i-1),f(3*i),i=1,nn)            !
15 format(3x,i3,1x,3e13.5)         !************************

k=0.0        !k����

do i=1,ne         !************************************************************                                                                     !
lg=0                                                                          !                                                                    !
lg(1)=(je(1,i)-1)*3+1                                                         !
lg(2)=lg(1)+1                                                                 !
lg(3)=lg(1)+2                                                                 !
lg(4)=(je(2,i)-1)*3+1                                                         !
lg(5)=lg(4)+1                                                                 !
lg(6)=lg(4)+2                                                                 !*****��װ�ܸ� 
call elxyz(lx,ly,lz,l,i,je,cx,cy,cz)  !�����ӳ���elxyz���㵥Ԫlx��ly��lz��l   !
call cdg(ke,lx,ly,lz,l,e(i)*a(i))     !�����ӳ���cdg���㵥��                  !
do m=1,6                                                                      !
do n=1,6                                                                      !
k(lg(m),lg(n))=k(lg(m),lg(n))+ke(m,n)                                         !
enddo                                                                         !
enddo                                                                         !
enddo             !************************************************************  
                                                         
do i=1,nfd         !*************
c=cfd(i)                        !
k(c,c)=(k(c,c)+1)*1.0e15        !*****�˴���������Լ������
enddo              !*************         

call gauss(k,f,nd) !����gauss�ӳ�������Է�����

do i=1,nfd !*****
c=cfd(i)        !
f(c)=0          !*****Լ��λ������
enddo      !*****  

write(6,*)     !���һ�пո�
 
write(6,16)                        !***********************
16 format(3x,'NODE      X-DISP      Y-DISP      Z-DISP')  !
write(6,17) (i,f(3*i-2),f(3*i-1),f(3*i),i=1,nn)           !*****����ڵ�λ��
17 format(3x,i3,1x,3e13.5)         !***********************
                        
do p=1,ne                                              !********************************                      
call elxyz(lx,ly,lz,l,p,je,cx,cy,cz)                                                   !
i=je(1,p)                                                                              !
j=je(2,p)                                                                              !*****���㵥Ԫ����
fne(p)=-e(p)*a(p)*(lx*(f(3*i-2)-f(3*j-2))+ly*(f(3*i-1)-f(3*j-1))+lz*(f(3*i)-f(3*j)))/l !
enddo                                                  !********************************

write(6,*)      !���һ�пո�

write(6,18)           !*****************
18 format(3x,'ELEMENT      FORCE')     !
write(6,19) (i,fne(i),i=1,ne)          !*****�����Ԫ����
19 format(3x,i4,3x,e13.5)              !
end                   !***************** 

!*********************************�ӳ��򲿷�*************************************************************************************
!�ӳ���1�����㵥Ԫ����l����Ԫlx,ly,lz
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
!�ӳ���2��3D truss���ռ����ӳ���
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
!�ӳ���3����˹��ȥ��
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
