!***************************��ά�Ľڵ�Ȳ�Ԫ*********************
program main
!jeΪ���浥Ԫ�ڵ������飬nodΪ��Ԫ�ڵ����ɶȱ������
integer,dimension(:,:),allocatable::je,nod             
!cfdΪԼ�����ɶ����飬pΪ�����һά��λ����
integer,dimension(:),allocatable::cfd,p        
!cx,cyΪ�ڵ��������飬fΪ����������غ����飬fqΪ�ڵ������飬skyΪ�նȾ��������ߴ�������
real,dimension(:),allocatable::cx,cy,f,fq,sky                  
!keΪ���գ�lgΪ��Ԫ�ڵ����ɶȶ�λ���飬deΪ���Ծ���niΪ�κ�����nidΪ�κ����ĵ�����nizΪ�κ�������������ϵ�µĵ���
!jacΪ�ſ˱Ⱦ���jac1Ϊ�ſ˱Ⱦ�����棬bΪ[B]����stΪ��ԪӦ�����飬uΪ��Ԫ�ڵ�λ������
dimension ke(8,8),lg(8),de(3,3),ni(8,2),nid(2,4),niz(2,4),jac(2,2),jac1(2,2),b(3,8),st(3,1),u(8,1)
real ke,ax,by,detj,st1,st2,sta,e,v,h,de,ni,nid,niz,jac,jac1,b,st,u,a1
!nnΪ�ڵ�����neΪ��Ԫ����ndΪ�����ɶ�����������Լ��ס�����ɶȣ���nfdΪ��Լ�������ɶ���,nfqΪ�������ڵ�����iiΪ���������ɶ�i������������ɶȵ���Сֵ��nwkΪ��������鳤��                             
integer nn,ne,nd,nfd,nfq,lg,ii,nwk
                    
open(5,file='shuju.txt')

!*****��������ڵ�������Ԫ���������ɶ�����ԭʼ���ɶ���������Լ��������������ڵ���������ģ�������ɱȣ����
read(5,*) nn,ne,nd,nfd,nfq,e,v,h                                           
open(6,file='result.txt')                          
write(6,7) nn,ne,nd,nfd,nfq,e,v,h                         
7 format(2x,'nn  ne  nd  nfd nfq           e           v           h'/5i4,3x,3e13.5)  

!�Զ�̬���鸳��ȷ����Χ���     
allocate(cx(nn),cy(nn),je(4,ne),nod(2,nn),cfd(nfd),f(nd-nfd),fq(nd),p(nd-nfd))

write(6,*)           !���һ�пո�

read(5,*) (cx(i),cy(i),i=1,nn)             !***************
write(6,8)                                                !
8 format(3x,'NODE   X-COORD   Y-COORD')                   !*****��������ڵ�����
write(6,9) (i,cx(i),cy(i),i=1,nn)                         !
9 format(2x,i4,1x,2f10.5)                  !*************** 
 
write(6,*)           !���һ�пո�

!******************���������Ԫ�ڵ���**************** 
read(5,*) ((je(i,j),i=1,4),j=1,ne)      
write(6,10)                                                             !
10 format(3x,'ELEMENT    NODE1  NODE2  NODE3  NODE4')  
write(6,11) (j,(je(i,j),i=1,4),j=1,ne)                  
11 format(3x,i4,3x,4i7) 
          
write(6,*)           !���һ�пո�

read(5,*) (cfd(i),i=1,nfd)   !********
write(6,12)                          !
12 format(5x,'Լ�����ɶȺ�:')        !*****���������Լ�������ɶ���ԭʼ���ɶ��еı��
write(6,13) (cfd(i),i=1,nfd)         !
13 format(2x,10i3)           !********

write(6,*)           !���һ�пո�

!***********��������ڵ��غ�***********************
fq=0
write(6,*) '                       ��Ԫ�ڵ��غ�'              
read(5,*) (j,fq(2*j-1),fq(2*j),i=1,nfq)                                                                       
write(6,14)                                                                                      
14 format(3x,'NODE      X-LOAD      Y-LOAD')   
write(6,15) (i,fq(2*i-1),fq(2*i),i=1,nn)                     
15 format(3x,i3,1x,2e13.5)                                     

!******������ڵ�������ɶȱ�ţ���Լ��ס�����ɶȸ��㣩**********
call cnod(nod,cfd,nfd,nn)


!**********�����ܸվ����������ǵ�ÿһ����С����Ԫ�ص��кź�ÿһ��������߾����е�Ԫ�ظ���
do i=1,nd-nfd
 ii=nd-nfd
 do j=1,ne
 call clg(lg,nod,je,j) 
  do m=1,8
  i1=lg(m)
   if(i1.eq.i) then
   call cbh(lg,i,ii)
   endif
  enddo
 enddo
p(i)=i-ii+1
enddo


!***********���������һά��λ���鲢ȷ����������鳤��************
p(1)=1
do i=2,nd-nfd
p(i)=p(i-1)+p(i)
enddo
nwk=p(nd-nfd)
allocate(sky(nwk))


!************�����ܸյ�����ߴ洢����**************
sky=0.0
call cde(e,v,de)

do i=1,ne

call clg(lg,nod,je,i)
!*********���㵥�գ�����ӳ������ѭ����**************
ke=0.0                                                                                  
  do j=1,2
   do l=1,2
   !******��ֵ��˹�㣨ȡ2*2�Ķ�ά��˹�㣩*****
    ax=((-1)**l)/sqrt(3.0)
    by=((-1)**j)/sqrt(3.0)
    call jbxhs(ax,by,ni,nid)
    call jacobian(jac,jac1,nid,i,je,cx,cy,detj)
    call ztxhs(niz,jac1,nid)
    call cb(niz,b)
    call ck(ke,b,detj,de,h)
   enddo
  enddo   
!***********ȷ��sky*************                                               
  do m=1,8
   if(lg(m).gt.0) then
   j=p(lg(m))
   sky(j)=sky(j)+ke(m,m)                                                                                  
    do n=1,8
	 if(lg(n).gt.0) then                                                                                 
      if(lg(n).lt.lg(m)) then
	  l=p(lg(m))+lg(n)-lg(m) 
	  sky(l)=sky(l)+ke(m,n)
	  endif 
	 endif                                                     
    enddo
   endif	                                                                                       
  enddo  
                                                                                  
enddo

!***********�����ӳ����������������f��ȥ�������ɶ�����Ӧ�����أ�*******
f=0.0
call cf(f,fq,cfd,nfd,nd)

!*********�����ӳ�������Է��������λ��*********
call skylineLDL(sky,p,f,nwk,nd-nfd)

!**********�����ӳ������������ĸ��ڵ�λ��*******
fq=0.0
call cfq(f,fq,cfd,nfd,nd)
                                                                      

!***************����ڵ�λ��*************
write(6,*) '                   �ڵ�λ��'        
write(6,20)                                                                                      
20 format(3x,'NODE      X-DISP       Y-DISP')  
write(6,21) (i,fq(2*i-1),fq(2*i),i=1,nn)                                                        
21 format(3x,i3,1x,2e13.5)                      


!*****************���㵥ԪӦ���������Ӧ����������ȡ��Ԫ�ֲ�����ϵ�µ��е���㣩***********
write(6,*) 
write(6,*) '                                           ��ԪӦ�� '
write(6,22)
22 format(3x,'ELEMENT        STX          STY         STXY          ST1          ST2          STA')
do i=1,ne
  do j=1,4
   j1=je(j,i)
   u(2*j-1,1)=fq(2*j1-1)
   u(2*j,1)=fq(2*j1)
  enddo
    ax=0.0
    by=0.0
    call jbxhs(ax,by,ni,nid)
    call jacobian(jac,jac1,nid,i,je,cx,cy,detj)
    call ztxhs(niz,jac1,nid)
    call cb(niz,b)
    call cst(st,st1,st2,sta,de,u,b)  
write(6,23) i,(st(j,1),j=1,3),st1,st2,sta
23 format(3x,i4,3x,6e13.5)
enddo
end
                        
!*********************************�ӳ��򲿷�*************************************************************************************
!�ӳ���1�����㵯�Ծ���
subroutine cde(e,v,de)
dimension de(3,3)
real e,v,de,a
a=e/(1-v*v)
de(1,1)=a
de(1,2)=a*v
de(1,3)=0.0
de(2,1)=de(1,2)
de(2,2)=a
de(2,3)=0.0
de(3,1)=0.0
de(3,2)=0.0
de(3,3)=a*(1-v)/2.0
end
!�ӳ���2�������κ��������ھֲ�����ϵ�µĵ���
subroutine jbxhs(ax,by,ni,nid)
dimension ni(8,2),nid(2,4)
real ax,by,ni,nid
ni=0.0
nid=0.0
!******����ni
ni(1,1)=0.25*(1-ax)*(1-by)
ni(2,2)=ni(1,1)
ni(3,1)=0.25*(1+ax)*(1-by)
ni(4,2)=ni(3,1)
ni(5,1)=0.25*(1+ax)*(1+by)
ni(6,2)=ni(5,1)
ni(7,1)=0.25*(1-ax)*(1+by)
ni(8,2)=ni(7,1)
!*******����nid
nid(1,1)=-0.25*(1-by)
nid(1,2)=-nid(1,1)
nid(1,3)=0.25*(1+by)
nid(1,4)=-nid(1,3)
nid(2,1)=-0.25*(1-ax)
nid(2,2)=-0.25*(1+ax)
nid(2,3)=-nid(2,2)
nid(2,4)=-nid(2,1)
end
!�ӳ���3�������ſɱȾ�������
subroutine jacobian(jac,jac1,nid,ie,je,cx,cy,detj)
dimension je(4,1),jac(2,2),jac1(2,2),nid(2,4),xy(4,2),cx(1),cy(1)
real jac,jac1,nid,xy,cx,cy,detj
integer je,ie
do i=1,4
	j=je(i,ie)
	xy(i,1)=cx(j)
	xy(i,2)=cy(j)
enddo
jac=matmul(nid,xy)
detj=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
jac1(1,1)=jac(2,2)/detj
jac1(1,2)=-jac(1,2)/detj
jac1(2,1)=-jac(2,1)/detj
jac1(2,2)=jac(1,1)/detj
end
!�ӳ���4�������κ�������������ϵ�µĵ���
subroutine ztxhs(niz,jac1,nid)
dimension niz(2,4),nid(2,4),jac1(2,2)
real niz,nid,jac1
niz=matmul(jac1,nid)
end
!�ӳ���5������[B]����
subroutine cb(niz,b)
dimension niz(2,4),b(3,8)
real niz,b
b=0.0
do i=1,4
	b(1,2*i-1)=niz(1,i)
	b(2,2*i)=niz(2,i)
	b(3,2*i-1)=niz(2,i)
	b(3,2*i)=niz(1,i)
enddo
end
!�ӳ���6�����㵥Ԫ�նȾ���[k](����������ѭ�������ܵõ�����)
subroutine ck(ke,b,detj,de,h)
dimension ke(8,8),b(3,8),b1(8,3),de(3,3),ke1(8,3),ke2(8,8)
real ke,ke1,ke2,b,b1,de,detj,h
b1=transpose(b)
ke1=matmul(b1,de)
ke2=matmul(ke1,b)
do i=1,8
 do j=1,8
   ke(i,j)=ke(i,j)+h*detj*ke2(i,j)
 enddo
enddo
end


!�ӳ���7�����㵥ԪӦ��������Ӧ����������
subroutine cst(st,st1,st2,sta,de,u,b)
dimension st(3,1),de(3,3),u(8,1),b(3,8),e(3,1)
real st,de,u,b,e,st1,st2,sta,h1,h2
e=matmul(b,u)
st=matmul(de,e)
h1=(st(1,1)+st(2,1))/2.0
h2=sqrt(st(3,1)*st(3,1)+((st(1,1)-st(2,1))**2)/4.0)
st1=h1+h2
st2=h1-h2
sta=-0.5*atan(2*st(3,1)/(st(1,1)-st(2,1)))
end


!�ӳ���8������ڵ����ɶȱ��
subroutine cnod(nod,cfd,nfd,nn)
dimension nod(2,nn),cfd(nfd)
integer nod,cfd,nfd,nn,a,b,nuf
nuf=0
do i=1,nn
a=1
b=1
 do j=1,nfd
 j1=cfd(j)
  if((2*i-1).eq.j1) then
  nod(1,i)=0
  a=0
  endif
  if((2*i).eq.j1) then
  nod(2,i)=0
  b=0
  endif
 enddo
 if(a) then
 nuf=nuf+1
 nod(1,i)=nuf
 endif
 if(b) then
 nuf=nuf+1
 nod(2,i)=nuf
 endif
enddo  
end


!�ӳ���9�����㵥Ԫ�ڵ����ɶ����������ɶȵĹ�������
subroutine clg(lg,nod,je,ie)
dimension lg(8),nod(2,1),je(4,1)
integer lg,nod,je,ie
do i=1,4
 i1=je(i,ie)
 lg(2*i-1)=nod(1,i1)
 lg(2*i)=nod(2,i1)
enddo
end

!�ӳ���10�������ܸ�ÿһ�е���С����Ԫ�ص��б��
subroutine cbh(lg,ie,ii)
dimension lg(8)
integer lg,ie,ii
do i=1,8
 if(lg(i)) then
  if(lg(i).lt.ii) then
  ii=lg(i)
  endif
 endif
enddo
end

!�ӳ���11������֪�ڵ��������������������Է������������������
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


!�ӳ���12������֪���Է��������������ڵ�λ��
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


!�ӳ���13������ߵ�LDL�ֽⷨ�������Է����飩
subroutine skylineLDL(sky,P,R,Nsky,NP)
real sky(Nsky),R(NP)
integer P(NP),nr(NP)
real sigma,term

nr(1)=1
do,i=2,NP
nr(i)=i-(P(i)-P(i-1))+1
enddo
!************���Ƿֽ⣬��L(ij)��d(ij)*******
do,j=2,NP
!==================��g(ij)================== 
 do,i=nr(j),j-1
 sigma=0
  nrmax=max(nr(i),nr(j))
   do,k=nrmax,i-1
   term=sky(P(i)-(i-k))*sky(P(j)-(j-k))
   sigma=sigma+term
   enddo
 sky(P(j)-(j-i))=sky(P(j)-(j-i))-sigma
 enddo
!=================��d(jj)===================
 sigma=0
  do,k=nr(j),j-1
  term=(sky(P(j)-(j-k)))**2/sky(P(k))
  sigma=sigma+term
  enddo
 sky(P(j))=sky(P(j))-sigma
!==================��L(ij)================== 
 do,i=nr(j),j-1
 sky(P(j)-(j-i))=sky(P(j)-(j-i))/sky(P(i))
 enddo
enddo
!************ǰ������{V}**********
do,j=2,NP
 sigma=0
  do,i=nr(j),j-1
  term=sky(P(j)-(j-i))*R(i)
  sigma=sigma+term
  enddo
 R(j)=R(j)-sigma
enddo
!************�ش�����{U}***********
!��{V-bar}��ֵ
do,i=1,NP
 R(i)=R(i)/sky(P(i))
enddo
!��{U}��ֵ
do,i=NP,2,-1
 do,j=nr(i),i-1
 R(j)=R(j)-sky(P(i)-(i-j))*R(i)
 enddo
enddo
end
