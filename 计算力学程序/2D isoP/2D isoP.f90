!***************************二维四节点等参元*********************
program main
!je为储存单元节点编号数组，nod为单元节点自由度编号数组
integer,dimension(:,:),allocatable::je,nod             
!cfd为约束自由度数组，p为天空线一维定位数组
integer,dimension(:),allocatable::cfd,p        
!cx,cy为节点坐标数组，f为求解所需外载荷数组，fq为节点力数组，sky为刚度矩阵的天空线储存数组
real,dimension(:),allocatable::cx,cy,f,fq,sky                  
!ke为单刚，lg为单元节点自由度定位数组，de为弹性矩阵，ni为形函数，nid为形函数的导数，niz为形函数在整体坐标系下的导数
!jac为雅克比矩阵，jac1为雅克比矩阵的逆，b为[B]矩阵，st为单元应力数组，u为单元节点位移数组
dimension ke(8,8),lg(8),de(3,3),ni(8,2),nid(2,4),niz(2,4),jac(2,2),jac1(2,2),b(3,8),st(3,1),u(8,1)
real ke,ax,by,detj,st1,st2,sta,e,v,h,de,ni,nid,niz,jac,jac1,b,st,u,a1
!nn为节点数，ne为单元数，nd为总自由度数（包含被约束住的自由度），nfd为被约束的自由度数,nfq为受外力节点数，ii为与总体自由度i的所以相关自由度的最小值，nwk为天空线数组长度                             
integer nn,ne,nd,nfd,nfq,lg,ii,nwk
                    
open(5,file='shuju.txt')

!*****输入输出节点数，单元数，总自由度数（原始自由度数，考虑约束情况），受力节点数，弹性模量，泊松比，厚度
read(5,*) nn,ne,nd,nfd,nfq,e,v,h                                           
open(6,file='result.txt')                          
write(6,7) nn,ne,nd,nfd,nfq,e,v,h                         
7 format(2x,'nn  ne  nd  nfd nfq           e           v           h'/5i4,3x,3e13.5)  

!对动态数组赋予确定范围宽度     
allocate(cx(nn),cy(nn),je(4,ne),nod(2,nn),cfd(nfd),f(nd-nfd),fq(nd),p(nd-nfd))

write(6,*)           !输出一行空格

read(5,*) (cx(i),cy(i),i=1,nn)             !***************
write(6,8)                                                !
8 format(3x,'NODE   X-COORD   Y-COORD')                   !*****输入输出节点坐标
write(6,9) (i,cx(i),cy(i),i=1,nn)                         !
9 format(2x,i4,1x,2f10.5)                  !*************** 
 
write(6,*)           !输出一行空格

!******************输入输出单元节点编号**************** 
read(5,*) ((je(i,j),i=1,4),j=1,ne)      
write(6,10)                                                             !
10 format(3x,'ELEMENT    NODE1  NODE2  NODE3  NODE4')  
write(6,11) (j,(je(i,j),i=1,4),j=1,ne)                  
11 format(3x,i4,3x,4i7) 
          
write(6,*)           !输出一行空格

read(5,*) (cfd(i),i=1,nfd)   !********
write(6,12)                          !
12 format(5x,'约束自由度号:')        !*****输入输出被约束的自由度在原始自由度中的编号
write(6,13) (cfd(i),i=1,nfd)         !
13 format(2x,10i3)           !********

write(6,*)           !输出一行空格

!***********输入输出节点载荷***********************
fq=0
write(6,*) '                       单元节点载荷'              
read(5,*) (j,fq(2*j-1),fq(2*j),i=1,nfq)                                                                       
write(6,14)                                                                                      
14 format(3x,'NODE      X-LOAD      Y-LOAD')   
write(6,15) (i,fq(2*i-1),fq(2*i),i=1,nn)                     
15 format(3x,i3,1x,2e13.5)                                     

!******计算各节点的新自由度编号（被约束住的自由度赋零）**********
call cnod(nod,cfd,nfd,nn)


!**********计算总刚矩阵处于上三角的每一列最小非零元素的行号和每一列在天空线矩阵中的元素个数
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


!***********计算天空线一维定位数组并确定天空线数组长度************
p(1)=1
do i=2,nd-nfd
p(i)=p(i-1)+p(i)
enddo
nwk=p(nd-nfd)
allocate(sky(nwk))


!************计算总刚的天空线存储数组**************
sky=0.0
call cde(e,v,de)

do i=1,ne

call clg(lg,nod,je,i)
!*********计算单刚（结合子程序加以循环）**************
ke=0.0                                                                                  
  do j=1,2
   do l=1,2
   !******赋值高斯点（取2*2的二维高斯点）*****
    ax=((-1)**l)/sqrt(3.0)
    by=((-1)**j)/sqrt(3.0)
    call jbxhs(ax,by,ni,nid)
    call jacobian(jac,jac1,nid,i,je,cx,cy,detj)
    call ztxhs(niz,jac1,nid)
    call cb(niz,b)
    call ck(ke,b,detj,de,h)
   enddo
  enddo   
!***********确定sky*************                                               
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

!***********调用子程序计算修正的外载f（去除零自由度所对应的外载）*******
f=0.0
call cf(f,fq,cfd,nfd,nd)

!*********调用子程序解线性方程组求解位移*********
call skylineLDL(sky,p,f,nwk,nd-nfd)

!**********调用子程序计算修正后的各节点位移*******
fq=0.0
call cfq(f,fq,cfd,nfd,nd)
                                                                      

!***************输出节点位移*************
write(6,*) '                   节点位移'        
write(6,20)                                                                                      
20 format(3x,'NODE      X-DISP       Y-DISP')  
write(6,21) (i,fq(2*i-1),fq(2*i),i=1,nn)                                                        
21 format(3x,i3,1x,2e13.5)                      


!*****************计算单元应力并求出主应力和主方向（取单元局部坐标系下的中点计算）***********
write(6,*) 
write(6,*) '                                           单元应力 '
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
                        
!*********************************子程序部分*************************************************************************************
!子程序1：计算弹性矩阵
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
!子程序2：计算形函数及其在局部坐标系下的导数
subroutine jbxhs(ax,by,ni,nid)
dimension ni(8,2),nid(2,4)
real ax,by,ni,nid
ni=0.0
nid=0.0
!******计算ni
ni(1,1)=0.25*(1-ax)*(1-by)
ni(2,2)=ni(1,1)
ni(3,1)=0.25*(1+ax)*(1-by)
ni(4,2)=ni(3,1)
ni(5,1)=0.25*(1+ax)*(1+by)
ni(6,2)=ni(5,1)
ni(7,1)=0.25*(1-ax)*(1+by)
ni(8,2)=ni(7,1)
!*******计算nid
nid(1,1)=-0.25*(1-by)
nid(1,2)=-nid(1,1)
nid(1,3)=0.25*(1+by)
nid(1,4)=-nid(1,3)
nid(2,1)=-0.25*(1-ax)
nid(2,2)=-0.25*(1+ax)
nid(2,3)=-nid(2,2)
nid(2,4)=-nid(2,1)
end
!子程序3：计算雅可比矩阵及其逆
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
!子程序4：计算形函数在整体坐标系下的导数
subroutine ztxhs(niz,jac1,nid)
dimension niz(2,4),nid(2,4),jac1(2,2)
real niz,nid,jac1
niz=matmul(jac1,nid)
end
!子程序5：计算[B]矩阵
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
!子程序6：计算单元刚度矩阵[k](结合主程序的循环语句才能得到单刚)
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


!子程序7：计算单元应力及其主应力和主方向
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


!子程序8：计算节点自由度编号
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


!子程序9：计算单元节点自由度与总体自由度的关联数组
subroutine clg(lg,nod,je,ie)
dimension lg(8),nod(2,1),je(4,1)
integer lg,nod,je,ie
do i=1,4
 i1=je(i,ie)
 lg(2*i-1)=nod(1,i1)
 lg(2*i)=nod(2,i1)
enddo
end

!子程序10：计算总刚每一列的最小非零元素的行编号
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

!子程序11：在已知节点力数组的情况下求解解线性方程组所需的外载数组
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


!子程序12：在已知线性方程组解的情况下求节点位移
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


!子程序13：天空线的LDL分解法（解线性方程组）
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
