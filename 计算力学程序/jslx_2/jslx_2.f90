!求平面桁架的节点位移
program main
dimension k(8,8),ke(4,4),x(2,3),y(2,3),lg(4),co(3),so(3),l(3)  !co,so为cos,sin
real k,ke,x,y,lg,co,so,l,a,b,c,e,v1,v2  !e代表EA,u1=u2=u3=u4=v3=v4=0,F1y=-1000N,F2y=0
!原始已知数据
e=69.0e7
co(1)=0
so(1)=1
l(1)=0.5
co(2)=-sqrt(5.0)/5.0
so(2)=2*sqrt(5.0)/5.0
l(2)=sqrt(5.0)/2.0
co(3)=-co(2)
so(3)=so(2)
l(3)=l(2)
x(1,1)=1
x(1,2)=3
x(1,3)=3
x(2,1)=3
x(2,2)=5
x(2,3)=7
do i=1,2
do j=1,3
y(i,j)=x(i,j)+1
enddo
enddo
!求总刚
do i=1,3
lg(1)=x(1,i)
lg(2)=y(1,i)
lg(3)=x(2,i)
lg(4)=y(2,i)
call f(ke,co(i),so(i),l(i),e)
do m=1,4
do n=1,4
k(lg(m),lg(n))=k(lg(m),lg(n))+ke(m,n)
enddo
enddo
enddo
open(5,file='jieguo.txt')
write(5,*)'总刚k='
write(5,'(1x,8(F14.2,1X))') ((k(i,j),j=1,8),i=1,8)
!求节点位移v1,v2
v2=-1000.0/(k(2,4)-k(2,2)*k(4,4)/k(4,2))
v1=-k(4,4)*v2/k(4,2)
write(5,*) 'v1=',v1,'v2=',v2
end
!求单刚矩阵的子程序
subroutine f(ke,c,s,l,e)
dimension ke(4,4)
real ke,c,s,l,e
ke(1,1)=e*c*c/l
ke(1,2)=e*c*s/l
ke(1,3)=-ke(1,1)
ke(1,4)=-ke(1,2)
ke(2,1)=ke(1,2)
ke(2,2)=e*s*s/l
ke(2,3)=-ke(1,2)
ke(2,4)=-ke(2,2)
ke(3,1)=-ke(1,1)
ke(3,2)=-ke(1,2)
ke(3,3)=ke(1,1)
ke(3,4)=ke(1,2)
ke(4,1)=-ke(1,2)
ke(4,2)=-ke(2,2)
ke(4,3)=ke(1,2)
ke(4,4)=ke(2,2)
end
