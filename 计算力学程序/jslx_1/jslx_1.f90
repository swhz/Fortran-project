!求平面桁架整体刚度阵
program main
dimension k(12,12),ke(4,4),x(2,9),y(2,9),lg(4),co(9),so(9),l(9)  !co,so为cos,sin
real k,ke,x,y,lg,co,so,l,a,b,c,r1,r2,e   !r1,r2为角度,e代表EA
read(*,*) a,b,c,r1,r2,e
!原始已知数据
co(1)=1
so(1)=0
l(1)=a
co(2)=cos(r1)
so(2)=sin(r1)
l(2)=b
co(3)=cos(-r2)
so(3)=sin(r2)
l(3)=c
co(4)=co(2)
so(4)=so(2)
l(4)=b
co(5)=1
co(6)=1
co(8)=1
so(5)=0
so(6)=0
so(8)=0
l(5)=a
l(6)=a
l(8)=a
co(7)=co(3)
so(7)=-so(3)
l(7)=c
co(9)=co(2)
so(9)=so(2)
l(9)=b
x(1,1)=1
x(1,2)=1
x(1,3)=3
x(1,4)=3
x(1,5)=5
x(1,6)=3
x(1,7)=7
x(1,8)=7
x(1,9)=9
x(2,1)=3
x(2,2)=5
x(2,3)=5
x(2,4)=7
x(2,5)=7
x(2,6)=9
x(2,7)=9
x(2,8)=11
x(2,9)=11
do i=1,2
do j=1,9
y(i,j)=x(i,j)+1
enddo
enddo
!求总刚
do i=1,9
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
write(5,'(1x,12(f9.3,2x))') ((k(i,j),j=1,12),i=1,12)
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
