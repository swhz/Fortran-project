!*************************�����޲�ַ�����Ĵ��Ӷ�����******************
!************�ó����������ı߶��Ǽ�֧�ľ��ΰ�*************
program main
!kΪ������Է������ϵ������wΪ�ӶȾ���faiΪfaiֵ�������
!l,p,r,s,zΪ�Ӷȵ�ϵ���������ı仯����
doubleprecision,dimension(:,:),allocatable::k,w,w1,fai,l,p,r,s,z
!bΪ���Է�������ұ������
doubleprecision,dimension(:),allocatable::b
!x��yΪ���ΰ�ĳ��ȺͿ�ȣ�wcmΪ�������õ��Ķ�Ӧ�ӶȵĲ�ֵ�����ֵ�ľ���ֵ��wmaxΪw�����ֵ��wbΪѭ���ж���
doubleprecision h,t,D,E,v,q,x,y,wmax,wcm,wb
integer m,n,cs
open(5,file='shuju.txt')
open(6,file='w.txt')
read(5,*) x,y
write(6,*) x,y
read(5,*) m,n
write(6,*) m,n
read(5,*) E,v,t,q
write(6,*) E,v,t,q 
h=x/m
!ȷ����̬����ķ�Χ���
allocate(k((m-1)*(n-1),(m-1)*(n-1)),w(m+1,n+1),w1(m+1,n+1),fai(m+1,n+1),b((m-1)*(n-1)),&
l((m-1),(n-1)),p((m-1),(n-1)),z((m-1),(n-1)),r((m-1),(n-1)),s((m-1),(n-1)))
D=E*t*t*t/(12.0*(1.0-v*v))
w=0.0
fai=0.0
wcm=1.0
wb=0
cs=0
do while(wcm.gt.wb)
w1=w
call cxs(l,p,r,s,z,fai,t,D,m,n)
call ckw(k,l,p,r,s,z,m,n)
b=q*(h**4)/D
call gauss(k,b,(m-1)*(n-1))
call cbwfai(w,b,m,n)
call cwcm(wcm,w,w1,m,n)
call cwmax(wmax,w,m,n)
wb=wmax*1.0e-6
call ckfai(k,m,n)
call cb(b,w,m,n,E)
call gauss(k,b,(m-1)*(n-1))
call cbwfai(fai,b,m,n)
cs=cs+1
write(*,*) cs,wmax,wcm    !���ѭ��������ÿ��ѭ����Ӧ������Ӷȣ�ǰ������ѭ����Ӧ�Ӷȵ�����ֵ
enddo
call cwmax(wmax,w,m,n)
write(6,10) w
10 format(1x,6e16.8)
write(6,*) wmax
end

!********************************�ӳ��򲿷�**********************************
!�ӳ���1����������Ӷ�w��ϵ������k
subroutine ckw(k,l,p,r,s,z,m,n)
dimension k((m-1)*(n-1),(m-1)*(n-1)),l(m-1,n-1),p(m-1,n-1),z(m-1,n-1),r(m-1,n-1),s(m-1,n-1)
doubleprecision k,l,p,z,r,s
integer m,n,a,b,c,d
k=0.0
a=0
b=0
do j=2,n
a=b
 do i=2,m
 b=a+i-1
 c=1+(m-1)*(j-2)
 d=(j-1)*(m-1)
 if((j.eq.2).or.(j.eq.n)) then
  if((i.eq.2).or.(i.eq.m)) then
  k(b,b)=z(i-1,j-1)-2.0
  else
  k(b,b)=z(i-1,j-1)-1.0
  endif
 else
  if((i.eq.2).or.(i.eq.m)) then
  k(b,b)=z(i-1,j-1)-1.0
  else
  k(b,b)=z(i-1,j-1)
  endif
 endif
 if(b.gt.c) then
 k(b,b-1)=s(i-1,j-1)
 endif
 if(b.lt.d) then
 k(b,b+1)=s(i-1,j-1)
 endif
 if(b.gt.(c+1)) then
 k(b,b-2)=1.0
 endif
 if(b.lt.(d-1)) then
 k(b,b+2)=1.0
 endif
 if(j.lt.n) then
 k(b,b+m-1)=p(i-1,j-1)
  if(b.gt.c) then
  k(b,b+m-2)=r(i-1,j-1)
  endif
  if(b.lt.d) then
  k(b,b+m)=l(i-1,j-1)
  endif
 endif
 if(j.gt.2) then
 k(b,b-m+1)=p(i-1,j-1)
  if(b.gt.c) then
  k(b,b-m)=l(i-1,j-1)
  endif
  if(b.lt.d) then
  k(b,b-m+2)=r(i-1,j-1)
  endif
 endif
 if(j.lt.(n-1)) then
 k(b,b+2*m-2)=1.0
 endif
 if(j.gt.3) then
 k(b,b-2*m+2)=1.0
 endif
 enddo
enddo
end
!�ӳ���2���������fai��ϵ������k
subroutine ckfai(k,m,n)
dimension k((m-1)*(n-1),(m-1)*(n-1))
doubleprecision k
integer m,n,a,b,c,d
k=0.0
a=0
b=0
do j=2,n
a=b
 do i=2,m
 b=a+i-1
 c=1+(m-1)*(j-2)
 d=(j-1)*(m-1)
 if((j.eq.2).or.(j.eq.n)) then
  if((i.eq.2).or.(i.eq.m)) then
  k(b,b)=22.0
  else
  k(b,b)=21.0
  endif
 else
  if((i.eq.2).or.(i.eq.m)) then
  k(b,b)=21.0
  else
  k(b,b)=20.0
  endif
 endif
 if(b.gt.c) then
 k(b,b-1)=-8.0
 endif
 if(b.lt.d) then
 k(b,b+1)=-8.0
 endif
 if(b.gt.(c+1)) then
 k(b,b-2)=1.0
 endif
 if(b.lt.(d-1)) then
 k(b,b+2)=1.0
 endif
 if(j.lt.n) then
 k(b,b+m-1)=-8.0
  if(b.gt.c) then
  k(b,b+m-2)=2.0
  endif
  if(b.lt.d) then
  k(b,b+m)=2.0
  endif
 endif
 if(j.gt.2) then
 k(b,b-m+1)=-8.0
  if(b.gt.c) then
  k(b,b-m)=2.0
  endif
  if(b.lt.d) then
  k(b,b-m+2)=2.0
  endif
 endif
 if(j.lt.(n-1)) then
 k(b,b+2*m-2)=1.0
 endif
 if(j.gt.3) then
 k(b,b-2*m+2)=1.0
 endif
 enddo
enddo
end
!�ӳ���3����������Ӷ�w��ϵ���������ı仯�����ֵ
subroutine cxs(l,p,r,s,z,fai,t,D,m,n)
dimension l((m-1),(n-1)),p((m-1),(n-1)),z((m-1),(n-1)),r((m-1),(n-1)),s((m-1),(n-1)),fai(m+1,n+1)
doubleprecision l,p,r,s,z,fai,t,D,a1,a2,a3
integer m,n
do j=2,n
 do i=2,m
 a1=t*(fai(i+1,j)-2.0*fai(i,j)+fai(i-1,j))/D
 a2=t*(fai(i,j+1)-2.0*fai(i,j)+fai(i,j-1))/D
 a3=-t*(fai(i+1,j+1)+fai(i-1,j-1)-fai(i+1,j-1)-fai(i-1,j+1))/(2.0*D)
 l(i-1,j-1)=2.0-a3/4.0
 p(i-1,j-1)=-a1-8.0
 r(i-1,j-1)=2.0+a3/4.0
 s(i-1,j-1)=-a2-8.0
 z(i-1,j-1)=20.0+2.0*a1+2.0*a2
 enddo
enddo
end
!�ӳ���4����������Ӷ�w�����Է������ұ������ֵ��������
subroutine cb(b,w,m,n,E)
dimension b((m-1)*(n-1)),w(m+1,n+1)
integer m,n
doubleprecision b,w,a,c,d,p,f,E
a=0
c=0
do j=2,n
a=c
 do i=2,m
  c=a+i-1
  d=w(i+1,j+1)+w(i-1,j-1)-w(i+1,j-1)-w(i-1,j+1)
  p=w(i+1,j)-2.0*w(i,j)+w(i-1,j)
  f=w(i,j+1)-2.0*w(i,j)+w(i,j-1)
  b(c)=((d*d/16.0)-p*f)*E
 enddo
enddo
end
!�ӳ���5���������fai�����Է������ұ������ֵ��������
subroutine cbwfai(a,b,m,n)
dimension a(m+1,n+1),b((m-1)*(n-1))
integer m,n,c,d
doubleprecision a,b
c=0
d=0
do j=2,n
c=d
 do i=2,m
 d=c+i-1
 a(i,j)=b(d)
 enddo
enddo
end
!�ӳ���6������ǰ�����μ�������w������ֵ
subroutine cwcm(wcm,w,w1,m,n)
dimension w(m+1,n+1),w1(m+1,n+1)
integer m,n
doubleprecision wcm,a,w,w1
wcm=0
do i=2,m
 do j=2,n
  a=abs(w(i,j)-w1(i,j))
  if(wcm.lt.a) then
  wcm=a
  endif
 enddo
enddo
end
!�ӳ���7������w�е�����Ӷ�wmax
subroutine cwmax(wmax,w,m,n)
dimension w(m+1,n+1),w1(m+1,n+1)
integer m,n
doubleprecision wmax,w,w1
wmax=0
do i=2,m
 do j=2,n
  if(w(i,j).gt.wmax) then
  wmax=w(i,j)
  endif
 enddo
enddo
end

!�ӳ���:8����˹��ȥ����������Է����飩
subroutine gauss(a,b,n)
dimension a(n,n),b(n)
doubleprecision a,b
integer n
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
end