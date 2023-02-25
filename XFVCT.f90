subroutine x0vct(x)
	use fsbr
	use nr
	use intface,only:pow,rad2deg,deg2rad
	implicit none
	real(dp), dimension(:) :: x
	real(dp) ::est,uist,uost
	real(dp)::pchi,pcho,omgb,omgo,omgc
	real(dp)::t1,t2
	real(dp)::pchang
	real(dp)::e,rdi,rdo
	real(dp)::alfst,appndi,appndo
	real(dp)::dy
	integer	 ::n,i

	n=size(x)

	!��brĳЩ�ڼ�������иı��һ�Ӱ�쵽��һ�ּ���ı������㣻
	!�˵���ڲ��ö����㷨�໥�ݽ��ĵ��÷ǳ���Ҫ��
!	call zerobr()  !�ô����ȡ��ǳ��ؼ���
	alfst=br%balno(1)%inn%alfst
	if(all(x==0.0_dp)) then
		!mzk ���Ƽ���
		if(br%absolute) then	  	
			!һ�㹤��
			!-3:mgx+mgz;-2:mgx;-1:mgz;0:no inert;1:fc only;2:fc+mgx;3:fc+mgz;4:fc+mgx+mgz
			if(br%inert>=1) then	
				!һ�㹤�������ǹ�����
				x(5)=br%balno(1)%inn%alfst
				x(6)=br%balno(1)%out%alfst
				call conangs_sym(x(5),x(6))
				omgc=br%balno(1)%omgc
				omgb=br%balno(1)%omgb
				pchang=br%pchang
				x(1)=omgc
				x(2)=-0.001_dp
				x(3)=omgb*cos(pchang)
				x(4)=omgb*sin(pchang)  
			else	
				!һ�㹤���������ǹ�����
				!1:mzk ���Ʒ�������ã�2���ڹ�����������,��Ȧ������3��������������ۣ���Ȧ����
				call ratint(br%alf,br%icom%appnd,alfst,appndi,dy)
				call ratint(br%alf,br%ocom%appnd,alfst,appndo,dy)
				uist=appndi*pow(br%balno(1)%inn%p,2.0_dp/3.0_dp)/2.0_dp
				uost=appndo*pow(br%balno(1)%out%p,2.0_dp/3.0_dp)/2.0_dp
				est=abs(br%icom%rx2)+abs(br%icom%ry2)-(br%ai+uist*2.0_dp)*cos(alfst)
				if(br%rcm==1) then
					omgc=((est-(br%r-uist)*cos(alfst))*br%omgi+(est+(br%r-uost)*cos(alfst))*br%omgo)/2.0_dp/est
					omgb=(omgc-br%omgi)*(est/(br%r-uist)-cos(alfst))
					pchang=alfst
				elseif(br%rcm==2) then
					!�ڹ�������!pchi �ڹ�������ʱ���������̬��
					pchi=atan(est*sin(alfst)/(est*cos(alfst)-(br%r-uist)))
					omgc=(est*(br%omgo+br%omgi)+(br%r-uist)*(br%omgo-br%omgi)*cos(alfst))/2.0_dp/est
					omgb=(omgc-br%omgi)*sin(alfst)/(sin(pchi)*cos(alfst)-cos(pchi)*sin(alfst))
					pchang=pchi
				elseif(br%rcm==3) then
					!���������
					!pcho ���������ʱ���������̬��
					pcho=atan(est*sin(alfst)/(est*cos(alfst)+(br%r-uost)))
					omgc=(est*(br%omgo+br%omgi)+(br%r-uist)*(br%omgo-br%omgi)*cos(alfst))/2.0_dp/est
					omgb=(omgc-br%omgo)*sin(alfst)/(sin(pcho)*cos(alfst)-cos(pcho)*sin(alfst))
					pchang=pcho
				endif
				x(1)=omgc
				x(2)=0.001_dp
				x(3)=omgb*cos(pchang)	
				x(4)=omgb*sin(pchang)
				x(5)=alfst
				x(6)=alfst
			endif
		else  
			!���⹤��
			!-3:mgx+mgz;-2:mgx;-1:mgz;0:no inert;1:fc only;2:fc+mgx;3:fc+mgz;4:fc+mgx+mgz
			if(br%inert>=1) then	
				!���´������ڼ���������
				i=br%rcm  !����
				!��������ƣ��ڽӴ�������
				br%rcm=3  !����Ϊ���
				x(5)=br%balno(1)%inn%alfst
				x(6)=br%balno(1)%out%alfst
				call conangs_sym(x(5),x(6))
				pcho=br%pchang	!����ʽ�ǳ��ӽ�
				omgb=br%balno(1)%omgb
!				t1=omgb*(cos(pcho)*sin(x(5))-sin(pcho)*cos(x(5)))+br%omgi*sin(x(5))
!				t2=omgb*(cos(pcho)*cos(x(5))+sin(pcho)*sin(x(5)))
!				br%s2rinn=t1/t2
				br%s2rinn=-tan(pcho-x(5))+sin(x(5))/(br%balno(1)%e/(br%r-br%balno(1)%inn%app*.5_dp)-cos(x(5)))
				br%balno(:)%e0_octl=br%balno(1)%e
				br%omgo_octl=br%omgo
				br%omgb_octl=omgb
				br%pchang_octl=pcho
				br%iconang_octl=x(5)
				br%oconang_octl=x(6)
				br%dlta0_octl=br%dlta0

				!�ڹ������ƣ���Ӵ�������
				br%rcm=2  !����Ϊ�ڿ�
				call conangs_sym(x(5),x(6))
				PCHI=br%pchang
				omgb=br%balno(1)%omgb
				omgo=br%omgo
!				t1=omgb*(sin(pchi)*cos(x(6))-cos(pchi)*sin(x(6)))+omgo*sin(x(6))
!				t2=omgb*(cos(pchi)*cos(x(6))+sin(pchi)*sin(x(6)))
!				br%s2rout=t1/t2
				br%s2rout=tan(pchi-x(6))+sin(x(6))/(br%balno(1)%e/(br%r-br%balno(1)%out%app*.5_dp)+cos(x(6)))
				br%balno(:)%e0_ictl=br%balno(1)%e
				br%omgo_ictl=br%omgo
				br%omgb_ictl=omgb
				br%pchang_ictl=pchi
				br%iconang_ictl=x(5)
				br%oconang_ictl=x(6)
				br%dlta0_ictl=br%dlta0

				!�ظ���ʼ״̬��������������
				br%rcm=i
				!���⹤�������ǹ�����
				x(5)=br%balno(1)%inn%alfst
				x(6)=br%balno(1)%out%alfst
				call conangs_sym(x(5),x(6))
				omgo=br%omgo
				omgb=br%balno(1)%omgb
				pchang=br%pchang
				x(1)=omgo
				x(2)=-0.001_dp
				x(3)=omgb*cos(pchang)
				x(4)=omgb*sin(pchang)  
				if(n>6) then
					x(7)=br%balno(1)%inn%p
					x(8)=br%balno(1)%out%p
				endif
			else  	
				!���⹤���������ǹ�����
				!1:mzk ���Ʒ�������ã�2���ڹ�����������,��Ȧ������3��������������ۣ���Ȧ����
				if(br%rcm==1) then
					uist=appndi*pow(br%balno(1)%inn%p,2.0_dp/3.0_dp)/2.0_dp
					uost=appndo*pow(br%balno(1)%out%p,2.0_dp/3.0_dp)/2.0_dp
					est=abs(br%icom%rx2)+abs(br%icom%ry2)-(br%ai+uist*2.0_dp)*cos(alfst)
					omgb=br%omgi*(est-(br%r-uist)*cos(alfst))/(br%r-uist)
					omgo=omgb*(br%r-uost)/(est+(br%r-uost)*cos(alfst))
					pchang=alfst
				elseif(br%rcm==2) then
					!�ڹ�������
					!pchi �ڹ�������ʱ���������̬��
					!���㹫ʽ�ο�harris������ע���������е��������ı˴˴���
					pchi=atan(br%dm/2.0*sin(br%balno(1)%inn%alfst)/(br%dm/2.0*cos(br%balno(1)%inn%alfst)-br%r))
					omgb=br%omgi*(-br%dm/2.0+br%r*cos(br%balno(1)%inn%alfst))/	&
											br%r/(cos(pchi)*cos(br%balno(1)%inn%alfst)+sin(pchi)*sin(br%balno(1)%inn%alfst))
					omgo=omgb/((br%dm/2.0+br%r*cos(br%balno(1)%inn%alfst))/	&
											br%r/(cos(br%balno(1)%inn%alfst)*cos(pchi)+sin(pchi)*sin(br%balno(1)%inn%alfst)))
					!�ڹ���������Ȧ���ٶ�
					omgo=abs(omgo)
					omgb=abs(omgb)
					pchang=pchi
				elseif(br%rcm==3) then
					!���������
					!pcho ���������ʱ���������̬��
					!���㹫ʽ�ο�harris������ע���������е��������ı˴˴���
					pcho=atan(br%dm/2.0*sin(br%balno(1)%inn%alfst)/(br%dm/2.0*cos(br%balno(1)%inn%alfst)+br%r))
					omgb=br%omgi*(-br%dm/2.0+br%r*cos(br%balno(1)%inn%alfst))/	&
											br%r/(cos(pcho)*cos(br%balno(1)%inn%alfst)+sin(pcho)*sin(br%balno(1)%inn%alfst))
					omgo=omgb/((br%dm/2.0+br%r*cos(br%balno(1)%inn%alfst))/	&
											br%r/(cos(br%balno(1)%inn%alfst)*cos(pcho)+sin(pcho)*sin(br%balno(1)%inn%alfst)))
					!�����������Ȧ���ٶ�'
					omgo=abs(omgo)
					omgb=abs(omgb)
					pchang=pcho
				endif
				x(1)=omgo
				x(2)=0.001_dp
				x(3)=omgb*cos(pchang)
				x(4)=omgb*sin(pchang)
				x(5)=alfst
				x(6)=alfst
				if(n>6) then
					x(7)=br%balno(1)%inn%p
					x(8)=br%balno(1)%out%p
				endif
			endif
		endif

		!���������⾲��ѧ�����������Ա����ã�
		br%x0(:)=x(:)
		br%x0(2)=0.0_dp	!�������۷�����ƫ���ǲ����ڣ�
	else
		x(:)=br%x(:)
	endif
end subroutine x0vct

subroutine printx0(x)
	use fsbr
	use intface,only:rad2deg
	implicit none
	integer :: i
	real(dp), dimension(br%n) :: x

	if(pntxo) then
		write(*,*) '    *************************��     ֵ***********************'
		write(*,*) '        index','          x'

		do i=1,size(x)
			if(i==5 .or. i==6) then
				write(*,*) i,x(i),'rad',rad2deg(x(i)),'deg'
			else
				write(*,*) i,x(i)
			endif
		enddo

		write(*,*) '    *********************************************************'
		pntxo=.false.  !ֻ�����һ������㷨�ĳ�ֵ��
	endif
endsubroutine printx0

subroutine printxf(x,f)
	use fsbr
	implicit none
	integer :: i
	real(dp), dimension(br%n) :: x
	real(dp), dimension(br%m) :: f
	real(dp) DNRM2

	write(*,*) '        index','          x','                     f'
	do i=1,size(x)
		write(*,*) i,x(i),f(i)
	enddo
	write(*,*) 'norm of F(x):',DNRM2(size(x),f,1)
	write(*,*) 'omgb',sqrt(x(3)**2+x(4)**2)
endsubroutine printxf

subroutine savexf(x,f)
	use fsbr
	implicit none
	integer :: i,j
	real(dp), dimension(br%n) :: x,f

	do i=1,br%n
		br%x(i)		=x(i)
		do j=1,br%z
			br%balno(j)%f(i)		=f(i)
		enddo
	enddo
endsubroutine savexf
