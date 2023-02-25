subroutine prefs(idx,outer)
	use fsbr
	implicit none
	integer,intent(in)		:: idx
	logical,intent(in)	  :: outer

	if(outer) then
		fs%ux =br%balno(idx)%out%ux
		fs%uh =br%balno(idx)%out%uh
		fs%uy =br%balno(idx)%out%uy
		fs%phx=br%balno(idx)%out%phx
		fs%phy=br%balno(idx)%out%phy
	else
		fs%ux =br%balno(idx)%inn%ux
		fs%uh =br%balno(idx)%inn%uh
		fs%uy =br%balno(idx)%inn%uy
		fs%phx=br%balno(idx)%inn%phx
		fs%phy=br%balno(idx)%inn%phy
	endif
endsubroutine prefs

subroutine postfs(idx,outer)
	use fsbr
	implicit none

	integer	,pointer :: nx,ny
	logical ,pointer :: simpmoment
	real(dp),pointer :: area(:)
	real(dp),pointer :: xc(:,:),yc(:,:),tx(:,:),ty(:,:)
	integer,intent(in)		:: idx
	logical,intent(in)	  :: outer
	real(dp)    :: tfx,tfy,mx,my
	real(dp)	:: temp,mzx,mzy,mz
	real(dp)	:: u
	integer		:: i,j

	nx			=>fs%nx
	ny			=>fs%ny
	simpmoment	=>br%simpmoment
	area		=>fs%area
	xc			=>fs%xc
	yc			=>fs%yc
	tx			=>fs%tx
	ty			=>fs%ty

	tfx	=0.0_dp
	tfy	=0.0_dp
	mx	=0.0_dp
	my	=0.0_dp
	mz	=0.0_dp
	mzx	=0.0_dp
	mzy	=0.0_dp
	
	if(outer) then
		u=br%balno(idx)%out%app/2.0_dp
	else
		u=br%balno(idx)%inn%app/2.0_dp
	endif

	do i=-ny,ny					!����j������y����
		if(i/=0) then				!����0������
			do j=nx,-nx,-1		!���ơ�����    !�˴��߽����Ҫ
				if(j/=0) then									!����0������
					!�����϶������϶����أ�����Ϊ�����弸�����ģ�����Ϊ�Ӵ����ʰ뾶
					temp=area(i)*tx(i,j)  
					tfx=tfx+temp
					mzx=mzx-temp*yc(i,j)	!mzx_f*b   mz=mzx_f*b+mzy_f*a
					temp=area(i)*ty(i,j)
					tfy=tfy+temp
					mzy=mzy+temp*xc(i,j)	!mzy_f*a		mz=mzx_f*b+mzy_f*a
!					mx=mx+y1*area(i)*rb   !���ٵ�����
!					my=my-(rb-yc(i,j)**2/2.0_dp/rb)*x1*area(i)  !y������ƴ���������������Ӱ��
				endif
			enddo
		endif
	enddo
	tfx=tfx/br%np			! longitudinal creep force!kalker simplied contact model pi/2 ��һ���������������ʱӦȥ���ñ任
	tfy=tfy/br%np			! lateral creep force !hertz�Ӵ� 2*pi/3

	if(simpmoment) then   !�ж��Ƿ���ü򻯷������������ڹ��������ĵ�����
		mx= tfy*(br%r-u)    !ע�����   ��������������������ϣ�������ͬ�ķ���ʹmx,my���ٻ���
		my=-tfx*(br%r-u)    !ע����� 1111111111111111111111111111111111
	endif

	fs%tfx=tfx
	fs%tfy=tfy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(br%sym) then
		if(outer) then
			temp=br%balno(1)%out%p*br%mu
			br%balno(:)%out%tx =tfx*temp
			br%balno(:)%out%ty =tfy*temp
			br%balno(:)%out%mx =mx*temp  !����
			br%balno(:)%out%my =my*temp
			mzx=mzx*temp*br%balno(1)%out%by
			mzy=mzy*temp*br%balno(1)%out%ax
			mz=mzx+mzy
			br%balno(:)%out%mz=mz
		else
			temp=br%balno(1)%inn%p*br%mu
			br%balno(:)%inn%tx =tfx*temp
			br%balno(:)%inn%ty =tfy*temp 
			br%balno(:)%inn%mx =mx*temp 
			br%balno(:)%inn%my =my*temp 
			mzx=			   mzx*temp*br%balno(1)%inn%by  !���� 
			mzy=			   mzy*temp*br%balno(1)%inn%ax  !����
			mz=mzx+mzy
			br%balno(:)%inn%mz=mz
		endif
	else
		if(outer) then
			br%balno(idx)%out%tx =tfx*br%balno(idx)%out%p*br%mu 
			br%balno(idx)%out%ty =tfy*br%balno(idx)%out%p*br%mu 
			br%balno(idx)%out%mx =mx*br%balno(idx)%out%p*br%mu 
			br%balno(idx)%out%my =my*br%balno(idx)%out%p*br%mu 
			mzx=			       mzx*br%balno(idx)%out%p*br%mu*br%balno(idx)%out%by
			mzy=			       mzy*br%balno(idx)%out%p*br%mu*br%balno(idx)%out%ax
			mz=mzx+mzy
			br%balno(idx)%out%mz=mz
		else
			br%balno(idx)%inn%tx =tfx*br%balno(idx)%inn%p*br%mu 
			br%balno(idx)%inn%ty =tfy*br%balno(idx)%inn%p*br%mu 
			br%balno(idx)%inn%mx =mx*br%balno(idx)%inn%p*br%mu 
			br%balno(idx)%inn%my =my*br%balno(idx)%inn%p*br%mu 
			mzx=			       mzx*br%balno(idx)%inn%p*br%mu*br%balno(idx)%inn%by
			mzy=			       mzy*br%balno(idx)%inn%p*br%mu*br%balno(idx)%inn%ax
			mz=mzx+mzy
			br%balno(idx)%inn%mz=mz
		endif
	endif
endsubroutine postfs