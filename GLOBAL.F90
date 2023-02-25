!��Ȼ����jj kalker�����񻮷ַ��������ڼ�����ԣ����ַ����Լ���Ҫ�������
!��������ı׶����ڱ߽�Ĵ�����Ϊ��������Ҫ���и���ļ��㡣
!�߽��Ĵ����Ƚ����⣬�������е�������
module fsbr
	use nrtype
	implicit none
	logical			:: pntxo
	logical			:: maxit
	logical			:: embed
	real(dp)		:: mprec,sqmprec
	type:: fastdata    !fastsim�㷨���ݽṹ
 		integer		:: nx,ny		!x,y�������ȷֶ�����nx=ny
		real(dp)	:: hy		!y��������߳�
		real(dp)	:: ux,uy,phx,phy   !�������们��������kalker���ٷ�������
		real(dp)	:: uh   !�������们��������kalker���ٷ�������
		real(dp)	:: spole(1,2)	!��ת����  
		logical		:: hertz						!hertzѹ���ֲ������֣�true:hertz��false:non-hertz
		logical		:: heathcote				!heathcote slip
		real(dp),allocatable:: hx(:)										!x��������߳�
		real(dp),allocatable:: area(:)									!��Ԫ���
		real(dp),allocatable:: xc(:,:),yc(:,:)					!����ѹ��������
		real(dp),allocatable:: tx(:,:),ty(:,:)	!�������϶�������
		real(dp),allocatable:: sx(:,:),sy(:,:)	!��Ի����ٶȣ������� relative velocity
		real(dp)	:: tfx,tfy
	end type fastdata

	!����Ϊ���ϵͳ�������塢��������ݽṹ
	type:: patchabk
		real(dp)	:: rx1,ry1,rx2,ry2
		real(dp)	:: rx,ry
		real(dp)	:: and,bnd,appnd !kΪ�Ӵ��ն�ϵ�� and,bnd Ϊ�������������İ��᳤�ȣ���������������
		real(dp)	:: ei1,ei2
		real(dp)	:: ke   !��Բ��
		real(dp)	:: kp   !��Բ�ʼ��㾫��
		real(dp)	:: e1
		real(dp)	:: e2
		real(dp)	:: nu1
		real(dp)	:: nu2
		logical		:: xlarger				!
	end type patchabk

	type comm  !���в���
		real(dp),allocatable	:: and(:),bnd(:),appnd(:)  !��ѹ�����ᡢ������
		real(dp),allocatable	:: ke(:)
		real(dp),allocatable	:: ei1(:),ei2(:)
		real(dp),allocatable	:: rx(:),ry(:)	!�ۺ����ʰ뾶�����ʺ�
		real(dp),allocatable	:: k(:)	!�Ӵ��ն�ϵ��
		real(dp)	:: rx2,ry2			!�������ʰ뾶
		logical		:: xlarger				!
		real(dp)	:: c11,c22,c23,c32,c33 !kalkerϵ��
		REAL(DP)	:: RC		!Common Radius of Curvature 
	endtype comm

	type contact
		real(dp)	:: ax,by,app,ab ! ���ᡢ������
		real(dp)	:: rdmd ! �Ӵ����ı��κ�뾶
		real(dp)	:: p		!�Ӵ�����ѹ��
		real(dp)	:: alfst !��̬�Ӵ���contact angle �����������������
		real(dp)	:: alf !contact angle
		real(dp)	:: xix,xiy,psi	!�们��
		real(dp)	:: ux,uh,uy,phx,phy   !�������们��������kalker���ٷ�������
		real(dp)	:: tx,ty !�������ܵ��ļ����϶���
		real(dp)	:: mx,my,mz !�������ܵ�������
		real(dp)	:: rd			!????????
		real(dp)	:: rbd		!��������κ������ٰ뾶
		real(dp)	:: spn2rol  !spin to roll ratio
	endtype contact

	type roller
		type(contact) :: inn,out
		real(dp)	:: e	!�������Բ�뾶
		real(dp)	:: e0_ictl,e0_octl	!�ڡ��������Ħ�������룩״̬�������Բ�뾶
		real(dp)	:: est	!��ƽ����Ħ��״̬�������Բ�뾶
		real(dp)	:: mgx,mgz !�������ܵ�����������
		real(dp)	:: fc	!�������ܵ���������
		real(dp)	:: posang !
		real(dp)	:: omgc		!ÿ�������幫ת���ٶ�
		real(dp)	:: omgb,omgx,omgy,omgz		!��������ת���ٶ�
		real(dp),allocatable	:: f(:)  !��������������
	endtype roller
	
	type:: bearsysdata
		type(comm)	:: icom,ocom
		logical		:: simpmoment			!�򻯼���
 		logical		:: sym			!�Ƿ�Գ��غ�
		logical		:: hertz						!hertzѹ���ֲ������֣�true:hertz��false:non-hertz
		logical		:: OPTIN
		logical		:: unscl	!unscale coordinates? TRUE:UNSCALE;FALSE:NO
		logical		:: rcm		!raceway control method
		logical		:: absolute		!�����˶�,���н��ٶȾ�����ڹ�������ϵ����
		logical		:: cmdcreep	!�����غ������µ��们ģ��
		integer		:: s2rmthd	!�����ȼ��㷽��
		integer		:: hcmthd  !heathcote���㷽��
		integer		:: z				!���������
		integer		:: cyc
		integer		:: inert	!0:no inert;1:fc only;2:fc+mgx;3:fc+mgz;4:fc+mgx+mgz
		real(dp),allocatable	:: kn(:)  !���ϸն�ϵ��������ƽ��ʱ
		real(dp)	:: d,r			!diameter of ballsֱ�����뾶
		real(dp)	:: ri,ro		!inner and outer raceway groove curvature radius
		real(dp)	:: d_i,d_o	!inner and outer raceway diameter
		real(dp)	:: alf0			!free contact angle
		real(dp)	:: a,ai,ao  !distance between raceway groove curvature centers
		real(dp)	:: rx1,ry1	!���������ʰ뾶
		real(dp)	:: dm				!mean diameter
		real(dp)	:: fo,fi		!raceway curvature/ball dia
		real(dp)	:: pd				!bearing diametral clearance
		real(dp)	:: kp,kpst	!��Բ�ʡ���̬�Ӵ��Ǽ��㾫��
		real(dp)	:: e1,g1		!�����嵯�ԡ�����ģ��
		real(dp)	:: e2,g2		!��Ȧ���ԡ�����ģ��
		real(dp)	:: m1				!����������
		real(dp)	:: ip				!������ת������
		real(dp)	:: mu				!Ħ��ϵ��
		real(dp)	:: nu1			!�����嵯��ģ��
		real(dp)	:: nu2			!�����嵯��ģ��
		real(dp)	:: rho1			!�ܶ�
		real(dp)	:: rho2			!�ܶ�
		real(dp)	:: np		!n prime  ���kalker fastsim�㷨����

		real(dp)	:: aload		!��������غ�
		real(dp)	:: rload		!��о����غ�
		real(dp)	:: mload		!��о����غ�
		real(dp)	:: xtran=0.0_dp	!ƽ��x���꣬����tecplot����ֵΪ��Բ����ı�����
		type(roller),allocatable :: balno(:)
		real(dp),allocatable	:: x(:),x0(:)
		real(dp)	:: dlta,dltr,theta
		real(dp)	:: f3(3)
		real(dp)	:: dltast,dlta0 !��ƽ�ⴿ����λ�ƣ�������Ħ���⾲��ѧ����λ��
		real(dp)	:: dlta0_ictl,dlta0_octl !�ڡ���ض�ƽ������λ��
		real(dp)	:: s2rinn	!��������ƣ��ڹ���������
		real(dp)	:: s2rout	!�ڹ������ƣ������������
		real(dp)	:: omgi		!��Ȧת��
		real(dp)	:: omgo		!��Ȧת��
		real(dp)	:: omgo_ictl,omgo_octl !�������Ȧת��
		real(dp)	:: omgc		!���ּܽ��ٶ�
		real(dp)	:: omgb_ictl,omgb_octl	!������ת�ٳ�ֵ
		real(dp)	:: iconang_ictl,oconang_ictl !�ڿؽӴ���
		real(dp)	:: iconang_octl,oconang_octl !�ڿؽӴ���
		real(dp)	:: startang !��ʼ������ķ�λ��
		real(dp)	:: pchang !�����帩����
		real(dp)	:: pchang_octl	!��������Ƹ�����
		real(dp)	:: pchang_ictl	!�ڹ������Ƹ�����

		real(dp)	:: aj,gama	!���ݲ���
		real(dp)	:: tmp1,tmp2	!���ݲ���

		integer		:: m=6,n=6  !���̡�δ֪���ĸ���
		
		integer		:: nsec
		real(dp),allocatable	:: alf(:)
		
		integer		::	index=0
		integer		::	intmthd	!1:�Ӵ���->������;2:������->�Ӵ��� ����Ħ����
		integer		::	crpmthd	!1:�Ӵ���->������;2:������->�Ӵ��ǡ�����Ħ����
		logical		::	init
	end type bearsysdata

	type(fastdata)	 ,target	:: fs
	type(bearsysdata),target	:: br
end module fsbr

subroutine freefs()
	use fsbr
	implicit none
	logical				::saveslip

	deallocate(fs%hx,fs%area,fs%xc,fs%yc,fs%tx,fs%ty,fs%sx,fs%sy)
	deallocate(br%balno)
end subroutine freefs

subroutine freemem()
	use fsbr
	implicit none
	logical				::saveslip

	deallocate(fs%hx,fs%area,fs%xc,fs%yc,fs%tx,fs%ty,fs%sx,fs%sy)
	deallocate(br%balno)
end subroutine freemem