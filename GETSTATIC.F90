!������Ħ���������غɾ�̬�Ӵ���
subroutine getstatic()
	use fsbr
	use nr,only:ratint
	use intface,only:rtsafe,funcd,pow
	implicit none
	real(dp)		:: contang,load,appnd,dy
	real(dp)	:: pcho,pchi,t1,t2,omgb,omgo
	
	!���㴿�����غ������µĽӴ��ǺͽӴ��غ�,���ں�������
	contang=rtsafe(funcd,0.0_dp,pio2_d,br%kpst)
	load=br%aload/br%z/sin(contang)
	br%balno(:)%inn%alfst=contang
	br%balno(:)%out%alfst=contang
	br%balno(:)%inn%p=load
	br%balno(:)%out%p=load
	!��ֵ�����õ���ǰ�Ӵ��Ǵ��ⵯ����������
	call ratint(br%alf,br%ocom%appnd,contang,appnd,dy)
	br%balno(:)%e=0.5_dp*br%d_o-br%ro+(br%ao+appnd*load**(2.0_dp/3.0_dp))*cos(contang)
	br%dltast=(0.5_dp*br%d_i+br%ri-0.5_dp*br%d_o+br%ro)*(tan(contang)-tan(br%alf0))

	br%balno(:)%est=br%balno(1)%e

	call ratint(br%alf,br%ocom%appnd,contang,appnd,dy)
	pcho=contang-atan(sin(contang)/(cos(contang)+br%balno(1)%e/(br%r-.5_dp*appnd*load**(2.0_dp/3.0_dp))))
	br%s2rinn=-tan(pcho-contang)+sin(contang)/(br%balno(1)%e/(br%r-.5_dp*appnd*load**(2.0_dp/3.0_dp))-cos(contang))
	br%balno(:)%e0_octl=br%balno(1)%e
	br%dlta0_octl=br%dltast
	br%pchang_octl=pcho
	br%iconang_octl=contang
	br%oconang_octl=contang
	call ratint(br%alf,br%icom%appnd,contang,appnd,dy)
	br%omgb_octl=br%omgi*(br%balno(1)%e/(br%r-.5_dp*appnd*load**(2.0_dp/3.0_dp))-cos(contang))/cos(pcho-contang)
	br%omgo_octl=-br%omgb_octl*sin(pcho-contang)/sin(contang)

	call ratint(br%alf,br%icom%appnd,contang,appnd,dy)
	pchi=contang-atan(sin(contang)/(cos(contang)-br%balno(1)%e/(br%r-.5_dp*appnd*load**(2.0_dp/3.0_dp))))
	br%s2rout=tan(pchi-contang)+sin(contang)/(br%balno(1)%e/(br%r-.5_dp*appnd*load**(2.0_dp/3.0_dp))+cos(contang))
	br%balno(:)%e0_ictl=br%balno(1)%e
	br%dlta0_ictl=br%dltast
	br%pchang_ictl=pchi
	br%iconang_ictl=contang
	br%oconang_ictl=contang
	br%omgb_ictl=br%omgi*sin(contang)/sin(pchi-contang)
	call ratint(br%alf,br%ocom%appnd,contang,appnd,dy)
	br%omgo_ictl=br%omgb_ictl*cos(pchi-contang)/(br%balno(1)%e/(br%r-.5_dp*appnd*load**(2.0_dp/3.0_dp))+cos(contang))
endsubroutine getstatic