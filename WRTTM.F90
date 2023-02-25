SUBROUTINE WRTTM(idx)
	USE FSBR
	USE INTFACE,ONLY:RAD2DEG,POW
	IMPLICIT NONE
	integer	::idx
	
	WRITE(*,*) "���ɽӴ���:",BR%ALf0 ,RAD2DEG(BR%ALf0)
	WRITE(*,*) "ʵ�ʽӴ���:",RAD2DEG(BR%balno(idx)%inn%alf),RAD2DEG(BR%balno(idx)%out%alf)
	WRITE(*,*) "��ƽ��Ӵ���:",BR%BALNO(idx)%INN%ALfST,RAD2DEG(BR%BALNO(idx)%INN%ALfST)
	!WRITE(*,*) "�޻�����Ȧ���ٶ�:",OMGO,BR%OMGI*(1.0_DP-GAMA)/(1.0_DP+GAMA)
	!WRITE(*,*) "�޻�����������ٶ�:",OMGB,BR%OMGI*BR%DM/BR%D*(1.0_DP-GAMA)

	WRITE(*,*) "TX=",BR%BALNO(idx)%OUT%TX
	WRITE(*,*) "TY=",BR%BALNO(idx)%OUT%TY
	WRITE(*,*) "PO=",BR%BALNO(idx)%OUT%P
	WRITE(*,*) "MX=",BR%BALNO(idx)%OUT%MX
	WRITE(*,*) "MY=",BR%BALNO(idx)%OUT%MY
	WRITE(*,*) "MZ=",BR%BALNO(idx)%OUT%MZ
	WRITE(*,*) "�϶�������:",SQRT(BR%BALNO(idx)%OUT%TX**2+BR%BALNO(idx)%OUT%TY**2)
	WRITE(*,*) "��Ӵ�����COULOMB��:",BR%MU*BR%BALNO(idx)%OUT%P
!	write(*,*) '��Բ��:',	BR%OCOM%KE,	SQRT(BR%OCOM%KE),1.0_DP/BR%OCOM%KE
	WRITE(*,*) 
	WRITE(*,*) "TX=",BR%BALNO(idx)%INN%TX
	WRITE(*,*) "TY=",BR%BALNO(idx)%INN%TY
	WRITE(*,*) "PI=",BR%BALNO(idx)%INN%P
	WRITE(*,*) "MX=",BR%BALNO(idx)%INN%MX
	WRITE(*,*) "MY=",BR%BALNO(idx)%INN%MY
	WRITE(*,*) "MZ=",BR%BALNO(idx)%INN%MZ
	write(*,*) "�϶�������:",sqrt(BR%BALNO(idx)%INN%TX**2+BR%BALNO(idx)%INN%TY**2)
	WRITE(*,*) "�ڽӴ�����COULOMB��:", BR%MU*BR%BALNO(idx)%INN%P
!	write(*,*) '��Բ��:',	BR%iCOM%KE,	SQRT(BR%iCOM%KE),1.0_DP/BR%iCOM%KE
	WRITE(*,*) 
	write(*,*) "������������:",br%balno(idx)%fc
	write(*,*) "������x���������ط���:",br%balno(idx)%mgx
	write(*,*) "������z���������ط���:",br%balno(idx)%mgz
ENDSUBROUTINE WRTTM

SUBROUTINE OPTSLU()
	USE FSBR
	USE INTFACE,ONLY:RAD2DEG,rad2rev
	IMPLICIT NONE
	REAL(DP)::PCHI,PCHO,omgb,omgo
	real(dp)::t1,t2
	character(len=1)			:: tab=char(9)
	integer	:: i
	real(dp) DNRM2
	
	!������ļ���������ȡ����
	write(12,*) '�����JONES�ļ��㷽���õ��ĳ���'
	write(12,*) BR%x0(1)*60_dp/(2.0_dp*PI_D),tab,"OMGO RPM"
	write(12,*) BR%x0(2)*60_dp/(2.0_dp*PI_D),tab,"OMGX RPM"
	write(12,*) BR%x0(3)*60_dp/(2.0_dp*PI_D),tab,'OMGY RPM'
	write(12,*) BR%x0(4)*60_dp/(2.0_dp*PI_D),tab,'OMGZ RPM'
	write(12,*) RAD2DEG(BR%x0(5)),tab,"ALPHAI"
	write(12,*) RAD2DEG(BR%x0(6)),tab,'ALPHAO'
	write(12,*) sqrt(BR%x0(3)**2+BR%x0(4)**2)*60_dp/(2.0_dp*PI_D),tab,'OMGB  RPM��������ת���ٶ�'	!
	write(12,*) RAD2DEG(BR%BALNO(1)%INN%ALfST),tab, 'ALPHA DEGREE'
	write(12,*) 0.0_dp,tab, 'beta DEGREE'
	write(12,*)	'������JONES�����Ƽ���'
	write(12,*)	RAD2DEG(br%pchang_octl),tab,'��������Ƹ�����HARRIS,��һ�����Ƽ���'
	write(12,*)	RAD2DEG(br%pchang_ictl),tab,'�ڹ������Ƹ�����HARRIS,��һ�����Ƽ���'
	write(12,*) br%omgb_octl,tab,'��������ƹ�����ת��RPM'
	write(12,*) br%omgo_octl,tab,'�����������Ȧת��RPM'
	write(12,*) br%s2rinn,tab,'��������ƣ��ڹ���������'
	write(12,*) br%omgb_ictl,tab,'�ڹ������ƹ�����ת��RPM'
	write(12,*) br%omgo_ictl,tab,'�ڹ���������Ȧת��RPM'
	write(12,*) br%s2rout,tab,'�ڹ������ƣ������������'

	write(12,*) '������̵����ս�'
	write(12,*) rad2rev(BR%x(2)),tab,'OMGX  RPM X�����(����ڱ��ּ�)'
	write(12,*) rad2rev(BR%x(3)),tab,'OMGY  RPM Y�����(����ڱ��ּ�)'
	write(12,*) rad2rev(BR%x(4)),tab,'OMGZ  RPM Z�����(����ڱ��ּ�)'
write(12,*) br%mu,tab,'br%mu'
write(12,*) br%aload,tab,'br%aload'
	if(br%absolute) then
		write(12,*) rad2rev(BR%x(1)),tab,'OMGC  RPM ���ּܽ��ٶ�'
	else
		write(12,*) rad2rev(BR%x(1)),tab,'OMGO  RPM ��Ȧ���ٶ�'
	endif
	write(12,*) RAD2DEG(BR%x(5)),tab,'ALPHAI DEGREE �ڽӴ���'
	write(12,*) RAD2DEG(BR%x(6)),tab,'ALPHAO DEGREE ��Ӵ���'
	write(12,*) RAD2DEG(atan(abs(BR%x(4))/SQRT(BR%x(2)**2+BR%x(3)**2))),tab,'������ALPHA DEGREE'
	write(12,*) RAD2DEG(atan(BR%x(2)/abs(BR%x(3)))),tab,'ƫ����BETA DEGREE'
	write(12,*) sqrt(rad2rev(BR%x(2))**2+rad2rev(BR%x(3))**2+rad2rev(BR%x(4))**2),tab, '��������ת���ٶ�'
	write(12,*) BR%BALNO(1)%out%spn2rol,tab,'��Ӵ�������'
	write(12,*) BR%BALNO(1)%inn%spn2rol,tab,'�ڽӴ�������'
	WRITE(12,*) BR%BALNO(1)%OUT%TX,	tab,"TX"
	WRITE(12,*) BR%BALNO(1)%OUT%TY,	tab,"TY"
	WRITE(12,*) BR%BALNO(1)%OUT%P,	tab,"PO"
	WRITE(12,*) BR%BALNO(1)%OUT%MX,	tab,"MX"
	WRITE(12,*) BR%BALNO(1)%OUT%MY,	tab,"MY"
	WRITE(12,*) BR%BALNO(1)%OUT%MZ,	tab,"MZ"
	WRITE(12,*) SQRT(BR%BALNO(1)%OUT%TX**2+BR%BALNO(1)%OUT%TY**2),tab,"�϶�������"
	WRITE(12,*) BR%BALNO(1)%INN%TX,	tab,	"TX"
	WRITE(12,*) BR%BALNO(1)%INN%TY,	tab,	"TY"
	WRITE(12,*) BR%BALNO(1)%INN%P,	tab,	"PI"
	WRITE(12,*) BR%BALNO(1)%INN%MX,	tab,	"MX"
	WRITE(12,*) BR%BALNO(1)%INN%MY,	tab,	"MY"
	WRITE(12,*) BR%BALNO(1)%INN%MZ,	tab,	"MZ"
	write(12,*) sqrt(BR%BALNO(1)%INN%TX**2+BR%BALNO(1)%INN%TY**2),tab,"�϶�������"
	write(12,*) '******************************************'
	WRITE(12,*) RAD2DEG(BR%ALf0),tab,"���ɽӴ���"
	WRITE(12,*) RAD2DEG(BR%X(5)),tab,"ʵ���ڽӴ���"
	WRITE(12,*) RAD2DEG(BR%X(6)),tab,"ʵ����Ӵ���"
	WRITE(12,*) RAD2DEG(BR%BALNO(1)%INN%ALfST),tab,"��ƽ��Ӵ���"
	!WRITE(*,*) "�޻�����Ȧ���ٶ�:",OMGO,BR%OMGI*(1.0_DP-GAMA)/(1.0_DP+GAMA)
	!WRITE(*,*) "�޻�����������ٶ�:",OMGB,BR%OMGI*BR%DM/BR%D*(1.0_DP-GAMA)
	WRITE(12,*) BR%MU*BR%BALNO(1)%OUT%P,tab,"��Ӵ�����COULOMB��"
!	write(12,*) '��Բ��:',	BR%OCOM%KE,	SQRT(BR%OCOM%KE),1.0_DP/BR%OCOM%KE
	WRITE(12,*) 
	WRITE(12,*) BR%MU*BR%BALNO(1)%INN%P,tab, "�ڽӴ�����COULOMB��"
!	write(12,*) '��Բ��:',	BR%iCOM%KE,	SQRT(BR%iCOM%KE),1.0_DP/BR%iCOM%KE
	WRITE(12,*) 
	write(12,*) br%balno(1)%fc,	tab,"������������"
	write(12,*) br%balno(1)%mgx,tab,"������x���������ط���"
	write(12,*) br%balno(1)%mgz,tab,"������z���������ط���"
	write(12,*) 
	do i=1,br%n
		write(12,"(e23.15,a1,a5,i1,a1)") br%balno(1)%f(i),tab,'br%f(',i,')'
	enddo
	write(12,*) DNRM2(size(br%balno(1)%f),br%balno(1)%f,1),tab,'norm of F(x)'

ENDSUBROUTINE OPTSLU

SUBROUTINE OPTLOG()	!�����������Ϣ
	USE FSBR
	IMPLICIT NONE
	character(len=1)			:: tab=char(9)
	
	!������������Ϣ���������ּ������
	write(12,*) '**************������������*********************'
	write(12,*)  pntxo,					tab,'pntxo'
	write(12,*)  maxit,					tab,'maxit'
	write(12,*)  embed,					tab,'embed'
	write(12,*)  mprec,					tab,'mprec'
	write(12,*)	 sqmprec,				tab,'sqmprec'
	write(12,*)  FS% HEATHCOTE,	tab,'FS%HEATHCOTE'
	write(12,*)  FS% HERTZ	,		tab,'FS%HERTZ'
	write(12,*)  FS%	NY,				tab,'FS%NY'
	write(12,*)  FS%	NX,				tab,'FS%NX'
	write(12,*)  br%simpmoment,	tab,'br%simpmoment'
 	write(12,*)  br%sym,				tab,'br%sym'
	write(12,*)  br%hertz,			tab,'br%hertz'						!hertzѹ���ֲ������֣�true:hertz��false:non-hertz
	write(12,*)  br%OPTIN,			tab,'br%OPTIN'
	write(12,*)  br%unscl,			tab,'br%unscl'	!unscale coordinates? TRUE:UNSCALE;FALSE:NO
	write(12,*)  br%rcm,				tab,'br%rcm'		!raceway control method
	write(12,*)  br%absolute,		tab,'br%absolute'		!�����˶�,���н��ٶȾ�����ڹ�������ϵ����
	write(12,*)  br%cmdcreep,		tab,'br%cmdcreep'	!�����غ������µ��们ģ��
	write(12,*)  br%hcmthd ,		tab,'br%hcmthd' !heathcote���㷽��
	write(12,*)  br%z		,				tab,'br%z'		!���������
	write(12,*)  br%cyc,				tab,'br%cyc'
	write(12,*)  br%inert,			tab,'br%inert'	!0:no inert;1:fc only;2:fc+mgx;3:fc+mgz;4:fc+mgx+mgz
	write(12,*)  br%d,					tab,'br%d'
	write(12,*)  br%r		,				tab,'br%r'	!diameter of ballsֱ�����뾶
	write(12,*)  br%ri,					tab,'br%ri'
	write(12,*)  br%ro	,				tab,'br%ro'	!inner and outer raceway groove curvature radius
	write(12,*)  br%d_i,				tab,'br%d_i'
	write(12,*)  br%d_o	,				tab,'br%d_o'!inner and outer raceway diameter
	write(12,*)  br%kp,					tab,'br%kp'
	write(12,*)  br%kpst,				tab,'br%kpst'	!��Բ�ʡ���̬�Ӵ��Ǽ��㾫��
	write(12,*)  br%e1,					tab,'br%e1'
	write(12,*)  br%g1	,				tab,'br%g1'	!�����嵯�ԡ�����ģ��
	write(12,*)  br%e2,					tab,'br%e2'
	write(12,*)  br%g2	,				tab,'br%g2'	!��Ȧ���ԡ�����ģ��
	write(12,*)  br%mu	,				tab,'br%mu'			!Ħ��ϵ��
	write(12,*)  br%nu1	,				tab,'br%nu1'		!�����嵯��ģ��
	write(12,*)  br%nu2	,				tab,'br%nu2'		!�����嵯��ģ��
	write(12,*)  br%rho1,				tab,'br%rho1'			!�ܶ�
	write(12,*)  br%rho2	,			tab,'br%rho2'		!�ܶ�
	write(12,*)  br%aload	,			tab,'br%aload'	!��������غ�
	write(12,*)  br%rload	,			tab,'br%rload'	!��о����غ�
	write(12,*)  br%mload	,			tab,'br%mload'	!��о����غ�
	write(12,*)  br%omgi,				tab,'br%omgi'		!��Ȧ���ٶ�
	write(12,*)  br%omgo,				tab,'br%omgo'		!��Ȧ���ٶ�
	write(12,*)  br%startang,		tab,'br%startang' !��ʼ������ķ�λ��
	write(12,*)  br%m,					tab,'br%m'
	write(12,*)  br%n,					tab,'br%n'
	write(12,*)  br%nsec,				tab,'br%nsec'
	write(12,*)  br%intmthd	,		tab,'br%intmthd'!1:�Ӵ���->������;2:������->�Ӵ��� ����Ħ����
	write(12,*)  br%crpmthd	,		tab,'br%crpmthd'!1:�Ӵ���->������;2:������->�Ӵ��ǡ�����Ħ����
	write(12,*)  br%init,				tab,'br%init'
	write(12,*) '***********************'
ENDSUBROUTINE OPTLOG