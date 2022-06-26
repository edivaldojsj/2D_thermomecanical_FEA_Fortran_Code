!**************************************************************************************************
!ANÁLISE TÉRMICA-ESTRUTURAL BIDIMENSIONAL DA BARRAGEM DE CONTRAFORTE DE ITAIPU UTILISANDO MEF

!ALGORITMO DESENVOLVIDO DE ACORDO COM O PLANO DE TRABAHO APRESENTADO À BANCA
!EXAMINADORA DO CURSO DE ENGENHARIA CIVIL DE INFRAESTRUTURA DA UNILA, COMO
!PARTE DOS REQUISITOS PARA OBTENÇÃO DO GRAU DE BACHAREL EM ENGENHARIA CIVIL
!
!DISCENTE:      EDIVALDO JOSÉ DA SILVA JUNIOR
!ORIENTADOR:    AREF KALILO LIMA KZAM
!
!PROGRAMA - CÁLCULO CAMPO DE TEMPERATURA, TENSÕES E DESLOCAMENTOS TÉRMICOS 2D



!******************************** PRÉ-PROCESSAMENTO *************************************

!****************************************************************************************************************************
! Variáveis para alocar valores do TXT  ****
!****************************************************************************************************************************
        program                 THERMAL_STRESS_2D
        implicit none

        INTEGER(4)                  i, j, k, m, n !CONTADORES
        INTEGER                  NODE_NUMBER,ELEM_NUMBER,MATERIAL_NUMBER,ELEM_NBC_NUMBER
        INTEGER                  NODE_RESTR_NUMBER, NODE_NBC_NUMBER, STANDARD_TEMP
        REAL,    ALLOCATABLE::   NODE_POSITION(:,:), MATERIAL_PROP(:,:)
        INTEGER, ALLOCATABLE::   ELEM_MATERIAL(:), ID_STRESS_STRAIN(:), ELEM_NODE(:,:), ID_ELEM_NBC(:,:)
        INTEGER, ALLOCATABLE::   ID_NODE_NBC(:,:)
        REAL,    ALLOCATABLE::   NBC_FORCA_X(:), NBC_FORCA_Y(:), NBC_DESLOC_X(:), NBC_DESLOC_Y(:), NBC_TEMP(:)
        REAL,    ALLOCATABLE::   NBC_FONTE_CALOR(:), NBC_FLUXO_X(:), NBC_FLUXO_Y(:), NBC_CONVEC_X(:), NBC_CONVEC_Y(:)
        INTEGER, ALLOCATABLE::   NODE_RESTR(:,:)


!****************************************************************************************************************************
!*** PROCESSAMENTO - VARIÁVEIS PARA O CÁLCULO DO CAMPO DE TEMPERATURA 2D  ****
!****************************************************************************************************************************
        REAL                     a, b, ky, kx, k_coef
        REAL,    ALLOCATABLE::   L_CONDUCTION_MATRIX(:,:,:),  L_CONVECTION_MATRIX(:,:,:), r_beta_t(:,:), RQ_VETOR(:,:)
        REAL,    ALLOCATABLE::   L_THERMAL_MATRIX(:,:,:), r_ELEM_VECTOR(:,:), G_THERMAL_MATRIX(:,:)
        REAL,    ALLOCATABLE::   G_THERMAL_MATRIX_BACKUP(:,:), V_GLOBAL(:), VETOR_TEMP_NODE(:), VETOR_TEMP_NODE_mult(:)
        REAL                     ALPHA, BETA, beta_fluxo, beta_convec, convec_media, temp_media
        REAL,    ALLOCATABLE::   s(:), t(:), N_T(:,:), VETOR_ELEM_TEMP(:,:), CENTER_ELEM_TEMP(:), VETOR_PIVOT(:)
        REAL                     TIE1, TIE2, TIE3

!****************************************************************************************************************************
!*** PROCESSAMENTO -  VARIÁVEIS PARA O CÁLCULO DOS DESLOCAMENTOS E TENSÕES TÉRMICAS 2D  ****
!****************************************************************************************************************************
        REAL,    ALLOCATABLE::   F1(:),F2(:),F3(:),B1(:),B2(:),B3(:)
        REAL,    ALLOCATABLE::   C1(:),C2(:),C3(:),ELEMENT_AREA(:)
        REAL                     X1, X2, X3, Y1, Y2, Y3
        REAL                     MOD_ELAST, POISSON,AREA, hA, CONST1, CONST2, CONST3, const4
        REAL                     DIV_CONST1, aldel_t, CONST5, CONST6, CONST7, CONST8, CONST9
        REAL,    ALLOCATABLE::   B_MATRIX_TRANSP(:,:,:),B_MATRIX(:,:,:), C_MATRIX(:,:,:), ELEM_THICKNESS(:)

        REAL,    ALLOCATABLE::   LOCAL_TS_MATRIX(:,:,:), G_TS_MATRIX(:,:), G_TS_MATRIX_BACKUP(:,:)
        REAL,    ALLOCATABLE::   G_FORCE_VECTOR(:), L_re_VETOR(:,:), G_DISP_VECTOR(:), ELEM_STRESS(:,:)
        INTEGER                  LIN_L, LIN_G, COL_L, COL_G


        REAL                     dx1, dy1, dx2, dy2, dx3, dy3
        REAL                     ELEM_AVE_TEMP, DELTA_TEMP
        REAL                     AVE_DISP_ELEM, AVE_THERMAL_ELEM

!****************************************************************************************************************************
!*** PROCESSAMENTO - VARIÁVEIS PARA INVETER MATRIZ LU=F  ****
!****************************************************************************************************************************
        REAL,    ALLOCATABLE::   THERMAL_INDX(:), THERMAL_STRESS_INDX(:) !Variável utilizada pela biblioteca de resolução de sistemas LU
        REAL                        ddd     !Variável utilizada pela biblioteca de resolução de sistemas LU

!****************************************************************************************************************************
!***  PÓS-PRECESSAMENTO - VARIÁVEIS PARA PLOTAR RESULTADOS ****
!****************************************************************************************************************************
        REAL                DISP_X_MAX, DISP_X_MIN, DISP_Y_MAX, DISP_Y_MIN

        REAL                X_MAX, Y_MAX
        REAL,                ALLOCATABLE::  MATRIX_VISUAL_2D(:,:), MATRIX_VISUAL_2D_2(:,:)
        INTEGER              X, Y
        REAL             STAND_X, STAND_Y

!***********************************************************************************************************************************
! PRÉ-PROCESSAMENTO:    LEITURA DE DADOS DO TXT
!***********************************************************************************************************************************
!open(unit=1,file="teste1_thermal_stress_Example 7_6_milimetros.txt")
!open(unit=1,file="ENTRADA1_EXEMPLO_VALIDACAO.txt")
!open(unit=1,file="ENTRADA_BLOCO_E6_THERMAL_STRESS_teste1.txt")
open(unit=1,file="ENTRADA_BLOCO_E6_THERMAL_STRESS_PLANE_STRAIN_WINTER.txt")
!open(unit=1,file="ENTRADA_BLOCO_E6_THERMAL_STRESS_PLANE_STRAIN_SUMMER_11_01_2010.txt")

    read(1,*)   !ESPAÇO RESERVADO PARA O CABEÇALHO DO PROBLEMA NO TXT
    read(1,*)
    read(1,*)   !   ''
    read(1,*)   !   ''
    read(1,*)   !   ''
    read(1,*)   !   ''
    read(1,*)   !   ''
    read(1,*)   !   ''
    read(1,*)   !   ''
    read(1,*)   !   ''
    read(1,*)   !   ''
    read(1,*)   !   ''
    read(1,*)   !   ''
    read(1,*)   !   ''
    read(1,*)   !   ''
    read(1,*)   !   ''
    read(1,*)   !   ''
    read(1,*)   !   ''
    read(1,*)   !   ''

    read(1,*)
!**************************************************************************************************
    write(*,*)"TOTAL NODAL, TOTAL ELEM, MATERIAL NUMBERS, ELEM.NBC NUMBER, NODE NBC NUMBER, NODE_RESTR_NUMBER, STANDARD TEMPETURE"
    read(1,*)  NODE_NUMBER,ELEM_NUMBER,MATERIAL_NUMBER,ELEM_NBC_NUMBER, NODE_NBC_NUMBER,NODE_RESTR_NUMBER, STANDARD_TEMP
    write(*,*) NODE_NUMBER,ELEM_NUMBER,MATERIAL_NUMBER,ELEM_NBC_NUMBER, NODE_NBC_NUMBER, NODE_RESTR_NUMBER, STANDARD_TEMP
    write(*,*)

    read(1,*)
    read(1,*)
!**************************************************************************************************
    STAND_X=0.D0; STAND_Y=0.D0

    !write(*,*)"ELEMENT STANDARD SIZE"
    read(1,*)STAND_X, STAND_Y
    !write(1,*)STAND_X, STAND_Y

    read(1,*)
    read(1,*)
!**************************************************************************************************
    allocate(NODE_POSITION(NODE_NUMBER ,2)); NODE_POSITION=0.D0
    !write(*,*)"NODE, X, Y (COORDENADAS)"

    do i=1, NODE_NUMBER
    read(1,*) n, (NODE_POSITION(n,j),j=1,2)
    !write(*,*) n, (NODE_POSITION(n,j),j=1,2)
    end do

    read(1,*)
    read(1,*)

    !NODE_POSITION=NODE_POSITION*1000! CONVERTENDO METROS EM MILÍMETROS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!**************************************************************************************************
allocate(MATERIAL_PROP(MATERIAL_NUMBER,4));MATERIAL_PROP=0.D0
    !write(*,*)
    !write(*,*)"MATERIAL, E(mod.Elast)[GPa], v(coef.Poisson), k(Cond.Termica), alpha(coef.dilat.Term.)[C]"

    do i=1, MATERIAL_NUMBER
    read(1,*) n,(MATERIAL_PROP(n,j),j=1,4)
    !write(*,*)n,(MATERIAL_PROP(n,j),j=1,4)
    end do

    read(1,*)
    read(1,*)
!**************************************************************************************************
allocate(ELEM_MATERIAL(ELEM_NUMBER),ELEM_THICKNESS(ELEM_NUMBER));ELEM_MATERIAL=0.D0; ELEM_THICKNESS=0.D0
allocate(ID_STRESS_STRAIN(ELEM_NUMBER), ELEM_NODE(ELEM_NUMBER,3));ID_STRESS_STRAIN=0.D0; ELEM_NODE=0.D0
    !write(*,*)
    !write(*,*)"ELEMENT, ID_STRESS_STRAIN, ELEM_THICKNESS´[m], PROP, Node 1, Node 2, Node 3"

    do i=1, ELEM_NUMBER
    read(1,*) n, ID_STRESS_STRAIN(n), ELEM_THICKNESS(n), ELEM_MATERIAL(n),(ELEM_NODE(n,j),j=1,3)
    !write(*,*)n, ID_STRESS_STRAIN(n), ELEM_THICKNESS(n), ELEM_MATERIAL(n),(ELEM_NODE(n,j),j=1,3)
    end do

    read(1,*)
    read(1,*)

!**************************************************************************************************
    allocate(ID_ELEM_NBC(ELEM_NUMBER,4)); ID_ELEM_NBC=0.D0

    !write(*,*)
    !write(*,*)"ELEMENT, CONDICAO NATURAL DE CONTORNO(NBC) NBC1, NBC2, NBC3 [0=NAO, 1=SIM]"

    do i=1, ELEM_NBC_NUMBER
    read(1,*) n, (ID_ELEM_NBC(n,j),j=1,3)
    !write(*,*) n,(ID_ELEM_NBC(n,j),j=1,3)
    end do

    read(1,*)
    read(1,*)
!**************************************************************************************************
    allocate(ID_NODE_NBC(NODE_NUMBER,6));ID_NODE_NBC=0.D0

    !write(*,*)
    !write(*,*)"NODE, FORCA, DESLOC, FONTE(Q), TEMPERATURA(T) FLUXO(q),CONVECÇÃO(h) [0=não, 1=sim]"

    do i=1, NODE_NBC_NUMBER
    read(1,*)n,(ID_NODE_NBC(n,j),j=1,6)
   ! write(*,*)n,(ID_NODE_NBC(n,j),j=1,6)
    end do

    read(1,*)
    read(1,*)
!**************************************************************************************************
    allocate(NBC_FORCA_X(NODE_NUMBER), NBC_FORCA_Y(NODE_NUMBER));NBC_FORCA_X=0.D0;NBC_FORCA_Y=0.D0
    allocate(NBC_DESLOC_X(NODE_NUMBER), NBC_DESLOC_Y(NODE_NUMBER));NBC_DESLOC_X=0.D0;NBC_DESLOC_Y=0.D0
    allocate(NBC_FONTE_CALOR(NODE_NUMBER),NBC_TEMP(NODE_NUMBER)); NBC_FONTE_CALOR=0.d0;  NBC_TEMP=0.D0
    allocate(NBC_FLUXO_X(NODE_NUMBER), NBC_FLUXO_Y(NODE_NUMBER))
    allocate(NBC_CONVEC_X(NODE_NUMBER), NBC_CONVEC_Y(NODE_NUMBER))
    NBC_FLUXO_X=0.D0; NBC_FLUXO_Y=0.D0; NBC_CONVEC_X=0.D0; NBC_CONVEC_Y=0.D0

    !write(*,*)
    !write(*,*)"Node, FORCA X, FORCA Y, DESLOC X, DESLOC Y"

    do j=1, NODE_NBC_NUMBER
    read(1,*) n, NBC_FORCA_X(n), NBC_FORCA_Y(n), NBC_DESLOC_X(n), NBC_DESLOC_Y(n)
    !write(*,*)n, NBC_FORCA_X(n), NBC_FORCA_Y(n), NBC_DESLOC_X(n), NBC_DESLOC_Y(n)
    end do

    read(1,*)
    read(1,*)
!**************************************************************************************************
    !write(*,*)
    !write(*,*)"Node, Q, T, qx, qy, hx, hy"

    do j=1, NODE_NBC_NUMBER
    read(1,*) n, NBC_FONTE_CALOR(n),NBC_TEMP(n), NBC_FLUXO_X(n), NBC_FLUXO_Y(n), NBC_CONVEC_X(n), NBC_CONVEC_Y(n)
    !write(*,*)n, NBC_FONTE_CALOR(n),NBC_TEMP(n), NBC_FLUXO_X(n), NBC_FLUXO_Y(n), NBC_CONVEC_X(n), NBC_CONVEC_Y(n)
    end do

    read(1,*)
    read(1,*)
!**************************************************************************************************
allocate(NODE_RESTR(NODE_NUMBER, 2)); NODE_RESTR=0.D0


   ! write(*,*)"Node,RESTR X, RESTR Y"

    do j=1, NODE_RESTR_NUMBER
    read(1,*) n, (NODE_RESTR(n,k),k=1,2)
    !write(*,*)n, (NODE_RESTR(n,k),k=1,2)
    end do


CLOSE(1)

WRITE(*,*)">> PRE-PROCESSAMENTO OK"
!***********************************************************************************************************************************



!******************************** PROCESSAMENTO *************************************

!**********************************************************************************
!**********************************************************************************
NODE_POSITION=NODE_POSITION*1000 ! convertendo as coordenadas nodais m --> mm

!**********************************************************************************


!***********************************************************************************************************************************
!PARTE 1 - THERMAL
!***********************************************************************************************************************************


!**********************************************************************************
! COMPONENTES DA MATRIZ DE RIGIDEZ DO ELEMENTO
!**********************************************************************************
allocate(F1(ELEM_NUMBER),F2(ELEM_NUMBER),F3(ELEM_NUMBER)); F1=0.D0;F2=0.D0;F3=0.D0
allocate(B1(ELEM_NUMBER),B2(ELEM_NUMBER),B3(ELEM_NUMBER)); B1=0.D0;B2=0.D0;B3=0.D0
allocate(C1(ELEM_NUMBER),C2(ELEM_NUMBER),C3(ELEM_NUMBER)); C1=0.D0;C2=0.D0;C3=0.D0
allocate(ELEMENT_AREA(ELEM_NUMBER)); ELEMENT_AREA=0.D0

allocate(L_CONDUCTION_MATRIX(ELEM_NUMBER,3,3));L_THERMAL_MATRIX=0.D0

do i=1, ELEM_NUMBER

    X1=NODE_POSITION(ELEM_NODE(i,1),1)
    X2=NODE_POSITION(ELEM_NODE(i,2),1)
    X3=NODE_POSITION(ELEM_NODE(i,3),1)

    Y1=NODE_POSITION(ELEM_NODE(i,1),2)
    Y2=NODE_POSITION(ELEM_NODE(i,2),2)
    Y3=NODE_POSITION(ELEM_NODE(i,3),2)

    F1(i)= X2*Y3-X3*Y2
    F2(i)= X3*Y1-X1*Y3
    F3(i)= X1*Y2-X2*Y1

    B1(i)= Y2-Y3
    B2(i)= Y3-Y1
    B3(i)= Y1-Y2

    C1(i)= X3-X2
    C2(i)= X1-X3
    C3(i)= X2-X1

    ELEMENT_AREA(i) = (F1(i)+F2(i)+F3(i))/2


!******************************************************************************************************************
! MATRIZ DE RIGIDEZ DO ELEMENTO
!******************************************************************************************************************
AREA=ELEMENT_AREA(i)
k_coef=MATERIAL_PROP(ELEM_MATERIAL(i),3)

L_CONDUCTION_MATRIX(i,1,1)= k_coef*(B1(i)**2+C1(i)**2)*(1/(4*AREA))
L_CONDUCTION_MATRIX(i,1,2)= k_coef*(B1(i)*B2(i)+C1(i)*C2(i))*(1/(4*AREA))
L_CONDUCTION_MATRIX(i,1,3)= k_coef*(B1(i)*B3(i)+C1(i)*C3(i))*(1/(4*AREA))

L_CONDUCTION_MATRIX(i,2,1)= k_coef*(B1(i)*B2(i)+C1(i)*C2(i))*(1/(4*AREA))
L_CONDUCTION_MATRIX(i,2,2)= k_coef*(B2(i)**2+C2(i)**2)*(1/(4*AREA))
L_CONDUCTION_MATRIX(i,2,3)= k_coef*(B2(i)*B3(i)+C2(i)*C3(i))*(1/(4*AREA))

L_CONDUCTION_MATRIX(i,3,1)= k_coef*(B1(i)*B3(i)+C1(i)*C3(i))*(1/(4*AREA))
L_CONDUCTION_MATRIX(i,3,2)= k_coef*(B2(i)*B3(i)+C2(i)*C3(i))*(1/(4*AREA))
L_CONDUCTION_MATRIX(i,3,3)= k_coef*(B3(i)**2+C3(i)**2)*(1/(4*AREA))

end do


!******************************************************************************************************************
!REALIZAR AS SOMAS (Kk+Kp+K_alpha)d=rq+r_beta (PROBLEMAS QUE ENVOLVEM CONVECÇÃO E FONTE DE CALOR)
!******************************************************************************************************************
allocate(L_THERMAL_MATRIX(ELEM_NUMBER,3,3), r_ELEM_VECTOR(ELEM_NUMBER,3)); L_THERMAL_MATRIX=0.D0; r_ELEM_VECTOR=0.D0

do i=1, ELEM_NUMBER
do j=1,3
 !   r_ELEM_VECTORi,j) = RQ_VETOR(i,j)+r_beta_t(i,j) !FONTE DE CALOR

do k=1,3

    L_THERMAL_MATRIX(i,j,k) = L_THERMAL_MATRIX(i,j,k) + L_CONDUCTION_MATRIX(i,j,k) ! + K_alpha(i,j,k)

end do
end do
end do



!******************************************************************************************************************
!***********MATRIZ GLOBAL DE RIGIDEZ ***********
!******************************************************************************************************************
allocate(G_THERMAL_MATRIX(NODE_NUMBER,NODE_NUMBER), V_GLOBAL(NODE_NUMBER)); G_THERMAL_MATRIX=0.D0; V_GLOBAL=0.D0
allocate(G_THERMAL_MATRIX_BACKUP(NODE_NUMBER,NODE_NUMBER)); G_THERMAL_MATRIX_BACKUP=0.D0


DO I=1, ELEM_NUMBER
do j=1,3

   V_GLOBAL(ELEM_NODE(i,j))= V_GLOBAL(ELEM_NODE(i,j)) !+ r_ELEM_VECTOR(i,j) !FONTE DE CALOR

do k=1, 3

    G_THERMAL_MATRIX(ELEM_NODE(i,j),ELEM_NODE(i,k))=G_THERMAL_MATRIX(ELEM_NODE(i,j),ELEM_NODE(i,k))+L_THERMAL_MATRIX(i,j,k)
    G_THERMAL_MATRIX_BACKUP(ELEM_NODE(i,j),ELEM_NODE(i,k))=G_THERMAL_MATRIX(ELEM_NODE(i,j),ELEM_NODE(i,k))

end do
end do
END DO

!******************************************************************************************************************
!***********APLICAR TEMPERATURAS DE CONTORNO NO VETOR GLOBAL DE TEMPERATURAS
!******************************************************************************************************************
ALLOCATE(VETOR_TEMP_NODE(NODE_NUMBER));VETOR_TEMP_NODE=0.D0

do i=1, NODE_NUMBER
    if(ID_NODE_NBC(i,4)==1)then
    VETOR_TEMP_NODE(i)=NBC_TEMP(i)
    endif
end do

!******************************************************************************************************************
!VISUALIZAR SE AS CONDIÇÕES FORAM ALOCADAS CORRETAMENTE

!    WRITE(*,*)"***************************************"
!    WRITE(*,*)"MATRIZ GLOBAL DE RIGIDEZ"
!    CALL MATRIXPRINT (G_THERMAL_MATRIX, NODE_NUMBER, NODE_NUMBER)

!    WRITE(*,*)"***************************************"
!    WRITE(*,*)"VETOR DE FLUXO"
!    do i=1,NODE_NUMBER
!    WRITE(*,*) i,  V_GLOBAL(i)
!    end do

!    WRITE(*,*)"***************************************"
!    WRITE(*,*)"CONDICAO DE CONTORNO: TEMPERATURA"
!    do i=1,NODE_NUMBER
!    WRITE(*,*) i,  VETOR_TEMP_NODE(i)
!    end do


!******************************************************************************************************************
! DELETAR TODAS AS LINHAS QUE REPRESENTAM OS NÓS QUE POSSUEM NBC DE TEMPERATURA
!******************************************************************************************************************
do i=1, NODE_NUMBER
    IF(ID_NODE_NBC(i,4)==1)THEN
       G_THERMAL_MATRIX(i,:)=0.D0
       G_THERMAL_MATRIX(i,i)=1.D0
    END IF
END DO

!******************************************************************************************************************
! MULTIPLICAR AS COLUNAS DOS NÓS COM NBC DE TEMPERATURA PELO VETOR NBC_TEMP
!******************************************************************************************************************!
allocate(VETOR_TEMP_NODE_mult(NODE_NUMBER));VETOR_TEMP_NODE_mult=1.D0
allocate(VETOR_PIVOT(NODE_NUMBER));VETOR_PIVOT=0.D0
    do i=1, NODE_NUMBER
        if(ID_NODE_NBC(i,4)==1)then
           VETOR_TEMP_NODE_mult(i)=VETOR_TEMP_NODE(i)
        end if
    end do


!******************************************************************************************************************
! SOMAR O (SOMATÓRIO DAS COLULAS * NBC TEMP) DE CADA NÓ COM O VETOR DE FLUXO. INVERTER SINAL!!!
!******************************************************************************************************************
do i=1, NODE_NUMBER
do j=1, NODE_NUMBER
    IF(ID_NODE_NBC(i,4)==1)THEN
       VETOR_PIVOT(j)= VETOR_PIVOT(j)+G_THERMAL_MATRIX(j,i)*VETOR_TEMP_NODE_mult(i)
    END IF
END DO
END DO


V_GLOBAL=V_GLOBAL+(-vetor_PIVOT) !AQUI INVERTE SINAL

do i=1, NODE_NUMBER
    if(ID_NODE_NBC(i,4)==1)then
    V_GLOBAL(i)=0.D0
    endif
end do



!********************************************************************
!VISUALIZAR SE AS CONDIÇÕES FORAM ALOCADAS CORRETAMENTE


!WRITE(*,*)"***************************************"
!    WRITE(*,*)"MATRIZ GLOBAL DE RIGIDEZ APOS CONDICOES DE CONTORNO"
!    CALL MATRIXPRINT (G_THERMAL_MATRIX, NODE_NUMBER, NODE_NUMBER)
!WRITE(*,*)"***************************************"
!WRITE(*,*)"vetor de temperatura de multiplicacao"
! do i=1,NODE_NUMBER
!    WRITE(*,*) i,  VETOR_TEMP_NODE_mult(i)
!    end do

!WRITE(*,*)"***************************************"
!WRITE(*,*)"VETOR DE SOMA DAS COLUNAS EXTRAIDAS"
! do i=1,NODE_NUMBER
!    WRITE(*,*) i,  VETOR_PIVOT(i)
!    end do

!WRITE(*,*)"***************************************"
!WRITE(*,*)"VETOR Global de fluxo"
! do i=1,NODE_NUMBER
!   WRITE(*,*) i,  V_GLOBAL(i)
!    end do

!**********************************************************************


!******************************************************************************************************************
        !SUBROTINAS QUE RESOLVEM O SISTEMA LU
!******************************************************************************************************************
allocate (THERMAL_INDX(NODE_NUMBER))

call ludcmp(G_THERMAL_MATRIX,NODE_NUMBER,NODE_NUMBER,THERMAL_INDX,ddd)
call lubksb(G_THERMAL_MATRIX,NODE_NUMBER,NODE_NUMBER,THERMAL_INDX,V_GLOBAL)
!******************************************************************************************************************




!******************************************************************************************************************
!   RESULTADO DAS TEMPERATURAS NODAIS
!******************************************************************************************************************
VETOR_TEMP_NODE=VETOR_TEMP_NODE+V_GLOBAL

!WRITE(*,*)"***************************************"
!    WRITE(*,*)"RESULTADO:"
!    WRITE(*,*)"NODE, TEMPETURE"
!    do i=1,NODE_NUMBER
!    WRITE(*,*) i,  VETOR_TEMP_NODE(i)
!    end do



!******************************************************************************************************************
!   COMPARAR AS TEMEPRATURAS DOS NÓS COM OS TERMÔMETROS INTERNOS DO BLOCO E6
!******************************************************************************************************************
!WRITE(*,*)"***************************************"
!WRITE(*,*)"TEMPERATURA DOS NOS DE VALIDACAO" COMPARAR COM OS TERMÔMETROS INTERNOS DA BARRAGEM

!    TIE1=(VETOR_TEMP_NODE(ELEM_NODE(370,1))+VETOR_TEMP_NODE(ELEM_NODE(370,2))+VETOR_TEMP_NODE(ELEM_NODE(370,3)))/3
!    TIE2=(VETOR_TEMP_NODE(ELEM_NODE(380,1))+VETOR_TEMP_NODE(ELEM_NODE(380,2))+VETOR_TEMP_NODE(ELEM_NODE(380,3)))/3
!    TIE3=(VETOR_TEMP_NODE(ELEM_NODE(436,1))+VETOR_TEMP_NODE(ELEM_NODE(436,2))+VETOR_TEMP_NODE(ELEM_NODE(436,3)))/3

!    WRITE(*,*) 'COMPARAR COM TI-E-1 -->', TIE1
!    WRITE(*,*) 'COMPARAR COM TI-E-2 -->', TIE2
!    WRITE(*,*) 'COMPARAR COM TI-E-2 -->', TIE3


!******************************************************************************************************************
WRITE(*,*)">> PROCESSAMENTO PARTE 1 - THERMAL OK"


!PROCESSAMENETO
!******************************************************************************************************************
!PART 2 - THERMAL STRESS
!******************************************************************************************************************


!**********************************************************************
! COMPONENTES DA MATRIZ DE RIGIDEZ DO ELEMENTO
!**********************************************************************
allocate(LOCAL_TS_MATRIX(ELEM_NUMBER,6,6));LOCAL_TS_MATRIX=0.D0
allocate(L_re_VETOR(ELEM_NUMBER,6));L_re_VETOR=0.D0

do i=1, ELEM_NUMBER

MOD_ELAST=MATERIAL_PROP(ELEM_MATERIAL(i),1)
POISSON=MATERIAL_PROP(ELEM_MATERIAL(i),2)

!***********************************************************************
! MATRIZ DE RIGIDEZ DO ELEMENTO
!***********************************************************************
AREA=ELEMENT_AREA(i)


    if(ID_STRESS_STRAIN(i)==1)then
!**********************************************
! ESTADO PLANO DE TENSÃO !!!!!!!
!**********************************************
    CONST1=(MOD_ELAST/(1-POISSON**2))
    CONST2=(MOD_ELAST*POISSON/(1-POISSON**2))
    CONST3= (MOD_ELAST*(1-POISSON)/(2*(1-POISSON**2)))
    DIV_CONST1=(4*AREA**2)
    hA=ELEM_THICKNESS(i)*ELEMENT_AREA(i)


    CONST4 = MATERIAL_PROP(ELEM_MATERIAL(i),4)

    else

!***************************************************
! ESTADO PLANO DE DEFORMAÇÃO !!!!!!!
!***************************************************
    CONST1=(MOD_ELAST*(1-POISSON)/((1+POISSON)*(1-2*POISSON)))
    CONST2=(MOD_ELAST*POISSON/((1+POISSON)*(1-2*POISSON)))
    CONST3= (MOD_ELAST*(1-2*POISSON)/(2*(1+POISSON)*(1-2*POISSON)))
    DIV_CONST1=(4*AREA**2)
    hA=ELEM_THICKNESS(i)*ELEMENT_AREA(i)

    CONST4 = (1+POISSON)*MATERIAL_PROP(ELEM_MATERIAL(i),4)

    end if


LOCAL_TS_MATRIX(i,1,1)= (B1(i)**2*CONST1/DIV_CONST1 + C1(i)**2*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,1,2)= ((C1(i)*B1(i))/DIV_CONST1*(CONST2 + CONST3))*hA
LOCAL_TS_MATRIX(i,1,3)= (B2(i)*B1(i)*CONST1/DIV_CONST1 + C2(i)*C1(i)*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,1,4)= (B1(i)*C2(i)*CONST2/DIV_CONST1 + B2(i)*C1(i)*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,1,5)= (B3(i)*B1(i)*CONST1/DIV_CONST1 + C3(i)*C1(i)*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,1,6)= (C3(i)*B1(i)*CONST2/DIV_CONST1 + B3(i)*C1(i)*CONST3/DIV_CONST1)*hA

LOCAL_TS_MATRIX(i,2,1)= ((C1(i)*B1(i))/DIV_CONST1*(CONST2 + CONST3))*hA
LOCAL_TS_MATRIX(i,2,2)= (C1(i)**2*CONST1/DIV_CONST1 + B1(i)**2*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,2,3)= (B2(i)*C1(i)*CONST2/DIV_CONST1 + C2(i)*B1(i)*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,2,4)= (C1(i)*C2(i)*CONST1/DIV_CONST1 + B1(i)*B2(i)*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,2,5)= (C1(i)*B3(i)*CONST2/DIV_CONST1 + B1(i)*C3(i)*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,2,6)= (C1(i)*C3(i)*CONST1/DIV_CONST1 + B1(i)*B3(i)*CONST3/DIV_CONST1)*hA

LOCAL_TS_MATRIX(i,3,1)= (B1(i)*B2(i)*CONST1/DIV_CONST1 + C2(i)*C1(i)*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,3,2)= (B2(i)*C1(i)*CONST2/DIV_CONST1 + B1(i)*C2(i)*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,3,3)= (B2(i)**2*CONST1/DIV_CONST1 + C2(i)**2*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,3,4)= (B2(i)*C2(i)*CONST2/DIV_CONST1 + B2(i)*C2(i)*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,3,5)= (B2(i)*B3(i)*CONST1/DIV_CONST1 + C3(i)*C2(i)*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,3,6)= (B2(i)*C3(i)*CONST2/DIV_CONST1 + B3(i)*C2(i)*CONST3/DIV_CONST1)*hA

LOCAL_TS_MATRIX(i,4,1)= (B1(i)*C2(i)*CONST2/DIV_CONST1 + B2(i)*C1(i)*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,4,2)= (C1(i)*C2(i)*CONST1/DIV_CONST1 + B2(i)*B1(i)*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,4,3)= (C2(i)*B2(i)/DIV_CONST1*(CONST2 + CONST3))*hA
LOCAL_TS_MATRIX(i,4,4)= (C2(i)**2*CONST1/DIV_CONST1 + B2(i)**2*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,4,5)= (B3(i)*C2(i)*CONST2/DIV_CONST1 + B2(i)*C3(i)*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,4,6)= (C3(i)*C2(i)*CONST1/DIV_CONST1 + B2(i)*B3(i)*CONST3/DIV_CONST1)*hA

LOCAL_TS_MATRIX(i,5,1)= (B1(i)*B3(i)*CONST1/DIV_CONST1 + C1(i)*C3(i)*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,5,2)= (C1(i)*B3(i)*CONST2/DIV_CONST1 + C3(i)*B1(i)*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,5,3)= (B2(i)*B3(i)*CONST1/DIV_CONST1 + C3(i)*C2(i)*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,5,4)= (C2(i)*B3(i)*CONST2/DIV_CONST1 + B2(i)*C3(i)*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,5,5)= (B3(i)**2*CONST1/DIV_CONST1 + C3(i)**2*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,5,6)= (B3(i)*C3(i)/DIV_CONST1*(CONST2 + CONST3))*hA

LOCAL_TS_MATRIX(i,6,1)= ((B1(i)*C3(i)*CONST2)/DIV_CONST1 + (B3(i)*C1(i)*CONST3)/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,6,2)= ((C1(i)*C3(i)*CONST1)/DIV_CONST1 + (B1(i)*B3(i)*CONST3)/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,6,3)= ((B2(i)*C3(i)*CONST2)/DIV_CONST1 + (C2(i)*B3(i)*CONST3)/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,6,4)= ((C2(i)*C3(i)*CONST1)/DIV_CONST1 + (B2(i)*B3(i))*CONST3/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,6,5)= (((B3(i)*C3(i))*(CONST2 + CONST3))/DIV_CONST1)*hA
LOCAL_TS_MATRIX(i,6,6)= ((C3(i)**2*CONST1)/DIV_CONST1 + (B3(i)**2*CONST3)/DIV_CONST1)*hA


!********************************************************************************************
! vetor rb --> The equivalent load vector due to body forces is (AINDA NÃO HABILITADO)
!********************************************************************************************

!do i=1, ELEM_NUMBER

!    rb(i,1)= (bx(i))*h(i)*ELEMENT_AREA(i)/3
!    rb(i,2)= (by(i))*h(i)*ELEMENT_AREA(i)/3
!    rb(i,3)= (bx(i))*h(i)*ELEMENT_AREA(i)/3
!    rb(i,4)= (by(i))*h(i)*ELEMENT_AREA(i)/3
!    rb(i,5)= (bx(i))*h(i)*ELEMENT_AREA(i)/3
!    rb(i,6)= (by(i))*h(i)*ELEMENT_AREA(i)/3
!end do

!******************************************************************************************************************
! vetor rq --> The equivalent load vector due to body forces is (AINDA NÃO HABILITADO)
!******************************************************************************************************************

!contribuição de cargas distribuídas aplicadas nos elementos.



!******************************************************************************************************************
! CONTRIBUIÇÃO DA TEMPERATURA --> The equivalent load vector from initial strains due to temperature change is
!******************************************************************************************************************
! AS TEMPERATURAS NODAIS SÃO UTILIZADAS COMO CONDIÇÃO DE FORÇA INICIAL DEVIDO A DILATAÇÃO TÉRMICA.
! VARIÁVEL DA TEMPERATURA NODAL VETOR_TEMP_NODE(ELEM_NODE(i,J)

L_re_VETOR(i,1)=const4*(VETOR_TEMP_NODE(ELEM_NODE(i,1))-STANDARD_TEMP)*B1(i)*(CONST1+CONST2)*hA/(2*ELEMENT_AREA(i))
L_re_VETOR(i,2)=const4*(VETOR_TEMP_NODE(ELEM_NODE(i,1))-STANDARD_TEMP)*C1(i)*(CONST1+CONST2)*hA/(2*ELEMENT_AREA(i))
L_re_VETOR(i,3)=const4*(VETOR_TEMP_NODE(ELEM_NODE(i,2))-STANDARD_TEMP)*B2(i)*(CONST1+CONST2)*hA/(2*ELEMENT_AREA(i))
L_re_VETOR(i,4)=const4*(VETOR_TEMP_NODE(ELEM_NODE(i,2))-STANDARD_TEMP)*C2(i)*(CONST1+CONST2)*hA/(2*ELEMENT_AREA(i))
L_re_VETOR(i,5)=const4*(VETOR_TEMP_NODE(ELEM_NODE(i,3))-STANDARD_TEMP)*B3(i)*(CONST1+CONST2)*hA/(2*ELEMENT_AREA(i))
L_re_VETOR(i,6)=const4*(VETOR_TEMP_NODE(ELEM_NODE(i,3))-STANDARD_TEMP)*C3(i)*(CONST1+CONST2)*hA/(2*ELEMENT_AREA(i))


!******************************************************************************************************************
end do




!********************************************************************
!VISUALIZAR SE AS CONDIÇÕES FORAM ALOCADAS CORRETAMENTE


!******************************************
! PRINT MATRIZ DE RIGIDEZ DO ELEMENTO
!******************************************
!do  i=1,ELEM_NUMBER
!    write(*,*)'***************************************************************************************************************'
!    write(*,*)'ELEMENT:', i
!    write(*,*)
!    write(*,*)'                  1                  2              3               4                 5                 6'
!do  j=1,6
!   write(*,*)j, (LOCAL_TS_MATRIX(i,j,k),k=1,6)
!end do
!end do
!******************************************


!write(*,*)'***************************************************************************************************************'

!do i=1, ELEM_NUMBER
!    write(*,*)'***************************************************************************************************************'
!    write(*,*)'ELEMENT:', i
!    WRITE(*,*)'ALPHA:', MATERIAL_PROP(ELEM_MATERIAL(i),4)
!    WRITE(*,*)'DELTA T NODE 1:', VETOR_TEMP_NODE(ELEM_NODE(i,1))-STANDARD_TEMP
!    WRITE(*,*)'DELTA T NODE 2:', VETOR_TEMP_NODE(ELEM_NODE(i,2))-STANDARD_TEMP
!    WRITE(*,*)'DELTA T NODE 3:', VETOR_TEMP_NODE(ELEM_NODE(i,3))-STANDARD_TEMP
 ! DO J=1, 6

!    WRITE(*,*)ELEM_NODE(i,1), L_re_VETOR(i,J)
!    end do
!end do


!****************************************************************************
        !CALCULA MATRIZ GLOBAL
!****************************************************************************
ALLOCATE (G_TS_MATRIX(2*NODE_NUMBER, 2*NODE_NUMBER)); G_TS_MATRIX=0.D0
ALLOCATE (G_TS_MATRIX_BACKUP(2*NODE_NUMBER, 2*NODE_NUMBER)); G_TS_MATRIX=0.D0

DO I=1, ELEM_NUMBER
DO J=1, 3 !NUMERO DE NOS DO ELEMENTO
DO K=1, 2 ! NUMERO DE GRAUS DE LIBERDADE DO NO

        DO M=1,3 !NUMERO DE NOS DO ELEMENTO
        DO N= 1,2 !NUMERO DE GRAUS DE LIBERDADE

        LIN_G=2*(ELEM_NODE(i,j)-1)+K
        COL_G=2*(ELEM_NODE(i,M)-1)+N

        LIN_L=2*(j-1)+k
        COL_L= 2*(M-1)+N

        G_TS_MATRIX(LIN_G,COL_G)=G_TS_MATRIX(LIN_G,COL_G)+LOCAL_TS_MATRIX(i,LIN_L,COL_L )
        G_TS_MATRIX_BACKUP(LIN_G,COL_G)=G_TS_MATRIX(LIN_G,COL_G)

        END DO
        END DO
END DO
END DO
END DO


!****************************************************************************
        !aplicacao das restricoes na matriz global
!****************************************************************************

    DO i=1, NODE_NUMBER
    do j=1, 2
        IF(NODE_RESTR(i,j)==1) THEN
        if(j==1)G_TS_MATRIX(2*(i-1)+1,:)=0.D0
        if(j==1)G_TS_MATRIX(:,2*(i-1)+1)=0.D0
        if(j==1)G_TS_MATRIX(2*(i-1)+1,2*(i-1)+1)=1.0D0

        if(j==2)G_TS_MATRIX(2*(i-1)+2,:)=0.D0
        if(j==2)G_TS_MATRIX(:,2*(i-1)+2)=0.D0
        if(j==2)G_TS_MATRIX(2*(i-1)+2,2*(i-1)+2)=1.0D0
        END IF
    END DO
    end do



!********************************************************************
!VISUALIZAR SE AS CONDIÇÕES FORAM ALOCADAS CORRETAMENTE

!WRITE(*,*)
!WRITE(*,*)'************************************************************************'
!write(*,*)'MATRI DE RIGIDEZ GLOBAL'
!CALL MATRIXPRINT  (G_TS_MATRIX_BACKUP, 2*NODE_NUMBER, 2*NODE_NUMBER)

!WRITE(*,*)
!WRITE(*,*)'************************************************************************'
!WRITE(*,*)'APOS RESTRICOES'

!CALL MATRIXPRINT  (G_TS_MATRIX, 2*NODE_NUMBER, 2*NODE_NUMBER)
!********************************************************************


!****************************************************************************
        !VETOR GLOBAL DE FORÇAS
!****************************************************************************
!SOMATÓRIO DE FORÇAS NODAIS, PRESSÃO NO ELEMENETO E FORÇA DEVIDO DILATAÇÃO TÉRMICA

ALLOCATE(G_FORCE_VECTOR(2*NODE_NUMBER));G_FORCE_VECTOR=0.D0
ALLOCATE(G_DISP_VECTOR(2*NODE_NUMBER));G_DISP_VECTOR=0.D0

    DO I=1,ELEM_NUMBER
    n=1

    do j=1,3
      do k=1,2

G_FORCE_VECTOR(2*(ELEM_NODE(i,j)-1)+k)=  L_re_VETOR(i,n) + G_FORCE_VECTOR(2*(ELEM_NODE(i,j)-1)+k)

        n=n+1

	END DO
	end do
    end do


!WRITE(*,*)'VETOR GLOBAL DE FORÇAS TERMICAS'
!do i=1, 2*node_number
!    write(*,*)i, G_FORCE_VECTOR(i)
!end do



!PREPARANDO O VETOR DE DESLOCAMENTO PARA A SUBROTINA QUE INVERTE MATRIZ

G_DISP_VECTOR=G_FORCE_VECTOR


 DO i=1, NODE_NUMBER
    do j=1, 2
        IF(NODE_RESTR(i,j)==1) THEN

        if(j==1)G_DISP_VECTOR(2*(i-1)+1)=0.D0 !APLICANDO RESTRIÇÃO
        if(j==2)G_DISP_VECTOR(2*(i-1)+2)=0.D0 !APLICANDO RESTRIÇÃO
        END IF
        end do
        end do



!****************************************************************************
WRITE(*,*)
WRITE(*,*)'VETOR GLOBAL DE FORCAS TERMICAS APOS RESTRICAO'

do i=1, 2*node_number
    write(*,*)i, G_DISP_VECTOR(i)
end do
!****************************************************************************




!****************************************************************************
        !SUBROTINAS QUE RESOLVEM O SISTEMA LU
!*****************************************************************************
allocate (THERMAL_STRESS_INDX(2*NODE_NUMBER))

call ludcmp(G_TS_MATRIX,2*NODE_NUMBER,2*NODE_NUMBER,THERMAL_STRESS_INDX,ddd)
call lubksb(G_TS_MATRIX,2*NODE_NUMBER,2*NODE_NUMBER,THERMAL_STRESS_INDX,G_DISP_VECTOR)



!****************************************************************************
! Cálculo das tensões normais em x, y e z
!*****************************************************************************
allocate(ELEM_STRESS(ELEM_NUMBER,3)); ELEM_STRESS=0.D0

do i=1, ELEM_NUMBER

       if(ID_STRESS_STRAIN(i)==1)then
!**********************************************
! ESTADO PLANO DE TENSÃO !!!!!!!
!**********************************************
    CONST4 = MATERIAL_PROP(ELEM_MATERIAL(i),4)
    CONST5=((MOD_ELAST)/(1-POISSON**2))
    CONST6=1.D0
    CONST7=(1-POISSON)/2
    else
!***************************************************
! ESTADO PLANO DE DEFORMAÇÃO !!!!!!!
!***************************************************
    CONST4 = (1+POISSON)*MATERIAL_PROP(ELEM_MATERIAL(i),4)
    CONST5=(MOD_ELAST/((1+POISSON)*(1-2*POISSON)))
    CONST6=(1-POISSON)
    CONST7=(1-2*POISSON)/2
    end if

dx1=G_DISP_VECTOR(2*(ELEM_NODE(i,1)-1)+1)
dy1=G_DISP_VECTOR(2*(ELEM_NODE(i,1)-1)+2)

dx2=G_DISP_VECTOR(2*(ELEM_NODE(i,2)-1)+1)
dy2=G_DISP_VECTOR(2*(ELEM_NODE(i,2)-1)+2)

dx3=G_DISP_VECTOR(2*(ELEM_NODE(i,3)-1)+1)
dy3=G_DISP_VECTOR(2*(ELEM_NODE(i,3)-1)+2)

ELEM_AVE_TEMP=(VETOR_TEMP_NODE(ELEM_NODE(i,1))+VETOR_TEMP_NODE(ELEM_NODE(i,2))+VETOR_TEMP_NODE(ELEM_NODE(i,3)))/3

DELTA_TEMP=ELEM_AVE_TEMP-STANDARD_TEMP
CONST8 = (B1(i)*dx1+B2(i)*dx2+B3(i)*dx3)/(2*ELEMENT_AREA(i))
CONST9 = (C1(i)*dy1+C2(i)*dy2+C3(i)*dy3)/(2*ELEMENT_AREA(i))

ELEM_STRESS(i,1)=CONST5*((CONST8-CONST4*DELTA_TEMP)*CONST6 + POISSON*CONST9-CONST4*DELTA_TEMP)        !TENSÃO EM X
ELEM_STRESS(i,2)=CONST5*((CONST8-CONST4*DELTA_TEMP)*POISSON + CONST6*CONST9-CONST4*DELTA_TEMP)        !TENSÃO EM Y
ELEM_STRESS(i,3)=CONST5*(CONST7*(C1(i)*dx1+B1(i)*dy1+C2(i)*dx2+B2(i)*dy2+C3(i)*dx3+B3(i)*dy3))  !TENSÃO EM Z

end do

write(*,*)
write(*,*)'constantes', const5
!********************************************************************
!VISUALIZAR RESULTADOS

!WRITE(*,*)
WRITE(*,*)"RESULTADO:"
WRITE(*,*)' NODE, NODE DISPLACEMENT:'

do i=1, 2*node_number
write(*,*)i, G_DISP_VECTOR(i)
end do


WRITE(*,*)
WRITE(*,*)"RESULTADO:"
WRITE(*,*)' ELEMENT, STRESS: X, Y, Z:'
do i=1, ELEM_NUMBER
WRITE(*,*)i, ELEM_STRESS(i,1), ELEM_STRESS(i,2), ELEM_STRESS(i,3)
end do

WRITE(*,*)">> PROCESSAMENTO PARTE 2 - THERMAL STRESS OK"
!*******************************************************************************************************



!******************************** POS-PROCESSAMENTO *************************************


if(STAND_X==0)then

    WRITE(*,*)">> POS-PROCESSAMENTO OK "

WRITE(*,*)
WRITE(*,*)"NAO HA MATRIZ 2D DE VISUALIZAÇÃO PORQUE OS ELEMENTOS NAO SAO UNIFORMES! =("

ELSE


!***************************************************************
!ROTINA PARA GERAR GRÁFICOS COLORIDOS NO GNOPLOT
!***************************************************************
X_MAX=0.D0;Y_MAX=0.D0

    do i=1, NODE_NUMBER
    if(X_MAX < NODE_POSITION(i,1))X_MAX=NODE_POSITION(i,1)
    if(Y_MAX < NODE_POSITION(i,2))Y_MAX=NODE_POSITION(i,2)
    end do

    X_MAX= X_MAX/STAND_X
    Y_MAX= Y_MAX/STAND_Y


!******************************************************************************************************************
!GERANDO GRÁFICOS - CAMPO DE TEMPERATURA
!******************************************************************************************************************
!allocate(MATRIX_VISUAL_2D(X_MAX,Y_MAX)); MATRIX_VISUAL_2D=100.0d0 !!!!!! !EXEMPLO DE VALIDAÇÃO
allocate(MATRIX_VISUAL_2D(X_MAX,Y_MAX)); MATRIX_VISUAL_2D=14.0d0 !!!!!! !BLOCO E6 INVERNO
!allocate(MATRIX_VISUAL_2D(X_MAX,Y_MAX)); MATRIX_VISUAL_2D=20.0d0 !!!!!! !BLOCO E6 VERÃO



do i=1, ELEM_NUMBER

    X= NODE_POSITION(ELEM_NODE(i,3),1)/STAND_X
    Y= NODE_POSITION(ELEM_NODE(i,3),2)/STAND_Y

AVE_THERMAL_ELEM=(VETOR_TEMP_NODE(ELEM_NODE(i,1))+VETOR_TEMP_NODE(ELEM_NODE(i,2))+VETOR_TEMP_NODE(ELEM_NODE(i,3)))/3

    MATRIX_VISUAL_2D(X,Y)=AVE_THERMAL_ELEM

end do

!**************************************************************************
!SAIDA_MATRIX_GRAPHIC.TXT - CADA POSIÇÃO MATRICIAL REPRESENTA UM ELEMENTO
!**************************************************************************
open(unit=2,file="SAIDA1_THERMAL_GRAPHIC_MATRIX.txt")
!PONTOS NODAIS
do i=1, Y_MAX
        write(2,*)(MATRIX_VISUAL_2D(j,i),j=1,X_MAX)
        !write(*,*)(MATRIX_VISUAL_2D(i,j),j=1,X_MAX)
end do

close(2)
!**************************************************************************



!******************************************************************************************************************
!GERANDO GRÁFICOS - CAMPO DE DESLOCAMENTOS EM X
!******************************************************************************************************************
MATRIX_VISUAL_2D=0.0005D0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!INVERNO
!MATRIX_VISUAL_2D=-0.0005D0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!VERÃO
!MATRIX_VISUAL_2D=0.d0       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!EXEMPLO DE VALIDAÇÃO

do i=1, ELEM_NUMBER

    X= NODE_POSITION(ELEM_NODE(i,3),1)/STAND_X
    Y= NODE_POSITION(ELEM_NODE(i,3),2)/STAND_Y

AVE_DISP_ELEM=(G_DISP_VECTOR(2*(ELEM_NODE(i,1)-1)+1)+G_DISP_VECTOR(2*(ELEM_NODE(i,2)-1)+1)+G_DISP_VECTOR(2*(ELEM_NODE(i,3)-1)+1))/3

    MATRIX_VISUAL_2D(X,Y)= AVE_DISP_ELEM

END DO
!**************************************************************************
!SAIDA_MATRIX_GRAPHIC.TXT - CADA POSIÇÃO MATRICIAL REPRESENTA UM ELEMENTO
!**************************************************************************
open(unit=3,file="SAIDA2_X_TS_DISP_GRAPHIC_MATRIX.txt")
!PONTOS NODAIS
do i=1, Y_MAX
        write(3,*)(MATRIX_VISUAL_2D(j,i),j=1,X_MAX)
end do

close(3)
!**************************************************************************



!******************************************************************************************************************
!GERANDO GRÁFICOS - CAMPO DE DESLOCAMENTOS EM Y
!******************************************************************************************************************
MATRIX_VISUAL_2D=2.0D0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!INVERNO
!MATRIX_VISUAL_2D=-0.01D0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!VERÃO
!MATRIX_VISUAL_2D=0.d0       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!EXEMPLO DE VALIDAÇÃO

do i=1, ELEM_NUMBER

        X= NODE_POSITION(ELEM_NODE(i,3),1)/STAND_X
        Y= NODE_POSITION(ELEM_NODE(i,3),2)/STAND_Y


AVE_DISP_ELEM=(G_DISP_VECTOR(2*(ELEM_NODE(i,1)-1)+2)+G_DISP_VECTOR(2*(ELEM_NODE(i,2)-1)+2)+G_DISP_VECTOR(2*(ELEM_NODE(i,3)-1)+2))/3

    MATRIX_VISUAL_2D(X,Y)= AVE_DISP_ELEM

END DO
!**************************************************************************
!SAIDA_MATRIX_GRAPHIC.TXT - CADA POSIÇÃO MATRICIAL REPRESENTA UM ELEMENTO
!**************************************************************************
open(unit=4,file="SAIDA3_Y_TS_DISP_GRAPHIC_MATRIX.txt")
!PONTOS NODAIS
do i=1, Y_MAX
        write(4,*)(MATRIX_VISUAL_2D(j,i),j=1,X_MAX)
end do

close(4)
!**************************************************************************





!******************************************************************************************************************
!GERANDO GRÁFICOS - CAMPO DE TENSÕES EM X
!******************************************************************************************************************
!MATRIX_VISUAL_2D=0.004D0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VALOR DEFINIDO DE ACORDO COM A AMPLITUDE DOS RESULTADOS
!do i=1, ELEM_NUMBER

 !    X= NODE_POSITION(ELEM_NODE(i,3),1)/STAND_X
 !    Y= NODE_POSITION(ELEM_NODE(i,3),2)/STAND_Y

!    MATRIX_VISUAL_2D(X,Y)=ELEM_STRESS(i,1)

!END DO
!**************************************************************************
!SAIDA_MATRIX_GRAPHIC.TXT - CADA POSIÇÃO MATRICIAL REPRESENTA UM ELEMENTO
!**************************************************************************
!open(unit=5,file="SAIDA4_X_TS__STRESS_GRAPHIC_MATRIX.txt")

!do i=1, Y_MAX
 !       write(5,*)(MATRIX_VISUAL_2D(j,i),j=1,X_MAX)
!end do

!close(5)
!**************************************************************************



!******************************************************************************************************************
!GERANDO GRÁFICOS - CAMPO DE TENSÕES EM Y
!******************************************************************************************************************
!MATRIX_VISUAL_2D=0.004D0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VALOR DEFINIDO DE ACORDO COM A AMPLITUDE DOS RESULTADOS
!do i=1, ELEM_NUMBER

 !    X= NODE_POSITION(ELEM_NODE(i,3),1)/STAND_X
 !   Y= NODE_POSITION(ELEM_NODE(i,3),2)/STAND_Y

  !  MATRIX_VISUAL_2D(X,Y)=ELEM_STRESS(i,2)

!END DO
!**************************************************************************
!SAIDA_MATRIX_GRAPHIC.TXT - CADA POSIÇÃO MATRICIAL REPRESENTA UM ELEMENTO
!**************************************************************************
!open(unit=6,file="SAIDA5_Y_TS__STRESS_GRAPHIC_MATRIX.txt")

!do i=1, Y_MAX
!write(6,*)(MATRIX_VISUAL_2D(j,i),j=1,X_MAX)
!end do

!close(6)
!**************************************************************************


!******************************************************************************************************************
!GERANDO GRÁFICOS - CAMPO DE TENSÕES EM Z
!******************************************************************************************************************
!MATRIX_VISUAL_2D=0.005D0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VALOR DEFINIDO DE ACORDO COM A AMPLITUDE DOS RESULTADOS
!do i=1, ELEM_NUMBER

 !X= NODE_POSITION(ELEM_NODE(i,3),1)/STAND_X
 !   Y= NODE_POSITION(ELEM_NODE(i,3),2)/STAND_Y

 !   MATRIX_VISUAL_2D(X,Y)=ELEM_STRESS(i,3)

!END DO
!**************************************************************************
!SAIDA_MATRIX_GRAPHIC.TXT - CADA POSIÇÃO MATRICIAL REPRESENTA UM ELEMENTO
!**************************************************************************
!open(unit=7,file="SAIDA6_Z_TS__STRESS_GRAPHIC_MATRIX.txt")

!do i=1, Y_MAX
!        write(7,*)(MATRIX_VISUAL_2D(j,i),j=1,X_MAX)
!end do

!close(7)
!**************************************************************************

WRITE(*,*)">> POS-PROCESSAMENTO OK"

WRITE(*,*)
WRITE(*,*)"UTILIZE O GNUPLOT PARA VISUALIZAR GRAFICAMENTE OS RESULTADOS =D"

end if





end program





! EXIBIR MATRIZES DE FORMA ORGANIZADA NO TERMINAL
subroutine   MATRIXPRINT  (a, m, n)
    integer           m, n, i, j, jref
      real a(m,n)
!
      do 2000  jref = 0,n-1,6
        print '(2X,6I12)',(j,j=jref+1,min(jref+6,n))
        do 1500  i = 1,m
          print '(I5,1P6E12.2)',i,(a(i,j),j=jref+1,min(jref+6,n))

 1500     continue
 2000   continue
!
      return
      end

SUBROUTINE ludcmp(a,n,np,indx,d)
        integer n,np,indx(n),nmax
        real d,a(np,np),TINY
        parameter (nmax=640,Tiny=1.0e-20)
        integer i,imax,j,k
        real aamax,dum,sum,vv(NMAX)
        d=1.
        do i=1,n
        aamax=0.
        do j=1,n
        if(abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        end do
        if(aamax.eq.0) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
        end do

        do j=1,n
        do i=1,j-1
        sum=a(i,j)
        do k=1,i-1
        sum=sum-a(i,k)*a(k,j)
        end do
        a(i,j)=sum
        end do
        aamax=0.
        do i=j,n
        sum=a(i,j)
        do k=1,j-1
        sum=sum-a(i,k)*a(k,j)
        enddo
        a(i,j)=sum
        dum=vv(i)*abs(sum)
        if(dum.ge.aamax)then
        imax=i
        aamax=dum
        endif
        enddo
        if(j.ne.imax)then
        do k=1,n
        dum=a(imax,k)
        a(imax,k)=a(j,k)
        a(j,k)=dum
        enddo
        d=-d
        vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=tiny
        if(j.ne.n)then
        dum=1./a(j,j)
        do i=j+1,n
        a(i,j)=a(i,j)*dum
        enddo
        endif
        enddo
        return
        end

subroutine lubksb(a,n,np,indx,b)
        integer n,np,indx(n)
        real a(np,np),b(n)
        integer i,ii,j,ll
        real sum
        ii=0
        do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if(ii.ne.0)then
        do j=ii,i-1
        sum=sum-a(i,j)*b(j)
        enddo
        else if(sum.ne.0.)then
        ii=i
        endif
        b(i)=sum
        enddo
        do i=n,1,-1
        sum=b(i)
        do j=i+1,n
        sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,i)
        enddo
        return
        end
