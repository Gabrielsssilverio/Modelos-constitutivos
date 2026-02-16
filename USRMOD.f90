! ==========================================================
! MODELO MOHR-COULOMB - VERSAO FINAL
! ==========================================================

SUBROUTINE USER_MOD(IDTask, iMod, IsUndr, iStep, iTer, Iel, Int, X, &
                    Y, Z, Time0, dTime, Props, Sigo, Swp0, StVaro, &
                    dEps, D, Bulk_W, Sig, Swp, StVar, ipl, nStat, &
                    NonSym, iStrsDep, iTimeDep, iTang, iPrjDir, &
                    iPrjLen, iAbort)
    IMPLICIT NONE
    !GCC$ ATTRIBUTES DLLEXPORT, STDCALL :: USER_MOD

    INTEGER, INTENT(IN) :: IDTask, iMod, IsUndr, iStep, iTer, Iel, Int, nStat, iPrjLen
    DOUBLE PRECISION, INTENT(IN) :: X, Y, Z, Time0, dTime, Bulk_W
    DOUBLE PRECISION, INTENT(INOUT) :: Props(50), Sigo(20), Swp0, StVaro(*) 
    DOUBLE PRECISION, INTENT(IN) :: dEps(12)     
    DOUBLE PRECISION, INTENT(INOUT) :: D(6,6), Sig(6), Swp, StVar(*)
    INTEGER, INTENT(INOUT) :: ipl, NonSym, iStrsDep, iTimeDep, iTang, iAbort
    CHARACTER(LEN=1), INTENT(IN) :: iPrjDir(*)

    DOUBLE PRECISION :: E, nu, G, fac, term1, term2, c, phi_rad
    DOUBLE PRECISION :: p_trial, q_trial, q_max, M, scaling, s(6)
    INTEGER :: i, j

    IF (IDTask == 2) THEN
        E = Props(1); nu = Props(2); c = Props(3)
        phi_rad = Props(4) * 3.1415926535D0 / 180.0D0
        
        DO i = 1, 6
            Sig(i) = Sigo(i)
            DO j = 1, 6
                Sig(i) = Sig(i) + D(i,j) * dEps(j)
            END DO
        END DO

        p_trial = -(Sig(1) + Sig(2) + Sig(3)) / 3.0D0
        s(1)=-Sig(1)-p_trial; s(2)=-Sig(2)-p_trial; s(3)=-Sig(3)-p_trial
        s(4)=-Sig(4); s(5)=-Sig(5); s(6)=-Sig(6)
        
        q_trial = SQRT(0.5D0*((s(1)-s(2))**2 + (s(2)-s(3))**2 + (s(3)-s(1))**2 + 6.0D0*(s(4)**2+s(5)**2+s(6)**2)))
        M = (6.0D0 * SIN(phi_rad)) / (3.0D0 - SIN(phi_rad))
        q_max = M * p_trial + (6.0D0 * c * COS(phi_rad)) / (3.0D0 - SIN(phi_rad))

        IF (p_trial < -1.0D-6) THEN
            Sig = 0.0D0
            ipl = 1
        ELSE IF (q_trial > q_max .AND. q_trial > 1.0D-6) THEN
            scaling = q_max / q_trial
            DO i = 1, 3
                Sig(i) = -(p_trial + s(i) * scaling)
            END DO
            DO i = 4, 6
                Sig(i) = -(s(i) * scaling)
            END DO
            ipl = 1 
        ELSE
            ipl = 0 
        END IF

    ELSEIF (IDTask == 3 .OR. IDTask == 6) THEN
        E = Props(1); nu = Props(2)
        G = 0.5D0 * E / (1.0D0 + nu)
        fac = 2.0D0 * G / (1.0D0 - 2.0D0 * nu)
        term1 = fac * (1.0D0 - nu); term2 = fac * nu
        D = 0.0D0
        D(1,1)=term1; D(2,2)=term1; D(3,3)=term1
        D(1,2)=term2; D(1,3)=term2; D(2,1)=term2; D(2,3)=term2; D(3,1)=term2; D(3,2)=term2
        D(4,4)=G; D(5,5)=G; D(6,6)=G

    ELSEIF (IDTask == 5) THEN
        NonSym = 1; iStrsDep = 0; iTimeDep = 0; iTang = 0
    END IF
END SUBROUTINE USER_MOD

SUBROUTINE GETMODELCOUNT(C)
    IMPLICIT NONE
    !GCC$ ATTRIBUTES DLLEXPORT, STDCALL :: GETMODELCOUNT
    INTEGER, INTENT(OUT) :: C
    C = 1
END SUBROUTINE GETMODELCOUNT

SUBROUTINE GETMODELNAME(iModel, Name)
    IMPLICIT NONE
    !GCC$ ATTRIBUTES DLLEXPORT, STDCALL :: GETMODELNAME
    INTEGER, INTENT(IN) :: iModel
    CHARACTER(LEN=256), INTENT(OUT) :: Name
    Name = CHAR(11) // "MC_UDSM"
END SUBROUTINE GETMODELNAME

SUBROUTINE GETPARAMCOUNT(iModel, C)
    IMPLICIT NONE
    !GCC$ ATTRIBUTES DLLEXPORT, STDCALL :: GETPARAMCOUNT
    INTEGER, INTENT(IN) :: iModel
    INTEGER, INTENT(OUT) :: C
    C = 4
END SUBROUTINE GETPARAMCOUNT

SUBROUTINE GETPARAMNAME(iModel, iParam, Name)
    IMPLICIT NONE
    !GCC$ ATTRIBUTES DLLEXPORT, STDCALL :: GETPARAMNAME
    INTEGER, INTENT(IN) :: iModel, iParam
    CHARACTER(LEN=256), INTENT(OUT) :: Name
    IF (iParam == 1) Name = CHAR(1) // "E"
    IF (iParam == 2) Name = CHAR(2) // "nu"
    IF (iParam == 3) Name = CHAR(1) // "c"
    IF (iParam == 4) Name = CHAR(3) // "phi"
END SUBROUTINE GETPARAMNAME

SUBROUTINE GETPARAMUNIT(iModel, iParam, Unit)
    IMPLICIT NONE
    !GCC$ ATTRIBUTES DLLEXPORT, STDCALL :: GETPARAMUNIT
    INTEGER, INTENT(IN) :: iModel, iParam
    CHARACTER(LEN=256), INTENT(OUT) :: Unit
    IF (iParam == 1 .OR. iParam == 3) Unit = CHAR(5) // "kN/m2"
    IF (iParam == 2 .OR. iParam == 4) Unit = CHAR(1) // "-"
END SUBROUTINE GETPARAMUNIT