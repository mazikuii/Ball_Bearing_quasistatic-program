SUBROUTINE NEFREE(CTRL)
	USE NEPARA
	IMPLICIT NONE
	TYPE(NECTRL)	::	CTRL

	DEALLOCATE(CTRL%TYPX)
	DEALLOCATE(CTRL%TYPY)
ENDSUBROUTINE NEFREE

