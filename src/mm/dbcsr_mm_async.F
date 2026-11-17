!--------------------------------------------------------------------------------------------------!
! Copyright (C) by the DBCSR developers group - All rights reserved                                !
! This file is part of the DBCSR library.                                                          !
!                                                                                                  !
! For information on the license, see the LICENSE file.                                            !
! For further information please visit https://dbcsr.cp2k.org                                      !
! SPDX-License-Identifier: GPL-2.0+                                                                !
!--------------------------------------------------------------------------------------------------!

MODULE dbcsr_mm_async
   !! Layer A asynchronous multiplication front-end.
   !!
   !! A multiplication started here performs all DBCSR scheduling and enqueues
   !! the accelerator work on the calling thread. Only the final completion
   !! (stream synchronization, device-to-host merge and product finalization)
   !! is deferred until dbcsr_multiply_wait_prv() is called. Consequently the
   !! caller may not touch any matrix involved in the multiplication (and, to
   !! be safe, should not issue further DBCSR operations) before the wait.
   !!
   !! Several multiplications may be in flight at the same time; each occupies
   !! one slot returned through dbcsr_multiply_async_slot_alloc().
   USE dbcsr_kinds, ONLY: real_8, int_8
   USE dbcsr_data_methods, ONLY: dbcsr_data_flush_retained, dbcsr_data_set_retain_slot
   USE dbcsr_mm, ONLY: dbcsr_multiply_complete_async_slot, &
                       dbcsr_multiply_generic, &
                       dbcsr_multiply_async_slot_alloc, &
                       dbcsr_multiply_async_slot_pending
   USE dbcsr_types, ONLY: dbcsr_scalar_type, dbcsr_type
#include "base/dbcsr_base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_async'

   LOGICAL, SAVE :: async_enabled = .TRUE.

   PUBLIC :: dbcsr_multiply_async_prv, dbcsr_multiply_wait_prv
   PUBLIC :: dbcsr_multiply_pending_prv, dbcsr_multiply_cancel_prv
   PUBLIC :: dbcsr_async_set_enabled, dbcsr_async_get_enabled

CONTAINS

   SUBROUTINE dbcsr_async_set_enabled(enabled)
      LOGICAL, INTENT(IN)                            :: enabled

      async_enabled = enabled
   END SUBROUTINE dbcsr_async_set_enabled

   FUNCTION dbcsr_async_get_enabled() RESULT(enabled)
      LOGICAL                                        :: enabled

      enabled = async_enabled
   END FUNCTION dbcsr_async_get_enabled

   SUBROUTINE dbcsr_multiply_async_prv(transa, transb, alpha, matrix_a, matrix_b, beta, matrix_c, &
                                       first_row, last_row, first_column, last_column, first_k, last_k, &
                                       retain_sparsity, filter_eps, flop, slot)
      CHARACTER(LEN=1), INTENT(IN)                  :: transa, transb
      TYPE(dbcsr_scalar_type), INTENT(IN)           :: alpha
      TYPE(dbcsr_type), INTENT(IN)                  :: matrix_a, matrix_b
      TYPE(dbcsr_scalar_type), INTENT(IN)           :: beta
      TYPE(dbcsr_type), INTENT(INOUT)               :: matrix_c
      INTEGER, INTENT(IN), OPTIONAL                 :: first_row, last_row, first_column, last_column, &
                                                       first_k, last_k
      LOGICAL, INTENT(IN), OPTIONAL                 :: retain_sparsity
      REAL(real_8), INTENT(IN), OPTIONAL            :: filter_eps
      INTEGER(int_8), INTENT(OUT), OPTIONAL         :: flop
      INTEGER, INTENT(OUT)                          :: slot

      IF (PRESENT(flop)) &
         DBCSR_ABORT("dbcsr_multiply_async_prv: flop is not available for asynchronous multiplies")

      IF (.NOT. async_enabled) THEN
         CALL dbcsr_multiply_generic(transa, transb, alpha, matrix_a, matrix_b, beta, matrix_c, &
                                     first_row, last_row, first_column, last_column, first_k, last_k, &
                                     retain_sparsity, filter_eps=filter_eps)
         slot = -1
         RETURN
      END IF

      slot = dbcsr_multiply_async_slot_alloc()
      IF (slot .LT. 1) &
         DBCSR_ABORT("dbcsr_multiply_async_prv: no free asynchronous multiplication slot")

      CALL dbcsr_multiply_generic(transa, transb, alpha, matrix_a, matrix_b, beta, matrix_c, &
                                  first_row, last_row, first_column, last_column, first_k, last_k, &
                                  retain_sparsity, filter_eps=filter_eps, &
                                  defer_completion=.TRUE., async_slot=slot)
      CALL dbcsr_data_set_retain_slot(-1)

      ! The generic routine falls back to a synchronous completion for
      ! unsupported configurations; report "no handle" in that case.
      IF (.NOT. dbcsr_multiply_async_slot_pending(slot)) THEN
         CALL dbcsr_data_flush_retained(slot)
         slot = -1
      END IF
   END SUBROUTINE dbcsr_multiply_async_prv

   SUBROUTINE dbcsr_multiply_wait_prv(slot)
      INTEGER, INTENT(IN)                            :: slot

      CALL dbcsr_multiply_complete_async_slot(slot)
   END SUBROUTINE dbcsr_multiply_wait_prv

   SUBROUTINE dbcsr_multiply_cancel_prv(slot)
      !! A running multiplication cannot be interrupted safely, so cancellation
      !! simply waits for it to finish.
      INTEGER, INTENT(IN)                            :: slot

      CALL dbcsr_multiply_complete_async_slot(slot)
   END SUBROUTINE dbcsr_multiply_cancel_prv

   FUNCTION dbcsr_multiply_pending_prv(slot) RESULT(pending)
      INTEGER, INTENT(IN)                            :: slot
      LOGICAL                                        :: pending

      pending = dbcsr_multiply_async_slot_pending(slot)
   END FUNCTION dbcsr_multiply_pending_prv

END MODULE dbcsr_mm_async
