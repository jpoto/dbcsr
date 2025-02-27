!--------------------------------------------------------------------------------------------------!
! Copyright (C) by the DBCSR developers group - All rights reserved                                !
! This file is part of the DBCSR library.                                                          !
!                                                                                                  !
! For information on the license, see the LICENSE file.                                            !
! For further information please visit https://dbcsr.cp2k.org                                      !
! SPDX-License-Identifier: GPL-2.0+                                                                !
!--------------------------------------------------------------------------------------------------!

MODULE dbcsr_mp_operations
   !! Wrappers to message passing calls.

   USE dbcsr_config, ONLY: has_MPI
   USE dbcsr_data_methods, ONLY: dbcsr_data_get_type
   USE dbcsr_mp_methods, ONLY: &
      dbcsr_mp_get_process, dbcsr_mp_grid_setup, dbcsr_mp_group, dbcsr_mp_has_subgroups, &
      dbcsr_mp_my_col_group, dbcsr_mp_my_row_group, dbcsr_mp_mynode, dbcsr_mp_mypcol, &
      dbcsr_mp_myprow, dbcsr_mp_npcols, dbcsr_mp_nprows, dbcsr_mp_numnodes, dbcsr_mp_pgrid
   USE dbcsr_ptr_util, ONLY: memory_copy
   USE dbcsr_types, ONLY: dbcsr_data_obj, &
                          dbcsr_mp_obj, &
                          dbcsr_type_complex_4, &
                          dbcsr_type_complex_8, &
                          dbcsr_type_int_4, &
                          dbcsr_type_real_4, &
                          dbcsr_type_real_8
   USE dbcsr_kinds, ONLY: real_4, &
                          real_8
   USE dbcsr_mpiwrap, ONLY: &
      mp_allgather, mp_alltoall, mp_gatherv, mp_ibcast, mp_irecv, mp_iscatter, mp_isend, &
      mp_isendrecv, mp_rget, mp_sendrecv, mp_type_descriptor_type, mp_type_indexed_make_c, &
      mp_type_indexed_make_d, mp_type_indexed_make_r, mp_type_indexed_make_z, mp_type_make, &
      mp_waitall, mp_win_create, mp_comm_type, mp_request_type, mp_win_type
#include "base/dbcsr_base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mp_operations'

   ! MP routines
   PUBLIC :: hybrid_alltoall_s1, hybrid_alltoall_d1, &
             hybrid_alltoall_c1, hybrid_alltoall_z1, &
             hybrid_alltoall_i1, hybrid_alltoall_any
   PUBLIC :: dbcsr_allgatherv
   PUBLIC :: dbcsr_sendrecv_any
   PUBLIC :: dbcsr_isend_any, dbcsr_irecv_any
   PUBLIC :: dbcsr_win_create_any, dbcsr_rget_any, dbcsr_ibcast_any
   PUBLIC :: dbcsr_iscatterv_any, dbcsr_gatherv_any
   PUBLIC :: dbcsr_isendrecv_any
   ! Type helpers
   PUBLIC :: dbcsr_mp_type_from_anytype

   INTERFACE dbcsr_hybrid_alltoall
      MODULE PROCEDURE hybrid_alltoall_s1, hybrid_alltoall_d1, &
         hybrid_alltoall_c1, hybrid_alltoall_z1
      MODULE PROCEDURE hybrid_alltoall_i1
      MODULE PROCEDURE hybrid_alltoall_any
   END INTERFACE

CONTAINS

   SUBROUTINE hybrid_alltoall_any(sb, scount, sdispl, &
                                  rb, rcount, rdispl, mp_env, most_ptp, remainder_ptp, no_hybrid)
      TYPE(dbcsr_data_obj), INTENT(IN)                   :: sb
      INTEGER, DIMENSION(:), CONTIGUOUS, INTENT(IN)      :: scount, sdispl
      TYPE(dbcsr_data_obj), INTENT(INOUT)                :: rb
      INTEGER, DIMENSION(:), CONTIGUOUS, INTENT(IN)      :: rcount, rdispl
      TYPE(dbcsr_mp_obj), INTENT(IN)                     :: mp_env
      LOGICAL, INTENT(in), OPTIONAL                      :: most_ptp, remainder_ptp, no_hybrid

      CHARACTER(len=*), PARAMETER :: routineN = 'hybrid_alltoall_any'

      INTEGER                                            :: error_handle

!   ---------------------------------------------------------------------------

      CALL timeset(routineN, error_handle)

      SELECT CASE (dbcsr_data_get_type(sb))
      CASE (dbcsr_type_real_4)
         CALL hybrid_alltoall_s1(sb%d%r_sp, scount, sdispl, &
                                 rb%d%r_sp, rcount, rdispl, mp_env, &
                                 most_ptp, remainder_ptp, no_hybrid)
      CASE (dbcsr_type_real_8)
         CALL hybrid_alltoall_d1(sb%d%r_dp, scount, sdispl, &
                                 rb%d%r_dp, rcount, rdispl, mp_env, &
                                 most_ptp, remainder_ptp, no_hybrid)
      CASE (dbcsr_type_complex_4)
         CALL hybrid_alltoall_c1(sb%d%c_sp, scount, sdispl, &
                                 rb%d%c_sp, rcount, rdispl, mp_env, &
                                 most_ptp, remainder_ptp, no_hybrid)
      CASE (dbcsr_type_complex_8)
         CALL hybrid_alltoall_z1(sb%d%c_dp, scount, sdispl, &
                                 rb%d%c_dp, rcount, rdispl, mp_env, &
                                 most_ptp, remainder_ptp, no_hybrid)
      CASE default
         DBCSR_ABORT("Invalid data type")
      END SELECT

      CALL timestop(error_handle)
   END SUBROUTINE hybrid_alltoall_any

   SUBROUTINE hybrid_alltoall_i1(sb, scount, sdispl, &
                                 rb, rcount, rdispl, mp_env, most_ptp, remainder_ptp, no_hybrid)
      !! Row/column and global all-to-all
      !!
      !! Communicator selection
      !! Uses row and column communicators for row/column
      !! sends. Remaining sends are performed using the global
      !! communicator.  Point-to-point isend/irecv are used if ptp is
      !! set, otherwise a alltoall collective call is issued.
      !! see mp_alltoall

      INTEGER, DIMENSION(:), CONTIGUOUS, INTENT(in), &
         TARGET                                 :: sb
      INTEGER, DIMENSION(:), CONTIGUOUS, INTENT(IN) :: scount, sdispl
      INTEGER, DIMENSION(:), CONTIGUOUS, &
         INTENT(INOUT), TARGET                  :: rb
      INTEGER, DIMENSION(:), CONTIGUOUS, INTENT(IN) :: rcount, rdispl
      TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
         !! MP Environment
      LOGICAL, INTENT(IN), OPTIONAL            :: most_ptp, remainder_ptp, &
                                                  no_hybrid
         !! Use point-to-point for row/column; default is no
         !! Use point-to-point for remaining; default is no
         !! Use regular global collective; default is no

      INTEGER :: mynode, mypcol, myprow, nall_rr, nall_sr, ncol_rr, &
                 ncol_sr, npcols, nprows, nrow_rr, nrow_sr, numnodes, dst, src, &
                 prow, pcol, send_cnt, recv_cnt, tag, i
      INTEGER, ALLOCATABLE, DIMENSION(:) :: new_rcount, new_rdispl, new_scount, new_sdispl
      INTEGER, DIMENSION(:, :), POINTER        :: pgrid
      LOGICAL                                  :: most_collective, &
                                                  remainder_collective, no_h
      INTEGER, DIMENSION(:), POINTER, CONTIGUOUS :: send_data_p, recv_data_p
      TYPE(dbcsr_mp_obj)                       :: mpe
      TYPE(mp_comm_type)                       :: all_group, grp
      TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:) :: all_rr, all_sr, col_rr, col_sr, row_rr, row_sr

      IF (.NOT. dbcsr_mp_has_subgroups(mp_env)) THEN
         mpe = mp_env
         CALL dbcsr_mp_grid_setup(mpe)
      END IF
      most_collective = .TRUE.
      remainder_collective = .TRUE.
      no_h = .FALSE.
      IF (PRESENT(most_ptp)) most_collective = .NOT. most_ptp
      IF (PRESENT(remainder_ptp)) remainder_collective = .NOT. remainder_ptp
      IF (PRESENT(no_hybrid)) no_h = no_hybrid
      all_group = dbcsr_mp_group(mp_env)
      ! Don't use subcommunicators if they're not defined.
      no_h = no_h .OR. .NOT. dbcsr_mp_has_subgroups(mp_env) .OR. .NOT. has_MPI
      subgrouped: IF (mp_env%mp%subgroups_defined .AND. .NOT. no_h) THEN
         mynode = dbcsr_mp_mynode(mp_env)
         numnodes = dbcsr_mp_numnodes(mp_env)
         nprows = dbcsr_mp_nprows(mp_env)
         npcols = dbcsr_mp_npcols(mp_env)
         myprow = dbcsr_mp_myprow(mp_env)
         mypcol = dbcsr_mp_mypcol(mp_env)
         pgrid => dbcsr_mp_pgrid(mp_env)
         ALLOCATE (row_sr(0:npcols - 1)); nrow_sr = 0
         ALLOCATE (row_rr(0:npcols - 1)); nrow_rr = 0
         ALLOCATE (col_sr(0:nprows - 1)); ncol_sr = 0
         ALLOCATE (col_rr(0:nprows - 1)); ncol_rr = 0
         ALLOCATE (all_sr(0:numnodes - 1)); nall_sr = 0
         ALLOCATE (all_rr(0:numnodes - 1)); nall_rr = 0
         ALLOCATE (new_scount(numnodes), new_rcount(numnodes))
         ALLOCATE (new_sdispl(numnodes), new_rdispl(numnodes))
         IF (.NOT. remainder_collective) THEN
            CALL remainder_point_to_point()
         END IF
         IF (.NOT. most_collective) THEN
            CALL most_point_to_point()
         ELSE
            CALL most_alltoall()
         END IF
         IF (remainder_collective) THEN
            CALL remainder_alltoall()
         END IF
         ! Wait for all issued sends and receives.
         IF (.NOT. most_collective) THEN
            CALL mp_waitall(row_sr(0:nrow_sr - 1))
            CALL mp_waitall(col_sr(0:ncol_sr - 1))
            CALL mp_waitall(row_rr(0:nrow_rr - 1))
            CALL mp_waitall(col_rr(0:ncol_rr - 1))
         END IF
         IF (.NOT. remainder_collective) THEN
            CALL mp_waitall(all_sr(1:nall_sr))
            CALL mp_waitall(all_rr(1:nall_rr))
         END IF
      ELSE
         CALL mp_alltoall(sb, scount, sdispl, &
                          rb, rcount, rdispl, &
                          all_group)
      END IF subgrouped
   CONTAINS
      SUBROUTINE most_alltoall()
         DO pcol = 0, npcols - 1
            new_scount(1 + pcol) = scount(1 + pgrid(myprow, pcol))
            new_rcount(1 + pcol) = rcount(1 + pgrid(myprow, pcol))
            new_sdispl(1 + pcol) = sdispl(1 + pgrid(myprow, pcol))
            new_rdispl(1 + pcol) = rdispl(1 + pgrid(myprow, pcol))
         END DO
         CALL mp_alltoall(sb, new_scount(1:npcols), new_sdispl(1:npcols), &
                          rb, new_rcount(1:npcols), new_rdispl(1:npcols), &
                          dbcsr_mp_my_row_group(mp_env))
         DO prow = 0, nprows - 1
            new_scount(1 + prow) = scount(1 + pgrid(prow, mypcol))
            new_rcount(1 + prow) = rcount(1 + pgrid(prow, mypcol))
            new_sdispl(1 + prow) = sdispl(1 + pgrid(prow, mypcol))
            new_rdispl(1 + prow) = rdispl(1 + pgrid(prow, mypcol))
         END DO
         CALL mp_alltoall(sb, new_scount(1:nprows), new_sdispl(1:nprows), &
                          rb, new_rcount(1:nprows), new_rdispl(1:nprows), &
                          dbcsr_mp_my_col_group(mp_env))
      END SUBROUTINE most_alltoall
      SUBROUTINE most_point_to_point()
         ! Go through my prow and exchange.
         DO i = 0, npcols - 1
            pcol = MOD(mypcol + i, npcols)
            grp = dbcsr_mp_my_row_group(mp_env)
            !
            dst = dbcsr_mp_get_process(mp_env, myprow, pcol)
            send_cnt = scount(dst + 1)
            send_data_p => sb(1 + sdispl(dst + 1):1 + sdispl(dst + 1) + send_cnt - 1)
            tag = 4*mypcol
            IF (send_cnt .GT. 0) THEN
               CALL mp_isend(send_data_p, pcol, grp, row_sr(nrow_sr), tag)
               nrow_sr = nrow_sr + 1
            END IF
            !
            pcol = MODULO(mypcol - i, npcols)
            src = dbcsr_mp_get_process(mp_env, myprow, pcol)
            recv_cnt = rcount(src + 1)
            recv_data_p => rb(1 + rdispl(src + 1):1 + rdispl(src + 1) + recv_cnt - 1)
            tag = 4*pcol
            IF (recv_cnt .GT. 0) THEN
               CALL mp_irecv(recv_data_p, pcol, grp, row_rr(nrow_rr), tag)
               nrow_rr = nrow_rr + 1
            END IF
         END DO
         ! go through my pcol and exchange
         DO i = 0, nprows - 1
            prow = MOD(myprow + i, nprows)
            grp = dbcsr_mp_my_col_group(mp_env)
            !
            dst = dbcsr_mp_get_process(mp_env, prow, mypcol)
            send_cnt = scount(dst + 1)
            IF (send_cnt .GT. 0) THEN
               send_data_p => sb(1 + sdispl(dst + 1):1 + sdispl(dst + 1) + send_cnt - 1)
               tag = 4*myprow + 1
               CALL mp_isend(send_data_p, prow, grp, col_sr(ncol_sr), tag)
               ncol_sr = ncol_sr + 1
            END IF
            !
            prow = MODULO(myprow - i, nprows)
            src = dbcsr_mp_get_process(mp_env, prow, mypcol)
            recv_cnt = rcount(src + 1)
            IF (recv_cnt .GT. 0) THEN
               recv_data_p => rb(1 + rdispl(src + 1):1 + rdispl(src + 1) + recv_cnt - 1)
               tag = 4*prow + 1
               CALL mp_irecv(recv_data_p, prow, grp, col_rr(ncol_rr), tag)
               ncol_rr = ncol_rr + 1
            END IF
         END DO
      END SUBROUTINE most_point_to_point
      SUBROUTINE remainder_alltoall()
         new_scount(:) = scount(:)
         new_rcount(:) = rcount(:)
         DO prow = 0, nprows - 1
            new_scount(1 + pgrid(prow, mypcol)) = 0
            new_rcount(1 + pgrid(prow, mypcol)) = 0
         END DO
         DO pcol = 0, npcols - 1
            new_scount(1 + pgrid(myprow, pcol)) = 0
            new_rcount(1 + pgrid(myprow, pcol)) = 0
         END DO
         CALL mp_alltoall(sb, new_scount, sdispl, &
                          rb, new_rcount, rdispl, all_group)
      END SUBROUTINE remainder_alltoall
      SUBROUTINE remainder_point_to_point()
         INTEGER                                            :: col, row

         DO row = 0, nprows - 1
            prow = MOD(row + myprow, nprows)
            IF (prow .EQ. myprow) CYCLE
            DO col = 0, npcols - 1
               pcol = MOD(col + mypcol, npcols)
               IF (pcol .EQ. mypcol) CYCLE
               dst = dbcsr_mp_get_process(mp_env, prow, pcol)
               send_cnt = scount(dst + 1)
               IF (send_cnt .GT. 0) THEN
                  tag = 4*mynode + 2
                  send_data_p => sb(1 + sdispl(dst + 1):1 + sdispl(dst + 1) + send_cnt - 1)
                  CALL mp_isend(send_data_p, dst, all_group, all_sr(nall_sr + 1), tag)
                  nall_sr = nall_sr + 1
               END IF
               !
               src = dbcsr_mp_get_process(mp_env, prow, pcol)
               recv_cnt = rcount(src + 1)
               IF (recv_cnt .GT. 0) THEN
                  recv_data_p => rb(1 + rdispl(src + 1):1 + rdispl(src + 1) + recv_cnt - 1)
                  tag = 4*src + 2
                  CALL mp_irecv(recv_data_p, src, all_group, all_rr(nall_rr + 1), tag)
                  nall_rr = nall_rr + 1
               END IF
            END DO
         END DO
      END SUBROUTINE remainder_point_to_point
   END SUBROUTINE hybrid_alltoall_i1

   FUNCTION dbcsr_mp_type_from_anytype(data_area) RESULT(mp_type)
      !! Creates an MPI combined type from the given anytype.

      TYPE(dbcsr_data_obj), INTENT(IN)                   :: data_area
         !! Data area of any type
      TYPE(mp_type_descriptor_type)                      :: mp_type
         !! Type descriptor

      SELECT CASE (data_area%d%data_type)
      CASE (dbcsr_type_int_4)
         mp_type = mp_type_make(data_area%d%i4)
      CASE (dbcsr_type_real_4)
         mp_type = mp_type_make(data_area%d%r_sp)
      CASE (dbcsr_type_real_8)
         mp_type = mp_type_make(data_area%d%r_dp)
      CASE (dbcsr_type_complex_4)
         mp_type = mp_type_make(data_area%d%c_sp)
      CASE (dbcsr_type_complex_8)
         mp_type = mp_type_make(data_area%d%c_dp)
      END SELECT
   END FUNCTION dbcsr_mp_type_from_anytype

   SUBROUTINE dbcsr_sendrecv_any(msgin, dest, msgout, source, comm)
      !! sendrecv of encapsulated data.
      !! @note see mp_sendrecv

      TYPE(dbcsr_data_obj), INTENT(IN)                   :: msgin
      INTEGER, INTENT(IN)                                :: dest
      TYPE(dbcsr_data_obj), INTENT(INOUT)                :: msgout
      INTEGER, INTENT(IN)                                :: source
      TYPE(mp_comm_type), INTENT(IN)                     :: comm

      IF (dbcsr_data_get_type(msgin) .NE. dbcsr_data_get_type(msgout)) &
         DBCSR_ABORT("Different data type for msgin and msgout")

      SELECT CASE (dbcsr_data_get_type(msgin))
      CASE (dbcsr_type_real_4)
         CALL mp_sendrecv(msgin%d%r_sp, dest, msgout%d%r_sp, source, comm)
      CASE (dbcsr_type_real_8)
         CALL mp_sendrecv(msgin%d%r_dp, dest, msgout%d%r_dp, source, comm)
      CASE (dbcsr_type_complex_4)
         CALL mp_sendrecv(msgin%d%c_sp, dest, msgout%d%c_sp, source, comm)
      CASE (dbcsr_type_complex_8)
         CALL mp_sendrecv(msgin%d%c_dp, dest, msgout%d%c_dp, source, comm)
      CASE default
         DBCSR_ABORT("Incorrect data type")
      END SELECT
   END SUBROUTINE dbcsr_sendrecv_any

   SUBROUTINE dbcsr_isend_any(msgin, dest, comm, request, tag)
      !! Non-blocking send of encapsulated data.
      !! @note see mp_isend_iv

      TYPE(dbcsr_data_obj), INTENT(IN)                   :: msgin
      INTEGER, INTENT(IN)                                :: dest
      TYPE(mp_comm_type), INTENT(IN)                     :: comm
      TYPE(mp_request_type), INTENT(OUT)                 :: request
      INTEGER, INTENT(IN), OPTIONAL                      :: tag

      SELECT CASE (dbcsr_data_get_type(msgin))
      CASE (dbcsr_type_real_4)
         CALL mp_isend(msgin%d%r_sp, dest, comm, request, tag)
      CASE (dbcsr_type_real_8)
         CALL mp_isend(msgin%d%r_dp, dest, comm, request, tag)
      CASE (dbcsr_type_complex_4)
         CALL mp_isend(msgin%d%c_sp, dest, comm, request, tag)
      CASE (dbcsr_type_complex_8)
         CALL mp_isend(msgin%d%c_dp, dest, comm, request, tag)
      CASE default
         DBCSR_ABORT("Incorrect data type")
      END SELECT
   END SUBROUTINE dbcsr_isend_any

   SUBROUTINE dbcsr_irecv_any(msgin, source, comm, request, tag)
      !! Non-blocking recv of encapsulated data.
      !! @note see mp_irecv_iv

      TYPE(dbcsr_data_obj), INTENT(IN)                   :: msgin
      INTEGER, INTENT(IN)                                :: source
      TYPE(mp_comm_type), INTENT(IN)                     :: comm
      TYPE(mp_request_type), INTENT(OUT)                 :: request
      INTEGER, INTENT(IN), OPTIONAL                      :: tag

      SELECT CASE (dbcsr_data_get_type(msgin))
      CASE (dbcsr_type_real_4)
         CALL mp_irecv(msgin%d%r_sp, source, comm, request, tag)
      CASE (dbcsr_type_real_8)
         CALL mp_irecv(msgin%d%r_dp, source, comm, request, tag)
      CASE (dbcsr_type_complex_4)
         CALL mp_irecv(msgin%d%c_sp, source, comm, request, tag)
      CASE (dbcsr_type_complex_8)
         CALL mp_irecv(msgin%d%c_dp, source, comm, request, tag)
      CASE default
         DBCSR_ABORT("Incorrect data type")
      END SELECT
   END SUBROUTINE dbcsr_irecv_any

   SUBROUTINE dbcsr_win_create_any(base, comm, win)
      !! Window initialization function of encapsulated data.
      TYPE(dbcsr_data_obj), INTENT(IN)                   :: base
      TYPE(mp_comm_type), INTENT(IN)                     :: comm
      TYPE(mp_win_type), INTENT(OUT)                     :: win

      SELECT CASE (dbcsr_data_get_type(base))
      CASE (dbcsr_type_real_4)
         CALL mp_win_create(base%d%r_sp, comm, win)
      CASE (dbcsr_type_real_8)
         CALL mp_win_create(base%d%r_dp, comm, win)
      CASE (dbcsr_type_complex_4)
         CALL mp_win_create(base%d%c_sp, comm, win)
      CASE (dbcsr_type_complex_8)
         CALL mp_win_create(base%d%c_dp, comm, win)
      CASE default
         DBCSR_ABORT("Incorrect data type")
      END SELECT
   END SUBROUTINE dbcsr_win_create_any

   SUBROUTINE dbcsr_rget_any(base, source, win, win_data, myproc, disp, request, &
      !! Single-sided Get function of encapsulated data.
                             origin_datatype, target_datatype)
      TYPE(dbcsr_data_obj), INTENT(IN)                   :: base
      INTEGER, INTENT(IN)                                :: source
      TYPE(mp_win_type), INTENT(IN)                      :: win
      TYPE(dbcsr_data_obj), INTENT(IN)                   :: win_data
      INTEGER, INTENT(IN), OPTIONAL                      :: myproc, disp
      TYPE(mp_request_type), INTENT(OUT)                 :: request
      TYPE(mp_type_descriptor_type), INTENT(IN), &
         OPTIONAL                                        :: origin_datatype, target_datatype

      IF (dbcsr_data_get_type(base) /= dbcsr_data_get_type(win_data)) &
         DBCSR_ABORT("Mismatch data type between buffer and window")

      SELECT CASE (dbcsr_data_get_type(base))
      CASE (dbcsr_type_real_4)
         CALL mp_rget(base%d%r_sp, source, win, win_data%d%r_sp, myproc, &
                      disp, request, origin_datatype, target_datatype)
      CASE (dbcsr_type_real_8)
         CALL mp_rget(base%d%r_dp, source, win, win_data%d%r_dp, myproc, &
                      disp, request, origin_datatype, target_datatype)
      CASE (dbcsr_type_complex_4)
         CALL mp_rget(base%d%c_sp, source, win, win_data%d%c_sp, myproc, &
                      disp, request, origin_datatype, target_datatype)
      CASE (dbcsr_type_complex_8)
         CALL mp_rget(base%d%c_dp, source, win, win_data%d%c_dp, myproc, &
                      disp, request, origin_datatype, target_datatype)
      CASE default
         DBCSR_ABORT("Incorrect data type")
      END SELECT
   END SUBROUTINE dbcsr_rget_any

   SUBROUTINE dbcsr_ibcast_any(base, source, grp, request)
      !! Bcast function of encapsulated data.
      TYPE(dbcsr_data_obj), INTENT(IN)                   :: base
      INTEGER, INTENT(IN)                                :: source
      TYPE(mp_comm_type), INTENT(IN)                     :: grp
      TYPE(mp_request_type), INTENT(INOUT)               :: request

      SELECT CASE (dbcsr_data_get_type(base))
      CASE (dbcsr_type_real_4)
         CALL mp_ibcast(base%d%r_sp, source, grp, request)
      CASE (dbcsr_type_real_8)
         CALL mp_ibcast(base%d%r_dp, source, grp, request)
      CASE (dbcsr_type_complex_4)
         CALL mp_ibcast(base%d%c_sp, source, grp, request)
      CASE (dbcsr_type_complex_8)
         CALL mp_ibcast(base%d%c_dp, source, grp, request)
      CASE default
         DBCSR_ABORT("Incorrect data type")
      END SELECT
   END SUBROUTINE dbcsr_ibcast_any

   SUBROUTINE dbcsr_iscatterv_any(base, counts, displs, msg, recvcount, root, grp, request)
      !! Scatter function of encapsulated data.
      TYPE(dbcsr_data_obj), INTENT(IN)                   :: base
      INTEGER, DIMENSION(:), INTENT(IN), CONTIGUOUS      :: counts, displs
      TYPE(dbcsr_data_obj), INTENT(INOUT)                :: msg
      INTEGER, INTENT(IN)                                :: recvcount, root
      TYPE(mp_comm_type), INTENT(IN)                     :: grp
      TYPE(mp_request_type), INTENT(INOUT)               :: request

      IF (dbcsr_data_get_type(base) .NE. dbcsr_data_get_type(msg)) &
         DBCSR_ABORT("Different data type for msgin and msgout")

      SELECT CASE (dbcsr_data_get_type(base))
      CASE (dbcsr_type_real_4)
         CALL mp_iscatter(base%d%r_sp, counts, displs, msg%d%r_sp, recvcount, root, grp, request)
      CASE (dbcsr_type_real_8)
         CALL mp_iscatter(base%d%r_dp, counts, displs, msg%d%r_dp, recvcount, root, grp, request)
      CASE (dbcsr_type_complex_4)
         CALL mp_iscatter(base%d%c_sp, counts, displs, msg%d%c_sp, recvcount, root, grp, request)
      CASE (dbcsr_type_complex_8)
         CALL mp_iscatter(base%d%c_dp, counts, displs, msg%d%c_dp, recvcount, root, grp, request)
      CASE default
         DBCSR_ABORT("Incorrect data type")
      END SELECT
   END SUBROUTINE dbcsr_iscatterv_any

   SUBROUTINE dbcsr_gatherv_any(base, ub_base, msg, counts, displs, root, grp)
      !! Gather function of encapsulated data.
      TYPE(dbcsr_data_obj), INTENT(IN)                   :: base
      INTEGER, INTENT(IN)                                :: ub_base
      TYPE(dbcsr_data_obj), INTENT(INOUT)                :: msg
      INTEGER, DIMENSION(:), CONTIGUOUS, INTENT(IN)      :: counts, displs
      INTEGER, INTENT(IN)                                :: root
      TYPE(mp_comm_type), INTENT(IN)                     :: grp

      IF (dbcsr_data_get_type(base) .NE. dbcsr_data_get_type(msg)) &
         DBCSR_ABORT("Different data type for msgin and msgout")

      SELECT CASE (dbcsr_data_get_type(base))
      CASE (dbcsr_type_real_4)
         CALL mp_gatherv(base%d%r_sp(:ub_base), msg%d%r_sp, counts, displs, root, grp)
      CASE (dbcsr_type_real_8)
         CALL mp_gatherv(base%d%r_dp(:ub_base), msg%d%r_dp, counts, displs, root, grp)
      CASE (dbcsr_type_complex_4)
         CALL mp_gatherv(base%d%c_sp(:ub_base), msg%d%c_sp, counts, displs, root, grp)
      CASE (dbcsr_type_complex_8)
         CALL mp_gatherv(base%d%c_dp(:ub_base), msg%d%c_dp, counts, displs, root, grp)
      CASE default
         DBCSR_ABORT("Incorrect data type")
      END SELECT
   END SUBROUTINE dbcsr_gatherv_any

   SUBROUTINE dbcsr_isendrecv_any(msgin, dest, msgout, source, grp, send_request, recv_request)
      !! Send/Recv function of encapsulated data.
      TYPE(dbcsr_data_obj), INTENT(IN)                   :: msgin
      INTEGER, INTENT(IN)                                :: dest
      TYPE(dbcsr_data_obj), INTENT(INOUT)                :: msgout
      INTEGER, INTENT(IN)                                :: source
      TYPE(mp_comm_type), INTENT(IN)                     :: grp
      TYPE(mp_request_type), INTENT(OUT)                 :: send_request, recv_request

      IF (dbcsr_data_get_type(msgin) .NE. dbcsr_data_get_type(msgout)) &
         DBCSR_ABORT("Different data type for msgin and msgout")

      SELECT CASE (dbcsr_data_get_type(msgin))
      CASE (dbcsr_type_real_4)
         CALL mp_isendrecv(msgin%d%r_sp, dest, &
                           msgout%d%r_sp, source, &
                           grp, send_request, recv_request)
      CASE (dbcsr_type_real_8)
         CALL mp_isendrecv(msgin%d%r_dp, dest, &
                           msgout%d%r_dp, source, &
                           grp, send_request, recv_request)
      CASE (dbcsr_type_complex_4)
         CALL mp_isendrecv(msgin%d%c_sp, dest, &
                           msgout%d%c_sp, source, &
                           grp, send_request, recv_request)
      CASE (dbcsr_type_complex_8)
         CALL mp_isendrecv(msgin%d%c_dp, dest, &
                           msgout%d%c_dp, source, &
                           grp, send_request, recv_request)
      CASE default
         DBCSR_ABORT("Incorrect data type")
      END SELECT
   END SUBROUTINE dbcsr_isendrecv_any

   SUBROUTINE dbcsr_allgatherv(send_data, scount, recv_data, recv_count, recv_displ, gid)
      !! Allgather of encapsulated data
      !! @note see mp_allgatherv_dv

      TYPE(dbcsr_data_obj), INTENT(IN)                   :: send_data
      INTEGER, INTENT(IN)                                :: scount
      TYPE(dbcsr_data_obj), INTENT(INOUT)                :: recv_data
      INTEGER, DIMENSION(:), CONTIGUOUS, INTENT(IN)      :: recv_count, recv_displ
      TYPE(mp_comm_type), INTENT(IN)                     :: gid

      IF (dbcsr_data_get_type(send_data) /= dbcsr_data_get_type(recv_data)) &
         DBCSR_ABORT("Data type mismatch")
      SELECT CASE (dbcsr_data_get_type(send_data))
      CASE (dbcsr_type_real_4)
         CALL mp_allgather(send_data%d%r_sp(1:scount), recv_data%d%r_sp, &
                           recv_count, recv_displ, gid)
      CASE (dbcsr_type_real_8)
         CALL mp_allgather(send_data%d%r_dp(1:scount), recv_data%d%r_dp, &
                           recv_count, recv_displ, gid)
      CASE (dbcsr_type_complex_4)
         CALL mp_allgather(send_data%d%c_sp(1:scount), recv_data%d%c_sp, &
                           recv_count, recv_displ, gid)
      CASE (dbcsr_type_complex_8)
         CALL mp_allgather(send_data%d%c_dp(1:scount), recv_data%d%c_dp, &
                           recv_count, recv_displ, gid)
      CASE default
         DBCSR_ABORT("Invalid data type")
      END SELECT
   END SUBROUTINE dbcsr_allgatherv

   #:include '../data/dbcsr.fypp'
   #:for n, nametype1, base1, prec1, kind1, type1, dkind1 in inst_params_float
      SUBROUTINE hybrid_alltoall_${nametype1}$1(sb, scount, sdispl, &
                                                rb, rcount, rdispl, mp_env, most_ptp, remainder_ptp, no_hybrid)
      !! Row/column and global all-to-all
      !!
      !! Communicator selection
      !! Uses row and column communicators for row/column
      !! sends. Remaining sends are performed using the global
      !! communicator.  Point-to-point isend/irecv are used if ptp is
      !! set, otherwise a alltoall collective call is issued.
      !! see mp_alltoall

         ${type1}$, DIMENSION(:), &
            CONTIGUOUS, INTENT(in), TARGET        :: sb
         INTEGER, DIMENSION(:), CONTIGUOUS, INTENT(IN) :: scount, sdispl
         ${type1}$, DIMENSION(:), &
            CONTIGUOUS, INTENT(INOUT), TARGET     :: rb
         INTEGER, DIMENSION(:), CONTIGUOUS, INTENT(IN) :: rcount, rdispl
         TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
         !! MP Environment
         LOGICAL, INTENT(in), OPTIONAL            :: most_ptp, remainder_ptp, &
                                                     no_hybrid
         !! Use point-to-point for row/column; default is no
         !! Use point-to-point for remaining; default is no
         !! Use regular global collective; default is no

         INTEGER :: mynode, mypcol, myprow, nall_rr, nall_sr, ncol_rr, &
                    ncol_sr, npcols, nprows, nrow_rr, nrow_sr, numnodes, dst, src, &
                    prow, pcol, send_cnt, recv_cnt, tag, i
         INTEGER, ALLOCATABLE, DIMENSION(:) :: new_rcount, new_rdispl, new_scount, new_sdispl
         INTEGER, DIMENSION(:, :), CONTIGUOUS, POINTER :: pgrid
         LOGICAL                                  :: most_collective, &
                                                     remainder_collective, no_h
         ${type1}$, DIMENSION(:), CONTIGUOUS, POINTER :: send_data_p, recv_data_p
         TYPE(dbcsr_mp_obj)                       :: mpe
         TYPE(mp_comm_type)                       :: all_group, grp
         TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:) :: all_rr, all_sr, col_rr, col_sr, row_rr, row_sr

         IF (.NOT. dbcsr_mp_has_subgroups(mp_env)) THEN
            mpe = mp_env
            CALL dbcsr_mp_grid_setup(mpe)
         END IF
         most_collective = .TRUE.
         remainder_collective = .TRUE.
         no_h = .FALSE.
         IF (PRESENT(most_ptp)) most_collective = .NOT. most_ptp
         IF (PRESENT(remainder_ptp)) remainder_collective = .NOT. remainder_ptp
         IF (PRESENT(no_hybrid)) no_h = no_hybrid
         all_group = dbcsr_mp_group(mp_env)
         ! Don't use subcommunicators if they're not defined.
         no_h = no_h .OR. .NOT. dbcsr_mp_has_subgroups(mp_env) .OR. .NOT. has_MPI
         subgrouped: IF (mp_env%mp%subgroups_defined .AND. .NOT. no_h) THEN
            mynode = dbcsr_mp_mynode(mp_env)
            numnodes = dbcsr_mp_numnodes(mp_env)
            nprows = dbcsr_mp_nprows(mp_env)
            npcols = dbcsr_mp_npcols(mp_env)
            myprow = dbcsr_mp_myprow(mp_env)
            mypcol = dbcsr_mp_mypcol(mp_env)
            pgrid => dbcsr_mp_pgrid(mp_env)
            ALLOCATE (row_sr(0:npcols - 1)); nrow_sr = 0
            ALLOCATE (row_rr(0:npcols - 1)); nrow_rr = 0
            ALLOCATE (col_sr(0:nprows - 1)); ncol_sr = 0
            ALLOCATE (col_rr(0:nprows - 1)); ncol_rr = 0
            ALLOCATE (all_sr(0:numnodes - 1)); nall_sr = 0
            ALLOCATE (all_rr(0:numnodes - 1)); nall_rr = 0
            ALLOCATE (new_scount(numnodes), new_rcount(numnodes))
            ALLOCATE (new_sdispl(numnodes), new_rdispl(numnodes))
            IF (.NOT. remainder_collective) THEN
               CALL remainder_point_to_point()
            END IF
            IF (.NOT. most_collective) THEN
               CALL most_point_to_point()
            ELSE
               CALL most_alltoall()
            END IF
            IF (remainder_collective) THEN
               CALL remainder_alltoall()
            END IF
            ! Wait for all issued sends and receives.
            IF (.NOT. most_collective) THEN
               CALL mp_waitall(row_sr(0:nrow_sr - 1))
               CALL mp_waitall(col_sr(0:ncol_sr - 1))
               CALL mp_waitall(row_rr(0:nrow_rr - 1))
               CALL mp_waitall(col_rr(0:ncol_rr - 1))
            END IF
            IF (.NOT. remainder_collective) THEN
               CALL mp_waitall(all_sr(1:nall_sr))
               CALL mp_waitall(all_rr(1:nall_rr))
            END IF
         ELSE
            CALL mp_alltoall(sb, scount, sdispl, &
                             rb, rcount, rdispl, &
                             all_group)
         END IF subgrouped
      CONTAINS
         SUBROUTINE most_alltoall()
            DO pcol = 0, npcols - 1
               new_scount(1 + pcol) = scount(1 + pgrid(myprow, pcol))
               new_rcount(1 + pcol) = rcount(1 + pgrid(myprow, pcol))
               new_sdispl(1 + pcol) = sdispl(1 + pgrid(myprow, pcol))
               new_rdispl(1 + pcol) = rdispl(1 + pgrid(myprow, pcol))
            END DO
            CALL mp_alltoall(sb, new_scount(1:npcols), new_sdispl(1:npcols), &
                             rb, new_rcount(1:npcols), new_rdispl(1:npcols), &
                             dbcsr_mp_my_row_group(mp_env))
            DO prow = 0, nprows - 1
               new_scount(1 + prow) = scount(1 + pgrid(prow, mypcol))
               new_rcount(1 + prow) = rcount(1 + pgrid(prow, mypcol))
               new_sdispl(1 + prow) = sdispl(1 + pgrid(prow, mypcol))
               new_rdispl(1 + prow) = rdispl(1 + pgrid(prow, mypcol))
            END DO
            CALL mp_alltoall(sb, new_scount(1:nprows), new_sdispl(1:nprows), &
                             rb, new_rcount(1:nprows), new_rdispl(1:nprows), &
                             dbcsr_mp_my_col_group(mp_env))
         END SUBROUTINE most_alltoall
         SUBROUTINE most_point_to_point()
            ! Go through my prow and exchange.
            DO i = 0, npcols - 1
               pcol = MOD(mypcol + i, npcols)
               grp = dbcsr_mp_my_row_group(mp_env)
               !
               dst = dbcsr_mp_get_process(mp_env, myprow, pcol)
               send_cnt = scount(dst + 1)
               IF (send_cnt .GT. 0) THEN
                  send_data_p => sb(1 + sdispl(dst + 1):1 + sdispl(dst + 1) + send_cnt - 1)
                  IF (pcol .NE. mypcol) THEN
                     tag = 4*mypcol
                     CALL mp_isend(send_data_p, pcol, grp, row_sr(nrow_sr), tag)
                     nrow_sr = nrow_sr + 1
                  END IF
               END IF
               !
               pcol = MODULO(mypcol - i, npcols)
               src = dbcsr_mp_get_process(mp_env, myprow, pcol)
               recv_cnt = rcount(src + 1)
               IF (recv_cnt .GT. 0) THEN
                  recv_data_p => rb(1 + rdispl(src + 1):1 + rdispl(src + 1) + recv_cnt - 1)
                  IF (pcol .NE. mypcol) THEN
                     tag = 4*pcol
                     CALL mp_irecv(recv_data_p, pcol, grp, row_rr(nrow_rr), tag)
                     nrow_rr = nrow_rr + 1
                  ELSE
                     CALL memory_copy(recv_data_p, send_data_p, recv_cnt)
                  END IF
               END IF
            END DO
            ! go through my pcol and exchange
            DO i = 0, nprows - 1
               prow = MOD(myprow + i, nprows)
               grp = dbcsr_mp_my_col_group(mp_env)
               !
               dst = dbcsr_mp_get_process(mp_env, prow, mypcol)
               send_cnt = scount(dst + 1)
               IF (send_cnt .GT. 0) THEN
                  send_data_p => sb(1 + sdispl(dst + 1):1 + sdispl(dst + 1) + send_cnt - 1)
                  IF (prow .NE. myprow) THEN
                     tag = 4*myprow + 1
                     CALL mp_isend(send_data_p, prow, grp, col_sr(ncol_sr), tag)
                     ncol_sr = ncol_sr + 1
                  END IF
               END IF
               !
               prow = MODULO(myprow - i, nprows)
               src = dbcsr_mp_get_process(mp_env, prow, mypcol)
               recv_cnt = rcount(src + 1)
               IF (recv_cnt .GT. 0) THEN
                  recv_data_p => rb(1 + rdispl(src + 1):1 + rdispl(src + 1) + recv_cnt - 1)
                  IF (prow .NE. myprow) THEN
                     tag = 4*prow + 1
                     CALL mp_irecv(recv_data_p, prow, grp, col_rr(ncol_rr), tag)
                     ncol_rr = ncol_rr + 1
                  ELSE
                     CALL memory_copy(recv_data_p, send_data_p, recv_cnt)
                  END IF
               END IF
            END DO
         END SUBROUTINE most_point_to_point
         SUBROUTINE remainder_alltoall()
            new_scount(:) = scount(:)
            new_rcount(:) = rcount(:)
            DO prow = 0, nprows - 1
               new_scount(1 + pgrid(prow, mypcol)) = 0
               new_rcount(1 + pgrid(prow, mypcol)) = 0
            END DO
            DO pcol = 0, npcols - 1
               new_scount(1 + pgrid(myprow, pcol)) = 0
               new_rcount(1 + pgrid(myprow, pcol)) = 0
            END DO
            CALL mp_alltoall(sb, new_scount, sdispl, &
                             rb, new_rcount, rdispl, all_group)
         END SUBROUTINE remainder_alltoall
         SUBROUTINE remainder_point_to_point()
            INTEGER                                  :: col, row

            DO row = 0, nprows - 1
               prow = MOD(row + myprow, nprows)
               IF (prow .EQ. myprow) CYCLE
               DO col = 0, npcols - 1
                  pcol = MOD(col + mypcol, npcols)
                  IF (pcol .EQ. mypcol) CYCLE
                  dst = dbcsr_mp_get_process(mp_env, prow, pcol)
                  send_cnt = scount(dst + 1)
                  IF (send_cnt .GT. 0) THEN
                     send_data_p => sb(1 + sdispl(dst + 1):1 + sdispl(dst + 1) + send_cnt - 1)
                     tag = 4*mynode + 2
                     CALL mp_isend(send_data_p, dst, all_group, all_sr(nall_sr + 1), tag)
                     nall_sr = nall_sr + 1
                  END IF
                  !
                  src = dbcsr_mp_get_process(mp_env, prow, pcol)
                  recv_cnt = rcount(src + 1)
                  IF (recv_cnt .GT. 0) THEN
                     recv_data_p => rb(1 + rdispl(src + 1):1 + rdispl(src + 1) + recv_cnt - 1)
                     tag = 4*src + 2
                     CALL mp_irecv(recv_data_p, src, all_group, all_rr(nall_rr + 1), tag)
                     nall_rr = nall_rr + 1
                  END IF
               END DO
            END DO
         END SUBROUTINE remainder_point_to_point
      END SUBROUTINE hybrid_alltoall_${nametype1}$1
   #:endfor

END MODULE dbcsr_mp_operations
