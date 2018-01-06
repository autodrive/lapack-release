#line 1 "dorm22.f"
/* dorm22.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

#line 1 "dorm22.f"
/* Table of constant values */

static doublereal c_b10 = 1.;

/* > \brief \b DORM22 multiplies a general matrix by a banded orthogonal matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORM22 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorm22.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorm22.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorm22.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*     SUBROUTINE DORM22( SIDE, TRANS, M, N, N1, N2, Q, LDQ, C, LDC, */
/*    $                   WORK, LWORK, INFO ) */

/*     .. Scalar Arguments .. */
/*     CHARACTER          SIDE, TRANS */
/*     INTEGER            M, N, N1, N2, LDQ, LDC, LWORK, INFO */
/*     .. */
/*     .. Array Arguments .. */
/*     DOUBLE PRECISION   Q( LDQ, * ), C( LDC, * ), WORK( * ) */
/*     .. */

/* > \par Purpose */
/*  ============ */
/* > */
/* > \verbatim */
/* > */
/* > */
/* >  DORM22 overwrites the general real M-by-N matrix C with */
/* > */
/* >                  SIDE = 'L'     SIDE = 'R' */
/* >  TRANS = 'N':      Q * C          C * Q */
/* >  TRANS = 'T':      Q**T * C       C * Q**T */
/* > */
/* >  where Q is a real orthogonal matrix of order NQ, with NQ = M if */
/* >  SIDE = 'L' and NQ = N if SIDE = 'R'. */
/* >  The orthogonal matrix Q processes a 2-by-2 block structure */
/* > */
/* >         [  Q11  Q12  ] */
/* >     Q = [            ] */
/* >         [  Q21  Q22  ], */
/* > */
/* >  where Q12 is an N1-by-N1 lower triangular matrix and Q21 is an */
/* >  N2-by-N2 upper triangular matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q**T from the Left; */
/* >          = 'R': apply Q or Q**T from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  apply Q (No transpose); */
/* >          = 'C':  apply Q**T (Conjugate transpose). */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix C. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix C. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N1 */
/* > \param[in] N2 */
/* > \verbatim */
/* >          N1 is INTEGER */
/* >          N2 is INTEGER */
/* >          The dimension of Q12 and Q21, respectively. N1, N2 >= 0. */
/* >          The following requirement must be satisfied: */
/* >          N1 + N2 = M if SIDE = 'L' and N1 + N2 = N if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] Q */
/* > \verbatim */
/* >          Q is DOUBLE PRECISION array, dimension */
/* >                                       (LDQ,M) if SIDE = 'L' */
/* >                                       (LDQ,N) if SIDE = 'R' */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q. */
/* >          LDQ >= max(1,M) if SIDE = 'L'; LDQ >= max(1,N) if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. LDC >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          If SIDE = 'L', LWORK >= max(1,N); */
/* >          if SIDE = 'R', LWORK >= max(1,M). */
/* >          For optimum performance LWORK >= M*N. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */


/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date January 2015 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dorm22_(char *side, char *trans, integer *m, integer *n, 
	integer *n1, integer *n2, doublereal *q, integer *ldq, doublereal *
	c__, integer *ldc, doublereal *work, integer *lwork, integer *info, 
	ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer q_dim1, q_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, nb, nq, nw, len;
    static logical left;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dlacpy_(
	    char *, integer *, integer *, doublereal *, integer *, doublereal 
	    *, integer *, ftnlen), xerbla_(char *, integer *, ftnlen);
    static logical notran;
    static integer ldwork, lwkopt;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     January 2015 */


/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 205 "dorm22.f"
    /* Parameter adjustments */
#line 205 "dorm22.f"
    q_dim1 = *ldq;
#line 205 "dorm22.f"
    q_offset = 1 + q_dim1;
#line 205 "dorm22.f"
    q -= q_offset;
#line 205 "dorm22.f"
    c_dim1 = *ldc;
#line 205 "dorm22.f"
    c_offset = 1 + c_dim1;
#line 205 "dorm22.f"
    c__ -= c_offset;
#line 205 "dorm22.f"
    --work;
#line 205 "dorm22.f"

#line 205 "dorm22.f"
    /* Function Body */
#line 205 "dorm22.f"
    *info = 0;
#line 206 "dorm22.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 207 "dorm22.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 208 "dorm22.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q; */
/*     NW is the minimum dimension of WORK. */

#line 213 "dorm22.f"
    if (left) {
#line 214 "dorm22.f"
	nq = *m;
#line 215 "dorm22.f"
    } else {
#line 216 "dorm22.f"
	nq = *n;
#line 217 "dorm22.f"
    }
#line 218 "dorm22.f"
    nw = nq;
#line 219 "dorm22.f"
    if (*n1 == 0 || *n2 == 0) {
#line 219 "dorm22.f"
	nw = 1;
#line 219 "dorm22.f"
    }
#line 220 "dorm22.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 221 "dorm22.f"
	*info = -1;
#line 222 "dorm22.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1)) {
#line 224 "dorm22.f"
	*info = -2;
#line 225 "dorm22.f"
    } else if (*m < 0) {
#line 226 "dorm22.f"
	*info = -3;
#line 227 "dorm22.f"
    } else if (*n < 0) {
#line 228 "dorm22.f"
	*info = -4;
#line 229 "dorm22.f"
    } else if (*n1 < 0 || *n1 + *n2 != nq) {
#line 230 "dorm22.f"
	*info = -5;
#line 231 "dorm22.f"
    } else if (*n2 < 0) {
#line 232 "dorm22.f"
	*info = -6;
#line 233 "dorm22.f"
    } else if (*ldq < max(1,nq)) {
#line 234 "dorm22.f"
	*info = -8;
#line 235 "dorm22.f"
    } else if (*ldc < max(1,*m)) {
#line 236 "dorm22.f"
	*info = -10;
#line 237 "dorm22.f"
    } else if (*lwork < nw && ! lquery) {
#line 238 "dorm22.f"
	*info = -12;
#line 239 "dorm22.f"
    }

#line 241 "dorm22.f"
    if (*info == 0) {
#line 242 "dorm22.f"
	lwkopt = *m * *n;
#line 243 "dorm22.f"
	work[1] = (doublereal) lwkopt;
#line 244 "dorm22.f"
    }

#line 246 "dorm22.f"
    if (*info != 0) {
#line 247 "dorm22.f"
	i__1 = -(*info);
#line 247 "dorm22.f"
	xerbla_("DORM22", &i__1, (ftnlen)6);
#line 248 "dorm22.f"
	return 0;
#line 249 "dorm22.f"
    } else if (lquery) {
#line 250 "dorm22.f"
	return 0;
#line 251 "dorm22.f"
    }

/*     Quick return if possible */

#line 255 "dorm22.f"
    if (*m == 0 || *n == 0) {
#line 256 "dorm22.f"
	work[1] = 1.;
#line 257 "dorm22.f"
	return 0;
#line 258 "dorm22.f"
    }

/*     Degenerate cases (N1 = 0 or N2 = 0) are handled using DTRMM. */

#line 262 "dorm22.f"
    if (*n1 == 0) {
#line 263 "dorm22.f"
	dtrmm_(side, "Upper", trans, "Non-Unit", m, n, &c_b10, &q[q_offset], 
		ldq, &c__[c_offset], ldc, (ftnlen)1, (ftnlen)5, (ftnlen)1, (
		ftnlen)8);
#line 265 "dorm22.f"
	work[1] = 1.;
#line 266 "dorm22.f"
	return 0;
#line 267 "dorm22.f"
    } else if (*n2 == 0) {
#line 268 "dorm22.f"
	dtrmm_(side, "Lower", trans, "Non-Unit", m, n, &c_b10, &q[q_offset], 
		ldq, &c__[c_offset], ldc, (ftnlen)1, (ftnlen)5, (ftnlen)1, (
		ftnlen)8);
#line 270 "dorm22.f"
	work[1] = 1.;
#line 271 "dorm22.f"
	return 0;
#line 272 "dorm22.f"
    }

/*     Compute the largest chunk size available from the workspace. */

/* Computing MAX */
#line 276 "dorm22.f"
    i__1 = 1, i__2 = min(*lwork,lwkopt) / nq;
#line 276 "dorm22.f"
    nb = max(i__1,i__2);

#line 278 "dorm22.f"
    if (left) {
#line 279 "dorm22.f"
	if (notran) {
#line 280 "dorm22.f"
	    i__1 = *n;
#line 280 "dorm22.f"
	    i__2 = nb;
#line 280 "dorm22.f"
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 281 "dorm22.f"
		i__3 = nb, i__4 = *n - i__ + 1;
#line 281 "dorm22.f"
		len = min(i__3,i__4);
#line 282 "dorm22.f"
		ldwork = *m;

/*              Multiply bottom part of C by Q12. */

#line 286 "dorm22.f"
		dlacpy_("All", n1, &len, &c__[*n2 + 1 + i__ * c_dim1], ldc, &
			work[1], &ldwork, (ftnlen)3);
#line 288 "dorm22.f"
		dtrmm_("Left", "Lower", "No Transpose", "Non-Unit", n1, &len, 
			&c_b10, &q[(*n2 + 1) * q_dim1 + 1], ldq, &work[1], &
			ldwork, (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*              Multiply top part of C by Q11. */

#line 294 "dorm22.f"
		dgemm_("No Transpose", "No Transpose", n1, &len, n2, &c_b10, &
			q[q_offset], ldq, &c__[i__ * c_dim1 + 1], ldc, &c_b10,
			 &work[1], &ldwork, (ftnlen)12, (ftnlen)12);

/*              Multiply top part of C by Q21. */

#line 300 "dorm22.f"
		dlacpy_("All", n2, &len, &c__[i__ * c_dim1 + 1], ldc, &work[*
			n1 + 1], &ldwork, (ftnlen)3);
#line 302 "dorm22.f"
		dtrmm_("Left", "Upper", "No Transpose", "Non-Unit", n2, &len, 
			&c_b10, &q[*n1 + 1 + q_dim1], ldq, &work[*n1 + 1], &
			ldwork, (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*              Multiply bottom part of C by Q22. */

#line 308 "dorm22.f"
		dgemm_("No Transpose", "No Transpose", n2, &len, n1, &c_b10, &
			q[*n1 + 1 + (*n2 + 1) * q_dim1], ldq, &c__[*n2 + 1 + 
			i__ * c_dim1], ldc, &c_b10, &work[*n1 + 1], &ldwork, (
			ftnlen)12, (ftnlen)12);

/*              Copy everything back. */

#line 314 "dorm22.f"
		dlacpy_("All", m, &len, &work[1], &ldwork, &c__[i__ * c_dim1 
			+ 1], ldc, (ftnlen)3);
#line 316 "dorm22.f"
	    }
#line 317 "dorm22.f"
	} else {
#line 318 "dorm22.f"
	    i__2 = *n;
#line 318 "dorm22.f"
	    i__1 = nb;
#line 318 "dorm22.f"
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 319 "dorm22.f"
		i__3 = nb, i__4 = *n - i__ + 1;
#line 319 "dorm22.f"
		len = min(i__3,i__4);
#line 320 "dorm22.f"
		ldwork = *m;

/*              Multiply bottom part of C by Q21**T. */

#line 324 "dorm22.f"
		dlacpy_("All", n2, &len, &c__[*n1 + 1 + i__ * c_dim1], ldc, &
			work[1], &ldwork, (ftnlen)3);
#line 326 "dorm22.f"
		dtrmm_("Left", "Upper", "Transpose", "Non-Unit", n2, &len, &
			c_b10, &q[*n1 + 1 + q_dim1], ldq, &work[1], &ldwork, (
			ftnlen)4, (ftnlen)5, (ftnlen)9, (ftnlen)8);

/*              Multiply top part of C by Q11**T. */

#line 332 "dorm22.f"
		dgemm_("Transpose", "No Transpose", n2, &len, n1, &c_b10, &q[
			q_offset], ldq, &c__[i__ * c_dim1 + 1], ldc, &c_b10, &
			work[1], &ldwork, (ftnlen)9, (ftnlen)12);

/*              Multiply top part of C by Q12**T. */

#line 338 "dorm22.f"
		dlacpy_("All", n1, &len, &c__[i__ * c_dim1 + 1], ldc, &work[*
			n2 + 1], &ldwork, (ftnlen)3);
#line 340 "dorm22.f"
		dtrmm_("Left", "Lower", "Transpose", "Non-Unit", n1, &len, &
			c_b10, &q[(*n2 + 1) * q_dim1 + 1], ldq, &work[*n2 + 1]
			, &ldwork, (ftnlen)4, (ftnlen)5, (ftnlen)9, (ftnlen)8)
			;

/*              Multiply bottom part of C by Q22**T. */

#line 346 "dorm22.f"
		dgemm_("Transpose", "No Transpose", n1, &len, n2, &c_b10, &q[*
			n1 + 1 + (*n2 + 1) * q_dim1], ldq, &c__[*n1 + 1 + i__ 
			* c_dim1], ldc, &c_b10, &work[*n2 + 1], &ldwork, (
			ftnlen)9, (ftnlen)12);

/*              Copy everything back. */

#line 352 "dorm22.f"
		dlacpy_("All", m, &len, &work[1], &ldwork, &c__[i__ * c_dim1 
			+ 1], ldc, (ftnlen)3);
#line 354 "dorm22.f"
	    }
#line 355 "dorm22.f"
	}
#line 356 "dorm22.f"
    } else {
#line 357 "dorm22.f"
	if (notran) {
#line 358 "dorm22.f"
	    i__1 = *m;
#line 358 "dorm22.f"
	    i__2 = nb;
#line 358 "dorm22.f"
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 359 "dorm22.f"
		i__3 = nb, i__4 = *m - i__ + 1;
#line 359 "dorm22.f"
		len = min(i__3,i__4);
#line 360 "dorm22.f"
		ldwork = len;

/*              Multiply right part of C by Q21. */

#line 364 "dorm22.f"
		dlacpy_("All", &len, n2, &c__[i__ + (*n1 + 1) * c_dim1], ldc, 
			&work[1], &ldwork, (ftnlen)3);
#line 366 "dorm22.f"
		dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", &len, n2,
			 &c_b10, &q[*n1 + 1 + q_dim1], ldq, &work[1], &ldwork,
			 (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*              Multiply left part of C by Q11. */

#line 372 "dorm22.f"
		dgemm_("No Transpose", "No Transpose", &len, n2, n1, &c_b10, &
			c__[i__ + c_dim1], ldc, &q[q_offset], ldq, &c_b10, &
			work[1], &ldwork, (ftnlen)12, (ftnlen)12);

/*              Multiply left part of C by Q12. */

#line 378 "dorm22.f"
		dlacpy_("All", &len, n1, &c__[i__ + c_dim1], ldc, &work[*n2 * 
			ldwork + 1], &ldwork, (ftnlen)3);
#line 380 "dorm22.f"
		dtrmm_("Right", "Lower", "No Transpose", "Non-Unit", &len, n1,
			 &c_b10, &q[(*n2 + 1) * q_dim1 + 1], ldq, &work[*n2 * 
			ldwork + 1], &ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			12, (ftnlen)8);

/*              Multiply right part of C by Q22. */

#line 386 "dorm22.f"
		dgemm_("No Transpose", "No Transpose", &len, n1, n2, &c_b10, &
			c__[i__ + (*n1 + 1) * c_dim1], ldc, &q[*n1 + 1 + (*n2 
			+ 1) * q_dim1], ldq, &c_b10, &work[*n2 * ldwork + 1], 
			&ldwork, (ftnlen)12, (ftnlen)12);

/*              Copy everything back. */

#line 392 "dorm22.f"
		dlacpy_("All", &len, n, &work[1], &ldwork, &c__[i__ + c_dim1],
			 ldc, (ftnlen)3);
#line 394 "dorm22.f"
	    }
#line 395 "dorm22.f"
	} else {
#line 396 "dorm22.f"
	    i__2 = *m;
#line 396 "dorm22.f"
	    i__1 = nb;
#line 396 "dorm22.f"
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 397 "dorm22.f"
		i__3 = nb, i__4 = *m - i__ + 1;
#line 397 "dorm22.f"
		len = min(i__3,i__4);
#line 398 "dorm22.f"
		ldwork = len;

/*              Multiply right part of C by Q12**T. */

#line 402 "dorm22.f"
		dlacpy_("All", &len, n1, &c__[i__ + (*n2 + 1) * c_dim1], ldc, 
			&work[1], &ldwork, (ftnlen)3);
#line 404 "dorm22.f"
		dtrmm_("Right", "Lower", "Transpose", "Non-Unit", &len, n1, &
			c_b10, &q[(*n2 + 1) * q_dim1 + 1], ldq, &work[1], &
			ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)8);

/*              Multiply left part of C by Q11**T. */

#line 410 "dorm22.f"
		dgemm_("No Transpose", "Transpose", &len, n1, n2, &c_b10, &
			c__[i__ + c_dim1], ldc, &q[q_offset], ldq, &c_b10, &
			work[1], &ldwork, (ftnlen)12, (ftnlen)9);

/*              Multiply left part of C by Q21**T. */

#line 416 "dorm22.f"
		dlacpy_("All", &len, n2, &c__[i__ + c_dim1], ldc, &work[*n1 * 
			ldwork + 1], &ldwork, (ftnlen)3);
#line 418 "dorm22.f"
		dtrmm_("Right", "Upper", "Transpose", "Non-Unit", &len, n2, &
			c_b10, &q[*n1 + 1 + q_dim1], ldq, &work[*n1 * ldwork 
			+ 1], &ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)9, (
			ftnlen)8);

/*              Multiply right part of C by Q22**T. */

#line 424 "dorm22.f"
		dgemm_("No Transpose", "Transpose", &len, n2, n1, &c_b10, &
			c__[i__ + (*n2 + 1) * c_dim1], ldc, &q[*n1 + 1 + (*n2 
			+ 1) * q_dim1], ldq, &c_b10, &work[*n1 * ldwork + 1], 
			&ldwork, (ftnlen)12, (ftnlen)9);

/*              Copy everything back. */

#line 430 "dorm22.f"
		dlacpy_("All", &len, n, &work[1], &ldwork, &c__[i__ + c_dim1],
			 ldc, (ftnlen)3);
#line 432 "dorm22.f"
	    }
#line 433 "dorm22.f"
	}
#line 434 "dorm22.f"
    }

#line 436 "dorm22.f"
    work[1] = (doublereal) lwkopt;
#line 437 "dorm22.f"
    return 0;

/*     End of DORM22 */

} /* dorm22_ */

