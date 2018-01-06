#line 1 "cunm22.f"
/* cunm22.f -- translated by f2c (version 20100827).
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

#line 1 "cunm22.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};

/* > \brief \b CUNM22 multiplies a general matrix by a banded unitary matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUNM22 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunm22.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunm22.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunm22.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*     SUBROUTINE CUNM22( SIDE, TRANS, M, N, N1, N2, Q, LDQ, C, LDC, */
/*    $                   WORK, LWORK, INFO ) */

/*     .. Scalar Arguments .. */
/*     CHARACTER          SIDE, TRANS */
/*     INTEGER            M, N, N1, N2, LDQ, LDC, LWORK, INFO */
/*     .. */
/*     .. Array Arguments .. */
/*     COMPLEX            Q( LDQ, * ), C( LDC, * ), WORK( * ) */
/*     .. */

/* > \par Purpose */
/*  ============ */
/* > */
/* > \verbatim */
/* > */
/* >  CUNM22 overwrites the general complex M-by-N matrix C with */
/* > */
/* >                  SIDE = 'L'     SIDE = 'R' */
/* >  TRANS = 'N':      Q * C          C * Q */
/* >  TRANS = 'C':      Q**H * C       C * Q**H */
/* > */
/* >  where Q is a complex unitary matrix of order NQ, with NQ = M if */
/* >  SIDE = 'L' and NQ = N if SIDE = 'R'. */
/* >  The unitary matrix Q processes a 2-by-2 block structure */
/* > */
/* >         [  Q11  Q12  ] */
/* >     Q = [            ] */
/* >         [  Q21  Q22  ], */
/* > */
/* >  where Q12 is an N1-by-N1 lower triangular matrix and Q21 is an */
/* >  N2-by-N2 upper triangular matrix. */
/* > \endverbatim */

/*  Arguments */
/*  ========= */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q**H from the Left; */
/* >          = 'R': apply Q or Q**H from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  apply Q (No transpose); */
/* >          = 'C':  apply Q**H (Conjugate transpose). */
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
/* >          Q is COMPLEX array, dimension */
/* >                              (LDQ,M) if SIDE = 'L' */
/* >                              (LDQ,N) if SIDE = 'R' */
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
/* >          C is COMPLEX array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q. */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
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
/* Subroutine */ int cunm22_(char *side, char *trans, integer *m, integer *n, 
	integer *n1, integer *n2, doublecomplex *q, integer *ldq, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork,
	 integer *info, ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer q_dim1, q_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, nb, nq, nw, len;
    static logical left;
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ctrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    clacpy_(char *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static logical notran;
    static integer ldwork, lwkopt;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.6.0) -- */
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

#line 204 "cunm22.f"
    /* Parameter adjustments */
#line 204 "cunm22.f"
    q_dim1 = *ldq;
#line 204 "cunm22.f"
    q_offset = 1 + q_dim1;
#line 204 "cunm22.f"
    q -= q_offset;
#line 204 "cunm22.f"
    c_dim1 = *ldc;
#line 204 "cunm22.f"
    c_offset = 1 + c_dim1;
#line 204 "cunm22.f"
    c__ -= c_offset;
#line 204 "cunm22.f"
    --work;
#line 204 "cunm22.f"

#line 204 "cunm22.f"
    /* Function Body */
#line 204 "cunm22.f"
    *info = 0;
#line 205 "cunm22.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 206 "cunm22.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 207 "cunm22.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q; */
/*     NW is the minimum dimension of WORK. */

#line 212 "cunm22.f"
    if (left) {
#line 213 "cunm22.f"
	nq = *m;
#line 214 "cunm22.f"
    } else {
#line 215 "cunm22.f"
	nq = *n;
#line 216 "cunm22.f"
    }
#line 217 "cunm22.f"
    nw = nq;
#line 218 "cunm22.f"
    if (*n1 == 0 || *n2 == 0) {
#line 218 "cunm22.f"
	nw = 1;
#line 218 "cunm22.f"
    }
#line 219 "cunm22.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 220 "cunm22.f"
	*info = -1;
#line 221 "cunm22.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "C", (ftnlen)1, (ftnlen)1)) {
#line 223 "cunm22.f"
	*info = -2;
#line 224 "cunm22.f"
    } else if (*m < 0) {
#line 225 "cunm22.f"
	*info = -3;
#line 226 "cunm22.f"
    } else if (*n < 0) {
#line 227 "cunm22.f"
	*info = -4;
#line 228 "cunm22.f"
    } else if (*n1 < 0 || *n1 + *n2 != nq) {
#line 229 "cunm22.f"
	*info = -5;
#line 230 "cunm22.f"
    } else if (*n2 < 0) {
#line 231 "cunm22.f"
	*info = -6;
#line 232 "cunm22.f"
    } else if (*ldq < max(1,nq)) {
#line 233 "cunm22.f"
	*info = -8;
#line 234 "cunm22.f"
    } else if (*ldc < max(1,*m)) {
#line 235 "cunm22.f"
	*info = -10;
#line 236 "cunm22.f"
    } else if (*lwork < nw && ! lquery) {
#line 237 "cunm22.f"
	*info = -12;
#line 238 "cunm22.f"
    }

#line 240 "cunm22.f"
    if (*info == 0) {
#line 241 "cunm22.f"
	lwkopt = *m * *n;
#line 242 "cunm22.f"
	z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 242 "cunm22.f"
	work[1].r = z__1.r, work[1].i = z__1.i;
#line 243 "cunm22.f"
    }

#line 245 "cunm22.f"
    if (*info != 0) {
#line 246 "cunm22.f"
	i__1 = -(*info);
#line 246 "cunm22.f"
	xerbla_("CUNM22", &i__1, (ftnlen)6);
#line 247 "cunm22.f"
	return 0;
#line 248 "cunm22.f"
    } else if (lquery) {
#line 249 "cunm22.f"
	return 0;
#line 250 "cunm22.f"
    }

/*     Quick return if possible */

#line 254 "cunm22.f"
    if (*m == 0 || *n == 0) {
#line 255 "cunm22.f"
	work[1].r = 1., work[1].i = 0.;
#line 256 "cunm22.f"
	return 0;
#line 257 "cunm22.f"
    }

/*     Degenerate cases (N1 = 0 or N2 = 0) are handled using CTRMM. */

#line 261 "cunm22.f"
    if (*n1 == 0) {
#line 262 "cunm22.f"
	ctrmm_(side, "Upper", trans, "Non-Unit", m, n, &c_b1, &q[q_offset], 
		ldq, &c__[c_offset], ldc, (ftnlen)1, (ftnlen)5, (ftnlen)1, (
		ftnlen)8);
#line 264 "cunm22.f"
	work[1].r = 1., work[1].i = 0.;
#line 265 "cunm22.f"
	return 0;
#line 266 "cunm22.f"
    } else if (*n2 == 0) {
#line 267 "cunm22.f"
	ctrmm_(side, "Lower", trans, "Non-Unit", m, n, &c_b1, &q[q_offset], 
		ldq, &c__[c_offset], ldc, (ftnlen)1, (ftnlen)5, (ftnlen)1, (
		ftnlen)8);
#line 269 "cunm22.f"
	work[1].r = 1., work[1].i = 0.;
#line 270 "cunm22.f"
	return 0;
#line 271 "cunm22.f"
    }

/*     Compute the largest chunk size available from the workspace. */

/* Computing MAX */
#line 275 "cunm22.f"
    i__1 = 1, i__2 = min(*lwork,lwkopt) / nq;
#line 275 "cunm22.f"
    nb = max(i__1,i__2);

#line 277 "cunm22.f"
    if (left) {
#line 278 "cunm22.f"
	if (notran) {
#line 279 "cunm22.f"
	    i__1 = *n;
#line 279 "cunm22.f"
	    i__2 = nb;
#line 279 "cunm22.f"
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 280 "cunm22.f"
		i__3 = nb, i__4 = *n - i__ + 1;
#line 280 "cunm22.f"
		len = min(i__3,i__4);
#line 281 "cunm22.f"
		ldwork = *m;

/*              Multiply bottom part of C by Q12. */

#line 285 "cunm22.f"
		clacpy_("All", n1, &len, &c__[*n2 + 1 + i__ * c_dim1], ldc, &
			work[1], &ldwork, (ftnlen)3);
#line 287 "cunm22.f"
		ctrmm_("Left", "Lower", "No Transpose", "Non-Unit", n1, &len, 
			&c_b1, &q[(*n2 + 1) * q_dim1 + 1], ldq, &work[1], &
			ldwork, (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*              Multiply top part of C by Q11. */

#line 293 "cunm22.f"
		cgemm_("No Transpose", "No Transpose", n1, &len, n2, &c_b1, &
			q[q_offset], ldq, &c__[i__ * c_dim1 + 1], ldc, &c_b1, 
			&work[1], &ldwork, (ftnlen)12, (ftnlen)12);

/*              Multiply top part of C by Q21. */

#line 299 "cunm22.f"
		clacpy_("All", n2, &len, &c__[i__ * c_dim1 + 1], ldc, &work[*
			n1 + 1], &ldwork, (ftnlen)3);
#line 301 "cunm22.f"
		ctrmm_("Left", "Upper", "No Transpose", "Non-Unit", n2, &len, 
			&c_b1, &q[*n1 + 1 + q_dim1], ldq, &work[*n1 + 1], &
			ldwork, (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*              Multiply bottom part of C by Q22. */

#line 307 "cunm22.f"
		cgemm_("No Transpose", "No Transpose", n2, &len, n1, &c_b1, &
			q[*n1 + 1 + (*n2 + 1) * q_dim1], ldq, &c__[*n2 + 1 + 
			i__ * c_dim1], ldc, &c_b1, &work[*n1 + 1], &ldwork, (
			ftnlen)12, (ftnlen)12);

/*              Copy everything back. */

#line 313 "cunm22.f"
		clacpy_("All", m, &len, &work[1], &ldwork, &c__[i__ * c_dim1 
			+ 1], ldc, (ftnlen)3);
#line 315 "cunm22.f"
	    }
#line 316 "cunm22.f"
	} else {
#line 317 "cunm22.f"
	    i__2 = *n;
#line 317 "cunm22.f"
	    i__1 = nb;
#line 317 "cunm22.f"
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 318 "cunm22.f"
		i__3 = nb, i__4 = *n - i__ + 1;
#line 318 "cunm22.f"
		len = min(i__3,i__4);
#line 319 "cunm22.f"
		ldwork = *m;

/*              Multiply bottom part of C by Q21**H. */

#line 323 "cunm22.f"
		clacpy_("All", n2, &len, &c__[*n1 + 1 + i__ * c_dim1], ldc, &
			work[1], &ldwork, (ftnlen)3);
#line 325 "cunm22.f"
		ctrmm_("Left", "Upper", "Conjugate", "Non-Unit", n2, &len, &
			c_b1, &q[*n1 + 1 + q_dim1], ldq, &work[1], &ldwork, (
			ftnlen)4, (ftnlen)5, (ftnlen)9, (ftnlen)8);

/*              Multiply top part of C by Q11**H. */

#line 331 "cunm22.f"
		cgemm_("Conjugate", "No Transpose", n2, &len, n1, &c_b1, &q[
			q_offset], ldq, &c__[i__ * c_dim1 + 1], ldc, &c_b1, &
			work[1], &ldwork, (ftnlen)9, (ftnlen)12);

/*              Multiply top part of C by Q12**H. */

#line 337 "cunm22.f"
		clacpy_("All", n1, &len, &c__[i__ * c_dim1 + 1], ldc, &work[*
			n2 + 1], &ldwork, (ftnlen)3);
#line 339 "cunm22.f"
		ctrmm_("Left", "Lower", "Conjugate", "Non-Unit", n1, &len, &
			c_b1, &q[(*n2 + 1) * q_dim1 + 1], ldq, &work[*n2 + 1],
			 &ldwork, (ftnlen)4, (ftnlen)5, (ftnlen)9, (ftnlen)8);

/*              Multiply bottom part of C by Q22**H. */

#line 345 "cunm22.f"
		cgemm_("Conjugate", "No Transpose", n1, &len, n2, &c_b1, &q[*
			n1 + 1 + (*n2 + 1) * q_dim1], ldq, &c__[*n1 + 1 + i__ 
			* c_dim1], ldc, &c_b1, &work[*n2 + 1], &ldwork, (
			ftnlen)9, (ftnlen)12);

/*              Copy everything back. */

#line 351 "cunm22.f"
		clacpy_("All", m, &len, &work[1], &ldwork, &c__[i__ * c_dim1 
			+ 1], ldc, (ftnlen)3);
#line 353 "cunm22.f"
	    }
#line 354 "cunm22.f"
	}
#line 355 "cunm22.f"
    } else {
#line 356 "cunm22.f"
	if (notran) {
#line 357 "cunm22.f"
	    i__1 = *m;
#line 357 "cunm22.f"
	    i__2 = nb;
#line 357 "cunm22.f"
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 358 "cunm22.f"
		i__3 = nb, i__4 = *m - i__ + 1;
#line 358 "cunm22.f"
		len = min(i__3,i__4);
#line 359 "cunm22.f"
		ldwork = len;

/*              Multiply right part of C by Q21. */

#line 363 "cunm22.f"
		clacpy_("All", &len, n2, &c__[i__ + (*n1 + 1) * c_dim1], ldc, 
			&work[1], &ldwork, (ftnlen)3);
#line 365 "cunm22.f"
		ctrmm_("Right", "Upper", "No Transpose", "Non-Unit", &len, n2,
			 &c_b1, &q[*n1 + 1 + q_dim1], ldq, &work[1], &ldwork, 
			(ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*              Multiply left part of C by Q11. */

#line 371 "cunm22.f"
		cgemm_("No Transpose", "No Transpose", &len, n2, n1, &c_b1, &
			c__[i__ + c_dim1], ldc, &q[q_offset], ldq, &c_b1, &
			work[1], &ldwork, (ftnlen)12, (ftnlen)12);

/*              Multiply left part of C by Q12. */

#line 377 "cunm22.f"
		clacpy_("All", &len, n1, &c__[i__ + c_dim1], ldc, &work[*n2 * 
			ldwork + 1], &ldwork, (ftnlen)3);
#line 379 "cunm22.f"
		ctrmm_("Right", "Lower", "No Transpose", "Non-Unit", &len, n1,
			 &c_b1, &q[(*n2 + 1) * q_dim1 + 1], ldq, &work[*n2 * 
			ldwork + 1], &ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			12, (ftnlen)8);

/*              Multiply right part of C by Q22. */

#line 385 "cunm22.f"
		cgemm_("No Transpose", "No Transpose", &len, n1, n2, &c_b1, &
			c__[i__ + (*n1 + 1) * c_dim1], ldc, &q[*n1 + 1 + (*n2 
			+ 1) * q_dim1], ldq, &c_b1, &work[*n2 * ldwork + 1], &
			ldwork, (ftnlen)12, (ftnlen)12);

/*              Copy everything back. */

#line 391 "cunm22.f"
		clacpy_("All", &len, n, &work[1], &ldwork, &c__[i__ + c_dim1],
			 ldc, (ftnlen)3);
#line 393 "cunm22.f"
	    }
#line 394 "cunm22.f"
	} else {
#line 395 "cunm22.f"
	    i__2 = *m;
#line 395 "cunm22.f"
	    i__1 = nb;
#line 395 "cunm22.f"
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 396 "cunm22.f"
		i__3 = nb, i__4 = *m - i__ + 1;
#line 396 "cunm22.f"
		len = min(i__3,i__4);
#line 397 "cunm22.f"
		ldwork = len;

/*              Multiply right part of C by Q12**H. */

#line 401 "cunm22.f"
		clacpy_("All", &len, n1, &c__[i__ + (*n2 + 1) * c_dim1], ldc, 
			&work[1], &ldwork, (ftnlen)3);
#line 403 "cunm22.f"
		ctrmm_("Right", "Lower", "Conjugate", "Non-Unit", &len, n1, &
			c_b1, &q[(*n2 + 1) * q_dim1 + 1], ldq, &work[1], &
			ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)8);

/*              Multiply left part of C by Q11**H. */

#line 409 "cunm22.f"
		cgemm_("No Transpose", "Conjugate", &len, n1, n2, &c_b1, &c__[
			i__ + c_dim1], ldc, &q[q_offset], ldq, &c_b1, &work[1]
			, &ldwork, (ftnlen)12, (ftnlen)9);

/*              Multiply left part of C by Q21**H. */

#line 415 "cunm22.f"
		clacpy_("All", &len, n2, &c__[i__ + c_dim1], ldc, &work[*n1 * 
			ldwork + 1], &ldwork, (ftnlen)3);
#line 417 "cunm22.f"
		ctrmm_("Right", "Upper", "Conjugate", "Non-Unit", &len, n2, &
			c_b1, &q[*n1 + 1 + q_dim1], ldq, &work[*n1 * ldwork + 
			1], &ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)
			8);

/*              Multiply right part of C by Q22**H. */

#line 423 "cunm22.f"
		cgemm_("No Transpose", "Conjugate", &len, n2, n1, &c_b1, &c__[
			i__ + (*n2 + 1) * c_dim1], ldc, &q[*n1 + 1 + (*n2 + 1)
			 * q_dim1], ldq, &c_b1, &work[*n1 * ldwork + 1], &
			ldwork, (ftnlen)12, (ftnlen)9);

/*              Copy everything back. */

#line 429 "cunm22.f"
		clacpy_("All", &len, n, &work[1], &ldwork, &c__[i__ + c_dim1],
			 ldc, (ftnlen)3);
#line 431 "cunm22.f"
	    }
#line 432 "cunm22.f"
	}
#line 433 "cunm22.f"
    }

#line 435 "cunm22.f"
    z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 435 "cunm22.f"
    work[1].r = z__1.r, work[1].i = z__1.i;
#line 436 "cunm22.f"
    return 0;

/*     End of CUNM22 */

} /* cunm22_ */

