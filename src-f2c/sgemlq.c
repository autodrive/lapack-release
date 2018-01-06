#line 1 "sgemlq.f"
/* sgemlq.f -- translated by f2c (version 20100827).
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

#line 1 "sgemlq.f"

/*  Definition: */
/*  =========== */

/*      SUBROUTINE SGEMLQ( SIDE, TRANS, M, N, K, A, LDA, T, */
/*     $                   TSIZE, C, LDC, WORK, LWORK, INFO ) */


/*     .. Scalar Arguments .. */
/*      CHARACTER          SIDE, TRANS */
/*      INTEGER            INFO, LDA, M, N, K, LDT, TSIZE, LWORK, LDC */
/*     .. */
/*     .. Array Arguments .. */
/*      REAL               A( LDA, * ), T( * ), C(LDC, * ), WORK( * ) */
/*     .. */

/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >     SGEMLQ overwrites the general real M-by-N matrix C with */
/* > */
/* >                    SIDE = 'L'     SIDE = 'R' */
/* >    TRANS = 'N':      Q * C          C * Q */
/* >    TRANS = 'T':      Q**T * C       C * Q**T */
/* >    where Q is a real orthogonal matrix defined as the product */
/* >    of blocked elementary reflectors computed by short wide LQ */
/* >    factorization (SGELQ) */
/* > */
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
/* >          = 'N':  No transpose, apply Q; */
/* >          = 'T':  Transpose, apply Q**T. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >=0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix C. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of elementary reflectors whose product defines */
/* >          the matrix Q. */
/* >          If SIDE = 'L', M >= K >= 0; */
/* >          if SIDE = 'R', N >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension */
/* >                               (LDA,M) if SIDE = 'L', */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          Part of the data structure to represent Q as returned by DGELQ. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,K). */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is REAL array, dimension (MAX(5,TSIZE)). */
/* >          Part of the data structure to represent Q as returned by SGELQ. */
/* > \endverbatim */
/* > */
/* > \param[in] TSIZE */
/* > \verbatim */
/* >          TSIZE is INTEGER */
/* >          The dimension of the array T. TSIZE >= 5. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array, dimension (LDC,N) */
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
/* >         (workspace) REAL array, dimension (MAX(1,LWORK)) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          If LWORK = -1, then a workspace query is assumed. The routine */
/* >          only calculates the size of the WORK array, returns this */
/* >          value as WORK(1), and no error message related to WORK */
/* >          is issued by XERBLA. */
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

/* > \par Further Details */
/*  ==================== */
/* > */
/* > \verbatim */
/* > */
/* > These details are particular for this LAPACK implementation. Users should not */
/* > take them for granted. These details may change in the future, and are unlikely not */
/* > true for another LAPACK implementation. These details are relevant if one wants */
/* > to try to understand the code. They are not part of the interface. */
/* > */
/* > In this version, */
/* > */
/* >          T(2): row block size (MB) */
/* >          T(3): column block size (NB) */
/* >          T(6:TSIZE): data structure needed for Q, computed by */
/* >                           SLASWLQ or SGELQT */
/* > */
/* >  Depending on the matrix dimensions M and N, and row and column */
/* >  block sizes MB and NB returned by ILAENV, SGELQ will use either */
/* >  SLASWLQ (if the matrix is wide-and-short) or SGELQT to compute */
/* >  the LQ factorization. */
/* >  This version of SGEMLQ will use either SLAMSWLQ or SGEMLQT to */
/* >  multiply matrix Q by another matrix. */
/* >  Further Details in SLAMSWLQ or SGEMLQT. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sgemlq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *t, integer *
	tsize, doublereal *c__, integer *ldc, doublereal *work, integer *
	lwork, integer *info, ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int slamswlq_(char *, char *, integer *, integer *
	    , integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer mb, nb, mn, lw;
    static logical left, tran;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    static integer nblcks;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical notran, lquery;
    extern /* Subroutine */ int sgemlqt_(char *, char *, integer *, integer *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. */
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

#line 204 "sgemlq.f"
    /* Parameter adjustments */
#line 204 "sgemlq.f"
    a_dim1 = *lda;
#line 204 "sgemlq.f"
    a_offset = 1 + a_dim1;
#line 204 "sgemlq.f"
    a -= a_offset;
#line 204 "sgemlq.f"
    --t;
#line 204 "sgemlq.f"
    c_dim1 = *ldc;
#line 204 "sgemlq.f"
    c_offset = 1 + c_dim1;
#line 204 "sgemlq.f"
    c__ -= c_offset;
#line 204 "sgemlq.f"
    --work;
#line 204 "sgemlq.f"

#line 204 "sgemlq.f"
    /* Function Body */
#line 204 "sgemlq.f"
    lquery = *lwork == -1;
#line 205 "sgemlq.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 206 "sgemlq.f"
    tran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 207 "sgemlq.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 208 "sgemlq.f"
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);

#line 210 "sgemlq.f"
    mb = (integer) t[2];
#line 211 "sgemlq.f"
    nb = (integer) t[3];
#line 212 "sgemlq.f"
    if (left) {
#line 213 "sgemlq.f"
	lw = *n * mb;
#line 214 "sgemlq.f"
	mn = *m;
#line 215 "sgemlq.f"
    } else {
#line 216 "sgemlq.f"
	lw = *m * mb;
#line 217 "sgemlq.f"
	mn = *n;
#line 218 "sgemlq.f"
    }

#line 220 "sgemlq.f"
    if (nb > *k && mn > *k) {
#line 221 "sgemlq.f"
	if ((mn - *k) % (nb - *k) == 0) {
#line 222 "sgemlq.f"
	    nblcks = (mn - *k) / (nb - *k);
#line 223 "sgemlq.f"
	} else {
#line 224 "sgemlq.f"
	    nblcks = (mn - *k) / (nb - *k) + 1;
#line 225 "sgemlq.f"
	}
#line 226 "sgemlq.f"
    } else {
#line 227 "sgemlq.f"
	nblcks = 1;
#line 228 "sgemlq.f"
    }

#line 230 "sgemlq.f"
    *info = 0;
#line 231 "sgemlq.f"
    if (! left && ! right) {
#line 232 "sgemlq.f"
	*info = -1;
#line 233 "sgemlq.f"
    } else if (! tran && ! notran) {
#line 234 "sgemlq.f"
	*info = -2;
#line 235 "sgemlq.f"
    } else if (*m < 0) {
#line 236 "sgemlq.f"
	*info = -3;
#line 237 "sgemlq.f"
    } else if (*n < 0) {
#line 238 "sgemlq.f"
	*info = -4;
#line 239 "sgemlq.f"
    } else if (*k < 0 || *k > mn) {
#line 240 "sgemlq.f"
	*info = -5;
#line 241 "sgemlq.f"
    } else if (*lda < max(1,*k)) {
#line 242 "sgemlq.f"
	*info = -7;
#line 243 "sgemlq.f"
    } else if (*tsize < 5) {
#line 244 "sgemlq.f"
	*info = -9;
#line 245 "sgemlq.f"
    } else if (*ldc < max(1,*m)) {
#line 246 "sgemlq.f"
	*info = -11;
#line 247 "sgemlq.f"
    } else if (*lwork < max(1,lw) && ! lquery) {
#line 248 "sgemlq.f"
	*info = -13;
#line 249 "sgemlq.f"
    }

#line 251 "sgemlq.f"
    if (*info == 0) {
#line 252 "sgemlq.f"
	work[1] = (doublereal) lw;
#line 253 "sgemlq.f"
    }

#line 255 "sgemlq.f"
    if (*info != 0) {
#line 256 "sgemlq.f"
	i__1 = -(*info);
#line 256 "sgemlq.f"
	xerbla_("SGEMLQ", &i__1, (ftnlen)6);
#line 257 "sgemlq.f"
	return 0;
#line 258 "sgemlq.f"
    } else if (lquery) {
#line 259 "sgemlq.f"
	return 0;
#line 260 "sgemlq.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 264 "sgemlq.f"
    i__1 = min(*m,*n);
#line 264 "sgemlq.f"
    if (min(i__1,*k) == 0) {
#line 265 "sgemlq.f"
	return 0;
#line 266 "sgemlq.f"
    }

/* Computing MAX */
#line 268 "sgemlq.f"
    i__1 = max(*m,*n);
#line 268 "sgemlq.f"
    if (left && *m <= *k || right && *n <= *k || nb <= *k || nb >= max(i__1,*
	    k)) {
#line 270 "sgemlq.f"
	sgemlqt_(side, trans, m, n, k, &mb, &a[a_offset], lda, &t[6], &mb, &
		c__[c_offset], ldc, &work[1], info, (ftnlen)1, (ftnlen)1);
#line 272 "sgemlq.f"
    } else {
#line 273 "sgemlq.f"
	slamswlq_(side, trans, m, n, k, &mb, &nb, &a[a_offset], lda, &t[6], &
		mb, &c__[c_offset], ldc, &work[1], lwork, info, (ftnlen)1, (
		ftnlen)1);
#line 275 "sgemlq.f"
    }

#line 277 "sgemlq.f"
    work[1] = (doublereal) lw;

#line 279 "sgemlq.f"
    return 0;

/*     End of SGEMLQ */

} /* sgemlq_ */

