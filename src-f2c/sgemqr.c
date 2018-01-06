#line 1 "sgemqr.f"
/* sgemqr.f -- translated by f2c (version 20100827).
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

#line 1 "sgemqr.f"

/*  Definition: */
/*  =========== */

/*      SUBROUTINE SGEMQR( SIDE, TRANS, M, N, K, A, LDA, T, */
/*     $                   TSIZE, C, LDC, WORK, LWORK, INFO ) */


/*     .. Scalar Arguments .. */
/*     CHARACTER         SIDE, TRANS */
/*     INTEGER           INFO, LDA, M, N, K, LDT, TSIZE, LWORK, LDC */
/*     .. */
/*     .. Array Arguments .. */
/*     REAL              A( LDA, * ), T( * ), C( LDC, * ), WORK( * ) */
/*     .. */

/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEMQR overwrites the general real M-by-N matrix C with */
/* > */
/* >                      SIDE = 'L'     SIDE = 'R' */
/* >     TRANS = 'N':      Q * C          C * Q */
/* >     TRANS = 'T':      Q**T * C       C * Q**T */
/* > */
/* > where Q is a real orthogonal matrix defined as the product */
/* > of blocked elementary reflectors computed by tall skinny */
/* > QR factorization (SGEQR) */
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
/* >          A is REAL array, dimension (LDA,K) */
/* >          Part of the data structure to represent Q as returned by SGEQR. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. */
/* >          If SIDE = 'L', LDA >= max(1,M); */
/* >          if SIDE = 'R', LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is REAL array, dimension (MAX(5,TSIZE)). */
/* >          Part of the data structure to represent Q as returned by SGEQR. */
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
/* >                           SLATSQR or SGEQRT */
/* > */
/* >  Depending on the matrix dimensions M and N, and row and column */
/* >  block sizes MB and NB returned by ILAENV, SGEQR will use either */
/* >  SLATSQR (if the matrix is tall-and-skinny) or SGEQRT to compute */
/* >  the QR factorization. */
/* >  This version of SGEMQR will use either SLAMTSQR or SGEMQRT to */
/* >  multiply matrix Q by another matrix. */
/* >  Further Details in SLAMTSQR or SGEMQRT. */
/* > */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sgemqr_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *t, integer *
	tsize, doublereal *c__, integer *ldc, doublereal *work, integer *
	lwork, integer *info, ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int slamtsqr_(char *, char *, integer *, integer *
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
    extern /* Subroutine */ int sgemqrt_(char *, char *, integer *, integer *,
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

#line 206 "sgemqr.f"
    /* Parameter adjustments */
#line 206 "sgemqr.f"
    a_dim1 = *lda;
#line 206 "sgemqr.f"
    a_offset = 1 + a_dim1;
#line 206 "sgemqr.f"
    a -= a_offset;
#line 206 "sgemqr.f"
    --t;
#line 206 "sgemqr.f"
    c_dim1 = *ldc;
#line 206 "sgemqr.f"
    c_offset = 1 + c_dim1;
#line 206 "sgemqr.f"
    c__ -= c_offset;
#line 206 "sgemqr.f"
    --work;
#line 206 "sgemqr.f"

#line 206 "sgemqr.f"
    /* Function Body */
#line 206 "sgemqr.f"
    lquery = *lwork == -1;
#line 207 "sgemqr.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 208 "sgemqr.f"
    tran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 209 "sgemqr.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 210 "sgemqr.f"
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);

#line 212 "sgemqr.f"
    mb = (integer) t[2];
#line 213 "sgemqr.f"
    nb = (integer) t[3];
#line 214 "sgemqr.f"
    if (left) {
#line 215 "sgemqr.f"
	lw = *n * nb;
#line 216 "sgemqr.f"
	mn = *m;
#line 217 "sgemqr.f"
    } else {
#line 218 "sgemqr.f"
	lw = mb * nb;
#line 219 "sgemqr.f"
	mn = *n;
#line 220 "sgemqr.f"
    }

#line 222 "sgemqr.f"
    if (mb > *k && mn > *k) {
#line 223 "sgemqr.f"
	if ((mn - *k) % (mb - *k) == 0) {
#line 224 "sgemqr.f"
	    nblcks = (mn - *k) / (mb - *k);
#line 225 "sgemqr.f"
	} else {
#line 226 "sgemqr.f"
	    nblcks = (mn - *k) / (mb - *k) + 1;
#line 227 "sgemqr.f"
	}
#line 228 "sgemqr.f"
    } else {
#line 229 "sgemqr.f"
	nblcks = 1;
#line 230 "sgemqr.f"
    }

#line 232 "sgemqr.f"
    *info = 0;
#line 233 "sgemqr.f"
    if (! left && ! right) {
#line 234 "sgemqr.f"
	*info = -1;
#line 235 "sgemqr.f"
    } else if (! tran && ! notran) {
#line 236 "sgemqr.f"
	*info = -2;
#line 237 "sgemqr.f"
    } else if (*m < 0) {
#line 238 "sgemqr.f"
	*info = -3;
#line 239 "sgemqr.f"
    } else if (*n < 0) {
#line 240 "sgemqr.f"
	*info = -4;
#line 241 "sgemqr.f"
    } else if (*k < 0 || *k > mn) {
#line 242 "sgemqr.f"
	*info = -5;
#line 243 "sgemqr.f"
    } else if (*lda < max(1,mn)) {
#line 244 "sgemqr.f"
	*info = -7;
#line 245 "sgemqr.f"
    } else if (*tsize < 5) {
#line 246 "sgemqr.f"
	*info = -9;
#line 247 "sgemqr.f"
    } else if (*ldc < max(1,*m)) {
#line 248 "sgemqr.f"
	*info = -11;
#line 249 "sgemqr.f"
    } else if (*lwork < max(1,lw) && ! lquery) {
#line 250 "sgemqr.f"
	*info = -13;
#line 251 "sgemqr.f"
    }

#line 253 "sgemqr.f"
    if (*info == 0) {
#line 254 "sgemqr.f"
	work[1] = (doublereal) lw;
#line 255 "sgemqr.f"
    }

#line 257 "sgemqr.f"
    if (*info != 0) {
#line 258 "sgemqr.f"
	i__1 = -(*info);
#line 258 "sgemqr.f"
	xerbla_("SGEMQR", &i__1, (ftnlen)6);
#line 259 "sgemqr.f"
	return 0;
#line 260 "sgemqr.f"
    } else if (lquery) {
#line 261 "sgemqr.f"
	return 0;
#line 262 "sgemqr.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 266 "sgemqr.f"
    i__1 = min(*m,*n);
#line 266 "sgemqr.f"
    if (min(i__1,*k) == 0) {
#line 267 "sgemqr.f"
	return 0;
#line 268 "sgemqr.f"
    }

/* Computing MAX */
#line 270 "sgemqr.f"
    i__1 = max(*m,*n);
#line 270 "sgemqr.f"
    if (left && *m <= *k || right && *n <= *k || mb <= *k || mb >= max(i__1,*
	    k)) {
#line 272 "sgemqr.f"
	sgemqrt_(side, trans, m, n, k, &nb, &a[a_offset], lda, &t[6], &nb, &
		c__[c_offset], ldc, &work[1], info, (ftnlen)1, (ftnlen)1);
#line 274 "sgemqr.f"
    } else {
#line 275 "sgemqr.f"
	slamtsqr_(side, trans, m, n, k, &mb, &nb, &a[a_offset], lda, &t[6], &
		nb, &c__[c_offset], ldc, &work[1], lwork, info, (ftnlen)1, (
		ftnlen)1);
#line 277 "sgemqr.f"
    }

#line 279 "sgemqr.f"
    work[1] = (doublereal) lw;

#line 281 "sgemqr.f"
    return 0;

/*     End of SGEMQR */

} /* sgemqr_ */

