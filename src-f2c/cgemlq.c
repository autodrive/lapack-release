#line 1 "cgemlq.f"
/* cgemlq.f -- translated by f2c (version 20100827).
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

#line 1 "cgemlq.f"

/*  Definition: */
/*  =========== */

/*      SUBROUTINE CGEMLQ( SIDE, TRANS, M, N, K, A, LDA, T, */
/*     $                   TSIZE, C, LDC, WORK, LWORK, INFO ) */


/*     .. Scalar Arguments .. */
/*     CHARACTER         SIDE, TRANS */
/*     INTEGER           INFO, LDA, M, N, K, LDT, TSIZE, LWORK, LDC */
/*     .. */
/*     .. Array Arguments .. */
/*     COMPLEX           A( LDA, * ), T( * ), C(LDC, * ), WORK( * ) */
/*     .. */

/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >     CGEMLQ overwrites the general real M-by-N matrix C with */
/* > */
/* >                      SIDE = 'L'     SIDE = 'R' */
/* >      TRANS = 'N':      Q * C          C * Q */
/* >      TRANS = 'C':      Q**H * C       C * Q**H */
/* >      where Q is a complex unitary matrix defined as the product */
/* >      of blocked elementary reflectors computed by short wide */
/* >      LQ factorization (CGELQ) */
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
/* >          A is COMPLEX array, dimension */
/* >                               (LDA,M) if SIDE = 'L', */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          Part of the data structure to represent Q as returned by CGELQ. */
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
/* >          T is COMPLEX array, dimension (MAX(5,TSIZE)). */
/* >          Part of the data structure to represent Q as returned by CGELQ. */
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
/* >          C is COMPLEX array, dimension (LDC,N) */
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
/* >         (workspace) COMPLEX array, dimension (MAX(1,LWORK)) */
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
/* >                           CLASWQR or CGELQT */
/* > */
/* >  Depending on the matrix dimensions M and N, and row and column */
/* >  block sizes MB and NB returned by ILAENV, CGELQ will use either */
/* >  CLASWLQ (if the matrix is wide-and-short) or CGELQT to compute */
/* >  the LQ factorization. */
/* >  This version of CGEMLQ will use either CLAMSWLQ or CGEMLQT to */
/* >  multiply matrix Q by another matrix. */
/* >  Further Details in CLAMSWLQ or CGEMLQT. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cgemlq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *t, integer 
	*tsize, doublecomplex *c__, integer *ldc, doublecomplex *work, 
	integer *lwork, integer *info, ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1;
    doublereal d__1;

    /* Local variables */
    extern /* Subroutine */ int clamswlq_(char *, char *, integer *, integer *
	    , integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen, ftnlen);
    static integer mb, nb, mn, lw;
    static logical left, tran;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    static integer nblcks;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical notran, lquery;
    extern /* Subroutine */ int cgemlqt_(char *, char *, integer *, integer *,
	     integer *, integer *, doublecomplex *, integer *, doublecomplex *
	    , integer *, doublecomplex *, integer *, doublecomplex *, integer 
	    *, ftnlen, ftnlen);


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

#line 204 "cgemlq.f"
    /* Parameter adjustments */
#line 204 "cgemlq.f"
    a_dim1 = *lda;
#line 204 "cgemlq.f"
    a_offset = 1 + a_dim1;
#line 204 "cgemlq.f"
    a -= a_offset;
#line 204 "cgemlq.f"
    --t;
#line 204 "cgemlq.f"
    c_dim1 = *ldc;
#line 204 "cgemlq.f"
    c_offset = 1 + c_dim1;
#line 204 "cgemlq.f"
    c__ -= c_offset;
#line 204 "cgemlq.f"
    --work;
#line 204 "cgemlq.f"

#line 204 "cgemlq.f"
    /* Function Body */
#line 204 "cgemlq.f"
    lquery = *lwork == -1;
#line 205 "cgemlq.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 206 "cgemlq.f"
    tran = lsame_(trans, "C", (ftnlen)1, (ftnlen)1);
#line 207 "cgemlq.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 208 "cgemlq.f"
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);

#line 210 "cgemlq.f"
    mb = (integer) t[2].r;
#line 211 "cgemlq.f"
    nb = (integer) t[3].r;
#line 212 "cgemlq.f"
    if (left) {
#line 213 "cgemlq.f"
	lw = *n * mb;
#line 214 "cgemlq.f"
	mn = *m;
#line 215 "cgemlq.f"
    } else {
#line 216 "cgemlq.f"
	lw = *m * mb;
#line 217 "cgemlq.f"
	mn = *n;
#line 218 "cgemlq.f"
    }

#line 220 "cgemlq.f"
    if (nb > *k && mn > *k) {
#line 221 "cgemlq.f"
	if ((mn - *k) % (nb - *k) == 0) {
#line 222 "cgemlq.f"
	    nblcks = (mn - *k) / (nb - *k);
#line 223 "cgemlq.f"
	} else {
#line 224 "cgemlq.f"
	    nblcks = (mn - *k) / (nb - *k) + 1;
#line 225 "cgemlq.f"
	}
#line 226 "cgemlq.f"
    } else {
#line 227 "cgemlq.f"
	nblcks = 1;
#line 228 "cgemlq.f"
    }

#line 230 "cgemlq.f"
    *info = 0;
#line 231 "cgemlq.f"
    if (! left && ! right) {
#line 232 "cgemlq.f"
	*info = -1;
#line 233 "cgemlq.f"
    } else if (! tran && ! notran) {
#line 234 "cgemlq.f"
	*info = -2;
#line 235 "cgemlq.f"
    } else if (*m < 0) {
#line 236 "cgemlq.f"
	*info = -3;
#line 237 "cgemlq.f"
    } else if (*n < 0) {
#line 238 "cgemlq.f"
	*info = -4;
#line 239 "cgemlq.f"
    } else if (*k < 0 || *k > mn) {
#line 240 "cgemlq.f"
	*info = -5;
#line 241 "cgemlq.f"
    } else if (*lda < max(1,*k)) {
#line 242 "cgemlq.f"
	*info = -7;
#line 243 "cgemlq.f"
    } else if (*tsize < 5) {
#line 244 "cgemlq.f"
	*info = -9;
#line 245 "cgemlq.f"
    } else if (*ldc < max(1,*m)) {
#line 246 "cgemlq.f"
	*info = -11;
#line 247 "cgemlq.f"
    } else if (*lwork < max(1,lw) && ! lquery) {
#line 248 "cgemlq.f"
	*info = -13;
#line 249 "cgemlq.f"
    }

#line 251 "cgemlq.f"
    if (*info == 0) {
#line 252 "cgemlq.f"
	d__1 = (doublereal) lw;
#line 252 "cgemlq.f"
	work[1].r = d__1, work[1].i = 0.;
#line 253 "cgemlq.f"
    }

#line 255 "cgemlq.f"
    if (*info != 0) {
#line 256 "cgemlq.f"
	i__1 = -(*info);
#line 256 "cgemlq.f"
	xerbla_("CGEMLQ", &i__1, (ftnlen)6);
#line 257 "cgemlq.f"
	return 0;
#line 258 "cgemlq.f"
    } else if (lquery) {
#line 259 "cgemlq.f"
	return 0;
#line 260 "cgemlq.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 264 "cgemlq.f"
    i__1 = min(*m,*n);
#line 264 "cgemlq.f"
    if (min(i__1,*k) == 0) {
#line 265 "cgemlq.f"
	return 0;
#line 266 "cgemlq.f"
    }

/* Computing MAX */
#line 268 "cgemlq.f"
    i__1 = max(*m,*n);
#line 268 "cgemlq.f"
    if (left && *m <= *k || right && *n <= *k || nb <= *k || nb >= max(i__1,*
	    k)) {
#line 270 "cgemlq.f"
	cgemlqt_(side, trans, m, n, k, &mb, &a[a_offset], lda, &t[6], &mb, &
		c__[c_offset], ldc, &work[1], info, (ftnlen)1, (ftnlen)1);
#line 272 "cgemlq.f"
    } else {
#line 273 "cgemlq.f"
	clamswlq_(side, trans, m, n, k, &mb, &nb, &a[a_offset], lda, &t[6], &
		mb, &c__[c_offset], ldc, &work[1], lwork, info, (ftnlen)1, (
		ftnlen)1);
#line 275 "cgemlq.f"
    }

#line 277 "cgemlq.f"
    d__1 = (doublereal) lw;
#line 277 "cgemlq.f"
    work[1].r = d__1, work[1].i = 0.;

#line 279 "cgemlq.f"
    return 0;

/*     End of CGEMLQ */

} /* cgemlq_ */

