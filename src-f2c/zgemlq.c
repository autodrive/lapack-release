#line 1 "zgemlq.f"
/* zgemlq.f -- translated by f2c (version 20100827).
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

#line 1 "zgemlq.f"

/*  Definition: */
/*  =========== */

/*      SUBROUTINE ZGEMLQ( SIDE, TRANS, M, N, K, A, LDA, T, */
/*     $                   TSIZE, C, LDC, WORK, LWORK, INFO ) */


/*     .. Scalar Arguments .. */
/*      CHARACTER          SIDE, TRANS */
/*      INTEGER            INFO, LDA, M, N, K, LDT, TSIZE, LWORK, LDC */
/*     .. */
/*     .. Array Arguments .. */
/*      COMPLEX*16         A( LDA, * ), T( * ), C(LDC, * ), WORK( * ) */
/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >     ZGEMLQ overwrites the general real M-by-N matrix C with */
/* > */
/* >                      SIDE = 'L'     SIDE = 'R' */
/* >      TRANS = 'N':      Q * C          C * Q */
/* >      TRANS = 'C':      Q**H * C       C * Q**H */
/* >      where Q is a complex unitary matrix defined as the product */
/* >      of blocked elementary reflectors computed by short wide */
/* >      LQ factorization (ZGELQ) */
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
/* > */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension */
/* >                               (LDA,M) if SIDE = 'L', */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          Part of the data structure to represent Q as returned by ZGELQ. */
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
/* >          T is COMPLEX*16 array, dimension (MAX(5,TSIZE)). */
/* >          Part of the data structure to represent Q as returned by ZGELQ. */
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
/* >          C is COMPLEX*16 array, dimension (LDC,N) */
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
/* >         (workspace) COMPLEX*16 array, dimension (MAX(1,LWORK)) */
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
/* >                           ZLASWLQ or ZGELQT */
/* > */
/* >  Depending on the matrix dimensions M and N, and row and column */
/* >  block sizes MB and NB returned by ILAENV, ZGELQ will use either */
/* >  ZLASWLQ (if the matrix is wide-and-short) or ZGELQT to compute */
/* >  the LQ factorization. */
/* >  This version of ZGEMLQ will use either ZLAMSWLQ or ZGEMLQT to */
/* >  multiply matrix Q by another matrix. */
/* >  Further Details in ZLAMSWLQ or ZGEMLQT. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgemlq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *t, integer 
	*tsize, doublecomplex *c__, integer *ldc, doublecomplex *work, 
	integer *lwork, integer *info, ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int zlamswlq_(char *, char *, integer *, integer *
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
    extern /* Subroutine */ int zgemlqt_(char *, char *, integer *, integer *,
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

#line 203 "zgemlq.f"
    /* Parameter adjustments */
#line 203 "zgemlq.f"
    a_dim1 = *lda;
#line 203 "zgemlq.f"
    a_offset = 1 + a_dim1;
#line 203 "zgemlq.f"
    a -= a_offset;
#line 203 "zgemlq.f"
    --t;
#line 203 "zgemlq.f"
    c_dim1 = *ldc;
#line 203 "zgemlq.f"
    c_offset = 1 + c_dim1;
#line 203 "zgemlq.f"
    c__ -= c_offset;
#line 203 "zgemlq.f"
    --work;
#line 203 "zgemlq.f"

#line 203 "zgemlq.f"
    /* Function Body */
#line 203 "zgemlq.f"
    lquery = *lwork == -1;
#line 204 "zgemlq.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 205 "zgemlq.f"
    tran = lsame_(trans, "C", (ftnlen)1, (ftnlen)1);
#line 206 "zgemlq.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 207 "zgemlq.f"
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);

#line 209 "zgemlq.f"
    mb = (integer) t[2].r;
#line 210 "zgemlq.f"
    nb = (integer) t[3].r;
#line 211 "zgemlq.f"
    if (left) {
#line 212 "zgemlq.f"
	lw = *n * mb;
#line 213 "zgemlq.f"
	mn = *m;
#line 214 "zgemlq.f"
    } else {
#line 215 "zgemlq.f"
	lw = *m * mb;
#line 216 "zgemlq.f"
	mn = *n;
#line 217 "zgemlq.f"
    }

#line 219 "zgemlq.f"
    if (nb > *k && mn > *k) {
#line 220 "zgemlq.f"
	if ((mn - *k) % (nb - *k) == 0) {
#line 221 "zgemlq.f"
	    nblcks = (mn - *k) / (nb - *k);
#line 222 "zgemlq.f"
	} else {
#line 223 "zgemlq.f"
	    nblcks = (mn - *k) / (nb - *k) + 1;
#line 224 "zgemlq.f"
	}
#line 225 "zgemlq.f"
    } else {
#line 226 "zgemlq.f"
	nblcks = 1;
#line 227 "zgemlq.f"
    }

#line 229 "zgemlq.f"
    *info = 0;
#line 230 "zgemlq.f"
    if (! left && ! right) {
#line 231 "zgemlq.f"
	*info = -1;
#line 232 "zgemlq.f"
    } else if (! tran && ! notran) {
#line 233 "zgemlq.f"
	*info = -2;
#line 234 "zgemlq.f"
    } else if (*m < 0) {
#line 235 "zgemlq.f"
	*info = -3;
#line 236 "zgemlq.f"
    } else if (*n < 0) {
#line 237 "zgemlq.f"
	*info = -4;
#line 238 "zgemlq.f"
    } else if (*k < 0 || *k > mn) {
#line 239 "zgemlq.f"
	*info = -5;
#line 240 "zgemlq.f"
    } else if (*lda < max(1,*k)) {
#line 241 "zgemlq.f"
	*info = -7;
#line 242 "zgemlq.f"
    } else if (*tsize < 5) {
#line 243 "zgemlq.f"
	*info = -9;
#line 244 "zgemlq.f"
    } else if (*ldc < max(1,*m)) {
#line 245 "zgemlq.f"
	*info = -11;
#line 246 "zgemlq.f"
    } else if (*lwork < max(1,lw) && ! lquery) {
#line 247 "zgemlq.f"
	*info = -13;
#line 248 "zgemlq.f"
    }

#line 250 "zgemlq.f"
    if (*info == 0) {
#line 251 "zgemlq.f"
	work[1].r = (doublereal) lw, work[1].i = 0.;
#line 252 "zgemlq.f"
    }

#line 254 "zgemlq.f"
    if (*info != 0) {
#line 255 "zgemlq.f"
	i__1 = -(*info);
#line 255 "zgemlq.f"
	xerbla_("ZGEMLQ", &i__1, (ftnlen)6);
#line 256 "zgemlq.f"
	return 0;
#line 257 "zgemlq.f"
    } else if (lquery) {
#line 258 "zgemlq.f"
	return 0;
#line 259 "zgemlq.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 263 "zgemlq.f"
    i__1 = min(*m,*n);
#line 263 "zgemlq.f"
    if (min(i__1,*k) == 0) {
#line 264 "zgemlq.f"
	return 0;
#line 265 "zgemlq.f"
    }

/* Computing MAX */
#line 267 "zgemlq.f"
    i__1 = max(*m,*n);
#line 267 "zgemlq.f"
    if (left && *m <= *k || right && *n <= *k || nb <= *k || nb >= max(i__1,*
	    k)) {
#line 269 "zgemlq.f"
	zgemlqt_(side, trans, m, n, k, &mb, &a[a_offset], lda, &t[6], &mb, &
		c__[c_offset], ldc, &work[1], info, (ftnlen)1, (ftnlen)1);
#line 271 "zgemlq.f"
    } else {
#line 272 "zgemlq.f"
	zlamswlq_(side, trans, m, n, k, &mb, &nb, &a[a_offset], lda, &t[6], &
		mb, &c__[c_offset], ldc, &work[1], lwork, info, (ftnlen)1, (
		ftnlen)1);
#line 274 "zgemlq.f"
    }

#line 276 "zgemlq.f"
    work[1].r = (doublereal) lw, work[1].i = 0.;

#line 278 "zgemlq.f"
    return 0;

/*     End of ZGEMLQ */

} /* zgemlq_ */

