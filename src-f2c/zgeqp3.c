#line 1 "zgeqp3.f"
/* zgeqp3.f -- translated by f2c (version 20100827).
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

#line 1 "zgeqp3.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b ZGEQP3 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEQP3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeqp3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeqp3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeqp3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, RWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            JPVT( * ) */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEQP3 computes a QR factorization with column pivoting of a */
/* > matrix A:  A*P = Q*R  using Level 3 BLAS. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the upper triangle of the array contains the */
/* >          min(M,N)-by-N upper trapezoidal matrix R; the elements below */
/* >          the diagonal, together with the array TAU, represent the */
/* >          unitary matrix Q as a product of min(M,N) elementary */
/* >          reflectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] JPVT */
/* > \verbatim */
/* >          JPVT is INTEGER array, dimension (N) */
/* >          On entry, if JPVT(J).ne.0, the J-th column of A is permuted */
/* >          to the front of A*P (a leading column); if JPVT(J)=0, */
/* >          the J-th column of A is a free column. */
/* >          On exit, if JPVT(J)=K, then the J-th column of A*P was the */
/* >          the K-th column of A. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX*16 array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO=0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. LWORK >= N+1. */
/* >          For optimal performance LWORK >= ( N+1 )*NB, where NB */
/* >          is the optimal blocksize. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit. */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

/* > \ingroup complex16GEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix Q is represented as a product of elementary reflectors */
/* > */
/* >     Q = H(1) H(2) . . . H(k), where k = min(m,n). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a real/complex vector */
/* >  with v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in */
/* >  A(i+1:m,i), and tau in TAU(i). */
/* > \endverbatim */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain */
/* >    X. Sun, Computer Science Dept., Duke University, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgeqp3_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, integer *jpvt, doublecomplex *tau, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Local variables */
    static integer j, jb, na, nb, sm, sn, nx, fjb, iws, nfxd, nbmin, minmn, 
	    minws;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zlaqp2_(integer *, integer *, 
	    integer *, doublecomplex *, integer *, integer *, doublecomplex *,
	     doublereal *, doublereal *, doublecomplex *);
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zgeqrf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    );
    static integer topbmn, sminmn;
    extern /* Subroutine */ int zlaqps_(integer *, integer *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, integer *, 
	    doublecomplex *, doublereal *, doublereal *, doublecomplex *, 
	    doublecomplex *, integer *);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int zunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test input arguments */
/*  ==================== */

#line 203 "zgeqp3.f"
    /* Parameter adjustments */
#line 203 "zgeqp3.f"
    a_dim1 = *lda;
#line 203 "zgeqp3.f"
    a_offset = 1 + a_dim1;
#line 203 "zgeqp3.f"
    a -= a_offset;
#line 203 "zgeqp3.f"
    --jpvt;
#line 203 "zgeqp3.f"
    --tau;
#line 203 "zgeqp3.f"
    --work;
#line 203 "zgeqp3.f"
    --rwork;
#line 203 "zgeqp3.f"

#line 203 "zgeqp3.f"
    /* Function Body */
#line 203 "zgeqp3.f"
    *info = 0;
#line 204 "zgeqp3.f"
    lquery = *lwork == -1;
#line 205 "zgeqp3.f"
    if (*m < 0) {
#line 206 "zgeqp3.f"
	*info = -1;
#line 207 "zgeqp3.f"
    } else if (*n < 0) {
#line 208 "zgeqp3.f"
	*info = -2;
#line 209 "zgeqp3.f"
    } else if (*lda < max(1,*m)) {
#line 210 "zgeqp3.f"
	*info = -4;
#line 211 "zgeqp3.f"
    }

#line 213 "zgeqp3.f"
    if (*info == 0) {
#line 214 "zgeqp3.f"
	minmn = min(*m,*n);
#line 215 "zgeqp3.f"
	if (minmn == 0) {
#line 216 "zgeqp3.f"
	    iws = 1;
#line 217 "zgeqp3.f"
	    lwkopt = 1;
#line 218 "zgeqp3.f"
	} else {
#line 219 "zgeqp3.f"
	    iws = *n + 1;
#line 220 "zgeqp3.f"
	    nb = ilaenv_(&c__1, "ZGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
#line 221 "zgeqp3.f"
	    lwkopt = (*n + 1) * nb;
#line 222 "zgeqp3.f"
	}
#line 223 "zgeqp3.f"
	z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 223 "zgeqp3.f"
	work[1].r = z__1.r, work[1].i = z__1.i;

#line 225 "zgeqp3.f"
	if (*lwork < iws && ! lquery) {
#line 226 "zgeqp3.f"
	    *info = -8;
#line 227 "zgeqp3.f"
	}
#line 228 "zgeqp3.f"
    }

#line 230 "zgeqp3.f"
    if (*info != 0) {
#line 231 "zgeqp3.f"
	i__1 = -(*info);
#line 231 "zgeqp3.f"
	xerbla_("ZGEQP3", &i__1, (ftnlen)6);
#line 232 "zgeqp3.f"
	return 0;
#line 233 "zgeqp3.f"
    } else if (lquery) {
#line 234 "zgeqp3.f"
	return 0;
#line 235 "zgeqp3.f"
    }

/*     Move initial columns up front. */

#line 239 "zgeqp3.f"
    nfxd = 1;
#line 240 "zgeqp3.f"
    i__1 = *n;
#line 240 "zgeqp3.f"
    for (j = 1; j <= i__1; ++j) {
#line 241 "zgeqp3.f"
	if (jpvt[j] != 0) {
#line 242 "zgeqp3.f"
	    if (j != nfxd) {
#line 243 "zgeqp3.f"
		zswap_(m, &a[j * a_dim1 + 1], &c__1, &a[nfxd * a_dim1 + 1], &
			c__1);
#line 244 "zgeqp3.f"
		jpvt[j] = jpvt[nfxd];
#line 245 "zgeqp3.f"
		jpvt[nfxd] = j;
#line 246 "zgeqp3.f"
	    } else {
#line 247 "zgeqp3.f"
		jpvt[j] = j;
#line 248 "zgeqp3.f"
	    }
#line 249 "zgeqp3.f"
	    ++nfxd;
#line 250 "zgeqp3.f"
	} else {
#line 251 "zgeqp3.f"
	    jpvt[j] = j;
#line 252 "zgeqp3.f"
	}
#line 253 "zgeqp3.f"
/* L10: */
#line 253 "zgeqp3.f"
    }
#line 254 "zgeqp3.f"
    --nfxd;

/*     Factorize fixed columns */
/*  ======================= */

/*     Compute the QR factorization of fixed columns and update */
/*     remaining columns. */

#line 262 "zgeqp3.f"
    if (nfxd > 0) {
#line 263 "zgeqp3.f"
	na = min(*m,nfxd);
/* CC      CALL ZGEQR2( M, NA, A, LDA, TAU, WORK, INFO ) */
#line 265 "zgeqp3.f"
	zgeqrf_(m, &na, &a[a_offset], lda, &tau[1], &work[1], lwork, info);
/* Computing MAX */
#line 266 "zgeqp3.f"
	i__1 = iws, i__2 = (integer) work[1].r;
#line 266 "zgeqp3.f"
	iws = max(i__1,i__2);
#line 267 "zgeqp3.f"
	if (na < *n) {
/* CC         CALL ZUNM2R( 'Left', 'Conjugate Transpose', M, N-NA, */
/* CC  $                   NA, A, LDA, TAU, A( 1, NA+1 ), LDA, WORK, */
/* CC  $                   INFO ) */
#line 271 "zgeqp3.f"
	    i__1 = *n - na;
#line 271 "zgeqp3.f"
	    zunmqr_("Left", "Conjugate Transpose", m, &i__1, &na, &a[a_offset]
		    , lda, &tau[1], &a[(na + 1) * a_dim1 + 1], lda, &work[1], 
		    lwork, info, (ftnlen)4, (ftnlen)19);
/* Computing MAX */
#line 274 "zgeqp3.f"
	    i__1 = iws, i__2 = (integer) work[1].r;
#line 274 "zgeqp3.f"
	    iws = max(i__1,i__2);
#line 275 "zgeqp3.f"
	}
#line 276 "zgeqp3.f"
    }

/*     Factorize free columns */
/*  ====================== */

#line 281 "zgeqp3.f"
    if (nfxd < minmn) {

#line 283 "zgeqp3.f"
	sm = *m - nfxd;
#line 284 "zgeqp3.f"
	sn = *n - nfxd;
#line 285 "zgeqp3.f"
	sminmn = minmn - nfxd;

/*        Determine the block size. */

#line 289 "zgeqp3.f"
	nb = ilaenv_(&c__1, "ZGEQRF", " ", &sm, &sn, &c_n1, &c_n1, (ftnlen)6, 
		(ftnlen)1);
#line 290 "zgeqp3.f"
	nbmin = 2;
#line 291 "zgeqp3.f"
	nx = 0;

#line 293 "zgeqp3.f"
	if (nb > 1 && nb < sminmn) {

/*           Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 297 "zgeqp3.f"
	    i__1 = 0, i__2 = ilaenv_(&c__3, "ZGEQRF", " ", &sm, &sn, &c_n1, &
		    c_n1, (ftnlen)6, (ftnlen)1);
#line 297 "zgeqp3.f"
	    nx = max(i__1,i__2);


#line 301 "zgeqp3.f"
	    if (nx < sminmn) {

/*              Determine if workspace is large enough for blocked code. */

#line 305 "zgeqp3.f"
		minws = (sn + 1) * nb;
#line 306 "zgeqp3.f"
		iws = max(iws,minws);
#line 307 "zgeqp3.f"
		if (*lwork < minws) {

/*                 Not enough workspace to use optimal NB: Reduce NB and */
/*                 determine the minimum value of NB. */

#line 312 "zgeqp3.f"
		    nb = *lwork / (sn + 1);
/* Computing MAX */
#line 313 "zgeqp3.f"
		    i__1 = 2, i__2 = ilaenv_(&c__2, "ZGEQRF", " ", &sm, &sn, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 313 "zgeqp3.f"
		    nbmin = max(i__1,i__2);


#line 317 "zgeqp3.f"
		}
#line 318 "zgeqp3.f"
	    }
#line 319 "zgeqp3.f"
	}

/*        Initialize partial column norms. The first N elements of work */
/*        store the exact column norms. */

#line 324 "zgeqp3.f"
	i__1 = *n;
#line 324 "zgeqp3.f"
	for (j = nfxd + 1; j <= i__1; ++j) {
#line 325 "zgeqp3.f"
	    rwork[j] = dznrm2_(&sm, &a[nfxd + 1 + j * a_dim1], &c__1);
#line 326 "zgeqp3.f"
	    rwork[*n + j] = rwork[j];
#line 327 "zgeqp3.f"
/* L20: */
#line 327 "zgeqp3.f"
	}

#line 329 "zgeqp3.f"
	if (nb >= nbmin && nb < sminmn && nx < sminmn) {

/*           Use blocked code initially. */

#line 334 "zgeqp3.f"
	    j = nfxd + 1;

/*           Compute factorization: while loop. */


#line 339 "zgeqp3.f"
	    topbmn = minmn - nx;
#line 340 "zgeqp3.f"
L30:
#line 341 "zgeqp3.f"
	    if (j <= topbmn) {
/* Computing MIN */
#line 342 "zgeqp3.f"
		i__1 = nb, i__2 = topbmn - j + 1;
#line 342 "zgeqp3.f"
		jb = min(i__1,i__2);

/*              Factorize JB columns among columns J:N. */

#line 346 "zgeqp3.f"
		i__1 = *n - j + 1;
#line 346 "zgeqp3.f"
		i__2 = j - 1;
#line 346 "zgeqp3.f"
		i__3 = *n - j + 1;
#line 346 "zgeqp3.f"
		zlaqps_(m, &i__1, &i__2, &jb, &fjb, &a[j * a_dim1 + 1], lda, &
			jpvt[j], &tau[j], &rwork[j], &rwork[*n + j], &work[1],
			 &work[jb + 1], &i__3);

#line 351 "zgeqp3.f"
		j += fjb;
#line 352 "zgeqp3.f"
		goto L30;
#line 353 "zgeqp3.f"
	    }
#line 354 "zgeqp3.f"
	} else {
#line 355 "zgeqp3.f"
	    j = nfxd + 1;
#line 356 "zgeqp3.f"
	}

/*        Use unblocked code to factor the last or only block. */


#line 361 "zgeqp3.f"
	if (j <= minmn) {
#line 361 "zgeqp3.f"
	    i__1 = *n - j + 1;
#line 361 "zgeqp3.f"
	    i__2 = j - 1;
#line 361 "zgeqp3.f"
	    zlaqp2_(m, &i__1, &i__2, &a[j * a_dim1 + 1], lda, &jpvt[j], &tau[
		    j], &rwork[j], &rwork[*n + j], &work[1]);
#line 361 "zgeqp3.f"
	}

#line 365 "zgeqp3.f"
    }

#line 367 "zgeqp3.f"
    z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 367 "zgeqp3.f"
    work[1].r = z__1.r, work[1].i = z__1.i;
#line 368 "zgeqp3.f"
    return 0;

/*     End of ZGEQP3 */

} /* zgeqp3_ */

