#line 1 "dgeqp3.f"
/* dgeqp3.f -- translated by f2c (version 20100827).
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

#line 1 "dgeqp3.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b DGEQP3 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEQP3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqp3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqp3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqp3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            JPVT( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEQP3 computes a QR factorization with column pivoting of a */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the upper triangle of the array contains the */
/* >          min(M,N)-by-N upper trapezoidal matrix R; the elements below */
/* >          the diagonal, together with the array TAU, represent the */
/* >          orthogonal matrix Q as a product of min(M,N) elementary */
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
/* >          TAU is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO=0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. LWORK >= 3*N+1. */
/* >          For optimal performance LWORK >= 2*N+( N+1 )*NB, where NB */
/* >          is the optimal blocksize. */
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
/* >          = 0: successful exit. */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup doubleGEcomputational */

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
/* >     H(i) = I - tau * v * v**T */
/* > */
/* >  where tau is a real scalar, and v is a real/complex vector */
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
/* Subroutine */ int dgeqp3_(integer *m, integer *n, doublereal *a, integer *
	lda, integer *jpvt, doublereal *tau, doublereal *work, integer *lwork,
	 integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer j, jb, na, nb, sm, sn, nx, fjb, iws, nfxd;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static integer nbmin, minmn;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer minws;
    extern /* Subroutine */ int dlaqp2_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), dgeqrf_(integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlaqps_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *);
    static integer topbmn, sminmn;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer lwkopt;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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

#line 194 "dgeqp3.f"
    /* Parameter adjustments */
#line 194 "dgeqp3.f"
    a_dim1 = *lda;
#line 194 "dgeqp3.f"
    a_offset = 1 + a_dim1;
#line 194 "dgeqp3.f"
    a -= a_offset;
#line 194 "dgeqp3.f"
    --jpvt;
#line 194 "dgeqp3.f"
    --tau;
#line 194 "dgeqp3.f"
    --work;
#line 194 "dgeqp3.f"

#line 194 "dgeqp3.f"
    /* Function Body */
#line 194 "dgeqp3.f"
    *info = 0;
#line 195 "dgeqp3.f"
    lquery = *lwork == -1;
#line 196 "dgeqp3.f"
    if (*m < 0) {
#line 197 "dgeqp3.f"
	*info = -1;
#line 198 "dgeqp3.f"
    } else if (*n < 0) {
#line 199 "dgeqp3.f"
	*info = -2;
#line 200 "dgeqp3.f"
    } else if (*lda < max(1,*m)) {
#line 201 "dgeqp3.f"
	*info = -4;
#line 202 "dgeqp3.f"
    }

#line 204 "dgeqp3.f"
    if (*info == 0) {
#line 205 "dgeqp3.f"
	minmn = min(*m,*n);
#line 206 "dgeqp3.f"
	if (minmn == 0) {
#line 207 "dgeqp3.f"
	    iws = 1;
#line 208 "dgeqp3.f"
	    lwkopt = 1;
#line 209 "dgeqp3.f"
	} else {
#line 210 "dgeqp3.f"
	    iws = *n * 3 + 1;
#line 211 "dgeqp3.f"
	    nb = ilaenv_(&c__1, "DGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
#line 212 "dgeqp3.f"
	    lwkopt = (*n << 1) + (*n + 1) * nb;
#line 213 "dgeqp3.f"
	}
#line 214 "dgeqp3.f"
	work[1] = (doublereal) lwkopt;

#line 216 "dgeqp3.f"
	if (*lwork < iws && ! lquery) {
#line 217 "dgeqp3.f"
	    *info = -8;
#line 218 "dgeqp3.f"
	}
#line 219 "dgeqp3.f"
    }

#line 221 "dgeqp3.f"
    if (*info != 0) {
#line 222 "dgeqp3.f"
	i__1 = -(*info);
#line 222 "dgeqp3.f"
	xerbla_("DGEQP3", &i__1, (ftnlen)6);
#line 223 "dgeqp3.f"
	return 0;
#line 224 "dgeqp3.f"
    } else if (lquery) {
#line 225 "dgeqp3.f"
	return 0;
#line 226 "dgeqp3.f"
    }

/*     Quick return if possible. */

#line 230 "dgeqp3.f"
    if (minmn == 0) {
#line 231 "dgeqp3.f"
	return 0;
#line 232 "dgeqp3.f"
    }

/*     Move initial columns up front. */

#line 236 "dgeqp3.f"
    nfxd = 1;
#line 237 "dgeqp3.f"
    i__1 = *n;
#line 237 "dgeqp3.f"
    for (j = 1; j <= i__1; ++j) {
#line 238 "dgeqp3.f"
	if (jpvt[j] != 0) {
#line 239 "dgeqp3.f"
	    if (j != nfxd) {
#line 240 "dgeqp3.f"
		dswap_(m, &a[j * a_dim1 + 1], &c__1, &a[nfxd * a_dim1 + 1], &
			c__1);
#line 241 "dgeqp3.f"
		jpvt[j] = jpvt[nfxd];
#line 242 "dgeqp3.f"
		jpvt[nfxd] = j;
#line 243 "dgeqp3.f"
	    } else {
#line 244 "dgeqp3.f"
		jpvt[j] = j;
#line 245 "dgeqp3.f"
	    }
#line 246 "dgeqp3.f"
	    ++nfxd;
#line 247 "dgeqp3.f"
	} else {
#line 248 "dgeqp3.f"
	    jpvt[j] = j;
#line 249 "dgeqp3.f"
	}
#line 250 "dgeqp3.f"
/* L10: */
#line 250 "dgeqp3.f"
    }
#line 251 "dgeqp3.f"
    --nfxd;

/*     Factorize fixed columns */
/*  ======================= */

/*     Compute the QR factorization of fixed columns and update */
/*     remaining columns. */

#line 259 "dgeqp3.f"
    if (nfxd > 0) {
#line 260 "dgeqp3.f"
	na = min(*m,nfxd);
/* CC      CALL DGEQR2( M, NA, A, LDA, TAU, WORK, INFO ) */
#line 262 "dgeqp3.f"
	dgeqrf_(m, &na, &a[a_offset], lda, &tau[1], &work[1], lwork, info);
/* Computing MAX */
#line 263 "dgeqp3.f"
	i__1 = iws, i__2 = (integer) work[1];
#line 263 "dgeqp3.f"
	iws = max(i__1,i__2);
#line 264 "dgeqp3.f"
	if (na < *n) {
/* CC         CALL DORM2R( 'Left', 'Transpose', M, N-NA, NA, A, LDA, */
/* CC  $                   TAU, A( 1, NA+1 ), LDA, WORK, INFO ) */
#line 267 "dgeqp3.f"
	    i__1 = *n - na;
#line 267 "dgeqp3.f"
	    dormqr_("Left", "Transpose", m, &i__1, &na, &a[a_offset], lda, &
		    tau[1], &a[(na + 1) * a_dim1 + 1], lda, &work[1], lwork, 
		    info, (ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 269 "dgeqp3.f"
	    i__1 = iws, i__2 = (integer) work[1];
#line 269 "dgeqp3.f"
	    iws = max(i__1,i__2);
#line 270 "dgeqp3.f"
	}
#line 271 "dgeqp3.f"
    }

/*     Factorize free columns */
/*  ====================== */

#line 276 "dgeqp3.f"
    if (nfxd < minmn) {

#line 278 "dgeqp3.f"
	sm = *m - nfxd;
#line 279 "dgeqp3.f"
	sn = *n - nfxd;
#line 280 "dgeqp3.f"
	sminmn = minmn - nfxd;

/*        Determine the block size. */

#line 284 "dgeqp3.f"
	nb = ilaenv_(&c__1, "DGEQRF", " ", &sm, &sn, &c_n1, &c_n1, (ftnlen)6, 
		(ftnlen)1);
#line 285 "dgeqp3.f"
	nbmin = 2;
#line 286 "dgeqp3.f"
	nx = 0;

#line 288 "dgeqp3.f"
	if (nb > 1 && nb < sminmn) {

/*           Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 292 "dgeqp3.f"
	    i__1 = 0, i__2 = ilaenv_(&c__3, "DGEQRF", " ", &sm, &sn, &c_n1, &
		    c_n1, (ftnlen)6, (ftnlen)1);
#line 292 "dgeqp3.f"
	    nx = max(i__1,i__2);


#line 296 "dgeqp3.f"
	    if (nx < sminmn) {

/*              Determine if workspace is large enough for blocked code. */

#line 300 "dgeqp3.f"
		minws = (sn << 1) + (sn + 1) * nb;
#line 301 "dgeqp3.f"
		iws = max(iws,minws);
#line 302 "dgeqp3.f"
		if (*lwork < minws) {

/*                 Not enough workspace to use optimal NB: Reduce NB and */
/*                 determine the minimum value of NB. */

#line 307 "dgeqp3.f"
		    nb = (*lwork - (sn << 1)) / (sn + 1);
/* Computing MAX */
#line 308 "dgeqp3.f"
		    i__1 = 2, i__2 = ilaenv_(&c__2, "DGEQRF", " ", &sm, &sn, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 308 "dgeqp3.f"
		    nbmin = max(i__1,i__2);


#line 312 "dgeqp3.f"
		}
#line 313 "dgeqp3.f"
	    }
#line 314 "dgeqp3.f"
	}

/*        Initialize partial column norms. The first N elements of work */
/*        store the exact column norms. */

#line 319 "dgeqp3.f"
	i__1 = *n;
#line 319 "dgeqp3.f"
	for (j = nfxd + 1; j <= i__1; ++j) {
#line 320 "dgeqp3.f"
	    work[j] = dnrm2_(&sm, &a[nfxd + 1 + j * a_dim1], &c__1);
#line 321 "dgeqp3.f"
	    work[*n + j] = work[j];
#line 322 "dgeqp3.f"
/* L20: */
#line 322 "dgeqp3.f"
	}

#line 324 "dgeqp3.f"
	if (nb >= nbmin && nb < sminmn && nx < sminmn) {

/*           Use blocked code initially. */

#line 329 "dgeqp3.f"
	    j = nfxd + 1;

/*           Compute factorization: while loop. */


#line 334 "dgeqp3.f"
	    topbmn = minmn - nx;
#line 335 "dgeqp3.f"
L30:
#line 336 "dgeqp3.f"
	    if (j <= topbmn) {
/* Computing MIN */
#line 337 "dgeqp3.f"
		i__1 = nb, i__2 = topbmn - j + 1;
#line 337 "dgeqp3.f"
		jb = min(i__1,i__2);

/*              Factorize JB columns among columns J:N. */

#line 341 "dgeqp3.f"
		i__1 = *n - j + 1;
#line 341 "dgeqp3.f"
		i__2 = j - 1;
#line 341 "dgeqp3.f"
		i__3 = *n - j + 1;
#line 341 "dgeqp3.f"
		dlaqps_(m, &i__1, &i__2, &jb, &fjb, &a[j * a_dim1 + 1], lda, &
			jpvt[j], &tau[j], &work[j], &work[*n + j], &work[(*n 
			<< 1) + 1], &work[(*n << 1) + jb + 1], &i__3);

#line 345 "dgeqp3.f"
		j += fjb;
#line 346 "dgeqp3.f"
		goto L30;
#line 347 "dgeqp3.f"
	    }
#line 348 "dgeqp3.f"
	} else {
#line 349 "dgeqp3.f"
	    j = nfxd + 1;
#line 350 "dgeqp3.f"
	}

/*        Use unblocked code to factor the last or only block. */


#line 355 "dgeqp3.f"
	if (j <= minmn) {
#line 355 "dgeqp3.f"
	    i__1 = *n - j + 1;
#line 355 "dgeqp3.f"
	    i__2 = j - 1;
#line 355 "dgeqp3.f"
	    dlaqp2_(m, &i__1, &i__2, &a[j * a_dim1 + 1], lda, &jpvt[j], &tau[
		    j], &work[j], &work[*n + j], &work[(*n << 1) + 1]);
#line 355 "dgeqp3.f"
	}

#line 360 "dgeqp3.f"
    }

#line 362 "dgeqp3.f"
    work[1] = (doublereal) iws;
#line 363 "dgeqp3.f"
    return 0;

/*     End of DGEQP3 */

} /* dgeqp3_ */

