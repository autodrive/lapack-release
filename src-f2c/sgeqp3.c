#line 1 "sgeqp3.f"
/* sgeqp3.f -- translated by f2c (version 20100827).
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

#line 1 "sgeqp3.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b SGEQP3 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGEQP3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeqp3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeqp3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeqp3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            JPVT( * ) */
/*       REAL               A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEQP3 computes a QR factorization with column pivoting of a */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          TAU is REAL array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
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

/* > \ingroup realGEcomputational */

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
/* Subroutine */ int sgeqp3_(integer *m, integer *n, doublereal *a, integer *
	lda, integer *jpvt, doublereal *tau, doublereal *work, integer *lwork,
	 integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer j, jb, na, nb, sm, sn, nx, fjb, iws, nfxd;
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    static integer nbmin, minmn, minws;
    extern /* Subroutine */ int sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), slaqp2_(integer *, integer *, integer *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int sgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
    static integer topbmn, sminmn;
    extern /* Subroutine */ int slaqps_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int sormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);


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

#line 191 "sgeqp3.f"
    /* Parameter adjustments */
#line 191 "sgeqp3.f"
    a_dim1 = *lda;
#line 191 "sgeqp3.f"
    a_offset = 1 + a_dim1;
#line 191 "sgeqp3.f"
    a -= a_offset;
#line 191 "sgeqp3.f"
    --jpvt;
#line 191 "sgeqp3.f"
    --tau;
#line 191 "sgeqp3.f"
    --work;
#line 191 "sgeqp3.f"

#line 191 "sgeqp3.f"
    /* Function Body */
#line 191 "sgeqp3.f"
    *info = 0;
#line 192 "sgeqp3.f"
    lquery = *lwork == -1;
#line 193 "sgeqp3.f"
    if (*m < 0) {
#line 194 "sgeqp3.f"
	*info = -1;
#line 195 "sgeqp3.f"
    } else if (*n < 0) {
#line 196 "sgeqp3.f"
	*info = -2;
#line 197 "sgeqp3.f"
    } else if (*lda < max(1,*m)) {
#line 198 "sgeqp3.f"
	*info = -4;
#line 199 "sgeqp3.f"
    }

#line 201 "sgeqp3.f"
    if (*info == 0) {
#line 202 "sgeqp3.f"
	minmn = min(*m,*n);
#line 203 "sgeqp3.f"
	if (minmn == 0) {
#line 204 "sgeqp3.f"
	    iws = 1;
#line 205 "sgeqp3.f"
	    lwkopt = 1;
#line 206 "sgeqp3.f"
	} else {
#line 207 "sgeqp3.f"
	    iws = *n * 3 + 1;
#line 208 "sgeqp3.f"
	    nb = ilaenv_(&c__1, "SGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
#line 209 "sgeqp3.f"
	    lwkopt = (*n << 1) + (*n + 1) * nb;
#line 210 "sgeqp3.f"
	}
#line 211 "sgeqp3.f"
	work[1] = (doublereal) lwkopt;

#line 213 "sgeqp3.f"
	if (*lwork < iws && ! lquery) {
#line 214 "sgeqp3.f"
	    *info = -8;
#line 215 "sgeqp3.f"
	}
#line 216 "sgeqp3.f"
    }

#line 218 "sgeqp3.f"
    if (*info != 0) {
#line 219 "sgeqp3.f"
	i__1 = -(*info);
#line 219 "sgeqp3.f"
	xerbla_("SGEQP3", &i__1, (ftnlen)6);
#line 220 "sgeqp3.f"
	return 0;
#line 221 "sgeqp3.f"
    } else if (lquery) {
#line 222 "sgeqp3.f"
	return 0;
#line 223 "sgeqp3.f"
    }

/*     Quick return if possible. */

#line 227 "sgeqp3.f"
    if (minmn == 0) {
#line 228 "sgeqp3.f"
	return 0;
#line 229 "sgeqp3.f"
    }

/*     Move initial columns up front. */

#line 233 "sgeqp3.f"
    nfxd = 1;
#line 234 "sgeqp3.f"
    i__1 = *n;
#line 234 "sgeqp3.f"
    for (j = 1; j <= i__1; ++j) {
#line 235 "sgeqp3.f"
	if (jpvt[j] != 0) {
#line 236 "sgeqp3.f"
	    if (j != nfxd) {
#line 237 "sgeqp3.f"
		sswap_(m, &a[j * a_dim1 + 1], &c__1, &a[nfxd * a_dim1 + 1], &
			c__1);
#line 238 "sgeqp3.f"
		jpvt[j] = jpvt[nfxd];
#line 239 "sgeqp3.f"
		jpvt[nfxd] = j;
#line 240 "sgeqp3.f"
	    } else {
#line 241 "sgeqp3.f"
		jpvt[j] = j;
#line 242 "sgeqp3.f"
	    }
#line 243 "sgeqp3.f"
	    ++nfxd;
#line 244 "sgeqp3.f"
	} else {
#line 245 "sgeqp3.f"
	    jpvt[j] = j;
#line 246 "sgeqp3.f"
	}
#line 247 "sgeqp3.f"
/* L10: */
#line 247 "sgeqp3.f"
    }
#line 248 "sgeqp3.f"
    --nfxd;

/*     Factorize fixed columns */
/*  ======================= */

/*     Compute the QR factorization of fixed columns and update */
/*     remaining columns. */

#line 256 "sgeqp3.f"
    if (nfxd > 0) {
#line 257 "sgeqp3.f"
	na = min(*m,nfxd);
/* CC      CALL SGEQR2( M, NA, A, LDA, TAU, WORK, INFO ) */
#line 259 "sgeqp3.f"
	sgeqrf_(m, &na, &a[a_offset], lda, &tau[1], &work[1], lwork, info);
/* Computing MAX */
#line 260 "sgeqp3.f"
	i__1 = iws, i__2 = (integer) work[1];
#line 260 "sgeqp3.f"
	iws = max(i__1,i__2);
#line 261 "sgeqp3.f"
	if (na < *n) {
/* CC         CALL SORM2R( 'Left', 'Transpose', M, N-NA, NA, A, LDA, */
/* CC  $                   TAU, A( 1, NA+1 ), LDA, WORK, INFO ) */
#line 264 "sgeqp3.f"
	    i__1 = *n - na;
#line 264 "sgeqp3.f"
	    sormqr_("Left", "Transpose", m, &i__1, &na, &a[a_offset], lda, &
		    tau[1], &a[(na + 1) * a_dim1 + 1], lda, &work[1], lwork, 
		    info, (ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 266 "sgeqp3.f"
	    i__1 = iws, i__2 = (integer) work[1];
#line 266 "sgeqp3.f"
	    iws = max(i__1,i__2);
#line 267 "sgeqp3.f"
	}
#line 268 "sgeqp3.f"
    }

/*     Factorize free columns */
/*  ====================== */

#line 273 "sgeqp3.f"
    if (nfxd < minmn) {

#line 275 "sgeqp3.f"
	sm = *m - nfxd;
#line 276 "sgeqp3.f"
	sn = *n - nfxd;
#line 277 "sgeqp3.f"
	sminmn = minmn - nfxd;

/*        Determine the block size. */

#line 281 "sgeqp3.f"
	nb = ilaenv_(&c__1, "SGEQRF", " ", &sm, &sn, &c_n1, &c_n1, (ftnlen)6, 
		(ftnlen)1);
#line 282 "sgeqp3.f"
	nbmin = 2;
#line 283 "sgeqp3.f"
	nx = 0;

#line 285 "sgeqp3.f"
	if (nb > 1 && nb < sminmn) {

/*           Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 289 "sgeqp3.f"
	    i__1 = 0, i__2 = ilaenv_(&c__3, "SGEQRF", " ", &sm, &sn, &c_n1, &
		    c_n1, (ftnlen)6, (ftnlen)1);
#line 289 "sgeqp3.f"
	    nx = max(i__1,i__2);


#line 293 "sgeqp3.f"
	    if (nx < sminmn) {

/*              Determine if workspace is large enough for blocked code. */

#line 297 "sgeqp3.f"
		minws = (sn << 1) + (sn + 1) * nb;
#line 298 "sgeqp3.f"
		iws = max(iws,minws);
#line 299 "sgeqp3.f"
		if (*lwork < minws) {

/*                 Not enough workspace to use optimal NB: Reduce NB and */
/*                 determine the minimum value of NB. */

#line 304 "sgeqp3.f"
		    nb = (*lwork - (sn << 1)) / (sn + 1);
/* Computing MAX */
#line 305 "sgeqp3.f"
		    i__1 = 2, i__2 = ilaenv_(&c__2, "SGEQRF", " ", &sm, &sn, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 305 "sgeqp3.f"
		    nbmin = max(i__1,i__2);


#line 309 "sgeqp3.f"
		}
#line 310 "sgeqp3.f"
	    }
#line 311 "sgeqp3.f"
	}

/*        Initialize partial column norms. The first N elements of work */
/*        store the exact column norms. */

#line 316 "sgeqp3.f"
	i__1 = *n;
#line 316 "sgeqp3.f"
	for (j = nfxd + 1; j <= i__1; ++j) {
#line 317 "sgeqp3.f"
	    work[j] = snrm2_(&sm, &a[nfxd + 1 + j * a_dim1], &c__1);
#line 318 "sgeqp3.f"
	    work[*n + j] = work[j];
#line 319 "sgeqp3.f"
/* L20: */
#line 319 "sgeqp3.f"
	}

#line 321 "sgeqp3.f"
	if (nb >= nbmin && nb < sminmn && nx < sminmn) {

/*           Use blocked code initially. */

#line 326 "sgeqp3.f"
	    j = nfxd + 1;

/*           Compute factorization: while loop. */


#line 331 "sgeqp3.f"
	    topbmn = minmn - nx;
#line 332 "sgeqp3.f"
L30:
#line 333 "sgeqp3.f"
	    if (j <= topbmn) {
/* Computing MIN */
#line 334 "sgeqp3.f"
		i__1 = nb, i__2 = topbmn - j + 1;
#line 334 "sgeqp3.f"
		jb = min(i__1,i__2);

/*              Factorize JB columns among columns J:N. */

#line 338 "sgeqp3.f"
		i__1 = *n - j + 1;
#line 338 "sgeqp3.f"
		i__2 = j - 1;
#line 338 "sgeqp3.f"
		i__3 = *n - j + 1;
#line 338 "sgeqp3.f"
		slaqps_(m, &i__1, &i__2, &jb, &fjb, &a[j * a_dim1 + 1], lda, &
			jpvt[j], &tau[j], &work[j], &work[*n + j], &work[(*n 
			<< 1) + 1], &work[(*n << 1) + jb + 1], &i__3);

#line 342 "sgeqp3.f"
		j += fjb;
#line 343 "sgeqp3.f"
		goto L30;
#line 344 "sgeqp3.f"
	    }
#line 345 "sgeqp3.f"
	} else {
#line 346 "sgeqp3.f"
	    j = nfxd + 1;
#line 347 "sgeqp3.f"
	}

/*        Use unblocked code to factor the last or only block. */


#line 352 "sgeqp3.f"
	if (j <= minmn) {
#line 352 "sgeqp3.f"
	    i__1 = *n - j + 1;
#line 352 "sgeqp3.f"
	    i__2 = j - 1;
#line 352 "sgeqp3.f"
	    slaqp2_(m, &i__1, &i__2, &a[j * a_dim1 + 1], lda, &jpvt[j], &tau[
		    j], &work[j], &work[*n + j], &work[(*n << 1) + 1]);
#line 352 "sgeqp3.f"
	}

#line 357 "sgeqp3.f"
    }

#line 359 "sgeqp3.f"
    work[1] = (doublereal) iws;
#line 360 "sgeqp3.f"
    return 0;

/*     End of SGEQP3 */

} /* sgeqp3_ */

