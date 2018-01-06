#line 1 "dtrsyl.f"
/* dtrsyl.f -- translated by f2c (version 20100827).
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

#line 1 "dtrsyl.f"
/* Table of constant values */

static integer c__1 = 1;
static logical c_false = FALSE_;
static integer c__2 = 2;
static doublereal c_b26 = 1.;
static doublereal c_b30 = 0.;
static logical c_true = TRUE_;

/* > \brief \b DTRSYL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTRSYL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrsyl.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrsyl.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrsyl.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, */
/*                          LDC, SCALE, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANA, TRANB */
/*       INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N */
/*       DOUBLE PRECISION   SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTRSYL solves the real Sylvester matrix equation: */
/* > */
/* >    op(A)*X + X*op(B) = scale*C or */
/* >    op(A)*X - X*op(B) = scale*C, */
/* > */
/* > where op(A) = A or A**T, and  A and B are both upper quasi- */
/* > triangular. A is M-by-M and B is N-by-N; the right hand side C and */
/* > the solution X are M-by-N; and scale is an output scale factor, set */
/* > <= 1 to avoid overflow in X. */
/* > */
/* > A and B must be in Schur canonical form (as returned by DHSEQR), that */
/* > is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; */
/* > each 2-by-2 diagonal block has its diagonal elements equal and its */
/* > off-diagonal elements of opposite sign. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANA */
/* > \verbatim */
/* >          TRANA is CHARACTER*1 */
/* >          Specifies the option op(A): */
/* >          = 'N': op(A) = A    (No transpose) */
/* >          = 'T': op(A) = A**T (Transpose) */
/* >          = 'C': op(A) = A**H (Conjugate transpose = Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] TRANB */
/* > \verbatim */
/* >          TRANB is CHARACTER*1 */
/* >          Specifies the option op(B): */
/* >          = 'N': op(B) = B    (No transpose) */
/* >          = 'T': op(B) = B**T (Transpose) */
/* >          = 'C': op(B) = B**H (Conjugate transpose = Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] ISGN */
/* > \verbatim */
/* >          ISGN is INTEGER */
/* >          Specifies the sign in the equation: */
/* >          = +1: solve op(A)*X + X*op(B) = scale*C */
/* >          = -1: solve op(A)*X - X*op(B) = scale*C */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The order of the matrix A, and the number of rows in the */
/* >          matrices X and C. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix B, and the number of columns in the */
/* >          matrices X and C. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,M) */
/* >          The upper quasi-triangular matrix A, in Schur canonical form. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,N) */
/* >          The upper quasi-triangular matrix B, in Schur canonical form. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (LDC,N) */
/* >          On entry, the M-by-N right hand side matrix C. */
/* >          On exit, C is overwritten by the solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. LDC >= max(1,M) */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is DOUBLE PRECISION */
/* >          The scale factor, scale, set <= 1 to avoid overflow in X. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* >          = 1: A and B have common or very close eigenvalues; perturbed */
/* >               values were used to solve the equation (but the matrices */
/* >               A and B are unchanged). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup doubleSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int dtrsyl_(char *trana, char *tranb, integer *isgn, integer 
	*m, integer *n, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *c__, integer *ldc, doublereal *scale, integer *info, 
	ftnlen trana_len, ftnlen tranb_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j, k, l;
    static doublereal x[4]	/* was [2][2] */;
    static integer k1, k2, l1, l2;
    static doublereal a11, db, da11, vec[4]	/* was [2][2] */, dum[1], eps,
	     sgn;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ierr;
    static doublereal smin, suml, sumr;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer knext, lnext;
    static doublereal xnorm;
    extern /* Subroutine */ int dlaln2_(logical *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, doublereal *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *),
	     dlasy2_(logical *, logical *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *), dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    static doublereal scaloc;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static logical notrna, notrnb;
    static doublereal smlnum;


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test input parameters */

#line 211 "dtrsyl.f"
    /* Parameter adjustments */
#line 211 "dtrsyl.f"
    a_dim1 = *lda;
#line 211 "dtrsyl.f"
    a_offset = 1 + a_dim1;
#line 211 "dtrsyl.f"
    a -= a_offset;
#line 211 "dtrsyl.f"
    b_dim1 = *ldb;
#line 211 "dtrsyl.f"
    b_offset = 1 + b_dim1;
#line 211 "dtrsyl.f"
    b -= b_offset;
#line 211 "dtrsyl.f"
    c_dim1 = *ldc;
#line 211 "dtrsyl.f"
    c_offset = 1 + c_dim1;
#line 211 "dtrsyl.f"
    c__ -= c_offset;
#line 211 "dtrsyl.f"

#line 211 "dtrsyl.f"
    /* Function Body */
#line 211 "dtrsyl.f"
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
#line 212 "dtrsyl.f"
    notrnb = lsame_(tranb, "N", (ftnlen)1, (ftnlen)1);

#line 214 "dtrsyl.f"
    *info = 0;
#line 215 "dtrsyl.f"
    if (! notrna && ! lsame_(trana, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trana, "C", (ftnlen)1, (ftnlen)1)) {
#line 217 "dtrsyl.f"
	*info = -1;
#line 218 "dtrsyl.f"
    } else if (! notrnb && ! lsame_(tranb, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(tranb, "C", (ftnlen)1, (ftnlen)1)) {
#line 220 "dtrsyl.f"
	*info = -2;
#line 221 "dtrsyl.f"
    } else if (*isgn != 1 && *isgn != -1) {
#line 222 "dtrsyl.f"
	*info = -3;
#line 223 "dtrsyl.f"
    } else if (*m < 0) {
#line 224 "dtrsyl.f"
	*info = -4;
#line 225 "dtrsyl.f"
    } else if (*n < 0) {
#line 226 "dtrsyl.f"
	*info = -5;
#line 227 "dtrsyl.f"
    } else if (*lda < max(1,*m)) {
#line 228 "dtrsyl.f"
	*info = -7;
#line 229 "dtrsyl.f"
    } else if (*ldb < max(1,*n)) {
#line 230 "dtrsyl.f"
	*info = -9;
#line 231 "dtrsyl.f"
    } else if (*ldc < max(1,*m)) {
#line 232 "dtrsyl.f"
	*info = -11;
#line 233 "dtrsyl.f"
    }
#line 234 "dtrsyl.f"
    if (*info != 0) {
#line 235 "dtrsyl.f"
	i__1 = -(*info);
#line 235 "dtrsyl.f"
	xerbla_("DTRSYL", &i__1, (ftnlen)6);
#line 236 "dtrsyl.f"
	return 0;
#line 237 "dtrsyl.f"
    }

/*     Quick return if possible */

#line 241 "dtrsyl.f"
    *scale = 1.;
#line 242 "dtrsyl.f"
    if (*m == 0 || *n == 0) {
#line 242 "dtrsyl.f"
	return 0;
#line 242 "dtrsyl.f"
    }

/*     Set constants to control overflow */

#line 247 "dtrsyl.f"
    eps = dlamch_("P", (ftnlen)1);
#line 248 "dtrsyl.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 249 "dtrsyl.f"
    bignum = 1. / smlnum;
#line 250 "dtrsyl.f"
    dlabad_(&smlnum, &bignum);
#line 251 "dtrsyl.f"
    smlnum = smlnum * (doublereal) (*m * *n) / eps;
#line 252 "dtrsyl.f"
    bignum = 1. / smlnum;

/* Computing MAX */
#line 254 "dtrsyl.f"
    d__1 = smlnum, d__2 = eps * dlange_("M", m, m, &a[a_offset], lda, dum, (
	    ftnlen)1), d__1 = max(d__1,d__2), d__2 = eps * dlange_("M", n, n, 
	    &b[b_offset], ldb, dum, (ftnlen)1);
#line 254 "dtrsyl.f"
    smin = max(d__1,d__2);

#line 257 "dtrsyl.f"
    sgn = (doublereal) (*isgn);

#line 259 "dtrsyl.f"
    if (notrna && notrnb) {

/*        Solve    A*X + ISGN*X*B = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        bottom-left corner column by column by */

/*         A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L) */

/*        Where */
/*                  M                         L-1 */
/*        R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)]. */
/*                I=K+1                       J=1 */

/*        Start column loop (index = L) */
/*        L1 (L2) : column index of the first (first) row of X(K,L). */

#line 276 "dtrsyl.f"
	lnext = 1;
#line 277 "dtrsyl.f"
	i__1 = *n;
#line 277 "dtrsyl.f"
	for (l = 1; l <= i__1; ++l) {
#line 278 "dtrsyl.f"
	    if (l < lnext) {
#line 278 "dtrsyl.f"
		goto L60;
#line 278 "dtrsyl.f"
	    }
#line 280 "dtrsyl.f"
	    if (l == *n) {
#line 281 "dtrsyl.f"
		l1 = l;
#line 282 "dtrsyl.f"
		l2 = l;
#line 283 "dtrsyl.f"
	    } else {
#line 284 "dtrsyl.f"
		if (b[l + 1 + l * b_dim1] != 0.) {
#line 285 "dtrsyl.f"
		    l1 = l;
#line 286 "dtrsyl.f"
		    l2 = l + 1;
#line 287 "dtrsyl.f"
		    lnext = l + 2;
#line 288 "dtrsyl.f"
		} else {
#line 289 "dtrsyl.f"
		    l1 = l;
#line 290 "dtrsyl.f"
		    l2 = l;
#line 291 "dtrsyl.f"
		    lnext = l + 1;
#line 292 "dtrsyl.f"
		}
#line 293 "dtrsyl.f"
	    }

/*           Start row loop (index = K) */
/*           K1 (K2): row index of the first (last) row of X(K,L). */

#line 298 "dtrsyl.f"
	    knext = *m;
#line 299 "dtrsyl.f"
	    for (k = *m; k >= 1; --k) {
#line 300 "dtrsyl.f"
		if (k > knext) {
#line 300 "dtrsyl.f"
		    goto L50;
#line 300 "dtrsyl.f"
		}
#line 302 "dtrsyl.f"
		if (k == 1) {
#line 303 "dtrsyl.f"
		    k1 = k;
#line 304 "dtrsyl.f"
		    k2 = k;
#line 305 "dtrsyl.f"
		} else {
#line 306 "dtrsyl.f"
		    if (a[k + (k - 1) * a_dim1] != 0.) {
#line 307 "dtrsyl.f"
			k1 = k - 1;
#line 308 "dtrsyl.f"
			k2 = k;
#line 309 "dtrsyl.f"
			knext = k - 2;
#line 310 "dtrsyl.f"
		    } else {
#line 311 "dtrsyl.f"
			k1 = k;
#line 312 "dtrsyl.f"
			k2 = k;
#line 313 "dtrsyl.f"
			knext = k - 1;
#line 314 "dtrsyl.f"
		    }
#line 315 "dtrsyl.f"
		}

#line 317 "dtrsyl.f"
		if (l1 == l2 && k1 == k2) {
#line 318 "dtrsyl.f"
		    i__2 = *m - k1;
/* Computing MIN */
#line 318 "dtrsyl.f"
		    i__3 = k1 + 1;
/* Computing MIN */
#line 318 "dtrsyl.f"
		    i__4 = k1 + 1;
#line 318 "dtrsyl.f"
		    suml = ddot_(&i__2, &a[k1 + min(i__3,*m) * a_dim1], lda, &
			    c__[min(i__4,*m) + l1 * c_dim1], &c__1);
#line 320 "dtrsyl.f"
		    i__2 = l1 - 1;
#line 320 "dtrsyl.f"
		    sumr = ddot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 321 "dtrsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
#line 322 "dtrsyl.f"
		    scaloc = 1.;

#line 324 "dtrsyl.f"
		    a11 = a[k1 + k1 * a_dim1] + sgn * b[l1 + l1 * b_dim1];
#line 325 "dtrsyl.f"
		    da11 = abs(a11);
#line 326 "dtrsyl.f"
		    if (da11 <= smin) {
#line 327 "dtrsyl.f"
			a11 = smin;
#line 328 "dtrsyl.f"
			da11 = smin;
#line 329 "dtrsyl.f"
			*info = 1;
#line 330 "dtrsyl.f"
		    }
#line 331 "dtrsyl.f"
		    db = abs(vec[0]);
#line 332 "dtrsyl.f"
		    if (da11 < 1. && db > 1.) {
#line 333 "dtrsyl.f"
			if (db > bignum * da11) {
#line 333 "dtrsyl.f"
			    scaloc = 1. / db;
#line 333 "dtrsyl.f"
			}
#line 335 "dtrsyl.f"
		    }
#line 336 "dtrsyl.f"
		    x[0] = vec[0] * scaloc / a11;

#line 338 "dtrsyl.f"
		    if (scaloc != 1.) {
#line 339 "dtrsyl.f"
			i__2 = *n;
#line 339 "dtrsyl.f"
			for (j = 1; j <= i__2; ++j) {
#line 340 "dtrsyl.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 341 "dtrsyl.f"
/* L10: */
#line 341 "dtrsyl.f"
			}
#line 342 "dtrsyl.f"
			*scale *= scaloc;
#line 343 "dtrsyl.f"
		    }
#line 344 "dtrsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];

#line 346 "dtrsyl.f"
		} else if (l1 == l2 && k1 != k2) {

#line 348 "dtrsyl.f"
		    i__2 = *m - k2;
/* Computing MIN */
#line 348 "dtrsyl.f"
		    i__3 = k2 + 1;
/* Computing MIN */
#line 348 "dtrsyl.f"
		    i__4 = k2 + 1;
#line 348 "dtrsyl.f"
		    suml = ddot_(&i__2, &a[k1 + min(i__3,*m) * a_dim1], lda, &
			    c__[min(i__4,*m) + l1 * c_dim1], &c__1);
#line 350 "dtrsyl.f"
		    i__2 = l1 - 1;
#line 350 "dtrsyl.f"
		    sumr = ddot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 351 "dtrsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

#line 353 "dtrsyl.f"
		    i__2 = *m - k2;
/* Computing MIN */
#line 353 "dtrsyl.f"
		    i__3 = k2 + 1;
/* Computing MIN */
#line 353 "dtrsyl.f"
		    i__4 = k2 + 1;
#line 353 "dtrsyl.f"
		    suml = ddot_(&i__2, &a[k2 + min(i__3,*m) * a_dim1], lda, &
			    c__[min(i__4,*m) + l1 * c_dim1], &c__1);
#line 355 "dtrsyl.f"
		    i__2 = l1 - 1;
#line 355 "dtrsyl.f"
		    sumr = ddot_(&i__2, &c__[k2 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 356 "dtrsyl.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

#line 358 "dtrsyl.f"
		    d__1 = -sgn * b[l1 + l1 * b_dim1];
#line 358 "dtrsyl.f"
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &a[k1 + k1 
			    * a_dim1], lda, &c_b26, &c_b26, vec, &c__2, &d__1,
			     &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
#line 361 "dtrsyl.f"
		    if (ierr != 0) {
#line 361 "dtrsyl.f"
			*info = 1;
#line 361 "dtrsyl.f"
		    }

#line 364 "dtrsyl.f"
		    if (scaloc != 1.) {
#line 365 "dtrsyl.f"
			i__2 = *n;
#line 365 "dtrsyl.f"
			for (j = 1; j <= i__2; ++j) {
#line 366 "dtrsyl.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 367 "dtrsyl.f"
/* L20: */
#line 367 "dtrsyl.f"
			}
#line 368 "dtrsyl.f"
			*scale *= scaloc;
#line 369 "dtrsyl.f"
		    }
#line 370 "dtrsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 371 "dtrsyl.f"
		    c__[k2 + l1 * c_dim1] = x[1];

#line 373 "dtrsyl.f"
		} else if (l1 != l2 && k1 == k2) {

#line 375 "dtrsyl.f"
		    i__2 = *m - k1;
/* Computing MIN */
#line 375 "dtrsyl.f"
		    i__3 = k1 + 1;
/* Computing MIN */
#line 375 "dtrsyl.f"
		    i__4 = k1 + 1;
#line 375 "dtrsyl.f"
		    suml = ddot_(&i__2, &a[k1 + min(i__3,*m) * a_dim1], lda, &
			    c__[min(i__4,*m) + l1 * c_dim1], &c__1);
#line 377 "dtrsyl.f"
		    i__2 = l1 - 1;
#line 377 "dtrsyl.f"
		    sumr = ddot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 378 "dtrsyl.f"
		    vec[0] = sgn * (c__[k1 + l1 * c_dim1] - (suml + sgn * 
			    sumr));

#line 380 "dtrsyl.f"
		    i__2 = *m - k1;
/* Computing MIN */
#line 380 "dtrsyl.f"
		    i__3 = k1 + 1;
/* Computing MIN */
#line 380 "dtrsyl.f"
		    i__4 = k1 + 1;
#line 380 "dtrsyl.f"
		    suml = ddot_(&i__2, &a[k1 + min(i__3,*m) * a_dim1], lda, &
			    c__[min(i__4,*m) + l2 * c_dim1], &c__1);
#line 382 "dtrsyl.f"
		    i__2 = l1 - 1;
#line 382 "dtrsyl.f"
		    sumr = ddot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l2 * 
			    b_dim1 + 1], &c__1);
#line 383 "dtrsyl.f"
		    vec[1] = sgn * (c__[k1 + l2 * c_dim1] - (suml + sgn * 
			    sumr));

#line 385 "dtrsyl.f"
		    d__1 = -sgn * a[k1 + k1 * a_dim1];
#line 385 "dtrsyl.f"
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &b[l1 + l1 *
			     b_dim1], ldb, &c_b26, &c_b26, vec, &c__2, &d__1, 
			    &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
#line 388 "dtrsyl.f"
		    if (ierr != 0) {
#line 388 "dtrsyl.f"
			*info = 1;
#line 388 "dtrsyl.f"
		    }

#line 391 "dtrsyl.f"
		    if (scaloc != 1.) {
#line 392 "dtrsyl.f"
			i__2 = *n;
#line 392 "dtrsyl.f"
			for (j = 1; j <= i__2; ++j) {
#line 393 "dtrsyl.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 394 "dtrsyl.f"
/* L30: */
#line 394 "dtrsyl.f"
			}
#line 395 "dtrsyl.f"
			*scale *= scaloc;
#line 396 "dtrsyl.f"
		    }
#line 397 "dtrsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 398 "dtrsyl.f"
		    c__[k1 + l2 * c_dim1] = x[1];

#line 400 "dtrsyl.f"
		} else if (l1 != l2 && k1 != k2) {

#line 402 "dtrsyl.f"
		    i__2 = *m - k2;
/* Computing MIN */
#line 402 "dtrsyl.f"
		    i__3 = k2 + 1;
/* Computing MIN */
#line 402 "dtrsyl.f"
		    i__4 = k2 + 1;
#line 402 "dtrsyl.f"
		    suml = ddot_(&i__2, &a[k1 + min(i__3,*m) * a_dim1], lda, &
			    c__[min(i__4,*m) + l1 * c_dim1], &c__1);
#line 404 "dtrsyl.f"
		    i__2 = l1 - 1;
#line 404 "dtrsyl.f"
		    sumr = ddot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 405 "dtrsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

#line 407 "dtrsyl.f"
		    i__2 = *m - k2;
/* Computing MIN */
#line 407 "dtrsyl.f"
		    i__3 = k2 + 1;
/* Computing MIN */
#line 407 "dtrsyl.f"
		    i__4 = k2 + 1;
#line 407 "dtrsyl.f"
		    suml = ddot_(&i__2, &a[k1 + min(i__3,*m) * a_dim1], lda, &
			    c__[min(i__4,*m) + l2 * c_dim1], &c__1);
#line 409 "dtrsyl.f"
		    i__2 = l1 - 1;
#line 409 "dtrsyl.f"
		    sumr = ddot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l2 * 
			    b_dim1 + 1], &c__1);
#line 410 "dtrsyl.f"
		    vec[2] = c__[k1 + l2 * c_dim1] - (suml + sgn * sumr);

#line 412 "dtrsyl.f"
		    i__2 = *m - k2;
/* Computing MIN */
#line 412 "dtrsyl.f"
		    i__3 = k2 + 1;
/* Computing MIN */
#line 412 "dtrsyl.f"
		    i__4 = k2 + 1;
#line 412 "dtrsyl.f"
		    suml = ddot_(&i__2, &a[k2 + min(i__3,*m) * a_dim1], lda, &
			    c__[min(i__4,*m) + l1 * c_dim1], &c__1);
#line 414 "dtrsyl.f"
		    i__2 = l1 - 1;
#line 414 "dtrsyl.f"
		    sumr = ddot_(&i__2, &c__[k2 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 415 "dtrsyl.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

#line 417 "dtrsyl.f"
		    i__2 = *m - k2;
/* Computing MIN */
#line 417 "dtrsyl.f"
		    i__3 = k2 + 1;
/* Computing MIN */
#line 417 "dtrsyl.f"
		    i__4 = k2 + 1;
#line 417 "dtrsyl.f"
		    suml = ddot_(&i__2, &a[k2 + min(i__3,*m) * a_dim1], lda, &
			    c__[min(i__4,*m) + l2 * c_dim1], &c__1);
#line 419 "dtrsyl.f"
		    i__2 = l1 - 1;
#line 419 "dtrsyl.f"
		    sumr = ddot_(&i__2, &c__[k2 + c_dim1], ldc, &b[l2 * 
			    b_dim1 + 1], &c__1);
#line 420 "dtrsyl.f"
		    vec[3] = c__[k2 + l2 * c_dim1] - (suml + sgn * sumr);

#line 422 "dtrsyl.f"
		    dlasy2_(&c_false, &c_false, isgn, &c__2, &c__2, &a[k1 + 
			    k1 * a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec,
			     &c__2, &scaloc, x, &c__2, &xnorm, &ierr);
#line 425 "dtrsyl.f"
		    if (ierr != 0) {
#line 425 "dtrsyl.f"
			*info = 1;
#line 425 "dtrsyl.f"
		    }

#line 428 "dtrsyl.f"
		    if (scaloc != 1.) {
#line 429 "dtrsyl.f"
			i__2 = *n;
#line 429 "dtrsyl.f"
			for (j = 1; j <= i__2; ++j) {
#line 430 "dtrsyl.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 431 "dtrsyl.f"
/* L40: */
#line 431 "dtrsyl.f"
			}
#line 432 "dtrsyl.f"
			*scale *= scaloc;
#line 433 "dtrsyl.f"
		    }
#line 434 "dtrsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 435 "dtrsyl.f"
		    c__[k1 + l2 * c_dim1] = x[2];
#line 436 "dtrsyl.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 437 "dtrsyl.f"
		    c__[k2 + l2 * c_dim1] = x[3];
#line 438 "dtrsyl.f"
		}

#line 440 "dtrsyl.f"
L50:
#line 440 "dtrsyl.f"
		;
#line 440 "dtrsyl.f"
	    }

#line 442 "dtrsyl.f"
L60:
#line 442 "dtrsyl.f"
	    ;
#line 442 "dtrsyl.f"
	}

#line 444 "dtrsyl.f"
    } else if (! notrna && notrnb) {

/*        Solve    A**T *X + ISGN*X*B = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        upper-left corner column by column by */

/*          A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L) */

/*        Where */
/*                   K-1        T                    L-1 */
/*          R(K,L) = SUM [A(I,K)**T*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)] */
/*                   I=1                          J=1 */

/*        Start column loop (index = L) */
/*        L1 (L2): column index of the first (last) row of X(K,L) */

#line 461 "dtrsyl.f"
	lnext = 1;
#line 462 "dtrsyl.f"
	i__1 = *n;
#line 462 "dtrsyl.f"
	for (l = 1; l <= i__1; ++l) {
#line 463 "dtrsyl.f"
	    if (l < lnext) {
#line 463 "dtrsyl.f"
		goto L120;
#line 463 "dtrsyl.f"
	    }
#line 465 "dtrsyl.f"
	    if (l == *n) {
#line 466 "dtrsyl.f"
		l1 = l;
#line 467 "dtrsyl.f"
		l2 = l;
#line 468 "dtrsyl.f"
	    } else {
#line 469 "dtrsyl.f"
		if (b[l + 1 + l * b_dim1] != 0.) {
#line 470 "dtrsyl.f"
		    l1 = l;
#line 471 "dtrsyl.f"
		    l2 = l + 1;
#line 472 "dtrsyl.f"
		    lnext = l + 2;
#line 473 "dtrsyl.f"
		} else {
#line 474 "dtrsyl.f"
		    l1 = l;
#line 475 "dtrsyl.f"
		    l2 = l;
#line 476 "dtrsyl.f"
		    lnext = l + 1;
#line 477 "dtrsyl.f"
		}
#line 478 "dtrsyl.f"
	    }

/*           Start row loop (index = K) */
/*           K1 (K2): row index of the first (last) row of X(K,L) */

#line 483 "dtrsyl.f"
	    knext = 1;
#line 484 "dtrsyl.f"
	    i__2 = *m;
#line 484 "dtrsyl.f"
	    for (k = 1; k <= i__2; ++k) {
#line 485 "dtrsyl.f"
		if (k < knext) {
#line 485 "dtrsyl.f"
		    goto L110;
#line 485 "dtrsyl.f"
		}
#line 487 "dtrsyl.f"
		if (k == *m) {
#line 488 "dtrsyl.f"
		    k1 = k;
#line 489 "dtrsyl.f"
		    k2 = k;
#line 490 "dtrsyl.f"
		} else {
#line 491 "dtrsyl.f"
		    if (a[k + 1 + k * a_dim1] != 0.) {
#line 492 "dtrsyl.f"
			k1 = k;
#line 493 "dtrsyl.f"
			k2 = k + 1;
#line 494 "dtrsyl.f"
			knext = k + 2;
#line 495 "dtrsyl.f"
		    } else {
#line 496 "dtrsyl.f"
			k1 = k;
#line 497 "dtrsyl.f"
			k2 = k;
#line 498 "dtrsyl.f"
			knext = k + 1;
#line 499 "dtrsyl.f"
		    }
#line 500 "dtrsyl.f"
		}

#line 502 "dtrsyl.f"
		if (l1 == l2 && k1 == k2) {
#line 503 "dtrsyl.f"
		    i__3 = k1 - 1;
#line 503 "dtrsyl.f"
		    suml = ddot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 504 "dtrsyl.f"
		    i__3 = l1 - 1;
#line 504 "dtrsyl.f"
		    sumr = ddot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 505 "dtrsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
#line 506 "dtrsyl.f"
		    scaloc = 1.;

#line 508 "dtrsyl.f"
		    a11 = a[k1 + k1 * a_dim1] + sgn * b[l1 + l1 * b_dim1];
#line 509 "dtrsyl.f"
		    da11 = abs(a11);
#line 510 "dtrsyl.f"
		    if (da11 <= smin) {
#line 511 "dtrsyl.f"
			a11 = smin;
#line 512 "dtrsyl.f"
			da11 = smin;
#line 513 "dtrsyl.f"
			*info = 1;
#line 514 "dtrsyl.f"
		    }
#line 515 "dtrsyl.f"
		    db = abs(vec[0]);
#line 516 "dtrsyl.f"
		    if (da11 < 1. && db > 1.) {
#line 517 "dtrsyl.f"
			if (db > bignum * da11) {
#line 517 "dtrsyl.f"
			    scaloc = 1. / db;
#line 517 "dtrsyl.f"
			}
#line 519 "dtrsyl.f"
		    }
#line 520 "dtrsyl.f"
		    x[0] = vec[0] * scaloc / a11;

#line 522 "dtrsyl.f"
		    if (scaloc != 1.) {
#line 523 "dtrsyl.f"
			i__3 = *n;
#line 523 "dtrsyl.f"
			for (j = 1; j <= i__3; ++j) {
#line 524 "dtrsyl.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 525 "dtrsyl.f"
/* L70: */
#line 525 "dtrsyl.f"
			}
#line 526 "dtrsyl.f"
			*scale *= scaloc;
#line 527 "dtrsyl.f"
		    }
#line 528 "dtrsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];

#line 530 "dtrsyl.f"
		} else if (l1 == l2 && k1 != k2) {

#line 532 "dtrsyl.f"
		    i__3 = k1 - 1;
#line 532 "dtrsyl.f"
		    suml = ddot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 533 "dtrsyl.f"
		    i__3 = l1 - 1;
#line 533 "dtrsyl.f"
		    sumr = ddot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 534 "dtrsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

#line 536 "dtrsyl.f"
		    i__3 = k1 - 1;
#line 536 "dtrsyl.f"
		    suml = ddot_(&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 537 "dtrsyl.f"
		    i__3 = l1 - 1;
#line 537 "dtrsyl.f"
		    sumr = ddot_(&i__3, &c__[k2 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 538 "dtrsyl.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

#line 540 "dtrsyl.f"
		    d__1 = -sgn * b[l1 + l1 * b_dim1];
#line 540 "dtrsyl.f"
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &a[k1 + k1 *
			     a_dim1], lda, &c_b26, &c_b26, vec, &c__2, &d__1, 
			    &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
#line 543 "dtrsyl.f"
		    if (ierr != 0) {
#line 543 "dtrsyl.f"
			*info = 1;
#line 543 "dtrsyl.f"
		    }

#line 546 "dtrsyl.f"
		    if (scaloc != 1.) {
#line 547 "dtrsyl.f"
			i__3 = *n;
#line 547 "dtrsyl.f"
			for (j = 1; j <= i__3; ++j) {
#line 548 "dtrsyl.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 549 "dtrsyl.f"
/* L80: */
#line 549 "dtrsyl.f"
			}
#line 550 "dtrsyl.f"
			*scale *= scaloc;
#line 551 "dtrsyl.f"
		    }
#line 552 "dtrsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 553 "dtrsyl.f"
		    c__[k2 + l1 * c_dim1] = x[1];

#line 555 "dtrsyl.f"
		} else if (l1 != l2 && k1 == k2) {

#line 557 "dtrsyl.f"
		    i__3 = k1 - 1;
#line 557 "dtrsyl.f"
		    suml = ddot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 558 "dtrsyl.f"
		    i__3 = l1 - 1;
#line 558 "dtrsyl.f"
		    sumr = ddot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 559 "dtrsyl.f"
		    vec[0] = sgn * (c__[k1 + l1 * c_dim1] - (suml + sgn * 
			    sumr));

#line 561 "dtrsyl.f"
		    i__3 = k1 - 1;
#line 561 "dtrsyl.f"
		    suml = ddot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
#line 562 "dtrsyl.f"
		    i__3 = l1 - 1;
#line 562 "dtrsyl.f"
		    sumr = ddot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l2 * 
			    b_dim1 + 1], &c__1);
#line 563 "dtrsyl.f"
		    vec[1] = sgn * (c__[k1 + l2 * c_dim1] - (suml + sgn * 
			    sumr));

#line 565 "dtrsyl.f"
		    d__1 = -sgn * a[k1 + k1 * a_dim1];
#line 565 "dtrsyl.f"
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &b[l1 + l1 *
			     b_dim1], ldb, &c_b26, &c_b26, vec, &c__2, &d__1, 
			    &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
#line 568 "dtrsyl.f"
		    if (ierr != 0) {
#line 568 "dtrsyl.f"
			*info = 1;
#line 568 "dtrsyl.f"
		    }

#line 571 "dtrsyl.f"
		    if (scaloc != 1.) {
#line 572 "dtrsyl.f"
			i__3 = *n;
#line 572 "dtrsyl.f"
			for (j = 1; j <= i__3; ++j) {
#line 573 "dtrsyl.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 574 "dtrsyl.f"
/* L90: */
#line 574 "dtrsyl.f"
			}
#line 575 "dtrsyl.f"
			*scale *= scaloc;
#line 576 "dtrsyl.f"
		    }
#line 577 "dtrsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 578 "dtrsyl.f"
		    c__[k1 + l2 * c_dim1] = x[1];

#line 580 "dtrsyl.f"
		} else if (l1 != l2 && k1 != k2) {

#line 582 "dtrsyl.f"
		    i__3 = k1 - 1;
#line 582 "dtrsyl.f"
		    suml = ddot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 583 "dtrsyl.f"
		    i__3 = l1 - 1;
#line 583 "dtrsyl.f"
		    sumr = ddot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 584 "dtrsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

#line 586 "dtrsyl.f"
		    i__3 = k1 - 1;
#line 586 "dtrsyl.f"
		    suml = ddot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
#line 587 "dtrsyl.f"
		    i__3 = l1 - 1;
#line 587 "dtrsyl.f"
		    sumr = ddot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l2 * 
			    b_dim1 + 1], &c__1);
#line 588 "dtrsyl.f"
		    vec[2] = c__[k1 + l2 * c_dim1] - (suml + sgn * sumr);

#line 590 "dtrsyl.f"
		    i__3 = k1 - 1;
#line 590 "dtrsyl.f"
		    suml = ddot_(&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 591 "dtrsyl.f"
		    i__3 = l1 - 1;
#line 591 "dtrsyl.f"
		    sumr = ddot_(&i__3, &c__[k2 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 592 "dtrsyl.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

#line 594 "dtrsyl.f"
		    i__3 = k1 - 1;
#line 594 "dtrsyl.f"
		    suml = ddot_(&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
#line 595 "dtrsyl.f"
		    i__3 = l1 - 1;
#line 595 "dtrsyl.f"
		    sumr = ddot_(&i__3, &c__[k2 + c_dim1], ldc, &b[l2 * 
			    b_dim1 + 1], &c__1);
#line 596 "dtrsyl.f"
		    vec[3] = c__[k2 + l2 * c_dim1] - (suml + sgn * sumr);

#line 598 "dtrsyl.f"
		    dlasy2_(&c_true, &c_false, isgn, &c__2, &c__2, &a[k1 + k1 
			    * a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec, &
			    c__2, &scaloc, x, &c__2, &xnorm, &ierr);
#line 601 "dtrsyl.f"
		    if (ierr != 0) {
#line 601 "dtrsyl.f"
			*info = 1;
#line 601 "dtrsyl.f"
		    }

#line 604 "dtrsyl.f"
		    if (scaloc != 1.) {
#line 605 "dtrsyl.f"
			i__3 = *n;
#line 605 "dtrsyl.f"
			for (j = 1; j <= i__3; ++j) {
#line 606 "dtrsyl.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 607 "dtrsyl.f"
/* L100: */
#line 607 "dtrsyl.f"
			}
#line 608 "dtrsyl.f"
			*scale *= scaloc;
#line 609 "dtrsyl.f"
		    }
#line 610 "dtrsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 611 "dtrsyl.f"
		    c__[k1 + l2 * c_dim1] = x[2];
#line 612 "dtrsyl.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 613 "dtrsyl.f"
		    c__[k2 + l2 * c_dim1] = x[3];
#line 614 "dtrsyl.f"
		}

#line 616 "dtrsyl.f"
L110:
#line 616 "dtrsyl.f"
		;
#line 616 "dtrsyl.f"
	    }
#line 617 "dtrsyl.f"
L120:
#line 617 "dtrsyl.f"
	    ;
#line 617 "dtrsyl.f"
	}

#line 619 "dtrsyl.f"
    } else if (! notrna && ! notrnb) {

/*        Solve    A**T*X + ISGN*X*B**T = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        top-right corner column by column by */

/*           A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L) */

/*        Where */
/*                     K-1                            N */
/*            R(K,L) = SUM [A(I,K)**T*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T]. */
/*                     I=1                          J=L+1 */

/*        Start column loop (index = L) */
/*        L1 (L2): column index of the first (last) row of X(K,L) */

#line 636 "dtrsyl.f"
	lnext = *n;
#line 637 "dtrsyl.f"
	for (l = *n; l >= 1; --l) {
#line 638 "dtrsyl.f"
	    if (l > lnext) {
#line 638 "dtrsyl.f"
		goto L180;
#line 638 "dtrsyl.f"
	    }
#line 640 "dtrsyl.f"
	    if (l == 1) {
#line 641 "dtrsyl.f"
		l1 = l;
#line 642 "dtrsyl.f"
		l2 = l;
#line 643 "dtrsyl.f"
	    } else {
#line 644 "dtrsyl.f"
		if (b[l + (l - 1) * b_dim1] != 0.) {
#line 645 "dtrsyl.f"
		    l1 = l - 1;
#line 646 "dtrsyl.f"
		    l2 = l;
#line 647 "dtrsyl.f"
		    lnext = l - 2;
#line 648 "dtrsyl.f"
		} else {
#line 649 "dtrsyl.f"
		    l1 = l;
#line 650 "dtrsyl.f"
		    l2 = l;
#line 651 "dtrsyl.f"
		    lnext = l - 1;
#line 652 "dtrsyl.f"
		}
#line 653 "dtrsyl.f"
	    }

/*           Start row loop (index = K) */
/*           K1 (K2): row index of the first (last) row of X(K,L) */

#line 658 "dtrsyl.f"
	    knext = 1;
#line 659 "dtrsyl.f"
	    i__1 = *m;
#line 659 "dtrsyl.f"
	    for (k = 1; k <= i__1; ++k) {
#line 660 "dtrsyl.f"
		if (k < knext) {
#line 660 "dtrsyl.f"
		    goto L170;
#line 660 "dtrsyl.f"
		}
#line 662 "dtrsyl.f"
		if (k == *m) {
#line 663 "dtrsyl.f"
		    k1 = k;
#line 664 "dtrsyl.f"
		    k2 = k;
#line 665 "dtrsyl.f"
		} else {
#line 666 "dtrsyl.f"
		    if (a[k + 1 + k * a_dim1] != 0.) {
#line 667 "dtrsyl.f"
			k1 = k;
#line 668 "dtrsyl.f"
			k2 = k + 1;
#line 669 "dtrsyl.f"
			knext = k + 2;
#line 670 "dtrsyl.f"
		    } else {
#line 671 "dtrsyl.f"
			k1 = k;
#line 672 "dtrsyl.f"
			k2 = k;
#line 673 "dtrsyl.f"
			knext = k + 1;
#line 674 "dtrsyl.f"
		    }
#line 675 "dtrsyl.f"
		}

#line 677 "dtrsyl.f"
		if (l1 == l2 && k1 == k2) {
#line 678 "dtrsyl.f"
		    i__2 = k1 - 1;
#line 678 "dtrsyl.f"
		    suml = ddot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 679 "dtrsyl.f"
		    i__2 = *n - l1;
/* Computing MIN */
#line 679 "dtrsyl.f"
		    i__3 = l1 + 1;
/* Computing MIN */
#line 679 "dtrsyl.f"
		    i__4 = l1 + 1;
#line 679 "dtrsyl.f"
		    sumr = ddot_(&i__2, &c__[k1 + min(i__3,*n) * c_dim1], ldc,
			     &b[l1 + min(i__4,*n) * b_dim1], ldb);
#line 681 "dtrsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
#line 682 "dtrsyl.f"
		    scaloc = 1.;

#line 684 "dtrsyl.f"
		    a11 = a[k1 + k1 * a_dim1] + sgn * b[l1 + l1 * b_dim1];
#line 685 "dtrsyl.f"
		    da11 = abs(a11);
#line 686 "dtrsyl.f"
		    if (da11 <= smin) {
#line 687 "dtrsyl.f"
			a11 = smin;
#line 688 "dtrsyl.f"
			da11 = smin;
#line 689 "dtrsyl.f"
			*info = 1;
#line 690 "dtrsyl.f"
		    }
#line 691 "dtrsyl.f"
		    db = abs(vec[0]);
#line 692 "dtrsyl.f"
		    if (da11 < 1. && db > 1.) {
#line 693 "dtrsyl.f"
			if (db > bignum * da11) {
#line 693 "dtrsyl.f"
			    scaloc = 1. / db;
#line 693 "dtrsyl.f"
			}
#line 695 "dtrsyl.f"
		    }
#line 696 "dtrsyl.f"
		    x[0] = vec[0] * scaloc / a11;

#line 698 "dtrsyl.f"
		    if (scaloc != 1.) {
#line 699 "dtrsyl.f"
			i__2 = *n;
#line 699 "dtrsyl.f"
			for (j = 1; j <= i__2; ++j) {
#line 700 "dtrsyl.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 701 "dtrsyl.f"
/* L130: */
#line 701 "dtrsyl.f"
			}
#line 702 "dtrsyl.f"
			*scale *= scaloc;
#line 703 "dtrsyl.f"
		    }
#line 704 "dtrsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];

#line 706 "dtrsyl.f"
		} else if (l1 == l2 && k1 != k2) {

#line 708 "dtrsyl.f"
		    i__2 = k1 - 1;
#line 708 "dtrsyl.f"
		    suml = ddot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 709 "dtrsyl.f"
		    i__2 = *n - l2;
/* Computing MIN */
#line 709 "dtrsyl.f"
		    i__3 = l2 + 1;
/* Computing MIN */
#line 709 "dtrsyl.f"
		    i__4 = l2 + 1;
#line 709 "dtrsyl.f"
		    sumr = ddot_(&i__2, &c__[k1 + min(i__3,*n) * c_dim1], ldc,
			     &b[l1 + min(i__4,*n) * b_dim1], ldb);
#line 711 "dtrsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

#line 713 "dtrsyl.f"
		    i__2 = k1 - 1;
#line 713 "dtrsyl.f"
		    suml = ddot_(&i__2, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 714 "dtrsyl.f"
		    i__2 = *n - l2;
/* Computing MIN */
#line 714 "dtrsyl.f"
		    i__3 = l2 + 1;
/* Computing MIN */
#line 714 "dtrsyl.f"
		    i__4 = l2 + 1;
#line 714 "dtrsyl.f"
		    sumr = ddot_(&i__2, &c__[k2 + min(i__3,*n) * c_dim1], ldc,
			     &b[l1 + min(i__4,*n) * b_dim1], ldb);
#line 716 "dtrsyl.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

#line 718 "dtrsyl.f"
		    d__1 = -sgn * b[l1 + l1 * b_dim1];
#line 718 "dtrsyl.f"
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &a[k1 + k1 *
			     a_dim1], lda, &c_b26, &c_b26, vec, &c__2, &d__1, 
			    &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
#line 721 "dtrsyl.f"
		    if (ierr != 0) {
#line 721 "dtrsyl.f"
			*info = 1;
#line 721 "dtrsyl.f"
		    }

#line 724 "dtrsyl.f"
		    if (scaloc != 1.) {
#line 725 "dtrsyl.f"
			i__2 = *n;
#line 725 "dtrsyl.f"
			for (j = 1; j <= i__2; ++j) {
#line 726 "dtrsyl.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 727 "dtrsyl.f"
/* L140: */
#line 727 "dtrsyl.f"
			}
#line 728 "dtrsyl.f"
			*scale *= scaloc;
#line 729 "dtrsyl.f"
		    }
#line 730 "dtrsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 731 "dtrsyl.f"
		    c__[k2 + l1 * c_dim1] = x[1];

#line 733 "dtrsyl.f"
		} else if (l1 != l2 && k1 == k2) {

#line 735 "dtrsyl.f"
		    i__2 = k1 - 1;
#line 735 "dtrsyl.f"
		    suml = ddot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 736 "dtrsyl.f"
		    i__2 = *n - l2;
/* Computing MIN */
#line 736 "dtrsyl.f"
		    i__3 = l2 + 1;
/* Computing MIN */
#line 736 "dtrsyl.f"
		    i__4 = l2 + 1;
#line 736 "dtrsyl.f"
		    sumr = ddot_(&i__2, &c__[k1 + min(i__3,*n) * c_dim1], ldc,
			     &b[l1 + min(i__4,*n) * b_dim1], ldb);
#line 738 "dtrsyl.f"
		    vec[0] = sgn * (c__[k1 + l1 * c_dim1] - (suml + sgn * 
			    sumr));

#line 740 "dtrsyl.f"
		    i__2 = k1 - 1;
#line 740 "dtrsyl.f"
		    suml = ddot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
#line 741 "dtrsyl.f"
		    i__2 = *n - l2;
/* Computing MIN */
#line 741 "dtrsyl.f"
		    i__3 = l2 + 1;
/* Computing MIN */
#line 741 "dtrsyl.f"
		    i__4 = l2 + 1;
#line 741 "dtrsyl.f"
		    sumr = ddot_(&i__2, &c__[k1 + min(i__3,*n) * c_dim1], ldc,
			     &b[l2 + min(i__4,*n) * b_dim1], ldb);
#line 743 "dtrsyl.f"
		    vec[1] = sgn * (c__[k1 + l2 * c_dim1] - (suml + sgn * 
			    sumr));

#line 745 "dtrsyl.f"
		    d__1 = -sgn * a[k1 + k1 * a_dim1];
#line 745 "dtrsyl.f"
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &b[l1 + l1 
			    * b_dim1], ldb, &c_b26, &c_b26, vec, &c__2, &d__1,
			     &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
#line 748 "dtrsyl.f"
		    if (ierr != 0) {
#line 748 "dtrsyl.f"
			*info = 1;
#line 748 "dtrsyl.f"
		    }

#line 751 "dtrsyl.f"
		    if (scaloc != 1.) {
#line 752 "dtrsyl.f"
			i__2 = *n;
#line 752 "dtrsyl.f"
			for (j = 1; j <= i__2; ++j) {
#line 753 "dtrsyl.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 754 "dtrsyl.f"
/* L150: */
#line 754 "dtrsyl.f"
			}
#line 755 "dtrsyl.f"
			*scale *= scaloc;
#line 756 "dtrsyl.f"
		    }
#line 757 "dtrsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 758 "dtrsyl.f"
		    c__[k1 + l2 * c_dim1] = x[1];

#line 760 "dtrsyl.f"
		} else if (l1 != l2 && k1 != k2) {

#line 762 "dtrsyl.f"
		    i__2 = k1 - 1;
#line 762 "dtrsyl.f"
		    suml = ddot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 763 "dtrsyl.f"
		    i__2 = *n - l2;
/* Computing MIN */
#line 763 "dtrsyl.f"
		    i__3 = l2 + 1;
/* Computing MIN */
#line 763 "dtrsyl.f"
		    i__4 = l2 + 1;
#line 763 "dtrsyl.f"
		    sumr = ddot_(&i__2, &c__[k1 + min(i__3,*n) * c_dim1], ldc,
			     &b[l1 + min(i__4,*n) * b_dim1], ldb);
#line 765 "dtrsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

#line 767 "dtrsyl.f"
		    i__2 = k1 - 1;
#line 767 "dtrsyl.f"
		    suml = ddot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
#line 768 "dtrsyl.f"
		    i__2 = *n - l2;
/* Computing MIN */
#line 768 "dtrsyl.f"
		    i__3 = l2 + 1;
/* Computing MIN */
#line 768 "dtrsyl.f"
		    i__4 = l2 + 1;
#line 768 "dtrsyl.f"
		    sumr = ddot_(&i__2, &c__[k1 + min(i__3,*n) * c_dim1], ldc,
			     &b[l2 + min(i__4,*n) * b_dim1], ldb);
#line 770 "dtrsyl.f"
		    vec[2] = c__[k1 + l2 * c_dim1] - (suml + sgn * sumr);

#line 772 "dtrsyl.f"
		    i__2 = k1 - 1;
#line 772 "dtrsyl.f"
		    suml = ddot_(&i__2, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 773 "dtrsyl.f"
		    i__2 = *n - l2;
/* Computing MIN */
#line 773 "dtrsyl.f"
		    i__3 = l2 + 1;
/* Computing MIN */
#line 773 "dtrsyl.f"
		    i__4 = l2 + 1;
#line 773 "dtrsyl.f"
		    sumr = ddot_(&i__2, &c__[k2 + min(i__3,*n) * c_dim1], ldc,
			     &b[l1 + min(i__4,*n) * b_dim1], ldb);
#line 775 "dtrsyl.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

#line 777 "dtrsyl.f"
		    i__2 = k1 - 1;
#line 777 "dtrsyl.f"
		    suml = ddot_(&i__2, &a[k2 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
#line 778 "dtrsyl.f"
		    i__2 = *n - l2;
/* Computing MIN */
#line 778 "dtrsyl.f"
		    i__3 = l2 + 1;
/* Computing MIN */
#line 778 "dtrsyl.f"
		    i__4 = l2 + 1;
#line 778 "dtrsyl.f"
		    sumr = ddot_(&i__2, &c__[k2 + min(i__3,*n) * c_dim1], ldc,
			     &b[l2 + min(i__4,*n) * b_dim1], ldb);
#line 780 "dtrsyl.f"
		    vec[3] = c__[k2 + l2 * c_dim1] - (suml + sgn * sumr);

#line 782 "dtrsyl.f"
		    dlasy2_(&c_true, &c_true, isgn, &c__2, &c__2, &a[k1 + k1 *
			     a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec, &
			    c__2, &scaloc, x, &c__2, &xnorm, &ierr);
#line 785 "dtrsyl.f"
		    if (ierr != 0) {
#line 785 "dtrsyl.f"
			*info = 1;
#line 785 "dtrsyl.f"
		    }

#line 788 "dtrsyl.f"
		    if (scaloc != 1.) {
#line 789 "dtrsyl.f"
			i__2 = *n;
#line 789 "dtrsyl.f"
			for (j = 1; j <= i__2; ++j) {
#line 790 "dtrsyl.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 791 "dtrsyl.f"
/* L160: */
#line 791 "dtrsyl.f"
			}
#line 792 "dtrsyl.f"
			*scale *= scaloc;
#line 793 "dtrsyl.f"
		    }
#line 794 "dtrsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 795 "dtrsyl.f"
		    c__[k1 + l2 * c_dim1] = x[2];
#line 796 "dtrsyl.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 797 "dtrsyl.f"
		    c__[k2 + l2 * c_dim1] = x[3];
#line 798 "dtrsyl.f"
		}

#line 800 "dtrsyl.f"
L170:
#line 800 "dtrsyl.f"
		;
#line 800 "dtrsyl.f"
	    }
#line 801 "dtrsyl.f"
L180:
#line 801 "dtrsyl.f"
	    ;
#line 801 "dtrsyl.f"
	}

#line 803 "dtrsyl.f"
    } else if (notrna && ! notrnb) {

/*        Solve    A*X + ISGN*X*B**T = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        bottom-right corner column by column by */

/*            A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L) */

/*        Where */
/*                      M                          N */
/*            R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T]. */
/*                    I=K+1                      J=L+1 */

/*        Start column loop (index = L) */
/*        L1 (L2): column index of the first (last) row of X(K,L) */

#line 820 "dtrsyl.f"
	lnext = *n;
#line 821 "dtrsyl.f"
	for (l = *n; l >= 1; --l) {
#line 822 "dtrsyl.f"
	    if (l > lnext) {
#line 822 "dtrsyl.f"
		goto L240;
#line 822 "dtrsyl.f"
	    }
#line 824 "dtrsyl.f"
	    if (l == 1) {
#line 825 "dtrsyl.f"
		l1 = l;
#line 826 "dtrsyl.f"
		l2 = l;
#line 827 "dtrsyl.f"
	    } else {
#line 828 "dtrsyl.f"
		if (b[l + (l - 1) * b_dim1] != 0.) {
#line 829 "dtrsyl.f"
		    l1 = l - 1;
#line 830 "dtrsyl.f"
		    l2 = l;
#line 831 "dtrsyl.f"
		    lnext = l - 2;
#line 832 "dtrsyl.f"
		} else {
#line 833 "dtrsyl.f"
		    l1 = l;
#line 834 "dtrsyl.f"
		    l2 = l;
#line 835 "dtrsyl.f"
		    lnext = l - 1;
#line 836 "dtrsyl.f"
		}
#line 837 "dtrsyl.f"
	    }

/*           Start row loop (index = K) */
/*           K1 (K2): row index of the first (last) row of X(K,L) */

#line 842 "dtrsyl.f"
	    knext = *m;
#line 843 "dtrsyl.f"
	    for (k = *m; k >= 1; --k) {
#line 844 "dtrsyl.f"
		if (k > knext) {
#line 844 "dtrsyl.f"
		    goto L230;
#line 844 "dtrsyl.f"
		}
#line 846 "dtrsyl.f"
		if (k == 1) {
#line 847 "dtrsyl.f"
		    k1 = k;
#line 848 "dtrsyl.f"
		    k2 = k;
#line 849 "dtrsyl.f"
		} else {
#line 850 "dtrsyl.f"
		    if (a[k + (k - 1) * a_dim1] != 0.) {
#line 851 "dtrsyl.f"
			k1 = k - 1;
#line 852 "dtrsyl.f"
			k2 = k;
#line 853 "dtrsyl.f"
			knext = k - 2;
#line 854 "dtrsyl.f"
		    } else {
#line 855 "dtrsyl.f"
			k1 = k;
#line 856 "dtrsyl.f"
			k2 = k;
#line 857 "dtrsyl.f"
			knext = k - 1;
#line 858 "dtrsyl.f"
		    }
#line 859 "dtrsyl.f"
		}

#line 861 "dtrsyl.f"
		if (l1 == l2 && k1 == k2) {
#line 862 "dtrsyl.f"
		    i__1 = *m - k1;
/* Computing MIN */
#line 862 "dtrsyl.f"
		    i__2 = k1 + 1;
/* Computing MIN */
#line 862 "dtrsyl.f"
		    i__3 = k1 + 1;
#line 862 "dtrsyl.f"
		    suml = ddot_(&i__1, &a[k1 + min(i__2,*m) * a_dim1], lda, &
			    c__[min(i__3,*m) + l1 * c_dim1], &c__1);
#line 864 "dtrsyl.f"
		    i__1 = *n - l1;
/* Computing MIN */
#line 864 "dtrsyl.f"
		    i__2 = l1 + 1;
/* Computing MIN */
#line 864 "dtrsyl.f"
		    i__3 = l1 + 1;
#line 864 "dtrsyl.f"
		    sumr = ddot_(&i__1, &c__[k1 + min(i__2,*n) * c_dim1], ldc,
			     &b[l1 + min(i__3,*n) * b_dim1], ldb);
#line 866 "dtrsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
#line 867 "dtrsyl.f"
		    scaloc = 1.;

#line 869 "dtrsyl.f"
		    a11 = a[k1 + k1 * a_dim1] + sgn * b[l1 + l1 * b_dim1];
#line 870 "dtrsyl.f"
		    da11 = abs(a11);
#line 871 "dtrsyl.f"
		    if (da11 <= smin) {
#line 872 "dtrsyl.f"
			a11 = smin;
#line 873 "dtrsyl.f"
			da11 = smin;
#line 874 "dtrsyl.f"
			*info = 1;
#line 875 "dtrsyl.f"
		    }
#line 876 "dtrsyl.f"
		    db = abs(vec[0]);
#line 877 "dtrsyl.f"
		    if (da11 < 1. && db > 1.) {
#line 878 "dtrsyl.f"
			if (db > bignum * da11) {
#line 878 "dtrsyl.f"
			    scaloc = 1. / db;
#line 878 "dtrsyl.f"
			}
#line 880 "dtrsyl.f"
		    }
#line 881 "dtrsyl.f"
		    x[0] = vec[0] * scaloc / a11;

#line 883 "dtrsyl.f"
		    if (scaloc != 1.) {
#line 884 "dtrsyl.f"
			i__1 = *n;
#line 884 "dtrsyl.f"
			for (j = 1; j <= i__1; ++j) {
#line 885 "dtrsyl.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 886 "dtrsyl.f"
/* L190: */
#line 886 "dtrsyl.f"
			}
#line 887 "dtrsyl.f"
			*scale *= scaloc;
#line 888 "dtrsyl.f"
		    }
#line 889 "dtrsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];

#line 891 "dtrsyl.f"
		} else if (l1 == l2 && k1 != k2) {

#line 893 "dtrsyl.f"
		    i__1 = *m - k2;
/* Computing MIN */
#line 893 "dtrsyl.f"
		    i__2 = k2 + 1;
/* Computing MIN */
#line 893 "dtrsyl.f"
		    i__3 = k2 + 1;
#line 893 "dtrsyl.f"
		    suml = ddot_(&i__1, &a[k1 + min(i__2,*m) * a_dim1], lda, &
			    c__[min(i__3,*m) + l1 * c_dim1], &c__1);
#line 895 "dtrsyl.f"
		    i__1 = *n - l2;
/* Computing MIN */
#line 895 "dtrsyl.f"
		    i__2 = l2 + 1;
/* Computing MIN */
#line 895 "dtrsyl.f"
		    i__3 = l2 + 1;
#line 895 "dtrsyl.f"
		    sumr = ddot_(&i__1, &c__[k1 + min(i__2,*n) * c_dim1], ldc,
			     &b[l1 + min(i__3,*n) * b_dim1], ldb);
#line 897 "dtrsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

#line 899 "dtrsyl.f"
		    i__1 = *m - k2;
/* Computing MIN */
#line 899 "dtrsyl.f"
		    i__2 = k2 + 1;
/* Computing MIN */
#line 899 "dtrsyl.f"
		    i__3 = k2 + 1;
#line 899 "dtrsyl.f"
		    suml = ddot_(&i__1, &a[k2 + min(i__2,*m) * a_dim1], lda, &
			    c__[min(i__3,*m) + l1 * c_dim1], &c__1);
#line 901 "dtrsyl.f"
		    i__1 = *n - l2;
/* Computing MIN */
#line 901 "dtrsyl.f"
		    i__2 = l2 + 1;
/* Computing MIN */
#line 901 "dtrsyl.f"
		    i__3 = l2 + 1;
#line 901 "dtrsyl.f"
		    sumr = ddot_(&i__1, &c__[k2 + min(i__2,*n) * c_dim1], ldc,
			     &b[l1 + min(i__3,*n) * b_dim1], ldb);
#line 903 "dtrsyl.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

#line 905 "dtrsyl.f"
		    d__1 = -sgn * b[l1 + l1 * b_dim1];
#line 905 "dtrsyl.f"
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &a[k1 + k1 
			    * a_dim1], lda, &c_b26, &c_b26, vec, &c__2, &d__1,
			     &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
#line 908 "dtrsyl.f"
		    if (ierr != 0) {
#line 908 "dtrsyl.f"
			*info = 1;
#line 908 "dtrsyl.f"
		    }

#line 911 "dtrsyl.f"
		    if (scaloc != 1.) {
#line 912 "dtrsyl.f"
			i__1 = *n;
#line 912 "dtrsyl.f"
			for (j = 1; j <= i__1; ++j) {
#line 913 "dtrsyl.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 914 "dtrsyl.f"
/* L200: */
#line 914 "dtrsyl.f"
			}
#line 915 "dtrsyl.f"
			*scale *= scaloc;
#line 916 "dtrsyl.f"
		    }
#line 917 "dtrsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 918 "dtrsyl.f"
		    c__[k2 + l1 * c_dim1] = x[1];

#line 920 "dtrsyl.f"
		} else if (l1 != l2 && k1 == k2) {

#line 922 "dtrsyl.f"
		    i__1 = *m - k1;
/* Computing MIN */
#line 922 "dtrsyl.f"
		    i__2 = k1 + 1;
/* Computing MIN */
#line 922 "dtrsyl.f"
		    i__3 = k1 + 1;
#line 922 "dtrsyl.f"
		    suml = ddot_(&i__1, &a[k1 + min(i__2,*m) * a_dim1], lda, &
			    c__[min(i__3,*m) + l1 * c_dim1], &c__1);
#line 924 "dtrsyl.f"
		    i__1 = *n - l2;
/* Computing MIN */
#line 924 "dtrsyl.f"
		    i__2 = l2 + 1;
/* Computing MIN */
#line 924 "dtrsyl.f"
		    i__3 = l2 + 1;
#line 924 "dtrsyl.f"
		    sumr = ddot_(&i__1, &c__[k1 + min(i__2,*n) * c_dim1], ldc,
			     &b[l1 + min(i__3,*n) * b_dim1], ldb);
#line 926 "dtrsyl.f"
		    vec[0] = sgn * (c__[k1 + l1 * c_dim1] - (suml + sgn * 
			    sumr));

#line 928 "dtrsyl.f"
		    i__1 = *m - k1;
/* Computing MIN */
#line 928 "dtrsyl.f"
		    i__2 = k1 + 1;
/* Computing MIN */
#line 928 "dtrsyl.f"
		    i__3 = k1 + 1;
#line 928 "dtrsyl.f"
		    suml = ddot_(&i__1, &a[k1 + min(i__2,*m) * a_dim1], lda, &
			    c__[min(i__3,*m) + l2 * c_dim1], &c__1);
#line 930 "dtrsyl.f"
		    i__1 = *n - l2;
/* Computing MIN */
#line 930 "dtrsyl.f"
		    i__2 = l2 + 1;
/* Computing MIN */
#line 930 "dtrsyl.f"
		    i__3 = l2 + 1;
#line 930 "dtrsyl.f"
		    sumr = ddot_(&i__1, &c__[k1 + min(i__2,*n) * c_dim1], ldc,
			     &b[l2 + min(i__3,*n) * b_dim1], ldb);
#line 932 "dtrsyl.f"
		    vec[1] = sgn * (c__[k1 + l2 * c_dim1] - (suml + sgn * 
			    sumr));

#line 934 "dtrsyl.f"
		    d__1 = -sgn * a[k1 + k1 * a_dim1];
#line 934 "dtrsyl.f"
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &b[l1 + l1 
			    * b_dim1], ldb, &c_b26, &c_b26, vec, &c__2, &d__1,
			     &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
#line 937 "dtrsyl.f"
		    if (ierr != 0) {
#line 937 "dtrsyl.f"
			*info = 1;
#line 937 "dtrsyl.f"
		    }

#line 940 "dtrsyl.f"
		    if (scaloc != 1.) {
#line 941 "dtrsyl.f"
			i__1 = *n;
#line 941 "dtrsyl.f"
			for (j = 1; j <= i__1; ++j) {
#line 942 "dtrsyl.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 943 "dtrsyl.f"
/* L210: */
#line 943 "dtrsyl.f"
			}
#line 944 "dtrsyl.f"
			*scale *= scaloc;
#line 945 "dtrsyl.f"
		    }
#line 946 "dtrsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 947 "dtrsyl.f"
		    c__[k1 + l2 * c_dim1] = x[1];

#line 949 "dtrsyl.f"
		} else if (l1 != l2 && k1 != k2) {

#line 951 "dtrsyl.f"
		    i__1 = *m - k2;
/* Computing MIN */
#line 951 "dtrsyl.f"
		    i__2 = k2 + 1;
/* Computing MIN */
#line 951 "dtrsyl.f"
		    i__3 = k2 + 1;
#line 951 "dtrsyl.f"
		    suml = ddot_(&i__1, &a[k1 + min(i__2,*m) * a_dim1], lda, &
			    c__[min(i__3,*m) + l1 * c_dim1], &c__1);
#line 953 "dtrsyl.f"
		    i__1 = *n - l2;
/* Computing MIN */
#line 953 "dtrsyl.f"
		    i__2 = l2 + 1;
/* Computing MIN */
#line 953 "dtrsyl.f"
		    i__3 = l2 + 1;
#line 953 "dtrsyl.f"
		    sumr = ddot_(&i__1, &c__[k1 + min(i__2,*n) * c_dim1], ldc,
			     &b[l1 + min(i__3,*n) * b_dim1], ldb);
#line 955 "dtrsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

#line 957 "dtrsyl.f"
		    i__1 = *m - k2;
/* Computing MIN */
#line 957 "dtrsyl.f"
		    i__2 = k2 + 1;
/* Computing MIN */
#line 957 "dtrsyl.f"
		    i__3 = k2 + 1;
#line 957 "dtrsyl.f"
		    suml = ddot_(&i__1, &a[k1 + min(i__2,*m) * a_dim1], lda, &
			    c__[min(i__3,*m) + l2 * c_dim1], &c__1);
#line 959 "dtrsyl.f"
		    i__1 = *n - l2;
/* Computing MIN */
#line 959 "dtrsyl.f"
		    i__2 = l2 + 1;
/* Computing MIN */
#line 959 "dtrsyl.f"
		    i__3 = l2 + 1;
#line 959 "dtrsyl.f"
		    sumr = ddot_(&i__1, &c__[k1 + min(i__2,*n) * c_dim1], ldc,
			     &b[l2 + min(i__3,*n) * b_dim1], ldb);
#line 961 "dtrsyl.f"
		    vec[2] = c__[k1 + l2 * c_dim1] - (suml + sgn * sumr);

#line 963 "dtrsyl.f"
		    i__1 = *m - k2;
/* Computing MIN */
#line 963 "dtrsyl.f"
		    i__2 = k2 + 1;
/* Computing MIN */
#line 963 "dtrsyl.f"
		    i__3 = k2 + 1;
#line 963 "dtrsyl.f"
		    suml = ddot_(&i__1, &a[k2 + min(i__2,*m) * a_dim1], lda, &
			    c__[min(i__3,*m) + l1 * c_dim1], &c__1);
#line 965 "dtrsyl.f"
		    i__1 = *n - l2;
/* Computing MIN */
#line 965 "dtrsyl.f"
		    i__2 = l2 + 1;
/* Computing MIN */
#line 965 "dtrsyl.f"
		    i__3 = l2 + 1;
#line 965 "dtrsyl.f"
		    sumr = ddot_(&i__1, &c__[k2 + min(i__2,*n) * c_dim1], ldc,
			     &b[l1 + min(i__3,*n) * b_dim1], ldb);
#line 967 "dtrsyl.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

#line 969 "dtrsyl.f"
		    i__1 = *m - k2;
/* Computing MIN */
#line 969 "dtrsyl.f"
		    i__2 = k2 + 1;
/* Computing MIN */
#line 969 "dtrsyl.f"
		    i__3 = k2 + 1;
#line 969 "dtrsyl.f"
		    suml = ddot_(&i__1, &a[k2 + min(i__2,*m) * a_dim1], lda, &
			    c__[min(i__3,*m) + l2 * c_dim1], &c__1);
#line 971 "dtrsyl.f"
		    i__1 = *n - l2;
/* Computing MIN */
#line 971 "dtrsyl.f"
		    i__2 = l2 + 1;
/* Computing MIN */
#line 971 "dtrsyl.f"
		    i__3 = l2 + 1;
#line 971 "dtrsyl.f"
		    sumr = ddot_(&i__1, &c__[k2 + min(i__2,*n) * c_dim1], ldc,
			     &b[l2 + min(i__3,*n) * b_dim1], ldb);
#line 973 "dtrsyl.f"
		    vec[3] = c__[k2 + l2 * c_dim1] - (suml + sgn * sumr);

#line 975 "dtrsyl.f"
		    dlasy2_(&c_false, &c_true, isgn, &c__2, &c__2, &a[k1 + k1 
			    * a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec, &
			    c__2, &scaloc, x, &c__2, &xnorm, &ierr);
#line 978 "dtrsyl.f"
		    if (ierr != 0) {
#line 978 "dtrsyl.f"
			*info = 1;
#line 978 "dtrsyl.f"
		    }

#line 981 "dtrsyl.f"
		    if (scaloc != 1.) {
#line 982 "dtrsyl.f"
			i__1 = *n;
#line 982 "dtrsyl.f"
			for (j = 1; j <= i__1; ++j) {
#line 983 "dtrsyl.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 984 "dtrsyl.f"
/* L220: */
#line 984 "dtrsyl.f"
			}
#line 985 "dtrsyl.f"
			*scale *= scaloc;
#line 986 "dtrsyl.f"
		    }
#line 987 "dtrsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 988 "dtrsyl.f"
		    c__[k1 + l2 * c_dim1] = x[2];
#line 989 "dtrsyl.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 990 "dtrsyl.f"
		    c__[k2 + l2 * c_dim1] = x[3];
#line 991 "dtrsyl.f"
		}

#line 993 "dtrsyl.f"
L230:
#line 993 "dtrsyl.f"
		;
#line 993 "dtrsyl.f"
	    }
#line 994 "dtrsyl.f"
L240:
#line 994 "dtrsyl.f"
	    ;
#line 994 "dtrsyl.f"
	}

#line 996 "dtrsyl.f"
    }

#line 998 "dtrsyl.f"
    return 0;

/*     End of DTRSYL */

} /* dtrsyl_ */

