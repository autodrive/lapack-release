#line 1 "strsyl.f"
/* strsyl.f -- translated by f2c (version 20100827).
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

#line 1 "strsyl.f"
/* Table of constant values */

static integer c__1 = 1;
static logical c_false = FALSE_;
static integer c__2 = 2;
static doublereal c_b26 = 1.;
static doublereal c_b30 = 0.;
static logical c_true = TRUE_;

/* > \brief \b STRSYL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STRSYL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strsyl.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strsyl.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strsyl.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, */
/*                          LDC, SCALE, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANA, TRANB */
/*       INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N */
/*       REAL               SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), B( LDB, * ), C( LDC, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STRSYL solves the real Sylvester matrix equation: */
/* > */
/* >    op(A)*X + X*op(B) = scale*C or */
/* >    op(A)*X - X*op(B) = scale*C, */
/* > */
/* > where op(A) = A or A**T, and  A and B are both upper quasi- */
/* > triangular. A is M-by-M and B is N-by-N; the right hand side C and */
/* > the solution X are M-by-N; and scale is an output scale factor, set */
/* > <= 1 to avoid overflow in X. */
/* > */
/* > A and B must be in Schur canonical form (as returned by SHSEQR), that */
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
/* >          A is REAL array, dimension (LDA,M) */
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
/* >          B is REAL array, dimension (LDB,N) */
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
/* >          C is REAL array, dimension (LDC,N) */
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
/* >          SCALE is REAL */
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

/* > \ingroup realSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int strsyl_(char *trana, char *tranb, integer *isgn, integer 
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
    static integer ierr;
    static doublereal smin;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal suml, sumr;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer knext, lnext;
    static doublereal xnorm;
    extern /* Subroutine */ int slaln2_(logical *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, doublereal *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *),
	     slasy2_(logical *, logical *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *), slabad_(doublereal *, doublereal *);
    static doublereal scaloc;
    extern doublereal slamch_(char *, ftnlen), slange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
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

#line 211 "strsyl.f"
    /* Parameter adjustments */
#line 211 "strsyl.f"
    a_dim1 = *lda;
#line 211 "strsyl.f"
    a_offset = 1 + a_dim1;
#line 211 "strsyl.f"
    a -= a_offset;
#line 211 "strsyl.f"
    b_dim1 = *ldb;
#line 211 "strsyl.f"
    b_offset = 1 + b_dim1;
#line 211 "strsyl.f"
    b -= b_offset;
#line 211 "strsyl.f"
    c_dim1 = *ldc;
#line 211 "strsyl.f"
    c_offset = 1 + c_dim1;
#line 211 "strsyl.f"
    c__ -= c_offset;
#line 211 "strsyl.f"

#line 211 "strsyl.f"
    /* Function Body */
#line 211 "strsyl.f"
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
#line 212 "strsyl.f"
    notrnb = lsame_(tranb, "N", (ftnlen)1, (ftnlen)1);

#line 214 "strsyl.f"
    *info = 0;
#line 215 "strsyl.f"
    if (! notrna && ! lsame_(trana, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trana, "C", (ftnlen)1, (ftnlen)1)) {
#line 217 "strsyl.f"
	*info = -1;
#line 218 "strsyl.f"
    } else if (! notrnb && ! lsame_(tranb, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(tranb, "C", (ftnlen)1, (ftnlen)1)) {
#line 220 "strsyl.f"
	*info = -2;
#line 221 "strsyl.f"
    } else if (*isgn != 1 && *isgn != -1) {
#line 222 "strsyl.f"
	*info = -3;
#line 223 "strsyl.f"
    } else if (*m < 0) {
#line 224 "strsyl.f"
	*info = -4;
#line 225 "strsyl.f"
    } else if (*n < 0) {
#line 226 "strsyl.f"
	*info = -5;
#line 227 "strsyl.f"
    } else if (*lda < max(1,*m)) {
#line 228 "strsyl.f"
	*info = -7;
#line 229 "strsyl.f"
    } else if (*ldb < max(1,*n)) {
#line 230 "strsyl.f"
	*info = -9;
#line 231 "strsyl.f"
    } else if (*ldc < max(1,*m)) {
#line 232 "strsyl.f"
	*info = -11;
#line 233 "strsyl.f"
    }
#line 234 "strsyl.f"
    if (*info != 0) {
#line 235 "strsyl.f"
	i__1 = -(*info);
#line 235 "strsyl.f"
	xerbla_("STRSYL", &i__1, (ftnlen)6);
#line 236 "strsyl.f"
	return 0;
#line 237 "strsyl.f"
    }

/*     Quick return if possible */

#line 241 "strsyl.f"
    *scale = 1.;
#line 242 "strsyl.f"
    if (*m == 0 || *n == 0) {
#line 242 "strsyl.f"
	return 0;
#line 242 "strsyl.f"
    }

/*     Set constants to control overflow */

#line 247 "strsyl.f"
    eps = slamch_("P", (ftnlen)1);
#line 248 "strsyl.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 249 "strsyl.f"
    bignum = 1. / smlnum;
#line 250 "strsyl.f"
    slabad_(&smlnum, &bignum);
#line 251 "strsyl.f"
    smlnum = smlnum * (doublereal) (*m * *n) / eps;
#line 252 "strsyl.f"
    bignum = 1. / smlnum;

/* Computing MAX */
#line 254 "strsyl.f"
    d__1 = smlnum, d__2 = eps * slange_("M", m, m, &a[a_offset], lda, dum, (
	    ftnlen)1), d__1 = max(d__1,d__2), d__2 = eps * slange_("M", n, n, 
	    &b[b_offset], ldb, dum, (ftnlen)1);
#line 254 "strsyl.f"
    smin = max(d__1,d__2);

#line 257 "strsyl.f"
    sgn = (doublereal) (*isgn);

#line 259 "strsyl.f"
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

#line 276 "strsyl.f"
	lnext = 1;
#line 277 "strsyl.f"
	i__1 = *n;
#line 277 "strsyl.f"
	for (l = 1; l <= i__1; ++l) {
#line 278 "strsyl.f"
	    if (l < lnext) {
#line 278 "strsyl.f"
		goto L70;
#line 278 "strsyl.f"
	    }
#line 280 "strsyl.f"
	    if (l == *n) {
#line 281 "strsyl.f"
		l1 = l;
#line 282 "strsyl.f"
		l2 = l;
#line 283 "strsyl.f"
	    } else {
#line 284 "strsyl.f"
		if (b[l + 1 + l * b_dim1] != 0.) {
#line 285 "strsyl.f"
		    l1 = l;
#line 286 "strsyl.f"
		    l2 = l + 1;
#line 287 "strsyl.f"
		    lnext = l + 2;
#line 288 "strsyl.f"
		} else {
#line 289 "strsyl.f"
		    l1 = l;
#line 290 "strsyl.f"
		    l2 = l;
#line 291 "strsyl.f"
		    lnext = l + 1;
#line 292 "strsyl.f"
		}
#line 293 "strsyl.f"
	    }

/*           Start row loop (index = K) */
/*           K1 (K2): row index of the first (last) row of X(K,L). */

#line 298 "strsyl.f"
	    knext = *m;
#line 299 "strsyl.f"
	    for (k = *m; k >= 1; --k) {
#line 300 "strsyl.f"
		if (k > knext) {
#line 300 "strsyl.f"
		    goto L60;
#line 300 "strsyl.f"
		}
#line 302 "strsyl.f"
		if (k == 1) {
#line 303 "strsyl.f"
		    k1 = k;
#line 304 "strsyl.f"
		    k2 = k;
#line 305 "strsyl.f"
		} else {
#line 306 "strsyl.f"
		    if (a[k + (k - 1) * a_dim1] != 0.) {
#line 307 "strsyl.f"
			k1 = k - 1;
#line 308 "strsyl.f"
			k2 = k;
#line 309 "strsyl.f"
			knext = k - 2;
#line 310 "strsyl.f"
		    } else {
#line 311 "strsyl.f"
			k1 = k;
#line 312 "strsyl.f"
			k2 = k;
#line 313 "strsyl.f"
			knext = k - 1;
#line 314 "strsyl.f"
		    }
#line 315 "strsyl.f"
		}

#line 317 "strsyl.f"
		if (l1 == l2 && k1 == k2) {
#line 318 "strsyl.f"
		    i__2 = *m - k1;
/* Computing MIN */
#line 318 "strsyl.f"
		    i__3 = k1 + 1;
/* Computing MIN */
#line 318 "strsyl.f"
		    i__4 = k1 + 1;
#line 318 "strsyl.f"
		    suml = sdot_(&i__2, &a[k1 + min(i__3,*m) * a_dim1], lda, &
			    c__[min(i__4,*m) + l1 * c_dim1], &c__1);
#line 320 "strsyl.f"
		    i__2 = l1 - 1;
#line 320 "strsyl.f"
		    sumr = sdot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 321 "strsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
#line 322 "strsyl.f"
		    scaloc = 1.;

#line 324 "strsyl.f"
		    a11 = a[k1 + k1 * a_dim1] + sgn * b[l1 + l1 * b_dim1];
#line 325 "strsyl.f"
		    da11 = abs(a11);
#line 326 "strsyl.f"
		    if (da11 <= smin) {
#line 327 "strsyl.f"
			a11 = smin;
#line 328 "strsyl.f"
			da11 = smin;
#line 329 "strsyl.f"
			*info = 1;
#line 330 "strsyl.f"
		    }
#line 331 "strsyl.f"
		    db = abs(vec[0]);
#line 332 "strsyl.f"
		    if (da11 < 1. && db > 1.) {
#line 333 "strsyl.f"
			if (db > bignum * da11) {
#line 333 "strsyl.f"
			    scaloc = 1. / db;
#line 333 "strsyl.f"
			}
#line 335 "strsyl.f"
		    }
#line 336 "strsyl.f"
		    x[0] = vec[0] * scaloc / a11;

#line 338 "strsyl.f"
		    if (scaloc != 1.) {
#line 339 "strsyl.f"
			i__2 = *n;
#line 339 "strsyl.f"
			for (j = 1; j <= i__2; ++j) {
#line 340 "strsyl.f"
			    sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 341 "strsyl.f"
/* L10: */
#line 341 "strsyl.f"
			}
#line 342 "strsyl.f"
			*scale *= scaloc;
#line 343 "strsyl.f"
		    }
#line 344 "strsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];

#line 346 "strsyl.f"
		} else if (l1 == l2 && k1 != k2) {

#line 348 "strsyl.f"
		    i__2 = *m - k2;
/* Computing MIN */
#line 348 "strsyl.f"
		    i__3 = k2 + 1;
/* Computing MIN */
#line 348 "strsyl.f"
		    i__4 = k2 + 1;
#line 348 "strsyl.f"
		    suml = sdot_(&i__2, &a[k1 + min(i__3,*m) * a_dim1], lda, &
			    c__[min(i__4,*m) + l1 * c_dim1], &c__1);
#line 350 "strsyl.f"
		    i__2 = l1 - 1;
#line 350 "strsyl.f"
		    sumr = sdot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 351 "strsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

#line 353 "strsyl.f"
		    i__2 = *m - k2;
/* Computing MIN */
#line 353 "strsyl.f"
		    i__3 = k2 + 1;
/* Computing MIN */
#line 353 "strsyl.f"
		    i__4 = k2 + 1;
#line 353 "strsyl.f"
		    suml = sdot_(&i__2, &a[k2 + min(i__3,*m) * a_dim1], lda, &
			    c__[min(i__4,*m) + l1 * c_dim1], &c__1);
#line 355 "strsyl.f"
		    i__2 = l1 - 1;
#line 355 "strsyl.f"
		    sumr = sdot_(&i__2, &c__[k2 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 356 "strsyl.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

#line 358 "strsyl.f"
		    d__1 = -sgn * b[l1 + l1 * b_dim1];
#line 358 "strsyl.f"
		    slaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &a[k1 + k1 
			    * a_dim1], lda, &c_b26, &c_b26, vec, &c__2, &d__1,
			     &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
#line 361 "strsyl.f"
		    if (ierr != 0) {
#line 361 "strsyl.f"
			*info = 1;
#line 361 "strsyl.f"
		    }

#line 364 "strsyl.f"
		    if (scaloc != 1.) {
#line 365 "strsyl.f"
			i__2 = *n;
#line 365 "strsyl.f"
			for (j = 1; j <= i__2; ++j) {
#line 366 "strsyl.f"
			    sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 367 "strsyl.f"
/* L20: */
#line 367 "strsyl.f"
			}
#line 368 "strsyl.f"
			*scale *= scaloc;
#line 369 "strsyl.f"
		    }
#line 370 "strsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 371 "strsyl.f"
		    c__[k2 + l1 * c_dim1] = x[1];

#line 373 "strsyl.f"
		} else if (l1 != l2 && k1 == k2) {

#line 375 "strsyl.f"
		    i__2 = *m - k1;
/* Computing MIN */
#line 375 "strsyl.f"
		    i__3 = k1 + 1;
/* Computing MIN */
#line 375 "strsyl.f"
		    i__4 = k1 + 1;
#line 375 "strsyl.f"
		    suml = sdot_(&i__2, &a[k1 + min(i__3,*m) * a_dim1], lda, &
			    c__[min(i__4,*m) + l1 * c_dim1], &c__1);
#line 377 "strsyl.f"
		    i__2 = l1 - 1;
#line 377 "strsyl.f"
		    sumr = sdot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 378 "strsyl.f"
		    vec[0] = sgn * (c__[k1 + l1 * c_dim1] - (suml + sgn * 
			    sumr));

#line 380 "strsyl.f"
		    i__2 = *m - k1;
/* Computing MIN */
#line 380 "strsyl.f"
		    i__3 = k1 + 1;
/* Computing MIN */
#line 380 "strsyl.f"
		    i__4 = k1 + 1;
#line 380 "strsyl.f"
		    suml = sdot_(&i__2, &a[k1 + min(i__3,*m) * a_dim1], lda, &
			    c__[min(i__4,*m) + l2 * c_dim1], &c__1);
#line 382 "strsyl.f"
		    i__2 = l1 - 1;
#line 382 "strsyl.f"
		    sumr = sdot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l2 * 
			    b_dim1 + 1], &c__1);
#line 383 "strsyl.f"
		    vec[1] = sgn * (c__[k1 + l2 * c_dim1] - (suml + sgn * 
			    sumr));

#line 385 "strsyl.f"
		    d__1 = -sgn * a[k1 + k1 * a_dim1];
#line 385 "strsyl.f"
		    slaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &b[l1 + l1 *
			     b_dim1], ldb, &c_b26, &c_b26, vec, &c__2, &d__1, 
			    &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
#line 388 "strsyl.f"
		    if (ierr != 0) {
#line 388 "strsyl.f"
			*info = 1;
#line 388 "strsyl.f"
		    }

#line 391 "strsyl.f"
		    if (scaloc != 1.) {
#line 392 "strsyl.f"
			i__2 = *n;
#line 392 "strsyl.f"
			for (j = 1; j <= i__2; ++j) {
#line 393 "strsyl.f"
			    sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 394 "strsyl.f"
/* L40: */
#line 394 "strsyl.f"
			}
#line 395 "strsyl.f"
			*scale *= scaloc;
#line 396 "strsyl.f"
		    }
#line 397 "strsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 398 "strsyl.f"
		    c__[k1 + l2 * c_dim1] = x[1];

#line 400 "strsyl.f"
		} else if (l1 != l2 && k1 != k2) {

#line 402 "strsyl.f"
		    i__2 = *m - k2;
/* Computing MIN */
#line 402 "strsyl.f"
		    i__3 = k2 + 1;
/* Computing MIN */
#line 402 "strsyl.f"
		    i__4 = k2 + 1;
#line 402 "strsyl.f"
		    suml = sdot_(&i__2, &a[k1 + min(i__3,*m) * a_dim1], lda, &
			    c__[min(i__4,*m) + l1 * c_dim1], &c__1);
#line 404 "strsyl.f"
		    i__2 = l1 - 1;
#line 404 "strsyl.f"
		    sumr = sdot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 405 "strsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

#line 407 "strsyl.f"
		    i__2 = *m - k2;
/* Computing MIN */
#line 407 "strsyl.f"
		    i__3 = k2 + 1;
/* Computing MIN */
#line 407 "strsyl.f"
		    i__4 = k2 + 1;
#line 407 "strsyl.f"
		    suml = sdot_(&i__2, &a[k1 + min(i__3,*m) * a_dim1], lda, &
			    c__[min(i__4,*m) + l2 * c_dim1], &c__1);
#line 409 "strsyl.f"
		    i__2 = l1 - 1;
#line 409 "strsyl.f"
		    sumr = sdot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l2 * 
			    b_dim1 + 1], &c__1);
#line 410 "strsyl.f"
		    vec[2] = c__[k1 + l2 * c_dim1] - (suml + sgn * sumr);

#line 412 "strsyl.f"
		    i__2 = *m - k2;
/* Computing MIN */
#line 412 "strsyl.f"
		    i__3 = k2 + 1;
/* Computing MIN */
#line 412 "strsyl.f"
		    i__4 = k2 + 1;
#line 412 "strsyl.f"
		    suml = sdot_(&i__2, &a[k2 + min(i__3,*m) * a_dim1], lda, &
			    c__[min(i__4,*m) + l1 * c_dim1], &c__1);
#line 414 "strsyl.f"
		    i__2 = l1 - 1;
#line 414 "strsyl.f"
		    sumr = sdot_(&i__2, &c__[k2 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 415 "strsyl.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

#line 417 "strsyl.f"
		    i__2 = *m - k2;
/* Computing MIN */
#line 417 "strsyl.f"
		    i__3 = k2 + 1;
/* Computing MIN */
#line 417 "strsyl.f"
		    i__4 = k2 + 1;
#line 417 "strsyl.f"
		    suml = sdot_(&i__2, &a[k2 + min(i__3,*m) * a_dim1], lda, &
			    c__[min(i__4,*m) + l2 * c_dim1], &c__1);
#line 419 "strsyl.f"
		    i__2 = l1 - 1;
#line 419 "strsyl.f"
		    sumr = sdot_(&i__2, &c__[k2 + c_dim1], ldc, &b[l2 * 
			    b_dim1 + 1], &c__1);
#line 420 "strsyl.f"
		    vec[3] = c__[k2 + l2 * c_dim1] - (suml + sgn * sumr);

#line 422 "strsyl.f"
		    slasy2_(&c_false, &c_false, isgn, &c__2, &c__2, &a[k1 + 
			    k1 * a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec,
			     &c__2, &scaloc, x, &c__2, &xnorm, &ierr);
#line 425 "strsyl.f"
		    if (ierr != 0) {
#line 425 "strsyl.f"
			*info = 1;
#line 425 "strsyl.f"
		    }

#line 428 "strsyl.f"
		    if (scaloc != 1.) {
#line 429 "strsyl.f"
			i__2 = *n;
#line 429 "strsyl.f"
			for (j = 1; j <= i__2; ++j) {
#line 430 "strsyl.f"
			    sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 431 "strsyl.f"
/* L50: */
#line 431 "strsyl.f"
			}
#line 432 "strsyl.f"
			*scale *= scaloc;
#line 433 "strsyl.f"
		    }
#line 434 "strsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 435 "strsyl.f"
		    c__[k1 + l2 * c_dim1] = x[2];
#line 436 "strsyl.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 437 "strsyl.f"
		    c__[k2 + l2 * c_dim1] = x[3];
#line 438 "strsyl.f"
		}

#line 440 "strsyl.f"
L60:
#line 440 "strsyl.f"
		;
#line 440 "strsyl.f"
	    }

#line 442 "strsyl.f"
L70:
#line 442 "strsyl.f"
	    ;
#line 442 "strsyl.f"
	}

#line 444 "strsyl.f"
    } else if (! notrna && notrnb) {

/*        Solve    A**T *X + ISGN*X*B = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        upper-left corner column by column by */

/*          A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L) */

/*        Where */
/*                   K-1                          L-1 */
/*          R(K,L) = SUM [A(I,K)**T*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)] */
/*                   I=1                          J=1 */

/*        Start column loop (index = L) */
/*        L1 (L2): column index of the first (last) row of X(K,L) */

#line 461 "strsyl.f"
	lnext = 1;
#line 462 "strsyl.f"
	i__1 = *n;
#line 462 "strsyl.f"
	for (l = 1; l <= i__1; ++l) {
#line 463 "strsyl.f"
	    if (l < lnext) {
#line 463 "strsyl.f"
		goto L130;
#line 463 "strsyl.f"
	    }
#line 465 "strsyl.f"
	    if (l == *n) {
#line 466 "strsyl.f"
		l1 = l;
#line 467 "strsyl.f"
		l2 = l;
#line 468 "strsyl.f"
	    } else {
#line 469 "strsyl.f"
		if (b[l + 1 + l * b_dim1] != 0.) {
#line 470 "strsyl.f"
		    l1 = l;
#line 471 "strsyl.f"
		    l2 = l + 1;
#line 472 "strsyl.f"
		    lnext = l + 2;
#line 473 "strsyl.f"
		} else {
#line 474 "strsyl.f"
		    l1 = l;
#line 475 "strsyl.f"
		    l2 = l;
#line 476 "strsyl.f"
		    lnext = l + 1;
#line 477 "strsyl.f"
		}
#line 478 "strsyl.f"
	    }

/*           Start row loop (index = K) */
/*           K1 (K2): row index of the first (last) row of X(K,L) */

#line 483 "strsyl.f"
	    knext = 1;
#line 484 "strsyl.f"
	    i__2 = *m;
#line 484 "strsyl.f"
	    for (k = 1; k <= i__2; ++k) {
#line 485 "strsyl.f"
		if (k < knext) {
#line 485 "strsyl.f"
		    goto L120;
#line 485 "strsyl.f"
		}
#line 487 "strsyl.f"
		if (k == *m) {
#line 488 "strsyl.f"
		    k1 = k;
#line 489 "strsyl.f"
		    k2 = k;
#line 490 "strsyl.f"
		} else {
#line 491 "strsyl.f"
		    if (a[k + 1 + k * a_dim1] != 0.) {
#line 492 "strsyl.f"
			k1 = k;
#line 493 "strsyl.f"
			k2 = k + 1;
#line 494 "strsyl.f"
			knext = k + 2;
#line 495 "strsyl.f"
		    } else {
#line 496 "strsyl.f"
			k1 = k;
#line 497 "strsyl.f"
			k2 = k;
#line 498 "strsyl.f"
			knext = k + 1;
#line 499 "strsyl.f"
		    }
#line 500 "strsyl.f"
		}

#line 502 "strsyl.f"
		if (l1 == l2 && k1 == k2) {
#line 503 "strsyl.f"
		    i__3 = k1 - 1;
#line 503 "strsyl.f"
		    suml = sdot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 504 "strsyl.f"
		    i__3 = l1 - 1;
#line 504 "strsyl.f"
		    sumr = sdot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 505 "strsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
#line 506 "strsyl.f"
		    scaloc = 1.;

#line 508 "strsyl.f"
		    a11 = a[k1 + k1 * a_dim1] + sgn * b[l1 + l1 * b_dim1];
#line 509 "strsyl.f"
		    da11 = abs(a11);
#line 510 "strsyl.f"
		    if (da11 <= smin) {
#line 511 "strsyl.f"
			a11 = smin;
#line 512 "strsyl.f"
			da11 = smin;
#line 513 "strsyl.f"
			*info = 1;
#line 514 "strsyl.f"
		    }
#line 515 "strsyl.f"
		    db = abs(vec[0]);
#line 516 "strsyl.f"
		    if (da11 < 1. && db > 1.) {
#line 517 "strsyl.f"
			if (db > bignum * da11) {
#line 517 "strsyl.f"
			    scaloc = 1. / db;
#line 517 "strsyl.f"
			}
#line 519 "strsyl.f"
		    }
#line 520 "strsyl.f"
		    x[0] = vec[0] * scaloc / a11;

#line 522 "strsyl.f"
		    if (scaloc != 1.) {
#line 523 "strsyl.f"
			i__3 = *n;
#line 523 "strsyl.f"
			for (j = 1; j <= i__3; ++j) {
#line 524 "strsyl.f"
			    sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 525 "strsyl.f"
/* L80: */
#line 525 "strsyl.f"
			}
#line 526 "strsyl.f"
			*scale *= scaloc;
#line 527 "strsyl.f"
		    }
#line 528 "strsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];

#line 530 "strsyl.f"
		} else if (l1 == l2 && k1 != k2) {

#line 532 "strsyl.f"
		    i__3 = k1 - 1;
#line 532 "strsyl.f"
		    suml = sdot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 533 "strsyl.f"
		    i__3 = l1 - 1;
#line 533 "strsyl.f"
		    sumr = sdot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 534 "strsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

#line 536 "strsyl.f"
		    i__3 = k1 - 1;
#line 536 "strsyl.f"
		    suml = sdot_(&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 537 "strsyl.f"
		    i__3 = l1 - 1;
#line 537 "strsyl.f"
		    sumr = sdot_(&i__3, &c__[k2 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 538 "strsyl.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

#line 540 "strsyl.f"
		    d__1 = -sgn * b[l1 + l1 * b_dim1];
#line 540 "strsyl.f"
		    slaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &a[k1 + k1 *
			     a_dim1], lda, &c_b26, &c_b26, vec, &c__2, &d__1, 
			    &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
#line 543 "strsyl.f"
		    if (ierr != 0) {
#line 543 "strsyl.f"
			*info = 1;
#line 543 "strsyl.f"
		    }

#line 546 "strsyl.f"
		    if (scaloc != 1.) {
#line 547 "strsyl.f"
			i__3 = *n;
#line 547 "strsyl.f"
			for (j = 1; j <= i__3; ++j) {
#line 548 "strsyl.f"
			    sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 549 "strsyl.f"
/* L90: */
#line 549 "strsyl.f"
			}
#line 550 "strsyl.f"
			*scale *= scaloc;
#line 551 "strsyl.f"
		    }
#line 552 "strsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 553 "strsyl.f"
		    c__[k2 + l1 * c_dim1] = x[1];

#line 555 "strsyl.f"
		} else if (l1 != l2 && k1 == k2) {

#line 557 "strsyl.f"
		    i__3 = k1 - 1;
#line 557 "strsyl.f"
		    suml = sdot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 558 "strsyl.f"
		    i__3 = l1 - 1;
#line 558 "strsyl.f"
		    sumr = sdot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 559 "strsyl.f"
		    vec[0] = sgn * (c__[k1 + l1 * c_dim1] - (suml + sgn * 
			    sumr));

#line 561 "strsyl.f"
		    i__3 = k1 - 1;
#line 561 "strsyl.f"
		    suml = sdot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
#line 562 "strsyl.f"
		    i__3 = l1 - 1;
#line 562 "strsyl.f"
		    sumr = sdot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l2 * 
			    b_dim1 + 1], &c__1);
#line 563 "strsyl.f"
		    vec[1] = sgn * (c__[k1 + l2 * c_dim1] - (suml + sgn * 
			    sumr));

#line 565 "strsyl.f"
		    d__1 = -sgn * a[k1 + k1 * a_dim1];
#line 565 "strsyl.f"
		    slaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &b[l1 + l1 *
			     b_dim1], ldb, &c_b26, &c_b26, vec, &c__2, &d__1, 
			    &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
#line 568 "strsyl.f"
		    if (ierr != 0) {
#line 568 "strsyl.f"
			*info = 1;
#line 568 "strsyl.f"
		    }

#line 571 "strsyl.f"
		    if (scaloc != 1.) {
#line 572 "strsyl.f"
			i__3 = *n;
#line 572 "strsyl.f"
			for (j = 1; j <= i__3; ++j) {
#line 573 "strsyl.f"
			    sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 574 "strsyl.f"
/* L100: */
#line 574 "strsyl.f"
			}
#line 575 "strsyl.f"
			*scale *= scaloc;
#line 576 "strsyl.f"
		    }
#line 577 "strsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 578 "strsyl.f"
		    c__[k1 + l2 * c_dim1] = x[1];

#line 580 "strsyl.f"
		} else if (l1 != l2 && k1 != k2) {

#line 582 "strsyl.f"
		    i__3 = k1 - 1;
#line 582 "strsyl.f"
		    suml = sdot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 583 "strsyl.f"
		    i__3 = l1 - 1;
#line 583 "strsyl.f"
		    sumr = sdot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 584 "strsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

#line 586 "strsyl.f"
		    i__3 = k1 - 1;
#line 586 "strsyl.f"
		    suml = sdot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
#line 587 "strsyl.f"
		    i__3 = l1 - 1;
#line 587 "strsyl.f"
		    sumr = sdot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l2 * 
			    b_dim1 + 1], &c__1);
#line 588 "strsyl.f"
		    vec[2] = c__[k1 + l2 * c_dim1] - (suml + sgn * sumr);

#line 590 "strsyl.f"
		    i__3 = k1 - 1;
#line 590 "strsyl.f"
		    suml = sdot_(&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 591 "strsyl.f"
		    i__3 = l1 - 1;
#line 591 "strsyl.f"
		    sumr = sdot_(&i__3, &c__[k2 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 592 "strsyl.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

#line 594 "strsyl.f"
		    i__3 = k1 - 1;
#line 594 "strsyl.f"
		    suml = sdot_(&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
#line 595 "strsyl.f"
		    i__3 = l1 - 1;
#line 595 "strsyl.f"
		    sumr = sdot_(&i__3, &c__[k2 + c_dim1], ldc, &b[l2 * 
			    b_dim1 + 1], &c__1);
#line 596 "strsyl.f"
		    vec[3] = c__[k2 + l2 * c_dim1] - (suml + sgn * sumr);

#line 598 "strsyl.f"
		    slasy2_(&c_true, &c_false, isgn, &c__2, &c__2, &a[k1 + k1 
			    * a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec, &
			    c__2, &scaloc, x, &c__2, &xnorm, &ierr);
#line 601 "strsyl.f"
		    if (ierr != 0) {
#line 601 "strsyl.f"
			*info = 1;
#line 601 "strsyl.f"
		    }

#line 604 "strsyl.f"
		    if (scaloc != 1.) {
#line 605 "strsyl.f"
			i__3 = *n;
#line 605 "strsyl.f"
			for (j = 1; j <= i__3; ++j) {
#line 606 "strsyl.f"
			    sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 607 "strsyl.f"
/* L110: */
#line 607 "strsyl.f"
			}
#line 608 "strsyl.f"
			*scale *= scaloc;
#line 609 "strsyl.f"
		    }
#line 610 "strsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 611 "strsyl.f"
		    c__[k1 + l2 * c_dim1] = x[2];
#line 612 "strsyl.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 613 "strsyl.f"
		    c__[k2 + l2 * c_dim1] = x[3];
#line 614 "strsyl.f"
		}

#line 616 "strsyl.f"
L120:
#line 616 "strsyl.f"
		;
#line 616 "strsyl.f"
	    }
#line 617 "strsyl.f"
L130:
#line 617 "strsyl.f"
	    ;
#line 617 "strsyl.f"
	}

#line 619 "strsyl.f"
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

#line 636 "strsyl.f"
	lnext = *n;
#line 637 "strsyl.f"
	for (l = *n; l >= 1; --l) {
#line 638 "strsyl.f"
	    if (l > lnext) {
#line 638 "strsyl.f"
		goto L190;
#line 638 "strsyl.f"
	    }
#line 640 "strsyl.f"
	    if (l == 1) {
#line 641 "strsyl.f"
		l1 = l;
#line 642 "strsyl.f"
		l2 = l;
#line 643 "strsyl.f"
	    } else {
#line 644 "strsyl.f"
		if (b[l + (l - 1) * b_dim1] != 0.) {
#line 645 "strsyl.f"
		    l1 = l - 1;
#line 646 "strsyl.f"
		    l2 = l;
#line 647 "strsyl.f"
		    lnext = l - 2;
#line 648 "strsyl.f"
		} else {
#line 649 "strsyl.f"
		    l1 = l;
#line 650 "strsyl.f"
		    l2 = l;
#line 651 "strsyl.f"
		    lnext = l - 1;
#line 652 "strsyl.f"
		}
#line 653 "strsyl.f"
	    }

/*           Start row loop (index = K) */
/*           K1 (K2): row index of the first (last) row of X(K,L) */

#line 658 "strsyl.f"
	    knext = 1;
#line 659 "strsyl.f"
	    i__1 = *m;
#line 659 "strsyl.f"
	    for (k = 1; k <= i__1; ++k) {
#line 660 "strsyl.f"
		if (k < knext) {
#line 660 "strsyl.f"
		    goto L180;
#line 660 "strsyl.f"
		}
#line 662 "strsyl.f"
		if (k == *m) {
#line 663 "strsyl.f"
		    k1 = k;
#line 664 "strsyl.f"
		    k2 = k;
#line 665 "strsyl.f"
		} else {
#line 666 "strsyl.f"
		    if (a[k + 1 + k * a_dim1] != 0.) {
#line 667 "strsyl.f"
			k1 = k;
#line 668 "strsyl.f"
			k2 = k + 1;
#line 669 "strsyl.f"
			knext = k + 2;
#line 670 "strsyl.f"
		    } else {
#line 671 "strsyl.f"
			k1 = k;
#line 672 "strsyl.f"
			k2 = k;
#line 673 "strsyl.f"
			knext = k + 1;
#line 674 "strsyl.f"
		    }
#line 675 "strsyl.f"
		}

#line 677 "strsyl.f"
		if (l1 == l2 && k1 == k2) {
#line 678 "strsyl.f"
		    i__2 = k1 - 1;
#line 678 "strsyl.f"
		    suml = sdot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 679 "strsyl.f"
		    i__2 = *n - l1;
/* Computing MIN */
#line 679 "strsyl.f"
		    i__3 = l1 + 1;
/* Computing MIN */
#line 679 "strsyl.f"
		    i__4 = l1 + 1;
#line 679 "strsyl.f"
		    sumr = sdot_(&i__2, &c__[k1 + min(i__3,*n) * c_dim1], ldc,
			     &b[l1 + min(i__4,*n) * b_dim1], ldb);
#line 681 "strsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
#line 682 "strsyl.f"
		    scaloc = 1.;

#line 684 "strsyl.f"
		    a11 = a[k1 + k1 * a_dim1] + sgn * b[l1 + l1 * b_dim1];
#line 685 "strsyl.f"
		    da11 = abs(a11);
#line 686 "strsyl.f"
		    if (da11 <= smin) {
#line 687 "strsyl.f"
			a11 = smin;
#line 688 "strsyl.f"
			da11 = smin;
#line 689 "strsyl.f"
			*info = 1;
#line 690 "strsyl.f"
		    }
#line 691 "strsyl.f"
		    db = abs(vec[0]);
#line 692 "strsyl.f"
		    if (da11 < 1. && db > 1.) {
#line 693 "strsyl.f"
			if (db > bignum * da11) {
#line 693 "strsyl.f"
			    scaloc = 1. / db;
#line 693 "strsyl.f"
			}
#line 695 "strsyl.f"
		    }
#line 696 "strsyl.f"
		    x[0] = vec[0] * scaloc / a11;

#line 698 "strsyl.f"
		    if (scaloc != 1.) {
#line 699 "strsyl.f"
			i__2 = *n;
#line 699 "strsyl.f"
			for (j = 1; j <= i__2; ++j) {
#line 700 "strsyl.f"
			    sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 701 "strsyl.f"
/* L140: */
#line 701 "strsyl.f"
			}
#line 702 "strsyl.f"
			*scale *= scaloc;
#line 703 "strsyl.f"
		    }
#line 704 "strsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];

#line 706 "strsyl.f"
		} else if (l1 == l2 && k1 != k2) {

#line 708 "strsyl.f"
		    i__2 = k1 - 1;
#line 708 "strsyl.f"
		    suml = sdot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 709 "strsyl.f"
		    i__2 = *n - l2;
/* Computing MIN */
#line 709 "strsyl.f"
		    i__3 = l2 + 1;
/* Computing MIN */
#line 709 "strsyl.f"
		    i__4 = l2 + 1;
#line 709 "strsyl.f"
		    sumr = sdot_(&i__2, &c__[k1 + min(i__3,*n) * c_dim1], ldc,
			     &b[l1 + min(i__4,*n) * b_dim1], ldb);
#line 711 "strsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

#line 713 "strsyl.f"
		    i__2 = k1 - 1;
#line 713 "strsyl.f"
		    suml = sdot_(&i__2, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 714 "strsyl.f"
		    i__2 = *n - l2;
/* Computing MIN */
#line 714 "strsyl.f"
		    i__3 = l2 + 1;
/* Computing MIN */
#line 714 "strsyl.f"
		    i__4 = l2 + 1;
#line 714 "strsyl.f"
		    sumr = sdot_(&i__2, &c__[k2 + min(i__3,*n) * c_dim1], ldc,
			     &b[l1 + min(i__4,*n) * b_dim1], ldb);
#line 716 "strsyl.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

#line 718 "strsyl.f"
		    d__1 = -sgn * b[l1 + l1 * b_dim1];
#line 718 "strsyl.f"
		    slaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &a[k1 + k1 *
			     a_dim1], lda, &c_b26, &c_b26, vec, &c__2, &d__1, 
			    &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
#line 721 "strsyl.f"
		    if (ierr != 0) {
#line 721 "strsyl.f"
			*info = 1;
#line 721 "strsyl.f"
		    }

#line 724 "strsyl.f"
		    if (scaloc != 1.) {
#line 725 "strsyl.f"
			i__2 = *n;
#line 725 "strsyl.f"
			for (j = 1; j <= i__2; ++j) {
#line 726 "strsyl.f"
			    sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 727 "strsyl.f"
/* L150: */
#line 727 "strsyl.f"
			}
#line 728 "strsyl.f"
			*scale *= scaloc;
#line 729 "strsyl.f"
		    }
#line 730 "strsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 731 "strsyl.f"
		    c__[k2 + l1 * c_dim1] = x[1];

#line 733 "strsyl.f"
		} else if (l1 != l2 && k1 == k2) {

#line 735 "strsyl.f"
		    i__2 = k1 - 1;
#line 735 "strsyl.f"
		    suml = sdot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 736 "strsyl.f"
		    i__2 = *n - l2;
/* Computing MIN */
#line 736 "strsyl.f"
		    i__3 = l2 + 1;
/* Computing MIN */
#line 736 "strsyl.f"
		    i__4 = l2 + 1;
#line 736 "strsyl.f"
		    sumr = sdot_(&i__2, &c__[k1 + min(i__3,*n) * c_dim1], ldc,
			     &b[l1 + min(i__4,*n) * b_dim1], ldb);
#line 738 "strsyl.f"
		    vec[0] = sgn * (c__[k1 + l1 * c_dim1] - (suml + sgn * 
			    sumr));

#line 740 "strsyl.f"
		    i__2 = k1 - 1;
#line 740 "strsyl.f"
		    suml = sdot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
#line 741 "strsyl.f"
		    i__2 = *n - l2;
/* Computing MIN */
#line 741 "strsyl.f"
		    i__3 = l2 + 1;
/* Computing MIN */
#line 741 "strsyl.f"
		    i__4 = l2 + 1;
#line 741 "strsyl.f"
		    sumr = sdot_(&i__2, &c__[k1 + min(i__3,*n) * c_dim1], ldc,
			     &b[l2 + min(i__4,*n) * b_dim1], ldb);
#line 743 "strsyl.f"
		    vec[1] = sgn * (c__[k1 + l2 * c_dim1] - (suml + sgn * 
			    sumr));

#line 745 "strsyl.f"
		    d__1 = -sgn * a[k1 + k1 * a_dim1];
#line 745 "strsyl.f"
		    slaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &b[l1 + l1 
			    * b_dim1], ldb, &c_b26, &c_b26, vec, &c__2, &d__1,
			     &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
#line 748 "strsyl.f"
		    if (ierr != 0) {
#line 748 "strsyl.f"
			*info = 1;
#line 748 "strsyl.f"
		    }

#line 751 "strsyl.f"
		    if (scaloc != 1.) {
#line 752 "strsyl.f"
			i__2 = *n;
#line 752 "strsyl.f"
			for (j = 1; j <= i__2; ++j) {
#line 753 "strsyl.f"
			    sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 754 "strsyl.f"
/* L160: */
#line 754 "strsyl.f"
			}
#line 755 "strsyl.f"
			*scale *= scaloc;
#line 756 "strsyl.f"
		    }
#line 757 "strsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 758 "strsyl.f"
		    c__[k1 + l2 * c_dim1] = x[1];

#line 760 "strsyl.f"
		} else if (l1 != l2 && k1 != k2) {

#line 762 "strsyl.f"
		    i__2 = k1 - 1;
#line 762 "strsyl.f"
		    suml = sdot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 763 "strsyl.f"
		    i__2 = *n - l2;
/* Computing MIN */
#line 763 "strsyl.f"
		    i__3 = l2 + 1;
/* Computing MIN */
#line 763 "strsyl.f"
		    i__4 = l2 + 1;
#line 763 "strsyl.f"
		    sumr = sdot_(&i__2, &c__[k1 + min(i__3,*n) * c_dim1], ldc,
			     &b[l1 + min(i__4,*n) * b_dim1], ldb);
#line 765 "strsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

#line 767 "strsyl.f"
		    i__2 = k1 - 1;
#line 767 "strsyl.f"
		    suml = sdot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
#line 768 "strsyl.f"
		    i__2 = *n - l2;
/* Computing MIN */
#line 768 "strsyl.f"
		    i__3 = l2 + 1;
/* Computing MIN */
#line 768 "strsyl.f"
		    i__4 = l2 + 1;
#line 768 "strsyl.f"
		    sumr = sdot_(&i__2, &c__[k1 + min(i__3,*n) * c_dim1], ldc,
			     &b[l2 + min(i__4,*n) * b_dim1], ldb);
#line 770 "strsyl.f"
		    vec[2] = c__[k1 + l2 * c_dim1] - (suml + sgn * sumr);

#line 772 "strsyl.f"
		    i__2 = k1 - 1;
#line 772 "strsyl.f"
		    suml = sdot_(&i__2, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 773 "strsyl.f"
		    i__2 = *n - l2;
/* Computing MIN */
#line 773 "strsyl.f"
		    i__3 = l2 + 1;
/* Computing MIN */
#line 773 "strsyl.f"
		    i__4 = l2 + 1;
#line 773 "strsyl.f"
		    sumr = sdot_(&i__2, &c__[k2 + min(i__3,*n) * c_dim1], ldc,
			     &b[l1 + min(i__4,*n) * b_dim1], ldb);
#line 775 "strsyl.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

#line 777 "strsyl.f"
		    i__2 = k1 - 1;
#line 777 "strsyl.f"
		    suml = sdot_(&i__2, &a[k2 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
#line 778 "strsyl.f"
		    i__2 = *n - l2;
/* Computing MIN */
#line 778 "strsyl.f"
		    i__3 = l2 + 1;
/* Computing MIN */
#line 778 "strsyl.f"
		    i__4 = l2 + 1;
#line 778 "strsyl.f"
		    sumr = sdot_(&i__2, &c__[k2 + min(i__3,*n) * c_dim1], ldc,
			     &b[l2 + min(i__4,*n) * b_dim1], ldb);
#line 780 "strsyl.f"
		    vec[3] = c__[k2 + l2 * c_dim1] - (suml + sgn * sumr);

#line 782 "strsyl.f"
		    slasy2_(&c_true, &c_true, isgn, &c__2, &c__2, &a[k1 + k1 *
			     a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec, &
			    c__2, &scaloc, x, &c__2, &xnorm, &ierr);
#line 785 "strsyl.f"
		    if (ierr != 0) {
#line 785 "strsyl.f"
			*info = 1;
#line 785 "strsyl.f"
		    }

#line 788 "strsyl.f"
		    if (scaloc != 1.) {
#line 789 "strsyl.f"
			i__2 = *n;
#line 789 "strsyl.f"
			for (j = 1; j <= i__2; ++j) {
#line 790 "strsyl.f"
			    sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 791 "strsyl.f"
/* L170: */
#line 791 "strsyl.f"
			}
#line 792 "strsyl.f"
			*scale *= scaloc;
#line 793 "strsyl.f"
		    }
#line 794 "strsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 795 "strsyl.f"
		    c__[k1 + l2 * c_dim1] = x[2];
#line 796 "strsyl.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 797 "strsyl.f"
		    c__[k2 + l2 * c_dim1] = x[3];
#line 798 "strsyl.f"
		}

#line 800 "strsyl.f"
L180:
#line 800 "strsyl.f"
		;
#line 800 "strsyl.f"
	    }
#line 801 "strsyl.f"
L190:
#line 801 "strsyl.f"
	    ;
#line 801 "strsyl.f"
	}

#line 803 "strsyl.f"
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

#line 820 "strsyl.f"
	lnext = *n;
#line 821 "strsyl.f"
	for (l = *n; l >= 1; --l) {
#line 822 "strsyl.f"
	    if (l > lnext) {
#line 822 "strsyl.f"
		goto L250;
#line 822 "strsyl.f"
	    }
#line 824 "strsyl.f"
	    if (l == 1) {
#line 825 "strsyl.f"
		l1 = l;
#line 826 "strsyl.f"
		l2 = l;
#line 827 "strsyl.f"
	    } else {
#line 828 "strsyl.f"
		if (b[l + (l - 1) * b_dim1] != 0.) {
#line 829 "strsyl.f"
		    l1 = l - 1;
#line 830 "strsyl.f"
		    l2 = l;
#line 831 "strsyl.f"
		    lnext = l - 2;
#line 832 "strsyl.f"
		} else {
#line 833 "strsyl.f"
		    l1 = l;
#line 834 "strsyl.f"
		    l2 = l;
#line 835 "strsyl.f"
		    lnext = l - 1;
#line 836 "strsyl.f"
		}
#line 837 "strsyl.f"
	    }

/*           Start row loop (index = K) */
/*           K1 (K2): row index of the first (last) row of X(K,L) */

#line 842 "strsyl.f"
	    knext = *m;
#line 843 "strsyl.f"
	    for (k = *m; k >= 1; --k) {
#line 844 "strsyl.f"
		if (k > knext) {
#line 844 "strsyl.f"
		    goto L240;
#line 844 "strsyl.f"
		}
#line 846 "strsyl.f"
		if (k == 1) {
#line 847 "strsyl.f"
		    k1 = k;
#line 848 "strsyl.f"
		    k2 = k;
#line 849 "strsyl.f"
		} else {
#line 850 "strsyl.f"
		    if (a[k + (k - 1) * a_dim1] != 0.) {
#line 851 "strsyl.f"
			k1 = k - 1;
#line 852 "strsyl.f"
			k2 = k;
#line 853 "strsyl.f"
			knext = k - 2;
#line 854 "strsyl.f"
		    } else {
#line 855 "strsyl.f"
			k1 = k;
#line 856 "strsyl.f"
			k2 = k;
#line 857 "strsyl.f"
			knext = k - 1;
#line 858 "strsyl.f"
		    }
#line 859 "strsyl.f"
		}

#line 861 "strsyl.f"
		if (l1 == l2 && k1 == k2) {
#line 862 "strsyl.f"
		    i__1 = *m - k1;
/* Computing MIN */
#line 862 "strsyl.f"
		    i__2 = k1 + 1;
/* Computing MIN */
#line 862 "strsyl.f"
		    i__3 = k1 + 1;
#line 862 "strsyl.f"
		    suml = sdot_(&i__1, &a[k1 + min(i__2,*m) * a_dim1], lda, &
			    c__[min(i__3,*m) + l1 * c_dim1], &c__1);
#line 864 "strsyl.f"
		    i__1 = *n - l1;
/* Computing MIN */
#line 864 "strsyl.f"
		    i__2 = l1 + 1;
/* Computing MIN */
#line 864 "strsyl.f"
		    i__3 = l1 + 1;
#line 864 "strsyl.f"
		    sumr = sdot_(&i__1, &c__[k1 + min(i__2,*n) * c_dim1], ldc,
			     &b[l1 + min(i__3,*n) * b_dim1], ldb);
#line 866 "strsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
#line 867 "strsyl.f"
		    scaloc = 1.;

#line 869 "strsyl.f"
		    a11 = a[k1 + k1 * a_dim1] + sgn * b[l1 + l1 * b_dim1];
#line 870 "strsyl.f"
		    da11 = abs(a11);
#line 871 "strsyl.f"
		    if (da11 <= smin) {
#line 872 "strsyl.f"
			a11 = smin;
#line 873 "strsyl.f"
			da11 = smin;
#line 874 "strsyl.f"
			*info = 1;
#line 875 "strsyl.f"
		    }
#line 876 "strsyl.f"
		    db = abs(vec[0]);
#line 877 "strsyl.f"
		    if (da11 < 1. && db > 1.) {
#line 878 "strsyl.f"
			if (db > bignum * da11) {
#line 878 "strsyl.f"
			    scaloc = 1. / db;
#line 878 "strsyl.f"
			}
#line 880 "strsyl.f"
		    }
#line 881 "strsyl.f"
		    x[0] = vec[0] * scaloc / a11;

#line 883 "strsyl.f"
		    if (scaloc != 1.) {
#line 884 "strsyl.f"
			i__1 = *n;
#line 884 "strsyl.f"
			for (j = 1; j <= i__1; ++j) {
#line 885 "strsyl.f"
			    sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 886 "strsyl.f"
/* L200: */
#line 886 "strsyl.f"
			}
#line 887 "strsyl.f"
			*scale *= scaloc;
#line 888 "strsyl.f"
		    }
#line 889 "strsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];

#line 891 "strsyl.f"
		} else if (l1 == l2 && k1 != k2) {

#line 893 "strsyl.f"
		    i__1 = *m - k2;
/* Computing MIN */
#line 893 "strsyl.f"
		    i__2 = k2 + 1;
/* Computing MIN */
#line 893 "strsyl.f"
		    i__3 = k2 + 1;
#line 893 "strsyl.f"
		    suml = sdot_(&i__1, &a[k1 + min(i__2,*m) * a_dim1], lda, &
			    c__[min(i__3,*m) + l1 * c_dim1], &c__1);
#line 895 "strsyl.f"
		    i__1 = *n - l2;
/* Computing MIN */
#line 895 "strsyl.f"
		    i__2 = l2 + 1;
/* Computing MIN */
#line 895 "strsyl.f"
		    i__3 = l2 + 1;
#line 895 "strsyl.f"
		    sumr = sdot_(&i__1, &c__[k1 + min(i__2,*n) * c_dim1], ldc,
			     &b[l1 + min(i__3,*n) * b_dim1], ldb);
#line 897 "strsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

#line 899 "strsyl.f"
		    i__1 = *m - k2;
/* Computing MIN */
#line 899 "strsyl.f"
		    i__2 = k2 + 1;
/* Computing MIN */
#line 899 "strsyl.f"
		    i__3 = k2 + 1;
#line 899 "strsyl.f"
		    suml = sdot_(&i__1, &a[k2 + min(i__2,*m) * a_dim1], lda, &
			    c__[min(i__3,*m) + l1 * c_dim1], &c__1);
#line 901 "strsyl.f"
		    i__1 = *n - l2;
/* Computing MIN */
#line 901 "strsyl.f"
		    i__2 = l2 + 1;
/* Computing MIN */
#line 901 "strsyl.f"
		    i__3 = l2 + 1;
#line 901 "strsyl.f"
		    sumr = sdot_(&i__1, &c__[k2 + min(i__2,*n) * c_dim1], ldc,
			     &b[l1 + min(i__3,*n) * b_dim1], ldb);
#line 903 "strsyl.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

#line 905 "strsyl.f"
		    d__1 = -sgn * b[l1 + l1 * b_dim1];
#line 905 "strsyl.f"
		    slaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &a[k1 + k1 
			    * a_dim1], lda, &c_b26, &c_b26, vec, &c__2, &d__1,
			     &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
#line 908 "strsyl.f"
		    if (ierr != 0) {
#line 908 "strsyl.f"
			*info = 1;
#line 908 "strsyl.f"
		    }

#line 911 "strsyl.f"
		    if (scaloc != 1.) {
#line 912 "strsyl.f"
			i__1 = *n;
#line 912 "strsyl.f"
			for (j = 1; j <= i__1; ++j) {
#line 913 "strsyl.f"
			    sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 914 "strsyl.f"
/* L210: */
#line 914 "strsyl.f"
			}
#line 915 "strsyl.f"
			*scale *= scaloc;
#line 916 "strsyl.f"
		    }
#line 917 "strsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 918 "strsyl.f"
		    c__[k2 + l1 * c_dim1] = x[1];

#line 920 "strsyl.f"
		} else if (l1 != l2 && k1 == k2) {

#line 922 "strsyl.f"
		    i__1 = *m - k1;
/* Computing MIN */
#line 922 "strsyl.f"
		    i__2 = k1 + 1;
/* Computing MIN */
#line 922 "strsyl.f"
		    i__3 = k1 + 1;
#line 922 "strsyl.f"
		    suml = sdot_(&i__1, &a[k1 + min(i__2,*m) * a_dim1], lda, &
			    c__[min(i__3,*m) + l1 * c_dim1], &c__1);
#line 924 "strsyl.f"
		    i__1 = *n - l2;
/* Computing MIN */
#line 924 "strsyl.f"
		    i__2 = l2 + 1;
/* Computing MIN */
#line 924 "strsyl.f"
		    i__3 = l2 + 1;
#line 924 "strsyl.f"
		    sumr = sdot_(&i__1, &c__[k1 + min(i__2,*n) * c_dim1], ldc,
			     &b[l1 + min(i__3,*n) * b_dim1], ldb);
#line 926 "strsyl.f"
		    vec[0] = sgn * (c__[k1 + l1 * c_dim1] - (suml + sgn * 
			    sumr));

#line 928 "strsyl.f"
		    i__1 = *m - k1;
/* Computing MIN */
#line 928 "strsyl.f"
		    i__2 = k1 + 1;
/* Computing MIN */
#line 928 "strsyl.f"
		    i__3 = k1 + 1;
#line 928 "strsyl.f"
		    suml = sdot_(&i__1, &a[k1 + min(i__2,*m) * a_dim1], lda, &
			    c__[min(i__3,*m) + l2 * c_dim1], &c__1);
#line 930 "strsyl.f"
		    i__1 = *n - l2;
/* Computing MIN */
#line 930 "strsyl.f"
		    i__2 = l2 + 1;
/* Computing MIN */
#line 930 "strsyl.f"
		    i__3 = l2 + 1;
#line 930 "strsyl.f"
		    sumr = sdot_(&i__1, &c__[k1 + min(i__2,*n) * c_dim1], ldc,
			     &b[l2 + min(i__3,*n) * b_dim1], ldb);
#line 932 "strsyl.f"
		    vec[1] = sgn * (c__[k1 + l2 * c_dim1] - (suml + sgn * 
			    sumr));

#line 934 "strsyl.f"
		    d__1 = -sgn * a[k1 + k1 * a_dim1];
#line 934 "strsyl.f"
		    slaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &b[l1 + l1 
			    * b_dim1], ldb, &c_b26, &c_b26, vec, &c__2, &d__1,
			     &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
#line 937 "strsyl.f"
		    if (ierr != 0) {
#line 937 "strsyl.f"
			*info = 1;
#line 937 "strsyl.f"
		    }

#line 940 "strsyl.f"
		    if (scaloc != 1.) {
#line 941 "strsyl.f"
			i__1 = *n;
#line 941 "strsyl.f"
			for (j = 1; j <= i__1; ++j) {
#line 942 "strsyl.f"
			    sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 943 "strsyl.f"
/* L220: */
#line 943 "strsyl.f"
			}
#line 944 "strsyl.f"
			*scale *= scaloc;
#line 945 "strsyl.f"
		    }
#line 946 "strsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 947 "strsyl.f"
		    c__[k1 + l2 * c_dim1] = x[1];

#line 949 "strsyl.f"
		} else if (l1 != l2 && k1 != k2) {

#line 951 "strsyl.f"
		    i__1 = *m - k2;
/* Computing MIN */
#line 951 "strsyl.f"
		    i__2 = k2 + 1;
/* Computing MIN */
#line 951 "strsyl.f"
		    i__3 = k2 + 1;
#line 951 "strsyl.f"
		    suml = sdot_(&i__1, &a[k1 + min(i__2,*m) * a_dim1], lda, &
			    c__[min(i__3,*m) + l1 * c_dim1], &c__1);
#line 953 "strsyl.f"
		    i__1 = *n - l2;
/* Computing MIN */
#line 953 "strsyl.f"
		    i__2 = l2 + 1;
/* Computing MIN */
#line 953 "strsyl.f"
		    i__3 = l2 + 1;
#line 953 "strsyl.f"
		    sumr = sdot_(&i__1, &c__[k1 + min(i__2,*n) * c_dim1], ldc,
			     &b[l1 + min(i__3,*n) * b_dim1], ldb);
#line 955 "strsyl.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

#line 957 "strsyl.f"
		    i__1 = *m - k2;
/* Computing MIN */
#line 957 "strsyl.f"
		    i__2 = k2 + 1;
/* Computing MIN */
#line 957 "strsyl.f"
		    i__3 = k2 + 1;
#line 957 "strsyl.f"
		    suml = sdot_(&i__1, &a[k1 + min(i__2,*m) * a_dim1], lda, &
			    c__[min(i__3,*m) + l2 * c_dim1], &c__1);
#line 959 "strsyl.f"
		    i__1 = *n - l2;
/* Computing MIN */
#line 959 "strsyl.f"
		    i__2 = l2 + 1;
/* Computing MIN */
#line 959 "strsyl.f"
		    i__3 = l2 + 1;
#line 959 "strsyl.f"
		    sumr = sdot_(&i__1, &c__[k1 + min(i__2,*n) * c_dim1], ldc,
			     &b[l2 + min(i__3,*n) * b_dim1], ldb);
#line 961 "strsyl.f"
		    vec[2] = c__[k1 + l2 * c_dim1] - (suml + sgn * sumr);

#line 963 "strsyl.f"
		    i__1 = *m - k2;
/* Computing MIN */
#line 963 "strsyl.f"
		    i__2 = k2 + 1;
/* Computing MIN */
#line 963 "strsyl.f"
		    i__3 = k2 + 1;
#line 963 "strsyl.f"
		    suml = sdot_(&i__1, &a[k2 + min(i__2,*m) * a_dim1], lda, &
			    c__[min(i__3,*m) + l1 * c_dim1], &c__1);
#line 965 "strsyl.f"
		    i__1 = *n - l2;
/* Computing MIN */
#line 965 "strsyl.f"
		    i__2 = l2 + 1;
/* Computing MIN */
#line 965 "strsyl.f"
		    i__3 = l2 + 1;
#line 965 "strsyl.f"
		    sumr = sdot_(&i__1, &c__[k2 + min(i__2,*n) * c_dim1], ldc,
			     &b[l1 + min(i__3,*n) * b_dim1], ldb);
#line 967 "strsyl.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

#line 969 "strsyl.f"
		    i__1 = *m - k2;
/* Computing MIN */
#line 969 "strsyl.f"
		    i__2 = k2 + 1;
/* Computing MIN */
#line 969 "strsyl.f"
		    i__3 = k2 + 1;
#line 969 "strsyl.f"
		    suml = sdot_(&i__1, &a[k2 + min(i__2,*m) * a_dim1], lda, &
			    c__[min(i__3,*m) + l2 * c_dim1], &c__1);
#line 971 "strsyl.f"
		    i__1 = *n - l2;
/* Computing MIN */
#line 971 "strsyl.f"
		    i__2 = l2 + 1;
/* Computing MIN */
#line 971 "strsyl.f"
		    i__3 = l2 + 1;
#line 971 "strsyl.f"
		    sumr = sdot_(&i__1, &c__[k2 + min(i__2,*n) * c_dim1], ldc,
			     &b[l2 + min(i__3,*n) * b_dim1], ldb);
#line 973 "strsyl.f"
		    vec[3] = c__[k2 + l2 * c_dim1] - (suml + sgn * sumr);

#line 975 "strsyl.f"
		    slasy2_(&c_false, &c_true, isgn, &c__2, &c__2, &a[k1 + k1 
			    * a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec, &
			    c__2, &scaloc, x, &c__2, &xnorm, &ierr);
#line 978 "strsyl.f"
		    if (ierr != 0) {
#line 978 "strsyl.f"
			*info = 1;
#line 978 "strsyl.f"
		    }

#line 981 "strsyl.f"
		    if (scaloc != 1.) {
#line 982 "strsyl.f"
			i__1 = *n;
#line 982 "strsyl.f"
			for (j = 1; j <= i__1; ++j) {
#line 983 "strsyl.f"
			    sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 984 "strsyl.f"
/* L230: */
#line 984 "strsyl.f"
			}
#line 985 "strsyl.f"
			*scale *= scaloc;
#line 986 "strsyl.f"
		    }
#line 987 "strsyl.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 988 "strsyl.f"
		    c__[k1 + l2 * c_dim1] = x[2];
#line 989 "strsyl.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 990 "strsyl.f"
		    c__[k2 + l2 * c_dim1] = x[3];
#line 991 "strsyl.f"
		}

#line 993 "strsyl.f"
L240:
#line 993 "strsyl.f"
		;
#line 993 "strsyl.f"
	    }
#line 994 "strsyl.f"
L250:
#line 994 "strsyl.f"
	    ;
#line 994 "strsyl.f"
	}

#line 996 "strsyl.f"
    }

#line 998 "strsyl.f"
    return 0;

/*     End of STRSYL */

} /* strsyl_ */

