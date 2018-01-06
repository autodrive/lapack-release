#line 1 "ctrsyl.f"
/* ctrsyl.f -- translated by f2c (version 20100827).
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

#line 1 "ctrsyl.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CTRSYL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTRSYL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrsyl.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrsyl.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrsyl.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, */
/*                          LDC, SCALE, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANA, TRANB */
/*       INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N */
/*       REAL               SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTRSYL solves the complex Sylvester matrix equation: */
/* > */
/* >    op(A)*X + X*op(B) = scale*C or */
/* >    op(A)*X - X*op(B) = scale*C, */
/* > */
/* > where op(A) = A or A**H, and A and B are both upper triangular. A is */
/* > M-by-M and B is N-by-N; the right hand side C and the solution X are */
/* > M-by-N; and scale is an output scale factor, set <= 1 to avoid */
/* > overflow in X. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANA */
/* > \verbatim */
/* >          TRANA is CHARACTER*1 */
/* >          Specifies the option op(A): */
/* >          = 'N': op(A) = A    (No transpose) */
/* >          = 'C': op(A) = A**H (Conjugate transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] TRANB */
/* > \verbatim */
/* >          TRANB is CHARACTER*1 */
/* >          Specifies the option op(B): */
/* >          = 'N': op(B) = B    (No transpose) */
/* >          = 'C': op(B) = B**H (Conjugate transpose) */
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
/* >          A is COMPLEX array, dimension (LDA,M) */
/* >          The upper triangular matrix A. */
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
/* >          B is COMPLEX array, dimension (LDB,N) */
/* >          The upper triangular matrix B. */
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
/* >          C is COMPLEX array, dimension (LDC,N) */
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

/* > \ingroup complexSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int ctrsyl_(char *trana, char *tranb, integer *isgn, integer 
	*m, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublecomplex *c__, integer *ldc, doublereal *scale, 
	integer *info, ftnlen trana_len, ftnlen tranb_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k, l;
    static doublecomplex a11;
    static doublereal db;
    static doublecomplex x11;
    static doublereal da11;
    static doublecomplex vec;
    static doublereal dum[1], eps, sgn, smin;
    static doublecomplex suml, sumr;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Double Complex */ VOID cdotu_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int slabad_(doublereal *, doublereal *);
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Double Complex */ VOID cladiv_(doublecomplex *, doublecomplex *,
	     doublecomplex *);
    static doublereal scaloc;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), xerbla_(char *, integer *, ftnlen);
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

#line 206 "ctrsyl.f"
    /* Parameter adjustments */
#line 206 "ctrsyl.f"
    a_dim1 = *lda;
#line 206 "ctrsyl.f"
    a_offset = 1 + a_dim1;
#line 206 "ctrsyl.f"
    a -= a_offset;
#line 206 "ctrsyl.f"
    b_dim1 = *ldb;
#line 206 "ctrsyl.f"
    b_offset = 1 + b_dim1;
#line 206 "ctrsyl.f"
    b -= b_offset;
#line 206 "ctrsyl.f"
    c_dim1 = *ldc;
#line 206 "ctrsyl.f"
    c_offset = 1 + c_dim1;
#line 206 "ctrsyl.f"
    c__ -= c_offset;
#line 206 "ctrsyl.f"

#line 206 "ctrsyl.f"
    /* Function Body */
#line 206 "ctrsyl.f"
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
#line 207 "ctrsyl.f"
    notrnb = lsame_(tranb, "N", (ftnlen)1, (ftnlen)1);

#line 209 "ctrsyl.f"
    *info = 0;
#line 210 "ctrsyl.f"
    if (! notrna && ! lsame_(trana, "C", (ftnlen)1, (ftnlen)1)) {
#line 211 "ctrsyl.f"
	*info = -1;
#line 212 "ctrsyl.f"
    } else if (! notrnb && ! lsame_(tranb, "C", (ftnlen)1, (ftnlen)1)) {
#line 213 "ctrsyl.f"
	*info = -2;
#line 214 "ctrsyl.f"
    } else if (*isgn != 1 && *isgn != -1) {
#line 215 "ctrsyl.f"
	*info = -3;
#line 216 "ctrsyl.f"
    } else if (*m < 0) {
#line 217 "ctrsyl.f"
	*info = -4;
#line 218 "ctrsyl.f"
    } else if (*n < 0) {
#line 219 "ctrsyl.f"
	*info = -5;
#line 220 "ctrsyl.f"
    } else if (*lda < max(1,*m)) {
#line 221 "ctrsyl.f"
	*info = -7;
#line 222 "ctrsyl.f"
    } else if (*ldb < max(1,*n)) {
#line 223 "ctrsyl.f"
	*info = -9;
#line 224 "ctrsyl.f"
    } else if (*ldc < max(1,*m)) {
#line 225 "ctrsyl.f"
	*info = -11;
#line 226 "ctrsyl.f"
    }
#line 227 "ctrsyl.f"
    if (*info != 0) {
#line 228 "ctrsyl.f"
	i__1 = -(*info);
#line 228 "ctrsyl.f"
	xerbla_("CTRSYL", &i__1, (ftnlen)6);
#line 229 "ctrsyl.f"
	return 0;
#line 230 "ctrsyl.f"
    }

/*     Quick return if possible */

#line 234 "ctrsyl.f"
    *scale = 1.;
#line 235 "ctrsyl.f"
    if (*m == 0 || *n == 0) {
#line 235 "ctrsyl.f"
	return 0;
#line 235 "ctrsyl.f"
    }

/*     Set constants to control overflow */

#line 240 "ctrsyl.f"
    eps = slamch_("P", (ftnlen)1);
#line 241 "ctrsyl.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 242 "ctrsyl.f"
    bignum = 1. / smlnum;
#line 243 "ctrsyl.f"
    slabad_(&smlnum, &bignum);
#line 244 "ctrsyl.f"
    smlnum = smlnum * (doublereal) (*m * *n) / eps;
#line 245 "ctrsyl.f"
    bignum = 1. / smlnum;
/* Computing MAX */
#line 246 "ctrsyl.f"
    d__1 = smlnum, d__2 = eps * clange_("M", m, m, &a[a_offset], lda, dum, (
	    ftnlen)1), d__1 = max(d__1,d__2), d__2 = eps * clange_("M", n, n, 
	    &b[b_offset], ldb, dum, (ftnlen)1);
#line 246 "ctrsyl.f"
    smin = max(d__1,d__2);
#line 248 "ctrsyl.f"
    sgn = (doublereal) (*isgn);

#line 250 "ctrsyl.f"
    if (notrna && notrnb) {

/*        Solve    A*X + ISGN*X*B = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        bottom-left corner column by column by */

/*            A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L) */

/*        Where */
/*                    M                        L-1 */
/*          R(K,L) = SUM [A(K,I)*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]. */
/*                  I=K+1                      J=1 */

#line 264 "ctrsyl.f"
	i__1 = *n;
#line 264 "ctrsyl.f"
	for (l = 1; l <= i__1; ++l) {
#line 265 "ctrsyl.f"
	    for (k = *m; k >= 1; --k) {

#line 267 "ctrsyl.f"
		i__2 = *m - k;
/* Computing MIN */
#line 267 "ctrsyl.f"
		i__3 = k + 1;
/* Computing MIN */
#line 267 "ctrsyl.f"
		i__4 = k + 1;
#line 267 "ctrsyl.f"
		cdotu_(&z__1, &i__2, &a[k + min(i__3,*m) * a_dim1], lda, &c__[
			min(i__4,*m) + l * c_dim1], &c__1);
#line 267 "ctrsyl.f"
		suml.r = z__1.r, suml.i = z__1.i;
#line 269 "ctrsyl.f"
		i__2 = l - 1;
#line 269 "ctrsyl.f"
		cdotu_(&z__1, &i__2, &c__[k + c_dim1], ldc, &b[l * b_dim1 + 1]
			, &c__1);
#line 269 "ctrsyl.f"
		sumr.r = z__1.r, sumr.i = z__1.i;
#line 270 "ctrsyl.f"
		i__2 = k + l * c_dim1;
#line 270 "ctrsyl.f"
		z__3.r = sgn * sumr.r, z__3.i = sgn * sumr.i;
#line 270 "ctrsyl.f"
		z__2.r = suml.r + z__3.r, z__2.i = suml.i + z__3.i;
#line 270 "ctrsyl.f"
		z__1.r = c__[i__2].r - z__2.r, z__1.i = c__[i__2].i - z__2.i;
#line 270 "ctrsyl.f"
		vec.r = z__1.r, vec.i = z__1.i;

#line 272 "ctrsyl.f"
		scaloc = 1.;
#line 273 "ctrsyl.f"
		i__2 = k + k * a_dim1;
#line 273 "ctrsyl.f"
		i__3 = l + l * b_dim1;
#line 273 "ctrsyl.f"
		z__2.r = sgn * b[i__3].r, z__2.i = sgn * b[i__3].i;
#line 273 "ctrsyl.f"
		z__1.r = a[i__2].r + z__2.r, z__1.i = a[i__2].i + z__2.i;
#line 273 "ctrsyl.f"
		a11.r = z__1.r, a11.i = z__1.i;
#line 274 "ctrsyl.f"
		da11 = (d__1 = a11.r, abs(d__1)) + (d__2 = d_imag(&a11), abs(
			d__2));
#line 275 "ctrsyl.f"
		if (da11 <= smin) {
#line 276 "ctrsyl.f"
		    a11.r = smin, a11.i = 0.;
#line 277 "ctrsyl.f"
		    da11 = smin;
#line 278 "ctrsyl.f"
		    *info = 1;
#line 279 "ctrsyl.f"
		}
#line 280 "ctrsyl.f"
		db = (d__1 = vec.r, abs(d__1)) + (d__2 = d_imag(&vec), abs(
			d__2));
#line 281 "ctrsyl.f"
		if (da11 < 1. && db > 1.) {
#line 282 "ctrsyl.f"
		    if (db > bignum * da11) {
#line 282 "ctrsyl.f"
			scaloc = 1. / db;
#line 282 "ctrsyl.f"
		    }
#line 284 "ctrsyl.f"
		}
#line 285 "ctrsyl.f"
		z__3.r = scaloc, z__3.i = 0.;
#line 285 "ctrsyl.f"
		z__2.r = vec.r * z__3.r - vec.i * z__3.i, z__2.i = vec.r * 
			z__3.i + vec.i * z__3.r;
#line 285 "ctrsyl.f"
		cladiv_(&z__1, &z__2, &a11);
#line 285 "ctrsyl.f"
		x11.r = z__1.r, x11.i = z__1.i;

#line 287 "ctrsyl.f"
		if (scaloc != 1.) {
#line 288 "ctrsyl.f"
		    i__2 = *n;
#line 288 "ctrsyl.f"
		    for (j = 1; j <= i__2; ++j) {
#line 289 "ctrsyl.f"
			csscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 290 "ctrsyl.f"
/* L10: */
#line 290 "ctrsyl.f"
		    }
#line 291 "ctrsyl.f"
		    *scale *= scaloc;
#line 292 "ctrsyl.f"
		}
#line 293 "ctrsyl.f"
		i__2 = k + l * c_dim1;
#line 293 "ctrsyl.f"
		c__[i__2].r = x11.r, c__[i__2].i = x11.i;

#line 295 "ctrsyl.f"
/* L20: */
#line 295 "ctrsyl.f"
	    }
#line 296 "ctrsyl.f"
/* L30: */
#line 296 "ctrsyl.f"
	}

#line 298 "ctrsyl.f"
    } else if (! notrna && notrnb) {

/*        Solve    A**H *X + ISGN*X*B = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        upper-left corner column by column by */

/*            A**H(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L) */

/*        Where */
/*                   K-1                           L-1 */
/*          R(K,L) = SUM [A**H(I,K)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)] */
/*                   I=1                           J=1 */

#line 312 "ctrsyl.f"
	i__1 = *n;
#line 312 "ctrsyl.f"
	for (l = 1; l <= i__1; ++l) {
#line 313 "ctrsyl.f"
	    i__2 = *m;
#line 313 "ctrsyl.f"
	    for (k = 1; k <= i__2; ++k) {

#line 315 "ctrsyl.f"
		i__3 = k - 1;
#line 315 "ctrsyl.f"
		cdotc_(&z__1, &i__3, &a[k * a_dim1 + 1], &c__1, &c__[l * 
			c_dim1 + 1], &c__1);
#line 315 "ctrsyl.f"
		suml.r = z__1.r, suml.i = z__1.i;
#line 316 "ctrsyl.f"
		i__3 = l - 1;
#line 316 "ctrsyl.f"
		cdotu_(&z__1, &i__3, &c__[k + c_dim1], ldc, &b[l * b_dim1 + 1]
			, &c__1);
#line 316 "ctrsyl.f"
		sumr.r = z__1.r, sumr.i = z__1.i;
#line 317 "ctrsyl.f"
		i__3 = k + l * c_dim1;
#line 317 "ctrsyl.f"
		z__3.r = sgn * sumr.r, z__3.i = sgn * sumr.i;
#line 317 "ctrsyl.f"
		z__2.r = suml.r + z__3.r, z__2.i = suml.i + z__3.i;
#line 317 "ctrsyl.f"
		z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 317 "ctrsyl.f"
		vec.r = z__1.r, vec.i = z__1.i;

#line 319 "ctrsyl.f"
		scaloc = 1.;
#line 320 "ctrsyl.f"
		d_cnjg(&z__2, &a[k + k * a_dim1]);
#line 320 "ctrsyl.f"
		i__3 = l + l * b_dim1;
#line 320 "ctrsyl.f"
		z__3.r = sgn * b[i__3].r, z__3.i = sgn * b[i__3].i;
#line 320 "ctrsyl.f"
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 320 "ctrsyl.f"
		a11.r = z__1.r, a11.i = z__1.i;
#line 321 "ctrsyl.f"
		da11 = (d__1 = a11.r, abs(d__1)) + (d__2 = d_imag(&a11), abs(
			d__2));
#line 322 "ctrsyl.f"
		if (da11 <= smin) {
#line 323 "ctrsyl.f"
		    a11.r = smin, a11.i = 0.;
#line 324 "ctrsyl.f"
		    da11 = smin;
#line 325 "ctrsyl.f"
		    *info = 1;
#line 326 "ctrsyl.f"
		}
#line 327 "ctrsyl.f"
		db = (d__1 = vec.r, abs(d__1)) + (d__2 = d_imag(&vec), abs(
			d__2));
#line 328 "ctrsyl.f"
		if (da11 < 1. && db > 1.) {
#line 329 "ctrsyl.f"
		    if (db > bignum * da11) {
#line 329 "ctrsyl.f"
			scaloc = 1. / db;
#line 329 "ctrsyl.f"
		    }
#line 331 "ctrsyl.f"
		}

#line 333 "ctrsyl.f"
		z__3.r = scaloc, z__3.i = 0.;
#line 333 "ctrsyl.f"
		z__2.r = vec.r * z__3.r - vec.i * z__3.i, z__2.i = vec.r * 
			z__3.i + vec.i * z__3.r;
#line 333 "ctrsyl.f"
		cladiv_(&z__1, &z__2, &a11);
#line 333 "ctrsyl.f"
		x11.r = z__1.r, x11.i = z__1.i;

#line 335 "ctrsyl.f"
		if (scaloc != 1.) {
#line 336 "ctrsyl.f"
		    i__3 = *n;
#line 336 "ctrsyl.f"
		    for (j = 1; j <= i__3; ++j) {
#line 337 "ctrsyl.f"
			csscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 338 "ctrsyl.f"
/* L40: */
#line 338 "ctrsyl.f"
		    }
#line 339 "ctrsyl.f"
		    *scale *= scaloc;
#line 340 "ctrsyl.f"
		}
#line 341 "ctrsyl.f"
		i__3 = k + l * c_dim1;
#line 341 "ctrsyl.f"
		c__[i__3].r = x11.r, c__[i__3].i = x11.i;

#line 343 "ctrsyl.f"
/* L50: */
#line 343 "ctrsyl.f"
	    }
#line 344 "ctrsyl.f"
/* L60: */
#line 344 "ctrsyl.f"
	}

#line 346 "ctrsyl.f"
    } else if (! notrna && ! notrnb) {

/*        Solve    A**H*X + ISGN*X*B**H = C. */

/*        The (K,L)th block of X is determined starting from */
/*        upper-right corner column by column by */

/*            A**H(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L) */

/*        Where */
/*                    K-1 */
/*           R(K,L) = SUM [A**H(I,K)*X(I,L)] + */
/*                    I=1 */
/*                           N */
/*                     ISGN*SUM [X(K,J)*B**H(L,J)]. */
/*                          J=L+1 */

#line 363 "ctrsyl.f"
	for (l = *n; l >= 1; --l) {
#line 364 "ctrsyl.f"
	    i__1 = *m;
#line 364 "ctrsyl.f"
	    for (k = 1; k <= i__1; ++k) {

#line 366 "ctrsyl.f"
		i__2 = k - 1;
#line 366 "ctrsyl.f"
		cdotc_(&z__1, &i__2, &a[k * a_dim1 + 1], &c__1, &c__[l * 
			c_dim1 + 1], &c__1);
#line 366 "ctrsyl.f"
		suml.r = z__1.r, suml.i = z__1.i;
#line 367 "ctrsyl.f"
		i__2 = *n - l;
/* Computing MIN */
#line 367 "ctrsyl.f"
		i__3 = l + 1;
/* Computing MIN */
#line 367 "ctrsyl.f"
		i__4 = l + 1;
#line 367 "ctrsyl.f"
		cdotc_(&z__1, &i__2, &c__[k + min(i__3,*n) * c_dim1], ldc, &b[
			l + min(i__4,*n) * b_dim1], ldb);
#line 367 "ctrsyl.f"
		sumr.r = z__1.r, sumr.i = z__1.i;
#line 369 "ctrsyl.f"
		i__2 = k + l * c_dim1;
#line 369 "ctrsyl.f"
		d_cnjg(&z__4, &sumr);
#line 369 "ctrsyl.f"
		z__3.r = sgn * z__4.r, z__3.i = sgn * z__4.i;
#line 369 "ctrsyl.f"
		z__2.r = suml.r + z__3.r, z__2.i = suml.i + z__3.i;
#line 369 "ctrsyl.f"
		z__1.r = c__[i__2].r - z__2.r, z__1.i = c__[i__2].i - z__2.i;
#line 369 "ctrsyl.f"
		vec.r = z__1.r, vec.i = z__1.i;

#line 371 "ctrsyl.f"
		scaloc = 1.;
#line 372 "ctrsyl.f"
		i__2 = k + k * a_dim1;
#line 372 "ctrsyl.f"
		i__3 = l + l * b_dim1;
#line 372 "ctrsyl.f"
		z__3.r = sgn * b[i__3].r, z__3.i = sgn * b[i__3].i;
#line 372 "ctrsyl.f"
		z__2.r = a[i__2].r + z__3.r, z__2.i = a[i__2].i + z__3.i;
#line 372 "ctrsyl.f"
		d_cnjg(&z__1, &z__2);
#line 372 "ctrsyl.f"
		a11.r = z__1.r, a11.i = z__1.i;
#line 373 "ctrsyl.f"
		da11 = (d__1 = a11.r, abs(d__1)) + (d__2 = d_imag(&a11), abs(
			d__2));
#line 374 "ctrsyl.f"
		if (da11 <= smin) {
#line 375 "ctrsyl.f"
		    a11.r = smin, a11.i = 0.;
#line 376 "ctrsyl.f"
		    da11 = smin;
#line 377 "ctrsyl.f"
		    *info = 1;
#line 378 "ctrsyl.f"
		}
#line 379 "ctrsyl.f"
		db = (d__1 = vec.r, abs(d__1)) + (d__2 = d_imag(&vec), abs(
			d__2));
#line 380 "ctrsyl.f"
		if (da11 < 1. && db > 1.) {
#line 381 "ctrsyl.f"
		    if (db > bignum * da11) {
#line 381 "ctrsyl.f"
			scaloc = 1. / db;
#line 381 "ctrsyl.f"
		    }
#line 383 "ctrsyl.f"
		}

#line 385 "ctrsyl.f"
		z__3.r = scaloc, z__3.i = 0.;
#line 385 "ctrsyl.f"
		z__2.r = vec.r * z__3.r - vec.i * z__3.i, z__2.i = vec.r * 
			z__3.i + vec.i * z__3.r;
#line 385 "ctrsyl.f"
		cladiv_(&z__1, &z__2, &a11);
#line 385 "ctrsyl.f"
		x11.r = z__1.r, x11.i = z__1.i;

#line 387 "ctrsyl.f"
		if (scaloc != 1.) {
#line 388 "ctrsyl.f"
		    i__2 = *n;
#line 388 "ctrsyl.f"
		    for (j = 1; j <= i__2; ++j) {
#line 389 "ctrsyl.f"
			csscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 390 "ctrsyl.f"
/* L70: */
#line 390 "ctrsyl.f"
		    }
#line 391 "ctrsyl.f"
		    *scale *= scaloc;
#line 392 "ctrsyl.f"
		}
#line 393 "ctrsyl.f"
		i__2 = k + l * c_dim1;
#line 393 "ctrsyl.f"
		c__[i__2].r = x11.r, c__[i__2].i = x11.i;

#line 395 "ctrsyl.f"
/* L80: */
#line 395 "ctrsyl.f"
	    }
#line 396 "ctrsyl.f"
/* L90: */
#line 396 "ctrsyl.f"
	}

#line 398 "ctrsyl.f"
    } else if (notrna && ! notrnb) {

/*        Solve    A*X + ISGN*X*B**H = C. */

/*        The (K,L)th block of X is determined starting from */
/*        bottom-left corner column by column by */

/*           A(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L) */

/*        Where */
/*                    M                          N */
/*          R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B**H(L,J)] */
/*                  I=K+1                      J=L+1 */

#line 412 "ctrsyl.f"
	for (l = *n; l >= 1; --l) {
#line 413 "ctrsyl.f"
	    for (k = *m; k >= 1; --k) {

#line 415 "ctrsyl.f"
		i__1 = *m - k;
/* Computing MIN */
#line 415 "ctrsyl.f"
		i__2 = k + 1;
/* Computing MIN */
#line 415 "ctrsyl.f"
		i__3 = k + 1;
#line 415 "ctrsyl.f"
		cdotu_(&z__1, &i__1, &a[k + min(i__2,*m) * a_dim1], lda, &c__[
			min(i__3,*m) + l * c_dim1], &c__1);
#line 415 "ctrsyl.f"
		suml.r = z__1.r, suml.i = z__1.i;
#line 417 "ctrsyl.f"
		i__1 = *n - l;
/* Computing MIN */
#line 417 "ctrsyl.f"
		i__2 = l + 1;
/* Computing MIN */
#line 417 "ctrsyl.f"
		i__3 = l + 1;
#line 417 "ctrsyl.f"
		cdotc_(&z__1, &i__1, &c__[k + min(i__2,*n) * c_dim1], ldc, &b[
			l + min(i__3,*n) * b_dim1], ldb);
#line 417 "ctrsyl.f"
		sumr.r = z__1.r, sumr.i = z__1.i;
#line 419 "ctrsyl.f"
		i__1 = k + l * c_dim1;
#line 419 "ctrsyl.f"
		d_cnjg(&z__4, &sumr);
#line 419 "ctrsyl.f"
		z__3.r = sgn * z__4.r, z__3.i = sgn * z__4.i;
#line 419 "ctrsyl.f"
		z__2.r = suml.r + z__3.r, z__2.i = suml.i + z__3.i;
#line 419 "ctrsyl.f"
		z__1.r = c__[i__1].r - z__2.r, z__1.i = c__[i__1].i - z__2.i;
#line 419 "ctrsyl.f"
		vec.r = z__1.r, vec.i = z__1.i;

#line 421 "ctrsyl.f"
		scaloc = 1.;
#line 422 "ctrsyl.f"
		i__1 = k + k * a_dim1;
#line 422 "ctrsyl.f"
		d_cnjg(&z__3, &b[l + l * b_dim1]);
#line 422 "ctrsyl.f"
		z__2.r = sgn * z__3.r, z__2.i = sgn * z__3.i;
#line 422 "ctrsyl.f"
		z__1.r = a[i__1].r + z__2.r, z__1.i = a[i__1].i + z__2.i;
#line 422 "ctrsyl.f"
		a11.r = z__1.r, a11.i = z__1.i;
#line 423 "ctrsyl.f"
		da11 = (d__1 = a11.r, abs(d__1)) + (d__2 = d_imag(&a11), abs(
			d__2));
#line 424 "ctrsyl.f"
		if (da11 <= smin) {
#line 425 "ctrsyl.f"
		    a11.r = smin, a11.i = 0.;
#line 426 "ctrsyl.f"
		    da11 = smin;
#line 427 "ctrsyl.f"
		    *info = 1;
#line 428 "ctrsyl.f"
		}
#line 429 "ctrsyl.f"
		db = (d__1 = vec.r, abs(d__1)) + (d__2 = d_imag(&vec), abs(
			d__2));
#line 430 "ctrsyl.f"
		if (da11 < 1. && db > 1.) {
#line 431 "ctrsyl.f"
		    if (db > bignum * da11) {
#line 431 "ctrsyl.f"
			scaloc = 1. / db;
#line 431 "ctrsyl.f"
		    }
#line 433 "ctrsyl.f"
		}

#line 435 "ctrsyl.f"
		z__3.r = scaloc, z__3.i = 0.;
#line 435 "ctrsyl.f"
		z__2.r = vec.r * z__3.r - vec.i * z__3.i, z__2.i = vec.r * 
			z__3.i + vec.i * z__3.r;
#line 435 "ctrsyl.f"
		cladiv_(&z__1, &z__2, &a11);
#line 435 "ctrsyl.f"
		x11.r = z__1.r, x11.i = z__1.i;

#line 437 "ctrsyl.f"
		if (scaloc != 1.) {
#line 438 "ctrsyl.f"
		    i__1 = *n;
#line 438 "ctrsyl.f"
		    for (j = 1; j <= i__1; ++j) {
#line 439 "ctrsyl.f"
			csscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 440 "ctrsyl.f"
/* L100: */
#line 440 "ctrsyl.f"
		    }
#line 441 "ctrsyl.f"
		    *scale *= scaloc;
#line 442 "ctrsyl.f"
		}
#line 443 "ctrsyl.f"
		i__1 = k + l * c_dim1;
#line 443 "ctrsyl.f"
		c__[i__1].r = x11.r, c__[i__1].i = x11.i;

#line 445 "ctrsyl.f"
/* L110: */
#line 445 "ctrsyl.f"
	    }
#line 446 "ctrsyl.f"
/* L120: */
#line 446 "ctrsyl.f"
	}

#line 448 "ctrsyl.f"
    }

#line 450 "ctrsyl.f"
    return 0;

/*     End of CTRSYL */

} /* ctrsyl_ */

