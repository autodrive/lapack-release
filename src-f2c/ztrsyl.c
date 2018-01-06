#line 1 "ztrsyl.f"
/* ztrsyl.f -- translated by f2c (version 20100827).
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

#line 1 "ztrsyl.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZTRSYL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTRSYL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrsyl.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrsyl.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrsyl.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, */
/*                          LDC, SCALE, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANA, TRANB */
/*       INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N */
/*       DOUBLE PRECISION   SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTRSYL solves the complex Sylvester matrix equation: */
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
/* >          A is COMPLEX*16 array, dimension (LDA,M) */
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
/* >          B is COMPLEX*16 array, dimension (LDB,N) */
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
/* >          C is COMPLEX*16 array, dimension (LDC,N) */
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

/* > \ingroup complex16SYcomputational */

/*  ===================================================================== */
/* Subroutine */ int ztrsyl_(char *trana, char *tranb, integer *isgn, integer 
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
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), zdotu_(
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal scaloc;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int zdscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    extern /* Double Complex */ VOID zladiv_(doublecomplex *, doublecomplex *,
	     doublecomplex *);
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

#line 206 "ztrsyl.f"
    /* Parameter adjustments */
#line 206 "ztrsyl.f"
    a_dim1 = *lda;
#line 206 "ztrsyl.f"
    a_offset = 1 + a_dim1;
#line 206 "ztrsyl.f"
    a -= a_offset;
#line 206 "ztrsyl.f"
    b_dim1 = *ldb;
#line 206 "ztrsyl.f"
    b_offset = 1 + b_dim1;
#line 206 "ztrsyl.f"
    b -= b_offset;
#line 206 "ztrsyl.f"
    c_dim1 = *ldc;
#line 206 "ztrsyl.f"
    c_offset = 1 + c_dim1;
#line 206 "ztrsyl.f"
    c__ -= c_offset;
#line 206 "ztrsyl.f"

#line 206 "ztrsyl.f"
    /* Function Body */
#line 206 "ztrsyl.f"
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
#line 207 "ztrsyl.f"
    notrnb = lsame_(tranb, "N", (ftnlen)1, (ftnlen)1);

#line 209 "ztrsyl.f"
    *info = 0;
#line 210 "ztrsyl.f"
    if (! notrna && ! lsame_(trana, "C", (ftnlen)1, (ftnlen)1)) {
#line 211 "ztrsyl.f"
	*info = -1;
#line 212 "ztrsyl.f"
    } else if (! notrnb && ! lsame_(tranb, "C", (ftnlen)1, (ftnlen)1)) {
#line 213 "ztrsyl.f"
	*info = -2;
#line 214 "ztrsyl.f"
    } else if (*isgn != 1 && *isgn != -1) {
#line 215 "ztrsyl.f"
	*info = -3;
#line 216 "ztrsyl.f"
    } else if (*m < 0) {
#line 217 "ztrsyl.f"
	*info = -4;
#line 218 "ztrsyl.f"
    } else if (*n < 0) {
#line 219 "ztrsyl.f"
	*info = -5;
#line 220 "ztrsyl.f"
    } else if (*lda < max(1,*m)) {
#line 221 "ztrsyl.f"
	*info = -7;
#line 222 "ztrsyl.f"
    } else if (*ldb < max(1,*n)) {
#line 223 "ztrsyl.f"
	*info = -9;
#line 224 "ztrsyl.f"
    } else if (*ldc < max(1,*m)) {
#line 225 "ztrsyl.f"
	*info = -11;
#line 226 "ztrsyl.f"
    }
#line 227 "ztrsyl.f"
    if (*info != 0) {
#line 228 "ztrsyl.f"
	i__1 = -(*info);
#line 228 "ztrsyl.f"
	xerbla_("ZTRSYL", &i__1, (ftnlen)6);
#line 229 "ztrsyl.f"
	return 0;
#line 230 "ztrsyl.f"
    }

/*     Quick return if possible */

#line 234 "ztrsyl.f"
    *scale = 1.;
#line 235 "ztrsyl.f"
    if (*m == 0 || *n == 0) {
#line 235 "ztrsyl.f"
	return 0;
#line 235 "ztrsyl.f"
    }

/*     Set constants to control overflow */

#line 240 "ztrsyl.f"
    eps = dlamch_("P", (ftnlen)1);
#line 241 "ztrsyl.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 242 "ztrsyl.f"
    bignum = 1. / smlnum;
#line 243 "ztrsyl.f"
    dlabad_(&smlnum, &bignum);
#line 244 "ztrsyl.f"
    smlnum = smlnum * (doublereal) (*m * *n) / eps;
#line 245 "ztrsyl.f"
    bignum = 1. / smlnum;
/* Computing MAX */
#line 246 "ztrsyl.f"
    d__1 = smlnum, d__2 = eps * zlange_("M", m, m, &a[a_offset], lda, dum, (
	    ftnlen)1), d__1 = max(d__1,d__2), d__2 = eps * zlange_("M", n, n, 
	    &b[b_offset], ldb, dum, (ftnlen)1);
#line 246 "ztrsyl.f"
    smin = max(d__1,d__2);
#line 248 "ztrsyl.f"
    sgn = (doublereal) (*isgn);

#line 250 "ztrsyl.f"
    if (notrna && notrnb) {

/*        Solve    A*X + ISGN*X*B = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        bottom-left corner column by column by */

/*            A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L) */

/*        Where */
/*                    M                        L-1 */
/*          R(K,L) = SUM [A(K,I)*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]. */
/*                  I=K+1                      J=1 */

#line 264 "ztrsyl.f"
	i__1 = *n;
#line 264 "ztrsyl.f"
	for (l = 1; l <= i__1; ++l) {
#line 265 "ztrsyl.f"
	    for (k = *m; k >= 1; --k) {

#line 267 "ztrsyl.f"
		i__2 = *m - k;
/* Computing MIN */
#line 267 "ztrsyl.f"
		i__3 = k + 1;
/* Computing MIN */
#line 267 "ztrsyl.f"
		i__4 = k + 1;
#line 267 "ztrsyl.f"
		zdotu_(&z__1, &i__2, &a[k + min(i__3,*m) * a_dim1], lda, &c__[
			min(i__4,*m) + l * c_dim1], &c__1);
#line 267 "ztrsyl.f"
		suml.r = z__1.r, suml.i = z__1.i;
#line 269 "ztrsyl.f"
		i__2 = l - 1;
#line 269 "ztrsyl.f"
		zdotu_(&z__1, &i__2, &c__[k + c_dim1], ldc, &b[l * b_dim1 + 1]
			, &c__1);
#line 269 "ztrsyl.f"
		sumr.r = z__1.r, sumr.i = z__1.i;
#line 270 "ztrsyl.f"
		i__2 = k + l * c_dim1;
#line 270 "ztrsyl.f"
		z__3.r = sgn * sumr.r, z__3.i = sgn * sumr.i;
#line 270 "ztrsyl.f"
		z__2.r = suml.r + z__3.r, z__2.i = suml.i + z__3.i;
#line 270 "ztrsyl.f"
		z__1.r = c__[i__2].r - z__2.r, z__1.i = c__[i__2].i - z__2.i;
#line 270 "ztrsyl.f"
		vec.r = z__1.r, vec.i = z__1.i;

#line 272 "ztrsyl.f"
		scaloc = 1.;
#line 273 "ztrsyl.f"
		i__2 = k + k * a_dim1;
#line 273 "ztrsyl.f"
		i__3 = l + l * b_dim1;
#line 273 "ztrsyl.f"
		z__2.r = sgn * b[i__3].r, z__2.i = sgn * b[i__3].i;
#line 273 "ztrsyl.f"
		z__1.r = a[i__2].r + z__2.r, z__1.i = a[i__2].i + z__2.i;
#line 273 "ztrsyl.f"
		a11.r = z__1.r, a11.i = z__1.i;
#line 274 "ztrsyl.f"
		da11 = (d__1 = a11.r, abs(d__1)) + (d__2 = d_imag(&a11), abs(
			d__2));
#line 275 "ztrsyl.f"
		if (da11 <= smin) {
#line 276 "ztrsyl.f"
		    a11.r = smin, a11.i = 0.;
#line 277 "ztrsyl.f"
		    da11 = smin;
#line 278 "ztrsyl.f"
		    *info = 1;
#line 279 "ztrsyl.f"
		}
#line 280 "ztrsyl.f"
		db = (d__1 = vec.r, abs(d__1)) + (d__2 = d_imag(&vec), abs(
			d__2));
#line 281 "ztrsyl.f"
		if (da11 < 1. && db > 1.) {
#line 282 "ztrsyl.f"
		    if (db > bignum * da11) {
#line 282 "ztrsyl.f"
			scaloc = 1. / db;
#line 282 "ztrsyl.f"
		    }
#line 284 "ztrsyl.f"
		}
#line 285 "ztrsyl.f"
		z__3.r = scaloc, z__3.i = 0.;
#line 285 "ztrsyl.f"
		z__2.r = vec.r * z__3.r - vec.i * z__3.i, z__2.i = vec.r * 
			z__3.i + vec.i * z__3.r;
#line 285 "ztrsyl.f"
		zladiv_(&z__1, &z__2, &a11);
#line 285 "ztrsyl.f"
		x11.r = z__1.r, x11.i = z__1.i;

#line 287 "ztrsyl.f"
		if (scaloc != 1.) {
#line 288 "ztrsyl.f"
		    i__2 = *n;
#line 288 "ztrsyl.f"
		    for (j = 1; j <= i__2; ++j) {
#line 289 "ztrsyl.f"
			zdscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 290 "ztrsyl.f"
/* L10: */
#line 290 "ztrsyl.f"
		    }
#line 291 "ztrsyl.f"
		    *scale *= scaloc;
#line 292 "ztrsyl.f"
		}
#line 293 "ztrsyl.f"
		i__2 = k + l * c_dim1;
#line 293 "ztrsyl.f"
		c__[i__2].r = x11.r, c__[i__2].i = x11.i;

#line 295 "ztrsyl.f"
/* L20: */
#line 295 "ztrsyl.f"
	    }
#line 296 "ztrsyl.f"
/* L30: */
#line 296 "ztrsyl.f"
	}

#line 298 "ztrsyl.f"
    } else if (! notrna && notrnb) {

/*        Solve    A**H *X + ISGN*X*B = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        upper-left corner column by column by */

/*            A**H(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L) */

/*        Where */
/*                   K-1                           L-1 */
/*          R(K,L) = SUM [A**H(I,K)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)] */
/*                   I=1                           J=1 */

#line 312 "ztrsyl.f"
	i__1 = *n;
#line 312 "ztrsyl.f"
	for (l = 1; l <= i__1; ++l) {
#line 313 "ztrsyl.f"
	    i__2 = *m;
#line 313 "ztrsyl.f"
	    for (k = 1; k <= i__2; ++k) {

#line 315 "ztrsyl.f"
		i__3 = k - 1;
#line 315 "ztrsyl.f"
		zdotc_(&z__1, &i__3, &a[k * a_dim1 + 1], &c__1, &c__[l * 
			c_dim1 + 1], &c__1);
#line 315 "ztrsyl.f"
		suml.r = z__1.r, suml.i = z__1.i;
#line 316 "ztrsyl.f"
		i__3 = l - 1;
#line 316 "ztrsyl.f"
		zdotu_(&z__1, &i__3, &c__[k + c_dim1], ldc, &b[l * b_dim1 + 1]
			, &c__1);
#line 316 "ztrsyl.f"
		sumr.r = z__1.r, sumr.i = z__1.i;
#line 317 "ztrsyl.f"
		i__3 = k + l * c_dim1;
#line 317 "ztrsyl.f"
		z__3.r = sgn * sumr.r, z__3.i = sgn * sumr.i;
#line 317 "ztrsyl.f"
		z__2.r = suml.r + z__3.r, z__2.i = suml.i + z__3.i;
#line 317 "ztrsyl.f"
		z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 317 "ztrsyl.f"
		vec.r = z__1.r, vec.i = z__1.i;

#line 319 "ztrsyl.f"
		scaloc = 1.;
#line 320 "ztrsyl.f"
		d_cnjg(&z__2, &a[k + k * a_dim1]);
#line 320 "ztrsyl.f"
		i__3 = l + l * b_dim1;
#line 320 "ztrsyl.f"
		z__3.r = sgn * b[i__3].r, z__3.i = sgn * b[i__3].i;
#line 320 "ztrsyl.f"
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 320 "ztrsyl.f"
		a11.r = z__1.r, a11.i = z__1.i;
#line 321 "ztrsyl.f"
		da11 = (d__1 = a11.r, abs(d__1)) + (d__2 = d_imag(&a11), abs(
			d__2));
#line 322 "ztrsyl.f"
		if (da11 <= smin) {
#line 323 "ztrsyl.f"
		    a11.r = smin, a11.i = 0.;
#line 324 "ztrsyl.f"
		    da11 = smin;
#line 325 "ztrsyl.f"
		    *info = 1;
#line 326 "ztrsyl.f"
		}
#line 327 "ztrsyl.f"
		db = (d__1 = vec.r, abs(d__1)) + (d__2 = d_imag(&vec), abs(
			d__2));
#line 328 "ztrsyl.f"
		if (da11 < 1. && db > 1.) {
#line 329 "ztrsyl.f"
		    if (db > bignum * da11) {
#line 329 "ztrsyl.f"
			scaloc = 1. / db;
#line 329 "ztrsyl.f"
		    }
#line 331 "ztrsyl.f"
		}

#line 333 "ztrsyl.f"
		z__3.r = scaloc, z__3.i = 0.;
#line 333 "ztrsyl.f"
		z__2.r = vec.r * z__3.r - vec.i * z__3.i, z__2.i = vec.r * 
			z__3.i + vec.i * z__3.r;
#line 333 "ztrsyl.f"
		zladiv_(&z__1, &z__2, &a11);
#line 333 "ztrsyl.f"
		x11.r = z__1.r, x11.i = z__1.i;

#line 335 "ztrsyl.f"
		if (scaloc != 1.) {
#line 336 "ztrsyl.f"
		    i__3 = *n;
#line 336 "ztrsyl.f"
		    for (j = 1; j <= i__3; ++j) {
#line 337 "ztrsyl.f"
			zdscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 338 "ztrsyl.f"
/* L40: */
#line 338 "ztrsyl.f"
		    }
#line 339 "ztrsyl.f"
		    *scale *= scaloc;
#line 340 "ztrsyl.f"
		}
#line 341 "ztrsyl.f"
		i__3 = k + l * c_dim1;
#line 341 "ztrsyl.f"
		c__[i__3].r = x11.r, c__[i__3].i = x11.i;

#line 343 "ztrsyl.f"
/* L50: */
#line 343 "ztrsyl.f"
	    }
#line 344 "ztrsyl.f"
/* L60: */
#line 344 "ztrsyl.f"
	}

#line 346 "ztrsyl.f"
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

#line 363 "ztrsyl.f"
	for (l = *n; l >= 1; --l) {
#line 364 "ztrsyl.f"
	    i__1 = *m;
#line 364 "ztrsyl.f"
	    for (k = 1; k <= i__1; ++k) {

#line 366 "ztrsyl.f"
		i__2 = k - 1;
#line 366 "ztrsyl.f"
		zdotc_(&z__1, &i__2, &a[k * a_dim1 + 1], &c__1, &c__[l * 
			c_dim1 + 1], &c__1);
#line 366 "ztrsyl.f"
		suml.r = z__1.r, suml.i = z__1.i;
#line 367 "ztrsyl.f"
		i__2 = *n - l;
/* Computing MIN */
#line 367 "ztrsyl.f"
		i__3 = l + 1;
/* Computing MIN */
#line 367 "ztrsyl.f"
		i__4 = l + 1;
#line 367 "ztrsyl.f"
		zdotc_(&z__1, &i__2, &c__[k + min(i__3,*n) * c_dim1], ldc, &b[
			l + min(i__4,*n) * b_dim1], ldb);
#line 367 "ztrsyl.f"
		sumr.r = z__1.r, sumr.i = z__1.i;
#line 369 "ztrsyl.f"
		i__2 = k + l * c_dim1;
#line 369 "ztrsyl.f"
		d_cnjg(&z__4, &sumr);
#line 369 "ztrsyl.f"
		z__3.r = sgn * z__4.r, z__3.i = sgn * z__4.i;
#line 369 "ztrsyl.f"
		z__2.r = suml.r + z__3.r, z__2.i = suml.i + z__3.i;
#line 369 "ztrsyl.f"
		z__1.r = c__[i__2].r - z__2.r, z__1.i = c__[i__2].i - z__2.i;
#line 369 "ztrsyl.f"
		vec.r = z__1.r, vec.i = z__1.i;

#line 371 "ztrsyl.f"
		scaloc = 1.;
#line 372 "ztrsyl.f"
		i__2 = k + k * a_dim1;
#line 372 "ztrsyl.f"
		i__3 = l + l * b_dim1;
#line 372 "ztrsyl.f"
		z__3.r = sgn * b[i__3].r, z__3.i = sgn * b[i__3].i;
#line 372 "ztrsyl.f"
		z__2.r = a[i__2].r + z__3.r, z__2.i = a[i__2].i + z__3.i;
#line 372 "ztrsyl.f"
		d_cnjg(&z__1, &z__2);
#line 372 "ztrsyl.f"
		a11.r = z__1.r, a11.i = z__1.i;
#line 373 "ztrsyl.f"
		da11 = (d__1 = a11.r, abs(d__1)) + (d__2 = d_imag(&a11), abs(
			d__2));
#line 374 "ztrsyl.f"
		if (da11 <= smin) {
#line 375 "ztrsyl.f"
		    a11.r = smin, a11.i = 0.;
#line 376 "ztrsyl.f"
		    da11 = smin;
#line 377 "ztrsyl.f"
		    *info = 1;
#line 378 "ztrsyl.f"
		}
#line 379 "ztrsyl.f"
		db = (d__1 = vec.r, abs(d__1)) + (d__2 = d_imag(&vec), abs(
			d__2));
#line 380 "ztrsyl.f"
		if (da11 < 1. && db > 1.) {
#line 381 "ztrsyl.f"
		    if (db > bignum * da11) {
#line 381 "ztrsyl.f"
			scaloc = 1. / db;
#line 381 "ztrsyl.f"
		    }
#line 383 "ztrsyl.f"
		}

#line 385 "ztrsyl.f"
		z__3.r = scaloc, z__3.i = 0.;
#line 385 "ztrsyl.f"
		z__2.r = vec.r * z__3.r - vec.i * z__3.i, z__2.i = vec.r * 
			z__3.i + vec.i * z__3.r;
#line 385 "ztrsyl.f"
		zladiv_(&z__1, &z__2, &a11);
#line 385 "ztrsyl.f"
		x11.r = z__1.r, x11.i = z__1.i;

#line 387 "ztrsyl.f"
		if (scaloc != 1.) {
#line 388 "ztrsyl.f"
		    i__2 = *n;
#line 388 "ztrsyl.f"
		    for (j = 1; j <= i__2; ++j) {
#line 389 "ztrsyl.f"
			zdscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 390 "ztrsyl.f"
/* L70: */
#line 390 "ztrsyl.f"
		    }
#line 391 "ztrsyl.f"
		    *scale *= scaloc;
#line 392 "ztrsyl.f"
		}
#line 393 "ztrsyl.f"
		i__2 = k + l * c_dim1;
#line 393 "ztrsyl.f"
		c__[i__2].r = x11.r, c__[i__2].i = x11.i;

#line 395 "ztrsyl.f"
/* L80: */
#line 395 "ztrsyl.f"
	    }
#line 396 "ztrsyl.f"
/* L90: */
#line 396 "ztrsyl.f"
	}

#line 398 "ztrsyl.f"
    } else if (notrna && ! notrnb) {

/*        Solve    A*X + ISGN*X*B**H = C. */

/*        The (K,L)th block of X is determined starting from */
/*        bottom-left corner column by column by */

/*           A(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L) */

/*        Where */
/*                    M                          N */
/*          R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B**H(L,J)] */
/*                  I=K+1                      J=L+1 */

#line 412 "ztrsyl.f"
	for (l = *n; l >= 1; --l) {
#line 413 "ztrsyl.f"
	    for (k = *m; k >= 1; --k) {

#line 415 "ztrsyl.f"
		i__1 = *m - k;
/* Computing MIN */
#line 415 "ztrsyl.f"
		i__2 = k + 1;
/* Computing MIN */
#line 415 "ztrsyl.f"
		i__3 = k + 1;
#line 415 "ztrsyl.f"
		zdotu_(&z__1, &i__1, &a[k + min(i__2,*m) * a_dim1], lda, &c__[
			min(i__3,*m) + l * c_dim1], &c__1);
#line 415 "ztrsyl.f"
		suml.r = z__1.r, suml.i = z__1.i;
#line 417 "ztrsyl.f"
		i__1 = *n - l;
/* Computing MIN */
#line 417 "ztrsyl.f"
		i__2 = l + 1;
/* Computing MIN */
#line 417 "ztrsyl.f"
		i__3 = l + 1;
#line 417 "ztrsyl.f"
		zdotc_(&z__1, &i__1, &c__[k + min(i__2,*n) * c_dim1], ldc, &b[
			l + min(i__3,*n) * b_dim1], ldb);
#line 417 "ztrsyl.f"
		sumr.r = z__1.r, sumr.i = z__1.i;
#line 419 "ztrsyl.f"
		i__1 = k + l * c_dim1;
#line 419 "ztrsyl.f"
		d_cnjg(&z__4, &sumr);
#line 419 "ztrsyl.f"
		z__3.r = sgn * z__4.r, z__3.i = sgn * z__4.i;
#line 419 "ztrsyl.f"
		z__2.r = suml.r + z__3.r, z__2.i = suml.i + z__3.i;
#line 419 "ztrsyl.f"
		z__1.r = c__[i__1].r - z__2.r, z__1.i = c__[i__1].i - z__2.i;
#line 419 "ztrsyl.f"
		vec.r = z__1.r, vec.i = z__1.i;

#line 421 "ztrsyl.f"
		scaloc = 1.;
#line 422 "ztrsyl.f"
		i__1 = k + k * a_dim1;
#line 422 "ztrsyl.f"
		d_cnjg(&z__3, &b[l + l * b_dim1]);
#line 422 "ztrsyl.f"
		z__2.r = sgn * z__3.r, z__2.i = sgn * z__3.i;
#line 422 "ztrsyl.f"
		z__1.r = a[i__1].r + z__2.r, z__1.i = a[i__1].i + z__2.i;
#line 422 "ztrsyl.f"
		a11.r = z__1.r, a11.i = z__1.i;
#line 423 "ztrsyl.f"
		da11 = (d__1 = a11.r, abs(d__1)) + (d__2 = d_imag(&a11), abs(
			d__2));
#line 424 "ztrsyl.f"
		if (da11 <= smin) {
#line 425 "ztrsyl.f"
		    a11.r = smin, a11.i = 0.;
#line 426 "ztrsyl.f"
		    da11 = smin;
#line 427 "ztrsyl.f"
		    *info = 1;
#line 428 "ztrsyl.f"
		}
#line 429 "ztrsyl.f"
		db = (d__1 = vec.r, abs(d__1)) + (d__2 = d_imag(&vec), abs(
			d__2));
#line 430 "ztrsyl.f"
		if (da11 < 1. && db > 1.) {
#line 431 "ztrsyl.f"
		    if (db > bignum * da11) {
#line 431 "ztrsyl.f"
			scaloc = 1. / db;
#line 431 "ztrsyl.f"
		    }
#line 433 "ztrsyl.f"
		}

#line 435 "ztrsyl.f"
		z__3.r = scaloc, z__3.i = 0.;
#line 435 "ztrsyl.f"
		z__2.r = vec.r * z__3.r - vec.i * z__3.i, z__2.i = vec.r * 
			z__3.i + vec.i * z__3.r;
#line 435 "ztrsyl.f"
		zladiv_(&z__1, &z__2, &a11);
#line 435 "ztrsyl.f"
		x11.r = z__1.r, x11.i = z__1.i;

#line 437 "ztrsyl.f"
		if (scaloc != 1.) {
#line 438 "ztrsyl.f"
		    i__1 = *n;
#line 438 "ztrsyl.f"
		    for (j = 1; j <= i__1; ++j) {
#line 439 "ztrsyl.f"
			zdscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 440 "ztrsyl.f"
/* L100: */
#line 440 "ztrsyl.f"
		    }
#line 441 "ztrsyl.f"
		    *scale *= scaloc;
#line 442 "ztrsyl.f"
		}
#line 443 "ztrsyl.f"
		i__1 = k + l * c_dim1;
#line 443 "ztrsyl.f"
		c__[i__1].r = x11.r, c__[i__1].i = x11.i;

#line 445 "ztrsyl.f"
/* L110: */
#line 445 "ztrsyl.f"
	    }
#line 446 "ztrsyl.f"
/* L120: */
#line 446 "ztrsyl.f"
	}

#line 448 "ztrsyl.f"
    }

#line 450 "ztrsyl.f"
    return 0;

/*     End of ZTRSYL */

} /* ztrsyl_ */

