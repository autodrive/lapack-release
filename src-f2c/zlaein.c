#line 1 "zlaein.f"
/* zlaein.f -- translated by f2c (version 20100827).
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

#line 1 "zlaein.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLAEIN computes a specified right or left eigenvector of an upper Hessenberg matrix by inverse 
iteration. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAEIN + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaein.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaein.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaein.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAEIN( RIGHTV, NOINIT, N, H, LDH, W, V, B, LDB, RWORK, */
/*                          EPS3, SMLNUM, INFO ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            NOINIT, RIGHTV */
/*       INTEGER            INFO, LDB, LDH, N */
/*       DOUBLE PRECISION   EPS3, SMLNUM */
/*       COMPLEX*16         W */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX*16         B( LDB, * ), H( LDH, * ), V( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAEIN uses inverse iteration to find a right or left eigenvector */
/* > corresponding to the eigenvalue W of a complex upper Hessenberg */
/* > matrix H. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] RIGHTV */
/* > \verbatim */
/* >          RIGHTV is LOGICAL */
/* >          = .TRUE. : compute right eigenvector; */
/* >          = .FALSE.: compute left eigenvector. */
/* > \endverbatim */
/* > */
/* > \param[in] NOINIT */
/* > \verbatim */
/* >          NOINIT is LOGICAL */
/* >          = .TRUE. : no initial vector supplied in V */
/* >          = .FALSE.: initial vector supplied in V. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix H.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] H */
/* > \verbatim */
/* >          H is COMPLEX*16 array, dimension (LDH,N) */
/* >          The upper Hessenberg matrix H. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* >          LDH is INTEGER */
/* >          The leading dimension of the array H.  LDH >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] W */
/* > \verbatim */
/* >          W is COMPLEX*16 */
/* >          The eigenvalue of H whose corresponding right or left */
/* >          eigenvector is to be computed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] V */
/* > \verbatim */
/* >          V is COMPLEX*16 array, dimension (N) */
/* >          On entry, if NOINIT = .FALSE., V must contain a starting */
/* >          vector for inverse iteration; otherwise V need not be set. */
/* >          On exit, V contains the computed eigenvector, normalized so */
/* >          that the component of largest magnitude has magnitude 1; here */
/* >          the magnitude of a complex number (x,y) is taken to be */
/* >          |x| + |y|. */
/* > \endverbatim */
/* > */
/* > \param[out] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB,N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[in] EPS3 */
/* > \verbatim */
/* >          EPS3 is DOUBLE PRECISION */
/* >          A small machine-dependent value which is used to perturb */
/* >          close eigenvalues, and to replace zero pivots. */
/* > \endverbatim */
/* > */
/* > \param[in] SMLNUM */
/* > \verbatim */
/* >          SMLNUM is DOUBLE PRECISION */
/* >          A machine-dependent value close to the underflow threshold. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          = 1:  inverse iteration did not converge; V is set to the */
/* >                last iterate. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zlaein_(logical *rightv, logical *noinit, integer *n, 
	doublecomplex *h__, integer *ldh, doublecomplex *w, doublecomplex *v, 
	doublecomplex *b, integer *ldb, doublereal *rwork, doublereal *eps3, 
	doublereal *smlnum, integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, h_dim1, h_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublecomplex x, ei, ej;
    static integer its, ierr;
    static doublecomplex temp;
    static doublereal scale;
    static char trans[1];
    static doublereal rtemp, rootn, vnorm;
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int zdscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    extern integer izamax_(integer *, doublecomplex *, integer *);
    extern /* Double Complex */ VOID zladiv_(doublecomplex *, doublecomplex *,
	     doublecomplex *);
    static char normin[1];
    extern doublereal dzasum_(integer *, doublecomplex *, integer *);
    static doublereal nrmsml;
    extern /* Subroutine */ int zlatrs_(char *, char *, char *, char *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    static doublereal growto;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 202 "zlaein.f"
    /* Parameter adjustments */
#line 202 "zlaein.f"
    h_dim1 = *ldh;
#line 202 "zlaein.f"
    h_offset = 1 + h_dim1;
#line 202 "zlaein.f"
    h__ -= h_offset;
#line 202 "zlaein.f"
    --v;
#line 202 "zlaein.f"
    b_dim1 = *ldb;
#line 202 "zlaein.f"
    b_offset = 1 + b_dim1;
#line 202 "zlaein.f"
    b -= b_offset;
#line 202 "zlaein.f"
    --rwork;
#line 202 "zlaein.f"

#line 202 "zlaein.f"
    /* Function Body */
#line 202 "zlaein.f"
    *info = 0;

/*     GROWTO is the threshold used in the acceptance test for an */
/*     eigenvector. */

#line 207 "zlaein.f"
    rootn = sqrt((doublereal) (*n));
#line 208 "zlaein.f"
    growto = .1 / rootn;
/* Computing MAX */
#line 209 "zlaein.f"
    d__1 = 1., d__2 = *eps3 * rootn;
#line 209 "zlaein.f"
    nrmsml = max(d__1,d__2) * *smlnum;

/*     Form B = H - W*I (except that the subdiagonal elements are not */
/*     stored). */

#line 214 "zlaein.f"
    i__1 = *n;
#line 214 "zlaein.f"
    for (j = 1; j <= i__1; ++j) {
#line 215 "zlaein.f"
	i__2 = j - 1;
#line 215 "zlaein.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 216 "zlaein.f"
	    i__3 = i__ + j * b_dim1;
#line 216 "zlaein.f"
	    i__4 = i__ + j * h_dim1;
#line 216 "zlaein.f"
	    b[i__3].r = h__[i__4].r, b[i__3].i = h__[i__4].i;
#line 217 "zlaein.f"
/* L10: */
#line 217 "zlaein.f"
	}
#line 218 "zlaein.f"
	i__2 = j + j * b_dim1;
#line 218 "zlaein.f"
	i__3 = j + j * h_dim1;
#line 218 "zlaein.f"
	z__1.r = h__[i__3].r - w->r, z__1.i = h__[i__3].i - w->i;
#line 218 "zlaein.f"
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 219 "zlaein.f"
/* L20: */
#line 219 "zlaein.f"
    }

#line 221 "zlaein.f"
    if (*noinit) {

/*        Initialize V. */

#line 225 "zlaein.f"
	i__1 = *n;
#line 225 "zlaein.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 226 "zlaein.f"
	    i__2 = i__;
#line 226 "zlaein.f"
	    v[i__2].r = *eps3, v[i__2].i = 0.;
#line 227 "zlaein.f"
/* L30: */
#line 227 "zlaein.f"
	}
#line 228 "zlaein.f"
    } else {

/*        Scale supplied initial vector. */

#line 232 "zlaein.f"
	vnorm = dznrm2_(n, &v[1], &c__1);
#line 233 "zlaein.f"
	d__1 = *eps3 * rootn / max(vnorm,nrmsml);
#line 233 "zlaein.f"
	zdscal_(n, &d__1, &v[1], &c__1);
#line 234 "zlaein.f"
    }

#line 236 "zlaein.f"
    if (*rightv) {

/*        LU decomposition with partial pivoting of B, replacing zero */
/*        pivots by EPS3. */

#line 241 "zlaein.f"
	i__1 = *n - 1;
#line 241 "zlaein.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 242 "zlaein.f"
	    i__2 = i__ + 1 + i__ * h_dim1;
#line 242 "zlaein.f"
	    ei.r = h__[i__2].r, ei.i = h__[i__2].i;
#line 243 "zlaein.f"
	    i__2 = i__ + i__ * b_dim1;
#line 243 "zlaein.f"
	    if ((d__1 = b[i__2].r, abs(d__1)) + (d__2 = d_imag(&b[i__ + i__ * 
		    b_dim1]), abs(d__2)) < (d__3 = ei.r, abs(d__3)) + (d__4 = 
		    d_imag(&ei), abs(d__4))) {

/*              Interchange rows and eliminate. */

#line 247 "zlaein.f"
		zladiv_(&z__1, &b[i__ + i__ * b_dim1], &ei);
#line 247 "zlaein.f"
		x.r = z__1.r, x.i = z__1.i;
#line 248 "zlaein.f"
		i__2 = i__ + i__ * b_dim1;
#line 248 "zlaein.f"
		b[i__2].r = ei.r, b[i__2].i = ei.i;
#line 249 "zlaein.f"
		i__2 = *n;
#line 249 "zlaein.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 250 "zlaein.f"
		    i__3 = i__ + 1 + j * b_dim1;
#line 250 "zlaein.f"
		    temp.r = b[i__3].r, temp.i = b[i__3].i;
#line 251 "zlaein.f"
		    i__3 = i__ + 1 + j * b_dim1;
#line 251 "zlaein.f"
		    i__4 = i__ + j * b_dim1;
#line 251 "zlaein.f"
		    z__2.r = x.r * temp.r - x.i * temp.i, z__2.i = x.r * 
			    temp.i + x.i * temp.r;
#line 251 "zlaein.f"
		    z__1.r = b[i__4].r - z__2.r, z__1.i = b[i__4].i - z__2.i;
#line 251 "zlaein.f"
		    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 252 "zlaein.f"
		    i__3 = i__ + j * b_dim1;
#line 252 "zlaein.f"
		    b[i__3].r = temp.r, b[i__3].i = temp.i;
#line 253 "zlaein.f"
/* L40: */
#line 253 "zlaein.f"
		}
#line 254 "zlaein.f"
	    } else {

/*              Eliminate without interchange. */

#line 258 "zlaein.f"
		i__2 = i__ + i__ * b_dim1;
#line 258 "zlaein.f"
		if (b[i__2].r == 0. && b[i__2].i == 0.) {
#line 258 "zlaein.f"
		    i__3 = i__ + i__ * b_dim1;
#line 258 "zlaein.f"
		    b[i__3].r = *eps3, b[i__3].i = 0.;
#line 258 "zlaein.f"
		}
#line 260 "zlaein.f"
		zladiv_(&z__1, &ei, &b[i__ + i__ * b_dim1]);
#line 260 "zlaein.f"
		x.r = z__1.r, x.i = z__1.i;
#line 261 "zlaein.f"
		if (x.r != 0. || x.i != 0.) {
#line 262 "zlaein.f"
		    i__2 = *n;
#line 262 "zlaein.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 263 "zlaein.f"
			i__3 = i__ + 1 + j * b_dim1;
#line 263 "zlaein.f"
			i__4 = i__ + 1 + j * b_dim1;
#line 263 "zlaein.f"
			i__5 = i__ + j * b_dim1;
#line 263 "zlaein.f"
			z__2.r = x.r * b[i__5].r - x.i * b[i__5].i, z__2.i = 
				x.r * b[i__5].i + x.i * b[i__5].r;
#line 263 "zlaein.f"
			z__1.r = b[i__4].r - z__2.r, z__1.i = b[i__4].i - 
				z__2.i;
#line 263 "zlaein.f"
			b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 264 "zlaein.f"
/* L50: */
#line 264 "zlaein.f"
		    }
#line 265 "zlaein.f"
		}
#line 266 "zlaein.f"
	    }
#line 267 "zlaein.f"
/* L60: */
#line 267 "zlaein.f"
	}
#line 268 "zlaein.f"
	i__1 = *n + *n * b_dim1;
#line 268 "zlaein.f"
	if (b[i__1].r == 0. && b[i__1].i == 0.) {
#line 268 "zlaein.f"
	    i__2 = *n + *n * b_dim1;
#line 268 "zlaein.f"
	    b[i__2].r = *eps3, b[i__2].i = 0.;
#line 268 "zlaein.f"
	}

#line 271 "zlaein.f"
	*(unsigned char *)trans = 'N';

#line 273 "zlaein.f"
    } else {

/*        UL decomposition with partial pivoting of B, replacing zero */
/*        pivots by EPS3. */

#line 278 "zlaein.f"
	for (j = *n; j >= 2; --j) {
#line 279 "zlaein.f"
	    i__1 = j + (j - 1) * h_dim1;
#line 279 "zlaein.f"
	    ej.r = h__[i__1].r, ej.i = h__[i__1].i;
#line 280 "zlaein.f"
	    i__1 = j + j * b_dim1;
#line 280 "zlaein.f"
	    if ((d__1 = b[i__1].r, abs(d__1)) + (d__2 = d_imag(&b[j + j * 
		    b_dim1]), abs(d__2)) < (d__3 = ej.r, abs(d__3)) + (d__4 = 
		    d_imag(&ej), abs(d__4))) {

/*              Interchange columns and eliminate. */

#line 284 "zlaein.f"
		zladiv_(&z__1, &b[j + j * b_dim1], &ej);
#line 284 "zlaein.f"
		x.r = z__1.r, x.i = z__1.i;
#line 285 "zlaein.f"
		i__1 = j + j * b_dim1;
#line 285 "zlaein.f"
		b[i__1].r = ej.r, b[i__1].i = ej.i;
#line 286 "zlaein.f"
		i__1 = j - 1;
#line 286 "zlaein.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 287 "zlaein.f"
		    i__2 = i__ + (j - 1) * b_dim1;
#line 287 "zlaein.f"
		    temp.r = b[i__2].r, temp.i = b[i__2].i;
#line 288 "zlaein.f"
		    i__2 = i__ + (j - 1) * b_dim1;
#line 288 "zlaein.f"
		    i__3 = i__ + j * b_dim1;
#line 288 "zlaein.f"
		    z__2.r = x.r * temp.r - x.i * temp.i, z__2.i = x.r * 
			    temp.i + x.i * temp.r;
#line 288 "zlaein.f"
		    z__1.r = b[i__3].r - z__2.r, z__1.i = b[i__3].i - z__2.i;
#line 288 "zlaein.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 289 "zlaein.f"
		    i__2 = i__ + j * b_dim1;
#line 289 "zlaein.f"
		    b[i__2].r = temp.r, b[i__2].i = temp.i;
#line 290 "zlaein.f"
/* L70: */
#line 290 "zlaein.f"
		}
#line 291 "zlaein.f"
	    } else {

/*              Eliminate without interchange. */

#line 295 "zlaein.f"
		i__1 = j + j * b_dim1;
#line 295 "zlaein.f"
		if (b[i__1].r == 0. && b[i__1].i == 0.) {
#line 295 "zlaein.f"
		    i__2 = j + j * b_dim1;
#line 295 "zlaein.f"
		    b[i__2].r = *eps3, b[i__2].i = 0.;
#line 295 "zlaein.f"
		}
#line 297 "zlaein.f"
		zladiv_(&z__1, &ej, &b[j + j * b_dim1]);
#line 297 "zlaein.f"
		x.r = z__1.r, x.i = z__1.i;
#line 298 "zlaein.f"
		if (x.r != 0. || x.i != 0.) {
#line 299 "zlaein.f"
		    i__1 = j - 1;
#line 299 "zlaein.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 300 "zlaein.f"
			i__2 = i__ + (j - 1) * b_dim1;
#line 300 "zlaein.f"
			i__3 = i__ + (j - 1) * b_dim1;
#line 300 "zlaein.f"
			i__4 = i__ + j * b_dim1;
#line 300 "zlaein.f"
			z__2.r = x.r * b[i__4].r - x.i * b[i__4].i, z__2.i = 
				x.r * b[i__4].i + x.i * b[i__4].r;
#line 300 "zlaein.f"
			z__1.r = b[i__3].r - z__2.r, z__1.i = b[i__3].i - 
				z__2.i;
#line 300 "zlaein.f"
			b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 301 "zlaein.f"
/* L80: */
#line 301 "zlaein.f"
		    }
#line 302 "zlaein.f"
		}
#line 303 "zlaein.f"
	    }
#line 304 "zlaein.f"
/* L90: */
#line 304 "zlaein.f"
	}
#line 305 "zlaein.f"
	i__1 = b_dim1 + 1;
#line 305 "zlaein.f"
	if (b[i__1].r == 0. && b[i__1].i == 0.) {
#line 305 "zlaein.f"
	    i__2 = b_dim1 + 1;
#line 305 "zlaein.f"
	    b[i__2].r = *eps3, b[i__2].i = 0.;
#line 305 "zlaein.f"
	}

#line 308 "zlaein.f"
	*(unsigned char *)trans = 'C';

#line 310 "zlaein.f"
    }

#line 312 "zlaein.f"
    *(unsigned char *)normin = 'N';
#line 313 "zlaein.f"
    i__1 = *n;
#line 313 "zlaein.f"
    for (its = 1; its <= i__1; ++its) {

/*        Solve U*x = scale*v for a right eigenvector */
/*          or U**H *x = scale*v for a left eigenvector, */
/*        overwriting x on v. */

#line 319 "zlaein.f"
	zlatrs_("Upper", trans, "Nonunit", normin, n, &b[b_offset], ldb, &v[1]
		, &scale, &rwork[1], &ierr, (ftnlen)5, (ftnlen)1, (ftnlen)7, (
		ftnlen)1);
#line 321 "zlaein.f"
	*(unsigned char *)normin = 'Y';

/*        Test for sufficient growth in the norm of v. */

#line 325 "zlaein.f"
	vnorm = dzasum_(n, &v[1], &c__1);
#line 326 "zlaein.f"
	if (vnorm >= growto * scale) {
#line 326 "zlaein.f"
	    goto L120;
#line 326 "zlaein.f"
	}

/*        Choose new orthogonal starting vector and try again. */

#line 331 "zlaein.f"
	rtemp = *eps3 / (rootn + 1.);
#line 332 "zlaein.f"
	v[1].r = *eps3, v[1].i = 0.;
#line 333 "zlaein.f"
	i__2 = *n;
#line 333 "zlaein.f"
	for (i__ = 2; i__ <= i__2; ++i__) {
#line 334 "zlaein.f"
	    i__3 = i__;
#line 334 "zlaein.f"
	    v[i__3].r = rtemp, v[i__3].i = 0.;
#line 335 "zlaein.f"
/* L100: */
#line 335 "zlaein.f"
	}
#line 336 "zlaein.f"
	i__2 = *n - its + 1;
#line 336 "zlaein.f"
	i__3 = *n - its + 1;
#line 336 "zlaein.f"
	d__1 = *eps3 * rootn;
#line 336 "zlaein.f"
	z__1.r = v[i__3].r - d__1, z__1.i = v[i__3].i;
#line 336 "zlaein.f"
	v[i__2].r = z__1.r, v[i__2].i = z__1.i;
#line 337 "zlaein.f"
/* L110: */
#line 337 "zlaein.f"
    }

/*     Failure to find eigenvector in N iterations. */

#line 341 "zlaein.f"
    *info = 1;

#line 343 "zlaein.f"
L120:

/*     Normalize eigenvector. */

#line 347 "zlaein.f"
    i__ = izamax_(n, &v[1], &c__1);
#line 348 "zlaein.f"
    i__1 = i__;
#line 348 "zlaein.f"
    d__3 = 1. / ((d__1 = v[i__1].r, abs(d__1)) + (d__2 = d_imag(&v[i__]), abs(
	    d__2)));
#line 348 "zlaein.f"
    zdscal_(n, &d__3, &v[1], &c__1);

#line 350 "zlaein.f"
    return 0;

/*     End of ZLAEIN */

} /* zlaein_ */

