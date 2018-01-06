#line 1 "dlasdq.f"
/* dlasdq.f -- translated by f2c (version 20100827).
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

#line 1 "dlasdq.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLASDQ computes the SVD of a real bidiagonal matrix with diagonal d and off-diagonal e. Used by
 sbdsdc. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASDQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasdq.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasdq.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasdq.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASDQ( UPLO, SQRE, N, NCVT, NRU, NCC, D, E, VT, LDVT, */
/*                          U, LDU, C, LDC, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU, SQRE */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   C( LDC, * ), D( * ), E( * ), U( LDU, * ), */
/*      $                   VT( LDVT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASDQ computes the singular value decomposition (SVD) of a real */
/* > (upper or lower) bidiagonal matrix with diagonal D and offdiagonal */
/* > E, accumulating the transformations if desired. Letting B denote */
/* > the input bidiagonal matrix, the algorithm computes orthogonal */
/* > matrices Q and P such that B = Q * S * P**T (P**T denotes the transpose */
/* > of P). The singular values S are overwritten on D. */
/* > */
/* > The input matrix U  is changed to U  * Q  if desired. */
/* > The input matrix VT is changed to P**T * VT if desired. */
/* > The input matrix C  is changed to Q**T * C  if desired. */
/* > */
/* > See "Computing  Small Singular Values of Bidiagonal Matrices With */
/* > Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan, */
/* > LAPACK Working Note #3, for a detailed description of the algorithm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >        On entry, UPLO specifies whether the input bidiagonal matrix */
/* >        is upper or lower bidiagonal, and whether it is square are */
/* >        not. */
/* >           UPLO = 'U' or 'u'   B is upper bidiagonal. */
/* >           UPLO = 'L' or 'l'   B is lower bidiagonal. */
/* > \endverbatim */
/* > */
/* > \param[in] SQRE */
/* > \verbatim */
/* >          SQRE is INTEGER */
/* >        = 0: then the input matrix is N-by-N. */
/* >        = 1: then the input matrix is N-by-(N+1) if UPLU = 'U' and */
/* >             (N+1)-by-N if UPLU = 'L'. */
/* > */
/* >        The bidiagonal matrix has */
/* >        N = NL + NR + 1 rows and */
/* >        M = N + SQRE >= N columns. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >        On entry, N specifies the number of rows and columns */
/* >        in the matrix. N must be at least 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NCVT */
/* > \verbatim */
/* >          NCVT is INTEGER */
/* >        On entry, NCVT specifies the number of columns of */
/* >        the matrix VT. NCVT must be at least 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRU */
/* > \verbatim */
/* >          NRU is INTEGER */
/* >        On entry, NRU specifies the number of rows of */
/* >        the matrix U. NRU must be at least 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NCC */
/* > \verbatim */
/* >          NCC is INTEGER */
/* >        On entry, NCC specifies the number of columns of */
/* >        the matrix C. NCC must be at least 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >        On entry, D contains the diagonal entries of the */
/* >        bidiagonal matrix whose SVD is desired. On normal exit, */
/* >        D contains the singular values in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array. */
/* >        dimension is (N-1) if SQRE = 0 and N if SQRE = 1. */
/* >        On entry, the entries of E contain the offdiagonal entries */
/* >        of the bidiagonal matrix whose SVD is desired. On normal */
/* >        exit, E will contain 0. If the algorithm does not converge, */
/* >        D and E will contain the diagonal and superdiagonal entries */
/* >        of a bidiagonal matrix orthogonally equivalent to the one */
/* >        given as input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VT */
/* > \verbatim */
/* >          VT is DOUBLE PRECISION array, dimension (LDVT, NCVT) */
/* >        On entry, contains a matrix which on exit has been */
/* >        premultiplied by P**T, dimension N-by-NCVT if SQRE = 0 */
/* >        and (N+1)-by-NCVT if SQRE = 1 (not referenced if NCVT=0). */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT */
/* > \verbatim */
/* >          LDVT is INTEGER */
/* >        On entry, LDVT specifies the leading dimension of VT as */
/* >        declared in the calling (sub) program. LDVT must be at */
/* >        least 1. If NCVT is nonzero LDVT must also be at least N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] U */
/* > \verbatim */
/* >          U is DOUBLE PRECISION array, dimension (LDU, N) */
/* >        On entry, contains a  matrix which on exit has been */
/* >        postmultiplied by Q, dimension NRU-by-N if SQRE = 0 */
/* >        and NRU-by-(N+1) if SQRE = 1 (not referenced if NRU=0). */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* >          LDU is INTEGER */
/* >        On entry, LDU  specifies the leading dimension of U as */
/* >        declared in the calling (sub) program. LDU must be at */
/* >        least max( 1, NRU ) . */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (LDC, NCC) */
/* >        On entry, contains an N-by-NCC matrix which on exit */
/* >        has been premultiplied by Q**T  dimension N-by-NCC if SQRE = 0 */
/* >        and (N+1)-by-NCC if SQRE = 1 (not referenced if NCC=0). */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >        On entry, LDC  specifies the leading dimension of C as */
/* >        declared in the calling (sub) program. LDC must be at */
/* >        least 1. If NCC is nonzero, LDC must also be at least N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (4*N) */
/* >        Workspace. Only referenced if one of NCVT, NRU, or NCC is */
/* >        nonzero, and if N is at least 2. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >        On exit, a value of 0 indicates a successful exit. */
/* >        If INFO < 0, argument number -INFO is illegal. */
/* >        If INFO > 0, the algorithm did not converge, and INFO */
/* >        specifies how many superdiagonals did not converge. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup OTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlasdq_(char *uplo, integer *sqre, integer *n, integer *
	ncvt, integer *nru, integer *ncc, doublereal *d__, doublereal *e, 
	doublereal *vt, integer *ldvt, doublereal *u, integer *ldu, 
	doublereal *c__, integer *ldc, doublereal *work, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer c_dim1, c_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, 
	    i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal r__, cs, sn;
    static integer np1, isub;
    static doublereal smin;
    static integer sqre1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlasr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen), dswap_(integer *, doublereal *, integer *
	    , doublereal *, integer *);
    static integer iuplo;
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), xerbla_(char *, 
	    integer *, ftnlen), dbdsqr_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen);
    static logical rotate;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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

/*     Test the input parameters. */

#line 253 "dlasdq.f"
    /* Parameter adjustments */
#line 253 "dlasdq.f"
    --d__;
#line 253 "dlasdq.f"
    --e;
#line 253 "dlasdq.f"
    vt_dim1 = *ldvt;
#line 253 "dlasdq.f"
    vt_offset = 1 + vt_dim1;
#line 253 "dlasdq.f"
    vt -= vt_offset;
#line 253 "dlasdq.f"
    u_dim1 = *ldu;
#line 253 "dlasdq.f"
    u_offset = 1 + u_dim1;
#line 253 "dlasdq.f"
    u -= u_offset;
#line 253 "dlasdq.f"
    c_dim1 = *ldc;
#line 253 "dlasdq.f"
    c_offset = 1 + c_dim1;
#line 253 "dlasdq.f"
    c__ -= c_offset;
#line 253 "dlasdq.f"
    --work;
#line 253 "dlasdq.f"

#line 253 "dlasdq.f"
    /* Function Body */
#line 253 "dlasdq.f"
    *info = 0;
#line 254 "dlasdq.f"
    iuplo = 0;
#line 255 "dlasdq.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 255 "dlasdq.f"
	iuplo = 1;
#line 255 "dlasdq.f"
    }
#line 257 "dlasdq.f"
    if (lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 257 "dlasdq.f"
	iuplo = 2;
#line 257 "dlasdq.f"
    }
#line 259 "dlasdq.f"
    if (iuplo == 0) {
#line 260 "dlasdq.f"
	*info = -1;
#line 261 "dlasdq.f"
    } else if (*sqre < 0 || *sqre > 1) {
#line 262 "dlasdq.f"
	*info = -2;
#line 263 "dlasdq.f"
    } else if (*n < 0) {
#line 264 "dlasdq.f"
	*info = -3;
#line 265 "dlasdq.f"
    } else if (*ncvt < 0) {
#line 266 "dlasdq.f"
	*info = -4;
#line 267 "dlasdq.f"
    } else if (*nru < 0) {
#line 268 "dlasdq.f"
	*info = -5;
#line 269 "dlasdq.f"
    } else if (*ncc < 0) {
#line 270 "dlasdq.f"
	*info = -6;
#line 271 "dlasdq.f"
    } else if (*ncvt == 0 && *ldvt < 1 || *ncvt > 0 && *ldvt < max(1,*n)) {
#line 273 "dlasdq.f"
	*info = -10;
#line 274 "dlasdq.f"
    } else if (*ldu < max(1,*nru)) {
#line 275 "dlasdq.f"
	*info = -12;
#line 276 "dlasdq.f"
    } else if (*ncc == 0 && *ldc < 1 || *ncc > 0 && *ldc < max(1,*n)) {
#line 278 "dlasdq.f"
	*info = -14;
#line 279 "dlasdq.f"
    }
#line 280 "dlasdq.f"
    if (*info != 0) {
#line 281 "dlasdq.f"
	i__1 = -(*info);
#line 281 "dlasdq.f"
	xerbla_("DLASDQ", &i__1, (ftnlen)6);
#line 282 "dlasdq.f"
	return 0;
#line 283 "dlasdq.f"
    }
#line 284 "dlasdq.f"
    if (*n == 0) {
#line 284 "dlasdq.f"
	return 0;
#line 284 "dlasdq.f"
    }

/*     ROTATE is true if any singular vectors desired, false otherwise */

#line 289 "dlasdq.f"
    rotate = *ncvt > 0 || *nru > 0 || *ncc > 0;
#line 290 "dlasdq.f"
    np1 = *n + 1;
#line 291 "dlasdq.f"
    sqre1 = *sqre;

/*     If matrix non-square upper bidiagonal, rotate to be lower */
/*     bidiagonal.  The rotations are on the right. */

#line 296 "dlasdq.f"
    if (iuplo == 1 && sqre1 == 1) {
#line 297 "dlasdq.f"
	i__1 = *n - 1;
#line 297 "dlasdq.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 298 "dlasdq.f"
	    dlartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
#line 299 "dlasdq.f"
	    d__[i__] = r__;
#line 300 "dlasdq.f"
	    e[i__] = sn * d__[i__ + 1];
#line 301 "dlasdq.f"
	    d__[i__ + 1] = cs * d__[i__ + 1];
#line 302 "dlasdq.f"
	    if (rotate) {
#line 303 "dlasdq.f"
		work[i__] = cs;
#line 304 "dlasdq.f"
		work[*n + i__] = sn;
#line 305 "dlasdq.f"
	    }
#line 306 "dlasdq.f"
/* L10: */
#line 306 "dlasdq.f"
	}
#line 307 "dlasdq.f"
	dlartg_(&d__[*n], &e[*n], &cs, &sn, &r__);
#line 308 "dlasdq.f"
	d__[*n] = r__;
#line 309 "dlasdq.f"
	e[*n] = 0.;
#line 310 "dlasdq.f"
	if (rotate) {
#line 311 "dlasdq.f"
	    work[*n] = cs;
#line 312 "dlasdq.f"
	    work[*n + *n] = sn;
#line 313 "dlasdq.f"
	}
#line 314 "dlasdq.f"
	iuplo = 2;
#line 315 "dlasdq.f"
	sqre1 = 0;

/*        Update singular vectors if desired. */

#line 319 "dlasdq.f"
	if (*ncvt > 0) {
#line 319 "dlasdq.f"
	    dlasr_("L", "V", "F", &np1, ncvt, &work[1], &work[np1], &vt[
		    vt_offset], ldvt, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 319 "dlasdq.f"
	}
#line 322 "dlasdq.f"
    }

/*     If matrix lower bidiagonal, rotate to be upper bidiagonal */
/*     by applying Givens rotations on the left. */

#line 327 "dlasdq.f"
    if (iuplo == 2) {
#line 328 "dlasdq.f"
	i__1 = *n - 1;
#line 328 "dlasdq.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 329 "dlasdq.f"
	    dlartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
#line 330 "dlasdq.f"
	    d__[i__] = r__;
#line 331 "dlasdq.f"
	    e[i__] = sn * d__[i__ + 1];
#line 332 "dlasdq.f"
	    d__[i__ + 1] = cs * d__[i__ + 1];
#line 333 "dlasdq.f"
	    if (rotate) {
#line 334 "dlasdq.f"
		work[i__] = cs;
#line 335 "dlasdq.f"
		work[*n + i__] = sn;
#line 336 "dlasdq.f"
	    }
#line 337 "dlasdq.f"
/* L20: */
#line 337 "dlasdq.f"
	}

/*        If matrix (N+1)-by-N lower bidiagonal, one additional */
/*        rotation is needed. */

#line 342 "dlasdq.f"
	if (sqre1 == 1) {
#line 343 "dlasdq.f"
	    dlartg_(&d__[*n], &e[*n], &cs, &sn, &r__);
#line 344 "dlasdq.f"
	    d__[*n] = r__;
#line 345 "dlasdq.f"
	    if (rotate) {
#line 346 "dlasdq.f"
		work[*n] = cs;
#line 347 "dlasdq.f"
		work[*n + *n] = sn;
#line 348 "dlasdq.f"
	    }
#line 349 "dlasdq.f"
	}

/*        Update singular vectors if desired. */

#line 353 "dlasdq.f"
	if (*nru > 0) {
#line 354 "dlasdq.f"
	    if (sqre1 == 0) {
#line 355 "dlasdq.f"
		dlasr_("R", "V", "F", nru, n, &work[1], &work[np1], &u[
			u_offset], ldu, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 357 "dlasdq.f"
	    } else {
#line 358 "dlasdq.f"
		dlasr_("R", "V", "F", nru, &np1, &work[1], &work[np1], &u[
			u_offset], ldu, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 360 "dlasdq.f"
	    }
#line 361 "dlasdq.f"
	}
#line 362 "dlasdq.f"
	if (*ncc > 0) {
#line 363 "dlasdq.f"
	    if (sqre1 == 0) {
#line 364 "dlasdq.f"
		dlasr_("L", "V", "F", n, ncc, &work[1], &work[np1], &c__[
			c_offset], ldc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 366 "dlasdq.f"
	    } else {
#line 367 "dlasdq.f"
		dlasr_("L", "V", "F", &np1, ncc, &work[1], &work[np1], &c__[
			c_offset], ldc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 369 "dlasdq.f"
	    }
#line 370 "dlasdq.f"
	}
#line 371 "dlasdq.f"
    }

/*     Call DBDSQR to compute the SVD of the reduced real */
/*     N-by-N upper bidiagonal matrix. */

#line 376 "dlasdq.f"
    dbdsqr_("U", n, ncvt, nru, ncc, &d__[1], &e[1], &vt[vt_offset], ldvt, &u[
	    u_offset], ldu, &c__[c_offset], ldc, &work[1], info, (ftnlen)1);

/*     Sort the singular values into ascending order (insertion sort on */
/*     singular values, but only one transposition per singular vector) */

#line 382 "dlasdq.f"
    i__1 = *n;
#line 382 "dlasdq.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Scan for smallest D(I). */

#line 386 "dlasdq.f"
	isub = i__;
#line 387 "dlasdq.f"
	smin = d__[i__];
#line 388 "dlasdq.f"
	i__2 = *n;
#line 388 "dlasdq.f"
	for (j = i__ + 1; j <= i__2; ++j) {
#line 389 "dlasdq.f"
	    if (d__[j] < smin) {
#line 390 "dlasdq.f"
		isub = j;
#line 391 "dlasdq.f"
		smin = d__[j];
#line 392 "dlasdq.f"
	    }
#line 393 "dlasdq.f"
/* L30: */
#line 393 "dlasdq.f"
	}
#line 394 "dlasdq.f"
	if (isub != i__) {

/*           Swap singular values and vectors. */

#line 398 "dlasdq.f"
	    d__[isub] = d__[i__];
#line 399 "dlasdq.f"
	    d__[i__] = smin;
#line 400 "dlasdq.f"
	    if (*ncvt > 0) {
#line 400 "dlasdq.f"
		dswap_(ncvt, &vt[isub + vt_dim1], ldvt, &vt[i__ + vt_dim1], 
			ldvt);
#line 400 "dlasdq.f"
	    }
#line 402 "dlasdq.f"
	    if (*nru > 0) {
#line 402 "dlasdq.f"
		dswap_(nru, &u[isub * u_dim1 + 1], &c__1, &u[i__ * u_dim1 + 1]
			, &c__1);
#line 402 "dlasdq.f"
	    }
#line 404 "dlasdq.f"
	    if (*ncc > 0) {
#line 404 "dlasdq.f"
		dswap_(ncc, &c__[isub + c_dim1], ldc, &c__[i__ + c_dim1], ldc)
			;
#line 404 "dlasdq.f"
	    }
#line 406 "dlasdq.f"
	}
#line 407 "dlasdq.f"
/* L40: */
#line 407 "dlasdq.f"
    }

#line 409 "dlasdq.f"
    return 0;

/*     End of DLASDQ */

} /* dlasdq_ */

