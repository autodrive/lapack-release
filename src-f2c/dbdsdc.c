#line 1 "dbdsdc.f"
/* dbdsdc.f -- translated by f2c (version 20100827).
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

#line 1 "dbdsdc.f"
/* Table of constant values */

static integer c__9 = 9;
static integer c__0 = 0;
static doublereal c_b15 = 1.;
static integer c__1 = 1;
static doublereal c_b29 = 0.;

/* > \brief \b DBDSDC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DBDSDC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dbdsdc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dbdsdc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dbdsdc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DBDSDC( UPLO, COMPQ, N, D, E, U, LDU, VT, LDVT, Q, IQ, */
/*                          WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPQ, UPLO */
/*       INTEGER            INFO, LDU, LDVT, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IQ( * ), IWORK( * ) */
/*       DOUBLE PRECISION   D( * ), E( * ), Q( * ), U( LDU, * ), */
/*      $                   VT( LDVT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DBDSDC computes the singular value decomposition (SVD) of a real */
/* > N-by-N (upper or lower) bidiagonal matrix B:  B = U * S * VT, */
/* > using a divide and conquer method, where S is a diagonal matrix */
/* > with non-negative diagonal elements (the singular values of B), and */
/* > U and VT are orthogonal matrices of left and right singular vectors, */
/* > respectively. DBDSDC can be used to compute all singular values, */
/* > and optionally, singular vectors or singular vectors in compact form. */
/* > */
/* > This code makes very mild assumptions about floating point */
/* > arithmetic. It will work on machines with a guard digit in */
/* > add/subtract, or on those binary machines without guard digits */
/* > which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2. */
/* > It could conceivably fail on hexadecimal or decimal machines */
/* > without guard digits, but we know of none.  See DLASD3 for details. */
/* > */
/* > The code currently calls DLASDQ if singular values only are desired. */
/* > However, it can be slightly modified to compute singular values */
/* > using the divide and conquer method. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  B is upper bidiagonal. */
/* >          = 'L':  B is lower bidiagonal. */
/* > \endverbatim */
/* > */
/* > \param[in] COMPQ */
/* > \verbatim */
/* >          COMPQ is CHARACTER*1 */
/* >          Specifies whether singular vectors are to be computed */
/* >          as follows: */
/* >          = 'N':  Compute singular values only; */
/* >          = 'P':  Compute singular values and compute singular */
/* >                  vectors in compact form; */
/* >          = 'I':  Compute singular values and singular vectors. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix B.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the n diagonal elements of the bidiagonal matrix B. */
/* >          On exit, if INFO=0, the singular values of B. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          On entry, the elements of E contain the offdiagonal */
/* >          elements of the bidiagonal matrix whose SVD is desired. */
/* >          On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is DOUBLE PRECISION array, dimension (LDU,N) */
/* >          If  COMPQ = 'I', then: */
/* >             On exit, if INFO = 0, U contains the left singular vectors */
/* >             of the bidiagonal matrix. */
/* >          For other values of COMPQ, U is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* >          LDU is INTEGER */
/* >          The leading dimension of the array U.  LDU >= 1. */
/* >          If singular vectors are desired, then LDU >= max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[out] VT */
/* > \verbatim */
/* >          VT is DOUBLE PRECISION array, dimension (LDVT,N) */
/* >          If  COMPQ = 'I', then: */
/* >             On exit, if INFO = 0, VT**T contains the right singular */
/* >             vectors of the bidiagonal matrix. */
/* >          For other values of COMPQ, VT is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT */
/* > \verbatim */
/* >          LDVT is INTEGER */
/* >          The leading dimension of the array VT.  LDVT >= 1. */
/* >          If singular vectors are desired, then LDVT >= max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* >          Q is DOUBLE PRECISION array, dimension (LDQ) */
/* >          If  COMPQ = 'P', then: */
/* >             On exit, if INFO = 0, Q and IQ contain the left */
/* >             and right singular vectors in a compact form, */
/* >             requiring O(N log N) space instead of 2*N**2. */
/* >             In particular, Q contains all the DOUBLE PRECISION data in */
/* >             LDQ >= N*(11 + 2*SMLSIZ + 8*INT(LOG_2(N/(SMLSIZ+1)))) */
/* >             words of memory, where SMLSIZ is returned by ILAENV and */
/* >             is equal to the maximum size of the subproblems at the */
/* >             bottom of the computation tree (usually about 25). */
/* >          For other values of COMPQ, Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] IQ */
/* > \verbatim */
/* >          IQ is INTEGER array, dimension (LDIQ) */
/* >          If  COMPQ = 'P', then: */
/* >             On exit, if INFO = 0, Q and IQ contain the left */
/* >             and right singular vectors in a compact form, */
/* >             requiring O(N log N) space instead of 2*N**2. */
/* >             In particular, IQ contains all INTEGER data in */
/* >             LDIQ >= N*(3 + 3*INT(LOG_2(N/(SMLSIZ+1)))) */
/* >             words of memory, where SMLSIZ is returned by ILAENV and */
/* >             is equal to the maximum size of the subproblems at the */
/* >             bottom of the computation tree (usually about 25). */
/* >          For other values of COMPQ, IQ is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          If COMPQ = 'N' then LWORK >= (4 * N). */
/* >          If COMPQ = 'P' then LWORK >= (6 * N). */
/* >          If COMPQ = 'I' then LWORK >= (3 * N**2 + 4 * N). */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (8*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  The algorithm failed to compute a singular value. */
/* >                The update process of divide and conquer failed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup auxOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dbdsdc_(char *uplo, char *compq, integer *n, doublereal *
	d__, doublereal *e, doublereal *u, integer *ldu, doublereal *vt, 
	integer *ldvt, doublereal *q, integer *iq, doublereal *work, integer *
	iwork, integer *info, ftnlen uplo_len, ftnlen compq_len)
{
    /* System generated locals */
    integer u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), log(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal p, r__;
    static integer z__, ic, ii, kk;
    static doublereal cs;
    static integer is, iu;
    static doublereal sn;
    static integer nm1;
    static doublereal eps;
    static integer ivt, difl, difr, ierr, perm, mlvl, sqre;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlasr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen), dcopy_(integer *, doublereal *, integer *
	    , doublereal *, integer *), dswap_(integer *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer poles, iuplo, nsize, start;
    extern /* Subroutine */ int dlasd0_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlasda_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *), dlascl_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen), dlasdq_(char *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen), dlartg_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer givcol;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    static integer icompq;
    static doublereal orgnrm;
    static integer givnum, givptr, qstart, smlsiz, wstart, smlszp;


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */
/*  Changed dimension statement in comment describing E from (N) to */
/*  (N-1).  Sven, 17 Feb 05. */
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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 256 "dbdsdc.f"
    /* Parameter adjustments */
#line 256 "dbdsdc.f"
    --d__;
#line 256 "dbdsdc.f"
    --e;
#line 256 "dbdsdc.f"
    u_dim1 = *ldu;
#line 256 "dbdsdc.f"
    u_offset = 1 + u_dim1;
#line 256 "dbdsdc.f"
    u -= u_offset;
#line 256 "dbdsdc.f"
    vt_dim1 = *ldvt;
#line 256 "dbdsdc.f"
    vt_offset = 1 + vt_dim1;
#line 256 "dbdsdc.f"
    vt -= vt_offset;
#line 256 "dbdsdc.f"
    --q;
#line 256 "dbdsdc.f"
    --iq;
#line 256 "dbdsdc.f"
    --work;
#line 256 "dbdsdc.f"
    --iwork;
#line 256 "dbdsdc.f"

#line 256 "dbdsdc.f"
    /* Function Body */
#line 256 "dbdsdc.f"
    *info = 0;

#line 258 "dbdsdc.f"
    iuplo = 0;
#line 259 "dbdsdc.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 259 "dbdsdc.f"
	iuplo = 1;
#line 259 "dbdsdc.f"
    }
#line 261 "dbdsdc.f"
    if (lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 261 "dbdsdc.f"
	iuplo = 2;
#line 261 "dbdsdc.f"
    }
#line 263 "dbdsdc.f"
    if (lsame_(compq, "N", (ftnlen)1, (ftnlen)1)) {
#line 264 "dbdsdc.f"
	icompq = 0;
#line 265 "dbdsdc.f"
    } else if (lsame_(compq, "P", (ftnlen)1, (ftnlen)1)) {
#line 266 "dbdsdc.f"
	icompq = 1;
#line 267 "dbdsdc.f"
    } else if (lsame_(compq, "I", (ftnlen)1, (ftnlen)1)) {
#line 268 "dbdsdc.f"
	icompq = 2;
#line 269 "dbdsdc.f"
    } else {
#line 270 "dbdsdc.f"
	icompq = -1;
#line 271 "dbdsdc.f"
    }
#line 272 "dbdsdc.f"
    if (iuplo == 0) {
#line 273 "dbdsdc.f"
	*info = -1;
#line 274 "dbdsdc.f"
    } else if (icompq < 0) {
#line 275 "dbdsdc.f"
	*info = -2;
#line 276 "dbdsdc.f"
    } else if (*n < 0) {
#line 277 "dbdsdc.f"
	*info = -3;
#line 278 "dbdsdc.f"
    } else if (*ldu < 1 || icompq == 2 && *ldu < *n) {
#line 280 "dbdsdc.f"
	*info = -7;
#line 281 "dbdsdc.f"
    } else if (*ldvt < 1 || icompq == 2 && *ldvt < *n) {
#line 283 "dbdsdc.f"
	*info = -9;
#line 284 "dbdsdc.f"
    }
#line 285 "dbdsdc.f"
    if (*info != 0) {
#line 286 "dbdsdc.f"
	i__1 = -(*info);
#line 286 "dbdsdc.f"
	xerbla_("DBDSDC", &i__1, (ftnlen)6);
#line 287 "dbdsdc.f"
	return 0;
#line 288 "dbdsdc.f"
    }

/*     Quick return if possible */

#line 292 "dbdsdc.f"
    if (*n == 0) {
#line 292 "dbdsdc.f"
	return 0;
#line 292 "dbdsdc.f"
    }
#line 294 "dbdsdc.f"
    smlsiz = ilaenv_(&c__9, "DBDSDC", " ", &c__0, &c__0, &c__0, &c__0, (
	    ftnlen)6, (ftnlen)1);
#line 295 "dbdsdc.f"
    if (*n == 1) {
#line 296 "dbdsdc.f"
	if (icompq == 1) {
#line 297 "dbdsdc.f"
	    q[1] = d_sign(&c_b15, &d__[1]);
#line 298 "dbdsdc.f"
	    q[smlsiz * *n + 1] = 1.;
#line 299 "dbdsdc.f"
	} else if (icompq == 2) {
#line 300 "dbdsdc.f"
	    u[u_dim1 + 1] = d_sign(&c_b15, &d__[1]);
#line 301 "dbdsdc.f"
	    vt[vt_dim1 + 1] = 1.;
#line 302 "dbdsdc.f"
	}
#line 303 "dbdsdc.f"
	d__[1] = abs(d__[1]);
#line 304 "dbdsdc.f"
	return 0;
#line 305 "dbdsdc.f"
    }
#line 306 "dbdsdc.f"
    nm1 = *n - 1;

/*     If matrix lower bidiagonal, rotate to be upper bidiagonal */
/*     by applying Givens rotations on the left */

#line 311 "dbdsdc.f"
    wstart = 1;
#line 312 "dbdsdc.f"
    qstart = 3;
#line 313 "dbdsdc.f"
    if (icompq == 1) {
#line 314 "dbdsdc.f"
	dcopy_(n, &d__[1], &c__1, &q[1], &c__1);
#line 315 "dbdsdc.f"
	i__1 = *n - 1;
#line 315 "dbdsdc.f"
	dcopy_(&i__1, &e[1], &c__1, &q[*n + 1], &c__1);
#line 316 "dbdsdc.f"
    }
#line 317 "dbdsdc.f"
    if (iuplo == 2) {
#line 318 "dbdsdc.f"
	qstart = 5;
#line 319 "dbdsdc.f"
	wstart = (*n << 1) - 1;
#line 320 "dbdsdc.f"
	i__1 = *n - 1;
#line 320 "dbdsdc.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 321 "dbdsdc.f"
	    dlartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
#line 322 "dbdsdc.f"
	    d__[i__] = r__;
#line 323 "dbdsdc.f"
	    e[i__] = sn * d__[i__ + 1];
#line 324 "dbdsdc.f"
	    d__[i__ + 1] = cs * d__[i__ + 1];
#line 325 "dbdsdc.f"
	    if (icompq == 1) {
#line 326 "dbdsdc.f"
		q[i__ + (*n << 1)] = cs;
#line 327 "dbdsdc.f"
		q[i__ + *n * 3] = sn;
#line 328 "dbdsdc.f"
	    } else if (icompq == 2) {
#line 329 "dbdsdc.f"
		work[i__] = cs;
#line 330 "dbdsdc.f"
		work[nm1 + i__] = -sn;
#line 331 "dbdsdc.f"
	    }
#line 332 "dbdsdc.f"
/* L10: */
#line 332 "dbdsdc.f"
	}
#line 333 "dbdsdc.f"
    }

/*     If ICOMPQ = 0, use DLASDQ to compute the singular values. */

#line 337 "dbdsdc.f"
    if (icompq == 0) {
#line 338 "dbdsdc.f"
	dlasdq_("U", &c__0, n, &c__0, &c__0, &c__0, &d__[1], &e[1], &vt[
		vt_offset], ldvt, &u[u_offset], ldu, &u[u_offset], ldu, &work[
		wstart], info, (ftnlen)1);
#line 340 "dbdsdc.f"
	goto L40;
#line 341 "dbdsdc.f"
    }

/*     If N is smaller than the minimum divide size SMLSIZ, then solve */
/*     the problem with another solver. */

#line 346 "dbdsdc.f"
    if (*n <= smlsiz) {
#line 347 "dbdsdc.f"
	if (icompq == 2) {
#line 348 "dbdsdc.f"
	    dlaset_("A", n, n, &c_b29, &c_b15, &u[u_offset], ldu, (ftnlen)1);
#line 349 "dbdsdc.f"
	    dlaset_("A", n, n, &c_b29, &c_b15, &vt[vt_offset], ldvt, (ftnlen)
		    1);
#line 350 "dbdsdc.f"
	    dlasdq_("U", &c__0, n, n, n, &c__0, &d__[1], &e[1], &vt[vt_offset]
		    , ldvt, &u[u_offset], ldu, &u[u_offset], ldu, &work[
		    wstart], info, (ftnlen)1);
#line 352 "dbdsdc.f"
	} else if (icompq == 1) {
#line 353 "dbdsdc.f"
	    iu = 1;
#line 354 "dbdsdc.f"
	    ivt = iu + *n;
#line 355 "dbdsdc.f"
	    dlaset_("A", n, n, &c_b29, &c_b15, &q[iu + (qstart - 1) * *n], n, 
		    (ftnlen)1);
#line 357 "dbdsdc.f"
	    dlaset_("A", n, n, &c_b29, &c_b15, &q[ivt + (qstart - 1) * *n], n,
		     (ftnlen)1);
#line 359 "dbdsdc.f"
	    dlasdq_("U", &c__0, n, n, n, &c__0, &d__[1], &e[1], &q[ivt + (
		    qstart - 1) * *n], n, &q[iu + (qstart - 1) * *n], n, &q[
		    iu + (qstart - 1) * *n], n, &work[wstart], info, (ftnlen)
		    1);
#line 364 "dbdsdc.f"
	}
#line 365 "dbdsdc.f"
	goto L40;
#line 366 "dbdsdc.f"
    }

#line 368 "dbdsdc.f"
    if (icompq == 2) {
#line 369 "dbdsdc.f"
	dlaset_("A", n, n, &c_b29, &c_b15, &u[u_offset], ldu, (ftnlen)1);
#line 370 "dbdsdc.f"
	dlaset_("A", n, n, &c_b29, &c_b15, &vt[vt_offset], ldvt, (ftnlen)1);
#line 371 "dbdsdc.f"
    }

/*     Scale. */

#line 375 "dbdsdc.f"
    orgnrm = dlanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 376 "dbdsdc.f"
    if (orgnrm == 0.) {
#line 376 "dbdsdc.f"
	return 0;
#line 376 "dbdsdc.f"
    }
#line 378 "dbdsdc.f"
    dlascl_("G", &c__0, &c__0, &orgnrm, &c_b15, n, &c__1, &d__[1], n, &ierr, (
	    ftnlen)1);
#line 379 "dbdsdc.f"
    dlascl_("G", &c__0, &c__0, &orgnrm, &c_b15, &nm1, &c__1, &e[1], &nm1, &
	    ierr, (ftnlen)1);

#line 381 "dbdsdc.f"
    eps = dlamch_("Epsilon", (ftnlen)7) * .9;

#line 383 "dbdsdc.f"
    mlvl = (integer) (log((doublereal) (*n) / (doublereal) (smlsiz + 1)) / 
	    log(2.)) + 1;
#line 384 "dbdsdc.f"
    smlszp = smlsiz + 1;

#line 386 "dbdsdc.f"
    if (icompq == 1) {
#line 387 "dbdsdc.f"
	iu = 1;
#line 388 "dbdsdc.f"
	ivt = smlsiz + 1;
#line 389 "dbdsdc.f"
	difl = ivt + smlszp;
#line 390 "dbdsdc.f"
	difr = difl + mlvl;
#line 391 "dbdsdc.f"
	z__ = difr + (mlvl << 1);
#line 392 "dbdsdc.f"
	ic = z__ + mlvl;
#line 393 "dbdsdc.f"
	is = ic + 1;
#line 394 "dbdsdc.f"
	poles = is + 1;
#line 395 "dbdsdc.f"
	givnum = poles + (mlvl << 1);

#line 397 "dbdsdc.f"
	k = 1;
#line 398 "dbdsdc.f"
	givptr = 2;
#line 399 "dbdsdc.f"
	perm = 3;
#line 400 "dbdsdc.f"
	givcol = perm + mlvl;
#line 401 "dbdsdc.f"
    }

#line 403 "dbdsdc.f"
    i__1 = *n;
#line 403 "dbdsdc.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 404 "dbdsdc.f"
	if ((d__1 = d__[i__], abs(d__1)) < eps) {
#line 405 "dbdsdc.f"
	    d__[i__] = d_sign(&eps, &d__[i__]);
#line 406 "dbdsdc.f"
	}
#line 407 "dbdsdc.f"
/* L20: */
#line 407 "dbdsdc.f"
    }

#line 409 "dbdsdc.f"
    start = 1;
#line 410 "dbdsdc.f"
    sqre = 0;

#line 412 "dbdsdc.f"
    i__1 = nm1;
#line 412 "dbdsdc.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 413 "dbdsdc.f"
	if ((d__1 = e[i__], abs(d__1)) < eps || i__ == nm1) {

/*        Subproblem found. First determine its size and then */
/*        apply divide and conquer on it. */

#line 418 "dbdsdc.f"
	    if (i__ < nm1) {

/*        A subproblem with E(I) small for I < NM1. */

#line 422 "dbdsdc.f"
		nsize = i__ - start + 1;
#line 423 "dbdsdc.f"
	    } else if ((d__1 = e[i__], abs(d__1)) >= eps) {

/*        A subproblem with E(NM1) not too small but I = NM1. */

#line 427 "dbdsdc.f"
		nsize = *n - start + 1;
#line 428 "dbdsdc.f"
	    } else {

/*        A subproblem with E(NM1) small. This implies an */
/*        1-by-1 subproblem at D(N). Solve this 1-by-1 problem */
/*        first. */

#line 434 "dbdsdc.f"
		nsize = i__ - start + 1;
#line 435 "dbdsdc.f"
		if (icompq == 2) {
#line 436 "dbdsdc.f"
		    u[*n + *n * u_dim1] = d_sign(&c_b15, &d__[*n]);
#line 437 "dbdsdc.f"
		    vt[*n + *n * vt_dim1] = 1.;
#line 438 "dbdsdc.f"
		} else if (icompq == 1) {
#line 439 "dbdsdc.f"
		    q[*n + (qstart - 1) * *n] = d_sign(&c_b15, &d__[*n]);
#line 440 "dbdsdc.f"
		    q[*n + (smlsiz + qstart - 1) * *n] = 1.;
#line 441 "dbdsdc.f"
		}
#line 442 "dbdsdc.f"
		d__[*n] = (d__1 = d__[*n], abs(d__1));
#line 443 "dbdsdc.f"
	    }
#line 444 "dbdsdc.f"
	    if (icompq == 2) {
#line 445 "dbdsdc.f"
		dlasd0_(&nsize, &sqre, &d__[start], &e[start], &u[start + 
			start * u_dim1], ldu, &vt[start + start * vt_dim1], 
			ldvt, &smlsiz, &iwork[1], &work[wstart], info);
#line 448 "dbdsdc.f"
	    } else {
#line 449 "dbdsdc.f"
		dlasda_(&icompq, &smlsiz, &nsize, &sqre, &d__[start], &e[
			start], &q[start + (iu + qstart - 2) * *n], n, &q[
			start + (ivt + qstart - 2) * *n], &iq[start + k * *n],
			 &q[start + (difl + qstart - 2) * *n], &q[start + (
			difr + qstart - 2) * *n], &q[start + (z__ + qstart - 
			2) * *n], &q[start + (poles + qstart - 2) * *n], &iq[
			start + givptr * *n], &iq[start + givcol * *n], n, &
			iq[start + perm * *n], &q[start + (givnum + qstart - 
			2) * *n], &q[start + (ic + qstart - 2) * *n], &q[
			start + (is + qstart - 2) * *n], &work[wstart], &
			iwork[1], info);
#line 462 "dbdsdc.f"
	    }
#line 463 "dbdsdc.f"
	    if (*info != 0) {
#line 464 "dbdsdc.f"
		return 0;
#line 465 "dbdsdc.f"
	    }
#line 466 "dbdsdc.f"
	    start = i__ + 1;
#line 467 "dbdsdc.f"
	}
#line 468 "dbdsdc.f"
/* L30: */
#line 468 "dbdsdc.f"
    }

/*     Unscale */

#line 472 "dbdsdc.f"
    dlascl_("G", &c__0, &c__0, &c_b15, &orgnrm, n, &c__1, &d__[1], n, &ierr, (
	    ftnlen)1);
#line 473 "dbdsdc.f"
L40:

/*     Use Selection Sort to minimize swaps of singular vectors */

#line 477 "dbdsdc.f"
    i__1 = *n;
#line 477 "dbdsdc.f"
    for (ii = 2; ii <= i__1; ++ii) {
#line 478 "dbdsdc.f"
	i__ = ii - 1;
#line 479 "dbdsdc.f"
	kk = i__;
#line 480 "dbdsdc.f"
	p = d__[i__];
#line 481 "dbdsdc.f"
	i__2 = *n;
#line 481 "dbdsdc.f"
	for (j = ii; j <= i__2; ++j) {
#line 482 "dbdsdc.f"
	    if (d__[j] > p) {
#line 483 "dbdsdc.f"
		kk = j;
#line 484 "dbdsdc.f"
		p = d__[j];
#line 485 "dbdsdc.f"
	    }
#line 486 "dbdsdc.f"
/* L50: */
#line 486 "dbdsdc.f"
	}
#line 487 "dbdsdc.f"
	if (kk != i__) {
#line 488 "dbdsdc.f"
	    d__[kk] = d__[i__];
#line 489 "dbdsdc.f"
	    d__[i__] = p;
#line 490 "dbdsdc.f"
	    if (icompq == 1) {
#line 491 "dbdsdc.f"
		iq[i__] = kk;
#line 492 "dbdsdc.f"
	    } else if (icompq == 2) {
#line 493 "dbdsdc.f"
		dswap_(n, &u[i__ * u_dim1 + 1], &c__1, &u[kk * u_dim1 + 1], &
			c__1);
#line 494 "dbdsdc.f"
		dswap_(n, &vt[i__ + vt_dim1], ldvt, &vt[kk + vt_dim1], ldvt);
#line 495 "dbdsdc.f"
	    }
#line 496 "dbdsdc.f"
	} else if (icompq == 1) {
#line 497 "dbdsdc.f"
	    iq[i__] = i__;
#line 498 "dbdsdc.f"
	}
#line 499 "dbdsdc.f"
/* L60: */
#line 499 "dbdsdc.f"
    }

/*     If ICOMPQ = 1, use IQ(N,1) as the indicator for UPLO */

#line 503 "dbdsdc.f"
    if (icompq == 1) {
#line 504 "dbdsdc.f"
	if (iuplo == 1) {
#line 505 "dbdsdc.f"
	    iq[*n] = 1;
#line 506 "dbdsdc.f"
	} else {
#line 507 "dbdsdc.f"
	    iq[*n] = 0;
#line 508 "dbdsdc.f"
	}
#line 509 "dbdsdc.f"
    }

/*     If B is lower bidiagonal, update U by those Givens rotations */
/*     which rotated B to be upper bidiagonal */

#line 514 "dbdsdc.f"
    if (iuplo == 2 && icompq == 2) {
#line 514 "dbdsdc.f"
	dlasr_("L", "V", "B", n, n, &work[1], &work[*n], &u[u_offset], ldu, (
		ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 514 "dbdsdc.f"
    }

#line 517 "dbdsdc.f"
    return 0;

/*     End of DBDSDC */

} /* dbdsdc_ */

