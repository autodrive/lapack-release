#line 1 "dlaexc.f"
/* dlaexc.f -- translated by f2c (version 20100827).
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

#line 1 "dlaexc.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__4 = 4;
static logical c_false = FALSE_;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__3 = 3;

/* > \brief \b DLAEXC swaps adjacent diagonal blocks of a real upper quasi-triangular matrix in Schur canonica
l form, by an orthogonal similarity transformation. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAEXC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaexc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaexc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaexc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAEXC( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            WANTQ */
/*       INTEGER            INFO, J1, LDQ, LDT, N, N1, N2 */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAEXC swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in */
/* > an upper quasi-triangular matrix T by an orthogonal similarity */
/* > transformation. */
/* > */
/* > T must be in Schur canonical form, that is, block upper triangular */
/* > with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block */
/* > has its diagonal elemnts equal and its off-diagonal elements of */
/* > opposite sign. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] WANTQ */
/* > \verbatim */
/* >          WANTQ is LOGICAL */
/* >          = .TRUE. : accumulate the transformation in the matrix Q; */
/* >          = .FALSE.: do not accumulate the transformation. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix T. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] T */
/* > \verbatim */
/* >          T is DOUBLE PRECISION array, dimension (LDT,N) */
/* >          On entry, the upper quasi-triangular matrix T, in Schur */
/* >          canonical form. */
/* >          On exit, the updated matrix T, again in Schur canonical form. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T. LDT >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is DOUBLE PRECISION array, dimension (LDQ,N) */
/* >          On entry, if WANTQ is .TRUE., the orthogonal matrix Q. */
/* >          On exit, if WANTQ is .TRUE., the updated matrix Q. */
/* >          If WANTQ is .FALSE., Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q. */
/* >          LDQ >= 1; and if WANTQ is .TRUE., LDQ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] J1 */
/* > \verbatim */
/* >          J1 is INTEGER */
/* >          The index of the first row of the first block T11. */
/* > \endverbatim */
/* > */
/* > \param[in] N1 */
/* > \verbatim */
/* >          N1 is INTEGER */
/* >          The order of the first block T11. N1 = 0, 1 or 2. */
/* > \endverbatim */
/* > */
/* > \param[in] N2 */
/* > \verbatim */
/* >          N2 is INTEGER */
/* >          The order of the second block T22. N2 = 0, 1 or 2. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          = 1: the transformed matrix T would be too far from Schur */
/* >               form; the blocks are not swapped and T and Q are */
/* >               unchanged. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlaexc_(logical *wantq, integer *n, doublereal *t, 
	integer *ldt, doublereal *q, integer *ldq, integer *j1, integer *n1, 
	integer *n2, doublereal *work, integer *info)
{
    /* System generated locals */
    integer q_dim1, q_offset, t_dim1, t_offset, i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static doublereal d__[16]	/* was [4][4] */;
    static integer k;
    static doublereal u[3], x[4]	/* was [2][2] */;
    static integer j2, j3, j4;
    static doublereal u1[3], u2[3];
    static integer nd;
    static doublereal cs, t11, t22, t33, sn, wi1, wi2, wr1, wr2, eps, tau, 
	    tau1, tau2;
    static integer ierr;
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static doublereal scale, dnorm, xnorm;
    extern /* Subroutine */ int dlanv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), dlasy2_(
	    logical *, logical *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), dlarfx_(char *, integer *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, ftnlen);
    static doublereal thresh, smlnum;


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
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 187 "dlaexc.f"
    /* Parameter adjustments */
#line 187 "dlaexc.f"
    t_dim1 = *ldt;
#line 187 "dlaexc.f"
    t_offset = 1 + t_dim1;
#line 187 "dlaexc.f"
    t -= t_offset;
#line 187 "dlaexc.f"
    q_dim1 = *ldq;
#line 187 "dlaexc.f"
    q_offset = 1 + q_dim1;
#line 187 "dlaexc.f"
    q -= q_offset;
#line 187 "dlaexc.f"
    --work;
#line 187 "dlaexc.f"

#line 187 "dlaexc.f"
    /* Function Body */
#line 187 "dlaexc.f"
    *info = 0;

/*     Quick return if possible */

#line 191 "dlaexc.f"
    if (*n == 0 || *n1 == 0 || *n2 == 0) {
#line 191 "dlaexc.f"
	return 0;
#line 191 "dlaexc.f"
    }
#line 193 "dlaexc.f"
    if (*j1 + *n1 > *n) {
#line 193 "dlaexc.f"
	return 0;
#line 193 "dlaexc.f"
    }

#line 196 "dlaexc.f"
    j2 = *j1 + 1;
#line 197 "dlaexc.f"
    j3 = *j1 + 2;
#line 198 "dlaexc.f"
    j4 = *j1 + 3;

#line 200 "dlaexc.f"
    if (*n1 == 1 && *n2 == 1) {

/*        Swap two 1-by-1 blocks. */

#line 204 "dlaexc.f"
	t11 = t[*j1 + *j1 * t_dim1];
#line 205 "dlaexc.f"
	t22 = t[j2 + j2 * t_dim1];

/*        Determine the transformation to perform the interchange. */

#line 209 "dlaexc.f"
	d__1 = t22 - t11;
#line 209 "dlaexc.f"
	dlartg_(&t[*j1 + j2 * t_dim1], &d__1, &cs, &sn, &temp);

/*        Apply transformation to the matrix T. */

#line 213 "dlaexc.f"
	if (j3 <= *n) {
#line 213 "dlaexc.f"
	    i__1 = *n - *j1 - 1;
#line 213 "dlaexc.f"
	    drot_(&i__1, &t[*j1 + j3 * t_dim1], ldt, &t[j2 + j3 * t_dim1], 
		    ldt, &cs, &sn);
#line 213 "dlaexc.f"
	}
#line 216 "dlaexc.f"
	i__1 = *j1 - 1;
#line 216 "dlaexc.f"
	drot_(&i__1, &t[*j1 * t_dim1 + 1], &c__1, &t[j2 * t_dim1 + 1], &c__1, 
		&cs, &sn);

#line 218 "dlaexc.f"
	t[*j1 + *j1 * t_dim1] = t22;
#line 219 "dlaexc.f"
	t[j2 + j2 * t_dim1] = t11;

#line 221 "dlaexc.f"
	if (*wantq) {

/*           Accumulate transformation in the matrix Q. */

#line 225 "dlaexc.f"
	    drot_(n, &q[*j1 * q_dim1 + 1], &c__1, &q[j2 * q_dim1 + 1], &c__1, 
		    &cs, &sn);
#line 226 "dlaexc.f"
	}

#line 228 "dlaexc.f"
    } else {

/*        Swapping involves at least one 2-by-2 block. */

/*        Copy the diagonal block of order N1+N2 to the local array D */
/*        and compute its norm. */

#line 235 "dlaexc.f"
	nd = *n1 + *n2;
#line 236 "dlaexc.f"
	dlacpy_("Full", &nd, &nd, &t[*j1 + *j1 * t_dim1], ldt, d__, &c__4, (
		ftnlen)4);
#line 237 "dlaexc.f"
	dnorm = dlange_("Max", &nd, &nd, d__, &c__4, &work[1], (ftnlen)3);

/*        Compute machine-dependent threshold for test for accepting */
/*        swap. */

#line 242 "dlaexc.f"
	eps = dlamch_("P", (ftnlen)1);
#line 243 "dlaexc.f"
	smlnum = dlamch_("S", (ftnlen)1) / eps;
/* Computing MAX */
#line 244 "dlaexc.f"
	d__1 = eps * 10. * dnorm;
#line 244 "dlaexc.f"
	thresh = max(d__1,smlnum);

/*        Solve T11*X - X*T22 = scale*T12 for X. */

#line 248 "dlaexc.f"
	dlasy2_(&c_false, &c_false, &c_n1, n1, n2, d__, &c__4, &d__[*n1 + 1 + 
		(*n1 + 1 << 2) - 5], &c__4, &d__[(*n1 + 1 << 2) - 4], &c__4, &
		scale, x, &c__2, &xnorm, &ierr);

/*        Swap the adjacent diagonal blocks. */

#line 254 "dlaexc.f"
	k = *n1 + *n1 + *n2 - 3;
#line 255 "dlaexc.f"
	switch (k) {
#line 255 "dlaexc.f"
	    case 1:  goto L10;
#line 255 "dlaexc.f"
	    case 2:  goto L20;
#line 255 "dlaexc.f"
	    case 3:  goto L30;
#line 255 "dlaexc.f"
	}

#line 257 "dlaexc.f"
L10:

/*        N1 = 1, N2 = 2: generate elementary reflector H so that: */

/*        ( scale, X11, X12 ) H = ( 0, 0, * ) */

#line 263 "dlaexc.f"
	u[0] = scale;
#line 264 "dlaexc.f"
	u[1] = x[0];
#line 265 "dlaexc.f"
	u[2] = x[2];
#line 266 "dlaexc.f"
	dlarfg_(&c__3, &u[2], u, &c__1, &tau);
#line 267 "dlaexc.f"
	u[2] = 1.;
#line 268 "dlaexc.f"
	t11 = t[*j1 + *j1 * t_dim1];

/*        Perform swap provisionally on diagonal block in D. */

#line 272 "dlaexc.f"
	dlarfx_("L", &c__3, &c__3, u, &tau, d__, &c__4, &work[1], (ftnlen)1);
#line 273 "dlaexc.f"
	dlarfx_("R", &c__3, &c__3, u, &tau, d__, &c__4, &work[1], (ftnlen)1);

/*        Test whether to reject swap. */

/* Computing MAX */
#line 277 "dlaexc.f"
	d__2 = abs(d__[2]), d__3 = abs(d__[6]), d__2 = max(d__2,d__3), d__3 = 
		(d__1 = d__[10] - t11, abs(d__1));
#line 277 "dlaexc.f"
	if (max(d__2,d__3) > thresh) {
#line 277 "dlaexc.f"
	    goto L50;
#line 277 "dlaexc.f"
	}

/*        Accept swap: apply transformation to the entire matrix T. */

#line 282 "dlaexc.f"
	i__1 = *n - *j1 + 1;
#line 282 "dlaexc.f"
	dlarfx_("L", &c__3, &i__1, u, &tau, &t[*j1 + *j1 * t_dim1], ldt, &
		work[1], (ftnlen)1);
#line 283 "dlaexc.f"
	dlarfx_("R", &j2, &c__3, u, &tau, &t[*j1 * t_dim1 + 1], ldt, &work[1],
		 (ftnlen)1);

#line 285 "dlaexc.f"
	t[j3 + *j1 * t_dim1] = 0.;
#line 286 "dlaexc.f"
	t[j3 + j2 * t_dim1] = 0.;
#line 287 "dlaexc.f"
	t[j3 + j3 * t_dim1] = t11;

#line 289 "dlaexc.f"
	if (*wantq) {

/*           Accumulate transformation in the matrix Q. */

#line 293 "dlaexc.f"
	    dlarfx_("R", n, &c__3, u, &tau, &q[*j1 * q_dim1 + 1], ldq, &work[
		    1], (ftnlen)1);
#line 294 "dlaexc.f"
	}
#line 295 "dlaexc.f"
	goto L40;

#line 297 "dlaexc.f"
L20:

/*        N1 = 2, N2 = 1: generate elementary reflector H so that: */

/*        H (  -X11 ) = ( * ) */
/*          (  -X21 ) = ( 0 ) */
/*          ( scale ) = ( 0 ) */

#line 305 "dlaexc.f"
	u[0] = -x[0];
#line 306 "dlaexc.f"
	u[1] = -x[1];
#line 307 "dlaexc.f"
	u[2] = scale;
#line 308 "dlaexc.f"
	dlarfg_(&c__3, u, &u[1], &c__1, &tau);
#line 309 "dlaexc.f"
	u[0] = 1.;
#line 310 "dlaexc.f"
	t33 = t[j3 + j3 * t_dim1];

/*        Perform swap provisionally on diagonal block in D. */

#line 314 "dlaexc.f"
	dlarfx_("L", &c__3, &c__3, u, &tau, d__, &c__4, &work[1], (ftnlen)1);
#line 315 "dlaexc.f"
	dlarfx_("R", &c__3, &c__3, u, &tau, d__, &c__4, &work[1], (ftnlen)1);

/*        Test whether to reject swap. */

/* Computing MAX */
#line 319 "dlaexc.f"
	d__2 = abs(d__[1]), d__3 = abs(d__[2]), d__2 = max(d__2,d__3), d__3 = 
		(d__1 = d__[0] - t33, abs(d__1));
#line 319 "dlaexc.f"
	if (max(d__2,d__3) > thresh) {
#line 319 "dlaexc.f"
	    goto L50;
#line 319 "dlaexc.f"
	}

/*        Accept swap: apply transformation to the entire matrix T. */

#line 324 "dlaexc.f"
	dlarfx_("R", &j3, &c__3, u, &tau, &t[*j1 * t_dim1 + 1], ldt, &work[1],
		 (ftnlen)1);
#line 325 "dlaexc.f"
	i__1 = *n - *j1;
#line 325 "dlaexc.f"
	dlarfx_("L", &c__3, &i__1, u, &tau, &t[*j1 + j2 * t_dim1], ldt, &work[
		1], (ftnlen)1);

#line 327 "dlaexc.f"
	t[*j1 + *j1 * t_dim1] = t33;
#line 328 "dlaexc.f"
	t[j2 + *j1 * t_dim1] = 0.;
#line 329 "dlaexc.f"
	t[j3 + *j1 * t_dim1] = 0.;

#line 331 "dlaexc.f"
	if (*wantq) {

/*           Accumulate transformation in the matrix Q. */

#line 335 "dlaexc.f"
	    dlarfx_("R", n, &c__3, u, &tau, &q[*j1 * q_dim1 + 1], ldq, &work[
		    1], (ftnlen)1);
#line 336 "dlaexc.f"
	}
#line 337 "dlaexc.f"
	goto L40;

#line 339 "dlaexc.f"
L30:

/*        N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so */
/*        that: */

/*        H(2) H(1) (  -X11  -X12 ) = (  *  * ) */
/*                  (  -X21  -X22 )   (  0  * ) */
/*                  ( scale    0  )   (  0  0 ) */
/*                  (    0  scale )   (  0  0 ) */

#line 349 "dlaexc.f"
	u1[0] = -x[0];
#line 350 "dlaexc.f"
	u1[1] = -x[1];
#line 351 "dlaexc.f"
	u1[2] = scale;
#line 352 "dlaexc.f"
	dlarfg_(&c__3, u1, &u1[1], &c__1, &tau1);
#line 353 "dlaexc.f"
	u1[0] = 1.;

#line 355 "dlaexc.f"
	temp = -tau1 * (x[2] + u1[1] * x[3]);
#line 356 "dlaexc.f"
	u2[0] = -temp * u1[1] - x[3];
#line 357 "dlaexc.f"
	u2[1] = -temp * u1[2];
#line 358 "dlaexc.f"
	u2[2] = scale;
#line 359 "dlaexc.f"
	dlarfg_(&c__3, u2, &u2[1], &c__1, &tau2);
#line 360 "dlaexc.f"
	u2[0] = 1.;

/*        Perform swap provisionally on diagonal block in D. */

#line 364 "dlaexc.f"
	dlarfx_("L", &c__3, &c__4, u1, &tau1, d__, &c__4, &work[1], (ftnlen)1)
		;
#line 365 "dlaexc.f"
	dlarfx_("R", &c__4, &c__3, u1, &tau1, d__, &c__4, &work[1], (ftnlen)1)
		;
#line 366 "dlaexc.f"
	dlarfx_("L", &c__3, &c__4, u2, &tau2, &d__[1], &c__4, &work[1], (
		ftnlen)1);
#line 367 "dlaexc.f"
	dlarfx_("R", &c__4, &c__3, u2, &tau2, &d__[4], &c__4, &work[1], (
		ftnlen)1);

/*        Test whether to reject swap. */

/* Computing MAX */
#line 371 "dlaexc.f"
	d__1 = abs(d__[2]), d__2 = abs(d__[6]), d__1 = max(d__1,d__2), d__2 = 
		abs(d__[3]), d__1 = max(d__1,d__2), d__2 = abs(d__[7]);
#line 371 "dlaexc.f"
	if (max(d__1,d__2) > thresh) {
#line 371 "dlaexc.f"
	    goto L50;
#line 371 "dlaexc.f"
	}

/*        Accept swap: apply transformation to the entire matrix T. */

#line 376 "dlaexc.f"
	i__1 = *n - *j1 + 1;
#line 376 "dlaexc.f"
	dlarfx_("L", &c__3, &i__1, u1, &tau1, &t[*j1 + *j1 * t_dim1], ldt, &
		work[1], (ftnlen)1);
#line 377 "dlaexc.f"
	dlarfx_("R", &j4, &c__3, u1, &tau1, &t[*j1 * t_dim1 + 1], ldt, &work[
		1], (ftnlen)1);
#line 378 "dlaexc.f"
	i__1 = *n - *j1 + 1;
#line 378 "dlaexc.f"
	dlarfx_("L", &c__3, &i__1, u2, &tau2, &t[j2 + *j1 * t_dim1], ldt, &
		work[1], (ftnlen)1);
#line 379 "dlaexc.f"
	dlarfx_("R", &j4, &c__3, u2, &tau2, &t[j2 * t_dim1 + 1], ldt, &work[1]
		, (ftnlen)1);

#line 381 "dlaexc.f"
	t[j3 + *j1 * t_dim1] = 0.;
#line 382 "dlaexc.f"
	t[j3 + j2 * t_dim1] = 0.;
#line 383 "dlaexc.f"
	t[j4 + *j1 * t_dim1] = 0.;
#line 384 "dlaexc.f"
	t[j4 + j2 * t_dim1] = 0.;

#line 386 "dlaexc.f"
	if (*wantq) {

/*           Accumulate transformation in the matrix Q. */

#line 390 "dlaexc.f"
	    dlarfx_("R", n, &c__3, u1, &tau1, &q[*j1 * q_dim1 + 1], ldq, &
		    work[1], (ftnlen)1);
#line 391 "dlaexc.f"
	    dlarfx_("R", n, &c__3, u2, &tau2, &q[j2 * q_dim1 + 1], ldq, &work[
		    1], (ftnlen)1);
#line 392 "dlaexc.f"
	}

#line 394 "dlaexc.f"
L40:

#line 396 "dlaexc.f"
	if (*n2 == 2) {

/*           Standardize new 2-by-2 block T11 */

#line 400 "dlaexc.f"
	    dlanv2_(&t[*j1 + *j1 * t_dim1], &t[*j1 + j2 * t_dim1], &t[j2 + *
		    j1 * t_dim1], &t[j2 + j2 * t_dim1], &wr1, &wi1, &wr2, &
		    wi2, &cs, &sn);
#line 402 "dlaexc.f"
	    i__1 = *n - *j1 - 1;
#line 402 "dlaexc.f"
	    drot_(&i__1, &t[*j1 + (*j1 + 2) * t_dim1], ldt, &t[j2 + (*j1 + 2) 
		    * t_dim1], ldt, &cs, &sn);
#line 404 "dlaexc.f"
	    i__1 = *j1 - 1;
#line 404 "dlaexc.f"
	    drot_(&i__1, &t[*j1 * t_dim1 + 1], &c__1, &t[j2 * t_dim1 + 1], &
		    c__1, &cs, &sn);
#line 405 "dlaexc.f"
	    if (*wantq) {
#line 405 "dlaexc.f"
		drot_(n, &q[*j1 * q_dim1 + 1], &c__1, &q[j2 * q_dim1 + 1], &
			c__1, &cs, &sn);
#line 405 "dlaexc.f"
	    }
#line 407 "dlaexc.f"
	}

#line 409 "dlaexc.f"
	if (*n1 == 2) {

/*           Standardize new 2-by-2 block T22 */

#line 413 "dlaexc.f"
	    j3 = *j1 + *n2;
#line 414 "dlaexc.f"
	    j4 = j3 + 1;
#line 415 "dlaexc.f"
	    dlanv2_(&t[j3 + j3 * t_dim1], &t[j3 + j4 * t_dim1], &t[j4 + j3 * 
		    t_dim1], &t[j4 + j4 * t_dim1], &wr1, &wi1, &wr2, &wi2, &
		    cs, &sn);
#line 417 "dlaexc.f"
	    if (j3 + 2 <= *n) {
#line 417 "dlaexc.f"
		i__1 = *n - j3 - 1;
#line 417 "dlaexc.f"
		drot_(&i__1, &t[j3 + (j3 + 2) * t_dim1], ldt, &t[j4 + (j3 + 2)
			 * t_dim1], ldt, &cs, &sn);
#line 417 "dlaexc.f"
	    }
#line 420 "dlaexc.f"
	    i__1 = j3 - 1;
#line 420 "dlaexc.f"
	    drot_(&i__1, &t[j3 * t_dim1 + 1], &c__1, &t[j4 * t_dim1 + 1], &
		    c__1, &cs, &sn);
#line 421 "dlaexc.f"
	    if (*wantq) {
#line 421 "dlaexc.f"
		drot_(n, &q[j3 * q_dim1 + 1], &c__1, &q[j4 * q_dim1 + 1], &
			c__1, &cs, &sn);
#line 421 "dlaexc.f"
	    }
#line 423 "dlaexc.f"
	}

#line 425 "dlaexc.f"
    }
#line 426 "dlaexc.f"
    return 0;

/*     Exit with INFO = 1 if swap was rejected. */

#line 430 "dlaexc.f"
L50:
#line 431 "dlaexc.f"
    *info = 1;
#line 432 "dlaexc.f"
    return 0;

/*     End of DLAEXC */

} /* dlaexc_ */

