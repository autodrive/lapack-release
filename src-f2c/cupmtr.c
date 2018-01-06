#line 1 "cupmtr.f"
/* cupmtr.f -- translated by f2c (version 20100827).
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

#line 1 "cupmtr.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CUPMTR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUPMTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cupmtr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cupmtr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cupmtr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CUPMTR( SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS, UPLO */
/*       INTEGER            INFO, LDC, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            AP( * ), C( LDC, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CUPMTR overwrites the general complex M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'C':      Q**H * C       C * Q**H */
/* > */
/* > where Q is a complex unitary matrix of order nq, with nq = m if */
/* > SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of */
/* > nq-1 elementary reflectors, as returned by CHPTRD using packed */
/* > storage: */
/* > */
/* > if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1); */
/* > */
/* > if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q**H from the Left; */
/* >          = 'R': apply Q or Q**H from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U': Upper triangular packed storage used in previous */
/* >                 call to CHPTRD; */
/* >          = 'L': Lower triangular packed storage used in previous */
/* >                 call to CHPTRD. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q; */
/* >          = 'C':  Conjugate transpose, apply Q**H. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix C. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix C. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX array, dimension */
/* >                               (M*(M+1)/2) if SIDE = 'L' */
/* >                               (N*(N+1)/2) if SIDE = 'R' */
/* >          The vectors which define the elementary reflectors, as */
/* >          returned by CHPTRD.  AP is modified by the routine but */
/* >          restored on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX array, dimension (M-1) if SIDE = 'L' */
/* >                                     or (N-1) if SIDE = 'R' */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by CHPTRD. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. LDC >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension */
/* >                                   (N) if SIDE = 'L' */
/* >                                   (M) if SIDE = 'R' */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cupmtr_(char *side, char *uplo, char *trans, integer *m, 
	integer *n, doublecomplex *ap, doublecomplex *tau, doublecomplex *c__,
	 integer *ldc, doublecomplex *work, integer *info, ftnlen side_len, 
	ftnlen uplo_len, ftnlen trans_len)
{
    /* System generated locals */
    integer c_dim1, c_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, i1, i2, i3, ic, jc, ii, mi, ni, nq;
    static doublecomplex aii;
    static logical left;
    static doublecomplex taui;
    extern /* Subroutine */ int clarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical notran, forwrd;


/*  -- LAPACK computational routine (version 3.7.0) -- */
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
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 191 "cupmtr.f"
    /* Parameter adjustments */
#line 191 "cupmtr.f"
    --ap;
#line 191 "cupmtr.f"
    --tau;
#line 191 "cupmtr.f"
    c_dim1 = *ldc;
#line 191 "cupmtr.f"
    c_offset = 1 + c_dim1;
#line 191 "cupmtr.f"
    c__ -= c_offset;
#line 191 "cupmtr.f"
    --work;
#line 191 "cupmtr.f"

#line 191 "cupmtr.f"
    /* Function Body */
#line 191 "cupmtr.f"
    *info = 0;
#line 192 "cupmtr.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 193 "cupmtr.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 194 "cupmtr.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

/*     NQ is the order of Q */

#line 198 "cupmtr.f"
    if (left) {
#line 199 "cupmtr.f"
	nq = *m;
#line 200 "cupmtr.f"
    } else {
#line 201 "cupmtr.f"
	nq = *n;
#line 202 "cupmtr.f"
    }
#line 203 "cupmtr.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 204 "cupmtr.f"
	*info = -1;
#line 205 "cupmtr.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 206 "cupmtr.f"
	*info = -2;
#line 207 "cupmtr.f"
    } else if (! notran && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 208 "cupmtr.f"
	*info = -3;
#line 209 "cupmtr.f"
    } else if (*m < 0) {
#line 210 "cupmtr.f"
	*info = -4;
#line 211 "cupmtr.f"
    } else if (*n < 0) {
#line 212 "cupmtr.f"
	*info = -5;
#line 213 "cupmtr.f"
    } else if (*ldc < max(1,*m)) {
#line 214 "cupmtr.f"
	*info = -9;
#line 215 "cupmtr.f"
    }
#line 216 "cupmtr.f"
    if (*info != 0) {
#line 217 "cupmtr.f"
	i__1 = -(*info);
#line 217 "cupmtr.f"
	xerbla_("CUPMTR", &i__1, (ftnlen)6);
#line 218 "cupmtr.f"
	return 0;
#line 219 "cupmtr.f"
    }

/*     Quick return if possible */

#line 223 "cupmtr.f"
    if (*m == 0 || *n == 0) {
#line 223 "cupmtr.f"
	return 0;
#line 223 "cupmtr.f"
    }

#line 226 "cupmtr.f"
    if (upper) {

/*        Q was determined by a call to CHPTRD with UPLO = 'U' */

#line 230 "cupmtr.f"
	forwrd = left && notran || ! left && ! notran;

#line 233 "cupmtr.f"
	if (forwrd) {
#line 234 "cupmtr.f"
	    i1 = 1;
#line 235 "cupmtr.f"
	    i2 = nq - 1;
#line 236 "cupmtr.f"
	    i3 = 1;
#line 237 "cupmtr.f"
	    ii = 2;
#line 238 "cupmtr.f"
	} else {
#line 239 "cupmtr.f"
	    i1 = nq - 1;
#line 240 "cupmtr.f"
	    i2 = 1;
#line 241 "cupmtr.f"
	    i3 = -1;
#line 242 "cupmtr.f"
	    ii = nq * (nq + 1) / 2 - 1;
#line 243 "cupmtr.f"
	}

#line 245 "cupmtr.f"
	if (left) {
#line 246 "cupmtr.f"
	    ni = *n;
#line 247 "cupmtr.f"
	} else {
#line 248 "cupmtr.f"
	    mi = *m;
#line 249 "cupmtr.f"
	}

#line 251 "cupmtr.f"
	i__1 = i2;
#line 251 "cupmtr.f"
	i__2 = i3;
#line 251 "cupmtr.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 252 "cupmtr.f"
	    if (left) {

/*              H(i) or H(i)**H is applied to C(1:i,1:n) */

#line 256 "cupmtr.f"
		mi = i__;
#line 257 "cupmtr.f"
	    } else {

/*              H(i) or H(i)**H is applied to C(1:m,1:i) */

#line 261 "cupmtr.f"
		ni = i__;
#line 262 "cupmtr.f"
	    }

/*           Apply H(i) or H(i)**H */

#line 266 "cupmtr.f"
	    if (notran) {
#line 267 "cupmtr.f"
		i__3 = i__;
#line 267 "cupmtr.f"
		taui.r = tau[i__3].r, taui.i = tau[i__3].i;
#line 268 "cupmtr.f"
	    } else {
#line 269 "cupmtr.f"
		d_cnjg(&z__1, &tau[i__]);
#line 269 "cupmtr.f"
		taui.r = z__1.r, taui.i = z__1.i;
#line 270 "cupmtr.f"
	    }
#line 271 "cupmtr.f"
	    i__3 = ii;
#line 271 "cupmtr.f"
	    aii.r = ap[i__3].r, aii.i = ap[i__3].i;
#line 272 "cupmtr.f"
	    i__3 = ii;
#line 272 "cupmtr.f"
	    ap[i__3].r = 1., ap[i__3].i = 0.;
#line 273 "cupmtr.f"
	    clarf_(side, &mi, &ni, &ap[ii - i__ + 1], &c__1, &taui, &c__[
		    c_offset], ldc, &work[1], (ftnlen)1);
#line 275 "cupmtr.f"
	    i__3 = ii;
#line 275 "cupmtr.f"
	    ap[i__3].r = aii.r, ap[i__3].i = aii.i;

#line 277 "cupmtr.f"
	    if (forwrd) {
#line 278 "cupmtr.f"
		ii = ii + i__ + 2;
#line 279 "cupmtr.f"
	    } else {
#line 280 "cupmtr.f"
		ii = ii - i__ - 1;
#line 281 "cupmtr.f"
	    }
#line 282 "cupmtr.f"
/* L10: */
#line 282 "cupmtr.f"
	}
#line 283 "cupmtr.f"
    } else {

/*        Q was determined by a call to CHPTRD with UPLO = 'L'. */

#line 287 "cupmtr.f"
	forwrd = left && ! notran || ! left && notran;

#line 290 "cupmtr.f"
	if (forwrd) {
#line 291 "cupmtr.f"
	    i1 = 1;
#line 292 "cupmtr.f"
	    i2 = nq - 1;
#line 293 "cupmtr.f"
	    i3 = 1;
#line 294 "cupmtr.f"
	    ii = 2;
#line 295 "cupmtr.f"
	} else {
#line 296 "cupmtr.f"
	    i1 = nq - 1;
#line 297 "cupmtr.f"
	    i2 = 1;
#line 298 "cupmtr.f"
	    i3 = -1;
#line 299 "cupmtr.f"
	    ii = nq * (nq + 1) / 2 - 1;
#line 300 "cupmtr.f"
	}

#line 302 "cupmtr.f"
	if (left) {
#line 303 "cupmtr.f"
	    ni = *n;
#line 304 "cupmtr.f"
	    jc = 1;
#line 305 "cupmtr.f"
	} else {
#line 306 "cupmtr.f"
	    mi = *m;
#line 307 "cupmtr.f"
	    ic = 1;
#line 308 "cupmtr.f"
	}

#line 310 "cupmtr.f"
	i__2 = i2;
#line 310 "cupmtr.f"
	i__1 = i3;
#line 310 "cupmtr.f"
	for (i__ = i1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
#line 311 "cupmtr.f"
	    i__3 = ii;
#line 311 "cupmtr.f"
	    aii.r = ap[i__3].r, aii.i = ap[i__3].i;
#line 312 "cupmtr.f"
	    i__3 = ii;
#line 312 "cupmtr.f"
	    ap[i__3].r = 1., ap[i__3].i = 0.;
#line 313 "cupmtr.f"
	    if (left) {

/*              H(i) or H(i)**H is applied to C(i+1:m,1:n) */

#line 317 "cupmtr.f"
		mi = *m - i__;
#line 318 "cupmtr.f"
		ic = i__ + 1;
#line 319 "cupmtr.f"
	    } else {

/*              H(i) or H(i)**H is applied to C(1:m,i+1:n) */

#line 323 "cupmtr.f"
		ni = *n - i__;
#line 324 "cupmtr.f"
		jc = i__ + 1;
#line 325 "cupmtr.f"
	    }

/*           Apply H(i) or H(i)**H */

#line 329 "cupmtr.f"
	    if (notran) {
#line 330 "cupmtr.f"
		i__3 = i__;
#line 330 "cupmtr.f"
		taui.r = tau[i__3].r, taui.i = tau[i__3].i;
#line 331 "cupmtr.f"
	    } else {
#line 332 "cupmtr.f"
		d_cnjg(&z__1, &tau[i__]);
#line 332 "cupmtr.f"
		taui.r = z__1.r, taui.i = z__1.i;
#line 333 "cupmtr.f"
	    }
#line 334 "cupmtr.f"
	    clarf_(side, &mi, &ni, &ap[ii], &c__1, &taui, &c__[ic + jc * 
		    c_dim1], ldc, &work[1], (ftnlen)1);
#line 336 "cupmtr.f"
	    i__3 = ii;
#line 336 "cupmtr.f"
	    ap[i__3].r = aii.r, ap[i__3].i = aii.i;

#line 338 "cupmtr.f"
	    if (forwrd) {
#line 339 "cupmtr.f"
		ii = ii + nq - i__ + 1;
#line 340 "cupmtr.f"
	    } else {
#line 341 "cupmtr.f"
		ii = ii - nq + i__ - 2;
#line 342 "cupmtr.f"
	    }
#line 343 "cupmtr.f"
/* L20: */
#line 343 "cupmtr.f"
	}
#line 344 "cupmtr.f"
    }
#line 345 "cupmtr.f"
    return 0;

/*     End of CUPMTR */

} /* cupmtr_ */

