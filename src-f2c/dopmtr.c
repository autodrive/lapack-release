#line 1 "dopmtr.f"
/* dopmtr.f -- translated by f2c (version 20100827).
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

#line 1 "dopmtr.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DOPMTR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DOPMTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dopmtr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dopmtr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dopmtr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DOPMTR( SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS, UPLO */
/*       INTEGER            INFO, LDC, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AP( * ), C( LDC, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DOPMTR overwrites the general real M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'T':      Q**T * C       C * Q**T */
/* > */
/* > where Q is a real orthogonal matrix of order nq, with nq = m if */
/* > SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of */
/* > nq-1 elementary reflectors, as returned by DSPTRD using packed */
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
/* >          = 'L': apply Q or Q**T from the Left; */
/* >          = 'R': apply Q or Q**T from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U': Upper triangular packed storage used in previous */
/* >                 call to DSPTRD; */
/* >          = 'L': Lower triangular packed storage used in previous */
/* >                 call to DSPTRD. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q; */
/* >          = 'T':  Transpose, apply Q**T. */
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
/* >          AP is DOUBLE PRECISION array, dimension */
/* >                               (M*(M+1)/2) if SIDE = 'L' */
/* >                               (N*(N+1)/2) if SIDE = 'R' */
/* >          The vectors which define the elementary reflectors, as */
/* >          returned by DSPTRD.  AP is modified by the routine but */
/* >          restored on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION array, dimension (M-1) if SIDE = 'L' */
/* >                                     or (N-1) if SIDE = 'R' */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by DSPTRD. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. */
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
/* >          WORK is DOUBLE PRECISION array, dimension */
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

/* > \date November 2011 */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dopmtr_(char *side, char *uplo, char *trans, integer *m, 
	integer *n, doublereal *ap, doublereal *tau, doublereal *c__, integer 
	*ldc, doublereal *work, integer *info, ftnlen side_len, ftnlen 
	uplo_len, ftnlen trans_len)
{
    /* System generated locals */
    integer c_dim1, c_offset, i__1, i__2;

    /* Local variables */
    static integer i__, i1, i2, i3, ic, jc, ii, mi, ni, nq;
    static doublereal aii;
    static logical left;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical notran, forwrd;


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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 191 "dopmtr.f"
    /* Parameter adjustments */
#line 191 "dopmtr.f"
    --ap;
#line 191 "dopmtr.f"
    --tau;
#line 191 "dopmtr.f"
    c_dim1 = *ldc;
#line 191 "dopmtr.f"
    c_offset = 1 + c_dim1;
#line 191 "dopmtr.f"
    c__ -= c_offset;
#line 191 "dopmtr.f"
    --work;
#line 191 "dopmtr.f"

#line 191 "dopmtr.f"
    /* Function Body */
#line 191 "dopmtr.f"
    *info = 0;
#line 192 "dopmtr.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 193 "dopmtr.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 194 "dopmtr.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

/*     NQ is the order of Q */

#line 198 "dopmtr.f"
    if (left) {
#line 199 "dopmtr.f"
	nq = *m;
#line 200 "dopmtr.f"
    } else {
#line 201 "dopmtr.f"
	nq = *n;
#line 202 "dopmtr.f"
    }
#line 203 "dopmtr.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 204 "dopmtr.f"
	*info = -1;
#line 205 "dopmtr.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 206 "dopmtr.f"
	*info = -2;
#line 207 "dopmtr.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 208 "dopmtr.f"
	*info = -3;
#line 209 "dopmtr.f"
    } else if (*m < 0) {
#line 210 "dopmtr.f"
	*info = -4;
#line 211 "dopmtr.f"
    } else if (*n < 0) {
#line 212 "dopmtr.f"
	*info = -5;
#line 213 "dopmtr.f"
    } else if (*ldc < max(1,*m)) {
#line 214 "dopmtr.f"
	*info = -9;
#line 215 "dopmtr.f"
    }
#line 216 "dopmtr.f"
    if (*info != 0) {
#line 217 "dopmtr.f"
	i__1 = -(*info);
#line 217 "dopmtr.f"
	xerbla_("DOPMTR", &i__1, (ftnlen)6);
#line 218 "dopmtr.f"
	return 0;
#line 219 "dopmtr.f"
    }

/*     Quick return if possible */

#line 223 "dopmtr.f"
    if (*m == 0 || *n == 0) {
#line 223 "dopmtr.f"
	return 0;
#line 223 "dopmtr.f"
    }

#line 226 "dopmtr.f"
    if (upper) {

/*        Q was determined by a call to DSPTRD with UPLO = 'U' */

#line 230 "dopmtr.f"
	forwrd = left && notran || ! left && ! notran;

#line 233 "dopmtr.f"
	if (forwrd) {
#line 234 "dopmtr.f"
	    i1 = 1;
#line 235 "dopmtr.f"
	    i2 = nq - 1;
#line 236 "dopmtr.f"
	    i3 = 1;
#line 237 "dopmtr.f"
	    ii = 2;
#line 238 "dopmtr.f"
	} else {
#line 239 "dopmtr.f"
	    i1 = nq - 1;
#line 240 "dopmtr.f"
	    i2 = 1;
#line 241 "dopmtr.f"
	    i3 = -1;
#line 242 "dopmtr.f"
	    ii = nq * (nq + 1) / 2 - 1;
#line 243 "dopmtr.f"
	}

#line 245 "dopmtr.f"
	if (left) {
#line 246 "dopmtr.f"
	    ni = *n;
#line 247 "dopmtr.f"
	} else {
#line 248 "dopmtr.f"
	    mi = *m;
#line 249 "dopmtr.f"
	}

#line 251 "dopmtr.f"
	i__1 = i2;
#line 251 "dopmtr.f"
	i__2 = i3;
#line 251 "dopmtr.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 252 "dopmtr.f"
	    if (left) {

/*              H(i) is applied to C(1:i,1:n) */

#line 256 "dopmtr.f"
		mi = i__;
#line 257 "dopmtr.f"
	    } else {

/*              H(i) is applied to C(1:m,1:i) */

#line 261 "dopmtr.f"
		ni = i__;
#line 262 "dopmtr.f"
	    }

/*           Apply H(i) */

#line 266 "dopmtr.f"
	    aii = ap[ii];
#line 267 "dopmtr.f"
	    ap[ii] = 1.;
#line 268 "dopmtr.f"
	    dlarf_(side, &mi, &ni, &ap[ii - i__ + 1], &c__1, &tau[i__], &c__[
		    c_offset], ldc, &work[1], (ftnlen)1);
#line 270 "dopmtr.f"
	    ap[ii] = aii;

#line 272 "dopmtr.f"
	    if (forwrd) {
#line 273 "dopmtr.f"
		ii = ii + i__ + 2;
#line 274 "dopmtr.f"
	    } else {
#line 275 "dopmtr.f"
		ii = ii - i__ - 1;
#line 276 "dopmtr.f"
	    }
#line 277 "dopmtr.f"
/* L10: */
#line 277 "dopmtr.f"
	}
#line 278 "dopmtr.f"
    } else {

/*        Q was determined by a call to DSPTRD with UPLO = 'L'. */

#line 282 "dopmtr.f"
	forwrd = left && ! notran || ! left && notran;

#line 285 "dopmtr.f"
	if (forwrd) {
#line 286 "dopmtr.f"
	    i1 = 1;
#line 287 "dopmtr.f"
	    i2 = nq - 1;
#line 288 "dopmtr.f"
	    i3 = 1;
#line 289 "dopmtr.f"
	    ii = 2;
#line 290 "dopmtr.f"
	} else {
#line 291 "dopmtr.f"
	    i1 = nq - 1;
#line 292 "dopmtr.f"
	    i2 = 1;
#line 293 "dopmtr.f"
	    i3 = -1;
#line 294 "dopmtr.f"
	    ii = nq * (nq + 1) / 2 - 1;
#line 295 "dopmtr.f"
	}

#line 297 "dopmtr.f"
	if (left) {
#line 298 "dopmtr.f"
	    ni = *n;
#line 299 "dopmtr.f"
	    jc = 1;
#line 300 "dopmtr.f"
	} else {
#line 301 "dopmtr.f"
	    mi = *m;
#line 302 "dopmtr.f"
	    ic = 1;
#line 303 "dopmtr.f"
	}

#line 305 "dopmtr.f"
	i__2 = i2;
#line 305 "dopmtr.f"
	i__1 = i3;
#line 305 "dopmtr.f"
	for (i__ = i1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
#line 306 "dopmtr.f"
	    aii = ap[ii];
#line 307 "dopmtr.f"
	    ap[ii] = 1.;
#line 308 "dopmtr.f"
	    if (left) {

/*              H(i) is applied to C(i+1:m,1:n) */

#line 312 "dopmtr.f"
		mi = *m - i__;
#line 313 "dopmtr.f"
		ic = i__ + 1;
#line 314 "dopmtr.f"
	    } else {

/*              H(i) is applied to C(1:m,i+1:n) */

#line 318 "dopmtr.f"
		ni = *n - i__;
#line 319 "dopmtr.f"
		jc = i__ + 1;
#line 320 "dopmtr.f"
	    }

/*           Apply H(i) */

#line 324 "dopmtr.f"
	    dlarf_(side, &mi, &ni, &ap[ii], &c__1, &tau[i__], &c__[ic + jc * 
		    c_dim1], ldc, &work[1], (ftnlen)1);
#line 326 "dopmtr.f"
	    ap[ii] = aii;

#line 328 "dopmtr.f"
	    if (forwrd) {
#line 329 "dopmtr.f"
		ii = ii + nq - i__ + 1;
#line 330 "dopmtr.f"
	    } else {
#line 331 "dopmtr.f"
		ii = ii - nq + i__ - 2;
#line 332 "dopmtr.f"
	    }
#line 333 "dopmtr.f"
/* L20: */
#line 333 "dopmtr.f"
	}
#line 334 "dopmtr.f"
    }
#line 335 "dopmtr.f"
    return 0;

/*     End of DOPMTR */

} /* dopmtr_ */

