#line 1 "sopmtr.f"
/* sopmtr.f -- translated by f2c (version 20100827).
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

#line 1 "sopmtr.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SOPMTR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SOPMTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sopmtr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sopmtr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sopmtr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SOPMTR( SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS, UPLO */
/*       INTEGER            INFO, LDC, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AP( * ), C( LDC, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SOPMTR overwrites the general real M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'T':      Q**T * C       C * Q**T */
/* > */
/* > where Q is a real orthogonal matrix of order nq, with nq = m if */
/* > SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of */
/* > nq-1 elementary reflectors, as returned by SSPTRD using packed */
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
/* >                 call to SSPTRD; */
/* >          = 'L': Lower triangular packed storage used in previous */
/* >                 call to SSPTRD. */
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
/* >          AP is REAL array, dimension */
/* >                               (M*(M+1)/2) if SIDE = 'L' */
/* >                               (N*(N+1)/2) if SIDE = 'R' */
/* >          The vectors which define the elementary reflectors, as */
/* >          returned by SSPTRD.  AP is modified by the routine but */
/* >          restored on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is REAL array, dimension (M-1) if SIDE = 'L' */
/* >                                     or (N-1) if SIDE = 'R' */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by SSPTRD. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array, dimension (LDC,N) */
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
/* >          WORK is REAL array, dimension */
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

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sopmtr_(char *side, char *uplo, char *trans, integer *m, 
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
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int slarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);
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

#line 191 "sopmtr.f"
    /* Parameter adjustments */
#line 191 "sopmtr.f"
    --ap;
#line 191 "sopmtr.f"
    --tau;
#line 191 "sopmtr.f"
    c_dim1 = *ldc;
#line 191 "sopmtr.f"
    c_offset = 1 + c_dim1;
#line 191 "sopmtr.f"
    c__ -= c_offset;
#line 191 "sopmtr.f"
    --work;
#line 191 "sopmtr.f"

#line 191 "sopmtr.f"
    /* Function Body */
#line 191 "sopmtr.f"
    *info = 0;
#line 192 "sopmtr.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 193 "sopmtr.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 194 "sopmtr.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

/*     NQ is the order of Q */

#line 198 "sopmtr.f"
    if (left) {
#line 199 "sopmtr.f"
	nq = *m;
#line 200 "sopmtr.f"
    } else {
#line 201 "sopmtr.f"
	nq = *n;
#line 202 "sopmtr.f"
    }
#line 203 "sopmtr.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 204 "sopmtr.f"
	*info = -1;
#line 205 "sopmtr.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 206 "sopmtr.f"
	*info = -2;
#line 207 "sopmtr.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 208 "sopmtr.f"
	*info = -3;
#line 209 "sopmtr.f"
    } else if (*m < 0) {
#line 210 "sopmtr.f"
	*info = -4;
#line 211 "sopmtr.f"
    } else if (*n < 0) {
#line 212 "sopmtr.f"
	*info = -5;
#line 213 "sopmtr.f"
    } else if (*ldc < max(1,*m)) {
#line 214 "sopmtr.f"
	*info = -9;
#line 215 "sopmtr.f"
    }
#line 216 "sopmtr.f"
    if (*info != 0) {
#line 217 "sopmtr.f"
	i__1 = -(*info);
#line 217 "sopmtr.f"
	xerbla_("SOPMTR", &i__1, (ftnlen)6);
#line 218 "sopmtr.f"
	return 0;
#line 219 "sopmtr.f"
    }

/*     Quick return if possible */

#line 223 "sopmtr.f"
    if (*m == 0 || *n == 0) {
#line 223 "sopmtr.f"
	return 0;
#line 223 "sopmtr.f"
    }

#line 226 "sopmtr.f"
    if (upper) {

/*        Q was determined by a call to SSPTRD with UPLO = 'U' */

#line 230 "sopmtr.f"
	forwrd = left && notran || ! left && ! notran;

#line 233 "sopmtr.f"
	if (forwrd) {
#line 234 "sopmtr.f"
	    i1 = 1;
#line 235 "sopmtr.f"
	    i2 = nq - 1;
#line 236 "sopmtr.f"
	    i3 = 1;
#line 237 "sopmtr.f"
	    ii = 2;
#line 238 "sopmtr.f"
	} else {
#line 239 "sopmtr.f"
	    i1 = nq - 1;
#line 240 "sopmtr.f"
	    i2 = 1;
#line 241 "sopmtr.f"
	    i3 = -1;
#line 242 "sopmtr.f"
	    ii = nq * (nq + 1) / 2 - 1;
#line 243 "sopmtr.f"
	}

#line 245 "sopmtr.f"
	if (left) {
#line 246 "sopmtr.f"
	    ni = *n;
#line 247 "sopmtr.f"
	} else {
#line 248 "sopmtr.f"
	    mi = *m;
#line 249 "sopmtr.f"
	}

#line 251 "sopmtr.f"
	i__1 = i2;
#line 251 "sopmtr.f"
	i__2 = i3;
#line 251 "sopmtr.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 252 "sopmtr.f"
	    if (left) {

/*              H(i) is applied to C(1:i,1:n) */

#line 256 "sopmtr.f"
		mi = i__;
#line 257 "sopmtr.f"
	    } else {

/*              H(i) is applied to C(1:m,1:i) */

#line 261 "sopmtr.f"
		ni = i__;
#line 262 "sopmtr.f"
	    }

/*           Apply H(i) */

#line 266 "sopmtr.f"
	    aii = ap[ii];
#line 267 "sopmtr.f"
	    ap[ii] = 1.;
#line 268 "sopmtr.f"
	    slarf_(side, &mi, &ni, &ap[ii - i__ + 1], &c__1, &tau[i__], &c__[
		    c_offset], ldc, &work[1], (ftnlen)1);
#line 270 "sopmtr.f"
	    ap[ii] = aii;

#line 272 "sopmtr.f"
	    if (forwrd) {
#line 273 "sopmtr.f"
		ii = ii + i__ + 2;
#line 274 "sopmtr.f"
	    } else {
#line 275 "sopmtr.f"
		ii = ii - i__ - 1;
#line 276 "sopmtr.f"
	    }
#line 277 "sopmtr.f"
/* L10: */
#line 277 "sopmtr.f"
	}
#line 278 "sopmtr.f"
    } else {

/*        Q was determined by a call to SSPTRD with UPLO = 'L'. */

#line 282 "sopmtr.f"
	forwrd = left && ! notran || ! left && notran;

#line 285 "sopmtr.f"
	if (forwrd) {
#line 286 "sopmtr.f"
	    i1 = 1;
#line 287 "sopmtr.f"
	    i2 = nq - 1;
#line 288 "sopmtr.f"
	    i3 = 1;
#line 289 "sopmtr.f"
	    ii = 2;
#line 290 "sopmtr.f"
	} else {
#line 291 "sopmtr.f"
	    i1 = nq - 1;
#line 292 "sopmtr.f"
	    i2 = 1;
#line 293 "sopmtr.f"
	    i3 = -1;
#line 294 "sopmtr.f"
	    ii = nq * (nq + 1) / 2 - 1;
#line 295 "sopmtr.f"
	}

#line 297 "sopmtr.f"
	if (left) {
#line 298 "sopmtr.f"
	    ni = *n;
#line 299 "sopmtr.f"
	    jc = 1;
#line 300 "sopmtr.f"
	} else {
#line 301 "sopmtr.f"
	    mi = *m;
#line 302 "sopmtr.f"
	    ic = 1;
#line 303 "sopmtr.f"
	}

#line 305 "sopmtr.f"
	i__2 = i2;
#line 305 "sopmtr.f"
	i__1 = i3;
#line 305 "sopmtr.f"
	for (i__ = i1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
#line 306 "sopmtr.f"
	    aii = ap[ii];
#line 307 "sopmtr.f"
	    ap[ii] = 1.;
#line 308 "sopmtr.f"
	    if (left) {

/*              H(i) is applied to C(i+1:m,1:n) */

#line 312 "sopmtr.f"
		mi = *m - i__;
#line 313 "sopmtr.f"
		ic = i__ + 1;
#line 314 "sopmtr.f"
	    } else {

/*              H(i) is applied to C(1:m,i+1:n) */

#line 318 "sopmtr.f"
		ni = *n - i__;
#line 319 "sopmtr.f"
		jc = i__ + 1;
#line 320 "sopmtr.f"
	    }

/*           Apply H(i) */

#line 324 "sopmtr.f"
	    slarf_(side, &mi, &ni, &ap[ii], &c__1, &tau[i__], &c__[ic + jc * 
		    c_dim1], ldc, &work[1], (ftnlen)1);
#line 326 "sopmtr.f"
	    ap[ii] = aii;

#line 328 "sopmtr.f"
	    if (forwrd) {
#line 329 "sopmtr.f"
		ii = ii + nq - i__ + 1;
#line 330 "sopmtr.f"
	    } else {
#line 331 "sopmtr.f"
		ii = ii - nq + i__ - 2;
#line 332 "sopmtr.f"
	    }
#line 333 "sopmtr.f"
/* L20: */
#line 333 "sopmtr.f"
	}
#line 334 "sopmtr.f"
    }
#line 335 "sopmtr.f"
    return 0;

/*     End of SOPMTR */

} /* sopmtr_ */

