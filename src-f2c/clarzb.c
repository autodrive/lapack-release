#line 1 "clarzb.f"
/* clarzb.f -- translated by f2c (version 20100827).
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

#line 1 "clarzb.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CLARZB applies a block reflector or its conjugate-transpose to a general matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLARZB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarzb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarzb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarzb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLARZB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, V, */
/*                          LDV, T, LDT, C, LDC, WORK, LDWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIRECT, SIDE, STOREV, TRANS */
/*       INTEGER            K, L, LDC, LDT, LDV, LDWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            C( LDC, * ), T( LDT, * ), V( LDV, * ), */
/*      $                   WORK( LDWORK, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARZB applies a complex block reflector H or its transpose H**H */
/* > to a complex distributed M-by-N  C from the left or the right. */
/* > */
/* > Currently, only STOREV = 'R' and DIRECT = 'B' are supported. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply H or H**H from the Left */
/* >          = 'R': apply H or H**H from the Right */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N': apply H (No transpose) */
/* >          = 'C': apply H**H (Conjugate transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] DIRECT */
/* > \verbatim */
/* >          DIRECT is CHARACTER*1 */
/* >          Indicates how H is formed from a product of elementary */
/* >          reflectors */
/* >          = 'F': H = H(1) H(2) . . . H(k) (Forward, not supported yet) */
/* >          = 'B': H = H(k) . . . H(2) H(1) (Backward) */
/* > \endverbatim */
/* > */
/* > \param[in] STOREV */
/* > \verbatim */
/* >          STOREV is CHARACTER*1 */
/* >          Indicates how the vectors which define the elementary */
/* >          reflectors are stored: */
/* >          = 'C': Columnwise                        (not supported yet) */
/* >          = 'R': Rowwise */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The order of the matrix T (= the number of elementary */
/* >          reflectors whose product defines the block reflector). */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* >          L is INTEGER */
/* >          The number of columns of the matrix V containing the */
/* >          meaningful part of the Householder reflectors. */
/* >          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is COMPLEX array, dimension (LDV,NV). */
/* >          If STOREV = 'C', NV = K; if STOREV = 'R', NV = L. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V. */
/* >          If STOREV = 'C', LDV >= L; if STOREV = 'R', LDV >= K. */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is COMPLEX array, dimension (LDT,K) */
/* >          The triangular K-by-K matrix T in the representation of the */
/* >          block reflector. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T. LDT >= K. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by H*C or H**H*C or C*H or C*H**H. */
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
/* >          WORK is COMPLEX array, dimension (LDWORK,K) */
/* > \endverbatim */
/* > */
/* > \param[in] LDWORK */
/* > \verbatim */
/* >          LDWORK is INTEGER */
/* >          The leading dimension of the array WORK. */
/* >          If SIDE = 'L', LDWORK >= max(1,N); */
/* >          if SIDE = 'R', LDWORK >= max(1,M). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complexOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int clarzb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, integer *l, doublecomplex 
	*v, integer *ldv, doublecomplex *t, integer *ldt, doublecomplex *c__, 
	integer *ldc, doublecomplex *work, integer *ldwork, ftnlen side_len, 
	ftnlen trans_len, ftnlen direct_len, ftnlen storev_len)
{
    /* System generated locals */
    integer c_dim1, c_offset, t_dim1, t_offset, v_dim1, v_offset, work_dim1, 
	    work_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j, info;
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), ctrmm_(char *, char *, char *, char *
	    , integer *, integer *, doublecomplex *, doublecomplex *, integer 
	    *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    clacgv_(integer *, doublecomplex *, integer *), xerbla_(char *, 
	    integer *, ftnlen);
    static char transt[1];


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 221 "clarzb.f"
    /* Parameter adjustments */
#line 221 "clarzb.f"
    v_dim1 = *ldv;
#line 221 "clarzb.f"
    v_offset = 1 + v_dim1;
#line 221 "clarzb.f"
    v -= v_offset;
#line 221 "clarzb.f"
    t_dim1 = *ldt;
#line 221 "clarzb.f"
    t_offset = 1 + t_dim1;
#line 221 "clarzb.f"
    t -= t_offset;
#line 221 "clarzb.f"
    c_dim1 = *ldc;
#line 221 "clarzb.f"
    c_offset = 1 + c_dim1;
#line 221 "clarzb.f"
    c__ -= c_offset;
#line 221 "clarzb.f"
    work_dim1 = *ldwork;
#line 221 "clarzb.f"
    work_offset = 1 + work_dim1;
#line 221 "clarzb.f"
    work -= work_offset;
#line 221 "clarzb.f"

#line 221 "clarzb.f"
    /* Function Body */
#line 221 "clarzb.f"
    if (*m <= 0 || *n <= 0) {
#line 221 "clarzb.f"
	return 0;
#line 221 "clarzb.f"
    }

/*     Check for currently supported options */

#line 226 "clarzb.f"
    info = 0;
#line 227 "clarzb.f"
    if (! lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 228 "clarzb.f"
	info = -3;
#line 229 "clarzb.f"
    } else if (! lsame_(storev, "R", (ftnlen)1, (ftnlen)1)) {
#line 230 "clarzb.f"
	info = -4;
#line 231 "clarzb.f"
    }
#line 232 "clarzb.f"
    if (info != 0) {
#line 233 "clarzb.f"
	i__1 = -info;
#line 233 "clarzb.f"
	xerbla_("CLARZB", &i__1, (ftnlen)6);
#line 234 "clarzb.f"
	return 0;
#line 235 "clarzb.f"
    }

#line 237 "clarzb.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 238 "clarzb.f"
	*(unsigned char *)transt = 'C';
#line 239 "clarzb.f"
    } else {
#line 240 "clarzb.f"
	*(unsigned char *)transt = 'N';
#line 241 "clarzb.f"
    }

#line 243 "clarzb.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*        Form  H * C  or  H**H * C */

/*        W( 1:n, 1:k ) = C( 1:k, 1:n )**H */

#line 249 "clarzb.f"
	i__1 = *k;
#line 249 "clarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 250 "clarzb.f"
	    ccopy_(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1], &c__1);
#line 251 "clarzb.f"
/* L10: */
#line 251 "clarzb.f"
	}

/*        W( 1:n, 1:k ) = W( 1:n, 1:k ) + ... */
/*                        C( m-l+1:m, 1:n )**H * V( 1:k, 1:l )**T */

#line 256 "clarzb.f"
	if (*l > 0) {
#line 256 "clarzb.f"
	    cgemm_("Transpose", "Conjugate transpose", n, k, l, &c_b1, &c__[*
		    m - *l + 1 + c_dim1], ldc, &v[v_offset], ldv, &c_b1, &
		    work[work_offset], ldwork, (ftnlen)9, (ftnlen)19);
#line 256 "clarzb.f"
	}

/*        W( 1:n, 1:k ) = W( 1:n, 1:k ) * T**T  or  W( 1:m, 1:k ) * T */

#line 263 "clarzb.f"
	ctrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b1, &t[t_offset]
		, ldt, &work[work_offset], ldwork, (ftnlen)5, (ftnlen)5, (
		ftnlen)1, (ftnlen)8);

/*        C( 1:k, 1:n ) = C( 1:k, 1:n ) - W( 1:n, 1:k )**H */

#line 268 "clarzb.f"
	i__1 = *n;
#line 268 "clarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 269 "clarzb.f"
	    i__2 = *k;
#line 269 "clarzb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 270 "clarzb.f"
		i__3 = i__ + j * c_dim1;
#line 270 "clarzb.f"
		i__4 = i__ + j * c_dim1;
#line 270 "clarzb.f"
		i__5 = j + i__ * work_dim1;
#line 270 "clarzb.f"
		z__1.r = c__[i__4].r - work[i__5].r, z__1.i = c__[i__4].i - 
			work[i__5].i;
#line 270 "clarzb.f"
		c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 271 "clarzb.f"
/* L20: */
#line 271 "clarzb.f"
	    }
#line 272 "clarzb.f"
/* L30: */
#line 272 "clarzb.f"
	}

/*        C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ... */
/*                            V( 1:k, 1:l )**H * W( 1:n, 1:k )**H */

#line 277 "clarzb.f"
	if (*l > 0) {
#line 277 "clarzb.f"
	    z__1.r = -1., z__1.i = -0.;
#line 277 "clarzb.f"
	    cgemm_("Transpose", "Transpose", l, n, k, &z__1, &v[v_offset], 
		    ldv, &work[work_offset], ldwork, &c_b1, &c__[*m - *l + 1 
		    + c_dim1], ldc, (ftnlen)9, (ftnlen)9);
#line 277 "clarzb.f"
	}

#line 281 "clarzb.f"
    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*        Form  C * H  or  C * H**H */

/*        W( 1:m, 1:k ) = C( 1:m, 1:k ) */

#line 287 "clarzb.f"
	i__1 = *k;
#line 287 "clarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 288 "clarzb.f"
	    ccopy_(m, &c__[j * c_dim1 + 1], &c__1, &work[j * work_dim1 + 1], &
		    c__1);
#line 289 "clarzb.f"
/* L40: */
#line 289 "clarzb.f"
	}

/*        W( 1:m, 1:k ) = W( 1:m, 1:k ) + ... */
/*                        C( 1:m, n-l+1:n ) * V( 1:k, 1:l )**H */

#line 294 "clarzb.f"
	if (*l > 0) {
#line 294 "clarzb.f"
	    cgemm_("No transpose", "Transpose", m, k, l, &c_b1, &c__[(*n - *l 
		    + 1) * c_dim1 + 1], ldc, &v[v_offset], ldv, &c_b1, &work[
		    work_offset], ldwork, (ftnlen)12, (ftnlen)9);
#line 294 "clarzb.f"
	}

/*        W( 1:m, 1:k ) = W( 1:m, 1:k ) * conjg( T )  or */
/*                        W( 1:m, 1:k ) * T**H */

#line 301 "clarzb.f"
	i__1 = *k;
#line 301 "clarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 302 "clarzb.f"
	    i__2 = *k - j + 1;
#line 302 "clarzb.f"
	    clacgv_(&i__2, &t[j + j * t_dim1], &c__1);
#line 303 "clarzb.f"
/* L50: */
#line 303 "clarzb.f"
	}
#line 304 "clarzb.f"
	ctrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b1, &t[t_offset],
		 ldt, &work[work_offset], ldwork, (ftnlen)5, (ftnlen)5, (
		ftnlen)1, (ftnlen)8);
#line 306 "clarzb.f"
	i__1 = *k;
#line 306 "clarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 307 "clarzb.f"
	    i__2 = *k - j + 1;
#line 307 "clarzb.f"
	    clacgv_(&i__2, &t[j + j * t_dim1], &c__1);
#line 308 "clarzb.f"
/* L60: */
#line 308 "clarzb.f"
	}

/*        C( 1:m, 1:k ) = C( 1:m, 1:k ) - W( 1:m, 1:k ) */

#line 312 "clarzb.f"
	i__1 = *k;
#line 312 "clarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 313 "clarzb.f"
	    i__2 = *m;
#line 313 "clarzb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 314 "clarzb.f"
		i__3 = i__ + j * c_dim1;
#line 314 "clarzb.f"
		i__4 = i__ + j * c_dim1;
#line 314 "clarzb.f"
		i__5 = i__ + j * work_dim1;
#line 314 "clarzb.f"
		z__1.r = c__[i__4].r - work[i__5].r, z__1.i = c__[i__4].i - 
			work[i__5].i;
#line 314 "clarzb.f"
		c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 315 "clarzb.f"
/* L70: */
#line 315 "clarzb.f"
	    }
#line 316 "clarzb.f"
/* L80: */
#line 316 "clarzb.f"
	}

/*        C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ... */
/*                            W( 1:m, 1:k ) * conjg( V( 1:k, 1:l ) ) */

#line 321 "clarzb.f"
	i__1 = *l;
#line 321 "clarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 322 "clarzb.f"
	    clacgv_(k, &v[j * v_dim1 + 1], &c__1);
#line 323 "clarzb.f"
/* L90: */
#line 323 "clarzb.f"
	}
#line 324 "clarzb.f"
	if (*l > 0) {
#line 324 "clarzb.f"
	    z__1.r = -1., z__1.i = -0.;
#line 324 "clarzb.f"
	    cgemm_("No transpose", "No transpose", m, l, k, &z__1, &work[
		    work_offset], ldwork, &v[v_offset], ldv, &c_b1, &c__[(*n 
		    - *l + 1) * c_dim1 + 1], ldc, (ftnlen)12, (ftnlen)12);
#line 324 "clarzb.f"
	}
#line 327 "clarzb.f"
	i__1 = *l;
#line 327 "clarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 328 "clarzb.f"
	    clacgv_(k, &v[j * v_dim1 + 1], &c__1);
#line 329 "clarzb.f"
/* L100: */
#line 329 "clarzb.f"
	}

#line 331 "clarzb.f"
    }

#line 333 "clarzb.f"
    return 0;

/*     End of CLARZB */

} /* clarzb_ */

