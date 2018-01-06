#line 1 "zlarzb.f"
/* zlarzb.f -- translated by f2c (version 20100827).
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

#line 1 "zlarzb.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZLARZB applies a block reflector or its conjugate-transpose to a general matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLARZB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarzb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarzb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarzb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLARZB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, V, */
/*                          LDV, T, LDT, C, LDC, WORK, LDWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIRECT, SIDE, STOREV, TRANS */
/*       INTEGER            K, L, LDC, LDT, LDV, LDWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         C( LDC, * ), T( LDT, * ), V( LDV, * ), */
/*      $                   WORK( LDWORK, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLARZB applies a complex block reflector H or its transpose H**H */
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
/* >          V is COMPLEX*16 array, dimension (LDV,NV). */
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
/* >          T is COMPLEX*16 array, dimension (LDT,K) */
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
/* >          C is COMPLEX*16 array, dimension (LDC,N) */
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
/* >          WORK is COMPLEX*16 array, dimension (LDWORK,K) */
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

/* > \date December 2016 */

/* > \ingroup complex16OTHERcomputational */

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
/* Subroutine */ int zlarzb_(char *side, char *trans, char *direct, char *
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
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), zcopy_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), ztrmm_(char *, char *, 
	    char *, char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen), 
	    zlacgv_(integer *, doublecomplex *, integer *);
    static char transt[1];


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
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 221 "zlarzb.f"
    /* Parameter adjustments */
#line 221 "zlarzb.f"
    v_dim1 = *ldv;
#line 221 "zlarzb.f"
    v_offset = 1 + v_dim1;
#line 221 "zlarzb.f"
    v -= v_offset;
#line 221 "zlarzb.f"
    t_dim1 = *ldt;
#line 221 "zlarzb.f"
    t_offset = 1 + t_dim1;
#line 221 "zlarzb.f"
    t -= t_offset;
#line 221 "zlarzb.f"
    c_dim1 = *ldc;
#line 221 "zlarzb.f"
    c_offset = 1 + c_dim1;
#line 221 "zlarzb.f"
    c__ -= c_offset;
#line 221 "zlarzb.f"
    work_dim1 = *ldwork;
#line 221 "zlarzb.f"
    work_offset = 1 + work_dim1;
#line 221 "zlarzb.f"
    work -= work_offset;
#line 221 "zlarzb.f"

#line 221 "zlarzb.f"
    /* Function Body */
#line 221 "zlarzb.f"
    if (*m <= 0 || *n <= 0) {
#line 221 "zlarzb.f"
	return 0;
#line 221 "zlarzb.f"
    }

/*     Check for currently supported options */

#line 226 "zlarzb.f"
    info = 0;
#line 227 "zlarzb.f"
    if (! lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 228 "zlarzb.f"
	info = -3;
#line 229 "zlarzb.f"
    } else if (! lsame_(storev, "R", (ftnlen)1, (ftnlen)1)) {
#line 230 "zlarzb.f"
	info = -4;
#line 231 "zlarzb.f"
    }
#line 232 "zlarzb.f"
    if (info != 0) {
#line 233 "zlarzb.f"
	i__1 = -info;
#line 233 "zlarzb.f"
	xerbla_("ZLARZB", &i__1, (ftnlen)6);
#line 234 "zlarzb.f"
	return 0;
#line 235 "zlarzb.f"
    }

#line 237 "zlarzb.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 238 "zlarzb.f"
	*(unsigned char *)transt = 'C';
#line 239 "zlarzb.f"
    } else {
#line 240 "zlarzb.f"
	*(unsigned char *)transt = 'N';
#line 241 "zlarzb.f"
    }

#line 243 "zlarzb.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*        Form  H * C  or  H**H * C */

/*        W( 1:n, 1:k ) = C( 1:k, 1:n )**H */

#line 249 "zlarzb.f"
	i__1 = *k;
#line 249 "zlarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 250 "zlarzb.f"
	    zcopy_(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1], &c__1);
#line 251 "zlarzb.f"
/* L10: */
#line 251 "zlarzb.f"
	}

/*        W( 1:n, 1:k ) = W( 1:n, 1:k ) + ... */
/*                        C( m-l+1:m, 1:n )**H * V( 1:k, 1:l )**T */

#line 256 "zlarzb.f"
	if (*l > 0) {
#line 256 "zlarzb.f"
	    zgemm_("Transpose", "Conjugate transpose", n, k, l, &c_b1, &c__[*
		    m - *l + 1 + c_dim1], ldc, &v[v_offset], ldv, &c_b1, &
		    work[work_offset], ldwork, (ftnlen)9, (ftnlen)19);
#line 256 "zlarzb.f"
	}

/*        W( 1:n, 1:k ) = W( 1:n, 1:k ) * T**T  or  W( 1:m, 1:k ) * T */

#line 263 "zlarzb.f"
	ztrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b1, &t[t_offset]
		, ldt, &work[work_offset], ldwork, (ftnlen)5, (ftnlen)5, (
		ftnlen)1, (ftnlen)8);

/*        C( 1:k, 1:n ) = C( 1:k, 1:n ) - W( 1:n, 1:k )**H */

#line 268 "zlarzb.f"
	i__1 = *n;
#line 268 "zlarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 269 "zlarzb.f"
	    i__2 = *k;
#line 269 "zlarzb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 270 "zlarzb.f"
		i__3 = i__ + j * c_dim1;
#line 270 "zlarzb.f"
		i__4 = i__ + j * c_dim1;
#line 270 "zlarzb.f"
		i__5 = j + i__ * work_dim1;
#line 270 "zlarzb.f"
		z__1.r = c__[i__4].r - work[i__5].r, z__1.i = c__[i__4].i - 
			work[i__5].i;
#line 270 "zlarzb.f"
		c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 271 "zlarzb.f"
/* L20: */
#line 271 "zlarzb.f"
	    }
#line 272 "zlarzb.f"
/* L30: */
#line 272 "zlarzb.f"
	}

/*        C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ... */
/*                            V( 1:k, 1:l )**H * W( 1:n, 1:k )**H */

#line 277 "zlarzb.f"
	if (*l > 0) {
#line 277 "zlarzb.f"
	    z__1.r = -1., z__1.i = -0.;
#line 277 "zlarzb.f"
	    zgemm_("Transpose", "Transpose", l, n, k, &z__1, &v[v_offset], 
		    ldv, &work[work_offset], ldwork, &c_b1, &c__[*m - *l + 1 
		    + c_dim1], ldc, (ftnlen)9, (ftnlen)9);
#line 277 "zlarzb.f"
	}

#line 281 "zlarzb.f"
    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*        Form  C * H  or  C * H**H */

/*        W( 1:m, 1:k ) = C( 1:m, 1:k ) */

#line 287 "zlarzb.f"
	i__1 = *k;
#line 287 "zlarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 288 "zlarzb.f"
	    zcopy_(m, &c__[j * c_dim1 + 1], &c__1, &work[j * work_dim1 + 1], &
		    c__1);
#line 289 "zlarzb.f"
/* L40: */
#line 289 "zlarzb.f"
	}

/*        W( 1:m, 1:k ) = W( 1:m, 1:k ) + ... */
/*                        C( 1:m, n-l+1:n ) * V( 1:k, 1:l )**H */

#line 294 "zlarzb.f"
	if (*l > 0) {
#line 294 "zlarzb.f"
	    zgemm_("No transpose", "Transpose", m, k, l, &c_b1, &c__[(*n - *l 
		    + 1) * c_dim1 + 1], ldc, &v[v_offset], ldv, &c_b1, &work[
		    work_offset], ldwork, (ftnlen)12, (ftnlen)9);
#line 294 "zlarzb.f"
	}

/*        W( 1:m, 1:k ) = W( 1:m, 1:k ) * conjg( T )  or */
/*                        W( 1:m, 1:k ) * T**H */

#line 301 "zlarzb.f"
	i__1 = *k;
#line 301 "zlarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 302 "zlarzb.f"
	    i__2 = *k - j + 1;
#line 302 "zlarzb.f"
	    zlacgv_(&i__2, &t[j + j * t_dim1], &c__1);
#line 303 "zlarzb.f"
/* L50: */
#line 303 "zlarzb.f"
	}
#line 304 "zlarzb.f"
	ztrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b1, &t[t_offset],
		 ldt, &work[work_offset], ldwork, (ftnlen)5, (ftnlen)5, (
		ftnlen)1, (ftnlen)8);
#line 306 "zlarzb.f"
	i__1 = *k;
#line 306 "zlarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 307 "zlarzb.f"
	    i__2 = *k - j + 1;
#line 307 "zlarzb.f"
	    zlacgv_(&i__2, &t[j + j * t_dim1], &c__1);
#line 308 "zlarzb.f"
/* L60: */
#line 308 "zlarzb.f"
	}

/*        C( 1:m, 1:k ) = C( 1:m, 1:k ) - W( 1:m, 1:k ) */

#line 312 "zlarzb.f"
	i__1 = *k;
#line 312 "zlarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 313 "zlarzb.f"
	    i__2 = *m;
#line 313 "zlarzb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 314 "zlarzb.f"
		i__3 = i__ + j * c_dim1;
#line 314 "zlarzb.f"
		i__4 = i__ + j * c_dim1;
#line 314 "zlarzb.f"
		i__5 = i__ + j * work_dim1;
#line 314 "zlarzb.f"
		z__1.r = c__[i__4].r - work[i__5].r, z__1.i = c__[i__4].i - 
			work[i__5].i;
#line 314 "zlarzb.f"
		c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 315 "zlarzb.f"
/* L70: */
#line 315 "zlarzb.f"
	    }
#line 316 "zlarzb.f"
/* L80: */
#line 316 "zlarzb.f"
	}

/*        C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ... */
/*                            W( 1:m, 1:k ) * conjg( V( 1:k, 1:l ) ) */

#line 321 "zlarzb.f"
	i__1 = *l;
#line 321 "zlarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 322 "zlarzb.f"
	    zlacgv_(k, &v[j * v_dim1 + 1], &c__1);
#line 323 "zlarzb.f"
/* L90: */
#line 323 "zlarzb.f"
	}
#line 324 "zlarzb.f"
	if (*l > 0) {
#line 324 "zlarzb.f"
	    z__1.r = -1., z__1.i = -0.;
#line 324 "zlarzb.f"
	    zgemm_("No transpose", "No transpose", m, l, k, &z__1, &work[
		    work_offset], ldwork, &v[v_offset], ldv, &c_b1, &c__[(*n 
		    - *l + 1) * c_dim1 + 1], ldc, (ftnlen)12, (ftnlen)12);
#line 324 "zlarzb.f"
	}
#line 327 "zlarzb.f"
	i__1 = *l;
#line 327 "zlarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 328 "zlarzb.f"
	    zlacgv_(k, &v[j * v_dim1 + 1], &c__1);
#line 329 "zlarzb.f"
/* L100: */
#line 329 "zlarzb.f"
	}

#line 331 "zlarzb.f"
    }

#line 333 "zlarzb.f"
    return 0;

/*     End of ZLARZB */

} /* zlarzb_ */

