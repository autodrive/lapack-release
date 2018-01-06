#line 1 "dlarzb.f"
/* dlarzb.f -- translated by f2c (version 20100827).
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

#line 1 "dlarzb.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b13 = 1.;
static doublereal c_b23 = -1.;

/* > \brief \b DLARZB applies a block reflector or its transpose to a general matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARZB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarzb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarzb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarzb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARZB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, V, */
/*                          LDV, T, LDT, C, LDC, WORK, LDWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIRECT, SIDE, STOREV, TRANS */
/*       INTEGER            K, L, LDC, LDT, LDV, LDWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ), */
/*      $                   WORK( LDWORK, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARZB applies a real block reflector H or its transpose H**T to */
/* > a real distributed M-by-N  C from the left or the right. */
/* > */
/* > Currently, only STOREV = 'R' and DIRECT = 'B' are supported. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply H or H**T from the Left */
/* >          = 'R': apply H or H**T from the Right */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N': apply H (No transpose) */
/* >          = 'C': apply H**T (Transpose) */
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
/* >          V is DOUBLE PRECISION array, dimension (LDV,NV). */
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
/* >          T is DOUBLE PRECISION array, dimension (LDT,K) */
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
/* >          C is DOUBLE PRECISION array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by H*C or H**T*C or C*H or C*H**T. */
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
/* >          WORK is DOUBLE PRECISION array, dimension (LDWORK,K) */
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

/* > \ingroup doubleOTHERcomputational */

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
/* Subroutine */ int dlarzb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, integer *l, doublereal *v,
	 integer *ldv, doublereal *t, integer *ldt, doublereal *c__, integer *
	ldc, doublereal *work, integer *ldwork, ftnlen side_len, ftnlen 
	trans_len, ftnlen direct_len, ftnlen storev_len)
{
    /* System generated locals */
    integer c_dim1, c_offset, t_dim1, t_offset, v_dim1, v_offset, work_dim1, 
	    work_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, info;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), xerbla_(
	    char *, integer *, ftnlen);
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

#line 221 "dlarzb.f"
    /* Parameter adjustments */
#line 221 "dlarzb.f"
    v_dim1 = *ldv;
#line 221 "dlarzb.f"
    v_offset = 1 + v_dim1;
#line 221 "dlarzb.f"
    v -= v_offset;
#line 221 "dlarzb.f"
    t_dim1 = *ldt;
#line 221 "dlarzb.f"
    t_offset = 1 + t_dim1;
#line 221 "dlarzb.f"
    t -= t_offset;
#line 221 "dlarzb.f"
    c_dim1 = *ldc;
#line 221 "dlarzb.f"
    c_offset = 1 + c_dim1;
#line 221 "dlarzb.f"
    c__ -= c_offset;
#line 221 "dlarzb.f"
    work_dim1 = *ldwork;
#line 221 "dlarzb.f"
    work_offset = 1 + work_dim1;
#line 221 "dlarzb.f"
    work -= work_offset;
#line 221 "dlarzb.f"

#line 221 "dlarzb.f"
    /* Function Body */
#line 221 "dlarzb.f"
    if (*m <= 0 || *n <= 0) {
#line 221 "dlarzb.f"
	return 0;
#line 221 "dlarzb.f"
    }

/*     Check for currently supported options */

#line 226 "dlarzb.f"
    info = 0;
#line 227 "dlarzb.f"
    if (! lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 228 "dlarzb.f"
	info = -3;
#line 229 "dlarzb.f"
    } else if (! lsame_(storev, "R", (ftnlen)1, (ftnlen)1)) {
#line 230 "dlarzb.f"
	info = -4;
#line 231 "dlarzb.f"
    }
#line 232 "dlarzb.f"
    if (info != 0) {
#line 233 "dlarzb.f"
	i__1 = -info;
#line 233 "dlarzb.f"
	xerbla_("DLARZB", &i__1, (ftnlen)6);
#line 234 "dlarzb.f"
	return 0;
#line 235 "dlarzb.f"
    }

#line 237 "dlarzb.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 238 "dlarzb.f"
	*(unsigned char *)transt = 'T';
#line 239 "dlarzb.f"
    } else {
#line 240 "dlarzb.f"
	*(unsigned char *)transt = 'N';
#line 241 "dlarzb.f"
    }

#line 243 "dlarzb.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*        Form  H * C  or  H**T * C */

/*        W( 1:n, 1:k ) = C( 1:k, 1:n )**T */

#line 249 "dlarzb.f"
	i__1 = *k;
#line 249 "dlarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 250 "dlarzb.f"
	    dcopy_(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1], &c__1);
#line 251 "dlarzb.f"
/* L10: */
#line 251 "dlarzb.f"
	}

/*        W( 1:n, 1:k ) = W( 1:n, 1:k ) + ... */
/*                        C( m-l+1:m, 1:n )**T * V( 1:k, 1:l )**T */

#line 256 "dlarzb.f"
	if (*l > 0) {
#line 256 "dlarzb.f"
	    dgemm_("Transpose", "Transpose", n, k, l, &c_b13, &c__[*m - *l + 
		    1 + c_dim1], ldc, &v[v_offset], ldv, &c_b13, &work[
		    work_offset], ldwork, (ftnlen)9, (ftnlen)9);
#line 256 "dlarzb.f"
	}

/*        W( 1:n, 1:k ) = W( 1:n, 1:k ) * T**T  or  W( 1:m, 1:k ) * T */

#line 262 "dlarzb.f"
	dtrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b13, &t[
		t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5, (
		ftnlen)5, (ftnlen)1, (ftnlen)8);

/*        C( 1:k, 1:n ) = C( 1:k, 1:n ) - W( 1:n, 1:k )**T */

#line 267 "dlarzb.f"
	i__1 = *n;
#line 267 "dlarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 268 "dlarzb.f"
	    i__2 = *k;
#line 268 "dlarzb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 269 "dlarzb.f"
		c__[i__ + j * c_dim1] -= work[j + i__ * work_dim1];
#line 270 "dlarzb.f"
/* L20: */
#line 270 "dlarzb.f"
	    }
#line 271 "dlarzb.f"
/* L30: */
#line 271 "dlarzb.f"
	}

/*        C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ... */
/*                            V( 1:k, 1:l )**T * W( 1:n, 1:k )**T */

#line 276 "dlarzb.f"
	if (*l > 0) {
#line 276 "dlarzb.f"
	    dgemm_("Transpose", "Transpose", l, n, k, &c_b23, &v[v_offset], 
		    ldv, &work[work_offset], ldwork, &c_b13, &c__[*m - *l + 1 
		    + c_dim1], ldc, (ftnlen)9, (ftnlen)9);
#line 276 "dlarzb.f"
	}

#line 280 "dlarzb.f"
    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*        Form  C * H  or  C * H**T */

/*        W( 1:m, 1:k ) = C( 1:m, 1:k ) */

#line 286 "dlarzb.f"
	i__1 = *k;
#line 286 "dlarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 287 "dlarzb.f"
	    dcopy_(m, &c__[j * c_dim1 + 1], &c__1, &work[j * work_dim1 + 1], &
		    c__1);
#line 288 "dlarzb.f"
/* L40: */
#line 288 "dlarzb.f"
	}

/*        W( 1:m, 1:k ) = W( 1:m, 1:k ) + ... */
/*                        C( 1:m, n-l+1:n ) * V( 1:k, 1:l )**T */

#line 293 "dlarzb.f"
	if (*l > 0) {
#line 293 "dlarzb.f"
	    dgemm_("No transpose", "Transpose", m, k, l, &c_b13, &c__[(*n - *
		    l + 1) * c_dim1 + 1], ldc, &v[v_offset], ldv, &c_b13, &
		    work[work_offset], ldwork, (ftnlen)12, (ftnlen)9);
#line 293 "dlarzb.f"
	}

/*        W( 1:m, 1:k ) = W( 1:m, 1:k ) * T  or  W( 1:m, 1:k ) * T**T */

#line 299 "dlarzb.f"
	dtrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b13, &t[t_offset]
		, ldt, &work[work_offset], ldwork, (ftnlen)5, (ftnlen)5, (
		ftnlen)1, (ftnlen)8);

/*        C( 1:m, 1:k ) = C( 1:m, 1:k ) - W( 1:m, 1:k ) */

#line 304 "dlarzb.f"
	i__1 = *k;
#line 304 "dlarzb.f"
	for (j = 1; j <= i__1; ++j) {
#line 305 "dlarzb.f"
	    i__2 = *m;
#line 305 "dlarzb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 306 "dlarzb.f"
		c__[i__ + j * c_dim1] -= work[i__ + j * work_dim1];
#line 307 "dlarzb.f"
/* L50: */
#line 307 "dlarzb.f"
	    }
#line 308 "dlarzb.f"
/* L60: */
#line 308 "dlarzb.f"
	}

/*        C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ... */
/*                            W( 1:m, 1:k ) * V( 1:k, 1:l ) */

#line 313 "dlarzb.f"
	if (*l > 0) {
#line 313 "dlarzb.f"
	    dgemm_("No transpose", "No transpose", m, l, k, &c_b23, &work[
		    work_offset], ldwork, &v[v_offset], ldv, &c_b13, &c__[(*n 
		    - *l + 1) * c_dim1 + 1], ldc, (ftnlen)12, (ftnlen)12);
#line 313 "dlarzb.f"
	}

#line 317 "dlarzb.f"
    }

#line 319 "dlarzb.f"
    return 0;

/*     End of DLARZB */

} /* dlarzb_ */

