#line 1 "sgbbrd.f"
/* sgbbrd.f -- translated by f2c (version 20100827).
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

#line 1 "sgbbrd.f"
/* Table of constant values */

static doublereal c_b8 = 0.;
static doublereal c_b9 = 1.;
static integer c__1 = 1;

/* > \brief \b SGBBRD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGBBRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgbbrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgbbrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgbbrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGBBRD( VECT, M, N, NCC, KL, KU, AB, LDAB, D, E, Q, */
/*                          LDQ, PT, LDPT, C, LDC, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          VECT */
/*       INTEGER            INFO, KL, KU, LDAB, LDC, LDPT, LDQ, M, N, NCC */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AB( LDAB, * ), C( LDC, * ), D( * ), E( * ), */
/*      $                   PT( LDPT, * ), Q( LDQ, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGBBRD reduces a real general m-by-n band matrix A to upper */
/* > bidiagonal form B by an orthogonal transformation: Q**T * A * P = B. */
/* > */
/* > The routine computes B, and optionally forms Q or P**T, or computes */
/* > Q**T*C for a given matrix C. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] VECT */
/* > \verbatim */
/* >          VECT is CHARACTER*1 */
/* >          Specifies whether or not the matrices Q and P**T are to be */
/* >          formed. */
/* >          = 'N': do not form Q or P**T; */
/* >          = 'Q': form Q only; */
/* >          = 'P': form P**T only; */
/* >          = 'B': form both. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NCC */
/* > \verbatim */
/* >          NCC is INTEGER */
/* >          The number of columns of the matrix C.  NCC >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KL */
/* > \verbatim */
/* >          KL is INTEGER */
/* >          The number of subdiagonals of the matrix A. KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* >          KU is INTEGER */
/* >          The number of superdiagonals of the matrix A. KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is REAL array, dimension (LDAB,N) */
/* >          On entry, the m-by-n band matrix A, stored in rows 1 to */
/* >          KL+KU+1. The j-th column of A is stored in the j-th column of */
/* >          the array AB as follows: */
/* >          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl). */
/* >          On exit, A is overwritten by values generated during the */
/* >          reduction. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array A. LDAB >= KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* >          D is REAL array, dimension (min(M,N)) */
/* >          The diagonal elements of the bidiagonal matrix B. */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (min(M,N)-1) */
/* >          The superdiagonal elements of the bidiagonal matrix B. */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* >          Q is REAL array, dimension (LDQ,M) */
/* >          If VECT = 'Q' or 'B', the m-by-m orthogonal matrix Q. */
/* >          If VECT = 'N' or 'P', the array Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q. */
/* >          LDQ >= max(1,M) if VECT = 'Q' or 'B'; LDQ >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[out] PT */
/* > \verbatim */
/* >          PT is REAL array, dimension (LDPT,N) */
/* >          If VECT = 'P' or 'B', the n-by-n orthogonal matrix P'. */
/* >          If VECT = 'N' or 'Q', the array PT is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDPT */
/* > \verbatim */
/* >          LDPT is INTEGER */
/* >          The leading dimension of the array PT. */
/* >          LDPT >= max(1,N) if VECT = 'P' or 'B'; LDPT >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array, dimension (LDC,NCC) */
/* >          On entry, an m-by-ncc matrix C. */
/* >          On exit, C is overwritten by Q**T*C. */
/* >          C is not referenced if NCC = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. */
/* >          LDC >= max(1,M) if NCC > 0; LDC >= 1 if NCC = 0. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (2*max(M,N)) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup realGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int sgbbrd_(char *vect, integer *m, integer *n, integer *ncc,
	 integer *kl, integer *ku, doublereal *ab, integer *ldab, doublereal *
	d__, doublereal *e, doublereal *q, integer *ldq, doublereal *pt, 
	integer *ldpt, doublereal *c__, integer *ldc, doublereal *work, 
	integer *info, ftnlen vect_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, c_dim1, c_offset, pt_dim1, pt_offset, q_dim1, 
	    q_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;

    /* Local variables */
    static integer i__, j, l, j1, j2, kb;
    static doublereal ra, rb, rc;
    static integer kk, ml, mn, nr, mu;
    static doublereal rs;
    static integer kb1, ml0, mu0, klm, kun, nrt, klu1, inca;
    extern /* Subroutine */ int srot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical wantb, wantc;
    static integer minmn;
    static logical wantq;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), slaset_(
	    char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), slartg_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), slargv_(
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *), slartv_(integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *);
    static logical wantpt;


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 230 "sgbbrd.f"
    /* Parameter adjustments */
#line 230 "sgbbrd.f"
    ab_dim1 = *ldab;
#line 230 "sgbbrd.f"
    ab_offset = 1 + ab_dim1;
#line 230 "sgbbrd.f"
    ab -= ab_offset;
#line 230 "sgbbrd.f"
    --d__;
#line 230 "sgbbrd.f"
    --e;
#line 230 "sgbbrd.f"
    q_dim1 = *ldq;
#line 230 "sgbbrd.f"
    q_offset = 1 + q_dim1;
#line 230 "sgbbrd.f"
    q -= q_offset;
#line 230 "sgbbrd.f"
    pt_dim1 = *ldpt;
#line 230 "sgbbrd.f"
    pt_offset = 1 + pt_dim1;
#line 230 "sgbbrd.f"
    pt -= pt_offset;
#line 230 "sgbbrd.f"
    c_dim1 = *ldc;
#line 230 "sgbbrd.f"
    c_offset = 1 + c_dim1;
#line 230 "sgbbrd.f"
    c__ -= c_offset;
#line 230 "sgbbrd.f"
    --work;
#line 230 "sgbbrd.f"

#line 230 "sgbbrd.f"
    /* Function Body */
#line 230 "sgbbrd.f"
    wantb = lsame_(vect, "B", (ftnlen)1, (ftnlen)1);
#line 231 "sgbbrd.f"
    wantq = lsame_(vect, "Q", (ftnlen)1, (ftnlen)1) || wantb;
#line 232 "sgbbrd.f"
    wantpt = lsame_(vect, "P", (ftnlen)1, (ftnlen)1) || wantb;
#line 233 "sgbbrd.f"
    wantc = *ncc > 0;
#line 234 "sgbbrd.f"
    klu1 = *kl + *ku + 1;
#line 235 "sgbbrd.f"
    *info = 0;
#line 236 "sgbbrd.f"
    if (! wantq && ! wantpt && ! lsame_(vect, "N", (ftnlen)1, (ftnlen)1)) {
#line 238 "sgbbrd.f"
	*info = -1;
#line 239 "sgbbrd.f"
    } else if (*m < 0) {
#line 240 "sgbbrd.f"
	*info = -2;
#line 241 "sgbbrd.f"
    } else if (*n < 0) {
#line 242 "sgbbrd.f"
	*info = -3;
#line 243 "sgbbrd.f"
    } else if (*ncc < 0) {
#line 244 "sgbbrd.f"
	*info = -4;
#line 245 "sgbbrd.f"
    } else if (*kl < 0) {
#line 246 "sgbbrd.f"
	*info = -5;
#line 247 "sgbbrd.f"
    } else if (*ku < 0) {
#line 248 "sgbbrd.f"
	*info = -6;
#line 249 "sgbbrd.f"
    } else if (*ldab < klu1) {
#line 250 "sgbbrd.f"
	*info = -8;
#line 251 "sgbbrd.f"
    } else if (*ldq < 1 || wantq && *ldq < max(1,*m)) {
#line 252 "sgbbrd.f"
	*info = -12;
#line 253 "sgbbrd.f"
    } else if (*ldpt < 1 || wantpt && *ldpt < max(1,*n)) {
#line 254 "sgbbrd.f"
	*info = -14;
#line 255 "sgbbrd.f"
    } else if (*ldc < 1 || wantc && *ldc < max(1,*m)) {
#line 256 "sgbbrd.f"
	*info = -16;
#line 257 "sgbbrd.f"
    }
#line 258 "sgbbrd.f"
    if (*info != 0) {
#line 259 "sgbbrd.f"
	i__1 = -(*info);
#line 259 "sgbbrd.f"
	xerbla_("SGBBRD", &i__1, (ftnlen)6);
#line 260 "sgbbrd.f"
	return 0;
#line 261 "sgbbrd.f"
    }

/*     Initialize Q and P**T to the unit matrix, if needed */

#line 265 "sgbbrd.f"
    if (wantq) {
#line 265 "sgbbrd.f"
	slaset_("Full", m, m, &c_b8, &c_b9, &q[q_offset], ldq, (ftnlen)4);
#line 265 "sgbbrd.f"
    }
#line 267 "sgbbrd.f"
    if (wantpt) {
#line 267 "sgbbrd.f"
	slaset_("Full", n, n, &c_b8, &c_b9, &pt[pt_offset], ldpt, (ftnlen)4);
#line 267 "sgbbrd.f"
    }

/*     Quick return if possible. */

#line 272 "sgbbrd.f"
    if (*m == 0 || *n == 0) {
#line 272 "sgbbrd.f"
	return 0;
#line 272 "sgbbrd.f"
    }

#line 275 "sgbbrd.f"
    minmn = min(*m,*n);

#line 277 "sgbbrd.f"
    if (*kl + *ku > 1) {

/*        Reduce to upper bidiagonal form if KU > 0; if KU = 0, reduce */
/*        first to lower bidiagonal form and then transform to upper */
/*        bidiagonal */

#line 283 "sgbbrd.f"
	if (*ku > 0) {
#line 284 "sgbbrd.f"
	    ml0 = 1;
#line 285 "sgbbrd.f"
	    mu0 = 2;
#line 286 "sgbbrd.f"
	} else {
#line 287 "sgbbrd.f"
	    ml0 = 2;
#line 288 "sgbbrd.f"
	    mu0 = 1;
#line 289 "sgbbrd.f"
	}

/*        Wherever possible, plane rotations are generated and applied in */
/*        vector operations of length NR over the index set J1:J2:KLU1. */

/*        The sines of the plane rotations are stored in WORK(1:max(m,n)) */
/*        and the cosines in WORK(max(m,n)+1:2*max(m,n)). */

#line 297 "sgbbrd.f"
	mn = max(*m,*n);
/* Computing MIN */
#line 298 "sgbbrd.f"
	i__1 = *m - 1;
#line 298 "sgbbrd.f"
	klm = min(i__1,*kl);
/* Computing MIN */
#line 299 "sgbbrd.f"
	i__1 = *n - 1;
#line 299 "sgbbrd.f"
	kun = min(i__1,*ku);
#line 300 "sgbbrd.f"
	kb = klm + kun;
#line 301 "sgbbrd.f"
	kb1 = kb + 1;
#line 302 "sgbbrd.f"
	inca = kb1 * *ldab;
#line 303 "sgbbrd.f"
	nr = 0;
#line 304 "sgbbrd.f"
	j1 = klm + 2;
#line 305 "sgbbrd.f"
	j2 = 1 - kun;

#line 307 "sgbbrd.f"
	i__1 = minmn;
#line 307 "sgbbrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Reduce i-th column and i-th row of matrix to bidiagonal form */

#line 311 "sgbbrd.f"
	    ml = klm + 1;
#line 312 "sgbbrd.f"
	    mu = kun + 1;
#line 313 "sgbbrd.f"
	    i__2 = kb;
#line 313 "sgbbrd.f"
	    for (kk = 1; kk <= i__2; ++kk) {
#line 314 "sgbbrd.f"
		j1 += kb;
#line 315 "sgbbrd.f"
		j2 += kb;

/*              generate plane rotations to annihilate nonzero elements */
/*              which have been created below the band */

#line 320 "sgbbrd.f"
		if (nr > 0) {
#line 320 "sgbbrd.f"
		    slargv_(&nr, &ab[klu1 + (j1 - klm - 1) * ab_dim1], &inca, 
			    &work[j1], &kb1, &work[mn + j1], &kb1);
#line 320 "sgbbrd.f"
		}

/*              apply plane rotations from the left */

#line 326 "sgbbrd.f"
		i__3 = kb;
#line 326 "sgbbrd.f"
		for (l = 1; l <= i__3; ++l) {
#line 327 "sgbbrd.f"
		    if (j2 - klm + l - 1 > *n) {
#line 328 "sgbbrd.f"
			nrt = nr - 1;
#line 329 "sgbbrd.f"
		    } else {
#line 330 "sgbbrd.f"
			nrt = nr;
#line 331 "sgbbrd.f"
		    }
#line 332 "sgbbrd.f"
		    if (nrt > 0) {
#line 332 "sgbbrd.f"
			slartv_(&nrt, &ab[klu1 - l + (j1 - klm + l - 1) * 
				ab_dim1], &inca, &ab[klu1 - l + 1 + (j1 - klm 
				+ l - 1) * ab_dim1], &inca, &work[mn + j1], &
				work[j1], &kb1);
#line 332 "sgbbrd.f"
		    }
#line 336 "sgbbrd.f"
/* L10: */
#line 336 "sgbbrd.f"
		}

#line 338 "sgbbrd.f"
		if (ml > ml0) {
#line 339 "sgbbrd.f"
		    if (ml <= *m - i__ + 1) {

/*                    generate plane rotation to annihilate a(i+ml-1,i) */
/*                    within the band, and apply rotation from the left */

#line 344 "sgbbrd.f"
			slartg_(&ab[*ku + ml - 1 + i__ * ab_dim1], &ab[*ku + 
				ml + i__ * ab_dim1], &work[mn + i__ + ml - 1],
				 &work[i__ + ml - 1], &ra);
#line 347 "sgbbrd.f"
			ab[*ku + ml - 1 + i__ * ab_dim1] = ra;
#line 348 "sgbbrd.f"
			if (i__ < *n) {
/* Computing MIN */
#line 348 "sgbbrd.f"
			    i__4 = *ku + ml - 2, i__5 = *n - i__;
#line 348 "sgbbrd.f"
			    i__3 = min(i__4,i__5);
#line 348 "sgbbrd.f"
			    i__6 = *ldab - 1;
#line 348 "sgbbrd.f"
			    i__7 = *ldab - 1;
#line 348 "sgbbrd.f"
			    srot_(&i__3, &ab[*ku + ml - 2 + (i__ + 1) * 
				    ab_dim1], &i__6, &ab[*ku + ml - 1 + (i__ 
				    + 1) * ab_dim1], &i__7, &work[mn + i__ + 
				    ml - 1], &work[i__ + ml - 1]);
#line 348 "sgbbrd.f"
			}
#line 353 "sgbbrd.f"
		    }
#line 354 "sgbbrd.f"
		    ++nr;
#line 355 "sgbbrd.f"
		    j1 -= kb1;
#line 356 "sgbbrd.f"
		}

#line 358 "sgbbrd.f"
		if (wantq) {

/*                 accumulate product of plane rotations in Q */

#line 362 "sgbbrd.f"
		    i__3 = j2;
#line 362 "sgbbrd.f"
		    i__4 = kb1;
#line 362 "sgbbrd.f"
		    for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) 
			    {
#line 363 "sgbbrd.f"
			srot_(m, &q[(j - 1) * q_dim1 + 1], &c__1, &q[j * 
				q_dim1 + 1], &c__1, &work[mn + j], &work[j]);
#line 365 "sgbbrd.f"
/* L20: */
#line 365 "sgbbrd.f"
		    }
#line 366 "sgbbrd.f"
		}

#line 368 "sgbbrd.f"
		if (wantc) {

/*                 apply plane rotations to C */

#line 372 "sgbbrd.f"
		    i__4 = j2;
#line 372 "sgbbrd.f"
		    i__3 = kb1;
#line 372 "sgbbrd.f"
		    for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) 
			    {
#line 373 "sgbbrd.f"
			srot_(ncc, &c__[j - 1 + c_dim1], ldc, &c__[j + c_dim1]
				, ldc, &work[mn + j], &work[j]);
#line 375 "sgbbrd.f"
/* L30: */
#line 375 "sgbbrd.f"
		    }
#line 376 "sgbbrd.f"
		}

#line 378 "sgbbrd.f"
		if (j2 + kun > *n) {

/*                 adjust J2 to keep within the bounds of the matrix */

#line 382 "sgbbrd.f"
		    --nr;
#line 383 "sgbbrd.f"
		    j2 -= kb1;
#line 384 "sgbbrd.f"
		}

#line 386 "sgbbrd.f"
		i__3 = j2;
#line 386 "sgbbrd.f"
		i__4 = kb1;
#line 386 "sgbbrd.f"
		for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {

/*                 create nonzero element a(j-1,j+ku) above the band */
/*                 and store it in WORK(n+1:2*n) */

#line 391 "sgbbrd.f"
		    work[j + kun] = work[j] * ab[(j + kun) * ab_dim1 + 1];
#line 392 "sgbbrd.f"
		    ab[(j + kun) * ab_dim1 + 1] = work[mn + j] * ab[(j + kun) 
			    * ab_dim1 + 1];
#line 393 "sgbbrd.f"
/* L40: */
#line 393 "sgbbrd.f"
		}

/*              generate plane rotations to annihilate nonzero elements */
/*              which have been generated above the band */

#line 398 "sgbbrd.f"
		if (nr > 0) {
#line 398 "sgbbrd.f"
		    slargv_(&nr, &ab[(j1 + kun - 1) * ab_dim1 + 1], &inca, &
			    work[j1 + kun], &kb1, &work[mn + j1 + kun], &kb1);
#line 398 "sgbbrd.f"
		}

/*              apply plane rotations from the right */

#line 405 "sgbbrd.f"
		i__4 = kb;
#line 405 "sgbbrd.f"
		for (l = 1; l <= i__4; ++l) {
#line 406 "sgbbrd.f"
		    if (j2 + l - 1 > *m) {
#line 407 "sgbbrd.f"
			nrt = nr - 1;
#line 408 "sgbbrd.f"
		    } else {
#line 409 "sgbbrd.f"
			nrt = nr;
#line 410 "sgbbrd.f"
		    }
#line 411 "sgbbrd.f"
		    if (nrt > 0) {
#line 411 "sgbbrd.f"
			slartv_(&nrt, &ab[l + 1 + (j1 + kun - 1) * ab_dim1], &
				inca, &ab[l + (j1 + kun) * ab_dim1], &inca, &
				work[mn + j1 + kun], &work[j1 + kun], &kb1);
#line 411 "sgbbrd.f"
		    }
#line 416 "sgbbrd.f"
/* L50: */
#line 416 "sgbbrd.f"
		}

#line 418 "sgbbrd.f"
		if (ml == ml0 && mu > mu0) {
#line 419 "sgbbrd.f"
		    if (mu <= *n - i__ + 1) {

/*                    generate plane rotation to annihilate a(i,i+mu-1) */
/*                    within the band, and apply rotation from the right */

#line 424 "sgbbrd.f"
			slartg_(&ab[*ku - mu + 3 + (i__ + mu - 2) * ab_dim1], 
				&ab[*ku - mu + 2 + (i__ + mu - 1) * ab_dim1], 
				&work[mn + i__ + mu - 1], &work[i__ + mu - 1],
				 &ra);
#line 428 "sgbbrd.f"
			ab[*ku - mu + 3 + (i__ + mu - 2) * ab_dim1] = ra;
/* Computing MIN */
#line 429 "sgbbrd.f"
			i__3 = *kl + mu - 2, i__5 = *m - i__;
#line 429 "sgbbrd.f"
			i__4 = min(i__3,i__5);
#line 429 "sgbbrd.f"
			srot_(&i__4, &ab[*ku - mu + 4 + (i__ + mu - 2) * 
				ab_dim1], &c__1, &ab[*ku - mu + 3 + (i__ + mu 
				- 1) * ab_dim1], &c__1, &work[mn + i__ + mu - 
				1], &work[i__ + mu - 1]);
#line 433 "sgbbrd.f"
		    }
#line 434 "sgbbrd.f"
		    ++nr;
#line 435 "sgbbrd.f"
		    j1 -= kb1;
#line 436 "sgbbrd.f"
		}

#line 438 "sgbbrd.f"
		if (wantpt) {

/*                 accumulate product of plane rotations in P**T */

#line 442 "sgbbrd.f"
		    i__4 = j2;
#line 442 "sgbbrd.f"
		    i__3 = kb1;
#line 442 "sgbbrd.f"
		    for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) 
			    {
#line 443 "sgbbrd.f"
			srot_(n, &pt[j + kun - 1 + pt_dim1], ldpt, &pt[j + 
				kun + pt_dim1], ldpt, &work[mn + j + kun], &
				work[j + kun]);
#line 446 "sgbbrd.f"
/* L60: */
#line 446 "sgbbrd.f"
		    }
#line 447 "sgbbrd.f"
		}

#line 449 "sgbbrd.f"
		if (j2 + kb > *m) {

/*                 adjust J2 to keep within the bounds of the matrix */

#line 453 "sgbbrd.f"
		    --nr;
#line 454 "sgbbrd.f"
		    j2 -= kb1;
#line 455 "sgbbrd.f"
		}

#line 457 "sgbbrd.f"
		i__3 = j2;
#line 457 "sgbbrd.f"
		i__4 = kb1;
#line 457 "sgbbrd.f"
		for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {

/*                 create nonzero element a(j+kl+ku,j+ku-1) below the */
/*                 band and store it in WORK(1:n) */

#line 462 "sgbbrd.f"
		    work[j + kb] = work[j + kun] * ab[klu1 + (j + kun) * 
			    ab_dim1];
#line 463 "sgbbrd.f"
		    ab[klu1 + (j + kun) * ab_dim1] = work[mn + j + kun] * ab[
			    klu1 + (j + kun) * ab_dim1];
#line 464 "sgbbrd.f"
/* L70: */
#line 464 "sgbbrd.f"
		}

#line 466 "sgbbrd.f"
		if (ml > ml0) {
#line 467 "sgbbrd.f"
		    --ml;
#line 468 "sgbbrd.f"
		} else {
#line 469 "sgbbrd.f"
		    --mu;
#line 470 "sgbbrd.f"
		}
#line 471 "sgbbrd.f"
/* L80: */
#line 471 "sgbbrd.f"
	    }
#line 472 "sgbbrd.f"
/* L90: */
#line 472 "sgbbrd.f"
	}
#line 473 "sgbbrd.f"
    }

#line 475 "sgbbrd.f"
    if (*ku == 0 && *kl > 0) {

/*        A has been reduced to lower bidiagonal form */

/*        Transform lower bidiagonal form to upper bidiagonal by applying */
/*        plane rotations from the left, storing diagonal elements in D */
/*        and off-diagonal elements in E */

/* Computing MIN */
#line 483 "sgbbrd.f"
	i__2 = *m - 1;
#line 483 "sgbbrd.f"
	i__1 = min(i__2,*n);
#line 483 "sgbbrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 484 "sgbbrd.f"
	    slartg_(&ab[i__ * ab_dim1 + 1], &ab[i__ * ab_dim1 + 2], &rc, &rs, 
		    &ra);
#line 485 "sgbbrd.f"
	    d__[i__] = ra;
#line 486 "sgbbrd.f"
	    if (i__ < *n) {
#line 487 "sgbbrd.f"
		e[i__] = rs * ab[(i__ + 1) * ab_dim1 + 1];
#line 488 "sgbbrd.f"
		ab[(i__ + 1) * ab_dim1 + 1] = rc * ab[(i__ + 1) * ab_dim1 + 1]
			;
#line 489 "sgbbrd.f"
	    }
#line 490 "sgbbrd.f"
	    if (wantq) {
#line 490 "sgbbrd.f"
		srot_(m, &q[i__ * q_dim1 + 1], &c__1, &q[(i__ + 1) * q_dim1 + 
			1], &c__1, &rc, &rs);
#line 490 "sgbbrd.f"
	    }
#line 492 "sgbbrd.f"
	    if (wantc) {
#line 492 "sgbbrd.f"
		srot_(ncc, &c__[i__ + c_dim1], ldc, &c__[i__ + 1 + c_dim1], 
			ldc, &rc, &rs);
#line 492 "sgbbrd.f"
	    }
#line 495 "sgbbrd.f"
/* L100: */
#line 495 "sgbbrd.f"
	}
#line 496 "sgbbrd.f"
	if (*m <= *n) {
#line 496 "sgbbrd.f"
	    d__[*m] = ab[*m * ab_dim1 + 1];
#line 496 "sgbbrd.f"
	}
#line 498 "sgbbrd.f"
    } else if (*ku > 0) {

/*        A has been reduced to upper bidiagonal form */

#line 502 "sgbbrd.f"
	if (*m < *n) {

/*           Annihilate a(m,m+1) by applying plane rotations from the */
/*           right, storing diagonal elements in D and off-diagonal */
/*           elements in E */

#line 508 "sgbbrd.f"
	    rb = ab[*ku + (*m + 1) * ab_dim1];
#line 509 "sgbbrd.f"
	    for (i__ = *m; i__ >= 1; --i__) {
#line 510 "sgbbrd.f"
		slartg_(&ab[*ku + 1 + i__ * ab_dim1], &rb, &rc, &rs, &ra);
#line 511 "sgbbrd.f"
		d__[i__] = ra;
#line 512 "sgbbrd.f"
		if (i__ > 1) {
#line 513 "sgbbrd.f"
		    rb = -rs * ab[*ku + i__ * ab_dim1];
#line 514 "sgbbrd.f"
		    e[i__ - 1] = rc * ab[*ku + i__ * ab_dim1];
#line 515 "sgbbrd.f"
		}
#line 516 "sgbbrd.f"
		if (wantpt) {
#line 516 "sgbbrd.f"
		    srot_(n, &pt[i__ + pt_dim1], ldpt, &pt[*m + 1 + pt_dim1], 
			    ldpt, &rc, &rs);
#line 516 "sgbbrd.f"
		}
#line 519 "sgbbrd.f"
/* L110: */
#line 519 "sgbbrd.f"
	    }
#line 520 "sgbbrd.f"
	} else {

/*           Copy off-diagonal elements to E and diagonal elements to D */

#line 524 "sgbbrd.f"
	    i__1 = minmn - 1;
#line 524 "sgbbrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 525 "sgbbrd.f"
		e[i__] = ab[*ku + (i__ + 1) * ab_dim1];
#line 526 "sgbbrd.f"
/* L120: */
#line 526 "sgbbrd.f"
	    }
#line 527 "sgbbrd.f"
	    i__1 = minmn;
#line 527 "sgbbrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 528 "sgbbrd.f"
		d__[i__] = ab[*ku + 1 + i__ * ab_dim1];
#line 529 "sgbbrd.f"
/* L130: */
#line 529 "sgbbrd.f"
	    }
#line 530 "sgbbrd.f"
	}
#line 531 "sgbbrd.f"
    } else {

/*        A is diagonal. Set elements of E to zero and copy diagonal */
/*        elements to D. */

#line 536 "sgbbrd.f"
	i__1 = minmn - 1;
#line 536 "sgbbrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 537 "sgbbrd.f"
	    e[i__] = 0.;
#line 538 "sgbbrd.f"
/* L140: */
#line 538 "sgbbrd.f"
	}
#line 539 "sgbbrd.f"
	i__1 = minmn;
#line 539 "sgbbrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 540 "sgbbrd.f"
	    d__[i__] = ab[i__ * ab_dim1 + 1];
#line 541 "sgbbrd.f"
/* L150: */
#line 541 "sgbbrd.f"
	}
#line 542 "sgbbrd.f"
    }
#line 543 "sgbbrd.f"
    return 0;

/*     End of SGBBRD */

} /* sgbbrd_ */

