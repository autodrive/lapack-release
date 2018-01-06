#line 1 "zgbbrd.f"
/* zgbbrd.f -- translated by f2c (version 20100827).
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

#line 1 "zgbbrd.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZGBBRD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGBBRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgbbrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgbbrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgbbrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGBBRD( VECT, M, N, NCC, KL, KU, AB, LDAB, D, E, Q, */
/*                          LDQ, PT, LDPT, C, LDC, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          VECT */
/*       INTEGER            INFO, KL, KU, LDAB, LDC, LDPT, LDQ, M, N, NCC */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), E( * ), RWORK( * ) */
/*       COMPLEX*16         AB( LDAB, * ), C( LDC, * ), PT( LDPT, * ), */
/*      $                   Q( LDQ, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGBBRD reduces a complex general m-by-n band matrix A to real upper */
/* > bidiagonal form B by a unitary transformation: Q**H * A * P = B. */
/* > */
/* > The routine computes B, and optionally forms Q or P**H, or computes */
/* > Q**H*C for a given matrix C. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] VECT */
/* > \verbatim */
/* >          VECT is CHARACTER*1 */
/* >          Specifies whether or not the matrices Q and P**H are to be */
/* >          formed. */
/* >          = 'N': do not form Q or P**H; */
/* >          = 'Q': form Q only; */
/* >          = 'P': form P**H only; */
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
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
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
/* >          D is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The diagonal elements of the bidiagonal matrix B. */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (min(M,N)-1) */
/* >          The superdiagonal elements of the bidiagonal matrix B. */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* >          Q is COMPLEX*16 array, dimension (LDQ,M) */
/* >          If VECT = 'Q' or 'B', the m-by-m unitary matrix Q. */
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
/* >          PT is COMPLEX*16 array, dimension (LDPT,N) */
/* >          If VECT = 'P' or 'B', the n-by-n unitary matrix P'. */
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
/* >          C is COMPLEX*16 array, dimension (LDC,NCC) */
/* >          On entry, an m-by-ncc matrix C. */
/* >          On exit, C is overwritten by Q**H*C. */
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
/* >          WORK is COMPLEX*16 array, dimension (max(M,N)) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (max(M,N)) */
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

/* > \date December 2016 */

/* > \ingroup complex16GBcomputational */

/*  ===================================================================== */
/* Subroutine */ int zgbbrd_(char *vect, integer *m, integer *n, integer *ncc,
	 integer *kl, integer *ku, doublecomplex *ab, integer *ldab, 
	doublereal *d__, doublereal *e, doublecomplex *q, integer *ldq, 
	doublecomplex *pt, integer *ldpt, doublecomplex *c__, integer *ldc, 
	doublecomplex *work, doublereal *rwork, integer *info, ftnlen 
	vect_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, c_dim1, c_offset, pt_dim1, pt_offset, q_dim1, 
	    q_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, j, l;
    static doublecomplex t;
    static integer j1, j2, kb;
    static doublecomplex ra, rb;
    static doublereal rc;
    static integer kk, ml, nr, mu;
    static doublecomplex rs;
    static integer kb1, ml0, mu0, klm, kun, nrt, klu1, inca;
    static doublereal abst;
    extern /* Subroutine */ int zrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical wantb, wantc;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static integer minmn;
    static logical wantq;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zlaset_(
	    char *, integer *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, ftnlen), zlartg_(doublecomplex *, 
	    doublecomplex *, doublereal *, doublecomplex *, doublecomplex *), 
	    zlargv_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublereal *, integer *);
    static logical wantpt;
    extern /* Subroutine */ int zlartv_(integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublereal *, doublecomplex *, 
	    integer *);


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 242 "zgbbrd.f"
    /* Parameter adjustments */
#line 242 "zgbbrd.f"
    ab_dim1 = *ldab;
#line 242 "zgbbrd.f"
    ab_offset = 1 + ab_dim1;
#line 242 "zgbbrd.f"
    ab -= ab_offset;
#line 242 "zgbbrd.f"
    --d__;
#line 242 "zgbbrd.f"
    --e;
#line 242 "zgbbrd.f"
    q_dim1 = *ldq;
#line 242 "zgbbrd.f"
    q_offset = 1 + q_dim1;
#line 242 "zgbbrd.f"
    q -= q_offset;
#line 242 "zgbbrd.f"
    pt_dim1 = *ldpt;
#line 242 "zgbbrd.f"
    pt_offset = 1 + pt_dim1;
#line 242 "zgbbrd.f"
    pt -= pt_offset;
#line 242 "zgbbrd.f"
    c_dim1 = *ldc;
#line 242 "zgbbrd.f"
    c_offset = 1 + c_dim1;
#line 242 "zgbbrd.f"
    c__ -= c_offset;
#line 242 "zgbbrd.f"
    --work;
#line 242 "zgbbrd.f"
    --rwork;
#line 242 "zgbbrd.f"

#line 242 "zgbbrd.f"
    /* Function Body */
#line 242 "zgbbrd.f"
    wantb = lsame_(vect, "B", (ftnlen)1, (ftnlen)1);
#line 243 "zgbbrd.f"
    wantq = lsame_(vect, "Q", (ftnlen)1, (ftnlen)1) || wantb;
#line 244 "zgbbrd.f"
    wantpt = lsame_(vect, "P", (ftnlen)1, (ftnlen)1) || wantb;
#line 245 "zgbbrd.f"
    wantc = *ncc > 0;
#line 246 "zgbbrd.f"
    klu1 = *kl + *ku + 1;
#line 247 "zgbbrd.f"
    *info = 0;
#line 248 "zgbbrd.f"
    if (! wantq && ! wantpt && ! lsame_(vect, "N", (ftnlen)1, (ftnlen)1)) {
#line 250 "zgbbrd.f"
	*info = -1;
#line 251 "zgbbrd.f"
    } else if (*m < 0) {
#line 252 "zgbbrd.f"
	*info = -2;
#line 253 "zgbbrd.f"
    } else if (*n < 0) {
#line 254 "zgbbrd.f"
	*info = -3;
#line 255 "zgbbrd.f"
    } else if (*ncc < 0) {
#line 256 "zgbbrd.f"
	*info = -4;
#line 257 "zgbbrd.f"
    } else if (*kl < 0) {
#line 258 "zgbbrd.f"
	*info = -5;
#line 259 "zgbbrd.f"
    } else if (*ku < 0) {
#line 260 "zgbbrd.f"
	*info = -6;
#line 261 "zgbbrd.f"
    } else if (*ldab < klu1) {
#line 262 "zgbbrd.f"
	*info = -8;
#line 263 "zgbbrd.f"
    } else if (*ldq < 1 || wantq && *ldq < max(1,*m)) {
#line 264 "zgbbrd.f"
	*info = -12;
#line 265 "zgbbrd.f"
    } else if (*ldpt < 1 || wantpt && *ldpt < max(1,*n)) {
#line 266 "zgbbrd.f"
	*info = -14;
#line 267 "zgbbrd.f"
    } else if (*ldc < 1 || wantc && *ldc < max(1,*m)) {
#line 268 "zgbbrd.f"
	*info = -16;
#line 269 "zgbbrd.f"
    }
#line 270 "zgbbrd.f"
    if (*info != 0) {
#line 271 "zgbbrd.f"
	i__1 = -(*info);
#line 271 "zgbbrd.f"
	xerbla_("ZGBBRD", &i__1, (ftnlen)6);
#line 272 "zgbbrd.f"
	return 0;
#line 273 "zgbbrd.f"
    }

/*     Initialize Q and P**H to the unit matrix, if needed */

#line 277 "zgbbrd.f"
    if (wantq) {
#line 277 "zgbbrd.f"
	zlaset_("Full", m, m, &c_b1, &c_b2, &q[q_offset], ldq, (ftnlen)4);
#line 277 "zgbbrd.f"
    }
#line 279 "zgbbrd.f"
    if (wantpt) {
#line 279 "zgbbrd.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &pt[pt_offset], ldpt, (ftnlen)4);
#line 279 "zgbbrd.f"
    }

/*     Quick return if possible. */

#line 284 "zgbbrd.f"
    if (*m == 0 || *n == 0) {
#line 284 "zgbbrd.f"
	return 0;
#line 284 "zgbbrd.f"
    }

#line 287 "zgbbrd.f"
    minmn = min(*m,*n);

#line 289 "zgbbrd.f"
    if (*kl + *ku > 1) {

/*        Reduce to upper bidiagonal form if KU > 0; if KU = 0, reduce */
/*        first to lower bidiagonal form and then transform to upper */
/*        bidiagonal */

#line 295 "zgbbrd.f"
	if (*ku > 0) {
#line 296 "zgbbrd.f"
	    ml0 = 1;
#line 297 "zgbbrd.f"
	    mu0 = 2;
#line 298 "zgbbrd.f"
	} else {
#line 299 "zgbbrd.f"
	    ml0 = 2;
#line 300 "zgbbrd.f"
	    mu0 = 1;
#line 301 "zgbbrd.f"
	}

/*        Wherever possible, plane rotations are generated and applied in */
/*        vector operations of length NR over the index set J1:J2:KLU1. */

/*        The complex sines of the plane rotations are stored in WORK, */
/*        and the real cosines in RWORK. */

/* Computing MIN */
#line 309 "zgbbrd.f"
	i__1 = *m - 1;
#line 309 "zgbbrd.f"
	klm = min(i__1,*kl);
/* Computing MIN */
#line 310 "zgbbrd.f"
	i__1 = *n - 1;
#line 310 "zgbbrd.f"
	kun = min(i__1,*ku);
#line 311 "zgbbrd.f"
	kb = klm + kun;
#line 312 "zgbbrd.f"
	kb1 = kb + 1;
#line 313 "zgbbrd.f"
	inca = kb1 * *ldab;
#line 314 "zgbbrd.f"
	nr = 0;
#line 315 "zgbbrd.f"
	j1 = klm + 2;
#line 316 "zgbbrd.f"
	j2 = 1 - kun;

#line 318 "zgbbrd.f"
	i__1 = minmn;
#line 318 "zgbbrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Reduce i-th column and i-th row of matrix to bidiagonal form */

#line 322 "zgbbrd.f"
	    ml = klm + 1;
#line 323 "zgbbrd.f"
	    mu = kun + 1;
#line 324 "zgbbrd.f"
	    i__2 = kb;
#line 324 "zgbbrd.f"
	    for (kk = 1; kk <= i__2; ++kk) {
#line 325 "zgbbrd.f"
		j1 += kb;
#line 326 "zgbbrd.f"
		j2 += kb;

/*              generate plane rotations to annihilate nonzero elements */
/*              which have been created below the band */

#line 331 "zgbbrd.f"
		if (nr > 0) {
#line 331 "zgbbrd.f"
		    zlargv_(&nr, &ab[klu1 + (j1 - klm - 1) * ab_dim1], &inca, 
			    &work[j1], &kb1, &rwork[j1], &kb1);
#line 331 "zgbbrd.f"
		}

/*              apply plane rotations from the left */

#line 337 "zgbbrd.f"
		i__3 = kb;
#line 337 "zgbbrd.f"
		for (l = 1; l <= i__3; ++l) {
#line 338 "zgbbrd.f"
		    if (j2 - klm + l - 1 > *n) {
#line 339 "zgbbrd.f"
			nrt = nr - 1;
#line 340 "zgbbrd.f"
		    } else {
#line 341 "zgbbrd.f"
			nrt = nr;
#line 342 "zgbbrd.f"
		    }
#line 343 "zgbbrd.f"
		    if (nrt > 0) {
#line 343 "zgbbrd.f"
			zlartv_(&nrt, &ab[klu1 - l + (j1 - klm + l - 1) * 
				ab_dim1], &inca, &ab[klu1 - l + 1 + (j1 - klm 
				+ l - 1) * ab_dim1], &inca, &rwork[j1], &work[
				j1], &kb1);
#line 343 "zgbbrd.f"
		    }
#line 347 "zgbbrd.f"
/* L10: */
#line 347 "zgbbrd.f"
		}

#line 349 "zgbbrd.f"
		if (ml > ml0) {
#line 350 "zgbbrd.f"
		    if (ml <= *m - i__ + 1) {

/*                    generate plane rotation to annihilate a(i+ml-1,i) */
/*                    within the band, and apply rotation from the left */

#line 355 "zgbbrd.f"
			zlartg_(&ab[*ku + ml - 1 + i__ * ab_dim1], &ab[*ku + 
				ml + i__ * ab_dim1], &rwork[i__ + ml - 1], &
				work[i__ + ml - 1], &ra);
#line 357 "zgbbrd.f"
			i__3 = *ku + ml - 1 + i__ * ab_dim1;
#line 357 "zgbbrd.f"
			ab[i__3].r = ra.r, ab[i__3].i = ra.i;
#line 358 "zgbbrd.f"
			if (i__ < *n) {
/* Computing MIN */
#line 358 "zgbbrd.f"
			    i__4 = *ku + ml - 2, i__5 = *n - i__;
#line 358 "zgbbrd.f"
			    i__3 = min(i__4,i__5);
#line 358 "zgbbrd.f"
			    i__6 = *ldab - 1;
#line 358 "zgbbrd.f"
			    i__7 = *ldab - 1;
#line 358 "zgbbrd.f"
			    zrot_(&i__3, &ab[*ku + ml - 2 + (i__ + 1) * 
				    ab_dim1], &i__6, &ab[*ku + ml - 1 + (i__ 
				    + 1) * ab_dim1], &i__7, &rwork[i__ + ml - 
				    1], &work[i__ + ml - 1]);
#line 358 "zgbbrd.f"
			}
#line 363 "zgbbrd.f"
		    }
#line 364 "zgbbrd.f"
		    ++nr;
#line 365 "zgbbrd.f"
		    j1 -= kb1;
#line 366 "zgbbrd.f"
		}

#line 368 "zgbbrd.f"
		if (wantq) {

/*                 accumulate product of plane rotations in Q */

#line 372 "zgbbrd.f"
		    i__3 = j2;
#line 372 "zgbbrd.f"
		    i__4 = kb1;
#line 372 "zgbbrd.f"
		    for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) 
			    {
#line 373 "zgbbrd.f"
			d_cnjg(&z__1, &work[j]);
#line 373 "zgbbrd.f"
			zrot_(m, &q[(j - 1) * q_dim1 + 1], &c__1, &q[j * 
				q_dim1 + 1], &c__1, &rwork[j], &z__1);
#line 375 "zgbbrd.f"
/* L20: */
#line 375 "zgbbrd.f"
		    }
#line 376 "zgbbrd.f"
		}

#line 378 "zgbbrd.f"
		if (wantc) {

/*                 apply plane rotations to C */

#line 382 "zgbbrd.f"
		    i__4 = j2;
#line 382 "zgbbrd.f"
		    i__3 = kb1;
#line 382 "zgbbrd.f"
		    for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) 
			    {
#line 383 "zgbbrd.f"
			zrot_(ncc, &c__[j - 1 + c_dim1], ldc, &c__[j + c_dim1]
				, ldc, &rwork[j], &work[j]);
#line 385 "zgbbrd.f"
/* L30: */
#line 385 "zgbbrd.f"
		    }
#line 386 "zgbbrd.f"
		}

#line 388 "zgbbrd.f"
		if (j2 + kun > *n) {

/*                 adjust J2 to keep within the bounds of the matrix */

#line 392 "zgbbrd.f"
		    --nr;
#line 393 "zgbbrd.f"
		    j2 -= kb1;
#line 394 "zgbbrd.f"
		}

#line 396 "zgbbrd.f"
		i__3 = j2;
#line 396 "zgbbrd.f"
		i__4 = kb1;
#line 396 "zgbbrd.f"
		for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {

/*                 create nonzero element a(j-1,j+ku) above the band */
/*                 and store it in WORK(n+1:2*n) */

#line 401 "zgbbrd.f"
		    i__5 = j + kun;
#line 401 "zgbbrd.f"
		    i__6 = j;
#line 401 "zgbbrd.f"
		    i__7 = (j + kun) * ab_dim1 + 1;
#line 401 "zgbbrd.f"
		    z__1.r = work[i__6].r * ab[i__7].r - work[i__6].i * ab[
			    i__7].i, z__1.i = work[i__6].r * ab[i__7].i + 
			    work[i__6].i * ab[i__7].r;
#line 401 "zgbbrd.f"
		    work[i__5].r = z__1.r, work[i__5].i = z__1.i;
#line 402 "zgbbrd.f"
		    i__5 = (j + kun) * ab_dim1 + 1;
#line 402 "zgbbrd.f"
		    i__6 = j;
#line 402 "zgbbrd.f"
		    i__7 = (j + kun) * ab_dim1 + 1;
#line 402 "zgbbrd.f"
		    z__1.r = rwork[i__6] * ab[i__7].r, z__1.i = rwork[i__6] * 
			    ab[i__7].i;
#line 402 "zgbbrd.f"
		    ab[i__5].r = z__1.r, ab[i__5].i = z__1.i;
#line 403 "zgbbrd.f"
/* L40: */
#line 403 "zgbbrd.f"
		}

/*              generate plane rotations to annihilate nonzero elements */
/*              which have been generated above the band */

#line 408 "zgbbrd.f"
		if (nr > 0) {
#line 408 "zgbbrd.f"
		    zlargv_(&nr, &ab[(j1 + kun - 1) * ab_dim1 + 1], &inca, &
			    work[j1 + kun], &kb1, &rwork[j1 + kun], &kb1);
#line 408 "zgbbrd.f"
		}

/*              apply plane rotations from the right */

#line 415 "zgbbrd.f"
		i__4 = kb;
#line 415 "zgbbrd.f"
		for (l = 1; l <= i__4; ++l) {
#line 416 "zgbbrd.f"
		    if (j2 + l - 1 > *m) {
#line 417 "zgbbrd.f"
			nrt = nr - 1;
#line 418 "zgbbrd.f"
		    } else {
#line 419 "zgbbrd.f"
			nrt = nr;
#line 420 "zgbbrd.f"
		    }
#line 421 "zgbbrd.f"
		    if (nrt > 0) {
#line 421 "zgbbrd.f"
			zlartv_(&nrt, &ab[l + 1 + (j1 + kun - 1) * ab_dim1], &
				inca, &ab[l + (j1 + kun) * ab_dim1], &inca, &
				rwork[j1 + kun], &work[j1 + kun], &kb1);
#line 421 "zgbbrd.f"
		    }
#line 425 "zgbbrd.f"
/* L50: */
#line 425 "zgbbrd.f"
		}

#line 427 "zgbbrd.f"
		if (ml == ml0 && mu > mu0) {
#line 428 "zgbbrd.f"
		    if (mu <= *n - i__ + 1) {

/*                    generate plane rotation to annihilate a(i,i+mu-1) */
/*                    within the band, and apply rotation from the right */

#line 433 "zgbbrd.f"
			zlartg_(&ab[*ku - mu + 3 + (i__ + mu - 2) * ab_dim1], 
				&ab[*ku - mu + 2 + (i__ + mu - 1) * ab_dim1], 
				&rwork[i__ + mu - 1], &work[i__ + mu - 1], &
				ra);
#line 436 "zgbbrd.f"
			i__4 = *ku - mu + 3 + (i__ + mu - 2) * ab_dim1;
#line 436 "zgbbrd.f"
			ab[i__4].r = ra.r, ab[i__4].i = ra.i;
/* Computing MIN */
#line 437 "zgbbrd.f"
			i__3 = *kl + mu - 2, i__5 = *m - i__;
#line 437 "zgbbrd.f"
			i__4 = min(i__3,i__5);
#line 437 "zgbbrd.f"
			zrot_(&i__4, &ab[*ku - mu + 4 + (i__ + mu - 2) * 
				ab_dim1], &c__1, &ab[*ku - mu + 3 + (i__ + mu 
				- 1) * ab_dim1], &c__1, &rwork[i__ + mu - 1], 
				&work[i__ + mu - 1]);
#line 441 "zgbbrd.f"
		    }
#line 442 "zgbbrd.f"
		    ++nr;
#line 443 "zgbbrd.f"
		    j1 -= kb1;
#line 444 "zgbbrd.f"
		}

#line 446 "zgbbrd.f"
		if (wantpt) {

/*                 accumulate product of plane rotations in P**H */

#line 450 "zgbbrd.f"
		    i__4 = j2;
#line 450 "zgbbrd.f"
		    i__3 = kb1;
#line 450 "zgbbrd.f"
		    for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) 
			    {
#line 451 "zgbbrd.f"
			d_cnjg(&z__1, &work[j + kun]);
#line 451 "zgbbrd.f"
			zrot_(n, &pt[j + kun - 1 + pt_dim1], ldpt, &pt[j + 
				kun + pt_dim1], ldpt, &rwork[j + kun], &z__1);
#line 454 "zgbbrd.f"
/* L60: */
#line 454 "zgbbrd.f"
		    }
#line 455 "zgbbrd.f"
		}

#line 457 "zgbbrd.f"
		if (j2 + kb > *m) {

/*                 adjust J2 to keep within the bounds of the matrix */

#line 461 "zgbbrd.f"
		    --nr;
#line 462 "zgbbrd.f"
		    j2 -= kb1;
#line 463 "zgbbrd.f"
		}

#line 465 "zgbbrd.f"
		i__3 = j2;
#line 465 "zgbbrd.f"
		i__4 = kb1;
#line 465 "zgbbrd.f"
		for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {

/*                 create nonzero element a(j+kl+ku,j+ku-1) below the */
/*                 band and store it in WORK(1:n) */

#line 470 "zgbbrd.f"
		    i__5 = j + kb;
#line 470 "zgbbrd.f"
		    i__6 = j + kun;
#line 470 "zgbbrd.f"
		    i__7 = klu1 + (j + kun) * ab_dim1;
#line 470 "zgbbrd.f"
		    z__1.r = work[i__6].r * ab[i__7].r - work[i__6].i * ab[
			    i__7].i, z__1.i = work[i__6].r * ab[i__7].i + 
			    work[i__6].i * ab[i__7].r;
#line 470 "zgbbrd.f"
		    work[i__5].r = z__1.r, work[i__5].i = z__1.i;
#line 471 "zgbbrd.f"
		    i__5 = klu1 + (j + kun) * ab_dim1;
#line 471 "zgbbrd.f"
		    i__6 = j + kun;
#line 471 "zgbbrd.f"
		    i__7 = klu1 + (j + kun) * ab_dim1;
#line 471 "zgbbrd.f"
		    z__1.r = rwork[i__6] * ab[i__7].r, z__1.i = rwork[i__6] * 
			    ab[i__7].i;
#line 471 "zgbbrd.f"
		    ab[i__5].r = z__1.r, ab[i__5].i = z__1.i;
#line 472 "zgbbrd.f"
/* L70: */
#line 472 "zgbbrd.f"
		}

#line 474 "zgbbrd.f"
		if (ml > ml0) {
#line 475 "zgbbrd.f"
		    --ml;
#line 476 "zgbbrd.f"
		} else {
#line 477 "zgbbrd.f"
		    --mu;
#line 478 "zgbbrd.f"
		}
#line 479 "zgbbrd.f"
/* L80: */
#line 479 "zgbbrd.f"
	    }
#line 480 "zgbbrd.f"
/* L90: */
#line 480 "zgbbrd.f"
	}
#line 481 "zgbbrd.f"
    }

#line 483 "zgbbrd.f"
    if (*ku == 0 && *kl > 0) {

/*        A has been reduced to complex lower bidiagonal form */

/*        Transform lower bidiagonal form to upper bidiagonal by applying */
/*        plane rotations from the left, overwriting superdiagonal */
/*        elements on subdiagonal elements */

/* Computing MIN */
#line 491 "zgbbrd.f"
	i__2 = *m - 1;
#line 491 "zgbbrd.f"
	i__1 = min(i__2,*n);
#line 491 "zgbbrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 492 "zgbbrd.f"
	    zlartg_(&ab[i__ * ab_dim1 + 1], &ab[i__ * ab_dim1 + 2], &rc, &rs, 
		    &ra);
#line 493 "zgbbrd.f"
	    i__2 = i__ * ab_dim1 + 1;
#line 493 "zgbbrd.f"
	    ab[i__2].r = ra.r, ab[i__2].i = ra.i;
#line 494 "zgbbrd.f"
	    if (i__ < *n) {
#line 495 "zgbbrd.f"
		i__2 = i__ * ab_dim1 + 2;
#line 495 "zgbbrd.f"
		i__4 = (i__ + 1) * ab_dim1 + 1;
#line 495 "zgbbrd.f"
		z__1.r = rs.r * ab[i__4].r - rs.i * ab[i__4].i, z__1.i = rs.r 
			* ab[i__4].i + rs.i * ab[i__4].r;
#line 495 "zgbbrd.f"
		ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 496 "zgbbrd.f"
		i__2 = (i__ + 1) * ab_dim1 + 1;
#line 496 "zgbbrd.f"
		i__4 = (i__ + 1) * ab_dim1 + 1;
#line 496 "zgbbrd.f"
		z__1.r = rc * ab[i__4].r, z__1.i = rc * ab[i__4].i;
#line 496 "zgbbrd.f"
		ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 497 "zgbbrd.f"
	    }
#line 498 "zgbbrd.f"
	    if (wantq) {
#line 498 "zgbbrd.f"
		d_cnjg(&z__1, &rs);
#line 498 "zgbbrd.f"
		zrot_(m, &q[i__ * q_dim1 + 1], &c__1, &q[(i__ + 1) * q_dim1 + 
			1], &c__1, &rc, &z__1);
#line 498 "zgbbrd.f"
	    }
#line 501 "zgbbrd.f"
	    if (wantc) {
#line 501 "zgbbrd.f"
		zrot_(ncc, &c__[i__ + c_dim1], ldc, &c__[i__ + 1 + c_dim1], 
			ldc, &rc, &rs);
#line 501 "zgbbrd.f"
	    }
#line 504 "zgbbrd.f"
/* L100: */
#line 504 "zgbbrd.f"
	}
#line 505 "zgbbrd.f"
    } else {

/*        A has been reduced to complex upper bidiagonal form or is */
/*        diagonal */

#line 510 "zgbbrd.f"
	if (*ku > 0 && *m < *n) {

/*           Annihilate a(m,m+1) by applying plane rotations from the */
/*           right */

#line 515 "zgbbrd.f"
	    i__1 = *ku + (*m + 1) * ab_dim1;
#line 515 "zgbbrd.f"
	    rb.r = ab[i__1].r, rb.i = ab[i__1].i;
#line 516 "zgbbrd.f"
	    for (i__ = *m; i__ >= 1; --i__) {
#line 517 "zgbbrd.f"
		zlartg_(&ab[*ku + 1 + i__ * ab_dim1], &rb, &rc, &rs, &ra);
#line 518 "zgbbrd.f"
		i__1 = *ku + 1 + i__ * ab_dim1;
#line 518 "zgbbrd.f"
		ab[i__1].r = ra.r, ab[i__1].i = ra.i;
#line 519 "zgbbrd.f"
		if (i__ > 1) {
#line 520 "zgbbrd.f"
		    d_cnjg(&z__3, &rs);
#line 520 "zgbbrd.f"
		    z__2.r = -z__3.r, z__2.i = -z__3.i;
#line 520 "zgbbrd.f"
		    i__1 = *ku + i__ * ab_dim1;
#line 520 "zgbbrd.f"
		    z__1.r = z__2.r * ab[i__1].r - z__2.i * ab[i__1].i, 
			    z__1.i = z__2.r * ab[i__1].i + z__2.i * ab[i__1]
			    .r;
#line 520 "zgbbrd.f"
		    rb.r = z__1.r, rb.i = z__1.i;
#line 521 "zgbbrd.f"
		    i__1 = *ku + i__ * ab_dim1;
#line 521 "zgbbrd.f"
		    i__2 = *ku + i__ * ab_dim1;
#line 521 "zgbbrd.f"
		    z__1.r = rc * ab[i__2].r, z__1.i = rc * ab[i__2].i;
#line 521 "zgbbrd.f"
		    ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 522 "zgbbrd.f"
		}
#line 523 "zgbbrd.f"
		if (wantpt) {
#line 523 "zgbbrd.f"
		    d_cnjg(&z__1, &rs);
#line 523 "zgbbrd.f"
		    zrot_(n, &pt[i__ + pt_dim1], ldpt, &pt[*m + 1 + pt_dim1], 
			    ldpt, &rc, &z__1);
#line 523 "zgbbrd.f"
		}
#line 526 "zgbbrd.f"
/* L110: */
#line 526 "zgbbrd.f"
	    }
#line 527 "zgbbrd.f"
	}
#line 528 "zgbbrd.f"
    }

/*     Make diagonal and superdiagonal elements real, storing them in D */
/*     and E */

#line 533 "zgbbrd.f"
    i__1 = *ku + 1 + ab_dim1;
#line 533 "zgbbrd.f"
    t.r = ab[i__1].r, t.i = ab[i__1].i;
#line 534 "zgbbrd.f"
    i__1 = minmn;
#line 534 "zgbbrd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 535 "zgbbrd.f"
	abst = z_abs(&t);
#line 536 "zgbbrd.f"
	d__[i__] = abst;
#line 537 "zgbbrd.f"
	if (abst != 0.) {
#line 538 "zgbbrd.f"
	    z__1.r = t.r / abst, z__1.i = t.i / abst;
#line 538 "zgbbrd.f"
	    t.r = z__1.r, t.i = z__1.i;
#line 539 "zgbbrd.f"
	} else {
#line 540 "zgbbrd.f"
	    t.r = 1., t.i = 0.;
#line 541 "zgbbrd.f"
	}
#line 542 "zgbbrd.f"
	if (wantq) {
#line 542 "zgbbrd.f"
	    zscal_(m, &t, &q[i__ * q_dim1 + 1], &c__1);
#line 542 "zgbbrd.f"
	}
#line 544 "zgbbrd.f"
	if (wantc) {
#line 544 "zgbbrd.f"
	    d_cnjg(&z__1, &t);
#line 544 "zgbbrd.f"
	    zscal_(ncc, &z__1, &c__[i__ + c_dim1], ldc);
#line 544 "zgbbrd.f"
	}
#line 546 "zgbbrd.f"
	if (i__ < minmn) {
#line 547 "zgbbrd.f"
	    if (*ku == 0 && *kl == 0) {
#line 548 "zgbbrd.f"
		e[i__] = 0.;
#line 549 "zgbbrd.f"
		i__2 = (i__ + 1) * ab_dim1 + 1;
#line 549 "zgbbrd.f"
		t.r = ab[i__2].r, t.i = ab[i__2].i;
#line 550 "zgbbrd.f"
	    } else {
#line 551 "zgbbrd.f"
		if (*ku == 0) {
#line 552 "zgbbrd.f"
		    i__2 = i__ * ab_dim1 + 2;
#line 552 "zgbbrd.f"
		    d_cnjg(&z__2, &t);
#line 552 "zgbbrd.f"
		    z__1.r = ab[i__2].r * z__2.r - ab[i__2].i * z__2.i, 
			    z__1.i = ab[i__2].r * z__2.i + ab[i__2].i * 
			    z__2.r;
#line 552 "zgbbrd.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 553 "zgbbrd.f"
		} else {
#line 554 "zgbbrd.f"
		    i__2 = *ku + (i__ + 1) * ab_dim1;
#line 554 "zgbbrd.f"
		    d_cnjg(&z__2, &t);
#line 554 "zgbbrd.f"
		    z__1.r = ab[i__2].r * z__2.r - ab[i__2].i * z__2.i, 
			    z__1.i = ab[i__2].r * z__2.i + ab[i__2].i * 
			    z__2.r;
#line 554 "zgbbrd.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 555 "zgbbrd.f"
		}
#line 556 "zgbbrd.f"
		abst = z_abs(&t);
#line 557 "zgbbrd.f"
		e[i__] = abst;
#line 558 "zgbbrd.f"
		if (abst != 0.) {
#line 559 "zgbbrd.f"
		    z__1.r = t.r / abst, z__1.i = t.i / abst;
#line 559 "zgbbrd.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 560 "zgbbrd.f"
		} else {
#line 561 "zgbbrd.f"
		    t.r = 1., t.i = 0.;
#line 562 "zgbbrd.f"
		}
#line 563 "zgbbrd.f"
		if (wantpt) {
#line 563 "zgbbrd.f"
		    zscal_(n, &t, &pt[i__ + 1 + pt_dim1], ldpt);
#line 563 "zgbbrd.f"
		}
#line 565 "zgbbrd.f"
		i__2 = *ku + 1 + (i__ + 1) * ab_dim1;
#line 565 "zgbbrd.f"
		d_cnjg(&z__2, &t);
#line 565 "zgbbrd.f"
		z__1.r = ab[i__2].r * z__2.r - ab[i__2].i * z__2.i, z__1.i = 
			ab[i__2].r * z__2.i + ab[i__2].i * z__2.r;
#line 565 "zgbbrd.f"
		t.r = z__1.r, t.i = z__1.i;
#line 566 "zgbbrd.f"
	    }
#line 567 "zgbbrd.f"
	}
#line 568 "zgbbrd.f"
/* L120: */
#line 568 "zgbbrd.f"
    }
#line 569 "zgbbrd.f"
    return 0;

/*     End of ZGBBRD */

} /* zgbbrd_ */

