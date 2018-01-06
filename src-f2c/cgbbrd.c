#line 1 "cgbbrd.f"
/* cgbbrd.f -- translated by f2c (version 20100827).
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

#line 1 "cgbbrd.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CGBBRD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGBBRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgbbrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgbbrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgbbrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGBBRD( VECT, M, N, NCC, KL, KU, AB, LDAB, D, E, Q, */
/*                          LDQ, PT, LDPT, C, LDC, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          VECT */
/*       INTEGER            INFO, KL, KU, LDAB, LDC, LDPT, LDQ, M, N, NCC */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), E( * ), RWORK( * ) */
/*       COMPLEX            AB( LDAB, * ), C( LDC, * ), PT( LDPT, * ), */
/*      $                   Q( LDQ, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGBBRD reduces a complex general m-by-n band matrix A to real upper */
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
/* >          AB is COMPLEX array, dimension (LDAB,N) */
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
/* >          Q is COMPLEX array, dimension (LDQ,M) */
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
/* >          PT is COMPLEX array, dimension (LDPT,N) */
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
/* >          C is COMPLEX array, dimension (LDC,NCC) */
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
/* >          WORK is COMPLEX array, dimension (max(M,N)) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (max(M,N)) */
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

/* > \ingroup complexGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int cgbbrd_(char *vect, integer *m, integer *n, integer *ncc,
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
    extern /* Subroutine */ int crot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *), 
	    cscal_(integer *, doublecomplex *, doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical wantb, wantc;
    static integer minmn;
    static logical wantq;
    extern /* Subroutine */ int claset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen), clartg_(doublecomplex *, doublecomplex *, doublereal *, 
	    doublecomplex *, doublecomplex *), xerbla_(char *, integer *, 
	    ftnlen), clargv_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, integer *), clartv_(
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublereal *, doublecomplex *, integer *);
    static logical wantpt;


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

#line 242 "cgbbrd.f"
    /* Parameter adjustments */
#line 242 "cgbbrd.f"
    ab_dim1 = *ldab;
#line 242 "cgbbrd.f"
    ab_offset = 1 + ab_dim1;
#line 242 "cgbbrd.f"
    ab -= ab_offset;
#line 242 "cgbbrd.f"
    --d__;
#line 242 "cgbbrd.f"
    --e;
#line 242 "cgbbrd.f"
    q_dim1 = *ldq;
#line 242 "cgbbrd.f"
    q_offset = 1 + q_dim1;
#line 242 "cgbbrd.f"
    q -= q_offset;
#line 242 "cgbbrd.f"
    pt_dim1 = *ldpt;
#line 242 "cgbbrd.f"
    pt_offset = 1 + pt_dim1;
#line 242 "cgbbrd.f"
    pt -= pt_offset;
#line 242 "cgbbrd.f"
    c_dim1 = *ldc;
#line 242 "cgbbrd.f"
    c_offset = 1 + c_dim1;
#line 242 "cgbbrd.f"
    c__ -= c_offset;
#line 242 "cgbbrd.f"
    --work;
#line 242 "cgbbrd.f"
    --rwork;
#line 242 "cgbbrd.f"

#line 242 "cgbbrd.f"
    /* Function Body */
#line 242 "cgbbrd.f"
    wantb = lsame_(vect, "B", (ftnlen)1, (ftnlen)1);
#line 243 "cgbbrd.f"
    wantq = lsame_(vect, "Q", (ftnlen)1, (ftnlen)1) || wantb;
#line 244 "cgbbrd.f"
    wantpt = lsame_(vect, "P", (ftnlen)1, (ftnlen)1) || wantb;
#line 245 "cgbbrd.f"
    wantc = *ncc > 0;
#line 246 "cgbbrd.f"
    klu1 = *kl + *ku + 1;
#line 247 "cgbbrd.f"
    *info = 0;
#line 248 "cgbbrd.f"
    if (! wantq && ! wantpt && ! lsame_(vect, "N", (ftnlen)1, (ftnlen)1)) {
#line 250 "cgbbrd.f"
	*info = -1;
#line 251 "cgbbrd.f"
    } else if (*m < 0) {
#line 252 "cgbbrd.f"
	*info = -2;
#line 253 "cgbbrd.f"
    } else if (*n < 0) {
#line 254 "cgbbrd.f"
	*info = -3;
#line 255 "cgbbrd.f"
    } else if (*ncc < 0) {
#line 256 "cgbbrd.f"
	*info = -4;
#line 257 "cgbbrd.f"
    } else if (*kl < 0) {
#line 258 "cgbbrd.f"
	*info = -5;
#line 259 "cgbbrd.f"
    } else if (*ku < 0) {
#line 260 "cgbbrd.f"
	*info = -6;
#line 261 "cgbbrd.f"
    } else if (*ldab < klu1) {
#line 262 "cgbbrd.f"
	*info = -8;
#line 263 "cgbbrd.f"
    } else if (*ldq < 1 || wantq && *ldq < max(1,*m)) {
#line 264 "cgbbrd.f"
	*info = -12;
#line 265 "cgbbrd.f"
    } else if (*ldpt < 1 || wantpt && *ldpt < max(1,*n)) {
#line 266 "cgbbrd.f"
	*info = -14;
#line 267 "cgbbrd.f"
    } else if (*ldc < 1 || wantc && *ldc < max(1,*m)) {
#line 268 "cgbbrd.f"
	*info = -16;
#line 269 "cgbbrd.f"
    }
#line 270 "cgbbrd.f"
    if (*info != 0) {
#line 271 "cgbbrd.f"
	i__1 = -(*info);
#line 271 "cgbbrd.f"
	xerbla_("CGBBRD", &i__1, (ftnlen)6);
#line 272 "cgbbrd.f"
	return 0;
#line 273 "cgbbrd.f"
    }

/*     Initialize Q and P**H to the unit matrix, if needed */

#line 277 "cgbbrd.f"
    if (wantq) {
#line 277 "cgbbrd.f"
	claset_("Full", m, m, &c_b1, &c_b2, &q[q_offset], ldq, (ftnlen)4);
#line 277 "cgbbrd.f"
    }
#line 279 "cgbbrd.f"
    if (wantpt) {
#line 279 "cgbbrd.f"
	claset_("Full", n, n, &c_b1, &c_b2, &pt[pt_offset], ldpt, (ftnlen)4);
#line 279 "cgbbrd.f"
    }

/*     Quick return if possible. */

#line 284 "cgbbrd.f"
    if (*m == 0 || *n == 0) {
#line 284 "cgbbrd.f"
	return 0;
#line 284 "cgbbrd.f"
    }

#line 287 "cgbbrd.f"
    minmn = min(*m,*n);

#line 289 "cgbbrd.f"
    if (*kl + *ku > 1) {

/*        Reduce to upper bidiagonal form if KU > 0; if KU = 0, reduce */
/*        first to lower bidiagonal form and then transform to upper */
/*        bidiagonal */

#line 295 "cgbbrd.f"
	if (*ku > 0) {
#line 296 "cgbbrd.f"
	    ml0 = 1;
#line 297 "cgbbrd.f"
	    mu0 = 2;
#line 298 "cgbbrd.f"
	} else {
#line 299 "cgbbrd.f"
	    ml0 = 2;
#line 300 "cgbbrd.f"
	    mu0 = 1;
#line 301 "cgbbrd.f"
	}

/*        Wherever possible, plane rotations are generated and applied in */
/*        vector operations of length NR over the index set J1:J2:KLU1. */

/*        The complex sines of the plane rotations are stored in WORK, */
/*        and the real cosines in RWORK. */

/* Computing MIN */
#line 309 "cgbbrd.f"
	i__1 = *m - 1;
#line 309 "cgbbrd.f"
	klm = min(i__1,*kl);
/* Computing MIN */
#line 310 "cgbbrd.f"
	i__1 = *n - 1;
#line 310 "cgbbrd.f"
	kun = min(i__1,*ku);
#line 311 "cgbbrd.f"
	kb = klm + kun;
#line 312 "cgbbrd.f"
	kb1 = kb + 1;
#line 313 "cgbbrd.f"
	inca = kb1 * *ldab;
#line 314 "cgbbrd.f"
	nr = 0;
#line 315 "cgbbrd.f"
	j1 = klm + 2;
#line 316 "cgbbrd.f"
	j2 = 1 - kun;

#line 318 "cgbbrd.f"
	i__1 = minmn;
#line 318 "cgbbrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Reduce i-th column and i-th row of matrix to bidiagonal form */

#line 322 "cgbbrd.f"
	    ml = klm + 1;
#line 323 "cgbbrd.f"
	    mu = kun + 1;
#line 324 "cgbbrd.f"
	    i__2 = kb;
#line 324 "cgbbrd.f"
	    for (kk = 1; kk <= i__2; ++kk) {
#line 325 "cgbbrd.f"
		j1 += kb;
#line 326 "cgbbrd.f"
		j2 += kb;

/*              generate plane rotations to annihilate nonzero elements */
/*              which have been created below the band */

#line 331 "cgbbrd.f"
		if (nr > 0) {
#line 331 "cgbbrd.f"
		    clargv_(&nr, &ab[klu1 + (j1 - klm - 1) * ab_dim1], &inca, 
			    &work[j1], &kb1, &rwork[j1], &kb1);
#line 331 "cgbbrd.f"
		}

/*              apply plane rotations from the left */

#line 337 "cgbbrd.f"
		i__3 = kb;
#line 337 "cgbbrd.f"
		for (l = 1; l <= i__3; ++l) {
#line 338 "cgbbrd.f"
		    if (j2 - klm + l - 1 > *n) {
#line 339 "cgbbrd.f"
			nrt = nr - 1;
#line 340 "cgbbrd.f"
		    } else {
#line 341 "cgbbrd.f"
			nrt = nr;
#line 342 "cgbbrd.f"
		    }
#line 343 "cgbbrd.f"
		    if (nrt > 0) {
#line 343 "cgbbrd.f"
			clartv_(&nrt, &ab[klu1 - l + (j1 - klm + l - 1) * 
				ab_dim1], &inca, &ab[klu1 - l + 1 + (j1 - klm 
				+ l - 1) * ab_dim1], &inca, &rwork[j1], &work[
				j1], &kb1);
#line 343 "cgbbrd.f"
		    }
#line 347 "cgbbrd.f"
/* L10: */
#line 347 "cgbbrd.f"
		}

#line 349 "cgbbrd.f"
		if (ml > ml0) {
#line 350 "cgbbrd.f"
		    if (ml <= *m - i__ + 1) {

/*                    generate plane rotation to annihilate a(i+ml-1,i) */
/*                    within the band, and apply rotation from the left */

#line 355 "cgbbrd.f"
			clartg_(&ab[*ku + ml - 1 + i__ * ab_dim1], &ab[*ku + 
				ml + i__ * ab_dim1], &rwork[i__ + ml - 1], &
				work[i__ + ml - 1], &ra);
#line 357 "cgbbrd.f"
			i__3 = *ku + ml - 1 + i__ * ab_dim1;
#line 357 "cgbbrd.f"
			ab[i__3].r = ra.r, ab[i__3].i = ra.i;
#line 358 "cgbbrd.f"
			if (i__ < *n) {
/* Computing MIN */
#line 358 "cgbbrd.f"
			    i__4 = *ku + ml - 2, i__5 = *n - i__;
#line 358 "cgbbrd.f"
			    i__3 = min(i__4,i__5);
#line 358 "cgbbrd.f"
			    i__6 = *ldab - 1;
#line 358 "cgbbrd.f"
			    i__7 = *ldab - 1;
#line 358 "cgbbrd.f"
			    crot_(&i__3, &ab[*ku + ml - 2 + (i__ + 1) * 
				    ab_dim1], &i__6, &ab[*ku + ml - 1 + (i__ 
				    + 1) * ab_dim1], &i__7, &rwork[i__ + ml - 
				    1], &work[i__ + ml - 1]);
#line 358 "cgbbrd.f"
			}
#line 363 "cgbbrd.f"
		    }
#line 364 "cgbbrd.f"
		    ++nr;
#line 365 "cgbbrd.f"
		    j1 -= kb1;
#line 366 "cgbbrd.f"
		}

#line 368 "cgbbrd.f"
		if (wantq) {

/*                 accumulate product of plane rotations in Q */

#line 372 "cgbbrd.f"
		    i__3 = j2;
#line 372 "cgbbrd.f"
		    i__4 = kb1;
#line 372 "cgbbrd.f"
		    for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) 
			    {
#line 373 "cgbbrd.f"
			d_cnjg(&z__1, &work[j]);
#line 373 "cgbbrd.f"
			crot_(m, &q[(j - 1) * q_dim1 + 1], &c__1, &q[j * 
				q_dim1 + 1], &c__1, &rwork[j], &z__1);
#line 375 "cgbbrd.f"
/* L20: */
#line 375 "cgbbrd.f"
		    }
#line 376 "cgbbrd.f"
		}

#line 378 "cgbbrd.f"
		if (wantc) {

/*                 apply plane rotations to C */

#line 382 "cgbbrd.f"
		    i__4 = j2;
#line 382 "cgbbrd.f"
		    i__3 = kb1;
#line 382 "cgbbrd.f"
		    for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) 
			    {
#line 383 "cgbbrd.f"
			crot_(ncc, &c__[j - 1 + c_dim1], ldc, &c__[j + c_dim1]
				, ldc, &rwork[j], &work[j]);
#line 385 "cgbbrd.f"
/* L30: */
#line 385 "cgbbrd.f"
		    }
#line 386 "cgbbrd.f"
		}

#line 388 "cgbbrd.f"
		if (j2 + kun > *n) {

/*                 adjust J2 to keep within the bounds of the matrix */

#line 392 "cgbbrd.f"
		    --nr;
#line 393 "cgbbrd.f"
		    j2 -= kb1;
#line 394 "cgbbrd.f"
		}

#line 396 "cgbbrd.f"
		i__3 = j2;
#line 396 "cgbbrd.f"
		i__4 = kb1;
#line 396 "cgbbrd.f"
		for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {

/*                 create nonzero element a(j-1,j+ku) above the band */
/*                 and store it in WORK(n+1:2*n) */

#line 401 "cgbbrd.f"
		    i__5 = j + kun;
#line 401 "cgbbrd.f"
		    i__6 = j;
#line 401 "cgbbrd.f"
		    i__7 = (j + kun) * ab_dim1 + 1;
#line 401 "cgbbrd.f"
		    z__1.r = work[i__6].r * ab[i__7].r - work[i__6].i * ab[
			    i__7].i, z__1.i = work[i__6].r * ab[i__7].i + 
			    work[i__6].i * ab[i__7].r;
#line 401 "cgbbrd.f"
		    work[i__5].r = z__1.r, work[i__5].i = z__1.i;
#line 402 "cgbbrd.f"
		    i__5 = (j + kun) * ab_dim1 + 1;
#line 402 "cgbbrd.f"
		    i__6 = j;
#line 402 "cgbbrd.f"
		    i__7 = (j + kun) * ab_dim1 + 1;
#line 402 "cgbbrd.f"
		    z__1.r = rwork[i__6] * ab[i__7].r, z__1.i = rwork[i__6] * 
			    ab[i__7].i;
#line 402 "cgbbrd.f"
		    ab[i__5].r = z__1.r, ab[i__5].i = z__1.i;
#line 403 "cgbbrd.f"
/* L40: */
#line 403 "cgbbrd.f"
		}

/*              generate plane rotations to annihilate nonzero elements */
/*              which have been generated above the band */

#line 408 "cgbbrd.f"
		if (nr > 0) {
#line 408 "cgbbrd.f"
		    clargv_(&nr, &ab[(j1 + kun - 1) * ab_dim1 + 1], &inca, &
			    work[j1 + kun], &kb1, &rwork[j1 + kun], &kb1);
#line 408 "cgbbrd.f"
		}

/*              apply plane rotations from the right */

#line 415 "cgbbrd.f"
		i__4 = kb;
#line 415 "cgbbrd.f"
		for (l = 1; l <= i__4; ++l) {
#line 416 "cgbbrd.f"
		    if (j2 + l - 1 > *m) {
#line 417 "cgbbrd.f"
			nrt = nr - 1;
#line 418 "cgbbrd.f"
		    } else {
#line 419 "cgbbrd.f"
			nrt = nr;
#line 420 "cgbbrd.f"
		    }
#line 421 "cgbbrd.f"
		    if (nrt > 0) {
#line 421 "cgbbrd.f"
			clartv_(&nrt, &ab[l + 1 + (j1 + kun - 1) * ab_dim1], &
				inca, &ab[l + (j1 + kun) * ab_dim1], &inca, &
				rwork[j1 + kun], &work[j1 + kun], &kb1);
#line 421 "cgbbrd.f"
		    }
#line 425 "cgbbrd.f"
/* L50: */
#line 425 "cgbbrd.f"
		}

#line 427 "cgbbrd.f"
		if (ml == ml0 && mu > mu0) {
#line 428 "cgbbrd.f"
		    if (mu <= *n - i__ + 1) {

/*                    generate plane rotation to annihilate a(i,i+mu-1) */
/*                    within the band, and apply rotation from the right */

#line 433 "cgbbrd.f"
			clartg_(&ab[*ku - mu + 3 + (i__ + mu - 2) * ab_dim1], 
				&ab[*ku - mu + 2 + (i__ + mu - 1) * ab_dim1], 
				&rwork[i__ + mu - 1], &work[i__ + mu - 1], &
				ra);
#line 436 "cgbbrd.f"
			i__4 = *ku - mu + 3 + (i__ + mu - 2) * ab_dim1;
#line 436 "cgbbrd.f"
			ab[i__4].r = ra.r, ab[i__4].i = ra.i;
/* Computing MIN */
#line 437 "cgbbrd.f"
			i__3 = *kl + mu - 2, i__5 = *m - i__;
#line 437 "cgbbrd.f"
			i__4 = min(i__3,i__5);
#line 437 "cgbbrd.f"
			crot_(&i__4, &ab[*ku - mu + 4 + (i__ + mu - 2) * 
				ab_dim1], &c__1, &ab[*ku - mu + 3 + (i__ + mu 
				- 1) * ab_dim1], &c__1, &rwork[i__ + mu - 1], 
				&work[i__ + mu - 1]);
#line 441 "cgbbrd.f"
		    }
#line 442 "cgbbrd.f"
		    ++nr;
#line 443 "cgbbrd.f"
		    j1 -= kb1;
#line 444 "cgbbrd.f"
		}

#line 446 "cgbbrd.f"
		if (wantpt) {

/*                 accumulate product of plane rotations in P**H */

#line 450 "cgbbrd.f"
		    i__4 = j2;
#line 450 "cgbbrd.f"
		    i__3 = kb1;
#line 450 "cgbbrd.f"
		    for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) 
			    {
#line 451 "cgbbrd.f"
			d_cnjg(&z__1, &work[j + kun]);
#line 451 "cgbbrd.f"
			crot_(n, &pt[j + kun - 1 + pt_dim1], ldpt, &pt[j + 
				kun + pt_dim1], ldpt, &rwork[j + kun], &z__1);
#line 454 "cgbbrd.f"
/* L60: */
#line 454 "cgbbrd.f"
		    }
#line 455 "cgbbrd.f"
		}

#line 457 "cgbbrd.f"
		if (j2 + kb > *m) {

/*                 adjust J2 to keep within the bounds of the matrix */

#line 461 "cgbbrd.f"
		    --nr;
#line 462 "cgbbrd.f"
		    j2 -= kb1;
#line 463 "cgbbrd.f"
		}

#line 465 "cgbbrd.f"
		i__3 = j2;
#line 465 "cgbbrd.f"
		i__4 = kb1;
#line 465 "cgbbrd.f"
		for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {

/*                 create nonzero element a(j+kl+ku,j+ku-1) below the */
/*                 band and store it in WORK(1:n) */

#line 470 "cgbbrd.f"
		    i__5 = j + kb;
#line 470 "cgbbrd.f"
		    i__6 = j + kun;
#line 470 "cgbbrd.f"
		    i__7 = klu1 + (j + kun) * ab_dim1;
#line 470 "cgbbrd.f"
		    z__1.r = work[i__6].r * ab[i__7].r - work[i__6].i * ab[
			    i__7].i, z__1.i = work[i__6].r * ab[i__7].i + 
			    work[i__6].i * ab[i__7].r;
#line 470 "cgbbrd.f"
		    work[i__5].r = z__1.r, work[i__5].i = z__1.i;
#line 471 "cgbbrd.f"
		    i__5 = klu1 + (j + kun) * ab_dim1;
#line 471 "cgbbrd.f"
		    i__6 = j + kun;
#line 471 "cgbbrd.f"
		    i__7 = klu1 + (j + kun) * ab_dim1;
#line 471 "cgbbrd.f"
		    z__1.r = rwork[i__6] * ab[i__7].r, z__1.i = rwork[i__6] * 
			    ab[i__7].i;
#line 471 "cgbbrd.f"
		    ab[i__5].r = z__1.r, ab[i__5].i = z__1.i;
#line 472 "cgbbrd.f"
/* L70: */
#line 472 "cgbbrd.f"
		}

#line 474 "cgbbrd.f"
		if (ml > ml0) {
#line 475 "cgbbrd.f"
		    --ml;
#line 476 "cgbbrd.f"
		} else {
#line 477 "cgbbrd.f"
		    --mu;
#line 478 "cgbbrd.f"
		}
#line 479 "cgbbrd.f"
/* L80: */
#line 479 "cgbbrd.f"
	    }
#line 480 "cgbbrd.f"
/* L90: */
#line 480 "cgbbrd.f"
	}
#line 481 "cgbbrd.f"
    }

#line 483 "cgbbrd.f"
    if (*ku == 0 && *kl > 0) {

/*        A has been reduced to complex lower bidiagonal form */

/*        Transform lower bidiagonal form to upper bidiagonal by applying */
/*        plane rotations from the left, overwriting superdiagonal */
/*        elements on subdiagonal elements */

/* Computing MIN */
#line 491 "cgbbrd.f"
	i__2 = *m - 1;
#line 491 "cgbbrd.f"
	i__1 = min(i__2,*n);
#line 491 "cgbbrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 492 "cgbbrd.f"
	    clartg_(&ab[i__ * ab_dim1 + 1], &ab[i__ * ab_dim1 + 2], &rc, &rs, 
		    &ra);
#line 493 "cgbbrd.f"
	    i__2 = i__ * ab_dim1 + 1;
#line 493 "cgbbrd.f"
	    ab[i__2].r = ra.r, ab[i__2].i = ra.i;
#line 494 "cgbbrd.f"
	    if (i__ < *n) {
#line 495 "cgbbrd.f"
		i__2 = i__ * ab_dim1 + 2;
#line 495 "cgbbrd.f"
		i__4 = (i__ + 1) * ab_dim1 + 1;
#line 495 "cgbbrd.f"
		z__1.r = rs.r * ab[i__4].r - rs.i * ab[i__4].i, z__1.i = rs.r 
			* ab[i__4].i + rs.i * ab[i__4].r;
#line 495 "cgbbrd.f"
		ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 496 "cgbbrd.f"
		i__2 = (i__ + 1) * ab_dim1 + 1;
#line 496 "cgbbrd.f"
		i__4 = (i__ + 1) * ab_dim1 + 1;
#line 496 "cgbbrd.f"
		z__1.r = rc * ab[i__4].r, z__1.i = rc * ab[i__4].i;
#line 496 "cgbbrd.f"
		ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 497 "cgbbrd.f"
	    }
#line 498 "cgbbrd.f"
	    if (wantq) {
#line 498 "cgbbrd.f"
		d_cnjg(&z__1, &rs);
#line 498 "cgbbrd.f"
		crot_(m, &q[i__ * q_dim1 + 1], &c__1, &q[(i__ + 1) * q_dim1 + 
			1], &c__1, &rc, &z__1);
#line 498 "cgbbrd.f"
	    }
#line 501 "cgbbrd.f"
	    if (wantc) {
#line 501 "cgbbrd.f"
		crot_(ncc, &c__[i__ + c_dim1], ldc, &c__[i__ + 1 + c_dim1], 
			ldc, &rc, &rs);
#line 501 "cgbbrd.f"
	    }
#line 504 "cgbbrd.f"
/* L100: */
#line 504 "cgbbrd.f"
	}
#line 505 "cgbbrd.f"
    } else {

/*        A has been reduced to complex upper bidiagonal form or is */
/*        diagonal */

#line 510 "cgbbrd.f"
	if (*ku > 0 && *m < *n) {

/*           Annihilate a(m,m+1) by applying plane rotations from the */
/*           right */

#line 515 "cgbbrd.f"
	    i__1 = *ku + (*m + 1) * ab_dim1;
#line 515 "cgbbrd.f"
	    rb.r = ab[i__1].r, rb.i = ab[i__1].i;
#line 516 "cgbbrd.f"
	    for (i__ = *m; i__ >= 1; --i__) {
#line 517 "cgbbrd.f"
		clartg_(&ab[*ku + 1 + i__ * ab_dim1], &rb, &rc, &rs, &ra);
#line 518 "cgbbrd.f"
		i__1 = *ku + 1 + i__ * ab_dim1;
#line 518 "cgbbrd.f"
		ab[i__1].r = ra.r, ab[i__1].i = ra.i;
#line 519 "cgbbrd.f"
		if (i__ > 1) {
#line 520 "cgbbrd.f"
		    d_cnjg(&z__3, &rs);
#line 520 "cgbbrd.f"
		    z__2.r = -z__3.r, z__2.i = -z__3.i;
#line 520 "cgbbrd.f"
		    i__1 = *ku + i__ * ab_dim1;
#line 520 "cgbbrd.f"
		    z__1.r = z__2.r * ab[i__1].r - z__2.i * ab[i__1].i, 
			    z__1.i = z__2.r * ab[i__1].i + z__2.i * ab[i__1]
			    .r;
#line 520 "cgbbrd.f"
		    rb.r = z__1.r, rb.i = z__1.i;
#line 521 "cgbbrd.f"
		    i__1 = *ku + i__ * ab_dim1;
#line 521 "cgbbrd.f"
		    i__2 = *ku + i__ * ab_dim1;
#line 521 "cgbbrd.f"
		    z__1.r = rc * ab[i__2].r, z__1.i = rc * ab[i__2].i;
#line 521 "cgbbrd.f"
		    ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 522 "cgbbrd.f"
		}
#line 523 "cgbbrd.f"
		if (wantpt) {
#line 523 "cgbbrd.f"
		    d_cnjg(&z__1, &rs);
#line 523 "cgbbrd.f"
		    crot_(n, &pt[i__ + pt_dim1], ldpt, &pt[*m + 1 + pt_dim1], 
			    ldpt, &rc, &z__1);
#line 523 "cgbbrd.f"
		}
#line 526 "cgbbrd.f"
/* L110: */
#line 526 "cgbbrd.f"
	    }
#line 527 "cgbbrd.f"
	}
#line 528 "cgbbrd.f"
    }

/*     Make diagonal and superdiagonal elements real, storing them in D */
/*     and E */

#line 533 "cgbbrd.f"
    i__1 = *ku + 1 + ab_dim1;
#line 533 "cgbbrd.f"
    t.r = ab[i__1].r, t.i = ab[i__1].i;
#line 534 "cgbbrd.f"
    i__1 = minmn;
#line 534 "cgbbrd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 535 "cgbbrd.f"
	abst = z_abs(&t);
#line 536 "cgbbrd.f"
	d__[i__] = abst;
#line 537 "cgbbrd.f"
	if (abst != 0.) {
#line 538 "cgbbrd.f"
	    z__1.r = t.r / abst, z__1.i = t.i / abst;
#line 538 "cgbbrd.f"
	    t.r = z__1.r, t.i = z__1.i;
#line 539 "cgbbrd.f"
	} else {
#line 540 "cgbbrd.f"
	    t.r = 1., t.i = 0.;
#line 541 "cgbbrd.f"
	}
#line 542 "cgbbrd.f"
	if (wantq) {
#line 542 "cgbbrd.f"
	    cscal_(m, &t, &q[i__ * q_dim1 + 1], &c__1);
#line 542 "cgbbrd.f"
	}
#line 544 "cgbbrd.f"
	if (wantc) {
#line 544 "cgbbrd.f"
	    d_cnjg(&z__1, &t);
#line 544 "cgbbrd.f"
	    cscal_(ncc, &z__1, &c__[i__ + c_dim1], ldc);
#line 544 "cgbbrd.f"
	}
#line 546 "cgbbrd.f"
	if (i__ < minmn) {
#line 547 "cgbbrd.f"
	    if (*ku == 0 && *kl == 0) {
#line 548 "cgbbrd.f"
		e[i__] = 0.;
#line 549 "cgbbrd.f"
		i__2 = (i__ + 1) * ab_dim1 + 1;
#line 549 "cgbbrd.f"
		t.r = ab[i__2].r, t.i = ab[i__2].i;
#line 550 "cgbbrd.f"
	    } else {
#line 551 "cgbbrd.f"
		if (*ku == 0) {
#line 552 "cgbbrd.f"
		    i__2 = i__ * ab_dim1 + 2;
#line 552 "cgbbrd.f"
		    d_cnjg(&z__2, &t);
#line 552 "cgbbrd.f"
		    z__1.r = ab[i__2].r * z__2.r - ab[i__2].i * z__2.i, 
			    z__1.i = ab[i__2].r * z__2.i + ab[i__2].i * 
			    z__2.r;
#line 552 "cgbbrd.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 553 "cgbbrd.f"
		} else {
#line 554 "cgbbrd.f"
		    i__2 = *ku + (i__ + 1) * ab_dim1;
#line 554 "cgbbrd.f"
		    d_cnjg(&z__2, &t);
#line 554 "cgbbrd.f"
		    z__1.r = ab[i__2].r * z__2.r - ab[i__2].i * z__2.i, 
			    z__1.i = ab[i__2].r * z__2.i + ab[i__2].i * 
			    z__2.r;
#line 554 "cgbbrd.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 555 "cgbbrd.f"
		}
#line 556 "cgbbrd.f"
		abst = z_abs(&t);
#line 557 "cgbbrd.f"
		e[i__] = abst;
#line 558 "cgbbrd.f"
		if (abst != 0.) {
#line 559 "cgbbrd.f"
		    z__1.r = t.r / abst, z__1.i = t.i / abst;
#line 559 "cgbbrd.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 560 "cgbbrd.f"
		} else {
#line 561 "cgbbrd.f"
		    t.r = 1., t.i = 0.;
#line 562 "cgbbrd.f"
		}
#line 563 "cgbbrd.f"
		if (wantpt) {
#line 563 "cgbbrd.f"
		    cscal_(n, &t, &pt[i__ + 1 + pt_dim1], ldpt);
#line 563 "cgbbrd.f"
		}
#line 565 "cgbbrd.f"
		i__2 = *ku + 1 + (i__ + 1) * ab_dim1;
#line 565 "cgbbrd.f"
		d_cnjg(&z__2, &t);
#line 565 "cgbbrd.f"
		z__1.r = ab[i__2].r * z__2.r - ab[i__2].i * z__2.i, z__1.i = 
			ab[i__2].r * z__2.i + ab[i__2].i * z__2.r;
#line 565 "cgbbrd.f"
		t.r = z__1.r, t.i = z__1.i;
#line 566 "cgbbrd.f"
	    }
#line 567 "cgbbrd.f"
	}
#line 568 "cgbbrd.f"
/* L120: */
#line 568 "cgbbrd.f"
    }
#line 569 "cgbbrd.f"
    return 0;

/*     End of CGBBRD */

} /* cgbbrd_ */

