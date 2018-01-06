#line 1 "zgebal.f"
/* zgebal.f -- translated by f2c (version 20100827).
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

#line 1 "zgebal.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZGEBAL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEBAL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgebal.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgebal.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgebal.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOB */
/*       INTEGER            IHI, ILO, INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   SCALE( * ) */
/*       COMPLEX*16         A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEBAL balances a general complex matrix A.  This involves, first, */
/* > permuting A by a similarity transformation to isolate eigenvalues */
/* > in the first 1 to ILO-1 and last IHI+1 to N elements on the */
/* > diagonal; and second, applying a diagonal similarity transformation */
/* > to rows and columns ILO to IHI to make the rows and columns as */
/* > close in norm as possible.  Both steps are optional. */
/* > */
/* > Balancing may reduce the 1-norm of the matrix, and improve the */
/* > accuracy of the computed eigenvalues and/or eigenvectors. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOB */
/* > \verbatim */
/* >          JOB is CHARACTER*1 */
/* >          Specifies the operations to be performed on A: */
/* >          = 'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0 */
/* >                  for i = 1,...,N; */
/* >          = 'P':  permute only; */
/* >          = 'S':  scale only; */
/* >          = 'B':  both permute and scale. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the input matrix A. */
/* >          On exit,  A is overwritten by the balanced matrix. */
/* >          If JOB = 'N', A is not referenced. */
/* >          See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] ILO */
/* > \verbatim */
/* > \endverbatim */
/* > */
/* > \param[out] IHI */
/* > \verbatim */
/* >          ILO and IHI are set to INTEGER such that on exit */
/* >          A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N. */
/* >          If JOB = 'N' or 'S', ILO = 1 and IHI = N. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is DOUBLE PRECISION array, dimension (N) */
/* >          Details of the permutations and scaling factors applied to */
/* >          A.  If P(j) is the index of the row and column interchanged */
/* >          with row and column j and D(j) is the scaling factor */
/* >          applied to row and column j, then */
/* >          SCALE(j) = P(j)    for j = 1,...,ILO-1 */
/* >                   = D(j)    for j = ILO,...,IHI */
/* >                   = P(j)    for j = IHI+1,...,N. */
/* >          The order in which the interchanges are made is N to IHI+1, */
/* >          then 1 to ILO-1. */
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

/* > \date November 2015 */

/* > \ingroup complex16GEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The permutations consist of row and column interchanges which put */
/* >  the matrix in the form */
/* > */
/* >             ( T1   X   Y  ) */
/* >     P A P = (  0   B   Z  ) */
/* >             (  0   0   T2 ) */
/* > */
/* >  where T1 and T2 are upper triangular matrices whose eigenvalues lie */
/* >  along the diagonal.  The column indices ILO and IHI mark the starting */
/* >  and ending columns of the submatrix B. Balancing consists of applying */
/* >  a diagonal similarity transformation inv(D) * B * D to make the */
/* >  1-norms of each row of B and its corresponding column nearly equal. */
/* >  The output matrix is */
/* > */
/* >     ( T1     X*D          Y    ) */
/* >     (  0  inv(D)*B*D  inv(D)*Z ). */
/* >     (  0      0           T2   ) */
/* > */
/* >  Information about the permutations P and the diagonal matrix D is */
/* >  returned in the vector SCALE. */
/* > */
/* >  This subroutine is based on the EISPACK routine CBAL. */
/* > */
/* >  Modified by Tzu-Yi Chen, Computer Science Division, University of */
/* >    California at Berkeley, USA */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgebal_(char *job, integer *n, doublecomplex *a, integer 
	*lda, integer *ilo, integer *ihi, doublereal *scale, integer *info, 
	ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_imag(doublecomplex *), z_abs(doublecomplex *);

    /* Local variables */
    static doublereal c__, f, g;
    static integer i__, j, k, l, m;
    static doublereal r__, s, ca, ra;
    static integer ica, ira, iexc;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static doublereal sfmin1, sfmin2, sfmax1, sfmax2;
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *), dlamch_(
	    char *, ftnlen);
    extern logical disnan_(doublereal *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    extern integer izamax_(integer *, doublecomplex *, integer *);
    static logical noconv;


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

/*     Test the input parameters */

#line 208 "zgebal.f"
    /* Parameter adjustments */
#line 208 "zgebal.f"
    a_dim1 = *lda;
#line 208 "zgebal.f"
    a_offset = 1 + a_dim1;
#line 208 "zgebal.f"
    a -= a_offset;
#line 208 "zgebal.f"
    --scale;
#line 208 "zgebal.f"

#line 208 "zgebal.f"
    /* Function Body */
#line 208 "zgebal.f"
    *info = 0;
#line 209 "zgebal.f"
    if (! lsame_(job, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(job, "P", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(job, "S", (ftnlen)1, (ftnlen)1) 
	    && ! lsame_(job, "B", (ftnlen)1, (ftnlen)1)) {
#line 211 "zgebal.f"
	*info = -1;
#line 212 "zgebal.f"
    } else if (*n < 0) {
#line 213 "zgebal.f"
	*info = -2;
#line 214 "zgebal.f"
    } else if (*lda < max(1,*n)) {
#line 215 "zgebal.f"
	*info = -4;
#line 216 "zgebal.f"
    }
#line 217 "zgebal.f"
    if (*info != 0) {
#line 218 "zgebal.f"
	i__1 = -(*info);
#line 218 "zgebal.f"
	xerbla_("ZGEBAL", &i__1, (ftnlen)6);
#line 219 "zgebal.f"
	return 0;
#line 220 "zgebal.f"
    }

#line 222 "zgebal.f"
    k = 1;
#line 223 "zgebal.f"
    l = *n;

#line 225 "zgebal.f"
    if (*n == 0) {
#line 225 "zgebal.f"
	goto L210;
#line 225 "zgebal.f"
    }

#line 228 "zgebal.f"
    if (lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 229 "zgebal.f"
	i__1 = *n;
#line 229 "zgebal.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 230 "zgebal.f"
	    scale[i__] = 1.;
#line 231 "zgebal.f"
/* L10: */
#line 231 "zgebal.f"
	}
#line 232 "zgebal.f"
	goto L210;
#line 233 "zgebal.f"
    }

#line 235 "zgebal.f"
    if (lsame_(job, "S", (ftnlen)1, (ftnlen)1)) {
#line 235 "zgebal.f"
	goto L120;
#line 235 "zgebal.f"
    }

/*     Permutation to isolate eigenvalues if possible */

#line 240 "zgebal.f"
    goto L50;

/*     Row and column exchange. */

#line 244 "zgebal.f"
L20:
#line 245 "zgebal.f"
    scale[m] = (doublereal) j;
#line 246 "zgebal.f"
    if (j == m) {
#line 246 "zgebal.f"
	goto L30;
#line 246 "zgebal.f"
    }

#line 249 "zgebal.f"
    zswap_(&l, &a[j * a_dim1 + 1], &c__1, &a[m * a_dim1 + 1], &c__1);
#line 250 "zgebal.f"
    i__1 = *n - k + 1;
#line 250 "zgebal.f"
    zswap_(&i__1, &a[j + k * a_dim1], lda, &a[m + k * a_dim1], lda);

#line 252 "zgebal.f"
L30:
#line 253 "zgebal.f"
    switch (iexc) {
#line 253 "zgebal.f"
	case 1:  goto L40;
#line 253 "zgebal.f"
	case 2:  goto L80;
#line 253 "zgebal.f"
    }

/*     Search for rows isolating an eigenvalue and push them down. */

#line 257 "zgebal.f"
L40:
#line 258 "zgebal.f"
    if (l == 1) {
#line 258 "zgebal.f"
	goto L210;
#line 258 "zgebal.f"
    }
#line 260 "zgebal.f"
    --l;

#line 262 "zgebal.f"
L50:
#line 263 "zgebal.f"
    for (j = l; j >= 1; --j) {

#line 265 "zgebal.f"
	i__1 = l;
#line 265 "zgebal.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 266 "zgebal.f"
	    if (i__ == j) {
#line 266 "zgebal.f"
		goto L60;
#line 266 "zgebal.f"
	    }
#line 268 "zgebal.f"
	    i__2 = j + i__ * a_dim1;
#line 268 "zgebal.f"
	    if (a[i__2].r != 0. || d_imag(&a[j + i__ * a_dim1]) != 0.) {
#line 268 "zgebal.f"
		goto L70;
#line 268 "zgebal.f"
	    }
#line 270 "zgebal.f"
L60:
#line 270 "zgebal.f"
	    ;
#line 270 "zgebal.f"
	}

#line 272 "zgebal.f"
	m = l;
#line 273 "zgebal.f"
	iexc = 1;
#line 274 "zgebal.f"
	goto L20;
#line 275 "zgebal.f"
L70:
#line 275 "zgebal.f"
	;
#line 275 "zgebal.f"
    }

#line 277 "zgebal.f"
    goto L90;

/*     Search for columns isolating an eigenvalue and push them left. */

#line 281 "zgebal.f"
L80:
#line 282 "zgebal.f"
    ++k;

#line 284 "zgebal.f"
L90:
#line 285 "zgebal.f"
    i__1 = l;
#line 285 "zgebal.f"
    for (j = k; j <= i__1; ++j) {

#line 287 "zgebal.f"
	i__2 = l;
#line 287 "zgebal.f"
	for (i__ = k; i__ <= i__2; ++i__) {
#line 288 "zgebal.f"
	    if (i__ == j) {
#line 288 "zgebal.f"
		goto L100;
#line 288 "zgebal.f"
	    }
#line 290 "zgebal.f"
	    i__3 = i__ + j * a_dim1;
#line 290 "zgebal.f"
	    if (a[i__3].r != 0. || d_imag(&a[i__ + j * a_dim1]) != 0.) {
#line 290 "zgebal.f"
		goto L110;
#line 290 "zgebal.f"
	    }
#line 292 "zgebal.f"
L100:
#line 292 "zgebal.f"
	    ;
#line 292 "zgebal.f"
	}

#line 294 "zgebal.f"
	m = k;
#line 295 "zgebal.f"
	iexc = 2;
#line 296 "zgebal.f"
	goto L20;
#line 297 "zgebal.f"
L110:
#line 297 "zgebal.f"
	;
#line 297 "zgebal.f"
    }

#line 299 "zgebal.f"
L120:
#line 300 "zgebal.f"
    i__1 = l;
#line 300 "zgebal.f"
    for (i__ = k; i__ <= i__1; ++i__) {
#line 301 "zgebal.f"
	scale[i__] = 1.;
#line 302 "zgebal.f"
/* L130: */
#line 302 "zgebal.f"
    }

#line 304 "zgebal.f"
    if (lsame_(job, "P", (ftnlen)1, (ftnlen)1)) {
#line 304 "zgebal.f"
	goto L210;
#line 304 "zgebal.f"
    }

/*     Balance the submatrix in rows K to L. */

/*     Iterative loop for norm reduction */

#line 311 "zgebal.f"
    sfmin1 = dlamch_("S", (ftnlen)1) / dlamch_("P", (ftnlen)1);
#line 312 "zgebal.f"
    sfmax1 = 1. / sfmin1;
#line 313 "zgebal.f"
    sfmin2 = sfmin1 * 2.;
#line 314 "zgebal.f"
    sfmax2 = 1. / sfmin2;
#line 315 "zgebal.f"
L140:
#line 316 "zgebal.f"
    noconv = FALSE_;

#line 318 "zgebal.f"
    i__1 = l;
#line 318 "zgebal.f"
    for (i__ = k; i__ <= i__1; ++i__) {

#line 320 "zgebal.f"
	i__2 = l - k + 1;
#line 320 "zgebal.f"
	c__ = dznrm2_(&i__2, &a[k + i__ * a_dim1], &c__1);
#line 321 "zgebal.f"
	i__2 = l - k + 1;
#line 321 "zgebal.f"
	r__ = dznrm2_(&i__2, &a[i__ + k * a_dim1], lda);
#line 322 "zgebal.f"
	ica = izamax_(&l, &a[i__ * a_dim1 + 1], &c__1);
#line 323 "zgebal.f"
	ca = z_abs(&a[ica + i__ * a_dim1]);
#line 324 "zgebal.f"
	i__2 = *n - k + 1;
#line 324 "zgebal.f"
	ira = izamax_(&i__2, &a[i__ + k * a_dim1], lda);
#line 325 "zgebal.f"
	ra = z_abs(&a[i__ + (ira + k - 1) * a_dim1]);

/*        Guard against zero C or R due to underflow. */

#line 329 "zgebal.f"
	if (c__ == 0. || r__ == 0.) {
#line 329 "zgebal.f"
	    goto L200;
#line 329 "zgebal.f"
	}
#line 331 "zgebal.f"
	g = r__ / 2.;
#line 332 "zgebal.f"
	f = 1.;
#line 333 "zgebal.f"
	s = c__ + r__;
#line 334 "zgebal.f"
L160:
/* Computing MAX */
#line 335 "zgebal.f"
	d__1 = max(f,c__);
/* Computing MIN */
#line 335 "zgebal.f"
	d__2 = min(r__,g);
#line 335 "zgebal.f"
	if (c__ >= g || max(d__1,ca) >= sfmax2 || min(d__2,ra) <= sfmin2) {
#line 335 "zgebal.f"
	    goto L170;
#line 335 "zgebal.f"
	}
#line 337 "zgebal.f"
	d__1 = c__ + f + ca + r__ + g + ra;
#line 337 "zgebal.f"
	if (disnan_(&d__1)) {

/*           Exit if NaN to avoid infinite loop */

#line 341 "zgebal.f"
	    *info = -3;
#line 342 "zgebal.f"
	    i__2 = -(*info);
#line 342 "zgebal.f"
	    xerbla_("ZGEBAL", &i__2, (ftnlen)6);
#line 343 "zgebal.f"
	    return 0;
#line 344 "zgebal.f"
	}
#line 345 "zgebal.f"
	f *= 2.;
#line 346 "zgebal.f"
	c__ *= 2.;
#line 347 "zgebal.f"
	ca *= 2.;
#line 348 "zgebal.f"
	r__ /= 2.;
#line 349 "zgebal.f"
	g /= 2.;
#line 350 "zgebal.f"
	ra /= 2.;
#line 351 "zgebal.f"
	goto L160;

#line 353 "zgebal.f"
L170:
#line 354 "zgebal.f"
	g = c__ / 2.;
#line 355 "zgebal.f"
L180:
/* Computing MIN */
#line 356 "zgebal.f"
	d__1 = min(f,c__), d__1 = min(d__1,g);
#line 356 "zgebal.f"
	if (g < r__ || max(r__,ra) >= sfmax2 || min(d__1,ca) <= sfmin2) {
#line 356 "zgebal.f"
	    goto L190;
#line 356 "zgebal.f"
	}
#line 358 "zgebal.f"
	f /= 2.;
#line 359 "zgebal.f"
	c__ /= 2.;
#line 360 "zgebal.f"
	g /= 2.;
#line 361 "zgebal.f"
	ca /= 2.;
#line 362 "zgebal.f"
	r__ *= 2.;
#line 363 "zgebal.f"
	ra *= 2.;
#line 364 "zgebal.f"
	goto L180;

/*        Now balance. */

#line 368 "zgebal.f"
L190:
#line 369 "zgebal.f"
	if (c__ + r__ >= s * .95) {
#line 369 "zgebal.f"
	    goto L200;
#line 369 "zgebal.f"
	}
#line 371 "zgebal.f"
	if (f < 1. && scale[i__] < 1.) {
#line 372 "zgebal.f"
	    if (f * scale[i__] <= sfmin1) {
#line 372 "zgebal.f"
		goto L200;
#line 372 "zgebal.f"
	    }
#line 374 "zgebal.f"
	}
#line 375 "zgebal.f"
	if (f > 1. && scale[i__] > 1.) {
#line 376 "zgebal.f"
	    if (scale[i__] >= sfmax1 / f) {
#line 376 "zgebal.f"
		goto L200;
#line 376 "zgebal.f"
	    }
#line 378 "zgebal.f"
	}
#line 379 "zgebal.f"
	g = 1. / f;
#line 380 "zgebal.f"
	scale[i__] *= f;
#line 381 "zgebal.f"
	noconv = TRUE_;

#line 383 "zgebal.f"
	i__2 = *n - k + 1;
#line 383 "zgebal.f"
	zdscal_(&i__2, &g, &a[i__ + k * a_dim1], lda);
#line 384 "zgebal.f"
	zdscal_(&l, &f, &a[i__ * a_dim1 + 1], &c__1);

#line 386 "zgebal.f"
L200:
#line 386 "zgebal.f"
	;
#line 386 "zgebal.f"
    }

#line 388 "zgebal.f"
    if (noconv) {
#line 388 "zgebal.f"
	goto L140;
#line 388 "zgebal.f"
    }

#line 391 "zgebal.f"
L210:
#line 392 "zgebal.f"
    *ilo = k;
#line 393 "zgebal.f"
    *ihi = l;

#line 395 "zgebal.f"
    return 0;

/*     End of ZGEBAL */

} /* zgebal_ */

