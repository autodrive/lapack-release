#line 1 "dgebal.f"
/* dgebal.f -- translated by f2c (version 20100827).
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

#line 1 "dgebal.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DGEBAL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEBAL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgebal.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgebal.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgebal.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOB */
/*       INTEGER            IHI, ILO, INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), SCALE( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEBAL balances a general real matrix A.  This involves, first, */
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
/* >          A is DOUBLE array, dimension (LDA,N) */
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
/* >          ILO is INTEGER */
/* > \endverbatim */
/* > \param[out] IHI */
/* > \verbatim */
/* >          IHI is INTEGER */
/* >          ILO and IHI are set to integers such that on exit */
/* >          A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N. */
/* >          If JOB = 'N' or 'S', ILO = 1 and IHI = N. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is DOUBLE array, dimension (N) */
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

/* > \date November 2013 */

/* > \ingroup doubleGEcomputational */

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
/* >  This subroutine is based on the EISPACK routine BALANC. */
/* > */
/* >  Modified by Tzu-Yi Chen, Computer Science Division, University of */
/* >    California at Berkeley, USA */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dgebal_(char *job, integer *n, doublereal *a, integer *
	lda, integer *ilo, integer *ihi, doublereal *scale, integer *info, 
	ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal c__, f, g;
    static integer i__, j, k, l, m;
    static doublereal r__, s, ca, ra;
    static integer ica, ira, iexc;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal sfmin1, sfmin2, sfmax1, sfmax2;
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern logical disnan_(doublereal *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical noconv;


/*  -- LAPACK computational routine (version 3.5.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2013 */

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

/*     Test the input parameters */

#line 208 "dgebal.f"
    /* Parameter adjustments */
#line 208 "dgebal.f"
    a_dim1 = *lda;
#line 208 "dgebal.f"
    a_offset = 1 + a_dim1;
#line 208 "dgebal.f"
    a -= a_offset;
#line 208 "dgebal.f"
    --scale;
#line 208 "dgebal.f"

#line 208 "dgebal.f"
    /* Function Body */
#line 208 "dgebal.f"
    *info = 0;
#line 209 "dgebal.f"
    if (! lsame_(job, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(job, "P", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(job, "S", (ftnlen)1, (ftnlen)1) 
	    && ! lsame_(job, "B", (ftnlen)1, (ftnlen)1)) {
#line 211 "dgebal.f"
	*info = -1;
#line 212 "dgebal.f"
    } else if (*n < 0) {
#line 213 "dgebal.f"
	*info = -2;
#line 214 "dgebal.f"
    } else if (*lda < max(1,*n)) {
#line 215 "dgebal.f"
	*info = -4;
#line 216 "dgebal.f"
    }
#line 217 "dgebal.f"
    if (*info != 0) {
#line 218 "dgebal.f"
	i__1 = -(*info);
#line 218 "dgebal.f"
	xerbla_("DGEBAL", &i__1, (ftnlen)6);
#line 219 "dgebal.f"
	return 0;
#line 220 "dgebal.f"
    }

#line 222 "dgebal.f"
    k = 1;
#line 223 "dgebal.f"
    l = *n;

#line 225 "dgebal.f"
    if (*n == 0) {
#line 225 "dgebal.f"
	goto L210;
#line 225 "dgebal.f"
    }

#line 228 "dgebal.f"
    if (lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 229 "dgebal.f"
	i__1 = *n;
#line 229 "dgebal.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 230 "dgebal.f"
	    scale[i__] = 1.;
#line 231 "dgebal.f"
/* L10: */
#line 231 "dgebal.f"
	}
#line 232 "dgebal.f"
	goto L210;
#line 233 "dgebal.f"
    }

#line 235 "dgebal.f"
    if (lsame_(job, "S", (ftnlen)1, (ftnlen)1)) {
#line 235 "dgebal.f"
	goto L120;
#line 235 "dgebal.f"
    }

/*     Permutation to isolate eigenvalues if possible */

#line 240 "dgebal.f"
    goto L50;

/*     Row and column exchange. */

#line 244 "dgebal.f"
L20:
#line 245 "dgebal.f"
    scale[m] = (doublereal) j;
#line 246 "dgebal.f"
    if (j == m) {
#line 246 "dgebal.f"
	goto L30;
#line 246 "dgebal.f"
    }

#line 249 "dgebal.f"
    dswap_(&l, &a[j * a_dim1 + 1], &c__1, &a[m * a_dim1 + 1], &c__1);
#line 250 "dgebal.f"
    i__1 = *n - k + 1;
#line 250 "dgebal.f"
    dswap_(&i__1, &a[j + k * a_dim1], lda, &a[m + k * a_dim1], lda);

#line 252 "dgebal.f"
L30:
#line 253 "dgebal.f"
    switch (iexc) {
#line 253 "dgebal.f"
	case 1:  goto L40;
#line 253 "dgebal.f"
	case 2:  goto L80;
#line 253 "dgebal.f"
    }

/*     Search for rows isolating an eigenvalue and push them down. */

#line 257 "dgebal.f"
L40:
#line 258 "dgebal.f"
    if (l == 1) {
#line 258 "dgebal.f"
	goto L210;
#line 258 "dgebal.f"
    }
#line 260 "dgebal.f"
    --l;

#line 262 "dgebal.f"
L50:
#line 263 "dgebal.f"
    for (j = l; j >= 1; --j) {

#line 265 "dgebal.f"
	i__1 = l;
#line 265 "dgebal.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 266 "dgebal.f"
	    if (i__ == j) {
#line 266 "dgebal.f"
		goto L60;
#line 266 "dgebal.f"
	    }
#line 268 "dgebal.f"
	    if (a[j + i__ * a_dim1] != 0.) {
#line 268 "dgebal.f"
		goto L70;
#line 268 "dgebal.f"
	    }
#line 270 "dgebal.f"
L60:
#line 270 "dgebal.f"
	    ;
#line 270 "dgebal.f"
	}

#line 272 "dgebal.f"
	m = l;
#line 273 "dgebal.f"
	iexc = 1;
#line 274 "dgebal.f"
	goto L20;
#line 275 "dgebal.f"
L70:
#line 275 "dgebal.f"
	;
#line 275 "dgebal.f"
    }

#line 277 "dgebal.f"
    goto L90;

/*     Search for columns isolating an eigenvalue and push them left. */

#line 281 "dgebal.f"
L80:
#line 282 "dgebal.f"
    ++k;

#line 284 "dgebal.f"
L90:
#line 285 "dgebal.f"
    i__1 = l;
#line 285 "dgebal.f"
    for (j = k; j <= i__1; ++j) {

#line 287 "dgebal.f"
	i__2 = l;
#line 287 "dgebal.f"
	for (i__ = k; i__ <= i__2; ++i__) {
#line 288 "dgebal.f"
	    if (i__ == j) {
#line 288 "dgebal.f"
		goto L100;
#line 288 "dgebal.f"
	    }
#line 290 "dgebal.f"
	    if (a[i__ + j * a_dim1] != 0.) {
#line 290 "dgebal.f"
		goto L110;
#line 290 "dgebal.f"
	    }
#line 292 "dgebal.f"
L100:
#line 292 "dgebal.f"
	    ;
#line 292 "dgebal.f"
	}

#line 294 "dgebal.f"
	m = k;
#line 295 "dgebal.f"
	iexc = 2;
#line 296 "dgebal.f"
	goto L20;
#line 297 "dgebal.f"
L110:
#line 297 "dgebal.f"
	;
#line 297 "dgebal.f"
    }

#line 299 "dgebal.f"
L120:
#line 300 "dgebal.f"
    i__1 = l;
#line 300 "dgebal.f"
    for (i__ = k; i__ <= i__1; ++i__) {
#line 301 "dgebal.f"
	scale[i__] = 1.;
#line 302 "dgebal.f"
/* L130: */
#line 302 "dgebal.f"
    }

#line 304 "dgebal.f"
    if (lsame_(job, "P", (ftnlen)1, (ftnlen)1)) {
#line 304 "dgebal.f"
	goto L210;
#line 304 "dgebal.f"
    }

/*     Balance the submatrix in rows K to L. */

/*     Iterative loop for norm reduction */

#line 311 "dgebal.f"
    sfmin1 = dlamch_("S", (ftnlen)1) / dlamch_("P", (ftnlen)1);
#line 312 "dgebal.f"
    sfmax1 = 1. / sfmin1;
#line 313 "dgebal.f"
    sfmin2 = sfmin1 * 2.;
#line 314 "dgebal.f"
    sfmax2 = 1. / sfmin2;

#line 316 "dgebal.f"
L140:
#line 317 "dgebal.f"
    noconv = FALSE_;

#line 319 "dgebal.f"
    i__1 = l;
#line 319 "dgebal.f"
    for (i__ = k; i__ <= i__1; ++i__) {

#line 321 "dgebal.f"
	i__2 = l - k + 1;
#line 321 "dgebal.f"
	c__ = dnrm2_(&i__2, &a[k + i__ * a_dim1], &c__1);
#line 322 "dgebal.f"
	i__2 = l - k + 1;
#line 322 "dgebal.f"
	r__ = dnrm2_(&i__2, &a[i__ + k * a_dim1], lda);
#line 323 "dgebal.f"
	ica = idamax_(&l, &a[i__ * a_dim1 + 1], &c__1);
#line 324 "dgebal.f"
	ca = (d__1 = a[ica + i__ * a_dim1], abs(d__1));
#line 325 "dgebal.f"
	i__2 = *n - k + 1;
#line 325 "dgebal.f"
	ira = idamax_(&i__2, &a[i__ + k * a_dim1], lda);
#line 326 "dgebal.f"
	ra = (d__1 = a[i__ + (ira + k - 1) * a_dim1], abs(d__1));

/*        Guard against zero C or R due to underflow. */

#line 330 "dgebal.f"
	if (c__ == 0. || r__ == 0.) {
#line 330 "dgebal.f"
	    goto L200;
#line 330 "dgebal.f"
	}
#line 332 "dgebal.f"
	g = r__ / 2.;
#line 333 "dgebal.f"
	f = 1.;
#line 334 "dgebal.f"
	s = c__ + r__;
#line 335 "dgebal.f"
L160:
/* Computing MAX */
#line 336 "dgebal.f"
	d__1 = max(f,c__);
/* Computing MIN */
#line 336 "dgebal.f"
	d__2 = min(r__,g);
#line 336 "dgebal.f"
	if (c__ >= g || max(d__1,ca) >= sfmax2 || min(d__2,ra) <= sfmin2) {
#line 336 "dgebal.f"
	    goto L170;
#line 336 "dgebal.f"
	}
#line 338 "dgebal.f"
	d__1 = c__ + f + ca + r__ + g + ra;
#line 338 "dgebal.f"
	if (disnan_(&d__1)) {

/*           Exit if NaN to avoid infinite loop */

#line 342 "dgebal.f"
	    *info = -3;
#line 343 "dgebal.f"
	    i__2 = -(*info);
#line 343 "dgebal.f"
	    xerbla_("DGEBAL", &i__2, (ftnlen)6);
#line 344 "dgebal.f"
	    return 0;
#line 345 "dgebal.f"
	}
#line 346 "dgebal.f"
	f *= 2.;
#line 347 "dgebal.f"
	c__ *= 2.;
#line 348 "dgebal.f"
	ca *= 2.;
#line 349 "dgebal.f"
	r__ /= 2.;
#line 350 "dgebal.f"
	g /= 2.;
#line 351 "dgebal.f"
	ra /= 2.;
#line 352 "dgebal.f"
	goto L160;

#line 354 "dgebal.f"
L170:
#line 355 "dgebal.f"
	g = c__ / 2.;
#line 356 "dgebal.f"
L180:
/* Computing MIN */
#line 357 "dgebal.f"
	d__1 = min(f,c__), d__1 = min(d__1,g);
#line 357 "dgebal.f"
	if (g < r__ || max(r__,ra) >= sfmax2 || min(d__1,ca) <= sfmin2) {
#line 357 "dgebal.f"
	    goto L190;
#line 357 "dgebal.f"
	}
#line 359 "dgebal.f"
	f /= 2.;
#line 360 "dgebal.f"
	c__ /= 2.;
#line 361 "dgebal.f"
	g /= 2.;
#line 362 "dgebal.f"
	ca /= 2.;
#line 363 "dgebal.f"
	r__ *= 2.;
#line 364 "dgebal.f"
	ra *= 2.;
#line 365 "dgebal.f"
	goto L180;

/*        Now balance. */

#line 369 "dgebal.f"
L190:
#line 370 "dgebal.f"
	if (c__ + r__ >= s * .95) {
#line 370 "dgebal.f"
	    goto L200;
#line 370 "dgebal.f"
	}
#line 372 "dgebal.f"
	if (f < 1. && scale[i__] < 1.) {
#line 373 "dgebal.f"
	    if (f * scale[i__] <= sfmin1) {
#line 373 "dgebal.f"
		goto L200;
#line 373 "dgebal.f"
	    }
#line 375 "dgebal.f"
	}
#line 376 "dgebal.f"
	if (f > 1. && scale[i__] > 1.) {
#line 377 "dgebal.f"
	    if (scale[i__] >= sfmax1 / f) {
#line 377 "dgebal.f"
		goto L200;
#line 377 "dgebal.f"
	    }
#line 379 "dgebal.f"
	}
#line 380 "dgebal.f"
	g = 1. / f;
#line 381 "dgebal.f"
	scale[i__] *= f;
#line 382 "dgebal.f"
	noconv = TRUE_;

#line 384 "dgebal.f"
	i__2 = *n - k + 1;
#line 384 "dgebal.f"
	dscal_(&i__2, &g, &a[i__ + k * a_dim1], lda);
#line 385 "dgebal.f"
	dscal_(&l, &f, &a[i__ * a_dim1 + 1], &c__1);

#line 387 "dgebal.f"
L200:
#line 387 "dgebal.f"
	;
#line 387 "dgebal.f"
    }

#line 389 "dgebal.f"
    if (noconv) {
#line 389 "dgebal.f"
	goto L140;
#line 389 "dgebal.f"
    }

#line 392 "dgebal.f"
L210:
#line 393 "dgebal.f"
    *ilo = k;
#line 394 "dgebal.f"
    *ihi = l;

#line 396 "dgebal.f"
    return 0;

/*     End of DGEBAL */

} /* dgebal_ */

