#line 1 "sgebal.f"
/* sgebal.f -- translated by f2c (version 20100827).
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

#line 1 "sgebal.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SGEBAL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGEBAL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgebal.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgebal.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgebal.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOB */
/*       INTEGER            IHI, ILO, INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), SCALE( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEBAL balances a general real matrix A.  This involves, first, */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          SCALE is REAL array, dimension (N) */
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

/* > \ingroup realGEcomputational */

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
/* Subroutine */ int sgebal_(char *job, integer *n, doublereal *a, integer *
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
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), sswap_(integer *, doublereal *, integer *, doublereal 
	    *, integer *);
    static doublereal sfmin1, sfmin2, sfmax1, sfmax2;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);
    extern logical sisnan_(doublereal *);
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

#line 206 "sgebal.f"
    /* Parameter adjustments */
#line 206 "sgebal.f"
    a_dim1 = *lda;
#line 206 "sgebal.f"
    a_offset = 1 + a_dim1;
#line 206 "sgebal.f"
    a -= a_offset;
#line 206 "sgebal.f"
    --scale;
#line 206 "sgebal.f"

#line 206 "sgebal.f"
    /* Function Body */
#line 206 "sgebal.f"
    *info = 0;
#line 207 "sgebal.f"
    if (! lsame_(job, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(job, "P", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(job, "S", (ftnlen)1, (ftnlen)1) 
	    && ! lsame_(job, "B", (ftnlen)1, (ftnlen)1)) {
#line 209 "sgebal.f"
	*info = -1;
#line 210 "sgebal.f"
    } else if (*n < 0) {
#line 211 "sgebal.f"
	*info = -2;
#line 212 "sgebal.f"
    } else if (*lda < max(1,*n)) {
#line 213 "sgebal.f"
	*info = -4;
#line 214 "sgebal.f"
    }
#line 215 "sgebal.f"
    if (*info != 0) {
#line 216 "sgebal.f"
	i__1 = -(*info);
#line 216 "sgebal.f"
	xerbla_("SGEBAL", &i__1, (ftnlen)6);
#line 217 "sgebal.f"
	return 0;
#line 218 "sgebal.f"
    }

#line 220 "sgebal.f"
    k = 1;
#line 221 "sgebal.f"
    l = *n;

#line 223 "sgebal.f"
    if (*n == 0) {
#line 223 "sgebal.f"
	goto L210;
#line 223 "sgebal.f"
    }

#line 226 "sgebal.f"
    if (lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 227 "sgebal.f"
	i__1 = *n;
#line 227 "sgebal.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 228 "sgebal.f"
	    scale[i__] = 1.;
#line 229 "sgebal.f"
/* L10: */
#line 229 "sgebal.f"
	}
#line 230 "sgebal.f"
	goto L210;
#line 231 "sgebal.f"
    }

#line 233 "sgebal.f"
    if (lsame_(job, "S", (ftnlen)1, (ftnlen)1)) {
#line 233 "sgebal.f"
	goto L120;
#line 233 "sgebal.f"
    }

/*     Permutation to isolate eigenvalues if possible */

#line 238 "sgebal.f"
    goto L50;

/*     Row and column exchange. */

#line 242 "sgebal.f"
L20:
#line 243 "sgebal.f"
    scale[m] = (doublereal) j;
#line 244 "sgebal.f"
    if (j == m) {
#line 244 "sgebal.f"
	goto L30;
#line 244 "sgebal.f"
    }

#line 247 "sgebal.f"
    sswap_(&l, &a[j * a_dim1 + 1], &c__1, &a[m * a_dim1 + 1], &c__1);
#line 248 "sgebal.f"
    i__1 = *n - k + 1;
#line 248 "sgebal.f"
    sswap_(&i__1, &a[j + k * a_dim1], lda, &a[m + k * a_dim1], lda);

#line 250 "sgebal.f"
L30:
#line 251 "sgebal.f"
    switch (iexc) {
#line 251 "sgebal.f"
	case 1:  goto L40;
#line 251 "sgebal.f"
	case 2:  goto L80;
#line 251 "sgebal.f"
    }

/*     Search for rows isolating an eigenvalue and push them down. */

#line 255 "sgebal.f"
L40:
#line 256 "sgebal.f"
    if (l == 1) {
#line 256 "sgebal.f"
	goto L210;
#line 256 "sgebal.f"
    }
#line 258 "sgebal.f"
    --l;

#line 260 "sgebal.f"
L50:
#line 261 "sgebal.f"
    for (j = l; j >= 1; --j) {

#line 263 "sgebal.f"
	i__1 = l;
#line 263 "sgebal.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 264 "sgebal.f"
	    if (i__ == j) {
#line 264 "sgebal.f"
		goto L60;
#line 264 "sgebal.f"
	    }
#line 266 "sgebal.f"
	    if (a[j + i__ * a_dim1] != 0.) {
#line 266 "sgebal.f"
		goto L70;
#line 266 "sgebal.f"
	    }
#line 268 "sgebal.f"
L60:
#line 268 "sgebal.f"
	    ;
#line 268 "sgebal.f"
	}

#line 270 "sgebal.f"
	m = l;
#line 271 "sgebal.f"
	iexc = 1;
#line 272 "sgebal.f"
	goto L20;
#line 273 "sgebal.f"
L70:
#line 273 "sgebal.f"
	;
#line 273 "sgebal.f"
    }

#line 275 "sgebal.f"
    goto L90;

/*     Search for columns isolating an eigenvalue and push them left. */

#line 279 "sgebal.f"
L80:
#line 280 "sgebal.f"
    ++k;

#line 282 "sgebal.f"
L90:
#line 283 "sgebal.f"
    i__1 = l;
#line 283 "sgebal.f"
    for (j = k; j <= i__1; ++j) {

#line 285 "sgebal.f"
	i__2 = l;
#line 285 "sgebal.f"
	for (i__ = k; i__ <= i__2; ++i__) {
#line 286 "sgebal.f"
	    if (i__ == j) {
#line 286 "sgebal.f"
		goto L100;
#line 286 "sgebal.f"
	    }
#line 288 "sgebal.f"
	    if (a[i__ + j * a_dim1] != 0.) {
#line 288 "sgebal.f"
		goto L110;
#line 288 "sgebal.f"
	    }
#line 290 "sgebal.f"
L100:
#line 290 "sgebal.f"
	    ;
#line 290 "sgebal.f"
	}

#line 292 "sgebal.f"
	m = k;
#line 293 "sgebal.f"
	iexc = 2;
#line 294 "sgebal.f"
	goto L20;
#line 295 "sgebal.f"
L110:
#line 295 "sgebal.f"
	;
#line 295 "sgebal.f"
    }

#line 297 "sgebal.f"
L120:
#line 298 "sgebal.f"
    i__1 = l;
#line 298 "sgebal.f"
    for (i__ = k; i__ <= i__1; ++i__) {
#line 299 "sgebal.f"
	scale[i__] = 1.;
#line 300 "sgebal.f"
/* L130: */
#line 300 "sgebal.f"
    }

#line 302 "sgebal.f"
    if (lsame_(job, "P", (ftnlen)1, (ftnlen)1)) {
#line 302 "sgebal.f"
	goto L210;
#line 302 "sgebal.f"
    }

/*     Balance the submatrix in rows K to L. */

/*     Iterative loop for norm reduction */

#line 309 "sgebal.f"
    sfmin1 = slamch_("S", (ftnlen)1) / slamch_("P", (ftnlen)1);
#line 310 "sgebal.f"
    sfmax1 = 1. / sfmin1;
#line 311 "sgebal.f"
    sfmin2 = sfmin1 * 2.;
#line 312 "sgebal.f"
    sfmax2 = 1. / sfmin2;
#line 313 "sgebal.f"
L140:
#line 314 "sgebal.f"
    noconv = FALSE_;

#line 316 "sgebal.f"
    i__1 = l;
#line 316 "sgebal.f"
    for (i__ = k; i__ <= i__1; ++i__) {

#line 318 "sgebal.f"
	i__2 = l - k + 1;
#line 318 "sgebal.f"
	c__ = snrm2_(&i__2, &a[k + i__ * a_dim1], &c__1);
#line 319 "sgebal.f"
	i__2 = l - k + 1;
#line 319 "sgebal.f"
	r__ = snrm2_(&i__2, &a[i__ + k * a_dim1], lda);
#line 320 "sgebal.f"
	ica = isamax_(&l, &a[i__ * a_dim1 + 1], &c__1);
#line 321 "sgebal.f"
	ca = (d__1 = a[ica + i__ * a_dim1], abs(d__1));
#line 322 "sgebal.f"
	i__2 = *n - k + 1;
#line 322 "sgebal.f"
	ira = isamax_(&i__2, &a[i__ + k * a_dim1], lda);
#line 323 "sgebal.f"
	ra = (d__1 = a[i__ + (ira + k - 1) * a_dim1], abs(d__1));

/*        Guard against zero C or R due to underflow. */

#line 327 "sgebal.f"
	if (c__ == 0. || r__ == 0.) {
#line 327 "sgebal.f"
	    goto L200;
#line 327 "sgebal.f"
	}
#line 329 "sgebal.f"
	g = r__ / 2.;
#line 330 "sgebal.f"
	f = 1.;
#line 331 "sgebal.f"
	s = c__ + r__;
#line 332 "sgebal.f"
L160:
/* Computing MAX */
#line 333 "sgebal.f"
	d__1 = max(f,c__);
/* Computing MIN */
#line 333 "sgebal.f"
	d__2 = min(r__,g);
#line 333 "sgebal.f"
	if (c__ >= g || max(d__1,ca) >= sfmax2 || min(d__2,ra) <= sfmin2) {
#line 333 "sgebal.f"
	    goto L170;
#line 333 "sgebal.f"
	}
#line 335 "sgebal.f"
	f *= 2.;
#line 336 "sgebal.f"
	c__ *= 2.;
#line 337 "sgebal.f"
	ca *= 2.;
#line 338 "sgebal.f"
	r__ /= 2.;
#line 339 "sgebal.f"
	g /= 2.;
#line 340 "sgebal.f"
	ra /= 2.;
#line 341 "sgebal.f"
	goto L160;

#line 343 "sgebal.f"
L170:
#line 344 "sgebal.f"
	g = c__ / 2.;
#line 345 "sgebal.f"
L180:
/* Computing MIN */
#line 346 "sgebal.f"
	d__1 = min(f,c__), d__1 = min(d__1,g);
#line 346 "sgebal.f"
	if (g < r__ || max(r__,ra) >= sfmax2 || min(d__1,ca) <= sfmin2) {
#line 346 "sgebal.f"
	    goto L190;
#line 346 "sgebal.f"
	}
#line 348 "sgebal.f"
	d__1 = c__ + f + ca + r__ + g + ra;
#line 348 "sgebal.f"
	if (sisnan_(&d__1)) {

/*           Exit if NaN to avoid infinite loop */

#line 352 "sgebal.f"
	    *info = -3;
#line 353 "sgebal.f"
	    i__2 = -(*info);
#line 353 "sgebal.f"
	    xerbla_("SGEBAL", &i__2, (ftnlen)6);
#line 354 "sgebal.f"
	    return 0;
#line 355 "sgebal.f"
	}
#line 356 "sgebal.f"
	f /= 2.;
#line 357 "sgebal.f"
	c__ /= 2.;
#line 358 "sgebal.f"
	g /= 2.;
#line 359 "sgebal.f"
	ca /= 2.;
#line 360 "sgebal.f"
	r__ *= 2.;
#line 361 "sgebal.f"
	ra *= 2.;
#line 362 "sgebal.f"
	goto L180;

/*        Now balance. */

#line 366 "sgebal.f"
L190:
#line 367 "sgebal.f"
	if (c__ + r__ >= s * .95) {
#line 367 "sgebal.f"
	    goto L200;
#line 367 "sgebal.f"
	}
#line 369 "sgebal.f"
	if (f < 1. && scale[i__] < 1.) {
#line 370 "sgebal.f"
	    if (f * scale[i__] <= sfmin1) {
#line 370 "sgebal.f"
		goto L200;
#line 370 "sgebal.f"
	    }
#line 372 "sgebal.f"
	}
#line 373 "sgebal.f"
	if (f > 1. && scale[i__] > 1.) {
#line 374 "sgebal.f"
	    if (scale[i__] >= sfmax1 / f) {
#line 374 "sgebal.f"
		goto L200;
#line 374 "sgebal.f"
	    }
#line 376 "sgebal.f"
	}
#line 377 "sgebal.f"
	g = 1. / f;
#line 378 "sgebal.f"
	scale[i__] *= f;
#line 379 "sgebal.f"
	noconv = TRUE_;

#line 381 "sgebal.f"
	i__2 = *n - k + 1;
#line 381 "sgebal.f"
	sscal_(&i__2, &g, &a[i__ + k * a_dim1], lda);
#line 382 "sgebal.f"
	sscal_(&l, &f, &a[i__ * a_dim1 + 1], &c__1);

#line 384 "sgebal.f"
L200:
#line 384 "sgebal.f"
	;
#line 384 "sgebal.f"
    }

#line 386 "sgebal.f"
    if (noconv) {
#line 386 "sgebal.f"
	goto L140;
#line 386 "sgebal.f"
    }

#line 389 "sgebal.f"
L210:
#line 390 "sgebal.f"
    *ilo = k;
#line 391 "sgebal.f"
    *ihi = l;

#line 393 "sgebal.f"
    return 0;

/*     End of SGEBAL */

} /* sgebal_ */

