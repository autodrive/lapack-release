#line 1 "cgebal.f"
/* cgebal.f -- translated by f2c (version 20100827).
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

#line 1 "cgebal.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CGEBAL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGEBAL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgebal.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgebal.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgebal.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOB */
/*       INTEGER            IHI, ILO, INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               SCALE( * ) */
/*       COMPLEX            A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEBAL balances a general complex matrix A.  This involves, first, */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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

/* > \date November 2013 */

/* > \ingroup complexGEcomputational */

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
/* Subroutine */ int cgebal_(char *job, integer *n, doublecomplex *a, integer 
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
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static doublereal sfmin1, sfmin2, sfmax1, sfmax2;
    extern doublereal scnrm2_(integer *, doublecomplex *, integer *);
    extern integer icamax_(integer *, doublecomplex *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), xerbla_(char *, integer *, ftnlen);
    extern logical sisnan_(doublereal *);
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 217 "cgebal.f"
    /* Parameter adjustments */
#line 217 "cgebal.f"
    a_dim1 = *lda;
#line 217 "cgebal.f"
    a_offset = 1 + a_dim1;
#line 217 "cgebal.f"
    a -= a_offset;
#line 217 "cgebal.f"
    --scale;
#line 217 "cgebal.f"

#line 217 "cgebal.f"
    /* Function Body */
#line 217 "cgebal.f"
    *info = 0;
#line 218 "cgebal.f"
    if (! lsame_(job, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(job, "P", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(job, "S", (ftnlen)1, (ftnlen)1) 
	    && ! lsame_(job, "B", (ftnlen)1, (ftnlen)1)) {
#line 220 "cgebal.f"
	*info = -1;
#line 221 "cgebal.f"
    } else if (*n < 0) {
#line 222 "cgebal.f"
	*info = -2;
#line 223 "cgebal.f"
    } else if (*lda < max(1,*n)) {
#line 224 "cgebal.f"
	*info = -4;
#line 225 "cgebal.f"
    }
#line 226 "cgebal.f"
    if (*info != 0) {
#line 227 "cgebal.f"
	i__1 = -(*info);
#line 227 "cgebal.f"
	xerbla_("CGEBAL", &i__1, (ftnlen)6);
#line 228 "cgebal.f"
	return 0;
#line 229 "cgebal.f"
    }

#line 231 "cgebal.f"
    k = 1;
#line 232 "cgebal.f"
    l = *n;

#line 234 "cgebal.f"
    if (*n == 0) {
#line 234 "cgebal.f"
	goto L210;
#line 234 "cgebal.f"
    }

#line 237 "cgebal.f"
    if (lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 238 "cgebal.f"
	i__1 = *n;
#line 238 "cgebal.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 239 "cgebal.f"
	    scale[i__] = 1.;
#line 240 "cgebal.f"
/* L10: */
#line 240 "cgebal.f"
	}
#line 241 "cgebal.f"
	goto L210;
#line 242 "cgebal.f"
    }

#line 244 "cgebal.f"
    if (lsame_(job, "S", (ftnlen)1, (ftnlen)1)) {
#line 244 "cgebal.f"
	goto L120;
#line 244 "cgebal.f"
    }

/*     Permutation to isolate eigenvalues if possible */

#line 249 "cgebal.f"
    goto L50;

/*     Row and column exchange. */

#line 253 "cgebal.f"
L20:
#line 254 "cgebal.f"
    scale[m] = (doublereal) j;
#line 255 "cgebal.f"
    if (j == m) {
#line 255 "cgebal.f"
	goto L30;
#line 255 "cgebal.f"
    }

#line 258 "cgebal.f"
    cswap_(&l, &a[j * a_dim1 + 1], &c__1, &a[m * a_dim1 + 1], &c__1);
#line 259 "cgebal.f"
    i__1 = *n - k + 1;
#line 259 "cgebal.f"
    cswap_(&i__1, &a[j + k * a_dim1], lda, &a[m + k * a_dim1], lda);

#line 261 "cgebal.f"
L30:
#line 262 "cgebal.f"
    switch (iexc) {
#line 262 "cgebal.f"
	case 1:  goto L40;
#line 262 "cgebal.f"
	case 2:  goto L80;
#line 262 "cgebal.f"
    }

/*     Search for rows isolating an eigenvalue and push them down. */

#line 266 "cgebal.f"
L40:
#line 267 "cgebal.f"
    if (l == 1) {
#line 267 "cgebal.f"
	goto L210;
#line 267 "cgebal.f"
    }
#line 269 "cgebal.f"
    --l;

#line 271 "cgebal.f"
L50:
#line 272 "cgebal.f"
    for (j = l; j >= 1; --j) {

#line 274 "cgebal.f"
	i__1 = l;
#line 274 "cgebal.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 275 "cgebal.f"
	    if (i__ == j) {
#line 275 "cgebal.f"
		goto L60;
#line 275 "cgebal.f"
	    }
#line 277 "cgebal.f"
	    i__2 = j + i__ * a_dim1;
#line 277 "cgebal.f"
	    if (a[i__2].r != 0. || d_imag(&a[j + i__ * a_dim1]) != 0.) {
#line 277 "cgebal.f"
		goto L70;
#line 277 "cgebal.f"
	    }
#line 279 "cgebal.f"
L60:
#line 279 "cgebal.f"
	    ;
#line 279 "cgebal.f"
	}

#line 281 "cgebal.f"
	m = l;
#line 282 "cgebal.f"
	iexc = 1;
#line 283 "cgebal.f"
	goto L20;
#line 284 "cgebal.f"
L70:
#line 284 "cgebal.f"
	;
#line 284 "cgebal.f"
    }

#line 286 "cgebal.f"
    goto L90;

/*     Search for columns isolating an eigenvalue and push them left. */

#line 290 "cgebal.f"
L80:
#line 291 "cgebal.f"
    ++k;

#line 293 "cgebal.f"
L90:
#line 294 "cgebal.f"
    i__1 = l;
#line 294 "cgebal.f"
    for (j = k; j <= i__1; ++j) {

#line 296 "cgebal.f"
	i__2 = l;
#line 296 "cgebal.f"
	for (i__ = k; i__ <= i__2; ++i__) {
#line 297 "cgebal.f"
	    if (i__ == j) {
#line 297 "cgebal.f"
		goto L100;
#line 297 "cgebal.f"
	    }
#line 299 "cgebal.f"
	    i__3 = i__ + j * a_dim1;
#line 299 "cgebal.f"
	    if (a[i__3].r != 0. || d_imag(&a[i__ + j * a_dim1]) != 0.) {
#line 299 "cgebal.f"
		goto L110;
#line 299 "cgebal.f"
	    }
#line 301 "cgebal.f"
L100:
#line 301 "cgebal.f"
	    ;
#line 301 "cgebal.f"
	}

#line 303 "cgebal.f"
	m = k;
#line 304 "cgebal.f"
	iexc = 2;
#line 305 "cgebal.f"
	goto L20;
#line 306 "cgebal.f"
L110:
#line 306 "cgebal.f"
	;
#line 306 "cgebal.f"
    }

#line 308 "cgebal.f"
L120:
#line 309 "cgebal.f"
    i__1 = l;
#line 309 "cgebal.f"
    for (i__ = k; i__ <= i__1; ++i__) {
#line 310 "cgebal.f"
	scale[i__] = 1.;
#line 311 "cgebal.f"
/* L130: */
#line 311 "cgebal.f"
    }

#line 313 "cgebal.f"
    if (lsame_(job, "P", (ftnlen)1, (ftnlen)1)) {
#line 313 "cgebal.f"
	goto L210;
#line 313 "cgebal.f"
    }

/*     Balance the submatrix in rows K to L. */

/*     Iterative loop for norm reduction */

#line 320 "cgebal.f"
    sfmin1 = slamch_("S", (ftnlen)1) / slamch_("P", (ftnlen)1);
#line 321 "cgebal.f"
    sfmax1 = 1. / sfmin1;
#line 322 "cgebal.f"
    sfmin2 = sfmin1 * 2.;
#line 323 "cgebal.f"
    sfmax2 = 1. / sfmin2;
#line 324 "cgebal.f"
L140:
#line 325 "cgebal.f"
    noconv = FALSE_;

#line 327 "cgebal.f"
    i__1 = l;
#line 327 "cgebal.f"
    for (i__ = k; i__ <= i__1; ++i__) {

#line 329 "cgebal.f"
	i__2 = l - k + 1;
#line 329 "cgebal.f"
	c__ = scnrm2_(&i__2, &a[k + i__ * a_dim1], &c__1);
#line 330 "cgebal.f"
	i__2 = l - k + 1;
#line 330 "cgebal.f"
	r__ = scnrm2_(&i__2, &a[i__ + k * a_dim1], lda);
#line 331 "cgebal.f"
	ica = icamax_(&l, &a[i__ * a_dim1 + 1], &c__1);
#line 332 "cgebal.f"
	ca = z_abs(&a[ica + i__ * a_dim1]);
#line 333 "cgebal.f"
	i__2 = *n - k + 1;
#line 333 "cgebal.f"
	ira = icamax_(&i__2, &a[i__ + k * a_dim1], lda);
#line 334 "cgebal.f"
	ra = z_abs(&a[i__ + (ira + k - 1) * a_dim1]);

/*        Guard against zero C or R due to underflow. */

#line 338 "cgebal.f"
	if (c__ == 0. || r__ == 0.) {
#line 338 "cgebal.f"
	    goto L200;
#line 338 "cgebal.f"
	}
#line 340 "cgebal.f"
	g = r__ / 2.;
#line 341 "cgebal.f"
	f = 1.;
#line 342 "cgebal.f"
	s = c__ + r__;
#line 343 "cgebal.f"
L160:
/* Computing MAX */
#line 344 "cgebal.f"
	d__1 = max(f,c__);
/* Computing MIN */
#line 344 "cgebal.f"
	d__2 = min(r__,g);
#line 344 "cgebal.f"
	if (c__ >= g || max(d__1,ca) >= sfmax2 || min(d__2,ra) <= sfmin2) {
#line 344 "cgebal.f"
	    goto L170;
#line 344 "cgebal.f"
	}
#line 346 "cgebal.f"
	d__1 = c__ + f + ca + r__ + g + ra;
#line 346 "cgebal.f"
	if (sisnan_(&d__1)) {

/*           Exit if NaN to avoid infinite loop */

#line 350 "cgebal.f"
	    *info = -3;
#line 351 "cgebal.f"
	    i__2 = -(*info);
#line 351 "cgebal.f"
	    xerbla_("CGEBAL", &i__2, (ftnlen)6);
#line 352 "cgebal.f"
	    return 0;
#line 353 "cgebal.f"
	}
#line 354 "cgebal.f"
	f *= 2.;
#line 355 "cgebal.f"
	c__ *= 2.;
#line 356 "cgebal.f"
	ca *= 2.;
#line 357 "cgebal.f"
	r__ /= 2.;
#line 358 "cgebal.f"
	g /= 2.;
#line 359 "cgebal.f"
	ra /= 2.;
#line 360 "cgebal.f"
	goto L160;

#line 362 "cgebal.f"
L170:
#line 363 "cgebal.f"
	g = c__ / 2.;
#line 364 "cgebal.f"
L180:
/* Computing MIN */
#line 365 "cgebal.f"
	d__1 = min(f,c__), d__1 = min(d__1,g);
#line 365 "cgebal.f"
	if (g < r__ || max(r__,ra) >= sfmax2 || min(d__1,ca) <= sfmin2) {
#line 365 "cgebal.f"
	    goto L190;
#line 365 "cgebal.f"
	}
#line 367 "cgebal.f"
	f /= 2.;
#line 368 "cgebal.f"
	c__ /= 2.;
#line 369 "cgebal.f"
	g /= 2.;
#line 370 "cgebal.f"
	ca /= 2.;
#line 371 "cgebal.f"
	r__ *= 2.;
#line 372 "cgebal.f"
	ra *= 2.;
#line 373 "cgebal.f"
	goto L180;

/*        Now balance. */

#line 377 "cgebal.f"
L190:
#line 378 "cgebal.f"
	if (c__ + r__ >= s * .95) {
#line 378 "cgebal.f"
	    goto L200;
#line 378 "cgebal.f"
	}
#line 380 "cgebal.f"
	if (f < 1. && scale[i__] < 1.) {
#line 381 "cgebal.f"
	    if (f * scale[i__] <= sfmin1) {
#line 381 "cgebal.f"
		goto L200;
#line 381 "cgebal.f"
	    }
#line 383 "cgebal.f"
	}
#line 384 "cgebal.f"
	if (f > 1. && scale[i__] > 1.) {
#line 385 "cgebal.f"
	    if (scale[i__] >= sfmax1 / f) {
#line 385 "cgebal.f"
		goto L200;
#line 385 "cgebal.f"
	    }
#line 387 "cgebal.f"
	}
#line 388 "cgebal.f"
	g = 1. / f;
#line 389 "cgebal.f"
	scale[i__] *= f;
#line 390 "cgebal.f"
	noconv = TRUE_;

#line 392 "cgebal.f"
	i__2 = *n - k + 1;
#line 392 "cgebal.f"
	csscal_(&i__2, &g, &a[i__ + k * a_dim1], lda);
#line 393 "cgebal.f"
	csscal_(&l, &f, &a[i__ * a_dim1 + 1], &c__1);

#line 395 "cgebal.f"
L200:
#line 395 "cgebal.f"
	;
#line 395 "cgebal.f"
    }

#line 397 "cgebal.f"
    if (noconv) {
#line 397 "cgebal.f"
	goto L140;
#line 397 "cgebal.f"
    }

#line 400 "cgebal.f"
L210:
#line 401 "cgebal.f"
    *ilo = k;
#line 402 "cgebal.f"
    *ihi = l;

#line 404 "cgebal.f"
    return 0;

/*     End of CGEBAL */

} /* cgebal_ */

