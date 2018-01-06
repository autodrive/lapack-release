#line 1 "cpteqr.f"
/* cpteqr.f -- translated by f2c (version 20100827).
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

#line 1 "cpteqr.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__0 = 0;
static integer c__1 = 1;

/* > \brief \b CPTEQR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPTEQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpteqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpteqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpteqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPZ */
/*       INTEGER            INFO, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), E( * ), WORK( * ) */
/*       COMPLEX            Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPTEQR computes all eigenvalues and, optionally, eigenvectors of a */
/* > symmetric positive definite tridiagonal matrix by first factoring the */
/* > matrix using SPTTRF and then calling CBDSQR to compute the singular */
/* > values of the bidiagonal factor. */
/* > */
/* > This routine computes the eigenvalues of the positive definite */
/* > tridiagonal matrix to high relative accuracy.  This means that if the */
/* > eigenvalues range over many orders of magnitude in size, then the */
/* > small eigenvalues and corresponding eigenvectors will be computed */
/* > more accurately than, for example, with the standard QR method. */
/* > */
/* > The eigenvectors of a full or band positive definite Hermitian matrix */
/* > can also be found if CHETRD, CHPTRD, or CHBTRD has been used to */
/* > reduce this matrix to tridiagonal form.  (The reduction to */
/* > tridiagonal form, however, may preclude the possibility of obtaining */
/* > high relative accuracy in the small eigenvalues of the original */
/* > matrix, if these eigenvalues range over many orders of magnitude.) */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] COMPZ */
/* > \verbatim */
/* >          COMPZ is CHARACTER*1 */
/* >          = 'N':  Compute eigenvalues only. */
/* >          = 'V':  Compute eigenvectors of original Hermitian */
/* >                  matrix also.  Array Z contains the unitary matrix */
/* >                  used to reduce the original matrix to tridiagonal */
/* >                  form. */
/* >          = 'I':  Compute eigenvectors of tridiagonal matrix also. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          On entry, the n diagonal elements of the tridiagonal matrix. */
/* >          On normal exit, D contains the eigenvalues, in descending */
/* >          order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N-1) */
/* >          On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* >          matrix. */
/* >          On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is COMPLEX array, dimension (LDZ, N) */
/* >          On entry, if COMPZ = 'V', the unitary matrix used in the */
/* >          reduction to tridiagonal form. */
/* >          On exit, if COMPZ = 'V', the orthonormal eigenvectors of the */
/* >          original Hermitian matrix; */
/* >          if COMPZ = 'I', the orthonormal eigenvectors of the */
/* >          tridiagonal matrix. */
/* >          If INFO > 0 on exit, Z contains the eigenvectors associated */
/* >          with only the stored eigenvalues. */
/* >          If  COMPZ = 'N', then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z.  LDZ >= 1, and if */
/* >          COMPZ = 'V' or 'I', LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (4*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if INFO = i, and i is: */
/* >                <= N  the Cholesky factorization of the matrix could */
/* >                      not be performed because the i-th principal minor */
/* >                      was not positive definite. */
/* >                > N   the SVD algorithm failed to converge; */
/* >                      if INFO = N+i, i off-diagonal elements of the */
/* >                      bidiagonal factor did not converge to zero. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complexPTcomputational */

/*  ===================================================================== */
/* Subroutine */ int cpteqr_(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublecomplex *z__, integer *ldz, doublereal *work, 
	integer *info, ftnlen compz_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublecomplex c__[1]	/* was [1][1] */;
    static integer i__;
    static doublecomplex vt[1]	/* was [1][1] */;
    static integer nru;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int claset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen), cbdsqr_(char *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublereal *, integer *, 
	    ftnlen);
    static integer icompz;
    extern /* Subroutine */ int spttrf_(integer *, doublereal *, doublereal *,
	     integer *);


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ==================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 189 "cpteqr.f"
    /* Parameter adjustments */
#line 189 "cpteqr.f"
    --d__;
#line 189 "cpteqr.f"
    --e;
#line 189 "cpteqr.f"
    z_dim1 = *ldz;
#line 189 "cpteqr.f"
    z_offset = 1 + z_dim1;
#line 189 "cpteqr.f"
    z__ -= z_offset;
#line 189 "cpteqr.f"
    --work;
#line 189 "cpteqr.f"

#line 189 "cpteqr.f"
    /* Function Body */
#line 189 "cpteqr.f"
    *info = 0;

#line 191 "cpteqr.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 192 "cpteqr.f"
	icompz = 0;
#line 193 "cpteqr.f"
    } else if (lsame_(compz, "V", (ftnlen)1, (ftnlen)1)) {
#line 194 "cpteqr.f"
	icompz = 1;
#line 195 "cpteqr.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 196 "cpteqr.f"
	icompz = 2;
#line 197 "cpteqr.f"
    } else {
#line 198 "cpteqr.f"
	icompz = -1;
#line 199 "cpteqr.f"
    }
#line 200 "cpteqr.f"
    if (icompz < 0) {
#line 201 "cpteqr.f"
	*info = -1;
#line 202 "cpteqr.f"
    } else if (*n < 0) {
#line 203 "cpteqr.f"
	*info = -2;
#line 204 "cpteqr.f"
    } else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n)) {
#line 206 "cpteqr.f"
	*info = -6;
#line 207 "cpteqr.f"
    }
#line 208 "cpteqr.f"
    if (*info != 0) {
#line 209 "cpteqr.f"
	i__1 = -(*info);
#line 209 "cpteqr.f"
	xerbla_("CPTEQR", &i__1, (ftnlen)6);
#line 210 "cpteqr.f"
	return 0;
#line 211 "cpteqr.f"
    }

/*     Quick return if possible */

#line 215 "cpteqr.f"
    if (*n == 0) {
#line 215 "cpteqr.f"
	return 0;
#line 215 "cpteqr.f"
    }

#line 218 "cpteqr.f"
    if (*n == 1) {
#line 219 "cpteqr.f"
	if (icompz > 0) {
#line 219 "cpteqr.f"
	    i__1 = z_dim1 + 1;
#line 219 "cpteqr.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 219 "cpteqr.f"
	}
#line 221 "cpteqr.f"
	return 0;
#line 222 "cpteqr.f"
    }
#line 223 "cpteqr.f"
    if (icompz == 2) {
#line 223 "cpteqr.f"
	claset_("Full", n, n, &c_b1, &c_b2, &z__[z_offset], ldz, (ftnlen)4);
#line 223 "cpteqr.f"
    }

/*     Call SPTTRF to factor the matrix. */

#line 228 "cpteqr.f"
    spttrf_(n, &d__[1], &e[1], info);
#line 229 "cpteqr.f"
    if (*info != 0) {
#line 229 "cpteqr.f"
	return 0;
#line 229 "cpteqr.f"
    }
#line 231 "cpteqr.f"
    i__1 = *n;
#line 231 "cpteqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 232 "cpteqr.f"
	d__[i__] = sqrt(d__[i__]);
#line 233 "cpteqr.f"
/* L10: */
#line 233 "cpteqr.f"
    }
#line 234 "cpteqr.f"
    i__1 = *n - 1;
#line 234 "cpteqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 235 "cpteqr.f"
	e[i__] *= d__[i__];
#line 236 "cpteqr.f"
/* L20: */
#line 236 "cpteqr.f"
    }

/*     Call CBDSQR to compute the singular values/vectors of the */
/*     bidiagonal factor. */

#line 241 "cpteqr.f"
    if (icompz > 0) {
#line 242 "cpteqr.f"
	nru = *n;
#line 243 "cpteqr.f"
    } else {
#line 244 "cpteqr.f"
	nru = 0;
#line 245 "cpteqr.f"
    }
#line 246 "cpteqr.f"
    cbdsqr_("Lower", n, &c__0, &nru, &c__0, &d__[1], &e[1], vt, &c__1, &z__[
	    z_offset], ldz, c__, &c__1, &work[1], info, (ftnlen)5);

/*     Square the singular values. */

#line 251 "cpteqr.f"
    if (*info == 0) {
#line 252 "cpteqr.f"
	i__1 = *n;
#line 252 "cpteqr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 253 "cpteqr.f"
	    d__[i__] *= d__[i__];
#line 254 "cpteqr.f"
/* L30: */
#line 254 "cpteqr.f"
	}
#line 255 "cpteqr.f"
    } else {
#line 256 "cpteqr.f"
	*info = *n + *info;
#line 257 "cpteqr.f"
    }

#line 259 "cpteqr.f"
    return 0;

/*     End of CPTEQR */

} /* cpteqr_ */

