#line 1 "spteqr.f"
/* spteqr.f -- translated by f2c (version 20100827).
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

#line 1 "spteqr.f"
/* Table of constant values */

static doublereal c_b7 = 0.;
static doublereal c_b8 = 1.;
static integer c__0 = 0;
static integer c__1 = 1;

/* > \brief \b SPTEQR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SPTEQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spteqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spteqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spteqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SPTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPZ */
/*       INTEGER            INFO, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPTEQR computes all eigenvalues and, optionally, eigenvectors of a */
/* > symmetric positive definite tridiagonal matrix by first factoring the */
/* > matrix using SPTTRF, and then calling SBDSQR to compute the singular */
/* > values of the bidiagonal factor. */
/* > */
/* > This routine computes the eigenvalues of the positive definite */
/* > tridiagonal matrix to high relative accuracy.  This means that if the */
/* > eigenvalues range over many orders of magnitude in size, then the */
/* > small eigenvalues and corresponding eigenvectors will be computed */
/* > more accurately than, for example, with the standard QR method. */
/* > */
/* > The eigenvectors of a full or band symmetric positive definite matrix */
/* > can also be found if SSYTRD, SSPTRD, or SSBTRD has been used to */
/* > reduce this matrix to tridiagonal form. (The reduction to tridiagonal */
/* > form, however, may preclude the possibility of obtaining high */
/* > relative accuracy in the small eigenvalues of the original matrix, if */
/* > these eigenvalues range over many orders of magnitude.) */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] COMPZ */
/* > \verbatim */
/* >          COMPZ is CHARACTER*1 */
/* >          = 'N':  Compute eigenvalues only. */
/* >          = 'V':  Compute eigenvectors of original symmetric */
/* >                  matrix also.  Array Z contains the orthogonal */
/* >                  matrix used to reduce the original matrix to */
/* >                  tridiagonal form. */
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
/* >          On entry, the n diagonal elements of the tridiagonal */
/* >          matrix. */
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
/* >          Z is REAL array, dimension (LDZ, N) */
/* >          On entry, if COMPZ = 'V', the orthogonal matrix used in the */
/* >          reduction to tridiagonal form. */
/* >          On exit, if COMPZ = 'V', the orthonormal eigenvectors of the */
/* >          original symmetric matrix; */
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

/* > \date December 2016 */

/* > \ingroup realPTcomputational */

/*  ===================================================================== */
/* Subroutine */ int spteqr_(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	integer *info, ftnlen compz_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c__[1]	/* was [1][1] */;
    static integer i__;
    static doublereal vt[1]	/* was [1][1] */;
    static integer nru;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), slaset_(
	    char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), sbdsqr_(char *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen);
    static integer icompz;
    extern /* Subroutine */ int spttrf_(integer *, doublereal *, doublereal *,
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

#line 187 "spteqr.f"
    /* Parameter adjustments */
#line 187 "spteqr.f"
    --d__;
#line 187 "spteqr.f"
    --e;
#line 187 "spteqr.f"
    z_dim1 = *ldz;
#line 187 "spteqr.f"
    z_offset = 1 + z_dim1;
#line 187 "spteqr.f"
    z__ -= z_offset;
#line 187 "spteqr.f"
    --work;
#line 187 "spteqr.f"

#line 187 "spteqr.f"
    /* Function Body */
#line 187 "spteqr.f"
    *info = 0;

#line 189 "spteqr.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 190 "spteqr.f"
	icompz = 0;
#line 191 "spteqr.f"
    } else if (lsame_(compz, "V", (ftnlen)1, (ftnlen)1)) {
#line 192 "spteqr.f"
	icompz = 1;
#line 193 "spteqr.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 194 "spteqr.f"
	icompz = 2;
#line 195 "spteqr.f"
    } else {
#line 196 "spteqr.f"
	icompz = -1;
#line 197 "spteqr.f"
    }
#line 198 "spteqr.f"
    if (icompz < 0) {
#line 199 "spteqr.f"
	*info = -1;
#line 200 "spteqr.f"
    } else if (*n < 0) {
#line 201 "spteqr.f"
	*info = -2;
#line 202 "spteqr.f"
    } else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n)) {
#line 204 "spteqr.f"
	*info = -6;
#line 205 "spteqr.f"
    }
#line 206 "spteqr.f"
    if (*info != 0) {
#line 207 "spteqr.f"
	i__1 = -(*info);
#line 207 "spteqr.f"
	xerbla_("SPTEQR", &i__1, (ftnlen)6);
#line 208 "spteqr.f"
	return 0;
#line 209 "spteqr.f"
    }

/*     Quick return if possible */

#line 213 "spteqr.f"
    if (*n == 0) {
#line 213 "spteqr.f"
	return 0;
#line 213 "spteqr.f"
    }

#line 216 "spteqr.f"
    if (*n == 1) {
#line 217 "spteqr.f"
	if (icompz > 0) {
#line 217 "spteqr.f"
	    z__[z_dim1 + 1] = 1.;
#line 217 "spteqr.f"
	}
#line 219 "spteqr.f"
	return 0;
#line 220 "spteqr.f"
    }
#line 221 "spteqr.f"
    if (icompz == 2) {
#line 221 "spteqr.f"
	slaset_("Full", n, n, &c_b7, &c_b8, &z__[z_offset], ldz, (ftnlen)4);
#line 221 "spteqr.f"
    }

/*     Call SPTTRF to factor the matrix. */

#line 226 "spteqr.f"
    spttrf_(n, &d__[1], &e[1], info);
#line 227 "spteqr.f"
    if (*info != 0) {
#line 227 "spteqr.f"
	return 0;
#line 227 "spteqr.f"
    }
#line 229 "spteqr.f"
    i__1 = *n;
#line 229 "spteqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 230 "spteqr.f"
	d__[i__] = sqrt(d__[i__]);
#line 231 "spteqr.f"
/* L10: */
#line 231 "spteqr.f"
    }
#line 232 "spteqr.f"
    i__1 = *n - 1;
#line 232 "spteqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 233 "spteqr.f"
	e[i__] *= d__[i__];
#line 234 "spteqr.f"
/* L20: */
#line 234 "spteqr.f"
    }

/*     Call SBDSQR to compute the singular values/vectors of the */
/*     bidiagonal factor. */

#line 239 "spteqr.f"
    if (icompz > 0) {
#line 240 "spteqr.f"
	nru = *n;
#line 241 "spteqr.f"
    } else {
#line 242 "spteqr.f"
	nru = 0;
#line 243 "spteqr.f"
    }
#line 244 "spteqr.f"
    sbdsqr_("Lower", n, &c__0, &nru, &c__0, &d__[1], &e[1], vt, &c__1, &z__[
	    z_offset], ldz, c__, &c__1, &work[1], info, (ftnlen)5);

/*     Square the singular values. */

#line 249 "spteqr.f"
    if (*info == 0) {
#line 250 "spteqr.f"
	i__1 = *n;
#line 250 "spteqr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 251 "spteqr.f"
	    d__[i__] *= d__[i__];
#line 252 "spteqr.f"
/* L30: */
#line 252 "spteqr.f"
	}
#line 253 "spteqr.f"
    } else {
#line 254 "spteqr.f"
	*info = *n + *info;
#line 255 "spteqr.f"
    }

#line 257 "spteqr.f"
    return 0;

/*     End of SPTEQR */

} /* spteqr_ */

