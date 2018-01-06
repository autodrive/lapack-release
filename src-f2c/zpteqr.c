#line 1 "zpteqr.f"
/* zpteqr.f -- translated by f2c (version 20100827).
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

#line 1 "zpteqr.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__0 = 0;
static integer c__1 = 1;

/* > \brief \b ZPTEQR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZPTEQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zpteqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zpteqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zpteqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZPTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPZ */
/*       INTEGER            INFO, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), E( * ), WORK( * ) */
/*       COMPLEX*16         Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPTEQR computes all eigenvalues and, optionally, eigenvectors of a */
/* > symmetric positive definite tridiagonal matrix by first factoring the */
/* > matrix using DPTTRF and then calling ZBDSQR to compute the singular */
/* > values of the bidiagonal factor. */
/* > */
/* > This routine computes the eigenvalues of the positive definite */
/* > tridiagonal matrix to high relative accuracy.  This means that if the */
/* > eigenvalues range over many orders of magnitude in size, then the */
/* > small eigenvalues and corresponding eigenvectors will be computed */
/* > more accurately than, for example, with the standard QR method. */
/* > */
/* > The eigenvectors of a full or band positive definite Hermitian matrix */
/* > can also be found if ZHETRD, ZHPTRD, or ZHBTRD has been used to */
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
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the n diagonal elements of the tridiagonal matrix. */
/* >          On normal exit, D contains the eigenvalues, in descending */
/* >          order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* >          matrix. */
/* >          On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is COMPLEX*16 array, dimension (LDZ, N) */
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
/* >          WORK is DOUBLE PRECISION array, dimension (4*N) */
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

/* > \ingroup complex16PTcomputational */

/*  ===================================================================== */
/* Subroutine */ int zpteqr_(char *compz, integer *n, doublereal *d__, 
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
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer icompz;
    extern /* Subroutine */ int zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen), dpttrf_(integer *, doublereal *, doublereal *, integer *)
	    , zbdsqr_(char *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 189 "zpteqr.f"
    /* Parameter adjustments */
#line 189 "zpteqr.f"
    --d__;
#line 189 "zpteqr.f"
    --e;
#line 189 "zpteqr.f"
    z_dim1 = *ldz;
#line 189 "zpteqr.f"
    z_offset = 1 + z_dim1;
#line 189 "zpteqr.f"
    z__ -= z_offset;
#line 189 "zpteqr.f"
    --work;
#line 189 "zpteqr.f"

#line 189 "zpteqr.f"
    /* Function Body */
#line 189 "zpteqr.f"
    *info = 0;

#line 191 "zpteqr.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 192 "zpteqr.f"
	icompz = 0;
#line 193 "zpteqr.f"
    } else if (lsame_(compz, "V", (ftnlen)1, (ftnlen)1)) {
#line 194 "zpteqr.f"
	icompz = 1;
#line 195 "zpteqr.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 196 "zpteqr.f"
	icompz = 2;
#line 197 "zpteqr.f"
    } else {
#line 198 "zpteqr.f"
	icompz = -1;
#line 199 "zpteqr.f"
    }
#line 200 "zpteqr.f"
    if (icompz < 0) {
#line 201 "zpteqr.f"
	*info = -1;
#line 202 "zpteqr.f"
    } else if (*n < 0) {
#line 203 "zpteqr.f"
	*info = -2;
#line 204 "zpteqr.f"
    } else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n)) {
#line 206 "zpteqr.f"
	*info = -6;
#line 207 "zpteqr.f"
    }
#line 208 "zpteqr.f"
    if (*info != 0) {
#line 209 "zpteqr.f"
	i__1 = -(*info);
#line 209 "zpteqr.f"
	xerbla_("ZPTEQR", &i__1, (ftnlen)6);
#line 210 "zpteqr.f"
	return 0;
#line 211 "zpteqr.f"
    }

/*     Quick return if possible */

#line 215 "zpteqr.f"
    if (*n == 0) {
#line 215 "zpteqr.f"
	return 0;
#line 215 "zpteqr.f"
    }

#line 218 "zpteqr.f"
    if (*n == 1) {
#line 219 "zpteqr.f"
	if (icompz > 0) {
#line 219 "zpteqr.f"
	    i__1 = z_dim1 + 1;
#line 219 "zpteqr.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 219 "zpteqr.f"
	}
#line 221 "zpteqr.f"
	return 0;
#line 222 "zpteqr.f"
    }
#line 223 "zpteqr.f"
    if (icompz == 2) {
#line 223 "zpteqr.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &z__[z_offset], ldz, (ftnlen)4);
#line 223 "zpteqr.f"
    }

/*     Call DPTTRF to factor the matrix. */

#line 228 "zpteqr.f"
    dpttrf_(n, &d__[1], &e[1], info);
#line 229 "zpteqr.f"
    if (*info != 0) {
#line 229 "zpteqr.f"
	return 0;
#line 229 "zpteqr.f"
    }
#line 231 "zpteqr.f"
    i__1 = *n;
#line 231 "zpteqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 232 "zpteqr.f"
	d__[i__] = sqrt(d__[i__]);
#line 233 "zpteqr.f"
/* L10: */
#line 233 "zpteqr.f"
    }
#line 234 "zpteqr.f"
    i__1 = *n - 1;
#line 234 "zpteqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 235 "zpteqr.f"
	e[i__] *= d__[i__];
#line 236 "zpteqr.f"
/* L20: */
#line 236 "zpteqr.f"
    }

/*     Call ZBDSQR to compute the singular values/vectors of the */
/*     bidiagonal factor. */

#line 241 "zpteqr.f"
    if (icompz > 0) {
#line 242 "zpteqr.f"
	nru = *n;
#line 243 "zpteqr.f"
    } else {
#line 244 "zpteqr.f"
	nru = 0;
#line 245 "zpteqr.f"
    }
#line 246 "zpteqr.f"
    zbdsqr_("Lower", n, &c__0, &nru, &c__0, &d__[1], &e[1], vt, &c__1, &z__[
	    z_offset], ldz, c__, &c__1, &work[1], info, (ftnlen)5);

/*     Square the singular values. */

#line 251 "zpteqr.f"
    if (*info == 0) {
#line 252 "zpteqr.f"
	i__1 = *n;
#line 252 "zpteqr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 253 "zpteqr.f"
	    d__[i__] *= d__[i__];
#line 254 "zpteqr.f"
/* L30: */
#line 254 "zpteqr.f"
	}
#line 255 "zpteqr.f"
    } else {
#line 256 "zpteqr.f"
	*info = *n + *info;
#line 257 "zpteqr.f"
    }

#line 259 "zpteqr.f"
    return 0;

/*     End of ZPTEQR */

} /* zpteqr_ */

