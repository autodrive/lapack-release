#line 1 "dstev.f"
/* dstev.f -- translated by f2c (version 20100827).
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

#line 1 "dstev.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> DSTEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER m
atrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSTEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dstev.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dstev.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dstev.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ */
/*       INTEGER            INFO, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSTEV computes all eigenvalues and, optionally, eigenvectors of a */
/* > real symmetric tridiagonal matrix A. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBZ */
/* > \verbatim */
/* >          JOBZ is CHARACTER*1 */
/* >          = 'N':  Compute eigenvalues only; */
/* >          = 'V':  Compute eigenvalues and eigenvectors. */
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
/* >          On entry, the n diagonal elements of the tridiagonal matrix */
/* >          A. */
/* >          On exit, if INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* >          matrix A, stored in elements 1 to N-1 of E. */
/* >          On exit, the contents of E are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, N) */
/* >          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal */
/* >          eigenvectors of the matrix A, with the i-th column of Z */
/* >          holding the eigenvector associated with D(i). */
/* >          If JOBZ = 'N', then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z.  LDZ >= 1, and if */
/* >          JOBZ = 'V', LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (max(1,2*N-2)) */
/* >          If JOBZ = 'N', WORK is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, the algorithm failed to converge; i */
/* >                off-diagonal elements of E did not converge to zero. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int dstev_(char *jobz, integer *n, doublereal *d__, 
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	integer *info, ftnlen jobz_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal eps;
    static integer imax;
    static doublereal rmin, rmax, tnrm;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical wantz;
    extern doublereal dlamch_(char *, ftnlen);
    static integer iscale;
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *), dsteqr_(char *, integer *, doublereal *, doublereal *
	    , doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal smlnum;


/*  -- LAPACK driver routine (version 3.7.0) -- */
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 159 "dstev.f"
    /* Parameter adjustments */
#line 159 "dstev.f"
    --d__;
#line 159 "dstev.f"
    --e;
#line 159 "dstev.f"
    z_dim1 = *ldz;
#line 159 "dstev.f"
    z_offset = 1 + z_dim1;
#line 159 "dstev.f"
    z__ -= z_offset;
#line 159 "dstev.f"
    --work;
#line 159 "dstev.f"

#line 159 "dstev.f"
    /* Function Body */
#line 159 "dstev.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);

#line 161 "dstev.f"
    *info = 0;
#line 162 "dstev.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 163 "dstev.f"
	*info = -1;
#line 164 "dstev.f"
    } else if (*n < 0) {
#line 165 "dstev.f"
	*info = -2;
#line 166 "dstev.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 167 "dstev.f"
	*info = -6;
#line 168 "dstev.f"
    }

#line 170 "dstev.f"
    if (*info != 0) {
#line 171 "dstev.f"
	i__1 = -(*info);
#line 171 "dstev.f"
	xerbla_("DSTEV ", &i__1, (ftnlen)6);
#line 172 "dstev.f"
	return 0;
#line 173 "dstev.f"
    }

/*     Quick return if possible */

#line 177 "dstev.f"
    if (*n == 0) {
#line 177 "dstev.f"
	return 0;
#line 177 "dstev.f"
    }

#line 180 "dstev.f"
    if (*n == 1) {
#line 181 "dstev.f"
	if (wantz) {
#line 181 "dstev.f"
	    z__[z_dim1 + 1] = 1.;
#line 181 "dstev.f"
	}
#line 183 "dstev.f"
	return 0;
#line 184 "dstev.f"
    }

/*     Get machine constants. */

#line 188 "dstev.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 189 "dstev.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 190 "dstev.f"
    smlnum = safmin / eps;
#line 191 "dstev.f"
    bignum = 1. / smlnum;
#line 192 "dstev.f"
    rmin = sqrt(smlnum);
#line 193 "dstev.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 197 "dstev.f"
    iscale = 0;
#line 198 "dstev.f"
    tnrm = dlanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 199 "dstev.f"
    if (tnrm > 0. && tnrm < rmin) {
#line 200 "dstev.f"
	iscale = 1;
#line 201 "dstev.f"
	sigma = rmin / tnrm;
#line 202 "dstev.f"
    } else if (tnrm > rmax) {
#line 203 "dstev.f"
	iscale = 1;
#line 204 "dstev.f"
	sigma = rmax / tnrm;
#line 205 "dstev.f"
    }
#line 206 "dstev.f"
    if (iscale == 1) {
#line 207 "dstev.f"
	dscal_(n, &sigma, &d__[1], &c__1);
#line 208 "dstev.f"
	i__1 = *n - 1;
#line 208 "dstev.f"
	dscal_(&i__1, &sigma, &e[1], &c__1);
#line 209 "dstev.f"
    }

/*     For eigenvalues only, call DSTERF.  For eigenvalues and */
/*     eigenvectors, call DSTEQR. */

#line 214 "dstev.f"
    if (! wantz) {
#line 215 "dstev.f"
	dsterf_(n, &d__[1], &e[1], info);
#line 216 "dstev.f"
    } else {
#line 217 "dstev.f"
	dsteqr_("I", n, &d__[1], &e[1], &z__[z_offset], ldz, &work[1], info, (
		ftnlen)1);
#line 218 "dstev.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 222 "dstev.f"
    if (iscale == 1) {
#line 223 "dstev.f"
	if (*info == 0) {
#line 224 "dstev.f"
	    imax = *n;
#line 225 "dstev.f"
	} else {
#line 226 "dstev.f"
	    imax = *info - 1;
#line 227 "dstev.f"
	}
#line 228 "dstev.f"
	d__1 = 1. / sigma;
#line 228 "dstev.f"
	dscal_(&imax, &d__1, &d__[1], &c__1);
#line 229 "dstev.f"
    }

#line 231 "dstev.f"
    return 0;

/*     End of DSTEV */

} /* dstev_ */

