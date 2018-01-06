#line 1 "sstev.f"
/* sstev.f -- translated by f2c (version 20100827).
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

#line 1 "sstev.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> SSTEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER m
atrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSTEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstev.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstev.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstev.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ */
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
/* > SSTEV computes all eigenvalues and, optionally, eigenvectors of a */
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
/* >          D is REAL array, dimension (N) */
/* >          On entry, the n diagonal elements of the tridiagonal matrix */
/* >          A. */
/* >          On exit, if INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N-1) */
/* >          On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* >          matrix A, stored in elements 1 to N-1 of E. */
/* >          On exit, the contents of E are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (LDZ, N) */
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
/* >          WORK is REAL array, dimension (max(1,2*N-2)) */
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

/* > \date November 2011 */

/* > \ingroup realOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int sstev_(char *jobz, integer *n, doublereal *d__, 
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
    static doublereal rmin, rmax, tnrm, sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical wantz;
    static integer iscale;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern doublereal slanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int ssterf_(integer *, doublereal *, doublereal *,
	     integer *);
    static doublereal smlnum;
    extern /* Subroutine */ int ssteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);


/*  -- LAPACK driver routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 159 "sstev.f"
    /* Parameter adjustments */
#line 159 "sstev.f"
    --d__;
#line 159 "sstev.f"
    --e;
#line 159 "sstev.f"
    z_dim1 = *ldz;
#line 159 "sstev.f"
    z_offset = 1 + z_dim1;
#line 159 "sstev.f"
    z__ -= z_offset;
#line 159 "sstev.f"
    --work;
#line 159 "sstev.f"

#line 159 "sstev.f"
    /* Function Body */
#line 159 "sstev.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);

#line 161 "sstev.f"
    *info = 0;
#line 162 "sstev.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 163 "sstev.f"
	*info = -1;
#line 164 "sstev.f"
    } else if (*n < 0) {
#line 165 "sstev.f"
	*info = -2;
#line 166 "sstev.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 167 "sstev.f"
	*info = -6;
#line 168 "sstev.f"
    }

#line 170 "sstev.f"
    if (*info != 0) {
#line 171 "sstev.f"
	i__1 = -(*info);
#line 171 "sstev.f"
	xerbla_("SSTEV ", &i__1, (ftnlen)6);
#line 172 "sstev.f"
	return 0;
#line 173 "sstev.f"
    }

/*     Quick return if possible */

#line 177 "sstev.f"
    if (*n == 0) {
#line 177 "sstev.f"
	return 0;
#line 177 "sstev.f"
    }

#line 180 "sstev.f"
    if (*n == 1) {
#line 181 "sstev.f"
	if (wantz) {
#line 181 "sstev.f"
	    z__[z_dim1 + 1] = 1.;
#line 181 "sstev.f"
	}
#line 183 "sstev.f"
	return 0;
#line 184 "sstev.f"
    }

/*     Get machine constants. */

#line 188 "sstev.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 189 "sstev.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 190 "sstev.f"
    smlnum = safmin / eps;
#line 191 "sstev.f"
    bignum = 1. / smlnum;
#line 192 "sstev.f"
    rmin = sqrt(smlnum);
#line 193 "sstev.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 197 "sstev.f"
    iscale = 0;
#line 198 "sstev.f"
    tnrm = slanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 199 "sstev.f"
    if (tnrm > 0. && tnrm < rmin) {
#line 200 "sstev.f"
	iscale = 1;
#line 201 "sstev.f"
	sigma = rmin / tnrm;
#line 202 "sstev.f"
    } else if (tnrm > rmax) {
#line 203 "sstev.f"
	iscale = 1;
#line 204 "sstev.f"
	sigma = rmax / tnrm;
#line 205 "sstev.f"
    }
#line 206 "sstev.f"
    if (iscale == 1) {
#line 207 "sstev.f"
	sscal_(n, &sigma, &d__[1], &c__1);
#line 208 "sstev.f"
	i__1 = *n - 1;
#line 208 "sstev.f"
	sscal_(&i__1, &sigma, &e[1], &c__1);
#line 209 "sstev.f"
    }

/*     For eigenvalues only, call SSTERF.  For eigenvalues and */
/*     eigenvectors, call SSTEQR. */

#line 214 "sstev.f"
    if (! wantz) {
#line 215 "sstev.f"
	ssterf_(n, &d__[1], &e[1], info);
#line 216 "sstev.f"
    } else {
#line 217 "sstev.f"
	ssteqr_("I", n, &d__[1], &e[1], &z__[z_offset], ldz, &work[1], info, (
		ftnlen)1);
#line 218 "sstev.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 222 "sstev.f"
    if (iscale == 1) {
#line 223 "sstev.f"
	if (*info == 0) {
#line 224 "sstev.f"
	    imax = *n;
#line 225 "sstev.f"
	} else {
#line 226 "sstev.f"
	    imax = *info - 1;
#line 227 "sstev.f"
	}
#line 228 "sstev.f"
	d__1 = 1. / sigma;
#line 228 "sstev.f"
	sscal_(&imax, &d__1, &d__[1], &c__1);
#line 229 "sstev.f"
    }

#line 231 "sstev.f"
    return 0;

/*     End of SSTEV */

} /* sstev_ */

