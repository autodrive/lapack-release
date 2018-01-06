#line 1 "sspev.f"
/* sspev.f -- translated by f2c (version 20100827).
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

#line 1 "sspev.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> SSPEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER m
atrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSPEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sspev.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sspev.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sspev.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AP( * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSPEV computes all the eigenvalues and, optionally, eigenvectors of a */
/* > real symmetric matrix A in packed storage. */
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
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is REAL array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the symmetric matrix */
/* >          A, packed columnwise in a linear array.  The j-th column of A */
/* >          is stored in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
/* > */
/* >          On exit, AP is overwritten by values generated during the */
/* >          reduction to tridiagonal form.  If UPLO = 'U', the diagonal */
/* >          and first superdiagonal of the tridiagonal matrix T overwrite */
/* >          the corresponding elements of A, and if UPLO = 'L', the */
/* >          diagonal and first subdiagonal of T overwrite the */
/* >          corresponding elements of A. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is REAL array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (LDZ, N) */
/* >          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal */
/* >          eigenvectors of the matrix A, with the i-th column of Z */
/* >          holding the eigenvector associated with W(i). */
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
/* >          WORK is REAL array, dimension (3*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if INFO = i, the algorithm failed to converge; i */
/* >                off-diagonal elements of an intermediate tridiagonal */
/* >                form did not converge to zero. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int sspev_(char *jobz, char *uplo, integer *n, doublereal *
	ap, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	integer *info, ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal eps;
    static integer inde;
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax, sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical wantz;
    static integer iscale;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static integer indtau, indwrk;
    extern doublereal slansp_(char *, char *, integer *, doublereal *, 
	    doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int ssterf_(integer *, doublereal *, doublereal *,
	     integer *);
    static doublereal smlnum;
    extern /* Subroutine */ int sopgtr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), ssptrd_(char *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen), ssteqr_(char *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen);


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

#line 173 "sspev.f"
    /* Parameter adjustments */
#line 173 "sspev.f"
    --ap;
#line 173 "sspev.f"
    --w;
#line 173 "sspev.f"
    z_dim1 = *ldz;
#line 173 "sspev.f"
    z_offset = 1 + z_dim1;
#line 173 "sspev.f"
    z__ -= z_offset;
#line 173 "sspev.f"
    --work;
#line 173 "sspev.f"

#line 173 "sspev.f"
    /* Function Body */
#line 173 "sspev.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);

#line 175 "sspev.f"
    *info = 0;
#line 176 "sspev.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 177 "sspev.f"
	*info = -1;
#line 178 "sspev.f"
    } else if (! (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) || lsame_(uplo, 
	    "L", (ftnlen)1, (ftnlen)1))) {
#line 180 "sspev.f"
	*info = -2;
#line 181 "sspev.f"
    } else if (*n < 0) {
#line 182 "sspev.f"
	*info = -3;
#line 183 "sspev.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 184 "sspev.f"
	*info = -7;
#line 185 "sspev.f"
    }

#line 187 "sspev.f"
    if (*info != 0) {
#line 188 "sspev.f"
	i__1 = -(*info);
#line 188 "sspev.f"
	xerbla_("SSPEV ", &i__1, (ftnlen)6);
#line 189 "sspev.f"
	return 0;
#line 190 "sspev.f"
    }

/*     Quick return if possible */

#line 194 "sspev.f"
    if (*n == 0) {
#line 194 "sspev.f"
	return 0;
#line 194 "sspev.f"
    }

#line 197 "sspev.f"
    if (*n == 1) {
#line 198 "sspev.f"
	w[1] = ap[1];
#line 199 "sspev.f"
	if (wantz) {
#line 199 "sspev.f"
	    z__[z_dim1 + 1] = 1.;
#line 199 "sspev.f"
	}
#line 201 "sspev.f"
	return 0;
#line 202 "sspev.f"
    }

/*     Get machine constants. */

#line 206 "sspev.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 207 "sspev.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 208 "sspev.f"
    smlnum = safmin / eps;
#line 209 "sspev.f"
    bignum = 1. / smlnum;
#line 210 "sspev.f"
    rmin = sqrt(smlnum);
#line 211 "sspev.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 215 "sspev.f"
    anrm = slansp_("M", uplo, n, &ap[1], &work[1], (ftnlen)1, (ftnlen)1);
#line 216 "sspev.f"
    iscale = 0;
#line 217 "sspev.f"
    if (anrm > 0. && anrm < rmin) {
#line 218 "sspev.f"
	iscale = 1;
#line 219 "sspev.f"
	sigma = rmin / anrm;
#line 220 "sspev.f"
    } else if (anrm > rmax) {
#line 221 "sspev.f"
	iscale = 1;
#line 222 "sspev.f"
	sigma = rmax / anrm;
#line 223 "sspev.f"
    }
#line 224 "sspev.f"
    if (iscale == 1) {
#line 225 "sspev.f"
	i__1 = *n * (*n + 1) / 2;
#line 225 "sspev.f"
	sscal_(&i__1, &sigma, &ap[1], &c__1);
#line 226 "sspev.f"
    }

/*     Call SSPTRD to reduce symmetric packed matrix to tridiagonal form. */

#line 230 "sspev.f"
    inde = 1;
#line 231 "sspev.f"
    indtau = inde + *n;
#line 232 "sspev.f"
    ssptrd_(uplo, n, &ap[1], &w[1], &work[inde], &work[indtau], &iinfo, (
	    ftnlen)1);

/*     For eigenvalues only, call SSTERF.  For eigenvectors, first call */
/*     SOPGTR to generate the orthogonal matrix, then call SSTEQR. */

#line 237 "sspev.f"
    if (! wantz) {
#line 238 "sspev.f"
	ssterf_(n, &w[1], &work[inde], info);
#line 239 "sspev.f"
    } else {
#line 240 "sspev.f"
	indwrk = indtau + *n;
#line 241 "sspev.f"
	sopgtr_(uplo, n, &ap[1], &work[indtau], &z__[z_offset], ldz, &work[
		indwrk], &iinfo, (ftnlen)1);
#line 243 "sspev.f"
	ssteqr_(jobz, n, &w[1], &work[inde], &z__[z_offset], ldz, &work[
		indtau], info, (ftnlen)1);
#line 245 "sspev.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 249 "sspev.f"
    if (iscale == 1) {
#line 250 "sspev.f"
	if (*info == 0) {
#line 251 "sspev.f"
	    imax = *n;
#line 252 "sspev.f"
	} else {
#line 253 "sspev.f"
	    imax = *info - 1;
#line 254 "sspev.f"
	}
#line 255 "sspev.f"
	d__1 = 1. / sigma;
#line 255 "sspev.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 256 "sspev.f"
    }

#line 258 "sspev.f"
    return 0;

/*     End of SSPEV */

} /* sspev_ */

