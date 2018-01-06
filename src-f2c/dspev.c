#line 1 "dspev.f"
/* dspev.f -- translated by f2c (version 20100827).
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

#line 1 "dspev.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> DSPEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER m
atrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSPEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dspev.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dspev.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dspev.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AP( * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSPEV computes all the eigenvalues and, optionally, eigenvectors of a */
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
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
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
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, N) */
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
/* >          WORK is DOUBLE PRECISION array, dimension (3*N) */
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

/* > \ingroup doubleOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int dspev_(char *jobz, char *uplo, integer *n, doublereal *
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
    static doublereal rmin, rmax;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static logical wantz;
    extern doublereal dlamch_(char *, ftnlen);
    static integer iscale;
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern doublereal dlansp_(char *, char *, integer *, doublereal *, 
	    doublereal *, ftnlen, ftnlen);
    static integer indtau;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *);
    static integer indwrk;
    extern /* Subroutine */ int dopgtr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), dsptrd_(char *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen), dsteqr_(char *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen);
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

#line 173 "dspev.f"
    /* Parameter adjustments */
#line 173 "dspev.f"
    --ap;
#line 173 "dspev.f"
    --w;
#line 173 "dspev.f"
    z_dim1 = *ldz;
#line 173 "dspev.f"
    z_offset = 1 + z_dim1;
#line 173 "dspev.f"
    z__ -= z_offset;
#line 173 "dspev.f"
    --work;
#line 173 "dspev.f"

#line 173 "dspev.f"
    /* Function Body */
#line 173 "dspev.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);

#line 175 "dspev.f"
    *info = 0;
#line 176 "dspev.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 177 "dspev.f"
	*info = -1;
#line 178 "dspev.f"
    } else if (! (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) || lsame_(uplo, 
	    "L", (ftnlen)1, (ftnlen)1))) {
#line 180 "dspev.f"
	*info = -2;
#line 181 "dspev.f"
    } else if (*n < 0) {
#line 182 "dspev.f"
	*info = -3;
#line 183 "dspev.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 184 "dspev.f"
	*info = -7;
#line 185 "dspev.f"
    }

#line 187 "dspev.f"
    if (*info != 0) {
#line 188 "dspev.f"
	i__1 = -(*info);
#line 188 "dspev.f"
	xerbla_("DSPEV ", &i__1, (ftnlen)6);
#line 189 "dspev.f"
	return 0;
#line 190 "dspev.f"
    }

/*     Quick return if possible */

#line 194 "dspev.f"
    if (*n == 0) {
#line 194 "dspev.f"
	return 0;
#line 194 "dspev.f"
    }

#line 197 "dspev.f"
    if (*n == 1) {
#line 198 "dspev.f"
	w[1] = ap[1];
#line 199 "dspev.f"
	if (wantz) {
#line 199 "dspev.f"
	    z__[z_dim1 + 1] = 1.;
#line 199 "dspev.f"
	}
#line 201 "dspev.f"
	return 0;
#line 202 "dspev.f"
    }

/*     Get machine constants. */

#line 206 "dspev.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 207 "dspev.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 208 "dspev.f"
    smlnum = safmin / eps;
#line 209 "dspev.f"
    bignum = 1. / smlnum;
#line 210 "dspev.f"
    rmin = sqrt(smlnum);
#line 211 "dspev.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 215 "dspev.f"
    anrm = dlansp_("M", uplo, n, &ap[1], &work[1], (ftnlen)1, (ftnlen)1);
#line 216 "dspev.f"
    iscale = 0;
#line 217 "dspev.f"
    if (anrm > 0. && anrm < rmin) {
#line 218 "dspev.f"
	iscale = 1;
#line 219 "dspev.f"
	sigma = rmin / anrm;
#line 220 "dspev.f"
    } else if (anrm > rmax) {
#line 221 "dspev.f"
	iscale = 1;
#line 222 "dspev.f"
	sigma = rmax / anrm;
#line 223 "dspev.f"
    }
#line 224 "dspev.f"
    if (iscale == 1) {
#line 225 "dspev.f"
	i__1 = *n * (*n + 1) / 2;
#line 225 "dspev.f"
	dscal_(&i__1, &sigma, &ap[1], &c__1);
#line 226 "dspev.f"
    }

/*     Call DSPTRD to reduce symmetric packed matrix to tridiagonal form. */

#line 230 "dspev.f"
    inde = 1;
#line 231 "dspev.f"
    indtau = inde + *n;
#line 232 "dspev.f"
    dsptrd_(uplo, n, &ap[1], &w[1], &work[inde], &work[indtau], &iinfo, (
	    ftnlen)1);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, first call */
/*     DOPGTR to generate the orthogonal matrix, then call DSTEQR. */

#line 237 "dspev.f"
    if (! wantz) {
#line 238 "dspev.f"
	dsterf_(n, &w[1], &work[inde], info);
#line 239 "dspev.f"
    } else {
#line 240 "dspev.f"
	indwrk = indtau + *n;
#line 241 "dspev.f"
	dopgtr_(uplo, n, &ap[1], &work[indtau], &z__[z_offset], ldz, &work[
		indwrk], &iinfo, (ftnlen)1);
#line 243 "dspev.f"
	dsteqr_(jobz, n, &w[1], &work[inde], &z__[z_offset], ldz, &work[
		indtau], info, (ftnlen)1);
#line 245 "dspev.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 249 "dspev.f"
    if (iscale == 1) {
#line 250 "dspev.f"
	if (*info == 0) {
#line 251 "dspev.f"
	    imax = *n;
#line 252 "dspev.f"
	} else {
#line 253 "dspev.f"
	    imax = *info - 1;
#line 254 "dspev.f"
	}
#line 255 "dspev.f"
	d__1 = 1. / sigma;
#line 255 "dspev.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 256 "dspev.f"
    }

#line 258 "dspev.f"
    return 0;

/*     End of DSPEV */

} /* dspev_ */

