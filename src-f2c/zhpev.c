#line 1 "zhpev.f"
/* zhpev.f -- translated by f2c (version 20100827).
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

#line 1 "zhpev.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> ZHPEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER m
atrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHPEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpev.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpev.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpev.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, RWORK, */
/*                         INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   RWORK( * ), W( * ) */
/*       COMPLEX*16         AP( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHPEV computes all the eigenvalues and, optionally, eigenvectors of a */
/* > complex Hermitian matrix in packed storage. */
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
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the Hermitian matrix */
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
/* >          Z is COMPLEX*16 array, dimension (LDZ, N) */
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
/* >          WORK is COMPLEX*16 array, dimension (max(1, 2*N-1)) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (max(1, 3*N-2)) */
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

/* > \date November 2011 */

/* > \ingroup complex16OTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int zhpev_(char *jobz, char *uplo, integer *n, doublecomplex 
	*ap, doublereal *w, doublecomplex *z__, integer *ldz, doublecomplex *
	work, doublereal *rwork, integer *info, ftnlen jobz_len, ftnlen 
	uplo_len)
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
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    static doublereal bignum;
    static integer indtau;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *);
    extern doublereal zlanhp_(char *, char *, integer *, doublecomplex *, 
	    doublereal *, ftnlen, ftnlen);
    static integer indrwk, indwrk;
    static doublereal smlnum;
    extern /* Subroutine */ int zhptrd_(char *, integer *, doublecomplex *, 
	    doublereal *, doublereal *, doublecomplex *, integer *, ftnlen), 
	    zsteqr_(char *, integer *, doublereal *, doublereal *, 
	    doublecomplex *, integer *, doublereal *, integer *, ftnlen), 
	    zupgtr_(char *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);


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

#line 184 "zhpev.f"
    /* Parameter adjustments */
#line 184 "zhpev.f"
    --ap;
#line 184 "zhpev.f"
    --w;
#line 184 "zhpev.f"
    z_dim1 = *ldz;
#line 184 "zhpev.f"
    z_offset = 1 + z_dim1;
#line 184 "zhpev.f"
    z__ -= z_offset;
#line 184 "zhpev.f"
    --work;
#line 184 "zhpev.f"
    --rwork;
#line 184 "zhpev.f"

#line 184 "zhpev.f"
    /* Function Body */
#line 184 "zhpev.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);

#line 186 "zhpev.f"
    *info = 0;
#line 187 "zhpev.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 188 "zhpev.f"
	*info = -1;
#line 189 "zhpev.f"
    } else if (! (lsame_(uplo, "L", (ftnlen)1, (ftnlen)1) || lsame_(uplo, 
	    "U", (ftnlen)1, (ftnlen)1))) {
#line 191 "zhpev.f"
	*info = -2;
#line 192 "zhpev.f"
    } else if (*n < 0) {
#line 193 "zhpev.f"
	*info = -3;
#line 194 "zhpev.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 195 "zhpev.f"
	*info = -7;
#line 196 "zhpev.f"
    }

#line 198 "zhpev.f"
    if (*info != 0) {
#line 199 "zhpev.f"
	i__1 = -(*info);
#line 199 "zhpev.f"
	xerbla_("ZHPEV ", &i__1, (ftnlen)6);
#line 200 "zhpev.f"
	return 0;
#line 201 "zhpev.f"
    }

/*     Quick return if possible */

#line 205 "zhpev.f"
    if (*n == 0) {
#line 205 "zhpev.f"
	return 0;
#line 205 "zhpev.f"
    }

#line 208 "zhpev.f"
    if (*n == 1) {
#line 209 "zhpev.f"
	w[1] = ap[1].r;
#line 210 "zhpev.f"
	rwork[1] = 1.;
#line 211 "zhpev.f"
	if (wantz) {
#line 211 "zhpev.f"
	    i__1 = z_dim1 + 1;
#line 211 "zhpev.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 211 "zhpev.f"
	}
#line 213 "zhpev.f"
	return 0;
#line 214 "zhpev.f"
    }

/*     Get machine constants. */

#line 218 "zhpev.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 219 "zhpev.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 220 "zhpev.f"
    smlnum = safmin / eps;
#line 221 "zhpev.f"
    bignum = 1. / smlnum;
#line 222 "zhpev.f"
    rmin = sqrt(smlnum);
#line 223 "zhpev.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 227 "zhpev.f"
    anrm = zlanhp_("M", uplo, n, &ap[1], &rwork[1], (ftnlen)1, (ftnlen)1);
#line 228 "zhpev.f"
    iscale = 0;
#line 229 "zhpev.f"
    if (anrm > 0. && anrm < rmin) {
#line 230 "zhpev.f"
	iscale = 1;
#line 231 "zhpev.f"
	sigma = rmin / anrm;
#line 232 "zhpev.f"
    } else if (anrm > rmax) {
#line 233 "zhpev.f"
	iscale = 1;
#line 234 "zhpev.f"
	sigma = rmax / anrm;
#line 235 "zhpev.f"
    }
#line 236 "zhpev.f"
    if (iscale == 1) {
#line 237 "zhpev.f"
	i__1 = *n * (*n + 1) / 2;
#line 237 "zhpev.f"
	zdscal_(&i__1, &sigma, &ap[1], &c__1);
#line 238 "zhpev.f"
    }

/*     Call ZHPTRD to reduce Hermitian packed matrix to tridiagonal form. */

#line 242 "zhpev.f"
    inde = 1;
#line 243 "zhpev.f"
    indtau = 1;
#line 244 "zhpev.f"
    zhptrd_(uplo, n, &ap[1], &w[1], &rwork[inde], &work[indtau], &iinfo, (
	    ftnlen)1);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, first call */
/*     ZUPGTR to generate the orthogonal matrix, then call ZSTEQR. */

#line 250 "zhpev.f"
    if (! wantz) {
#line 251 "zhpev.f"
	dsterf_(n, &w[1], &rwork[inde], info);
#line 252 "zhpev.f"
    } else {
#line 253 "zhpev.f"
	indwrk = indtau + *n;
#line 254 "zhpev.f"
	zupgtr_(uplo, n, &ap[1], &work[indtau], &z__[z_offset], ldz, &work[
		indwrk], &iinfo, (ftnlen)1);
#line 256 "zhpev.f"
	indrwk = inde + *n;
#line 257 "zhpev.f"
	zsteqr_(jobz, n, &w[1], &rwork[inde], &z__[z_offset], ldz, &rwork[
		indrwk], info, (ftnlen)1);
#line 259 "zhpev.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 263 "zhpev.f"
    if (iscale == 1) {
#line 264 "zhpev.f"
	if (*info == 0) {
#line 265 "zhpev.f"
	    imax = *n;
#line 266 "zhpev.f"
	} else {
#line 267 "zhpev.f"
	    imax = *info - 1;
#line 268 "zhpev.f"
	}
#line 269 "zhpev.f"
	d__1 = 1. / sigma;
#line 269 "zhpev.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 270 "zhpev.f"
    }

#line 272 "zhpev.f"
    return 0;

/*     End of ZHPEV */

} /* zhpev_ */

