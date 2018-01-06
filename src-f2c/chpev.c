#line 1 "chpev.f"
/* chpev.f -- translated by f2c (version 20100827).
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

#line 1 "chpev.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> CHPEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER m
atrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHPEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chpev.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chpev.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chpev.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, RWORK, */
/*                         INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            AP( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHPEV computes all the eigenvalues and, optionally, eigenvectors of a */
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
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
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
/* >          W is REAL array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is COMPLEX array, dimension (LDZ, N) */
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
/* >          WORK is COMPLEX array, dimension (max(1, 2*N-1)) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (max(1, 3*N-2)) */
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

/* > \ingroup complexOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int chpev_(char *jobz, char *uplo, integer *n, doublecomplex 
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
    static doublereal rmin, rmax, sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical wantz;
    static integer iscale;
    extern doublereal clanhp_(char *, char *, integer *, doublecomplex *, 
	    doublereal *, ftnlen, ftnlen), slamch_(char *, ftnlen);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static integer indtau;
    extern /* Subroutine */ int chptrd_(char *, integer *, doublecomplex *, 
	    doublereal *, doublereal *, doublecomplex *, integer *, ftnlen);
    static integer indrwk, indwrk;
    extern /* Subroutine */ int csteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *,
	     ftnlen), cupgtr_(char *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, ftnlen), ssterf_(integer *, doublereal *, doublereal *,
	     integer *);
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

#line 184 "chpev.f"
    /* Parameter adjustments */
#line 184 "chpev.f"
    --ap;
#line 184 "chpev.f"
    --w;
#line 184 "chpev.f"
    z_dim1 = *ldz;
#line 184 "chpev.f"
    z_offset = 1 + z_dim1;
#line 184 "chpev.f"
    z__ -= z_offset;
#line 184 "chpev.f"
    --work;
#line 184 "chpev.f"
    --rwork;
#line 184 "chpev.f"

#line 184 "chpev.f"
    /* Function Body */
#line 184 "chpev.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);

#line 186 "chpev.f"
    *info = 0;
#line 187 "chpev.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 188 "chpev.f"
	*info = -1;
#line 189 "chpev.f"
    } else if (! (lsame_(uplo, "L", (ftnlen)1, (ftnlen)1) || lsame_(uplo, 
	    "U", (ftnlen)1, (ftnlen)1))) {
#line 191 "chpev.f"
	*info = -2;
#line 192 "chpev.f"
    } else if (*n < 0) {
#line 193 "chpev.f"
	*info = -3;
#line 194 "chpev.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 195 "chpev.f"
	*info = -7;
#line 196 "chpev.f"
    }

#line 198 "chpev.f"
    if (*info != 0) {
#line 199 "chpev.f"
	i__1 = -(*info);
#line 199 "chpev.f"
	xerbla_("CHPEV ", &i__1, (ftnlen)6);
#line 200 "chpev.f"
	return 0;
#line 201 "chpev.f"
    }

/*     Quick return if possible */

#line 205 "chpev.f"
    if (*n == 0) {
#line 205 "chpev.f"
	return 0;
#line 205 "chpev.f"
    }

#line 208 "chpev.f"
    if (*n == 1) {
#line 209 "chpev.f"
	w[1] = ap[1].r;
#line 210 "chpev.f"
	rwork[1] = 1.;
#line 211 "chpev.f"
	if (wantz) {
#line 211 "chpev.f"
	    i__1 = z_dim1 + 1;
#line 211 "chpev.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 211 "chpev.f"
	}
#line 213 "chpev.f"
	return 0;
#line 214 "chpev.f"
    }

/*     Get machine constants. */

#line 218 "chpev.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 219 "chpev.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 220 "chpev.f"
    smlnum = safmin / eps;
#line 221 "chpev.f"
    bignum = 1. / smlnum;
#line 222 "chpev.f"
    rmin = sqrt(smlnum);
#line 223 "chpev.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 227 "chpev.f"
    anrm = clanhp_("M", uplo, n, &ap[1], &rwork[1], (ftnlen)1, (ftnlen)1);
#line 228 "chpev.f"
    iscale = 0;
#line 229 "chpev.f"
    if (anrm > 0. && anrm < rmin) {
#line 230 "chpev.f"
	iscale = 1;
#line 231 "chpev.f"
	sigma = rmin / anrm;
#line 232 "chpev.f"
    } else if (anrm > rmax) {
#line 233 "chpev.f"
	iscale = 1;
#line 234 "chpev.f"
	sigma = rmax / anrm;
#line 235 "chpev.f"
    }
#line 236 "chpev.f"
    if (iscale == 1) {
#line 237 "chpev.f"
	i__1 = *n * (*n + 1) / 2;
#line 237 "chpev.f"
	csscal_(&i__1, &sigma, &ap[1], &c__1);
#line 238 "chpev.f"
    }

/*     Call CHPTRD to reduce Hermitian packed matrix to tridiagonal form. */

#line 242 "chpev.f"
    inde = 1;
#line 243 "chpev.f"
    indtau = 1;
#line 244 "chpev.f"
    chptrd_(uplo, n, &ap[1], &w[1], &rwork[inde], &work[indtau], &iinfo, (
	    ftnlen)1);

/*     For eigenvalues only, call SSTERF.  For eigenvectors, first call */
/*     CUPGTR to generate the orthogonal matrix, then call CSTEQR. */

#line 250 "chpev.f"
    if (! wantz) {
#line 251 "chpev.f"
	ssterf_(n, &w[1], &rwork[inde], info);
#line 252 "chpev.f"
    } else {
#line 253 "chpev.f"
	indwrk = indtau + *n;
#line 254 "chpev.f"
	cupgtr_(uplo, n, &ap[1], &work[indtau], &z__[z_offset], ldz, &work[
		indwrk], &iinfo, (ftnlen)1);
#line 256 "chpev.f"
	indrwk = inde + *n;
#line 257 "chpev.f"
	csteqr_(jobz, n, &w[1], &rwork[inde], &z__[z_offset], ldz, &rwork[
		indrwk], info, (ftnlen)1);
#line 259 "chpev.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 263 "chpev.f"
    if (iscale == 1) {
#line 264 "chpev.f"
	if (*info == 0) {
#line 265 "chpev.f"
	    imax = *n;
#line 266 "chpev.f"
	} else {
#line 267 "chpev.f"
	    imax = *info - 1;
#line 268 "chpev.f"
	}
#line 269 "chpev.f"
	d__1 = 1. / sigma;
#line 269 "chpev.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 270 "chpev.f"
    }

#line 272 "chpev.f"
    return 0;

/*     End of CHPEV */

} /* chpev_ */

