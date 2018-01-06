#line 1 "chbev.f"
/* chbev.f -- translated by f2c (version 20100827).
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

#line 1 "chbev.f"
/* Table of constant values */

static doublereal c_b11 = 1.;
static integer c__1 = 1;

/* > \brief <b> CHBEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER m
atrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHBEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chbev.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chbev.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chbev.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHBEV( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, */
/*                         RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, KD, LDAB, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            AB( LDAB, * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHBEV computes all the eigenvalues and, optionally, eigenvectors of */
/* > a complex Hermitian band matrix A. */
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
/* > \param[in] KD */
/* > \verbatim */
/* >          KD is INTEGER */
/* >          The number of superdiagonals of the matrix A if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'.  KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is COMPLEX array, dimension (LDAB, N) */
/* >          On entry, the upper or lower triangle of the Hermitian band */
/* >          matrix A, stored in the first KD+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* > */
/* >          On exit, AB is overwritten by values generated during the */
/* >          reduction to tridiagonal form.  If UPLO = 'U', the first */
/* >          superdiagonal and the diagonal of the tridiagonal matrix T */
/* >          are returned in rows KD and KD+1 of AB, and if UPLO = 'L', */
/* >          the diagonal and first subdiagonal of T are returned in the */
/* >          first two rows of AB. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KD + 1. */
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
/* >          WORK is COMPLEX array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (max(1,3*N-2)) */
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
/* Subroutine */ int chbev_(char *jobz, char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *w, doublecomplex *z__, 
	integer *ldz, doublecomplex *work, doublereal *rwork, integer *info, 
	ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, z_dim1, z_offset, i__1;
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
    static logical lower, wantz;
    extern doublereal clanhb_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, ftnlen, ftnlen);
    static integer iscale;
    extern /* Subroutine */ int clascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen), chbtrd_(char *, char *, integer *,
	     integer *, doublecomplex *, integer *, doublereal *, doublereal *
	    , doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static integer indrwk;
    extern /* Subroutine */ int csteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *,
	     ftnlen), ssterf_(integer *, doublereal *, doublereal *, integer *
	    );
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

#line 196 "chbev.f"
    /* Parameter adjustments */
#line 196 "chbev.f"
    ab_dim1 = *ldab;
#line 196 "chbev.f"
    ab_offset = 1 + ab_dim1;
#line 196 "chbev.f"
    ab -= ab_offset;
#line 196 "chbev.f"
    --w;
#line 196 "chbev.f"
    z_dim1 = *ldz;
#line 196 "chbev.f"
    z_offset = 1 + z_dim1;
#line 196 "chbev.f"
    z__ -= z_offset;
#line 196 "chbev.f"
    --work;
#line 196 "chbev.f"
    --rwork;
#line 196 "chbev.f"

#line 196 "chbev.f"
    /* Function Body */
#line 196 "chbev.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 197 "chbev.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);

#line 199 "chbev.f"
    *info = 0;
#line 200 "chbev.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 201 "chbev.f"
	*info = -1;
#line 202 "chbev.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 203 "chbev.f"
	*info = -2;
#line 204 "chbev.f"
    } else if (*n < 0) {
#line 205 "chbev.f"
	*info = -3;
#line 206 "chbev.f"
    } else if (*kd < 0) {
#line 207 "chbev.f"
	*info = -4;
#line 208 "chbev.f"
    } else if (*ldab < *kd + 1) {
#line 209 "chbev.f"
	*info = -6;
#line 210 "chbev.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 211 "chbev.f"
	*info = -9;
#line 212 "chbev.f"
    }

#line 214 "chbev.f"
    if (*info != 0) {
#line 215 "chbev.f"
	i__1 = -(*info);
#line 215 "chbev.f"
	xerbla_("CHBEV ", &i__1, (ftnlen)6);
#line 216 "chbev.f"
	return 0;
#line 217 "chbev.f"
    }

/*     Quick return if possible */

#line 221 "chbev.f"
    if (*n == 0) {
#line 221 "chbev.f"
	return 0;
#line 221 "chbev.f"
    }

#line 224 "chbev.f"
    if (*n == 1) {
#line 225 "chbev.f"
	if (lower) {
#line 226 "chbev.f"
	    i__1 = ab_dim1 + 1;
#line 226 "chbev.f"
	    w[1] = ab[i__1].r;
#line 227 "chbev.f"
	} else {
#line 228 "chbev.f"
	    i__1 = *kd + 1 + ab_dim1;
#line 228 "chbev.f"
	    w[1] = ab[i__1].r;
#line 229 "chbev.f"
	}
#line 230 "chbev.f"
	if (wantz) {
#line 230 "chbev.f"
	    i__1 = z_dim1 + 1;
#line 230 "chbev.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 230 "chbev.f"
	}
#line 232 "chbev.f"
	return 0;
#line 233 "chbev.f"
    }

/*     Get machine constants. */

#line 237 "chbev.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 238 "chbev.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 239 "chbev.f"
    smlnum = safmin / eps;
#line 240 "chbev.f"
    bignum = 1. / smlnum;
#line 241 "chbev.f"
    rmin = sqrt(smlnum);
#line 242 "chbev.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 246 "chbev.f"
    anrm = clanhb_("M", uplo, n, kd, &ab[ab_offset], ldab, &rwork[1], (ftnlen)
	    1, (ftnlen)1);
#line 247 "chbev.f"
    iscale = 0;
#line 248 "chbev.f"
    if (anrm > 0. && anrm < rmin) {
#line 249 "chbev.f"
	iscale = 1;
#line 250 "chbev.f"
	sigma = rmin / anrm;
#line 251 "chbev.f"
    } else if (anrm > rmax) {
#line 252 "chbev.f"
	iscale = 1;
#line 253 "chbev.f"
	sigma = rmax / anrm;
#line 254 "chbev.f"
    }
#line 255 "chbev.f"
    if (iscale == 1) {
#line 256 "chbev.f"
	if (lower) {
#line 257 "chbev.f"
	    clascl_("B", kd, kd, &c_b11, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 258 "chbev.f"
	} else {
#line 259 "chbev.f"
	    clascl_("Q", kd, kd, &c_b11, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 260 "chbev.f"
	}
#line 261 "chbev.f"
    }

/*     Call CHBTRD to reduce Hermitian band matrix to tridiagonal form. */

#line 265 "chbev.f"
    inde = 1;
#line 266 "chbev.f"
    chbtrd_(jobz, uplo, n, kd, &ab[ab_offset], ldab, &w[1], &rwork[inde], &
	    z__[z_offset], ldz, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);

/*     For eigenvalues only, call SSTERF.  For eigenvectors, call CSTEQR. */

#line 271 "chbev.f"
    if (! wantz) {
#line 272 "chbev.f"
	ssterf_(n, &w[1], &rwork[inde], info);
#line 273 "chbev.f"
    } else {
#line 274 "chbev.f"
	indrwk = inde + *n;
#line 275 "chbev.f"
	csteqr_(jobz, n, &w[1], &rwork[inde], &z__[z_offset], ldz, &rwork[
		indrwk], info, (ftnlen)1);
#line 277 "chbev.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 281 "chbev.f"
    if (iscale == 1) {
#line 282 "chbev.f"
	if (*info == 0) {
#line 283 "chbev.f"
	    imax = *n;
#line 284 "chbev.f"
	} else {
#line 285 "chbev.f"
	    imax = *info - 1;
#line 286 "chbev.f"
	}
#line 287 "chbev.f"
	d__1 = 1. / sigma;
#line 287 "chbev.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 288 "chbev.f"
    }

#line 290 "chbev.f"
    return 0;

/*     End of CHBEV */

} /* chbev_ */

