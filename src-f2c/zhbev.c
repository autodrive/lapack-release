#line 1 "zhbev.f"
/* zhbev.f -- translated by f2c (version 20100827).
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

#line 1 "zhbev.f"
/* Table of constant values */

static doublereal c_b11 = 1.;
static integer c__1 = 1;

/* > \brief <b> ZHBEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER m
atrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHBEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhbev.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhbev.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhbev.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHBEV( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, */
/*                         RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, KD, LDAB, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   RWORK( * ), W( * ) */
/*       COMPLEX*16         AB( LDAB, * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHBEV computes all the eigenvalues and, optionally, eigenvectors of */
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
/* >          AB is COMPLEX*16 array, dimension (LDAB, N) */
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
/* >          WORK is COMPLEX*16 array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (max(1,3*N-2)) */
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

/* > \ingroup complex16OTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int zhbev_(char *jobz, char *uplo, integer *n, integer *kd, 
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
    static doublereal rmin, rmax;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static logical lower, wantz;
    extern doublereal dlamch_(char *, ftnlen);
    static integer iscale;
    static doublereal safmin;
    extern doublereal zlanhb_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *), zlascl_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublecomplex *, integer *, 
	    integer *, ftnlen), zhbtrd_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen);
    static integer indrwk;
    static doublereal smlnum;
    extern /* Subroutine */ int zsteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *,
	     ftnlen);


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

#line 196 "zhbev.f"
    /* Parameter adjustments */
#line 196 "zhbev.f"
    ab_dim1 = *ldab;
#line 196 "zhbev.f"
    ab_offset = 1 + ab_dim1;
#line 196 "zhbev.f"
    ab -= ab_offset;
#line 196 "zhbev.f"
    --w;
#line 196 "zhbev.f"
    z_dim1 = *ldz;
#line 196 "zhbev.f"
    z_offset = 1 + z_dim1;
#line 196 "zhbev.f"
    z__ -= z_offset;
#line 196 "zhbev.f"
    --work;
#line 196 "zhbev.f"
    --rwork;
#line 196 "zhbev.f"

#line 196 "zhbev.f"
    /* Function Body */
#line 196 "zhbev.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 197 "zhbev.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);

#line 199 "zhbev.f"
    *info = 0;
#line 200 "zhbev.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 201 "zhbev.f"
	*info = -1;
#line 202 "zhbev.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 203 "zhbev.f"
	*info = -2;
#line 204 "zhbev.f"
    } else if (*n < 0) {
#line 205 "zhbev.f"
	*info = -3;
#line 206 "zhbev.f"
    } else if (*kd < 0) {
#line 207 "zhbev.f"
	*info = -4;
#line 208 "zhbev.f"
    } else if (*ldab < *kd + 1) {
#line 209 "zhbev.f"
	*info = -6;
#line 210 "zhbev.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 211 "zhbev.f"
	*info = -9;
#line 212 "zhbev.f"
    }

#line 214 "zhbev.f"
    if (*info != 0) {
#line 215 "zhbev.f"
	i__1 = -(*info);
#line 215 "zhbev.f"
	xerbla_("ZHBEV ", &i__1, (ftnlen)6);
#line 216 "zhbev.f"
	return 0;
#line 217 "zhbev.f"
    }

/*     Quick return if possible */

#line 221 "zhbev.f"
    if (*n == 0) {
#line 221 "zhbev.f"
	return 0;
#line 221 "zhbev.f"
    }

#line 224 "zhbev.f"
    if (*n == 1) {
#line 225 "zhbev.f"
	if (lower) {
#line 226 "zhbev.f"
	    i__1 = ab_dim1 + 1;
#line 226 "zhbev.f"
	    w[1] = ab[i__1].r;
#line 227 "zhbev.f"
	} else {
#line 228 "zhbev.f"
	    i__1 = *kd + 1 + ab_dim1;
#line 228 "zhbev.f"
	    w[1] = ab[i__1].r;
#line 229 "zhbev.f"
	}
#line 230 "zhbev.f"
	if (wantz) {
#line 230 "zhbev.f"
	    i__1 = z_dim1 + 1;
#line 230 "zhbev.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 230 "zhbev.f"
	}
#line 232 "zhbev.f"
	return 0;
#line 233 "zhbev.f"
    }

/*     Get machine constants. */

#line 237 "zhbev.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 238 "zhbev.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 239 "zhbev.f"
    smlnum = safmin / eps;
#line 240 "zhbev.f"
    bignum = 1. / smlnum;
#line 241 "zhbev.f"
    rmin = sqrt(smlnum);
#line 242 "zhbev.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 246 "zhbev.f"
    anrm = zlanhb_("M", uplo, n, kd, &ab[ab_offset], ldab, &rwork[1], (ftnlen)
	    1, (ftnlen)1);
#line 247 "zhbev.f"
    iscale = 0;
#line 248 "zhbev.f"
    if (anrm > 0. && anrm < rmin) {
#line 249 "zhbev.f"
	iscale = 1;
#line 250 "zhbev.f"
	sigma = rmin / anrm;
#line 251 "zhbev.f"
    } else if (anrm > rmax) {
#line 252 "zhbev.f"
	iscale = 1;
#line 253 "zhbev.f"
	sigma = rmax / anrm;
#line 254 "zhbev.f"
    }
#line 255 "zhbev.f"
    if (iscale == 1) {
#line 256 "zhbev.f"
	if (lower) {
#line 257 "zhbev.f"
	    zlascl_("B", kd, kd, &c_b11, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 258 "zhbev.f"
	} else {
#line 259 "zhbev.f"
	    zlascl_("Q", kd, kd, &c_b11, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 260 "zhbev.f"
	}
#line 261 "zhbev.f"
    }

/*     Call ZHBTRD to reduce Hermitian band matrix to tridiagonal form. */

#line 265 "zhbev.f"
    inde = 1;
#line 266 "zhbev.f"
    zhbtrd_(jobz, uplo, n, kd, &ab[ab_offset], ldab, &w[1], &rwork[inde], &
	    z__[z_offset], ldz, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, call ZSTEQR. */

#line 271 "zhbev.f"
    if (! wantz) {
#line 272 "zhbev.f"
	dsterf_(n, &w[1], &rwork[inde], info);
#line 273 "zhbev.f"
    } else {
#line 274 "zhbev.f"
	indrwk = inde + *n;
#line 275 "zhbev.f"
	zsteqr_(jobz, n, &w[1], &rwork[inde], &z__[z_offset], ldz, &rwork[
		indrwk], info, (ftnlen)1);
#line 277 "zhbev.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 281 "zhbev.f"
    if (iscale == 1) {
#line 282 "zhbev.f"
	if (*info == 0) {
#line 283 "zhbev.f"
	    imax = *n;
#line 284 "zhbev.f"
	} else {
#line 285 "zhbev.f"
	    imax = *info - 1;
#line 286 "zhbev.f"
	}
#line 287 "zhbev.f"
	d__1 = 1. / sigma;
#line 287 "zhbev.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 288 "zhbev.f"
    }

#line 290 "zhbev.f"
    return 0;

/*     End of ZHBEV */

} /* zhbev_ */

