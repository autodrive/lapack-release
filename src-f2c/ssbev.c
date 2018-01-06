#line 1 "ssbev.f"
/* ssbev.f -- translated by f2c (version 20100827).
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

#line 1 "ssbev.f"
/* Table of constant values */

static doublereal c_b11 = 1.;
static integer c__1 = 1;

/* > \brief <b> SSBEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER m
atrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSBEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssbev.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssbev.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssbev.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSBEV( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, */
/*                         INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, KD, LDAB, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSBEV computes all the eigenvalues and, optionally, eigenvectors of */
/* > a real symmetric band matrix A. */
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
/* >          AB is REAL array, dimension (LDAB, N) */
/* >          On entry, the upper or lower triangle of the symmetric band */
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
/* >          WORK is REAL array, dimension (max(1,3*N-2)) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
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

/* > \ingroup realOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int ssbev_(char *jobz, char *uplo, integer *n, integer *kd, 
	doublereal *ab, integer *ldab, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *info, ftnlen jobz_len, 
	ftnlen uplo_len)
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
    static integer iscale;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern doublereal slansb_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer indwrk;
    extern /* Subroutine */ int ssbtrd_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, ftnlen, ftnlen), ssterf_(
	    integer *, doublereal *, doublereal *, integer *);
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

#line 189 "ssbev.f"
    /* Parameter adjustments */
#line 189 "ssbev.f"
    ab_dim1 = *ldab;
#line 189 "ssbev.f"
    ab_offset = 1 + ab_dim1;
#line 189 "ssbev.f"
    ab -= ab_offset;
#line 189 "ssbev.f"
    --w;
#line 189 "ssbev.f"
    z_dim1 = *ldz;
#line 189 "ssbev.f"
    z_offset = 1 + z_dim1;
#line 189 "ssbev.f"
    z__ -= z_offset;
#line 189 "ssbev.f"
    --work;
#line 189 "ssbev.f"

#line 189 "ssbev.f"
    /* Function Body */
#line 189 "ssbev.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 190 "ssbev.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);

#line 192 "ssbev.f"
    *info = 0;
#line 193 "ssbev.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 194 "ssbev.f"
	*info = -1;
#line 195 "ssbev.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 196 "ssbev.f"
	*info = -2;
#line 197 "ssbev.f"
    } else if (*n < 0) {
#line 198 "ssbev.f"
	*info = -3;
#line 199 "ssbev.f"
    } else if (*kd < 0) {
#line 200 "ssbev.f"
	*info = -4;
#line 201 "ssbev.f"
    } else if (*ldab < *kd + 1) {
#line 202 "ssbev.f"
	*info = -6;
#line 203 "ssbev.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 204 "ssbev.f"
	*info = -9;
#line 205 "ssbev.f"
    }

#line 207 "ssbev.f"
    if (*info != 0) {
#line 208 "ssbev.f"
	i__1 = -(*info);
#line 208 "ssbev.f"
	xerbla_("SSBEV ", &i__1, (ftnlen)6);
#line 209 "ssbev.f"
	return 0;
#line 210 "ssbev.f"
    }

/*     Quick return if possible */

#line 214 "ssbev.f"
    if (*n == 0) {
#line 214 "ssbev.f"
	return 0;
#line 214 "ssbev.f"
    }

#line 217 "ssbev.f"
    if (*n == 1) {
#line 218 "ssbev.f"
	if (lower) {
#line 219 "ssbev.f"
	    w[1] = ab[ab_dim1 + 1];
#line 220 "ssbev.f"
	} else {
#line 221 "ssbev.f"
	    w[1] = ab[*kd + 1 + ab_dim1];
#line 222 "ssbev.f"
	}
#line 223 "ssbev.f"
	if (wantz) {
#line 223 "ssbev.f"
	    z__[z_dim1 + 1] = 1.;
#line 223 "ssbev.f"
	}
#line 225 "ssbev.f"
	return 0;
#line 226 "ssbev.f"
    }

/*     Get machine constants. */

#line 230 "ssbev.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 231 "ssbev.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 232 "ssbev.f"
    smlnum = safmin / eps;
#line 233 "ssbev.f"
    bignum = 1. / smlnum;
#line 234 "ssbev.f"
    rmin = sqrt(smlnum);
#line 235 "ssbev.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 239 "ssbev.f"
    anrm = slansb_("M", uplo, n, kd, &ab[ab_offset], ldab, &work[1], (ftnlen)
	    1, (ftnlen)1);
#line 240 "ssbev.f"
    iscale = 0;
#line 241 "ssbev.f"
    if (anrm > 0. && anrm < rmin) {
#line 242 "ssbev.f"
	iscale = 1;
#line 243 "ssbev.f"
	sigma = rmin / anrm;
#line 244 "ssbev.f"
    } else if (anrm > rmax) {
#line 245 "ssbev.f"
	iscale = 1;
#line 246 "ssbev.f"
	sigma = rmax / anrm;
#line 247 "ssbev.f"
    }
#line 248 "ssbev.f"
    if (iscale == 1) {
#line 249 "ssbev.f"
	if (lower) {
#line 250 "ssbev.f"
	    slascl_("B", kd, kd, &c_b11, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 251 "ssbev.f"
	} else {
#line 252 "ssbev.f"
	    slascl_("Q", kd, kd, &c_b11, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 253 "ssbev.f"
	}
#line 254 "ssbev.f"
    }

/*     Call SSBTRD to reduce symmetric band matrix to tridiagonal form. */

#line 258 "ssbev.f"
    inde = 1;
#line 259 "ssbev.f"
    indwrk = inde + *n;
#line 260 "ssbev.f"
    ssbtrd_(jobz, uplo, n, kd, &ab[ab_offset], ldab, &w[1], &work[inde], &z__[
	    z_offset], ldz, &work[indwrk], &iinfo, (ftnlen)1, (ftnlen)1);

/*     For eigenvalues only, call SSTERF.  For eigenvectors, call SSTEQR. */

#line 265 "ssbev.f"
    if (! wantz) {
#line 266 "ssbev.f"
	ssterf_(n, &w[1], &work[inde], info);
#line 267 "ssbev.f"
    } else {
#line 268 "ssbev.f"
	ssteqr_(jobz, n, &w[1], &work[inde], &z__[z_offset], ldz, &work[
		indwrk], info, (ftnlen)1);
#line 270 "ssbev.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 274 "ssbev.f"
    if (iscale == 1) {
#line 275 "ssbev.f"
	if (*info == 0) {
#line 276 "ssbev.f"
	    imax = *n;
#line 277 "ssbev.f"
	} else {
#line 278 "ssbev.f"
	    imax = *info - 1;
#line 279 "ssbev.f"
	}
#line 280 "ssbev.f"
	d__1 = 1. / sigma;
#line 280 "ssbev.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 281 "ssbev.f"
    }

#line 283 "ssbev.f"
    return 0;

/*     End of SSBEV */

} /* ssbev_ */

