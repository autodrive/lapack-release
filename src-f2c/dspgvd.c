#line 1 "dspgvd.f"
/* dspgvd.f -- translated by f2c (version 20100827).
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

#line 1 "dspgvd.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DSPGVD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSPGVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dspgvd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dspgvd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dspgvd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSPGVD( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK, */
/*                          LWORK, IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, ITYPE, LDZ, LIWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   AP( * ), BP( * ), W( * ), WORK( * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSPGVD computes all the eigenvalues, and optionally, the eigenvectors */
/* > of a real generalized symmetric-definite eigenproblem, of the form */
/* > A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and */
/* > B are assumed to be symmetric, stored in packed format, and B is also */
/* > positive definite. */
/* > If eigenvectors are desired, it uses a divide and conquer algorithm. */
/* > */
/* > The divide and conquer algorithm makes very mild assumptions about */
/* > floating point arithmetic. It will work on machines with a guard */
/* > digit in add/subtract, or on those binary machines without guard */
/* > digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or */
/* > Cray-2. It could conceivably fail on hexadecimal or decimal machines */
/* > without guard digits, but we know of none. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ITYPE */
/* > \verbatim */
/* >          ITYPE is INTEGER */
/* >          Specifies the problem type to be solved: */
/* >          = 1:  A*x = (lambda)*B*x */
/* >          = 2:  A*B*x = (lambda)*x */
/* >          = 3:  B*A*x = (lambda)*x */
/* > \endverbatim */
/* > */
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
/* >          = 'U':  Upper triangles of A and B are stored; */
/* >          = 'L':  Lower triangles of A and B are stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices A and B.  N >= 0. */
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
/* >          On exit, the contents of AP are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] BP */
/* > \verbatim */
/* >          BP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the symmetric matrix */
/* >          B, packed columnwise in a linear array.  The j-th column of B */
/* >          is stored in the array BP as follows: */
/* >          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n. */
/* > */
/* >          On exit, the triangular factor U or L from the Cholesky */
/* >          factorization B = U**T*U or B = L*L**T, in the same storage */
/* >          format as B. */
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
/* >          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of */
/* >          eigenvectors.  The eigenvectors are normalized as follows: */
/* >          if ITYPE = 1 or 2, Z**T*B*Z = I; */
/* >          if ITYPE = 3, Z**T*inv(B)*Z = I. */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the required LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          If N <= 1,               LWORK >= 1. */
/* >          If JOBZ = 'N' and N > 1, LWORK >= 2*N. */
/* >          If JOBZ = 'V' and N > 1, LWORK >= 1 + 6*N + 2*N**2. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the required sizes of the WORK and IWORK */
/* >          arrays, returns these values as the first entries of the WORK */
/* >          and IWORK arrays, and no error message related to LWORK or */
/* >          LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* >          On exit, if INFO = 0, IWORK(1) returns the required LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* >          LIWORK is INTEGER */
/* >          The dimension of the array IWORK. */
/* >          If JOBZ  = 'N' or N <= 1, LIWORK >= 1. */
/* >          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N. */
/* > */
/* >          If LIWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the required sizes of the WORK and */
/* >          IWORK arrays, returns these values as the first entries of */
/* >          the WORK and IWORK arrays, and no error message related to */
/* >          LWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  DPPTRF or DSPEVD returned an error code: */
/* >             <= N:  if INFO = i, DSPEVD failed to converge; */
/* >                    i off-diagonal elements of an intermediate */
/* >                    tridiagonal form did not converge to zero; */
/* >             > N:   if INFO = N + i, for 1 <= i <= N, then the leading */
/* >                    minor of order i of B is not positive definite. */
/* >                    The factorization of B could not be completed and */
/* >                    no eigenvalues or eigenvectors were computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleOTHEReigen */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */

/*  ===================================================================== */
/* Subroutine */ int dspgvd_(integer *itype, char *jobz, char *uplo, integer *
	n, doublereal *ap, doublereal *bp, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *lwork, integer *iwork, 
	integer *liwork, integer *info, ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j, neig;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer lwmin;
    static char trans[1];
    static logical upper;
    extern /* Subroutine */ int dtpmv_(char *, char *, char *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen), 
	    dtpsv_(char *, char *, char *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    static logical wantz;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dspevd_(
	    char *, char *, integer *, doublereal *, doublereal *, doublereal 
	    *, integer *, doublereal *, integer *, integer *, integer *, 
	    integer *, ftnlen, ftnlen);
    static integer liwmin;
    extern /* Subroutine */ int dpptrf_(char *, integer *, doublereal *, 
	    integer *, ftnlen), dspgst_(integer *, char *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static logical lquery;


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

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

#line 249 "dspgvd.f"
    /* Parameter adjustments */
#line 249 "dspgvd.f"
    --ap;
#line 249 "dspgvd.f"
    --bp;
#line 249 "dspgvd.f"
    --w;
#line 249 "dspgvd.f"
    z_dim1 = *ldz;
#line 249 "dspgvd.f"
    z_offset = 1 + z_dim1;
#line 249 "dspgvd.f"
    z__ -= z_offset;
#line 249 "dspgvd.f"
    --work;
#line 249 "dspgvd.f"
    --iwork;
#line 249 "dspgvd.f"

#line 249 "dspgvd.f"
    /* Function Body */
#line 249 "dspgvd.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 250 "dspgvd.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 251 "dspgvd.f"
    lquery = *lwork == -1 || *liwork == -1;

#line 253 "dspgvd.f"
    *info = 0;
#line 254 "dspgvd.f"
    if (*itype < 1 || *itype > 3) {
#line 255 "dspgvd.f"
	*info = -1;
#line 256 "dspgvd.f"
    } else if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 257 "dspgvd.f"
	*info = -2;
#line 258 "dspgvd.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 259 "dspgvd.f"
	*info = -3;
#line 260 "dspgvd.f"
    } else if (*n < 0) {
#line 261 "dspgvd.f"
	*info = -4;
#line 262 "dspgvd.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 263 "dspgvd.f"
	*info = -9;
#line 264 "dspgvd.f"
    }

#line 266 "dspgvd.f"
    if (*info == 0) {
#line 267 "dspgvd.f"
	if (*n <= 1) {
#line 268 "dspgvd.f"
	    liwmin = 1;
#line 269 "dspgvd.f"
	    lwmin = 1;
#line 270 "dspgvd.f"
	} else {
#line 271 "dspgvd.f"
	    if (wantz) {
#line 272 "dspgvd.f"
		liwmin = *n * 5 + 3;
/* Computing 2nd power */
#line 273 "dspgvd.f"
		i__1 = *n;
#line 273 "dspgvd.f"
		lwmin = *n * 6 + 1 + (i__1 * i__1 << 1);
#line 274 "dspgvd.f"
	    } else {
#line 275 "dspgvd.f"
		liwmin = 1;
#line 276 "dspgvd.f"
		lwmin = *n << 1;
#line 277 "dspgvd.f"
	    }
#line 278 "dspgvd.f"
	}
#line 279 "dspgvd.f"
	work[1] = (doublereal) lwmin;
#line 280 "dspgvd.f"
	iwork[1] = liwmin;
#line 281 "dspgvd.f"
	if (*lwork < lwmin && ! lquery) {
#line 282 "dspgvd.f"
	    *info = -11;
#line 283 "dspgvd.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 284 "dspgvd.f"
	    *info = -13;
#line 285 "dspgvd.f"
	}
#line 286 "dspgvd.f"
    }

#line 288 "dspgvd.f"
    if (*info != 0) {
#line 289 "dspgvd.f"
	i__1 = -(*info);
#line 289 "dspgvd.f"
	xerbla_("DSPGVD", &i__1, (ftnlen)6);
#line 290 "dspgvd.f"
	return 0;
#line 291 "dspgvd.f"
    } else if (lquery) {
#line 292 "dspgvd.f"
	return 0;
#line 293 "dspgvd.f"
    }

/*     Quick return if possible */

#line 297 "dspgvd.f"
    if (*n == 0) {
#line 297 "dspgvd.f"
	return 0;
#line 297 "dspgvd.f"
    }

/*     Form a Cholesky factorization of BP. */

#line 302 "dspgvd.f"
    dpptrf_(uplo, n, &bp[1], info, (ftnlen)1);
#line 303 "dspgvd.f"
    if (*info != 0) {
#line 304 "dspgvd.f"
	*info = *n + *info;
#line 305 "dspgvd.f"
	return 0;
#line 306 "dspgvd.f"
    }

/*     Transform problem to standard eigenvalue problem and solve. */

#line 310 "dspgvd.f"
    dspgst_(itype, uplo, n, &ap[1], &bp[1], info, (ftnlen)1);
#line 311 "dspgvd.f"
    dspevd_(jobz, uplo, n, &ap[1], &w[1], &z__[z_offset], ldz, &work[1], 
	    lwork, &iwork[1], liwork, info, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 313 "dspgvd.f"
    d__1 = (doublereal) lwmin;
#line 313 "dspgvd.f"
    lwmin = (integer) max(d__1,work[1]);
/* Computing MAX */
#line 314 "dspgvd.f"
    d__1 = (doublereal) liwmin, d__2 = (doublereal) iwork[1];
#line 314 "dspgvd.f"
    liwmin = (integer) max(d__1,d__2);

#line 316 "dspgvd.f"
    if (wantz) {

/*        Backtransform eigenvectors to the original problem. */

#line 320 "dspgvd.f"
	neig = *n;
#line 321 "dspgvd.f"
	if (*info > 0) {
#line 321 "dspgvd.f"
	    neig = *info - 1;
#line 321 "dspgvd.f"
	}
#line 323 "dspgvd.f"
	if (*itype == 1 || *itype == 2) {

/*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x; */
/*           backtransform eigenvectors: x = inv(L)**T *y or inv(U)*y */

#line 328 "dspgvd.f"
	    if (upper) {
#line 329 "dspgvd.f"
		*(unsigned char *)trans = 'N';
#line 330 "dspgvd.f"
	    } else {
#line 331 "dspgvd.f"
		*(unsigned char *)trans = 'T';
#line 332 "dspgvd.f"
	    }

#line 334 "dspgvd.f"
	    i__1 = neig;
#line 334 "dspgvd.f"
	    for (j = 1; j <= i__1; ++j) {
#line 335 "dspgvd.f"
		dtpsv_(uplo, trans, "Non-unit", n, &bp[1], &z__[j * z_dim1 + 
			1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)8);
#line 337 "dspgvd.f"
/* L10: */
#line 337 "dspgvd.f"
	    }

#line 339 "dspgvd.f"
	} else if (*itype == 3) {

/*           For B*A*x=(lambda)*x; */
/*           backtransform eigenvectors: x = L*y or U**T *y */

#line 344 "dspgvd.f"
	    if (upper) {
#line 345 "dspgvd.f"
		*(unsigned char *)trans = 'T';
#line 346 "dspgvd.f"
	    } else {
#line 347 "dspgvd.f"
		*(unsigned char *)trans = 'N';
#line 348 "dspgvd.f"
	    }

#line 350 "dspgvd.f"
	    i__1 = neig;
#line 350 "dspgvd.f"
	    for (j = 1; j <= i__1; ++j) {
#line 351 "dspgvd.f"
		dtpmv_(uplo, trans, "Non-unit", n, &bp[1], &z__[j * z_dim1 + 
			1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)8);
#line 353 "dspgvd.f"
/* L20: */
#line 353 "dspgvd.f"
	    }
#line 354 "dspgvd.f"
	}
#line 355 "dspgvd.f"
    }

#line 357 "dspgvd.f"
    work[1] = (doublereal) lwmin;
#line 358 "dspgvd.f"
    iwork[1] = liwmin;

#line 360 "dspgvd.f"
    return 0;

/*     End of DSPGVD */

} /* dspgvd_ */

