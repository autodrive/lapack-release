#line 1 "slaed0.f"
/* slaed0.f -- translated by f2c (version 20100827).
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

#line 1 "slaed0.f"
/* Table of constant values */

static integer c__9 = 9;
static integer c__0 = 0;
static integer c__2 = 2;
static doublereal c_b23 = 1.;
static doublereal c_b24 = 0.;
static integer c__1 = 1;

/* > \brief \b SLAED0 used by sstedc. Computes all eigenvalues and corresponding eigenvectors of an unreduced 
symmetric tridiagonal matrix using the divide and conquer method. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAED0 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed0.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed0.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed0.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAED0( ICOMPQ, QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS, */
/*                          WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            ICOMPQ, INFO, LDQ, LDQS, N, QSIZ */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               D( * ), E( * ), Q( LDQ, * ), QSTORE( LDQS, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAED0 computes all eigenvalues and corresponding eigenvectors of a */
/* > symmetric tridiagonal matrix using the divide and conquer method. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ICOMPQ */
/* > \verbatim */
/* >          ICOMPQ is INTEGER */
/* >          = 0:  Compute eigenvalues only. */
/* >          = 1:  Compute eigenvectors of original dense symmetric matrix */
/* >                also.  On entry, Q contains the orthogonal matrix used */
/* >                to reduce the original matrix to tridiagonal form. */
/* >          = 2:  Compute eigenvalues and eigenvectors of tridiagonal */
/* >                matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] QSIZ */
/* > \verbatim */
/* >          QSIZ is INTEGER */
/* >         The dimension of the orthogonal matrix used to reduce */
/* >         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         The dimension of the symmetric tridiagonal matrix.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >         On entry, the main diagonal of the tridiagonal matrix. */
/* >         On exit, its eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N-1) */
/* >         The off-diagonal elements of the tridiagonal matrix. */
/* >         On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is REAL array, dimension (LDQ, N) */
/* >         On entry, Q must contain an N-by-N orthogonal matrix. */
/* >         If ICOMPQ = 0    Q is not referenced. */
/* >         If ICOMPQ = 1    On entry, Q is a subset of the columns of the */
/* >                          orthogonal matrix used to reduce the full */
/* >                          matrix to tridiagonal form corresponding to */
/* >                          the subset of the full matrix which is being */
/* >                          decomposed at this time. */
/* >         If ICOMPQ = 2    On entry, Q will be the identity matrix. */
/* >                          On exit, Q contains the eigenvectors of the */
/* >                          tridiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >         The leading dimension of the array Q.  If eigenvectors are */
/* >         desired, then  LDQ >= max(1,N).  In any case,  LDQ >= 1. */
/* > \endverbatim */
/* > */
/* > \param[out] QSTORE */
/* > \verbatim */
/* >          QSTORE is REAL array, dimension (LDQS, N) */
/* >         Referenced only when ICOMPQ = 1.  Used to store parts of */
/* >         the eigenvector matrix when the updating matrix multiplies */
/* >         take place. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQS */
/* > \verbatim */
/* >          LDQS is INTEGER */
/* >         The leading dimension of the array QSTORE.  If ICOMPQ = 1, */
/* >         then  LDQS >= max(1,N).  In any case,  LDQS >= 1. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, */
/* >         If ICOMPQ = 0 or 1, the dimension of WORK must be at least */
/* >                     1 + 3*N + 2*N*lg N + 3*N**2 */
/* >                     ( lg( N ) = smallest integer k */
/* >                                 such that 2^k >= N ) */
/* >         If ICOMPQ = 2, the dimension of WORK must be at least */
/* >                     4*N + N**2. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, */
/* >         If ICOMPQ = 0 or 1, the dimension of IWORK must be at least */
/* >                        6 + 6*N + 5*N*lg N. */
/* >                        ( lg( N ) = smallest integer k */
/* >                                    such that 2^k >= N ) */
/* >         If ICOMPQ = 2, the dimension of IWORK must be at least */
/* >                        3 + 5*N. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  The algorithm failed to compute an eigenvalue while */
/* >                working on the submatrix lying in rows and columns */
/* >                INFO/(N+1) through mod(INFO,N+1). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup auxOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int slaed0_(integer *icompq, integer *qsiz, integer *n, 
	doublereal *d__, doublereal *e, doublereal *q, integer *ldq, 
	doublereal *qstore, integer *ldqs, doublereal *work, integer *iwork, 
	integer *info)
{
    /* System generated locals */
    integer q_dim1, q_offset, qstore_dim1, qstore_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal);
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static integer i__, j, k, iq, lgn, msd2, smm1, spm1, spm2;
    static doublereal temp;
    static integer curr;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer iperm, indxq, iwrem;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer iqptr, tlvls;
    extern /* Subroutine */ int slaed1_(integer *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *), slaed7_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, integer *, integer *);
    static integer igivcl;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer igivnm, submat;
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static integer curprb, subpbs, igivpt, curlvl, matsiz, iprmpt, smlsiz;
    extern /* Subroutine */ int ssteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 217 "slaed0.f"
    /* Parameter adjustments */
#line 217 "slaed0.f"
    --d__;
#line 217 "slaed0.f"
    --e;
#line 217 "slaed0.f"
    q_dim1 = *ldq;
#line 217 "slaed0.f"
    q_offset = 1 + q_dim1;
#line 217 "slaed0.f"
    q -= q_offset;
#line 217 "slaed0.f"
    qstore_dim1 = *ldqs;
#line 217 "slaed0.f"
    qstore_offset = 1 + qstore_dim1;
#line 217 "slaed0.f"
    qstore -= qstore_offset;
#line 217 "slaed0.f"
    --work;
#line 217 "slaed0.f"
    --iwork;
#line 217 "slaed0.f"

#line 217 "slaed0.f"
    /* Function Body */
#line 217 "slaed0.f"
    *info = 0;

#line 219 "slaed0.f"
    if (*icompq < 0 || *icompq > 2) {
#line 220 "slaed0.f"
	*info = -1;
#line 221 "slaed0.f"
    } else if (*icompq == 1 && *qsiz < max(0,*n)) {
#line 222 "slaed0.f"
	*info = -2;
#line 223 "slaed0.f"
    } else if (*n < 0) {
#line 224 "slaed0.f"
	*info = -3;
#line 225 "slaed0.f"
    } else if (*ldq < max(1,*n)) {
#line 226 "slaed0.f"
	*info = -7;
#line 227 "slaed0.f"
    } else if (*ldqs < max(1,*n)) {
#line 228 "slaed0.f"
	*info = -9;
#line 229 "slaed0.f"
    }
#line 230 "slaed0.f"
    if (*info != 0) {
#line 231 "slaed0.f"
	i__1 = -(*info);
#line 231 "slaed0.f"
	xerbla_("SLAED0", &i__1, (ftnlen)6);
#line 232 "slaed0.f"
	return 0;
#line 233 "slaed0.f"
    }

/*     Quick return if possible */

#line 237 "slaed0.f"
    if (*n == 0) {
#line 237 "slaed0.f"
	return 0;
#line 237 "slaed0.f"
    }

#line 240 "slaed0.f"
    smlsiz = ilaenv_(&c__9, "SLAED0", " ", &c__0, &c__0, &c__0, &c__0, (
	    ftnlen)6, (ftnlen)1);

/*     Determine the size and placement of the submatrices, and save in */
/*     the leading elements of IWORK. */

#line 245 "slaed0.f"
    iwork[1] = *n;
#line 246 "slaed0.f"
    subpbs = 1;
#line 247 "slaed0.f"
    tlvls = 0;
#line 248 "slaed0.f"
L10:
#line 249 "slaed0.f"
    if (iwork[subpbs] > smlsiz) {
#line 250 "slaed0.f"
	for (j = subpbs; j >= 1; --j) {
#line 251 "slaed0.f"
	    iwork[j * 2] = (iwork[j] + 1) / 2;
#line 252 "slaed0.f"
	    iwork[(j << 1) - 1] = iwork[j] / 2;
#line 253 "slaed0.f"
/* L20: */
#line 253 "slaed0.f"
	}
#line 254 "slaed0.f"
	++tlvls;
#line 255 "slaed0.f"
	subpbs <<= 1;
#line 256 "slaed0.f"
	goto L10;
#line 257 "slaed0.f"
    }
#line 258 "slaed0.f"
    i__1 = subpbs;
#line 258 "slaed0.f"
    for (j = 2; j <= i__1; ++j) {
#line 259 "slaed0.f"
	iwork[j] += iwork[j - 1];
#line 260 "slaed0.f"
/* L30: */
#line 260 "slaed0.f"
    }

/*     Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1 */
/*     using rank-1 modifications (cuts). */

#line 265 "slaed0.f"
    spm1 = subpbs - 1;
#line 266 "slaed0.f"
    i__1 = spm1;
#line 266 "slaed0.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 267 "slaed0.f"
	submat = iwork[i__] + 1;
#line 268 "slaed0.f"
	smm1 = submat - 1;
#line 269 "slaed0.f"
	d__[smm1] -= (d__1 = e[smm1], abs(d__1));
#line 270 "slaed0.f"
	d__[submat] -= (d__1 = e[smm1], abs(d__1));
#line 271 "slaed0.f"
/* L40: */
#line 271 "slaed0.f"
    }

#line 273 "slaed0.f"
    indxq = (*n << 2) + 3;
#line 274 "slaed0.f"
    if (*icompq != 2) {

/*        Set up workspaces for eigenvalues only/accumulate new vectors */
/*        routine */

#line 279 "slaed0.f"
	temp = log((doublereal) (*n)) / log(2.);
#line 280 "slaed0.f"
	lgn = (integer) temp;
#line 281 "slaed0.f"
	if (pow_ii(&c__2, &lgn) < *n) {
#line 281 "slaed0.f"
	    ++lgn;
#line 281 "slaed0.f"
	}
#line 283 "slaed0.f"
	if (pow_ii(&c__2, &lgn) < *n) {
#line 283 "slaed0.f"
	    ++lgn;
#line 283 "slaed0.f"
	}
#line 285 "slaed0.f"
	iprmpt = indxq + *n + 1;
#line 286 "slaed0.f"
	iperm = iprmpt + *n * lgn;
#line 287 "slaed0.f"
	iqptr = iperm + *n * lgn;
#line 288 "slaed0.f"
	igivpt = iqptr + *n + 2;
#line 289 "slaed0.f"
	igivcl = igivpt + *n * lgn;

#line 291 "slaed0.f"
	igivnm = 1;
#line 292 "slaed0.f"
	iq = igivnm + (*n << 1) * lgn;
/* Computing 2nd power */
#line 293 "slaed0.f"
	i__1 = *n;
#line 293 "slaed0.f"
	iwrem = iq + i__1 * i__1 + 1;

/*        Initialize pointers */

#line 297 "slaed0.f"
	i__1 = subpbs;
#line 297 "slaed0.f"
	for (i__ = 0; i__ <= i__1; ++i__) {
#line 298 "slaed0.f"
	    iwork[iprmpt + i__] = 1;
#line 299 "slaed0.f"
	    iwork[igivpt + i__] = 1;
#line 300 "slaed0.f"
/* L50: */
#line 300 "slaed0.f"
	}
#line 301 "slaed0.f"
	iwork[iqptr] = 1;
#line 302 "slaed0.f"
    }

/*     Solve each submatrix eigenproblem at the bottom of the divide and */
/*     conquer tree. */

#line 307 "slaed0.f"
    curr = 0;
#line 308 "slaed0.f"
    i__1 = spm1;
#line 308 "slaed0.f"
    for (i__ = 0; i__ <= i__1; ++i__) {
#line 309 "slaed0.f"
	if (i__ == 0) {
#line 310 "slaed0.f"
	    submat = 1;
#line 311 "slaed0.f"
	    matsiz = iwork[1];
#line 312 "slaed0.f"
	} else {
#line 313 "slaed0.f"
	    submat = iwork[i__] + 1;
#line 314 "slaed0.f"
	    matsiz = iwork[i__ + 1] - iwork[i__];
#line 315 "slaed0.f"
	}
#line 316 "slaed0.f"
	if (*icompq == 2) {
#line 317 "slaed0.f"
	    ssteqr_("I", &matsiz, &d__[submat], &e[submat], &q[submat + 
		    submat * q_dim1], ldq, &work[1], info, (ftnlen)1);
#line 319 "slaed0.f"
	    if (*info != 0) {
#line 319 "slaed0.f"
		goto L130;
#line 319 "slaed0.f"
	    }
#line 321 "slaed0.f"
	} else {
#line 322 "slaed0.f"
	    ssteqr_("I", &matsiz, &d__[submat], &e[submat], &work[iq - 1 + 
		    iwork[iqptr + curr]], &matsiz, &work[1], info, (ftnlen)1);
#line 325 "slaed0.f"
	    if (*info != 0) {
#line 325 "slaed0.f"
		goto L130;
#line 325 "slaed0.f"
	    }
#line 327 "slaed0.f"
	    if (*icompq == 1) {
#line 328 "slaed0.f"
		sgemm_("N", "N", qsiz, &matsiz, &matsiz, &c_b23, &q[submat * 
			q_dim1 + 1], ldq, &work[iq - 1 + iwork[iqptr + curr]],
			 &matsiz, &c_b24, &qstore[submat * qstore_dim1 + 1], 
			ldqs, (ftnlen)1, (ftnlen)1);
#line 332 "slaed0.f"
	    }
/* Computing 2nd power */
#line 333 "slaed0.f"
	    i__2 = matsiz;
#line 333 "slaed0.f"
	    iwork[iqptr + curr + 1] = iwork[iqptr + curr] + i__2 * i__2;
#line 334 "slaed0.f"
	    ++curr;
#line 335 "slaed0.f"
	}
#line 336 "slaed0.f"
	k = 1;
#line 337 "slaed0.f"
	i__2 = iwork[i__ + 1];
#line 337 "slaed0.f"
	for (j = submat; j <= i__2; ++j) {
#line 338 "slaed0.f"
	    iwork[indxq + j] = k;
#line 339 "slaed0.f"
	    ++k;
#line 340 "slaed0.f"
/* L60: */
#line 340 "slaed0.f"
	}
#line 341 "slaed0.f"
/* L70: */
#line 341 "slaed0.f"
    }

/*     Successively merge eigensystems of adjacent submatrices */
/*     into eigensystem for the corresponding larger matrix. */

/*     while ( SUBPBS > 1 ) */

#line 348 "slaed0.f"
    curlvl = 1;
#line 349 "slaed0.f"
L80:
#line 350 "slaed0.f"
    if (subpbs > 1) {
#line 351 "slaed0.f"
	spm2 = subpbs - 2;
#line 352 "slaed0.f"
	i__1 = spm2;
#line 352 "slaed0.f"
	for (i__ = 0; i__ <= i__1; i__ += 2) {
#line 353 "slaed0.f"
	    if (i__ == 0) {
#line 354 "slaed0.f"
		submat = 1;
#line 355 "slaed0.f"
		matsiz = iwork[2];
#line 356 "slaed0.f"
		msd2 = iwork[1];
#line 357 "slaed0.f"
		curprb = 0;
#line 358 "slaed0.f"
	    } else {
#line 359 "slaed0.f"
		submat = iwork[i__] + 1;
#line 360 "slaed0.f"
		matsiz = iwork[i__ + 2] - iwork[i__];
#line 361 "slaed0.f"
		msd2 = matsiz / 2;
#line 362 "slaed0.f"
		++curprb;
#line 363 "slaed0.f"
	    }

/*     Merge lower order eigensystems (of size MSD2 and MATSIZ - MSD2) */
/*     into an eigensystem of size MATSIZ. */
/*     SLAED1 is used only for the full eigensystem of a tridiagonal */
/*     matrix. */
/*     SLAED7 handles the cases in which eigenvalues only or eigenvalues */
/*     and eigenvectors of a full symmetric matrix (which was reduced to */
/*     tridiagonal form) are desired. */

#line 373 "slaed0.f"
	    if (*icompq == 2) {
#line 374 "slaed0.f"
		slaed1_(&matsiz, &d__[submat], &q[submat + submat * q_dim1], 
			ldq, &iwork[indxq + submat], &e[submat + msd2 - 1], &
			msd2, &work[1], &iwork[subpbs + 1], info);
#line 378 "slaed0.f"
	    } else {
#line 379 "slaed0.f"
		slaed7_(icompq, &matsiz, qsiz, &tlvls, &curlvl, &curprb, &d__[
			submat], &qstore[submat * qstore_dim1 + 1], ldqs, &
			iwork[indxq + submat], &e[submat + msd2 - 1], &msd2, &
			work[iq], &iwork[iqptr], &iwork[iprmpt], &iwork[iperm]
			, &iwork[igivpt], &iwork[igivcl], &work[igivnm], &
			work[iwrem], &iwork[subpbs + 1], info);
#line 387 "slaed0.f"
	    }
#line 388 "slaed0.f"
	    if (*info != 0) {
#line 388 "slaed0.f"
		goto L130;
#line 388 "slaed0.f"
	    }
#line 390 "slaed0.f"
	    iwork[i__ / 2 + 1] = iwork[i__ + 2];
#line 391 "slaed0.f"
/* L90: */
#line 391 "slaed0.f"
	}
#line 392 "slaed0.f"
	subpbs /= 2;
#line 393 "slaed0.f"
	++curlvl;
#line 394 "slaed0.f"
	goto L80;
#line 395 "slaed0.f"
    }

/*     end while */

/*     Re-merge the eigenvalues/vectors which were deflated at the final */
/*     merge step. */

#line 402 "slaed0.f"
    if (*icompq == 1) {
#line 403 "slaed0.f"
	i__1 = *n;
#line 403 "slaed0.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 404 "slaed0.f"
	    j = iwork[indxq + i__];
#line 405 "slaed0.f"
	    work[i__] = d__[j];
#line 406 "slaed0.f"
	    scopy_(qsiz, &qstore[j * qstore_dim1 + 1], &c__1, &q[i__ * q_dim1 
		    + 1], &c__1);
#line 407 "slaed0.f"
/* L100: */
#line 407 "slaed0.f"
	}
#line 408 "slaed0.f"
	scopy_(n, &work[1], &c__1, &d__[1], &c__1);
#line 409 "slaed0.f"
    } else if (*icompq == 2) {
#line 410 "slaed0.f"
	i__1 = *n;
#line 410 "slaed0.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 411 "slaed0.f"
	    j = iwork[indxq + i__];
#line 412 "slaed0.f"
	    work[i__] = d__[j];
#line 413 "slaed0.f"
	    scopy_(n, &q[j * q_dim1 + 1], &c__1, &work[*n * i__ + 1], &c__1);
#line 414 "slaed0.f"
/* L110: */
#line 414 "slaed0.f"
	}
#line 415 "slaed0.f"
	scopy_(n, &work[1], &c__1, &d__[1], &c__1);
#line 416 "slaed0.f"
	slacpy_("A", n, n, &work[*n + 1], n, &q[q_offset], ldq, (ftnlen)1);
#line 417 "slaed0.f"
    } else {
#line 418 "slaed0.f"
	i__1 = *n;
#line 418 "slaed0.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 419 "slaed0.f"
	    j = iwork[indxq + i__];
#line 420 "slaed0.f"
	    work[i__] = d__[j];
#line 421 "slaed0.f"
/* L120: */
#line 421 "slaed0.f"
	}
#line 422 "slaed0.f"
	scopy_(n, &work[1], &c__1, &d__[1], &c__1);
#line 423 "slaed0.f"
    }
#line 424 "slaed0.f"
    goto L140;

#line 426 "slaed0.f"
L130:
#line 427 "slaed0.f"
    *info = submat * (*n + 1) + submat + matsiz - 1;

#line 429 "slaed0.f"
L140:
#line 430 "slaed0.f"
    return 0;

/*     End of SLAED0 */

} /* slaed0_ */

