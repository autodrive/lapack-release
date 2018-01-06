#line 1 "zlaed0.f"
/* zlaed0.f -- translated by f2c (version 20100827).
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

#line 1 "zlaed0.f"
/* Table of constant values */

static integer c__9 = 9;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c__1 = 1;

/* > \brief \b ZLAED0 used by sstedc. Computes all eigenvalues and corresponding eigenvectors of an unreduced 
symmetric tridiagonal matrix using the divide and conquer method. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAED0 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaed0.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaed0.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaed0.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAED0( QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS, RWORK, */
/*                          IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDQ, LDQS, N, QSIZ */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   D( * ), E( * ), RWORK( * ) */
/*       COMPLEX*16         Q( LDQ, * ), QSTORE( LDQS, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Using the divide and conquer method, ZLAED0 computes all eigenvalues */
/* > of a symmetric tridiagonal matrix which is one diagonal block of */
/* > those from reducing a dense or band Hermitian matrix and */
/* > corresponding eigenvectors of the dense or band matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] QSIZ */
/* > \verbatim */
/* >          QSIZ is INTEGER */
/* >         The dimension of the unitary matrix used to reduce */
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
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >         On entry, the diagonal elements of the tridiagonal matrix. */
/* >         On exit, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >         On entry, the off-diagonal elements of the tridiagonal matrix. */
/* >         On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is COMPLEX*16 array, dimension (LDQ,N) */
/* >         On entry, Q must contain an QSIZ x N matrix whose columns */
/* >         unitarily orthonormal. It is a part of the unitary matrix */
/* >         that reduces the full dense Hermitian matrix to a */
/* >         (reducible) symmetric tridiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >         The leading dimension of the array Q.  LDQ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, */
/* >         the dimension of IWORK must be at least */
/* >                      6 + 6*N + 5*N*lg N */
/* >                      ( lg( N ) = smallest integer k */
/* >                                  such that 2^k >= N ) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, */
/* >                               dimension (1 + 3*N + 2*N*lg N + 3*N**2) */
/* >                        ( lg( N ) = smallest integer k */
/* >                                    such that 2^k >= N ) */
/* > \endverbatim */
/* > */
/* > \param[out] QSTORE */
/* > \verbatim */
/* >          QSTORE is COMPLEX*16 array, dimension (LDQS, N) */
/* >         Used to store parts of */
/* >         the eigenvector matrix when the updating matrix multiplies */
/* >         take place. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQS */
/* > \verbatim */
/* >          LDQS is INTEGER */
/* >         The leading dimension of the array QSTORE. */
/* >         LDQS >= max(1,N). */
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

/* > \date September 2012 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zlaed0_(integer *qsiz, integer *n, doublereal *d__, 
	doublereal *e, doublecomplex *q, integer *ldq, doublecomplex *qstore, 
	integer *ldqs, doublereal *rwork, integer *iwork, integer *info)
{
    /* System generated locals */
    integer q_dim1, q_offset, qstore_dim1, qstore_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal);
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static integer i__, j, k, ll, iq, lgn, msd2, smm1, spm1, spm2;
    static doublereal temp;
    static integer curr, iperm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer indxq, iwrem, iqptr, tlvls;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zlaed7_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublecomplex *, integer *, doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublecomplex *, doublereal *, integer *, integer *)
	    ;
    static integer igivcl;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zlacrm_(integer *, integer *, doublecomplex *,
	     integer *, doublereal *, integer *, doublecomplex *, integer *, 
	    doublereal *);
    static integer igivnm, submat, curprb, subpbs, igivpt;
    extern /* Subroutine */ int dsteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);
    static integer curlvl, matsiz, iprmpt, smlsiz;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*  Warning:      N could be as big as QSIZ! */

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

#line 191 "zlaed0.f"
    /* Parameter adjustments */
#line 191 "zlaed0.f"
    --d__;
#line 191 "zlaed0.f"
    --e;
#line 191 "zlaed0.f"
    q_dim1 = *ldq;
#line 191 "zlaed0.f"
    q_offset = 1 + q_dim1;
#line 191 "zlaed0.f"
    q -= q_offset;
#line 191 "zlaed0.f"
    qstore_dim1 = *ldqs;
#line 191 "zlaed0.f"
    qstore_offset = 1 + qstore_dim1;
#line 191 "zlaed0.f"
    qstore -= qstore_offset;
#line 191 "zlaed0.f"
    --rwork;
#line 191 "zlaed0.f"
    --iwork;
#line 191 "zlaed0.f"

#line 191 "zlaed0.f"
    /* Function Body */
#line 191 "zlaed0.f"
    *info = 0;

/*     IF( ICOMPQ .LT. 0 .OR. ICOMPQ .GT. 2 ) THEN */
/*        INFO = -1 */
/*     ELSE IF( ( ICOMPQ .EQ. 1 ) .AND. ( QSIZ .LT. MAX( 0, N ) ) ) */
/*    $        THEN */
#line 197 "zlaed0.f"
    if (*qsiz < max(0,*n)) {
#line 198 "zlaed0.f"
	*info = -1;
#line 199 "zlaed0.f"
    } else if (*n < 0) {
#line 200 "zlaed0.f"
	*info = -2;
#line 201 "zlaed0.f"
    } else if (*ldq < max(1,*n)) {
#line 202 "zlaed0.f"
	*info = -6;
#line 203 "zlaed0.f"
    } else if (*ldqs < max(1,*n)) {
#line 204 "zlaed0.f"
	*info = -8;
#line 205 "zlaed0.f"
    }
#line 206 "zlaed0.f"
    if (*info != 0) {
#line 207 "zlaed0.f"
	i__1 = -(*info);
#line 207 "zlaed0.f"
	xerbla_("ZLAED0", &i__1, (ftnlen)6);
#line 208 "zlaed0.f"
	return 0;
#line 209 "zlaed0.f"
    }

/*     Quick return if possible */

#line 213 "zlaed0.f"
    if (*n == 0) {
#line 213 "zlaed0.f"
	return 0;
#line 213 "zlaed0.f"
    }

#line 216 "zlaed0.f"
    smlsiz = ilaenv_(&c__9, "ZLAED0", " ", &c__0, &c__0, &c__0, &c__0, (
	    ftnlen)6, (ftnlen)1);

/*     Determine the size and placement of the submatrices, and save in */
/*     the leading elements of IWORK. */

#line 221 "zlaed0.f"
    iwork[1] = *n;
#line 222 "zlaed0.f"
    subpbs = 1;
#line 223 "zlaed0.f"
    tlvls = 0;
#line 224 "zlaed0.f"
L10:
#line 225 "zlaed0.f"
    if (iwork[subpbs] > smlsiz) {
#line 226 "zlaed0.f"
	for (j = subpbs; j >= 1; --j) {
#line 227 "zlaed0.f"
	    iwork[j * 2] = (iwork[j] + 1) / 2;
#line 228 "zlaed0.f"
	    iwork[(j << 1) - 1] = iwork[j] / 2;
#line 229 "zlaed0.f"
/* L20: */
#line 229 "zlaed0.f"
	}
#line 230 "zlaed0.f"
	++tlvls;
#line 231 "zlaed0.f"
	subpbs <<= 1;
#line 232 "zlaed0.f"
	goto L10;
#line 233 "zlaed0.f"
    }
#line 234 "zlaed0.f"
    i__1 = subpbs;
#line 234 "zlaed0.f"
    for (j = 2; j <= i__1; ++j) {
#line 235 "zlaed0.f"
	iwork[j] += iwork[j - 1];
#line 236 "zlaed0.f"
/* L30: */
#line 236 "zlaed0.f"
    }

/*     Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1 */
/*     using rank-1 modifications (cuts). */

#line 241 "zlaed0.f"
    spm1 = subpbs - 1;
#line 242 "zlaed0.f"
    i__1 = spm1;
#line 242 "zlaed0.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 243 "zlaed0.f"
	submat = iwork[i__] + 1;
#line 244 "zlaed0.f"
	smm1 = submat - 1;
#line 245 "zlaed0.f"
	d__[smm1] -= (d__1 = e[smm1], abs(d__1));
#line 246 "zlaed0.f"
	d__[submat] -= (d__1 = e[smm1], abs(d__1));
#line 247 "zlaed0.f"
/* L40: */
#line 247 "zlaed0.f"
    }

#line 249 "zlaed0.f"
    indxq = (*n << 2) + 3;

/*     Set up workspaces for eigenvalues only/accumulate new vectors */
/*     routine */

#line 254 "zlaed0.f"
    temp = log((doublereal) (*n)) / log(2.);
#line 255 "zlaed0.f"
    lgn = (integer) temp;
#line 256 "zlaed0.f"
    if (pow_ii(&c__2, &lgn) < *n) {
#line 256 "zlaed0.f"
	++lgn;
#line 256 "zlaed0.f"
    }
#line 258 "zlaed0.f"
    if (pow_ii(&c__2, &lgn) < *n) {
#line 258 "zlaed0.f"
	++lgn;
#line 258 "zlaed0.f"
    }
#line 260 "zlaed0.f"
    iprmpt = indxq + *n + 1;
#line 261 "zlaed0.f"
    iperm = iprmpt + *n * lgn;
#line 262 "zlaed0.f"
    iqptr = iperm + *n * lgn;
#line 263 "zlaed0.f"
    igivpt = iqptr + *n + 2;
#line 264 "zlaed0.f"
    igivcl = igivpt + *n * lgn;

#line 266 "zlaed0.f"
    igivnm = 1;
#line 267 "zlaed0.f"
    iq = igivnm + (*n << 1) * lgn;
/* Computing 2nd power */
#line 268 "zlaed0.f"
    i__1 = *n;
#line 268 "zlaed0.f"
    iwrem = iq + i__1 * i__1 + 1;
/*     Initialize pointers */
#line 270 "zlaed0.f"
    i__1 = subpbs;
#line 270 "zlaed0.f"
    for (i__ = 0; i__ <= i__1; ++i__) {
#line 271 "zlaed0.f"
	iwork[iprmpt + i__] = 1;
#line 272 "zlaed0.f"
	iwork[igivpt + i__] = 1;
#line 273 "zlaed0.f"
/* L50: */
#line 273 "zlaed0.f"
    }
#line 274 "zlaed0.f"
    iwork[iqptr] = 1;

/*     Solve each submatrix eigenproblem at the bottom of the divide and */
/*     conquer tree. */

#line 279 "zlaed0.f"
    curr = 0;
#line 280 "zlaed0.f"
    i__1 = spm1;
#line 280 "zlaed0.f"
    for (i__ = 0; i__ <= i__1; ++i__) {
#line 281 "zlaed0.f"
	if (i__ == 0) {
#line 282 "zlaed0.f"
	    submat = 1;
#line 283 "zlaed0.f"
	    matsiz = iwork[1];
#line 284 "zlaed0.f"
	} else {
#line 285 "zlaed0.f"
	    submat = iwork[i__] + 1;
#line 286 "zlaed0.f"
	    matsiz = iwork[i__ + 1] - iwork[i__];
#line 287 "zlaed0.f"
	}
#line 288 "zlaed0.f"
	ll = iq - 1 + iwork[iqptr + curr];
#line 289 "zlaed0.f"
	dsteqr_("I", &matsiz, &d__[submat], &e[submat], &rwork[ll], &matsiz, &
		rwork[1], info, (ftnlen)1);
#line 291 "zlaed0.f"
	zlacrm_(qsiz, &matsiz, &q[submat * q_dim1 + 1], ldq, &rwork[ll], &
		matsiz, &qstore[submat * qstore_dim1 + 1], ldqs, &rwork[iwrem]
		);
/* Computing 2nd power */
#line 294 "zlaed0.f"
	i__2 = matsiz;
#line 294 "zlaed0.f"
	iwork[iqptr + curr + 1] = iwork[iqptr + curr] + i__2 * i__2;
#line 295 "zlaed0.f"
	++curr;
#line 296 "zlaed0.f"
	if (*info > 0) {
#line 297 "zlaed0.f"
	    *info = submat * (*n + 1) + submat + matsiz - 1;
#line 298 "zlaed0.f"
	    return 0;
#line 299 "zlaed0.f"
	}
#line 300 "zlaed0.f"
	k = 1;
#line 301 "zlaed0.f"
	i__2 = iwork[i__ + 1];
#line 301 "zlaed0.f"
	for (j = submat; j <= i__2; ++j) {
#line 302 "zlaed0.f"
	    iwork[indxq + j] = k;
#line 303 "zlaed0.f"
	    ++k;
#line 304 "zlaed0.f"
/* L60: */
#line 304 "zlaed0.f"
	}
#line 305 "zlaed0.f"
/* L70: */
#line 305 "zlaed0.f"
    }

/*     Successively merge eigensystems of adjacent submatrices */
/*     into eigensystem for the corresponding larger matrix. */

/*     while ( SUBPBS > 1 ) */

#line 312 "zlaed0.f"
    curlvl = 1;
#line 313 "zlaed0.f"
L80:
#line 314 "zlaed0.f"
    if (subpbs > 1) {
#line 315 "zlaed0.f"
	spm2 = subpbs - 2;
#line 316 "zlaed0.f"
	i__1 = spm2;
#line 316 "zlaed0.f"
	for (i__ = 0; i__ <= i__1; i__ += 2) {
#line 317 "zlaed0.f"
	    if (i__ == 0) {
#line 318 "zlaed0.f"
		submat = 1;
#line 319 "zlaed0.f"
		matsiz = iwork[2];
#line 320 "zlaed0.f"
		msd2 = iwork[1];
#line 321 "zlaed0.f"
		curprb = 0;
#line 322 "zlaed0.f"
	    } else {
#line 323 "zlaed0.f"
		submat = iwork[i__] + 1;
#line 324 "zlaed0.f"
		matsiz = iwork[i__ + 2] - iwork[i__];
#line 325 "zlaed0.f"
		msd2 = matsiz / 2;
#line 326 "zlaed0.f"
		++curprb;
#line 327 "zlaed0.f"
	    }

/*     Merge lower order eigensystems (of size MSD2 and MATSIZ - MSD2) */
/*     into an eigensystem of size MATSIZ.  ZLAED7 handles the case */
/*     when the eigenvectors of a full or band Hermitian matrix (which */
/*     was reduced to tridiagonal form) are desired. */

/*     I am free to use Q as a valuable working space until Loop 150. */

#line 336 "zlaed0.f"
	    zlaed7_(&matsiz, &msd2, qsiz, &tlvls, &curlvl, &curprb, &d__[
		    submat], &qstore[submat * qstore_dim1 + 1], ldqs, &e[
		    submat + msd2 - 1], &iwork[indxq + submat], &rwork[iq], &
		    iwork[iqptr], &iwork[iprmpt], &iwork[iperm], &iwork[
		    igivpt], &iwork[igivcl], &rwork[igivnm], &q[submat * 
		    q_dim1 + 1], &rwork[iwrem], &iwork[subpbs + 1], info);
#line 344 "zlaed0.f"
	    if (*info > 0) {
#line 345 "zlaed0.f"
		*info = submat * (*n + 1) + submat + matsiz - 1;
#line 346 "zlaed0.f"
		return 0;
#line 347 "zlaed0.f"
	    }
#line 348 "zlaed0.f"
	    iwork[i__ / 2 + 1] = iwork[i__ + 2];
#line 349 "zlaed0.f"
/* L90: */
#line 349 "zlaed0.f"
	}
#line 350 "zlaed0.f"
	subpbs /= 2;
#line 351 "zlaed0.f"
	++curlvl;
#line 352 "zlaed0.f"
	goto L80;
#line 353 "zlaed0.f"
    }

/*     end while */

/*     Re-merge the eigenvalues/vectors which were deflated at the final */
/*     merge step. */

#line 360 "zlaed0.f"
    i__1 = *n;
#line 360 "zlaed0.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 361 "zlaed0.f"
	j = iwork[indxq + i__];
#line 362 "zlaed0.f"
	rwork[i__] = d__[j];
#line 363 "zlaed0.f"
	zcopy_(qsiz, &qstore[j * qstore_dim1 + 1], &c__1, &q[i__ * q_dim1 + 1]
		, &c__1);
#line 364 "zlaed0.f"
/* L100: */
#line 364 "zlaed0.f"
    }
#line 365 "zlaed0.f"
    dcopy_(n, &rwork[1], &c__1, &d__[1], &c__1);

#line 367 "zlaed0.f"
    return 0;

/*     End of ZLAED0 */

} /* zlaed0_ */

