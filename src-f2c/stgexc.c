#line 1 "stgexc.f"
/* stgexc.f -- translated by f2c (version 20100827).
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

#line 1 "stgexc.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b STGEXC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STGEXC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgexc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgexc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgexc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, */
/*                          LDZ, IFST, ILST, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            WANTQ, WANTZ */
/*       INTEGER            IFST, ILST, INFO, LDA, LDB, LDQ, LDZ, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/*      $                   WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STGEXC reorders the generalized real Schur decomposition of a real */
/* > matrix pair (A,B) using an orthogonal equivalence transformation */
/* > */
/* >                (A, B) = Q * (A, B) * Z**T, */
/* > */
/* > so that the diagonal block of (A, B) with row index IFST is moved */
/* > to row ILST. */
/* > */
/* > (A, B) must be in generalized real Schur canonical form (as returned */
/* > by SGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2 */
/* > diagonal blocks. B is upper triangular. */
/* > */
/* > Optionally, the matrices Q and Z of generalized Schur vectors are */
/* > updated. */
/* > */
/* >        Q(in) * A(in) * Z(in)**T = Q(out) * A(out) * Z(out)**T */
/* >        Q(in) * B(in) * Z(in)**T = Q(out) * B(out) * Z(out)**T */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] WANTQ */
/* > \verbatim */
/* >          WANTQ is LOGICAL */
/* >          .TRUE. : update the left transformation matrix Q; */
/* >          .FALSE.: do not update Q. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* >          WANTZ is LOGICAL */
/* >          .TRUE. : update the right transformation matrix Z; */
/* >          .FALSE.: do not update Z. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices A and B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the matrix A in generalized real Schur canonical */
/* >          form. */
/* >          On exit, the updated matrix A, again in generalized */
/* >          real Schur canonical form. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,N) */
/* >          On entry, the matrix B in generalized real Schur canonical */
/* >          form (A,B). */
/* >          On exit, the updated matrix B, again in generalized */
/* >          real Schur canonical form (A,B). */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is REAL array, dimension (LDZ,N) */
/* >          On entry, if WANTQ = .TRUE., the orthogonal matrix Q. */
/* >          On exit, the updated matrix Q. */
/* >          If WANTQ = .FALSE., Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q. LDQ >= 1. */
/* >          If WANTQ = .TRUE., LDQ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (LDZ,N) */
/* >          On entry, if WANTZ = .TRUE., the orthogonal matrix Z. */
/* >          On exit, the updated matrix Z. */
/* >          If WANTZ = .FALSE., Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z. LDZ >= 1. */
/* >          If WANTZ = .TRUE., LDZ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] IFST */
/* > \verbatim */
/* >          IFST is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in,out] ILST */
/* > \verbatim */
/* >          ILST is INTEGER */
/* >          Specify the reordering of the diagonal blocks of (A, B). */
/* >          The block with row index IFST is moved to row ILST, by a */
/* >          sequence of swapping between adjacent blocks. */
/* >          On exit, if IFST pointed on entry to the second row of */
/* >          a 2-by-2 block, it is changed to point to the first row; */
/* >          ILST always points to the first row of the block in its */
/* >          final position (which may differ from its input value by */
/* >          +1 or -1). 1 <= IFST, ILST <= N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          LWORK >= 1 when N <= 1, otherwise LWORK >= 4*N + 16. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >           =0:  successful exit. */
/* >           <0:  if INFO = -i, the i-th argument had an illegal value. */
/* >           =1:  The transformed matrix pair (A, B) would be too far */
/* >                from generalized Schur form; the problem is ill- */
/* >                conditioned. (A, B) may have been partially reordered, */
/* >                and ILST points to the first row of the current */
/* >                position of the block being moved. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realGEcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* >     Umea University, S-901 87 Umea, Sweden. */

/* > \par References: */
/*  ================ */
/* > */
/* > \verbatim */
/* > */
/* >  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the */
/* >      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in */
/* >      M.S. Moonen et al (eds), Linear Algebra for Large Scale and */
/* >      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int stgexc_(logical *wantq, logical *wantz, integer *n, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	q, integer *ldq, doublereal *z__, integer *ldz, integer *ifst, 
	integer *ilst, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1;

    /* Local variables */
    static integer nbf, nbl, here, lwmin;
    extern /* Subroutine */ int stgex2_(logical *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, doublereal *, integer *, integer *), xerbla_(char *, integer *,
	     ftnlen);
    static integer nbnext;
    static logical lquery;


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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and test input arguments. */

#line 257 "stgexc.f"
    /* Parameter adjustments */
#line 257 "stgexc.f"
    a_dim1 = *lda;
#line 257 "stgexc.f"
    a_offset = 1 + a_dim1;
#line 257 "stgexc.f"
    a -= a_offset;
#line 257 "stgexc.f"
    b_dim1 = *ldb;
#line 257 "stgexc.f"
    b_offset = 1 + b_dim1;
#line 257 "stgexc.f"
    b -= b_offset;
#line 257 "stgexc.f"
    q_dim1 = *ldq;
#line 257 "stgexc.f"
    q_offset = 1 + q_dim1;
#line 257 "stgexc.f"
    q -= q_offset;
#line 257 "stgexc.f"
    z_dim1 = *ldz;
#line 257 "stgexc.f"
    z_offset = 1 + z_dim1;
#line 257 "stgexc.f"
    z__ -= z_offset;
#line 257 "stgexc.f"
    --work;
#line 257 "stgexc.f"

#line 257 "stgexc.f"
    /* Function Body */
#line 257 "stgexc.f"
    *info = 0;
#line 258 "stgexc.f"
    lquery = *lwork == -1;
#line 259 "stgexc.f"
    if (*n < 0) {
#line 260 "stgexc.f"
	*info = -3;
#line 261 "stgexc.f"
    } else if (*lda < max(1,*n)) {
#line 262 "stgexc.f"
	*info = -5;
#line 263 "stgexc.f"
    } else if (*ldb < max(1,*n)) {
#line 264 "stgexc.f"
	*info = -7;
#line 265 "stgexc.f"
    } else if (*ldq < 1 || *wantq && *ldq < max(1,*n)) {
#line 266 "stgexc.f"
	*info = -9;
#line 267 "stgexc.f"
    } else if (*ldz < 1 || *wantz && *ldz < max(1,*n)) {
#line 268 "stgexc.f"
	*info = -11;
#line 269 "stgexc.f"
    } else if (*ifst < 1 || *ifst > *n) {
#line 270 "stgexc.f"
	*info = -12;
#line 271 "stgexc.f"
    } else if (*ilst < 1 || *ilst > *n) {
#line 272 "stgexc.f"
	*info = -13;
#line 273 "stgexc.f"
    }

#line 275 "stgexc.f"
    if (*info == 0) {
#line 276 "stgexc.f"
	if (*n <= 1) {
#line 277 "stgexc.f"
	    lwmin = 1;
#line 278 "stgexc.f"
	} else {
#line 279 "stgexc.f"
	    lwmin = (*n << 2) + 16;
#line 280 "stgexc.f"
	}
#line 281 "stgexc.f"
	work[1] = (doublereal) lwmin;

#line 283 "stgexc.f"
	if (*lwork < lwmin && ! lquery) {
#line 284 "stgexc.f"
	    *info = -15;
#line 285 "stgexc.f"
	}
#line 286 "stgexc.f"
    }

#line 288 "stgexc.f"
    if (*info != 0) {
#line 289 "stgexc.f"
	i__1 = -(*info);
#line 289 "stgexc.f"
	xerbla_("STGEXC", &i__1, (ftnlen)6);
#line 290 "stgexc.f"
	return 0;
#line 291 "stgexc.f"
    } else if (lquery) {
#line 292 "stgexc.f"
	return 0;
#line 293 "stgexc.f"
    }

/*     Quick return if possible */

#line 297 "stgexc.f"
    if (*n <= 1) {
#line 297 "stgexc.f"
	return 0;
#line 297 "stgexc.f"
    }

/*     Determine the first row of the specified block and find out */
/*     if it is 1-by-1 or 2-by-2. */

#line 303 "stgexc.f"
    if (*ifst > 1) {
#line 304 "stgexc.f"
	if (a[*ifst + (*ifst - 1) * a_dim1] != 0.) {
#line 304 "stgexc.f"
	    --(*ifst);
#line 304 "stgexc.f"
	}
#line 306 "stgexc.f"
    }
#line 307 "stgexc.f"
    nbf = 1;
#line 308 "stgexc.f"
    if (*ifst < *n) {
#line 309 "stgexc.f"
	if (a[*ifst + 1 + *ifst * a_dim1] != 0.) {
#line 309 "stgexc.f"
	    nbf = 2;
#line 309 "stgexc.f"
	}
#line 311 "stgexc.f"
    }

/*     Determine the first row of the final block */
/*     and find out if it is 1-by-1 or 2-by-2. */

#line 316 "stgexc.f"
    if (*ilst > 1) {
#line 317 "stgexc.f"
	if (a[*ilst + (*ilst - 1) * a_dim1] != 0.) {
#line 317 "stgexc.f"
	    --(*ilst);
#line 317 "stgexc.f"
	}
#line 319 "stgexc.f"
    }
#line 320 "stgexc.f"
    nbl = 1;
#line 321 "stgexc.f"
    if (*ilst < *n) {
#line 322 "stgexc.f"
	if (a[*ilst + 1 + *ilst * a_dim1] != 0.) {
#line 322 "stgexc.f"
	    nbl = 2;
#line 322 "stgexc.f"
	}
#line 324 "stgexc.f"
    }
#line 325 "stgexc.f"
    if (*ifst == *ilst) {
#line 325 "stgexc.f"
	return 0;
#line 325 "stgexc.f"
    }

#line 328 "stgexc.f"
    if (*ifst < *ilst) {

/*        Update ILST. */

#line 332 "stgexc.f"
	if (nbf == 2 && nbl == 1) {
#line 332 "stgexc.f"
	    --(*ilst);
#line 332 "stgexc.f"
	}
#line 334 "stgexc.f"
	if (nbf == 1 && nbl == 2) {
#line 334 "stgexc.f"
	    ++(*ilst);
#line 334 "stgexc.f"
	}

#line 337 "stgexc.f"
	here = *ifst;

#line 339 "stgexc.f"
L10:

/*        Swap with next one below. */

#line 343 "stgexc.f"
	if (nbf == 1 || nbf == 2) {

/*           Current block either 1-by-1 or 2-by-2. */

#line 347 "stgexc.f"
	    nbnext = 1;
#line 348 "stgexc.f"
	    if (here + nbf + 1 <= *n) {
#line 349 "stgexc.f"
		if (a[here + nbf + 1 + (here + nbf) * a_dim1] != 0.) {
#line 349 "stgexc.f"
		    nbnext = 2;
#line 349 "stgexc.f"
		}
#line 351 "stgexc.f"
	    }
#line 352 "stgexc.f"
	    stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[
		    q_offset], ldq, &z__[z_offset], ldz, &here, &nbf, &nbnext,
		     &work[1], lwork, info);
#line 354 "stgexc.f"
	    if (*info != 0) {
#line 355 "stgexc.f"
		*ilst = here;
#line 356 "stgexc.f"
		return 0;
#line 357 "stgexc.f"
	    }
#line 358 "stgexc.f"
	    here += nbnext;

/*           Test if 2-by-2 block breaks into two 1-by-1 blocks. */

#line 362 "stgexc.f"
	    if (nbf == 2) {
#line 363 "stgexc.f"
		if (a[here + 1 + here * a_dim1] == 0.) {
#line 363 "stgexc.f"
		    nbf = 3;
#line 363 "stgexc.f"
		}
#line 365 "stgexc.f"
	    }

#line 367 "stgexc.f"
	} else {

/*           Current block consists of two 1-by-1 blocks, each of which */
/*           must be swapped individually. */

#line 372 "stgexc.f"
	    nbnext = 1;
#line 373 "stgexc.f"
	    if (here + 3 <= *n) {
#line 374 "stgexc.f"
		if (a[here + 3 + (here + 2) * a_dim1] != 0.) {
#line 374 "stgexc.f"
		    nbnext = 2;
#line 374 "stgexc.f"
		}
#line 376 "stgexc.f"
	    }
#line 377 "stgexc.f"
	    i__1 = here + 1;
#line 377 "stgexc.f"
	    stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[
		    q_offset], ldq, &z__[z_offset], ldz, &i__1, &c__1, &
		    nbnext, &work[1], lwork, info);
#line 379 "stgexc.f"
	    if (*info != 0) {
#line 380 "stgexc.f"
		*ilst = here;
#line 381 "stgexc.f"
		return 0;
#line 382 "stgexc.f"
	    }
#line 383 "stgexc.f"
	    if (nbnext == 1) {

/*              Swap two 1-by-1 blocks. */

#line 387 "stgexc.f"
		stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb,
			 &q[q_offset], ldq, &z__[z_offset], ldz, &here, &c__1,
			 &c__1, &work[1], lwork, info);
#line 389 "stgexc.f"
		if (*info != 0) {
#line 390 "stgexc.f"
		    *ilst = here;
#line 391 "stgexc.f"
		    return 0;
#line 392 "stgexc.f"
		}
#line 393 "stgexc.f"
		++here;

#line 395 "stgexc.f"
	    } else {

/*              Recompute NBNEXT in case of 2-by-2 split. */

#line 399 "stgexc.f"
		if (a[here + 2 + (here + 1) * a_dim1] == 0.) {
#line 399 "stgexc.f"
		    nbnext = 1;
#line 399 "stgexc.f"
		}
#line 401 "stgexc.f"
		if (nbnext == 2) {

/*                 2-by-2 block did not split. */

#line 405 "stgexc.f"
		    stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], 
			    ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &
			    here, &c__1, &nbnext, &work[1], lwork, info);
#line 408 "stgexc.f"
		    if (*info != 0) {
#line 409 "stgexc.f"
			*ilst = here;
#line 410 "stgexc.f"
			return 0;
#line 411 "stgexc.f"
		    }
#line 412 "stgexc.f"
		    here += 2;
#line 413 "stgexc.f"
		} else {

/*                 2-by-2 block did split. */

#line 417 "stgexc.f"
		    stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], 
			    ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &
			    here, &c__1, &c__1, &work[1], lwork, info);
#line 419 "stgexc.f"
		    if (*info != 0) {
#line 420 "stgexc.f"
			*ilst = here;
#line 421 "stgexc.f"
			return 0;
#line 422 "stgexc.f"
		    }
#line 423 "stgexc.f"
		    ++here;
#line 424 "stgexc.f"
		    stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], 
			    ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &
			    here, &c__1, &c__1, &work[1], lwork, info);
#line 426 "stgexc.f"
		    if (*info != 0) {
#line 427 "stgexc.f"
			*ilst = here;
#line 428 "stgexc.f"
			return 0;
#line 429 "stgexc.f"
		    }
#line 430 "stgexc.f"
		    ++here;
#line 431 "stgexc.f"
		}

#line 433 "stgexc.f"
	    }
#line 434 "stgexc.f"
	}
#line 435 "stgexc.f"
	if (here < *ilst) {
#line 435 "stgexc.f"
	    goto L10;
#line 435 "stgexc.f"
	}
#line 437 "stgexc.f"
    } else {
#line 438 "stgexc.f"
	here = *ifst;

#line 440 "stgexc.f"
L20:

/*        Swap with next one below. */

#line 444 "stgexc.f"
	if (nbf == 1 || nbf == 2) {

/*           Current block either 1-by-1 or 2-by-2. */

#line 448 "stgexc.f"
	    nbnext = 1;
#line 449 "stgexc.f"
	    if (here >= 3) {
#line 450 "stgexc.f"
		if (a[here - 1 + (here - 2) * a_dim1] != 0.) {
#line 450 "stgexc.f"
		    nbnext = 2;
#line 450 "stgexc.f"
		}
#line 452 "stgexc.f"
	    }
#line 453 "stgexc.f"
	    i__1 = here - nbnext;
#line 453 "stgexc.f"
	    stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[
		    q_offset], ldq, &z__[z_offset], ldz, &i__1, &nbnext, &nbf,
		     &work[1], lwork, info);
#line 456 "stgexc.f"
	    if (*info != 0) {
#line 457 "stgexc.f"
		*ilst = here;
#line 458 "stgexc.f"
		return 0;
#line 459 "stgexc.f"
	    }
#line 460 "stgexc.f"
	    here -= nbnext;

/*           Test if 2-by-2 block breaks into two 1-by-1 blocks. */

#line 464 "stgexc.f"
	    if (nbf == 2) {
#line 465 "stgexc.f"
		if (a[here + 1 + here * a_dim1] == 0.) {
#line 465 "stgexc.f"
		    nbf = 3;
#line 465 "stgexc.f"
		}
#line 467 "stgexc.f"
	    }

#line 469 "stgexc.f"
	} else {

/*           Current block consists of two 1-by-1 blocks, each of which */
/*           must be swapped individually. */

#line 474 "stgexc.f"
	    nbnext = 1;
#line 475 "stgexc.f"
	    if (here >= 3) {
#line 476 "stgexc.f"
		if (a[here - 1 + (here - 2) * a_dim1] != 0.) {
#line 476 "stgexc.f"
		    nbnext = 2;
#line 476 "stgexc.f"
		}
#line 478 "stgexc.f"
	    }
#line 479 "stgexc.f"
	    i__1 = here - nbnext;
#line 479 "stgexc.f"
	    stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[
		    q_offset], ldq, &z__[z_offset], ldz, &i__1, &nbnext, &
		    c__1, &work[1], lwork, info);
#line 482 "stgexc.f"
	    if (*info != 0) {
#line 483 "stgexc.f"
		*ilst = here;
#line 484 "stgexc.f"
		return 0;
#line 485 "stgexc.f"
	    }
#line 486 "stgexc.f"
	    if (nbnext == 1) {

/*              Swap two 1-by-1 blocks. */

#line 490 "stgexc.f"
		stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb,
			 &q[q_offset], ldq, &z__[z_offset], ldz, &here, &
			nbnext, &c__1, &work[1], lwork, info);
#line 492 "stgexc.f"
		if (*info != 0) {
#line 493 "stgexc.f"
		    *ilst = here;
#line 494 "stgexc.f"
		    return 0;
#line 495 "stgexc.f"
		}
#line 496 "stgexc.f"
		--here;
#line 497 "stgexc.f"
	    } else {

/*             Recompute NBNEXT in case of 2-by-2 split. */

#line 501 "stgexc.f"
		if (a[here + (here - 1) * a_dim1] == 0.) {
#line 501 "stgexc.f"
		    nbnext = 1;
#line 501 "stgexc.f"
		}
#line 503 "stgexc.f"
		if (nbnext == 2) {

/*                 2-by-2 block did not split. */

#line 507 "stgexc.f"
		    i__1 = here - 1;
#line 507 "stgexc.f"
		    stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], 
			    ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &
			    i__1, &c__2, &c__1, &work[1], lwork, info);
#line 509 "stgexc.f"
		    if (*info != 0) {
#line 510 "stgexc.f"
			*ilst = here;
#line 511 "stgexc.f"
			return 0;
#line 512 "stgexc.f"
		    }
#line 513 "stgexc.f"
		    here += -2;
#line 514 "stgexc.f"
		} else {

/*                 2-by-2 block did split. */

#line 518 "stgexc.f"
		    stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], 
			    ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &
			    here, &c__1, &c__1, &work[1], lwork, info);
#line 520 "stgexc.f"
		    if (*info != 0) {
#line 521 "stgexc.f"
			*ilst = here;
#line 522 "stgexc.f"
			return 0;
#line 523 "stgexc.f"
		    }
#line 524 "stgexc.f"
		    --here;
#line 525 "stgexc.f"
		    stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], 
			    ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &
			    here, &c__1, &c__1, &work[1], lwork, info);
#line 527 "stgexc.f"
		    if (*info != 0) {
#line 528 "stgexc.f"
			*ilst = here;
#line 529 "stgexc.f"
			return 0;
#line 530 "stgexc.f"
		    }
#line 531 "stgexc.f"
		    --here;
#line 532 "stgexc.f"
		}
#line 533 "stgexc.f"
	    }
#line 534 "stgexc.f"
	}
#line 535 "stgexc.f"
	if (here > *ilst) {
#line 535 "stgexc.f"
	    goto L20;
#line 535 "stgexc.f"
	}
#line 537 "stgexc.f"
    }
#line 538 "stgexc.f"
    *ilst = here;
#line 539 "stgexc.f"
    work[1] = (doublereal) lwmin;
#line 540 "stgexc.f"
    return 0;

/*     End of STGEXC */

} /* stgexc_ */

