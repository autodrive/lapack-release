#line 1 "dtrexc.f"
/* dtrexc.f -- translated by f2c (version 20100827).
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

#line 1 "dtrexc.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b DTREXC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTREXC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrexc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrexc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrexc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPQ */
/*       INTEGER            IFST, ILST, INFO, LDQ, LDT, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTREXC reorders the real Schur factorization of a real matrix */
/* > A = Q*T*Q**T, so that the diagonal block of T with row index IFST is */
/* > moved to row ILST. */
/* > */
/* > The real Schur form T is reordered by an orthogonal similarity */
/* > transformation Z**T*T*Z, and optionally the matrix Q of Schur vectors */
/* > is updated by postmultiplying it with Z. */
/* > */
/* > T must be in Schur canonical form (as returned by DHSEQR), that is, */
/* > block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each */
/* > 2-by-2 diagonal block has its diagonal elements equal and its */
/* > off-diagonal elements of opposite sign. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] COMPQ */
/* > \verbatim */
/* >          COMPQ is CHARACTER*1 */
/* >          = 'V':  update the matrix Q of Schur vectors; */
/* >          = 'N':  do not update Q. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix T. N >= 0. */
/* >          If N == 0 arguments ILST and IFST may be any value. */
/* > \endverbatim */
/* > */
/* > \param[in,out] T */
/* > \verbatim */
/* >          T is DOUBLE PRECISION array, dimension (LDT,N) */
/* >          On entry, the upper quasi-triangular matrix T, in Schur */
/* >          Schur canonical form. */
/* >          On exit, the reordered upper quasi-triangular matrix, again */
/* >          in Schur canonical form. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T. LDT >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is DOUBLE PRECISION array, dimension (LDQ,N) */
/* >          On entry, if COMPQ = 'V', the matrix Q of Schur vectors. */
/* >          On exit, if COMPQ = 'V', Q has been postmultiplied by the */
/* >          orthogonal transformation matrix Z which reorders T. */
/* >          If COMPQ = 'N', Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q.  LDQ >= 1, and if */
/* >          COMPQ = 'V', LDQ >= max(1,N). */
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
/* > */
/* >          Specify the reordering of the diagonal blocks of T. */
/* >          The block with row index IFST is moved to row ILST, by a */
/* >          sequence of transpositions between adjacent blocks. */
/* >          On exit, if IFST pointed on entry to the second row of a */
/* >          2-by-2 block, it is changed to point to the first row; ILST */
/* >          always points to the first row of the block in its final */
/* >          position (which may differ from its input value by +1 or -1). */
/* >          1 <= IFST <= N; 1 <= ILST <= N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          = 1:  two adjacent blocks were too close to swap (the problem */
/* >                is very ill-conditioned); T may have been partially */
/* >                reordered, and ILST points to the first row of the */
/* >                current position of the block being moved. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dtrexc_(char *compq, integer *n, doublereal *t, integer *
	ldt, doublereal *q, integer *ldq, integer *ifst, integer *ilst, 
	doublereal *work, integer *info, ftnlen compq_len)
{
    /* System generated locals */
    integer q_dim1, q_offset, t_dim1, t_offset, i__1;

    /* Local variables */
    static integer nbf, nbl, here;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical wantq;
    extern /* Subroutine */ int dlaexc_(logical *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, doublereal *, integer *), xerbla_(char *, integer *, ftnlen);
    static integer nbnext;


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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and test the input arguments. */

#line 188 "dtrexc.f"
    /* Parameter adjustments */
#line 188 "dtrexc.f"
    t_dim1 = *ldt;
#line 188 "dtrexc.f"
    t_offset = 1 + t_dim1;
#line 188 "dtrexc.f"
    t -= t_offset;
#line 188 "dtrexc.f"
    q_dim1 = *ldq;
#line 188 "dtrexc.f"
    q_offset = 1 + q_dim1;
#line 188 "dtrexc.f"
    q -= q_offset;
#line 188 "dtrexc.f"
    --work;
#line 188 "dtrexc.f"

#line 188 "dtrexc.f"
    /* Function Body */
#line 188 "dtrexc.f"
    *info = 0;
#line 189 "dtrexc.f"
    wantq = lsame_(compq, "V", (ftnlen)1, (ftnlen)1);
#line 190 "dtrexc.f"
    if (! wantq && ! lsame_(compq, "N", (ftnlen)1, (ftnlen)1)) {
#line 191 "dtrexc.f"
	*info = -1;
#line 192 "dtrexc.f"
    } else if (*n < 0) {
#line 193 "dtrexc.f"
	*info = -2;
#line 194 "dtrexc.f"
    } else if (*ldt < max(1,*n)) {
#line 195 "dtrexc.f"
	*info = -4;
#line 196 "dtrexc.f"
    } else if (*ldq < 1 || wantq && *ldq < max(1,*n)) {
#line 197 "dtrexc.f"
	*info = -6;
#line 198 "dtrexc.f"
    } else if ((*ifst < 1 || *ifst > *n) && *n > 0) {
#line 199 "dtrexc.f"
	*info = -7;
#line 200 "dtrexc.f"
    } else if ((*ilst < 1 || *ilst > *n) && *n > 0) {
#line 201 "dtrexc.f"
	*info = -8;
#line 202 "dtrexc.f"
    }
#line 203 "dtrexc.f"
    if (*info != 0) {
#line 204 "dtrexc.f"
	i__1 = -(*info);
#line 204 "dtrexc.f"
	xerbla_("DTREXC", &i__1, (ftnlen)6);
#line 205 "dtrexc.f"
	return 0;
#line 206 "dtrexc.f"
    }

/*     Quick return if possible */

#line 210 "dtrexc.f"
    if (*n <= 1) {
#line 210 "dtrexc.f"
	return 0;
#line 210 "dtrexc.f"
    }

/*     Determine the first row of specified block */
/*     and find out it is 1 by 1 or 2 by 2. */

#line 216 "dtrexc.f"
    if (*ifst > 1) {
#line 217 "dtrexc.f"
	if (t[*ifst + (*ifst - 1) * t_dim1] != 0.) {
#line 217 "dtrexc.f"
	    --(*ifst);
#line 217 "dtrexc.f"
	}
#line 219 "dtrexc.f"
    }
#line 220 "dtrexc.f"
    nbf = 1;
#line 221 "dtrexc.f"
    if (*ifst < *n) {
#line 222 "dtrexc.f"
	if (t[*ifst + 1 + *ifst * t_dim1] != 0.) {
#line 222 "dtrexc.f"
	    nbf = 2;
#line 222 "dtrexc.f"
	}
#line 224 "dtrexc.f"
    }

/*     Determine the first row of the final block */
/*     and find out it is 1 by 1 or 2 by 2. */

#line 229 "dtrexc.f"
    if (*ilst > 1) {
#line 230 "dtrexc.f"
	if (t[*ilst + (*ilst - 1) * t_dim1] != 0.) {
#line 230 "dtrexc.f"
	    --(*ilst);
#line 230 "dtrexc.f"
	}
#line 232 "dtrexc.f"
    }
#line 233 "dtrexc.f"
    nbl = 1;
#line 234 "dtrexc.f"
    if (*ilst < *n) {
#line 235 "dtrexc.f"
	if (t[*ilst + 1 + *ilst * t_dim1] != 0.) {
#line 235 "dtrexc.f"
	    nbl = 2;
#line 235 "dtrexc.f"
	}
#line 237 "dtrexc.f"
    }

#line 239 "dtrexc.f"
    if (*ifst == *ilst) {
#line 239 "dtrexc.f"
	return 0;
#line 239 "dtrexc.f"
    }

#line 242 "dtrexc.f"
    if (*ifst < *ilst) {

/*        Update ILST */

#line 246 "dtrexc.f"
	if (nbf == 2 && nbl == 1) {
#line 246 "dtrexc.f"
	    --(*ilst);
#line 246 "dtrexc.f"
	}
#line 248 "dtrexc.f"
	if (nbf == 1 && nbl == 2) {
#line 248 "dtrexc.f"
	    ++(*ilst);
#line 248 "dtrexc.f"
	}

#line 251 "dtrexc.f"
	here = *ifst;

#line 253 "dtrexc.f"
L10:

/*        Swap block with next one below */

#line 257 "dtrexc.f"
	if (nbf == 1 || nbf == 2) {

/*           Current block either 1 by 1 or 2 by 2 */

#line 261 "dtrexc.f"
	    nbnext = 1;
#line 262 "dtrexc.f"
	    if (here + nbf + 1 <= *n) {
#line 263 "dtrexc.f"
		if (t[here + nbf + 1 + (here + nbf) * t_dim1] != 0.) {
#line 263 "dtrexc.f"
		    nbnext = 2;
#line 263 "dtrexc.f"
		}
#line 265 "dtrexc.f"
	    }
#line 266 "dtrexc.f"
	    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &here, &
		    nbf, &nbnext, &work[1], info);
#line 268 "dtrexc.f"
	    if (*info != 0) {
#line 269 "dtrexc.f"
		*ilst = here;
#line 270 "dtrexc.f"
		return 0;
#line 271 "dtrexc.f"
	    }
#line 272 "dtrexc.f"
	    here += nbnext;

/*           Test if 2 by 2 block breaks into two 1 by 1 blocks */

#line 276 "dtrexc.f"
	    if (nbf == 2) {
#line 277 "dtrexc.f"
		if (t[here + 1 + here * t_dim1] == 0.) {
#line 277 "dtrexc.f"
		    nbf = 3;
#line 277 "dtrexc.f"
		}
#line 279 "dtrexc.f"
	    }

#line 281 "dtrexc.f"
	} else {

/*           Current block consists of two 1 by 1 blocks each of which */
/*           must be swapped individually */

#line 286 "dtrexc.f"
	    nbnext = 1;
#line 287 "dtrexc.f"
	    if (here + 3 <= *n) {
#line 288 "dtrexc.f"
		if (t[here + 3 + (here + 2) * t_dim1] != 0.) {
#line 288 "dtrexc.f"
		    nbnext = 2;
#line 288 "dtrexc.f"
		}
#line 290 "dtrexc.f"
	    }
#line 291 "dtrexc.f"
	    i__1 = here + 1;
#line 291 "dtrexc.f"
	    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &
		    c__1, &nbnext, &work[1], info);
#line 293 "dtrexc.f"
	    if (*info != 0) {
#line 294 "dtrexc.f"
		*ilst = here;
#line 295 "dtrexc.f"
		return 0;
#line 296 "dtrexc.f"
	    }
#line 297 "dtrexc.f"
	    if (nbnext == 1) {

/*              Swap two 1 by 1 blocks, no problems possible */

#line 301 "dtrexc.f"
		dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			here, &c__1, &nbnext, &work[1], info);
#line 303 "dtrexc.f"
		++here;
#line 304 "dtrexc.f"
	    } else {

/*              Recompute NBNEXT in case 2 by 2 split */

#line 308 "dtrexc.f"
		if (t[here + 2 + (here + 1) * t_dim1] == 0.) {
#line 308 "dtrexc.f"
		    nbnext = 1;
#line 308 "dtrexc.f"
		}
#line 310 "dtrexc.f"
		if (nbnext == 2) {

/*                 2 by 2 Block did not split */

#line 314 "dtrexc.f"
		    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			    here, &c__1, &nbnext, &work[1], info);
#line 316 "dtrexc.f"
		    if (*info != 0) {
#line 317 "dtrexc.f"
			*ilst = here;
#line 318 "dtrexc.f"
			return 0;
#line 319 "dtrexc.f"
		    }
#line 320 "dtrexc.f"
		    here += 2;
#line 321 "dtrexc.f"
		} else {

/*                 2 by 2 Block did split */

#line 325 "dtrexc.f"
		    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			    here, &c__1, &c__1, &work[1], info);
#line 327 "dtrexc.f"
		    i__1 = here + 1;
#line 327 "dtrexc.f"
		    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			    i__1, &c__1, &c__1, &work[1], info);
#line 329 "dtrexc.f"
		    here += 2;
#line 330 "dtrexc.f"
		}
#line 331 "dtrexc.f"
	    }
#line 332 "dtrexc.f"
	}
#line 333 "dtrexc.f"
	if (here < *ilst) {
#line 333 "dtrexc.f"
	    goto L10;
#line 333 "dtrexc.f"
	}

#line 336 "dtrexc.f"
    } else {

#line 338 "dtrexc.f"
	here = *ifst;
#line 339 "dtrexc.f"
L20:

/*        Swap block with next one above */

#line 343 "dtrexc.f"
	if (nbf == 1 || nbf == 2) {

/*           Current block either 1 by 1 or 2 by 2 */

#line 347 "dtrexc.f"
	    nbnext = 1;
#line 348 "dtrexc.f"
	    if (here >= 3) {
#line 349 "dtrexc.f"
		if (t[here - 1 + (here - 2) * t_dim1] != 0.) {
#line 349 "dtrexc.f"
		    nbnext = 2;
#line 349 "dtrexc.f"
		}
#line 351 "dtrexc.f"
	    }
#line 352 "dtrexc.f"
	    i__1 = here - nbnext;
#line 352 "dtrexc.f"
	    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &
		    nbnext, &nbf, &work[1], info);
#line 354 "dtrexc.f"
	    if (*info != 0) {
#line 355 "dtrexc.f"
		*ilst = here;
#line 356 "dtrexc.f"
		return 0;
#line 357 "dtrexc.f"
	    }
#line 358 "dtrexc.f"
	    here -= nbnext;

/*           Test if 2 by 2 block breaks into two 1 by 1 blocks */

#line 362 "dtrexc.f"
	    if (nbf == 2) {
#line 363 "dtrexc.f"
		if (t[here + 1 + here * t_dim1] == 0.) {
#line 363 "dtrexc.f"
		    nbf = 3;
#line 363 "dtrexc.f"
		}
#line 365 "dtrexc.f"
	    }

#line 367 "dtrexc.f"
	} else {

/*           Current block consists of two 1 by 1 blocks each of which */
/*           must be swapped individually */

#line 372 "dtrexc.f"
	    nbnext = 1;
#line 373 "dtrexc.f"
	    if (here >= 3) {
#line 374 "dtrexc.f"
		if (t[here - 1 + (here - 2) * t_dim1] != 0.) {
#line 374 "dtrexc.f"
		    nbnext = 2;
#line 374 "dtrexc.f"
		}
#line 376 "dtrexc.f"
	    }
#line 377 "dtrexc.f"
	    i__1 = here - nbnext;
#line 377 "dtrexc.f"
	    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &
		    nbnext, &c__1, &work[1], info);
#line 379 "dtrexc.f"
	    if (*info != 0) {
#line 380 "dtrexc.f"
		*ilst = here;
#line 381 "dtrexc.f"
		return 0;
#line 382 "dtrexc.f"
	    }
#line 383 "dtrexc.f"
	    if (nbnext == 1) {

/*              Swap two 1 by 1 blocks, no problems possible */

#line 387 "dtrexc.f"
		dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			here, &nbnext, &c__1, &work[1], info);
#line 389 "dtrexc.f"
		--here;
#line 390 "dtrexc.f"
	    } else {

/*              Recompute NBNEXT in case 2 by 2 split */

#line 394 "dtrexc.f"
		if (t[here + (here - 1) * t_dim1] == 0.) {
#line 394 "dtrexc.f"
		    nbnext = 1;
#line 394 "dtrexc.f"
		}
#line 396 "dtrexc.f"
		if (nbnext == 2) {

/*                 2 by 2 Block did not split */

#line 400 "dtrexc.f"
		    i__1 = here - 1;
#line 400 "dtrexc.f"
		    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			    i__1, &c__2, &c__1, &work[1], info);
#line 402 "dtrexc.f"
		    if (*info != 0) {
#line 403 "dtrexc.f"
			*ilst = here;
#line 404 "dtrexc.f"
			return 0;
#line 405 "dtrexc.f"
		    }
#line 406 "dtrexc.f"
		    here += -2;
#line 407 "dtrexc.f"
		} else {

/*                 2 by 2 Block did split */

#line 411 "dtrexc.f"
		    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			    here, &c__1, &c__1, &work[1], info);
#line 413 "dtrexc.f"
		    i__1 = here - 1;
#line 413 "dtrexc.f"
		    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			    i__1, &c__1, &c__1, &work[1], info);
#line 415 "dtrexc.f"
		    here += -2;
#line 416 "dtrexc.f"
		}
#line 417 "dtrexc.f"
	    }
#line 418 "dtrexc.f"
	}
#line 419 "dtrexc.f"
	if (here > *ilst) {
#line 419 "dtrexc.f"
	    goto L20;
#line 419 "dtrexc.f"
	}
#line 421 "dtrexc.f"
    }
#line 422 "dtrexc.f"
    *ilst = here;

#line 424 "dtrexc.f"
    return 0;

/*     End of DTREXC */

} /* dtrexc_ */

