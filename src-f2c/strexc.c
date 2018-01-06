#line 1 "strexc.f"
/* strexc.f -- translated by f2c (version 20100827).
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

#line 1 "strexc.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b STREXC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STREXC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strexc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strexc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strexc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPQ */
/*       INTEGER            IFST, ILST, INFO, LDQ, LDT, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               Q( LDQ, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STREXC reorders the real Schur factorization of a real matrix */
/* > A = Q*T*Q**T, so that the diagonal block of T with row index IFST is */
/* > moved to row ILST. */
/* > */
/* > The real Schur form T is reordered by an orthogonal similarity */
/* > transformation Z**T*T*Z, and optionally the matrix Q of Schur vectors */
/* > is updated by postmultiplying it with Z. */
/* > */
/* > T must be in Schur canonical form (as returned by SHSEQR), that is, */
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
/* >          T is REAL array, dimension (LDT,N) */
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
/* >          Q is REAL array, dimension (LDQ,N) */
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
/* >          WORK is REAL array, dimension (N) */
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

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int strexc_(char *compq, integer *n, doublereal *t, integer *
	ldt, doublereal *q, integer *ldq, integer *ifst, integer *ilst, 
	doublereal *work, integer *info, ftnlen compq_len)
{
    /* System generated locals */
    integer q_dim1, q_offset, t_dim1, t_offset, i__1;

    /* Local variables */
    static integer nbf, nbl, here;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical wantq;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), slaexc_(
	    logical *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *);
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

#line 188 "strexc.f"
    /* Parameter adjustments */
#line 188 "strexc.f"
    t_dim1 = *ldt;
#line 188 "strexc.f"
    t_offset = 1 + t_dim1;
#line 188 "strexc.f"
    t -= t_offset;
#line 188 "strexc.f"
    q_dim1 = *ldq;
#line 188 "strexc.f"
    q_offset = 1 + q_dim1;
#line 188 "strexc.f"
    q -= q_offset;
#line 188 "strexc.f"
    --work;
#line 188 "strexc.f"

#line 188 "strexc.f"
    /* Function Body */
#line 188 "strexc.f"
    *info = 0;
#line 189 "strexc.f"
    wantq = lsame_(compq, "V", (ftnlen)1, (ftnlen)1);
#line 190 "strexc.f"
    if (! wantq && ! lsame_(compq, "N", (ftnlen)1, (ftnlen)1)) {
#line 191 "strexc.f"
	*info = -1;
#line 192 "strexc.f"
    } else if (*n < 0) {
#line 193 "strexc.f"
	*info = -2;
#line 194 "strexc.f"
    } else if (*ldt < max(1,*n)) {
#line 195 "strexc.f"
	*info = -4;
#line 196 "strexc.f"
    } else if (*ldq < 1 || wantq && *ldq < max(1,*n)) {
#line 197 "strexc.f"
	*info = -6;
#line 198 "strexc.f"
    } else if ((*ifst < 1 || *ifst > *n) && *n > 0) {
#line 199 "strexc.f"
	*info = -7;
#line 200 "strexc.f"
    } else if ((*ilst < 1 || *ilst > *n) && *n > 0) {
#line 201 "strexc.f"
	*info = -8;
#line 202 "strexc.f"
    }
#line 203 "strexc.f"
    if (*info != 0) {
#line 204 "strexc.f"
	i__1 = -(*info);
#line 204 "strexc.f"
	xerbla_("STREXC", &i__1, (ftnlen)6);
#line 205 "strexc.f"
	return 0;
#line 206 "strexc.f"
    }

/*     Quick return if possible */

#line 210 "strexc.f"
    if (*n <= 1) {
#line 210 "strexc.f"
	return 0;
#line 210 "strexc.f"
    }

/*     Determine the first row of specified block */
/*     and find out it is 1 by 1 or 2 by 2. */

#line 216 "strexc.f"
    if (*ifst > 1) {
#line 217 "strexc.f"
	if (t[*ifst + (*ifst - 1) * t_dim1] != 0.) {
#line 217 "strexc.f"
	    --(*ifst);
#line 217 "strexc.f"
	}
#line 219 "strexc.f"
    }
#line 220 "strexc.f"
    nbf = 1;
#line 221 "strexc.f"
    if (*ifst < *n) {
#line 222 "strexc.f"
	if (t[*ifst + 1 + *ifst * t_dim1] != 0.) {
#line 222 "strexc.f"
	    nbf = 2;
#line 222 "strexc.f"
	}
#line 224 "strexc.f"
    }

/*     Determine the first row of the final block */
/*     and find out it is 1 by 1 or 2 by 2. */

#line 229 "strexc.f"
    if (*ilst > 1) {
#line 230 "strexc.f"
	if (t[*ilst + (*ilst - 1) * t_dim1] != 0.) {
#line 230 "strexc.f"
	    --(*ilst);
#line 230 "strexc.f"
	}
#line 232 "strexc.f"
    }
#line 233 "strexc.f"
    nbl = 1;
#line 234 "strexc.f"
    if (*ilst < *n) {
#line 235 "strexc.f"
	if (t[*ilst + 1 + *ilst * t_dim1] != 0.) {
#line 235 "strexc.f"
	    nbl = 2;
#line 235 "strexc.f"
	}
#line 237 "strexc.f"
    }

#line 239 "strexc.f"
    if (*ifst == *ilst) {
#line 239 "strexc.f"
	return 0;
#line 239 "strexc.f"
    }

#line 242 "strexc.f"
    if (*ifst < *ilst) {

/*        Update ILST */

#line 246 "strexc.f"
	if (nbf == 2 && nbl == 1) {
#line 246 "strexc.f"
	    --(*ilst);
#line 246 "strexc.f"
	}
#line 248 "strexc.f"
	if (nbf == 1 && nbl == 2) {
#line 248 "strexc.f"
	    ++(*ilst);
#line 248 "strexc.f"
	}

#line 251 "strexc.f"
	here = *ifst;

#line 253 "strexc.f"
L10:

/*        Swap block with next one below */

#line 257 "strexc.f"
	if (nbf == 1 || nbf == 2) {

/*           Current block either 1 by 1 or 2 by 2 */

#line 261 "strexc.f"
	    nbnext = 1;
#line 262 "strexc.f"
	    if (here + nbf + 1 <= *n) {
#line 263 "strexc.f"
		if (t[here + nbf + 1 + (here + nbf) * t_dim1] != 0.) {
#line 263 "strexc.f"
		    nbnext = 2;
#line 263 "strexc.f"
		}
#line 265 "strexc.f"
	    }
#line 266 "strexc.f"
	    slaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &here, &
		    nbf, &nbnext, &work[1], info);
#line 268 "strexc.f"
	    if (*info != 0) {
#line 269 "strexc.f"
		*ilst = here;
#line 270 "strexc.f"
		return 0;
#line 271 "strexc.f"
	    }
#line 272 "strexc.f"
	    here += nbnext;

/*           Test if 2 by 2 block breaks into two 1 by 1 blocks */

#line 276 "strexc.f"
	    if (nbf == 2) {
#line 277 "strexc.f"
		if (t[here + 1 + here * t_dim1] == 0.) {
#line 277 "strexc.f"
		    nbf = 3;
#line 277 "strexc.f"
		}
#line 279 "strexc.f"
	    }

#line 281 "strexc.f"
	} else {

/*           Current block consists of two 1 by 1 blocks each of which */
/*           must be swapped individually */

#line 286 "strexc.f"
	    nbnext = 1;
#line 287 "strexc.f"
	    if (here + 3 <= *n) {
#line 288 "strexc.f"
		if (t[here + 3 + (here + 2) * t_dim1] != 0.) {
#line 288 "strexc.f"
		    nbnext = 2;
#line 288 "strexc.f"
		}
#line 290 "strexc.f"
	    }
#line 291 "strexc.f"
	    i__1 = here + 1;
#line 291 "strexc.f"
	    slaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &
		    c__1, &nbnext, &work[1], info);
#line 293 "strexc.f"
	    if (*info != 0) {
#line 294 "strexc.f"
		*ilst = here;
#line 295 "strexc.f"
		return 0;
#line 296 "strexc.f"
	    }
#line 297 "strexc.f"
	    if (nbnext == 1) {

/*              Swap two 1 by 1 blocks, no problems possible */

#line 301 "strexc.f"
		slaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			here, &c__1, &nbnext, &work[1], info);
#line 303 "strexc.f"
		++here;
#line 304 "strexc.f"
	    } else {

/*              Recompute NBNEXT in case 2 by 2 split */

#line 308 "strexc.f"
		if (t[here + 2 + (here + 1) * t_dim1] == 0.) {
#line 308 "strexc.f"
		    nbnext = 1;
#line 308 "strexc.f"
		}
#line 310 "strexc.f"
		if (nbnext == 2) {

/*                 2 by 2 Block did not split */

#line 314 "strexc.f"
		    slaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			    here, &c__1, &nbnext, &work[1], info);
#line 316 "strexc.f"
		    if (*info != 0) {
#line 317 "strexc.f"
			*ilst = here;
#line 318 "strexc.f"
			return 0;
#line 319 "strexc.f"
		    }
#line 320 "strexc.f"
		    here += 2;
#line 321 "strexc.f"
		} else {

/*                 2 by 2 Block did split */

#line 325 "strexc.f"
		    slaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			    here, &c__1, &c__1, &work[1], info);
#line 327 "strexc.f"
		    i__1 = here + 1;
#line 327 "strexc.f"
		    slaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			    i__1, &c__1, &c__1, &work[1], info);
#line 329 "strexc.f"
		    here += 2;
#line 330 "strexc.f"
		}
#line 331 "strexc.f"
	    }
#line 332 "strexc.f"
	}
#line 333 "strexc.f"
	if (here < *ilst) {
#line 333 "strexc.f"
	    goto L10;
#line 333 "strexc.f"
	}

#line 336 "strexc.f"
    } else {

#line 338 "strexc.f"
	here = *ifst;
#line 339 "strexc.f"
L20:

/*        Swap block with next one above */

#line 343 "strexc.f"
	if (nbf == 1 || nbf == 2) {

/*           Current block either 1 by 1 or 2 by 2 */

#line 347 "strexc.f"
	    nbnext = 1;
#line 348 "strexc.f"
	    if (here >= 3) {
#line 349 "strexc.f"
		if (t[here - 1 + (here - 2) * t_dim1] != 0.) {
#line 349 "strexc.f"
		    nbnext = 2;
#line 349 "strexc.f"
		}
#line 351 "strexc.f"
	    }
#line 352 "strexc.f"
	    i__1 = here - nbnext;
#line 352 "strexc.f"
	    slaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &
		    nbnext, &nbf, &work[1], info);
#line 354 "strexc.f"
	    if (*info != 0) {
#line 355 "strexc.f"
		*ilst = here;
#line 356 "strexc.f"
		return 0;
#line 357 "strexc.f"
	    }
#line 358 "strexc.f"
	    here -= nbnext;

/*           Test if 2 by 2 block breaks into two 1 by 1 blocks */

#line 362 "strexc.f"
	    if (nbf == 2) {
#line 363 "strexc.f"
		if (t[here + 1 + here * t_dim1] == 0.) {
#line 363 "strexc.f"
		    nbf = 3;
#line 363 "strexc.f"
		}
#line 365 "strexc.f"
	    }

#line 367 "strexc.f"
	} else {

/*           Current block consists of two 1 by 1 blocks each of which */
/*           must be swapped individually */

#line 372 "strexc.f"
	    nbnext = 1;
#line 373 "strexc.f"
	    if (here >= 3) {
#line 374 "strexc.f"
		if (t[here - 1 + (here - 2) * t_dim1] != 0.) {
#line 374 "strexc.f"
		    nbnext = 2;
#line 374 "strexc.f"
		}
#line 376 "strexc.f"
	    }
#line 377 "strexc.f"
	    i__1 = here - nbnext;
#line 377 "strexc.f"
	    slaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &
		    nbnext, &c__1, &work[1], info);
#line 379 "strexc.f"
	    if (*info != 0) {
#line 380 "strexc.f"
		*ilst = here;
#line 381 "strexc.f"
		return 0;
#line 382 "strexc.f"
	    }
#line 383 "strexc.f"
	    if (nbnext == 1) {

/*              Swap two 1 by 1 blocks, no problems possible */

#line 387 "strexc.f"
		slaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			here, &nbnext, &c__1, &work[1], info);
#line 389 "strexc.f"
		--here;
#line 390 "strexc.f"
	    } else {

/*              Recompute NBNEXT in case 2 by 2 split */

#line 394 "strexc.f"
		if (t[here + (here - 1) * t_dim1] == 0.) {
#line 394 "strexc.f"
		    nbnext = 1;
#line 394 "strexc.f"
		}
#line 396 "strexc.f"
		if (nbnext == 2) {

/*                 2 by 2 Block did not split */

#line 400 "strexc.f"
		    i__1 = here - 1;
#line 400 "strexc.f"
		    slaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			    i__1, &c__2, &c__1, &work[1], info);
#line 402 "strexc.f"
		    if (*info != 0) {
#line 403 "strexc.f"
			*ilst = here;
#line 404 "strexc.f"
			return 0;
#line 405 "strexc.f"
		    }
#line 406 "strexc.f"
		    here += -2;
#line 407 "strexc.f"
		} else {

/*                 2 by 2 Block did split */

#line 411 "strexc.f"
		    slaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			    here, &c__1, &c__1, &work[1], info);
#line 413 "strexc.f"
		    i__1 = here - 1;
#line 413 "strexc.f"
		    slaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			    i__1, &c__1, &c__1, &work[1], info);
#line 415 "strexc.f"
		    here += -2;
#line 416 "strexc.f"
		}
#line 417 "strexc.f"
	    }
#line 418 "strexc.f"
	}
#line 419 "strexc.f"
	if (here > *ilst) {
#line 419 "strexc.f"
	    goto L20;
#line 419 "strexc.f"
	}
#line 421 "strexc.f"
    }
#line 422 "strexc.f"
    *ilst = here;

#line 424 "strexc.f"
    return 0;

/*     End of STREXC */

} /* strexc_ */

