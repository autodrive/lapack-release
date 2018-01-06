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
/* >          The leading dimension of the array Q.  LDQ >= max(1,N). */
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

/* > \date November 2011 */

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


/*  -- LAPACK computational routine (version 3.4.0) -- */
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

/*     Decode and test the input arguments. */

#line 186 "dtrexc.f"
    /* Parameter adjustments */
#line 186 "dtrexc.f"
    t_dim1 = *ldt;
#line 186 "dtrexc.f"
    t_offset = 1 + t_dim1;
#line 186 "dtrexc.f"
    t -= t_offset;
#line 186 "dtrexc.f"
    q_dim1 = *ldq;
#line 186 "dtrexc.f"
    q_offset = 1 + q_dim1;
#line 186 "dtrexc.f"
    q -= q_offset;
#line 186 "dtrexc.f"
    --work;
#line 186 "dtrexc.f"

#line 186 "dtrexc.f"
    /* Function Body */
#line 186 "dtrexc.f"
    *info = 0;
#line 187 "dtrexc.f"
    wantq = lsame_(compq, "V", (ftnlen)1, (ftnlen)1);
#line 188 "dtrexc.f"
    if (! wantq && ! lsame_(compq, "N", (ftnlen)1, (ftnlen)1)) {
#line 189 "dtrexc.f"
	*info = -1;
#line 190 "dtrexc.f"
    } else if (*n < 0) {
#line 191 "dtrexc.f"
	*info = -2;
#line 192 "dtrexc.f"
    } else if (*ldt < max(1,*n)) {
#line 193 "dtrexc.f"
	*info = -4;
#line 194 "dtrexc.f"
    } else if (*ldq < 1 || wantq && *ldq < max(1,*n)) {
#line 195 "dtrexc.f"
	*info = -6;
#line 196 "dtrexc.f"
    } else if (*ifst < 1 || *ifst > *n) {
#line 197 "dtrexc.f"
	*info = -7;
#line 198 "dtrexc.f"
    } else if (*ilst < 1 || *ilst > *n) {
#line 199 "dtrexc.f"
	*info = -8;
#line 200 "dtrexc.f"
    }
#line 201 "dtrexc.f"
    if (*info != 0) {
#line 202 "dtrexc.f"
	i__1 = -(*info);
#line 202 "dtrexc.f"
	xerbla_("DTREXC", &i__1, (ftnlen)6);
#line 203 "dtrexc.f"
	return 0;
#line 204 "dtrexc.f"
    }

/*     Quick return if possible */

#line 208 "dtrexc.f"
    if (*n <= 1) {
#line 208 "dtrexc.f"
	return 0;
#line 208 "dtrexc.f"
    }

/*     Determine the first row of specified block */
/*     and find out it is 1 by 1 or 2 by 2. */

#line 214 "dtrexc.f"
    if (*ifst > 1) {
#line 215 "dtrexc.f"
	if (t[*ifst + (*ifst - 1) * t_dim1] != 0.) {
#line 215 "dtrexc.f"
	    --(*ifst);
#line 215 "dtrexc.f"
	}
#line 217 "dtrexc.f"
    }
#line 218 "dtrexc.f"
    nbf = 1;
#line 219 "dtrexc.f"
    if (*ifst < *n) {
#line 220 "dtrexc.f"
	if (t[*ifst + 1 + *ifst * t_dim1] != 0.) {
#line 220 "dtrexc.f"
	    nbf = 2;
#line 220 "dtrexc.f"
	}
#line 222 "dtrexc.f"
    }

/*     Determine the first row of the final block */
/*     and find out it is 1 by 1 or 2 by 2. */

#line 227 "dtrexc.f"
    if (*ilst > 1) {
#line 228 "dtrexc.f"
	if (t[*ilst + (*ilst - 1) * t_dim1] != 0.) {
#line 228 "dtrexc.f"
	    --(*ilst);
#line 228 "dtrexc.f"
	}
#line 230 "dtrexc.f"
    }
#line 231 "dtrexc.f"
    nbl = 1;
#line 232 "dtrexc.f"
    if (*ilst < *n) {
#line 233 "dtrexc.f"
	if (t[*ilst + 1 + *ilst * t_dim1] != 0.) {
#line 233 "dtrexc.f"
	    nbl = 2;
#line 233 "dtrexc.f"
	}
#line 235 "dtrexc.f"
    }

#line 237 "dtrexc.f"
    if (*ifst == *ilst) {
#line 237 "dtrexc.f"
	return 0;
#line 237 "dtrexc.f"
    }

#line 240 "dtrexc.f"
    if (*ifst < *ilst) {

/*        Update ILST */

#line 244 "dtrexc.f"
	if (nbf == 2 && nbl == 1) {
#line 244 "dtrexc.f"
	    --(*ilst);
#line 244 "dtrexc.f"
	}
#line 246 "dtrexc.f"
	if (nbf == 1 && nbl == 2) {
#line 246 "dtrexc.f"
	    ++(*ilst);
#line 246 "dtrexc.f"
	}

#line 249 "dtrexc.f"
	here = *ifst;

#line 251 "dtrexc.f"
L10:

/*        Swap block with next one below */

#line 255 "dtrexc.f"
	if (nbf == 1 || nbf == 2) {

/*           Current block either 1 by 1 or 2 by 2 */

#line 259 "dtrexc.f"
	    nbnext = 1;
#line 260 "dtrexc.f"
	    if (here + nbf + 1 <= *n) {
#line 261 "dtrexc.f"
		if (t[here + nbf + 1 + (here + nbf) * t_dim1] != 0.) {
#line 261 "dtrexc.f"
		    nbnext = 2;
#line 261 "dtrexc.f"
		}
#line 263 "dtrexc.f"
	    }
#line 264 "dtrexc.f"
	    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &here, &
		    nbf, &nbnext, &work[1], info);
#line 266 "dtrexc.f"
	    if (*info != 0) {
#line 267 "dtrexc.f"
		*ilst = here;
#line 268 "dtrexc.f"
		return 0;
#line 269 "dtrexc.f"
	    }
#line 270 "dtrexc.f"
	    here += nbnext;

/*           Test if 2 by 2 block breaks into two 1 by 1 blocks */

#line 274 "dtrexc.f"
	    if (nbf == 2) {
#line 275 "dtrexc.f"
		if (t[here + 1 + here * t_dim1] == 0.) {
#line 275 "dtrexc.f"
		    nbf = 3;
#line 275 "dtrexc.f"
		}
#line 277 "dtrexc.f"
	    }

#line 279 "dtrexc.f"
	} else {

/*           Current block consists of two 1 by 1 blocks each of which */
/*           must be swapped individually */

#line 284 "dtrexc.f"
	    nbnext = 1;
#line 285 "dtrexc.f"
	    if (here + 3 <= *n) {
#line 286 "dtrexc.f"
		if (t[here + 3 + (here + 2) * t_dim1] != 0.) {
#line 286 "dtrexc.f"
		    nbnext = 2;
#line 286 "dtrexc.f"
		}
#line 288 "dtrexc.f"
	    }
#line 289 "dtrexc.f"
	    i__1 = here + 1;
#line 289 "dtrexc.f"
	    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &
		    c__1, &nbnext, &work[1], info);
#line 291 "dtrexc.f"
	    if (*info != 0) {
#line 292 "dtrexc.f"
		*ilst = here;
#line 293 "dtrexc.f"
		return 0;
#line 294 "dtrexc.f"
	    }
#line 295 "dtrexc.f"
	    if (nbnext == 1) {

/*              Swap two 1 by 1 blocks, no problems possible */

#line 299 "dtrexc.f"
		dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			here, &c__1, &nbnext, &work[1], info);
#line 301 "dtrexc.f"
		++here;
#line 302 "dtrexc.f"
	    } else {

/*              Recompute NBNEXT in case 2 by 2 split */

#line 306 "dtrexc.f"
		if (t[here + 2 + (here + 1) * t_dim1] == 0.) {
#line 306 "dtrexc.f"
		    nbnext = 1;
#line 306 "dtrexc.f"
		}
#line 308 "dtrexc.f"
		if (nbnext == 2) {

/*                 2 by 2 Block did not split */

#line 312 "dtrexc.f"
		    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			    here, &c__1, &nbnext, &work[1], info);
#line 314 "dtrexc.f"
		    if (*info != 0) {
#line 315 "dtrexc.f"
			*ilst = here;
#line 316 "dtrexc.f"
			return 0;
#line 317 "dtrexc.f"
		    }
#line 318 "dtrexc.f"
		    here += 2;
#line 319 "dtrexc.f"
		} else {

/*                 2 by 2 Block did split */

#line 323 "dtrexc.f"
		    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			    here, &c__1, &c__1, &work[1], info);
#line 325 "dtrexc.f"
		    i__1 = here + 1;
#line 325 "dtrexc.f"
		    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			    i__1, &c__1, &c__1, &work[1], info);
#line 327 "dtrexc.f"
		    here += 2;
#line 328 "dtrexc.f"
		}
#line 329 "dtrexc.f"
	    }
#line 330 "dtrexc.f"
	}
#line 331 "dtrexc.f"
	if (here < *ilst) {
#line 331 "dtrexc.f"
	    goto L10;
#line 331 "dtrexc.f"
	}

#line 334 "dtrexc.f"
    } else {

#line 336 "dtrexc.f"
	here = *ifst;
#line 337 "dtrexc.f"
L20:

/*        Swap block with next one above */

#line 341 "dtrexc.f"
	if (nbf == 1 || nbf == 2) {

/*           Current block either 1 by 1 or 2 by 2 */

#line 345 "dtrexc.f"
	    nbnext = 1;
#line 346 "dtrexc.f"
	    if (here >= 3) {
#line 347 "dtrexc.f"
		if (t[here - 1 + (here - 2) * t_dim1] != 0.) {
#line 347 "dtrexc.f"
		    nbnext = 2;
#line 347 "dtrexc.f"
		}
#line 349 "dtrexc.f"
	    }
#line 350 "dtrexc.f"
	    i__1 = here - nbnext;
#line 350 "dtrexc.f"
	    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &
		    nbnext, &nbf, &work[1], info);
#line 352 "dtrexc.f"
	    if (*info != 0) {
#line 353 "dtrexc.f"
		*ilst = here;
#line 354 "dtrexc.f"
		return 0;
#line 355 "dtrexc.f"
	    }
#line 356 "dtrexc.f"
	    here -= nbnext;

/*           Test if 2 by 2 block breaks into two 1 by 1 blocks */

#line 360 "dtrexc.f"
	    if (nbf == 2) {
#line 361 "dtrexc.f"
		if (t[here + 1 + here * t_dim1] == 0.) {
#line 361 "dtrexc.f"
		    nbf = 3;
#line 361 "dtrexc.f"
		}
#line 363 "dtrexc.f"
	    }

#line 365 "dtrexc.f"
	} else {

/*           Current block consists of two 1 by 1 blocks each of which */
/*           must be swapped individually */

#line 370 "dtrexc.f"
	    nbnext = 1;
#line 371 "dtrexc.f"
	    if (here >= 3) {
#line 372 "dtrexc.f"
		if (t[here - 1 + (here - 2) * t_dim1] != 0.) {
#line 372 "dtrexc.f"
		    nbnext = 2;
#line 372 "dtrexc.f"
		}
#line 374 "dtrexc.f"
	    }
#line 375 "dtrexc.f"
	    i__1 = here - nbnext;
#line 375 "dtrexc.f"
	    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &
		    nbnext, &c__1, &work[1], info);
#line 377 "dtrexc.f"
	    if (*info != 0) {
#line 378 "dtrexc.f"
		*ilst = here;
#line 379 "dtrexc.f"
		return 0;
#line 380 "dtrexc.f"
	    }
#line 381 "dtrexc.f"
	    if (nbnext == 1) {

/*              Swap two 1 by 1 blocks, no problems possible */

#line 385 "dtrexc.f"
		dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			here, &nbnext, &c__1, &work[1], info);
#line 387 "dtrexc.f"
		--here;
#line 388 "dtrexc.f"
	    } else {

/*              Recompute NBNEXT in case 2 by 2 split */

#line 392 "dtrexc.f"
		if (t[here + (here - 1) * t_dim1] == 0.) {
#line 392 "dtrexc.f"
		    nbnext = 1;
#line 392 "dtrexc.f"
		}
#line 394 "dtrexc.f"
		if (nbnext == 2) {

/*                 2 by 2 Block did not split */

#line 398 "dtrexc.f"
		    i__1 = here - 1;
#line 398 "dtrexc.f"
		    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			    i__1, &c__2, &c__1, &work[1], info);
#line 400 "dtrexc.f"
		    if (*info != 0) {
#line 401 "dtrexc.f"
			*ilst = here;
#line 402 "dtrexc.f"
			return 0;
#line 403 "dtrexc.f"
		    }
#line 404 "dtrexc.f"
		    here += -2;
#line 405 "dtrexc.f"
		} else {

/*                 2 by 2 Block did split */

#line 409 "dtrexc.f"
		    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			    here, &c__1, &c__1, &work[1], info);
#line 411 "dtrexc.f"
		    i__1 = here - 1;
#line 411 "dtrexc.f"
		    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			    i__1, &c__1, &c__1, &work[1], info);
#line 413 "dtrexc.f"
		    here += -2;
#line 414 "dtrexc.f"
		}
#line 415 "dtrexc.f"
	    }
#line 416 "dtrexc.f"
	}
#line 417 "dtrexc.f"
	if (here > *ilst) {
#line 417 "dtrexc.f"
	    goto L20;
#line 417 "dtrexc.f"
	}
#line 419 "dtrexc.f"
    }
#line 420 "dtrexc.f"
    *ilst = here;

#line 422 "dtrexc.f"
    return 0;

/*     End of DTREXC */

} /* dtrexc_ */

