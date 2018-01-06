#line 1 "slascl.f"
/* slascl.f -- translated by f2c (version 20100827).
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

#line 1 "slascl.f"
/* > \brief \b SLASCL multiplies a general rectangular matrix by a real scalar defined as cto/cfrom. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLASCL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slascl.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slascl.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slascl.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TYPE */
/*       INTEGER            INFO, KL, KU, LDA, M, N */
/*       REAL               CFROM, CTO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASCL multiplies the M by N real matrix A by the real scalar */
/* > CTO/CFROM.  This is done without over/underflow as long as the final */
/* > result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that */
/* > A may be full, upper triangular, lower triangular, upper Hessenberg, */
/* > or banded. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TYPE */
/* > \verbatim */
/* >          TYPE is CHARACTER*1 */
/* >          TYPE indices the storage type of the input matrix. */
/* >          = 'G':  A is a full matrix. */
/* >          = 'L':  A is a lower triangular matrix. */
/* >          = 'U':  A is an upper triangular matrix. */
/* >          = 'H':  A is an upper Hessenberg matrix. */
/* >          = 'B':  A is a symmetric band matrix with lower bandwidth KL */
/* >                  and upper bandwidth KU and with the only the lower */
/* >                  half stored. */
/* >          = 'Q':  A is a symmetric band matrix with lower bandwidth KL */
/* >                  and upper bandwidth KU and with the only the upper */
/* >                  half stored. */
/* >          = 'Z':  A is a band matrix with lower bandwidth KL and upper */
/* >                  bandwidth KU. See SGBTRF for storage details. */
/* > \endverbatim */
/* > */
/* > \param[in] KL */
/* > \verbatim */
/* >          KL is INTEGER */
/* >          The lower bandwidth of A.  Referenced only if TYPE = 'B', */
/* >          'Q' or 'Z'. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* >          KU is INTEGER */
/* >          The upper bandwidth of A.  Referenced only if TYPE = 'B', */
/* >          'Q' or 'Z'. */
/* > \endverbatim */
/* > */
/* > \param[in] CFROM */
/* > \verbatim */
/* >          CFROM is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] CTO */
/* > \verbatim */
/* >          CTO is REAL */
/* > */
/* >          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed */
/* >          without over/underflow if the final result CTO*A(I,J)/CFROM */
/* >          can be represented without over/underflow.  CFROM must be */
/* >          nonzero. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          The matrix to be multiplied by CTO/CFROM.  See TYPE for the */
/* >          storage type. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. */
/* >          If TYPE = 'G', 'L', 'U', 'H', LDA >= max(1,M); */
/* >             TYPE = 'B', LDA >= KL+1; */
/* >             TYPE = 'Q', LDA >= KU+1; */
/* >             TYPE = 'Z', LDA >= 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          0  - successful exit */
/* >          <0 - if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slascl_(char *type__, integer *kl, integer *ku, 
	doublereal *cfrom, doublereal *cto, integer *m, integer *n, 
	doublereal *a, integer *lda, integer *info, ftnlen type_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, j, k1, k2, k3, k4;
    static doublereal mul, cto1;
    static logical done;
    static doublereal ctoc;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer itype;
    static doublereal cfrom1;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal cfromc;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern logical sisnan_(doublereal *);
    static doublereal smlnum;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 186 "slascl.f"
    /* Parameter adjustments */
#line 186 "slascl.f"
    a_dim1 = *lda;
#line 186 "slascl.f"
    a_offset = 1 + a_dim1;
#line 186 "slascl.f"
    a -= a_offset;
#line 186 "slascl.f"

#line 186 "slascl.f"
    /* Function Body */
#line 186 "slascl.f"
    *info = 0;

#line 188 "slascl.f"
    if (lsame_(type__, "G", (ftnlen)1, (ftnlen)1)) {
#line 189 "slascl.f"
	itype = 0;
#line 190 "slascl.f"
    } else if (lsame_(type__, "L", (ftnlen)1, (ftnlen)1)) {
#line 191 "slascl.f"
	itype = 1;
#line 192 "slascl.f"
    } else if (lsame_(type__, "U", (ftnlen)1, (ftnlen)1)) {
#line 193 "slascl.f"
	itype = 2;
#line 194 "slascl.f"
    } else if (lsame_(type__, "H", (ftnlen)1, (ftnlen)1)) {
#line 195 "slascl.f"
	itype = 3;
#line 196 "slascl.f"
    } else if (lsame_(type__, "B", (ftnlen)1, (ftnlen)1)) {
#line 197 "slascl.f"
	itype = 4;
#line 198 "slascl.f"
    } else if (lsame_(type__, "Q", (ftnlen)1, (ftnlen)1)) {
#line 199 "slascl.f"
	itype = 5;
#line 200 "slascl.f"
    } else if (lsame_(type__, "Z", (ftnlen)1, (ftnlen)1)) {
#line 201 "slascl.f"
	itype = 6;
#line 202 "slascl.f"
    } else {
#line 203 "slascl.f"
	itype = -1;
#line 204 "slascl.f"
    }

#line 206 "slascl.f"
    if (itype == -1) {
#line 207 "slascl.f"
	*info = -1;
#line 208 "slascl.f"
    } else if (*cfrom == 0. || sisnan_(cfrom)) {
#line 209 "slascl.f"
	*info = -4;
#line 210 "slascl.f"
    } else if (sisnan_(cto)) {
#line 211 "slascl.f"
	*info = -5;
#line 212 "slascl.f"
    } else if (*m < 0) {
#line 213 "slascl.f"
	*info = -6;
#line 214 "slascl.f"
    } else if (*n < 0 || itype == 4 && *n != *m || itype == 5 && *n != *m) {
#line 216 "slascl.f"
	*info = -7;
#line 217 "slascl.f"
    } else if (itype <= 3 && *lda < max(1,*m)) {
#line 218 "slascl.f"
	*info = -9;
#line 219 "slascl.f"
    } else if (itype >= 4) {
/* Computing MAX */
#line 220 "slascl.f"
	i__1 = *m - 1;
#line 220 "slascl.f"
	if (*kl < 0 || *kl > max(i__1,0)) {
#line 221 "slascl.f"
	    *info = -2;
#line 222 "slascl.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 222 "slascl.f"
	    i__1 = *n - 1;
#line 222 "slascl.f"
	    if (*ku < 0 || *ku > max(i__1,0) || (itype == 4 || itype == 5) && 
		    *kl != *ku) {
#line 225 "slascl.f"
		*info = -3;
#line 226 "slascl.f"
	    } else if (itype == 4 && *lda < *kl + 1 || itype == 5 && *lda < *
		    ku + 1 || itype == 6 && *lda < (*kl << 1) + *ku + 1) {
#line 229 "slascl.f"
		*info = -9;
#line 230 "slascl.f"
	    }
#line 230 "slascl.f"
	}
#line 231 "slascl.f"
    }

#line 233 "slascl.f"
    if (*info != 0) {
#line 234 "slascl.f"
	i__1 = -(*info);
#line 234 "slascl.f"
	xerbla_("SLASCL", &i__1, (ftnlen)6);
#line 235 "slascl.f"
	return 0;
#line 236 "slascl.f"
    }

/*     Quick return if possible */

#line 240 "slascl.f"
    if (*n == 0 || *m == 0) {
#line 240 "slascl.f"
	return 0;
#line 240 "slascl.f"
    }

/*     Get machine parameters */

#line 245 "slascl.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 246 "slascl.f"
    bignum = 1. / smlnum;

#line 248 "slascl.f"
    cfromc = *cfrom;
#line 249 "slascl.f"
    ctoc = *cto;

#line 251 "slascl.f"
L10:
#line 252 "slascl.f"
    cfrom1 = cfromc * smlnum;
#line 253 "slascl.f"
    if (cfrom1 == cfromc) {
/*        CFROMC is an inf.  Multiply by a correctly signed zero for */
/*        finite CTOC, or a NaN if CTOC is infinite. */
#line 256 "slascl.f"
	mul = ctoc / cfromc;
#line 257 "slascl.f"
	done = TRUE_;
#line 258 "slascl.f"
	cto1 = ctoc;
#line 259 "slascl.f"
    } else {
#line 260 "slascl.f"
	cto1 = ctoc / bignum;
#line 261 "slascl.f"
	if (cto1 == ctoc) {
/*           CTOC is either 0 or an inf.  In both cases, CTOC itself */
/*           serves as the correct multiplication factor. */
#line 264 "slascl.f"
	    mul = ctoc;
#line 265 "slascl.f"
	    done = TRUE_;
#line 266 "slascl.f"
	    cfromc = 1.;
#line 267 "slascl.f"
	} else if (abs(cfrom1) > abs(ctoc) && ctoc != 0.) {
#line 268 "slascl.f"
	    mul = smlnum;
#line 269 "slascl.f"
	    done = FALSE_;
#line 270 "slascl.f"
	    cfromc = cfrom1;
#line 271 "slascl.f"
	} else if (abs(cto1) > abs(cfromc)) {
#line 272 "slascl.f"
	    mul = bignum;
#line 273 "slascl.f"
	    done = FALSE_;
#line 274 "slascl.f"
	    ctoc = cto1;
#line 275 "slascl.f"
	} else {
#line 276 "slascl.f"
	    mul = ctoc / cfromc;
#line 277 "slascl.f"
	    done = TRUE_;
#line 278 "slascl.f"
	}
#line 279 "slascl.f"
    }

#line 281 "slascl.f"
    if (itype == 0) {

/*        Full matrix */

#line 285 "slascl.f"
	i__1 = *n;
#line 285 "slascl.f"
	for (j = 1; j <= i__1; ++j) {
#line 286 "slascl.f"
	    i__2 = *m;
#line 286 "slascl.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 287 "slascl.f"
		a[i__ + j * a_dim1] *= mul;
#line 288 "slascl.f"
/* L20: */
#line 288 "slascl.f"
	    }
#line 289 "slascl.f"
/* L30: */
#line 289 "slascl.f"
	}

#line 291 "slascl.f"
    } else if (itype == 1) {

/*        Lower triangular matrix */

#line 295 "slascl.f"
	i__1 = *n;
#line 295 "slascl.f"
	for (j = 1; j <= i__1; ++j) {
#line 296 "slascl.f"
	    i__2 = *m;
#line 296 "slascl.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
#line 297 "slascl.f"
		a[i__ + j * a_dim1] *= mul;
#line 298 "slascl.f"
/* L40: */
#line 298 "slascl.f"
	    }
#line 299 "slascl.f"
/* L50: */
#line 299 "slascl.f"
	}

#line 301 "slascl.f"
    } else if (itype == 2) {

/*        Upper triangular matrix */

#line 305 "slascl.f"
	i__1 = *n;
#line 305 "slascl.f"
	for (j = 1; j <= i__1; ++j) {
#line 306 "slascl.f"
	    i__2 = min(j,*m);
#line 306 "slascl.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 307 "slascl.f"
		a[i__ + j * a_dim1] *= mul;
#line 308 "slascl.f"
/* L60: */
#line 308 "slascl.f"
	    }
#line 309 "slascl.f"
/* L70: */
#line 309 "slascl.f"
	}

#line 311 "slascl.f"
    } else if (itype == 3) {

/*        Upper Hessenberg matrix */

#line 315 "slascl.f"
	i__1 = *n;
#line 315 "slascl.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 316 "slascl.f"
	    i__3 = j + 1;
#line 316 "slascl.f"
	    i__2 = min(i__3,*m);
#line 316 "slascl.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 317 "slascl.f"
		a[i__ + j * a_dim1] *= mul;
#line 318 "slascl.f"
/* L80: */
#line 318 "slascl.f"
	    }
#line 319 "slascl.f"
/* L90: */
#line 319 "slascl.f"
	}

#line 321 "slascl.f"
    } else if (itype == 4) {

/*        Lower half of a symmetric band matrix */

#line 325 "slascl.f"
	k3 = *kl + 1;
#line 326 "slascl.f"
	k4 = *n + 1;
#line 327 "slascl.f"
	i__1 = *n;
#line 327 "slascl.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 328 "slascl.f"
	    i__3 = k3, i__4 = k4 - j;
#line 328 "slascl.f"
	    i__2 = min(i__3,i__4);
#line 328 "slascl.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 329 "slascl.f"
		a[i__ + j * a_dim1] *= mul;
#line 330 "slascl.f"
/* L100: */
#line 330 "slascl.f"
	    }
#line 331 "slascl.f"
/* L110: */
#line 331 "slascl.f"
	}

#line 333 "slascl.f"
    } else if (itype == 5) {

/*        Upper half of a symmetric band matrix */

#line 337 "slascl.f"
	k1 = *ku + 2;
#line 338 "slascl.f"
	k3 = *ku + 1;
#line 339 "slascl.f"
	i__1 = *n;
#line 339 "slascl.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 340 "slascl.f"
	    i__2 = k1 - j;
#line 340 "slascl.f"
	    i__3 = k3;
#line 340 "slascl.f"
	    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 341 "slascl.f"
		a[i__ + j * a_dim1] *= mul;
#line 342 "slascl.f"
/* L120: */
#line 342 "slascl.f"
	    }
#line 343 "slascl.f"
/* L130: */
#line 343 "slascl.f"
	}

#line 345 "slascl.f"
    } else if (itype == 6) {

/*        Band matrix */

#line 349 "slascl.f"
	k1 = *kl + *ku + 2;
#line 350 "slascl.f"
	k2 = *kl + 1;
#line 351 "slascl.f"
	k3 = (*kl << 1) + *ku + 1;
#line 352 "slascl.f"
	k4 = *kl + *ku + 1 + *m;
#line 353 "slascl.f"
	i__1 = *n;
#line 353 "slascl.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 354 "slascl.f"
	    i__3 = k1 - j;
/* Computing MIN */
#line 354 "slascl.f"
	    i__4 = k3, i__5 = k4 - j;
#line 354 "slascl.f"
	    i__2 = min(i__4,i__5);
#line 354 "slascl.f"
	    for (i__ = max(i__3,k2); i__ <= i__2; ++i__) {
#line 355 "slascl.f"
		a[i__ + j * a_dim1] *= mul;
#line 356 "slascl.f"
/* L140: */
#line 356 "slascl.f"
	    }
#line 357 "slascl.f"
/* L150: */
#line 357 "slascl.f"
	}

#line 359 "slascl.f"
    }

#line 361 "slascl.f"
    if (! done) {
#line 361 "slascl.f"
	goto L10;
#line 361 "slascl.f"
    }

#line 364 "slascl.f"
    return 0;

/*     End of SLASCL */

} /* slascl_ */

