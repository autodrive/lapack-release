#line 1 "clascl.f"
/* clascl.f -- translated by f2c (version 20100827).
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

#line 1 "clascl.f"
/* > \brief \b CLASCL multiplies a general rectangular matrix by a real scalar defined as cto/cfrom. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLASCL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clascl.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clascl.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clascl.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TYPE */
/*       INTEGER            INFO, KL, KU, LDA, M, N */
/*       REAL               CFROM, CTO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLASCL multiplies the M by N complex matrix A by the real scalar */
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
/* >                  bandwidth KU. See CGBTRF for storage details. */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int clascl_(char *type__, integer *kl, integer *ku, 
	doublereal *cfrom, doublereal *cto, integer *m, integer *n, 
	doublecomplex *a, integer *lda, integer *info, ftnlen type_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1;

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

#line 186 "clascl.f"
    /* Parameter adjustments */
#line 186 "clascl.f"
    a_dim1 = *lda;
#line 186 "clascl.f"
    a_offset = 1 + a_dim1;
#line 186 "clascl.f"
    a -= a_offset;
#line 186 "clascl.f"

#line 186 "clascl.f"
    /* Function Body */
#line 186 "clascl.f"
    *info = 0;

#line 188 "clascl.f"
    if (lsame_(type__, "G", (ftnlen)1, (ftnlen)1)) {
#line 189 "clascl.f"
	itype = 0;
#line 190 "clascl.f"
    } else if (lsame_(type__, "L", (ftnlen)1, (ftnlen)1)) {
#line 191 "clascl.f"
	itype = 1;
#line 192 "clascl.f"
    } else if (lsame_(type__, "U", (ftnlen)1, (ftnlen)1)) {
#line 193 "clascl.f"
	itype = 2;
#line 194 "clascl.f"
    } else if (lsame_(type__, "H", (ftnlen)1, (ftnlen)1)) {
#line 195 "clascl.f"
	itype = 3;
#line 196 "clascl.f"
    } else if (lsame_(type__, "B", (ftnlen)1, (ftnlen)1)) {
#line 197 "clascl.f"
	itype = 4;
#line 198 "clascl.f"
    } else if (lsame_(type__, "Q", (ftnlen)1, (ftnlen)1)) {
#line 199 "clascl.f"
	itype = 5;
#line 200 "clascl.f"
    } else if (lsame_(type__, "Z", (ftnlen)1, (ftnlen)1)) {
#line 201 "clascl.f"
	itype = 6;
#line 202 "clascl.f"
    } else {
#line 203 "clascl.f"
	itype = -1;
#line 204 "clascl.f"
    }

#line 206 "clascl.f"
    if (itype == -1) {
#line 207 "clascl.f"
	*info = -1;
#line 208 "clascl.f"
    } else if (*cfrom == 0. || sisnan_(cfrom)) {
#line 209 "clascl.f"
	*info = -4;
#line 210 "clascl.f"
    } else if (sisnan_(cto)) {
#line 211 "clascl.f"
	*info = -5;
#line 212 "clascl.f"
    } else if (*m < 0) {
#line 213 "clascl.f"
	*info = -6;
#line 214 "clascl.f"
    } else if (*n < 0 || itype == 4 && *n != *m || itype == 5 && *n != *m) {
#line 216 "clascl.f"
	*info = -7;
#line 217 "clascl.f"
    } else if (itype <= 3 && *lda < max(1,*m)) {
#line 218 "clascl.f"
	*info = -9;
#line 219 "clascl.f"
    } else if (itype >= 4) {
/* Computing MAX */
#line 220 "clascl.f"
	i__1 = *m - 1;
#line 220 "clascl.f"
	if (*kl < 0 || *kl > max(i__1,0)) {
#line 221 "clascl.f"
	    *info = -2;
#line 222 "clascl.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 222 "clascl.f"
	    i__1 = *n - 1;
#line 222 "clascl.f"
	    if (*ku < 0 || *ku > max(i__1,0) || (itype == 4 || itype == 5) && 
		    *kl != *ku) {
#line 225 "clascl.f"
		*info = -3;
#line 226 "clascl.f"
	    } else if (itype == 4 && *lda < *kl + 1 || itype == 5 && *lda < *
		    ku + 1 || itype == 6 && *lda < (*kl << 1) + *ku + 1) {
#line 229 "clascl.f"
		*info = -9;
#line 230 "clascl.f"
	    }
#line 230 "clascl.f"
	}
#line 231 "clascl.f"
    }

#line 233 "clascl.f"
    if (*info != 0) {
#line 234 "clascl.f"
	i__1 = -(*info);
#line 234 "clascl.f"
	xerbla_("CLASCL", &i__1, (ftnlen)6);
#line 235 "clascl.f"
	return 0;
#line 236 "clascl.f"
    }

/*     Quick return if possible */

#line 240 "clascl.f"
    if (*n == 0 || *m == 0) {
#line 240 "clascl.f"
	return 0;
#line 240 "clascl.f"
    }

/*     Get machine parameters */

#line 245 "clascl.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 246 "clascl.f"
    bignum = 1. / smlnum;

#line 248 "clascl.f"
    cfromc = *cfrom;
#line 249 "clascl.f"
    ctoc = *cto;

#line 251 "clascl.f"
L10:
#line 252 "clascl.f"
    cfrom1 = cfromc * smlnum;
#line 253 "clascl.f"
    if (cfrom1 == cfromc) {
/*        CFROMC is an inf.  Multiply by a correctly signed zero for */
/*        finite CTOC, or a NaN if CTOC is infinite. */
#line 256 "clascl.f"
	mul = ctoc / cfromc;
#line 257 "clascl.f"
	done = TRUE_;
#line 258 "clascl.f"
	cto1 = ctoc;
#line 259 "clascl.f"
    } else {
#line 260 "clascl.f"
	cto1 = ctoc / bignum;
#line 261 "clascl.f"
	if (cto1 == ctoc) {
/*           CTOC is either 0 or an inf.  In both cases, CTOC itself */
/*           serves as the correct multiplication factor. */
#line 264 "clascl.f"
	    mul = ctoc;
#line 265 "clascl.f"
	    done = TRUE_;
#line 266 "clascl.f"
	    cfromc = 1.;
#line 267 "clascl.f"
	} else if (abs(cfrom1) > abs(ctoc) && ctoc != 0.) {
#line 268 "clascl.f"
	    mul = smlnum;
#line 269 "clascl.f"
	    done = FALSE_;
#line 270 "clascl.f"
	    cfromc = cfrom1;
#line 271 "clascl.f"
	} else if (abs(cto1) > abs(cfromc)) {
#line 272 "clascl.f"
	    mul = bignum;
#line 273 "clascl.f"
	    done = FALSE_;
#line 274 "clascl.f"
	    ctoc = cto1;
#line 275 "clascl.f"
	} else {
#line 276 "clascl.f"
	    mul = ctoc / cfromc;
#line 277 "clascl.f"
	    done = TRUE_;
#line 278 "clascl.f"
	}
#line 279 "clascl.f"
    }

#line 281 "clascl.f"
    if (itype == 0) {

/*        Full matrix */

#line 285 "clascl.f"
	i__1 = *n;
#line 285 "clascl.f"
	for (j = 1; j <= i__1; ++j) {
#line 286 "clascl.f"
	    i__2 = *m;
#line 286 "clascl.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 287 "clascl.f"
		i__3 = i__ + j * a_dim1;
#line 287 "clascl.f"
		i__4 = i__ + j * a_dim1;
#line 287 "clascl.f"
		z__1.r = mul * a[i__4].r, z__1.i = mul * a[i__4].i;
#line 287 "clascl.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 288 "clascl.f"
/* L20: */
#line 288 "clascl.f"
	    }
#line 289 "clascl.f"
/* L30: */
#line 289 "clascl.f"
	}

#line 291 "clascl.f"
    } else if (itype == 1) {

/*        Lower triangular matrix */

#line 295 "clascl.f"
	i__1 = *n;
#line 295 "clascl.f"
	for (j = 1; j <= i__1; ++j) {
#line 296 "clascl.f"
	    i__2 = *m;
#line 296 "clascl.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
#line 297 "clascl.f"
		i__3 = i__ + j * a_dim1;
#line 297 "clascl.f"
		i__4 = i__ + j * a_dim1;
#line 297 "clascl.f"
		z__1.r = mul * a[i__4].r, z__1.i = mul * a[i__4].i;
#line 297 "clascl.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 298 "clascl.f"
/* L40: */
#line 298 "clascl.f"
	    }
#line 299 "clascl.f"
/* L50: */
#line 299 "clascl.f"
	}

#line 301 "clascl.f"
    } else if (itype == 2) {

/*        Upper triangular matrix */

#line 305 "clascl.f"
	i__1 = *n;
#line 305 "clascl.f"
	for (j = 1; j <= i__1; ++j) {
#line 306 "clascl.f"
	    i__2 = min(j,*m);
#line 306 "clascl.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 307 "clascl.f"
		i__3 = i__ + j * a_dim1;
#line 307 "clascl.f"
		i__4 = i__ + j * a_dim1;
#line 307 "clascl.f"
		z__1.r = mul * a[i__4].r, z__1.i = mul * a[i__4].i;
#line 307 "clascl.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 308 "clascl.f"
/* L60: */
#line 308 "clascl.f"
	    }
#line 309 "clascl.f"
/* L70: */
#line 309 "clascl.f"
	}

#line 311 "clascl.f"
    } else if (itype == 3) {

/*        Upper Hessenberg matrix */

#line 315 "clascl.f"
	i__1 = *n;
#line 315 "clascl.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 316 "clascl.f"
	    i__3 = j + 1;
#line 316 "clascl.f"
	    i__2 = min(i__3,*m);
#line 316 "clascl.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 317 "clascl.f"
		i__3 = i__ + j * a_dim1;
#line 317 "clascl.f"
		i__4 = i__ + j * a_dim1;
#line 317 "clascl.f"
		z__1.r = mul * a[i__4].r, z__1.i = mul * a[i__4].i;
#line 317 "clascl.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 318 "clascl.f"
/* L80: */
#line 318 "clascl.f"
	    }
#line 319 "clascl.f"
/* L90: */
#line 319 "clascl.f"
	}

#line 321 "clascl.f"
    } else if (itype == 4) {

/*        Lower half of a symmetric band matrix */

#line 325 "clascl.f"
	k3 = *kl + 1;
#line 326 "clascl.f"
	k4 = *n + 1;
#line 327 "clascl.f"
	i__1 = *n;
#line 327 "clascl.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 328 "clascl.f"
	    i__3 = k3, i__4 = k4 - j;
#line 328 "clascl.f"
	    i__2 = min(i__3,i__4);
#line 328 "clascl.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 329 "clascl.f"
		i__3 = i__ + j * a_dim1;
#line 329 "clascl.f"
		i__4 = i__ + j * a_dim1;
#line 329 "clascl.f"
		z__1.r = mul * a[i__4].r, z__1.i = mul * a[i__4].i;
#line 329 "clascl.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 330 "clascl.f"
/* L100: */
#line 330 "clascl.f"
	    }
#line 331 "clascl.f"
/* L110: */
#line 331 "clascl.f"
	}

#line 333 "clascl.f"
    } else if (itype == 5) {

/*        Upper half of a symmetric band matrix */

#line 337 "clascl.f"
	k1 = *ku + 2;
#line 338 "clascl.f"
	k3 = *ku + 1;
#line 339 "clascl.f"
	i__1 = *n;
#line 339 "clascl.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 340 "clascl.f"
	    i__2 = k1 - j;
#line 340 "clascl.f"
	    i__3 = k3;
#line 340 "clascl.f"
	    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 341 "clascl.f"
		i__2 = i__ + j * a_dim1;
#line 341 "clascl.f"
		i__4 = i__ + j * a_dim1;
#line 341 "clascl.f"
		z__1.r = mul * a[i__4].r, z__1.i = mul * a[i__4].i;
#line 341 "clascl.f"
		a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 342 "clascl.f"
/* L120: */
#line 342 "clascl.f"
	    }
#line 343 "clascl.f"
/* L130: */
#line 343 "clascl.f"
	}

#line 345 "clascl.f"
    } else if (itype == 6) {

/*        Band matrix */

#line 349 "clascl.f"
	k1 = *kl + *ku + 2;
#line 350 "clascl.f"
	k2 = *kl + 1;
#line 351 "clascl.f"
	k3 = (*kl << 1) + *ku + 1;
#line 352 "clascl.f"
	k4 = *kl + *ku + 1 + *m;
#line 353 "clascl.f"
	i__1 = *n;
#line 353 "clascl.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 354 "clascl.f"
	    i__3 = k1 - j;
/* Computing MIN */
#line 354 "clascl.f"
	    i__4 = k3, i__5 = k4 - j;
#line 354 "clascl.f"
	    i__2 = min(i__4,i__5);
#line 354 "clascl.f"
	    for (i__ = max(i__3,k2); i__ <= i__2; ++i__) {
#line 355 "clascl.f"
		i__3 = i__ + j * a_dim1;
#line 355 "clascl.f"
		i__4 = i__ + j * a_dim1;
#line 355 "clascl.f"
		z__1.r = mul * a[i__4].r, z__1.i = mul * a[i__4].i;
#line 355 "clascl.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 356 "clascl.f"
/* L140: */
#line 356 "clascl.f"
	    }
#line 357 "clascl.f"
/* L150: */
#line 357 "clascl.f"
	}

#line 359 "clascl.f"
    }

#line 361 "clascl.f"
    if (! done) {
#line 361 "clascl.f"
	goto L10;
#line 361 "clascl.f"
    }

#line 364 "clascl.f"
    return 0;

/*     End of CLASCL */

} /* clascl_ */

