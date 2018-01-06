#line 1 "zlascl.f"
/* zlascl.f -- translated by f2c (version 20100827).
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

#line 1 "zlascl.f"
/* > \brief \b ZLASCL multiplies a general rectangular matrix by a real scalar defined as cto/cfrom. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLASCL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlascl.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlascl.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlascl.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TYPE */
/*       INTEGER            INFO, KL, KU, LDA, M, N */
/*       DOUBLE PRECISION   CFROM, CTO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLASCL multiplies the M by N complex matrix A by the real scalar */
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
/* >                  bandwidth KU. See ZGBTRF for storage details. */
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
/* >          CFROM is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] CTO */
/* > \verbatim */
/* >          CTO is DOUBLE PRECISION */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          The matrix to be multiplied by CTO/CFROM.  See TYPE for the */
/* >          storage type. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
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

/* > \date September 2012 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zlascl_(char *type__, integer *kl, integer *ku, 
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
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal cfromc;
    extern logical disnan_(doublereal *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum, smlnum;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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

#line 182 "zlascl.f"
    /* Parameter adjustments */
#line 182 "zlascl.f"
    a_dim1 = *lda;
#line 182 "zlascl.f"
    a_offset = 1 + a_dim1;
#line 182 "zlascl.f"
    a -= a_offset;
#line 182 "zlascl.f"

#line 182 "zlascl.f"
    /* Function Body */
#line 182 "zlascl.f"
    *info = 0;

#line 184 "zlascl.f"
    if (lsame_(type__, "G", (ftnlen)1, (ftnlen)1)) {
#line 185 "zlascl.f"
	itype = 0;
#line 186 "zlascl.f"
    } else if (lsame_(type__, "L", (ftnlen)1, (ftnlen)1)) {
#line 187 "zlascl.f"
	itype = 1;
#line 188 "zlascl.f"
    } else if (lsame_(type__, "U", (ftnlen)1, (ftnlen)1)) {
#line 189 "zlascl.f"
	itype = 2;
#line 190 "zlascl.f"
    } else if (lsame_(type__, "H", (ftnlen)1, (ftnlen)1)) {
#line 191 "zlascl.f"
	itype = 3;
#line 192 "zlascl.f"
    } else if (lsame_(type__, "B", (ftnlen)1, (ftnlen)1)) {
#line 193 "zlascl.f"
	itype = 4;
#line 194 "zlascl.f"
    } else if (lsame_(type__, "Q", (ftnlen)1, (ftnlen)1)) {
#line 195 "zlascl.f"
	itype = 5;
#line 196 "zlascl.f"
    } else if (lsame_(type__, "Z", (ftnlen)1, (ftnlen)1)) {
#line 197 "zlascl.f"
	itype = 6;
#line 198 "zlascl.f"
    } else {
#line 199 "zlascl.f"
	itype = -1;
#line 200 "zlascl.f"
    }

#line 202 "zlascl.f"
    if (itype == -1) {
#line 203 "zlascl.f"
	*info = -1;
#line 204 "zlascl.f"
    } else if (*cfrom == 0. || disnan_(cfrom)) {
#line 205 "zlascl.f"
	*info = -4;
#line 206 "zlascl.f"
    } else if (disnan_(cto)) {
#line 207 "zlascl.f"
	*info = -5;
#line 208 "zlascl.f"
    } else if (*m < 0) {
#line 209 "zlascl.f"
	*info = -6;
#line 210 "zlascl.f"
    } else if (*n < 0 || itype == 4 && *n != *m || itype == 5 && *n != *m) {
#line 212 "zlascl.f"
	*info = -7;
#line 213 "zlascl.f"
    } else if (itype <= 3 && *lda < max(1,*m)) {
#line 214 "zlascl.f"
	*info = -9;
#line 215 "zlascl.f"
    } else if (itype >= 4) {
/* Computing MAX */
#line 216 "zlascl.f"
	i__1 = *m - 1;
#line 216 "zlascl.f"
	if (*kl < 0 || *kl > max(i__1,0)) {
#line 217 "zlascl.f"
	    *info = -2;
#line 218 "zlascl.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 218 "zlascl.f"
	    i__1 = *n - 1;
#line 218 "zlascl.f"
	    if (*ku < 0 || *ku > max(i__1,0) || (itype == 4 || itype == 5) && 
		    *kl != *ku) {
#line 221 "zlascl.f"
		*info = -3;
#line 222 "zlascl.f"
	    } else if (itype == 4 && *lda < *kl + 1 || itype == 5 && *lda < *
		    ku + 1 || itype == 6 && *lda < (*kl << 1) + *ku + 1) {
#line 225 "zlascl.f"
		*info = -9;
#line 226 "zlascl.f"
	    }
#line 226 "zlascl.f"
	}
#line 227 "zlascl.f"
    }

#line 229 "zlascl.f"
    if (*info != 0) {
#line 230 "zlascl.f"
	i__1 = -(*info);
#line 230 "zlascl.f"
	xerbla_("ZLASCL", &i__1, (ftnlen)6);
#line 231 "zlascl.f"
	return 0;
#line 232 "zlascl.f"
    }

/*     Quick return if possible */

#line 236 "zlascl.f"
    if (*n == 0 || *m == 0) {
#line 236 "zlascl.f"
	return 0;
#line 236 "zlascl.f"
    }

/*     Get machine parameters */

#line 241 "zlascl.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 242 "zlascl.f"
    bignum = 1. / smlnum;

#line 244 "zlascl.f"
    cfromc = *cfrom;
#line 245 "zlascl.f"
    ctoc = *cto;

#line 247 "zlascl.f"
L10:
#line 248 "zlascl.f"
    cfrom1 = cfromc * smlnum;
#line 249 "zlascl.f"
    if (cfrom1 == cfromc) {
/*        CFROMC is an inf.  Multiply by a correctly signed zero for */
/*        finite CTOC, or a NaN if CTOC is infinite. */
#line 252 "zlascl.f"
	mul = ctoc / cfromc;
#line 253 "zlascl.f"
	done = TRUE_;
#line 254 "zlascl.f"
	cto1 = ctoc;
#line 255 "zlascl.f"
    } else {
#line 256 "zlascl.f"
	cto1 = ctoc / bignum;
#line 257 "zlascl.f"
	if (cto1 == ctoc) {
/*           CTOC is either 0 or an inf.  In both cases, CTOC itself */
/*           serves as the correct multiplication factor. */
#line 260 "zlascl.f"
	    mul = ctoc;
#line 261 "zlascl.f"
	    done = TRUE_;
#line 262 "zlascl.f"
	    cfromc = 1.;
#line 263 "zlascl.f"
	} else if (abs(cfrom1) > abs(ctoc) && ctoc != 0.) {
#line 264 "zlascl.f"
	    mul = smlnum;
#line 265 "zlascl.f"
	    done = FALSE_;
#line 266 "zlascl.f"
	    cfromc = cfrom1;
#line 267 "zlascl.f"
	} else if (abs(cto1) > abs(cfromc)) {
#line 268 "zlascl.f"
	    mul = bignum;
#line 269 "zlascl.f"
	    done = FALSE_;
#line 270 "zlascl.f"
	    ctoc = cto1;
#line 271 "zlascl.f"
	} else {
#line 272 "zlascl.f"
	    mul = ctoc / cfromc;
#line 273 "zlascl.f"
	    done = TRUE_;
#line 274 "zlascl.f"
	}
#line 275 "zlascl.f"
    }

#line 277 "zlascl.f"
    if (itype == 0) {

/*        Full matrix */

#line 281 "zlascl.f"
	i__1 = *n;
#line 281 "zlascl.f"
	for (j = 1; j <= i__1; ++j) {
#line 282 "zlascl.f"
	    i__2 = *m;
#line 282 "zlascl.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 283 "zlascl.f"
		i__3 = i__ + j * a_dim1;
#line 283 "zlascl.f"
		i__4 = i__ + j * a_dim1;
#line 283 "zlascl.f"
		z__1.r = mul * a[i__4].r, z__1.i = mul * a[i__4].i;
#line 283 "zlascl.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 284 "zlascl.f"
/* L20: */
#line 284 "zlascl.f"
	    }
#line 285 "zlascl.f"
/* L30: */
#line 285 "zlascl.f"
	}

#line 287 "zlascl.f"
    } else if (itype == 1) {

/*        Lower triangular matrix */

#line 291 "zlascl.f"
	i__1 = *n;
#line 291 "zlascl.f"
	for (j = 1; j <= i__1; ++j) {
#line 292 "zlascl.f"
	    i__2 = *m;
#line 292 "zlascl.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
#line 293 "zlascl.f"
		i__3 = i__ + j * a_dim1;
#line 293 "zlascl.f"
		i__4 = i__ + j * a_dim1;
#line 293 "zlascl.f"
		z__1.r = mul * a[i__4].r, z__1.i = mul * a[i__4].i;
#line 293 "zlascl.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 294 "zlascl.f"
/* L40: */
#line 294 "zlascl.f"
	    }
#line 295 "zlascl.f"
/* L50: */
#line 295 "zlascl.f"
	}

#line 297 "zlascl.f"
    } else if (itype == 2) {

/*        Upper triangular matrix */

#line 301 "zlascl.f"
	i__1 = *n;
#line 301 "zlascl.f"
	for (j = 1; j <= i__1; ++j) {
#line 302 "zlascl.f"
	    i__2 = min(j,*m);
#line 302 "zlascl.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 303 "zlascl.f"
		i__3 = i__ + j * a_dim1;
#line 303 "zlascl.f"
		i__4 = i__ + j * a_dim1;
#line 303 "zlascl.f"
		z__1.r = mul * a[i__4].r, z__1.i = mul * a[i__4].i;
#line 303 "zlascl.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 304 "zlascl.f"
/* L60: */
#line 304 "zlascl.f"
	    }
#line 305 "zlascl.f"
/* L70: */
#line 305 "zlascl.f"
	}

#line 307 "zlascl.f"
    } else if (itype == 3) {

/*        Upper Hessenberg matrix */

#line 311 "zlascl.f"
	i__1 = *n;
#line 311 "zlascl.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 312 "zlascl.f"
	    i__3 = j + 1;
#line 312 "zlascl.f"
	    i__2 = min(i__3,*m);
#line 312 "zlascl.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 313 "zlascl.f"
		i__3 = i__ + j * a_dim1;
#line 313 "zlascl.f"
		i__4 = i__ + j * a_dim1;
#line 313 "zlascl.f"
		z__1.r = mul * a[i__4].r, z__1.i = mul * a[i__4].i;
#line 313 "zlascl.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 314 "zlascl.f"
/* L80: */
#line 314 "zlascl.f"
	    }
#line 315 "zlascl.f"
/* L90: */
#line 315 "zlascl.f"
	}

#line 317 "zlascl.f"
    } else if (itype == 4) {

/*        Lower half of a symmetric band matrix */

#line 321 "zlascl.f"
	k3 = *kl + 1;
#line 322 "zlascl.f"
	k4 = *n + 1;
#line 323 "zlascl.f"
	i__1 = *n;
#line 323 "zlascl.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 324 "zlascl.f"
	    i__3 = k3, i__4 = k4 - j;
#line 324 "zlascl.f"
	    i__2 = min(i__3,i__4);
#line 324 "zlascl.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 325 "zlascl.f"
		i__3 = i__ + j * a_dim1;
#line 325 "zlascl.f"
		i__4 = i__ + j * a_dim1;
#line 325 "zlascl.f"
		z__1.r = mul * a[i__4].r, z__1.i = mul * a[i__4].i;
#line 325 "zlascl.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 326 "zlascl.f"
/* L100: */
#line 326 "zlascl.f"
	    }
#line 327 "zlascl.f"
/* L110: */
#line 327 "zlascl.f"
	}

#line 329 "zlascl.f"
    } else if (itype == 5) {

/*        Upper half of a symmetric band matrix */

#line 333 "zlascl.f"
	k1 = *ku + 2;
#line 334 "zlascl.f"
	k3 = *ku + 1;
#line 335 "zlascl.f"
	i__1 = *n;
#line 335 "zlascl.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 336 "zlascl.f"
	    i__2 = k1 - j;
#line 336 "zlascl.f"
	    i__3 = k3;
#line 336 "zlascl.f"
	    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 337 "zlascl.f"
		i__2 = i__ + j * a_dim1;
#line 337 "zlascl.f"
		i__4 = i__ + j * a_dim1;
#line 337 "zlascl.f"
		z__1.r = mul * a[i__4].r, z__1.i = mul * a[i__4].i;
#line 337 "zlascl.f"
		a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 338 "zlascl.f"
/* L120: */
#line 338 "zlascl.f"
	    }
#line 339 "zlascl.f"
/* L130: */
#line 339 "zlascl.f"
	}

#line 341 "zlascl.f"
    } else if (itype == 6) {

/*        Band matrix */

#line 345 "zlascl.f"
	k1 = *kl + *ku + 2;
#line 346 "zlascl.f"
	k2 = *kl + 1;
#line 347 "zlascl.f"
	k3 = (*kl << 1) + *ku + 1;
#line 348 "zlascl.f"
	k4 = *kl + *ku + 1 + *m;
#line 349 "zlascl.f"
	i__1 = *n;
#line 349 "zlascl.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 350 "zlascl.f"
	    i__3 = k1 - j;
/* Computing MIN */
#line 350 "zlascl.f"
	    i__4 = k3, i__5 = k4 - j;
#line 350 "zlascl.f"
	    i__2 = min(i__4,i__5);
#line 350 "zlascl.f"
	    for (i__ = max(i__3,k2); i__ <= i__2; ++i__) {
#line 351 "zlascl.f"
		i__3 = i__ + j * a_dim1;
#line 351 "zlascl.f"
		i__4 = i__ + j * a_dim1;
#line 351 "zlascl.f"
		z__1.r = mul * a[i__4].r, z__1.i = mul * a[i__4].i;
#line 351 "zlascl.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 352 "zlascl.f"
/* L140: */
#line 352 "zlascl.f"
	    }
#line 353 "zlascl.f"
/* L150: */
#line 353 "zlascl.f"
	}

#line 355 "zlascl.f"
    }

#line 357 "zlascl.f"
    if (! done) {
#line 357 "zlascl.f"
	goto L10;
#line 357 "zlascl.f"
    }

#line 360 "zlascl.f"
    return 0;

/*     End of ZLASCL */

} /* zlascl_ */

