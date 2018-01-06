#line 1 "dspgst.f"
/* dspgst.f -- translated by f2c (version 20100827).
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

#line 1 "dspgst.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b9 = -1.;
static doublereal c_b11 = 1.;

/* > \brief \b DSPGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSPGST + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dspgst.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dspgst.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dspgst.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSPGST( ITYPE, UPLO, N, AP, BP, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, ITYPE, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AP( * ), BP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSPGST reduces a real symmetric-definite generalized eigenproblem */
/* > to standard form, using packed storage. */
/* > */
/* > If ITYPE = 1, the problem is A*x = lambda*B*x, */
/* > and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T) */
/* > */
/* > If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or */
/* > B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L. */
/* > */
/* > B must have been previously factorized as U**T*U or L*L**T by DPPTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ITYPE */
/* > \verbatim */
/* >          ITYPE is INTEGER */
/* >          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T); */
/* >          = 2 or 3: compute U*A*U**T or L**T*A*L. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored and B is factored as */
/* >                  U**T*U; */
/* >          = 'L':  Lower triangle of A is stored and B is factored as */
/* >                  L*L**T. */
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
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > */
/* >          On exit, if INFO = 0, the transformed matrix, stored in the */
/* >          same format as A. */
/* > \endverbatim */
/* > */
/* > \param[in] BP */
/* > \verbatim */
/* >          BP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* >          The triangular factor from the Cholesky factorization of B, */
/* >          stored in the same format as A, as returned by DPPTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
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
/* Subroutine */ int dspgst_(integer *itype, char *uplo, integer *n, 
	doublereal *ap, doublereal *bp, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer j, k, j1, k1, jj, kk;
    static doublereal ct, ajj;
    static integer j1j1;
    static doublereal akk;
    static integer k1k1;
    static doublereal bjj, bkk;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int dspr2_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    ftnlen), dscal_(integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dspmv_(char *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *, ftnlen);
    static logical upper;
    extern /* Subroutine */ int dtpmv_(char *, char *, char *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen), 
	    dtpsv_(char *, char *, char *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 153 "dspgst.f"
    /* Parameter adjustments */
#line 153 "dspgst.f"
    --bp;
#line 153 "dspgst.f"
    --ap;
#line 153 "dspgst.f"

#line 153 "dspgst.f"
    /* Function Body */
#line 153 "dspgst.f"
    *info = 0;
#line 154 "dspgst.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 155 "dspgst.f"
    if (*itype < 1 || *itype > 3) {
#line 156 "dspgst.f"
	*info = -1;
#line 157 "dspgst.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 158 "dspgst.f"
	*info = -2;
#line 159 "dspgst.f"
    } else if (*n < 0) {
#line 160 "dspgst.f"
	*info = -3;
#line 161 "dspgst.f"
    }
#line 162 "dspgst.f"
    if (*info != 0) {
#line 163 "dspgst.f"
	i__1 = -(*info);
#line 163 "dspgst.f"
	xerbla_("DSPGST", &i__1, (ftnlen)6);
#line 164 "dspgst.f"
	return 0;
#line 165 "dspgst.f"
    }

#line 167 "dspgst.f"
    if (*itype == 1) {
#line 168 "dspgst.f"
	if (upper) {

/*           Compute inv(U**T)*A*inv(U) */

/*           J1 and JJ are the indices of A(1,j) and A(j,j) */

#line 174 "dspgst.f"
	    jj = 0;
#line 175 "dspgst.f"
	    i__1 = *n;
#line 175 "dspgst.f"
	    for (j = 1; j <= i__1; ++j) {
#line 176 "dspgst.f"
		j1 = jj + 1;
#line 177 "dspgst.f"
		jj += j;

/*              Compute the j-th column of the upper triangle of A */

#line 181 "dspgst.f"
		bjj = bp[jj];
#line 182 "dspgst.f"
		dtpsv_(uplo, "Transpose", "Nonunit", &j, &bp[1], &ap[j1], &
			c__1, (ftnlen)1, (ftnlen)9, (ftnlen)7);
#line 184 "dspgst.f"
		i__2 = j - 1;
#line 184 "dspgst.f"
		dspmv_(uplo, &i__2, &c_b9, &ap[1], &bp[j1], &c__1, &c_b11, &
			ap[j1], &c__1, (ftnlen)1);
#line 186 "dspgst.f"
		i__2 = j - 1;
#line 186 "dspgst.f"
		d__1 = 1. / bjj;
#line 186 "dspgst.f"
		dscal_(&i__2, &d__1, &ap[j1], &c__1);
#line 187 "dspgst.f"
		i__2 = j - 1;
#line 187 "dspgst.f"
		ap[jj] = (ap[jj] - ddot_(&i__2, &ap[j1], &c__1, &bp[j1], &
			c__1)) / bjj;
#line 189 "dspgst.f"
/* L10: */
#line 189 "dspgst.f"
	    }
#line 190 "dspgst.f"
	} else {

/*           Compute inv(L)*A*inv(L**T) */

/*           KK and K1K1 are the indices of A(k,k) and A(k+1,k+1) */

#line 196 "dspgst.f"
	    kk = 1;
#line 197 "dspgst.f"
	    i__1 = *n;
#line 197 "dspgst.f"
	    for (k = 1; k <= i__1; ++k) {
#line 198 "dspgst.f"
		k1k1 = kk + *n - k + 1;

/*              Update the lower triangle of A(k:n,k:n) */

#line 202 "dspgst.f"
		akk = ap[kk];
#line 203 "dspgst.f"
		bkk = bp[kk];
/* Computing 2nd power */
#line 204 "dspgst.f"
		d__1 = bkk;
#line 204 "dspgst.f"
		akk /= d__1 * d__1;
#line 205 "dspgst.f"
		ap[kk] = akk;
#line 206 "dspgst.f"
		if (k < *n) {
#line 207 "dspgst.f"
		    i__2 = *n - k;
#line 207 "dspgst.f"
		    d__1 = 1. / bkk;
#line 207 "dspgst.f"
		    dscal_(&i__2, &d__1, &ap[kk + 1], &c__1);
#line 208 "dspgst.f"
		    ct = akk * -.5;
#line 209 "dspgst.f"
		    i__2 = *n - k;
#line 209 "dspgst.f"
		    daxpy_(&i__2, &ct, &bp[kk + 1], &c__1, &ap[kk + 1], &c__1)
			    ;
#line 210 "dspgst.f"
		    i__2 = *n - k;
#line 210 "dspgst.f"
		    dspr2_(uplo, &i__2, &c_b9, &ap[kk + 1], &c__1, &bp[kk + 1]
			    , &c__1, &ap[k1k1], (ftnlen)1);
#line 212 "dspgst.f"
		    i__2 = *n - k;
#line 212 "dspgst.f"
		    daxpy_(&i__2, &ct, &bp[kk + 1], &c__1, &ap[kk + 1], &c__1)
			    ;
#line 213 "dspgst.f"
		    i__2 = *n - k;
#line 213 "dspgst.f"
		    dtpsv_(uplo, "No transpose", "Non-unit", &i__2, &bp[k1k1],
			     &ap[kk + 1], &c__1, (ftnlen)1, (ftnlen)12, (
			    ftnlen)8);
#line 215 "dspgst.f"
		}
#line 216 "dspgst.f"
		kk = k1k1;
#line 217 "dspgst.f"
/* L20: */
#line 217 "dspgst.f"
	    }
#line 218 "dspgst.f"
	}
#line 219 "dspgst.f"
    } else {
#line 220 "dspgst.f"
	if (upper) {

/*           Compute U*A*U**T */

/*           K1 and KK are the indices of A(1,k) and A(k,k) */

#line 226 "dspgst.f"
	    kk = 0;
#line 227 "dspgst.f"
	    i__1 = *n;
#line 227 "dspgst.f"
	    for (k = 1; k <= i__1; ++k) {
#line 228 "dspgst.f"
		k1 = kk + 1;
#line 229 "dspgst.f"
		kk += k;

/*              Update the upper triangle of A(1:k,1:k) */

#line 233 "dspgst.f"
		akk = ap[kk];
#line 234 "dspgst.f"
		bkk = bp[kk];
#line 235 "dspgst.f"
		i__2 = k - 1;
#line 235 "dspgst.f"
		dtpmv_(uplo, "No transpose", "Non-unit", &i__2, &bp[1], &ap[
			k1], &c__1, (ftnlen)1, (ftnlen)12, (ftnlen)8);
#line 237 "dspgst.f"
		ct = akk * .5;
#line 238 "dspgst.f"
		i__2 = k - 1;
#line 238 "dspgst.f"
		daxpy_(&i__2, &ct, &bp[k1], &c__1, &ap[k1], &c__1);
#line 239 "dspgst.f"
		i__2 = k - 1;
#line 239 "dspgst.f"
		dspr2_(uplo, &i__2, &c_b11, &ap[k1], &c__1, &bp[k1], &c__1, &
			ap[1], (ftnlen)1);
#line 241 "dspgst.f"
		i__2 = k - 1;
#line 241 "dspgst.f"
		daxpy_(&i__2, &ct, &bp[k1], &c__1, &ap[k1], &c__1);
#line 242 "dspgst.f"
		i__2 = k - 1;
#line 242 "dspgst.f"
		dscal_(&i__2, &bkk, &ap[k1], &c__1);
/* Computing 2nd power */
#line 243 "dspgst.f"
		d__1 = bkk;
#line 243 "dspgst.f"
		ap[kk] = akk * (d__1 * d__1);
#line 244 "dspgst.f"
/* L30: */
#line 244 "dspgst.f"
	    }
#line 245 "dspgst.f"
	} else {

/*           Compute L**T *A*L */

/*           JJ and J1J1 are the indices of A(j,j) and A(j+1,j+1) */

#line 251 "dspgst.f"
	    jj = 1;
#line 252 "dspgst.f"
	    i__1 = *n;
#line 252 "dspgst.f"
	    for (j = 1; j <= i__1; ++j) {
#line 253 "dspgst.f"
		j1j1 = jj + *n - j + 1;

/*              Compute the j-th column of the lower triangle of A */

#line 257 "dspgst.f"
		ajj = ap[jj];
#line 258 "dspgst.f"
		bjj = bp[jj];
#line 259 "dspgst.f"
		i__2 = *n - j;
#line 259 "dspgst.f"
		ap[jj] = ajj * bjj + ddot_(&i__2, &ap[jj + 1], &c__1, &bp[jj 
			+ 1], &c__1);
#line 261 "dspgst.f"
		i__2 = *n - j;
#line 261 "dspgst.f"
		dscal_(&i__2, &bjj, &ap[jj + 1], &c__1);
#line 262 "dspgst.f"
		i__2 = *n - j;
#line 262 "dspgst.f"
		dspmv_(uplo, &i__2, &c_b11, &ap[j1j1], &bp[jj + 1], &c__1, &
			c_b11, &ap[jj + 1], &c__1, (ftnlen)1);
#line 264 "dspgst.f"
		i__2 = *n - j + 1;
#line 264 "dspgst.f"
		dtpmv_(uplo, "Transpose", "Non-unit", &i__2, &bp[jj], &ap[jj],
			 &c__1, (ftnlen)1, (ftnlen)9, (ftnlen)8);
#line 266 "dspgst.f"
		jj = j1j1;
#line 267 "dspgst.f"
/* L40: */
#line 267 "dspgst.f"
	    }
#line 268 "dspgst.f"
	}
#line 269 "dspgst.f"
    }
#line 270 "dspgst.f"
    return 0;

/*     End of DSPGST */

} /* dspgst_ */

