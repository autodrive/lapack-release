#line 1 "chpgst.f"
/* chpgst.f -- translated by f2c (version 20100827).
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

#line 1 "chpgst.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CHPGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHPGST + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chpgst.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chpgst.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chpgst.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHPGST( ITYPE, UPLO, N, AP, BP, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, ITYPE, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            AP( * ), BP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHPGST reduces a complex Hermitian-definite generalized */
/* > eigenproblem to standard form, using packed storage. */
/* > */
/* > If ITYPE = 1, the problem is A*x = lambda*B*x, */
/* > and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H) */
/* > */
/* > If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or */
/* > B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L. */
/* > */
/* > B must have been previously factorized as U**H*U or L*L**H by CPPTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ITYPE */
/* > \verbatim */
/* >          ITYPE is INTEGER */
/* >          = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H); */
/* >          = 2 or 3: compute U*A*U**H or L**H*A*L. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored and B is factored as */
/* >                  U**H*U; */
/* >          = 'L':  Lower triangle of A is stored and B is factored as */
/* >                  L*L**H. */
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
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the Hermitian matrix */
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
/* >          BP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          The triangular factor from the Cholesky factorization of B, */
/* >          stored in the same format as A, as returned by CPPTRF. */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int chpgst_(integer *itype, char *uplo, integer *n, 
	doublecomplex *ap, doublecomplex *bp, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;

    /* Local variables */
    static integer j, k, j1, k1, jj, kk;
    static doublecomplex ct;
    static doublereal ajj;
    static integer j1j1;
    static doublereal akk;
    static integer k1k1;
    static doublereal bjj, bkk;
    extern /* Subroutine */ int chpr2_(char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, ftnlen);
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int chpmv_(char *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, ftnlen), caxpy_(integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *), ctpmv_(char *, char *, char *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static logical upper;
    extern /* Subroutine */ int ctpsv_(char *, char *, char *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen, 
	    ftnlen), csscal_(integer *, doublereal *, doublecomplex *, 
	    integer *), xerbla_(char *, integer *, ftnlen);


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
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 159 "chpgst.f"
    /* Parameter adjustments */
#line 159 "chpgst.f"
    --bp;
#line 159 "chpgst.f"
    --ap;
#line 159 "chpgst.f"

#line 159 "chpgst.f"
    /* Function Body */
#line 159 "chpgst.f"
    *info = 0;
#line 160 "chpgst.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 161 "chpgst.f"
    if (*itype < 1 || *itype > 3) {
#line 162 "chpgst.f"
	*info = -1;
#line 163 "chpgst.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 164 "chpgst.f"
	*info = -2;
#line 165 "chpgst.f"
    } else if (*n < 0) {
#line 166 "chpgst.f"
	*info = -3;
#line 167 "chpgst.f"
    }
#line 168 "chpgst.f"
    if (*info != 0) {
#line 169 "chpgst.f"
	i__1 = -(*info);
#line 169 "chpgst.f"
	xerbla_("CHPGST", &i__1, (ftnlen)6);
#line 170 "chpgst.f"
	return 0;
#line 171 "chpgst.f"
    }

#line 173 "chpgst.f"
    if (*itype == 1) {
#line 174 "chpgst.f"
	if (upper) {

/*           Compute inv(U**H)*A*inv(U) */

/*           J1 and JJ are the indices of A(1,j) and A(j,j) */

#line 180 "chpgst.f"
	    jj = 0;
#line 181 "chpgst.f"
	    i__1 = *n;
#line 181 "chpgst.f"
	    for (j = 1; j <= i__1; ++j) {
#line 182 "chpgst.f"
		j1 = jj + 1;
#line 183 "chpgst.f"
		jj += j;

/*              Compute the j-th column of the upper triangle of A */

#line 187 "chpgst.f"
		i__2 = jj;
#line 187 "chpgst.f"
		i__3 = jj;
#line 187 "chpgst.f"
		d__1 = ap[i__3].r;
#line 187 "chpgst.f"
		ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 188 "chpgst.f"
		i__2 = jj;
#line 188 "chpgst.f"
		bjj = bp[i__2].r;
#line 189 "chpgst.f"
		ctpsv_(uplo, "Conjugate transpose", "Non-unit", &j, &bp[1], &
			ap[j1], &c__1, (ftnlen)1, (ftnlen)19, (ftnlen)8);
#line 191 "chpgst.f"
		i__2 = j - 1;
#line 191 "chpgst.f"
		z__1.r = -1., z__1.i = -0.;
#line 191 "chpgst.f"
		chpmv_(uplo, &i__2, &z__1, &ap[1], &bp[j1], &c__1, &c_b1, &ap[
			j1], &c__1, (ftnlen)1);
#line 193 "chpgst.f"
		i__2 = j - 1;
#line 193 "chpgst.f"
		d__1 = 1. / bjj;
#line 193 "chpgst.f"
		csscal_(&i__2, &d__1, &ap[j1], &c__1);
#line 194 "chpgst.f"
		i__2 = jj;
#line 194 "chpgst.f"
		i__3 = jj;
#line 194 "chpgst.f"
		i__4 = j - 1;
#line 194 "chpgst.f"
		cdotc_(&z__3, &i__4, &ap[j1], &c__1, &bp[j1], &c__1);
#line 194 "chpgst.f"
		z__2.r = ap[i__3].r - z__3.r, z__2.i = ap[i__3].i - z__3.i;
#line 194 "chpgst.f"
		z__1.r = z__2.r / bjj, z__1.i = z__2.i / bjj;
#line 194 "chpgst.f"
		ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
#line 196 "chpgst.f"
/* L10: */
#line 196 "chpgst.f"
	    }
#line 197 "chpgst.f"
	} else {

/*           Compute inv(L)*A*inv(L**H) */

/*           KK and K1K1 are the indices of A(k,k) and A(k+1,k+1) */

#line 203 "chpgst.f"
	    kk = 1;
#line 204 "chpgst.f"
	    i__1 = *n;
#line 204 "chpgst.f"
	    for (k = 1; k <= i__1; ++k) {
#line 205 "chpgst.f"
		k1k1 = kk + *n - k + 1;

/*              Update the lower triangle of A(k:n,k:n) */

#line 209 "chpgst.f"
		i__2 = kk;
#line 209 "chpgst.f"
		akk = ap[i__2].r;
#line 210 "chpgst.f"
		i__2 = kk;
#line 210 "chpgst.f"
		bkk = bp[i__2].r;
/* Computing 2nd power */
#line 211 "chpgst.f"
		d__1 = bkk;
#line 211 "chpgst.f"
		akk /= d__1 * d__1;
#line 212 "chpgst.f"
		i__2 = kk;
#line 212 "chpgst.f"
		ap[i__2].r = akk, ap[i__2].i = 0.;
#line 213 "chpgst.f"
		if (k < *n) {
#line 214 "chpgst.f"
		    i__2 = *n - k;
#line 214 "chpgst.f"
		    d__1 = 1. / bkk;
#line 214 "chpgst.f"
		    csscal_(&i__2, &d__1, &ap[kk + 1], &c__1);
#line 215 "chpgst.f"
		    d__1 = akk * -.5;
#line 215 "chpgst.f"
		    ct.r = d__1, ct.i = 0.;
#line 216 "chpgst.f"
		    i__2 = *n - k;
#line 216 "chpgst.f"
		    caxpy_(&i__2, &ct, &bp[kk + 1], &c__1, &ap[kk + 1], &c__1)
			    ;
#line 217 "chpgst.f"
		    i__2 = *n - k;
#line 217 "chpgst.f"
		    z__1.r = -1., z__1.i = -0.;
#line 217 "chpgst.f"
		    chpr2_(uplo, &i__2, &z__1, &ap[kk + 1], &c__1, &bp[kk + 1]
			    , &c__1, &ap[k1k1], (ftnlen)1);
#line 219 "chpgst.f"
		    i__2 = *n - k;
#line 219 "chpgst.f"
		    caxpy_(&i__2, &ct, &bp[kk + 1], &c__1, &ap[kk + 1], &c__1)
			    ;
#line 220 "chpgst.f"
		    i__2 = *n - k;
#line 220 "chpgst.f"
		    ctpsv_(uplo, "No transpose", "Non-unit", &i__2, &bp[k1k1],
			     &ap[kk + 1], &c__1, (ftnlen)1, (ftnlen)12, (
			    ftnlen)8);
#line 222 "chpgst.f"
		}
#line 223 "chpgst.f"
		kk = k1k1;
#line 224 "chpgst.f"
/* L20: */
#line 224 "chpgst.f"
	    }
#line 225 "chpgst.f"
	}
#line 226 "chpgst.f"
    } else {
#line 227 "chpgst.f"
	if (upper) {

/*           Compute U*A*U**H */

/*           K1 and KK are the indices of A(1,k) and A(k,k) */

#line 233 "chpgst.f"
	    kk = 0;
#line 234 "chpgst.f"
	    i__1 = *n;
#line 234 "chpgst.f"
	    for (k = 1; k <= i__1; ++k) {
#line 235 "chpgst.f"
		k1 = kk + 1;
#line 236 "chpgst.f"
		kk += k;

/*              Update the upper triangle of A(1:k,1:k) */

#line 240 "chpgst.f"
		i__2 = kk;
#line 240 "chpgst.f"
		akk = ap[i__2].r;
#line 241 "chpgst.f"
		i__2 = kk;
#line 241 "chpgst.f"
		bkk = bp[i__2].r;
#line 242 "chpgst.f"
		i__2 = k - 1;
#line 242 "chpgst.f"
		ctpmv_(uplo, "No transpose", "Non-unit", &i__2, &bp[1], &ap[
			k1], &c__1, (ftnlen)1, (ftnlen)12, (ftnlen)8);
#line 244 "chpgst.f"
		d__1 = akk * .5;
#line 244 "chpgst.f"
		ct.r = d__1, ct.i = 0.;
#line 245 "chpgst.f"
		i__2 = k - 1;
#line 245 "chpgst.f"
		caxpy_(&i__2, &ct, &bp[k1], &c__1, &ap[k1], &c__1);
#line 246 "chpgst.f"
		i__2 = k - 1;
#line 246 "chpgst.f"
		chpr2_(uplo, &i__2, &c_b1, &ap[k1], &c__1, &bp[k1], &c__1, &
			ap[1], (ftnlen)1);
#line 248 "chpgst.f"
		i__2 = k - 1;
#line 248 "chpgst.f"
		caxpy_(&i__2, &ct, &bp[k1], &c__1, &ap[k1], &c__1);
#line 249 "chpgst.f"
		i__2 = k - 1;
#line 249 "chpgst.f"
		csscal_(&i__2, &bkk, &ap[k1], &c__1);
#line 250 "chpgst.f"
		i__2 = kk;
/* Computing 2nd power */
#line 250 "chpgst.f"
		d__2 = bkk;
#line 250 "chpgst.f"
		d__1 = akk * (d__2 * d__2);
#line 250 "chpgst.f"
		ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 251 "chpgst.f"
/* L30: */
#line 251 "chpgst.f"
	    }
#line 252 "chpgst.f"
	} else {

/*           Compute L**H *A*L */

/*           JJ and J1J1 are the indices of A(j,j) and A(j+1,j+1) */

#line 258 "chpgst.f"
	    jj = 1;
#line 259 "chpgst.f"
	    i__1 = *n;
#line 259 "chpgst.f"
	    for (j = 1; j <= i__1; ++j) {
#line 260 "chpgst.f"
		j1j1 = jj + *n - j + 1;

/*              Compute the j-th column of the lower triangle of A */

#line 264 "chpgst.f"
		i__2 = jj;
#line 264 "chpgst.f"
		ajj = ap[i__2].r;
#line 265 "chpgst.f"
		i__2 = jj;
#line 265 "chpgst.f"
		bjj = bp[i__2].r;
#line 266 "chpgst.f"
		i__2 = jj;
#line 266 "chpgst.f"
		d__1 = ajj * bjj;
#line 266 "chpgst.f"
		i__3 = *n - j;
#line 266 "chpgst.f"
		cdotc_(&z__2, &i__3, &ap[jj + 1], &c__1, &bp[jj + 1], &c__1);
#line 266 "chpgst.f"
		z__1.r = d__1 + z__2.r, z__1.i = z__2.i;
#line 266 "chpgst.f"
		ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
#line 268 "chpgst.f"
		i__2 = *n - j;
#line 268 "chpgst.f"
		csscal_(&i__2, &bjj, &ap[jj + 1], &c__1);
#line 269 "chpgst.f"
		i__2 = *n - j;
#line 269 "chpgst.f"
		chpmv_(uplo, &i__2, &c_b1, &ap[j1j1], &bp[jj + 1], &c__1, &
			c_b1, &ap[jj + 1], &c__1, (ftnlen)1);
#line 271 "chpgst.f"
		i__2 = *n - j + 1;
#line 271 "chpgst.f"
		ctpmv_(uplo, "Conjugate transpose", "Non-unit", &i__2, &bp[jj]
			, &ap[jj], &c__1, (ftnlen)1, (ftnlen)19, (ftnlen)8);
#line 273 "chpgst.f"
		jj = j1j1;
#line 274 "chpgst.f"
/* L40: */
#line 274 "chpgst.f"
	    }
#line 275 "chpgst.f"
	}
#line 276 "chpgst.f"
    }
#line 277 "chpgst.f"
    return 0;

/*     End of CHPGST */

} /* chpgst_ */

