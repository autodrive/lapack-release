#line 1 "zlaqps.f"
/* zlaqps.f -- translated by f2c (version 20100827).
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

#line 1 "zlaqps.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZLAQPS computes a step of QR factorization with column pivoting of a real m-by-n matrix A by us
ing BLAS level 3. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAQPS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqps.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqps.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqps.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAQPS( M, N, OFFSET, NB, KB, A, LDA, JPVT, TAU, VN1, */
/*                          VN2, AUXV, F, LDF ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            KB, LDA, LDF, M, N, NB, OFFSET */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            JPVT( * ) */
/*       DOUBLE PRECISION   VN1( * ), VN2( * ) */
/*       COMPLEX*16         A( LDA, * ), AUXV( * ), F( LDF, * ), TAU( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAQPS computes a step of QR factorization with column pivoting */
/* > of a complex M-by-N matrix A by using Blas-3.  It tries to factorize */
/* > NB columns from A starting from the row OFFSET+1, and updates all */
/* > of the matrix with Blas-3 xGEMM. */
/* > */
/* > In some cases, due to catastrophic cancellations, it cannot */
/* > factorize NB columns.  Hence, the actual number of factorized */
/* > columns is returned in KB. */
/* > */
/* > Block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A. N >= 0 */
/* > \endverbatim */
/* > */
/* > \param[in] OFFSET */
/* > \verbatim */
/* >          OFFSET is INTEGER */
/* >          The number of rows of A that have been factorized in */
/* >          previous steps. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The number of columns to factorize. */
/* > \endverbatim */
/* > */
/* > \param[out] KB */
/* > \verbatim */
/* >          KB is INTEGER */
/* >          The number of columns actually factorized. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, block A(OFFSET+1:M,1:KB) is the triangular */
/* >          factor obtained and block A(1:OFFSET,1:N) has been */
/* >          accordingly pivoted, but no factorized. */
/* >          The rest of the matrix, block A(OFFSET+1:M,KB+1:N) has */
/* >          been updated. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] JPVT */
/* > \verbatim */
/* >          JPVT is INTEGER array, dimension (N) */
/* >          JPVT(I) = K <==> Column K of the full matrix A has been */
/* >          permuted into position I in AP. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX*16 array, dimension (KB) */
/* >          The scalar factors of the elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VN1 */
/* > \verbatim */
/* >          VN1 is DOUBLE PRECISION array, dimension (N) */
/* >          The vector with the partial column norms. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VN2 */
/* > \verbatim */
/* >          VN2 is DOUBLE PRECISION array, dimension (N) */
/* >          The vector with the exact column norms. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AUXV */
/* > \verbatim */
/* >          AUXV is COMPLEX*16 array, dimension (NB) */
/* >          Auxiliar vector. */
/* > \endverbatim */
/* > */
/* > \param[in,out] F */
/* > \verbatim */
/* >          F is COMPLEX*16 array, dimension (LDF,NB) */
/* >          Matrix F**H = L * Y**H * A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDF */
/* > \verbatim */
/* >          LDF is INTEGER */
/* >          The leading dimension of the array F. LDF >= max(1,N). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complex16OTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain */
/* >    X. Sun, Computer Science Dept., Duke University, USA */
/* > \n */
/* >  Partial column norm updating strategy modified on April 2011 */
/* >    Z. Drmac and Z. Bujanovic, Dept. of Mathematics, */
/* >    University of Zagreb, Croatia. */

/* > \par References: */
/*  ================ */
/* > */
/* > LAPACK Working Note 176 */

/* > \htmlonly */
/* > <a href="http://www.netlib.org/lapack/lawnspdf/lawn176.pdf">[PDF]</a> */
/* > \endhtmlonly */

/*  ===================================================================== */
/* Subroutine */ int zlaqps_(integer *m, integer *n, integer *offset, integer 
	*nb, integer *kb, doublecomplex *a, integer *lda, integer *jpvt, 
	doublecomplex *tau, doublereal *vn1, doublereal *vn2, doublecomplex *
	auxv, doublecomplex *f, integer *ldf)
{
    /* System generated locals */
    integer a_dim1, a_offset, f_dim1, f_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    double sqrt(doublereal);
    void d_cnjg(doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *);
    integer i_dnnt(doublereal *);

    /* Local variables */
    static integer j, k, rk;
    static doublecomplex akk;
    static integer pvt;
    static doublereal temp, temp2, tol3z;
    static integer itemp;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    zswap_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *);
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *), dlamch_(
	    char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    static integer lsticc;
    extern /* Subroutine */ int zlarfg_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *);
    static integer lastrk;


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 221 "zlaqps.f"
    /* Parameter adjustments */
#line 221 "zlaqps.f"
    a_dim1 = *lda;
#line 221 "zlaqps.f"
    a_offset = 1 + a_dim1;
#line 221 "zlaqps.f"
    a -= a_offset;
#line 221 "zlaqps.f"
    --jpvt;
#line 221 "zlaqps.f"
    --tau;
#line 221 "zlaqps.f"
    --vn1;
#line 221 "zlaqps.f"
    --vn2;
#line 221 "zlaqps.f"
    --auxv;
#line 221 "zlaqps.f"
    f_dim1 = *ldf;
#line 221 "zlaqps.f"
    f_offset = 1 + f_dim1;
#line 221 "zlaqps.f"
    f -= f_offset;
#line 221 "zlaqps.f"

#line 221 "zlaqps.f"
    /* Function Body */
/* Computing MIN */
#line 221 "zlaqps.f"
    i__1 = *m, i__2 = *n + *offset;
#line 221 "zlaqps.f"
    lastrk = min(i__1,i__2);
#line 222 "zlaqps.f"
    lsticc = 0;
#line 223 "zlaqps.f"
    k = 0;
#line 224 "zlaqps.f"
    tol3z = sqrt(dlamch_("Epsilon", (ftnlen)7));

/*     Beginning of while loop. */

#line 228 "zlaqps.f"
L10:
#line 229 "zlaqps.f"
    if (k < *nb && lsticc == 0) {
#line 230 "zlaqps.f"
	++k;
#line 231 "zlaqps.f"
	rk = *offset + k;

/*        Determine ith pivot column and swap if necessary */

#line 235 "zlaqps.f"
	i__1 = *n - k + 1;
#line 235 "zlaqps.f"
	pvt = k - 1 + idamax_(&i__1, &vn1[k], &c__1);
#line 236 "zlaqps.f"
	if (pvt != k) {
#line 237 "zlaqps.f"
	    zswap_(m, &a[pvt * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
#line 238 "zlaqps.f"
	    i__1 = k - 1;
#line 238 "zlaqps.f"
	    zswap_(&i__1, &f[pvt + f_dim1], ldf, &f[k + f_dim1], ldf);
#line 239 "zlaqps.f"
	    itemp = jpvt[pvt];
#line 240 "zlaqps.f"
	    jpvt[pvt] = jpvt[k];
#line 241 "zlaqps.f"
	    jpvt[k] = itemp;
#line 242 "zlaqps.f"
	    vn1[pvt] = vn1[k];
#line 243 "zlaqps.f"
	    vn2[pvt] = vn2[k];
#line 244 "zlaqps.f"
	}

/*        Apply previous Householder reflectors to column K: */
/*        A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)**H. */

#line 249 "zlaqps.f"
	if (k > 1) {
#line 250 "zlaqps.f"
	    i__1 = k - 1;
#line 250 "zlaqps.f"
	    for (j = 1; j <= i__1; ++j) {
#line 251 "zlaqps.f"
		i__2 = k + j * f_dim1;
#line 251 "zlaqps.f"
		d_cnjg(&z__1, &f[k + j * f_dim1]);
#line 251 "zlaqps.f"
		f[i__2].r = z__1.r, f[i__2].i = z__1.i;
#line 252 "zlaqps.f"
/* L20: */
#line 252 "zlaqps.f"
	    }
#line 253 "zlaqps.f"
	    i__1 = *m - rk + 1;
#line 253 "zlaqps.f"
	    i__2 = k - 1;
#line 253 "zlaqps.f"
	    z__1.r = -1., z__1.i = -0.;
#line 253 "zlaqps.f"
	    zgemv_("No transpose", &i__1, &i__2, &z__1, &a[rk + a_dim1], lda, 
		    &f[k + f_dim1], ldf, &c_b2, &a[rk + k * a_dim1], &c__1, (
		    ftnlen)12);
#line 255 "zlaqps.f"
	    i__1 = k - 1;
#line 255 "zlaqps.f"
	    for (j = 1; j <= i__1; ++j) {
#line 256 "zlaqps.f"
		i__2 = k + j * f_dim1;
#line 256 "zlaqps.f"
		d_cnjg(&z__1, &f[k + j * f_dim1]);
#line 256 "zlaqps.f"
		f[i__2].r = z__1.r, f[i__2].i = z__1.i;
#line 257 "zlaqps.f"
/* L30: */
#line 257 "zlaqps.f"
	    }
#line 258 "zlaqps.f"
	}

/*        Generate elementary reflector H(k). */

#line 262 "zlaqps.f"
	if (rk < *m) {
#line 263 "zlaqps.f"
	    i__1 = *m - rk + 1;
#line 263 "zlaqps.f"
	    zlarfg_(&i__1, &a[rk + k * a_dim1], &a[rk + 1 + k * a_dim1], &
		    c__1, &tau[k]);
#line 264 "zlaqps.f"
	} else {
#line 265 "zlaqps.f"
	    zlarfg_(&c__1, &a[rk + k * a_dim1], &a[rk + k * a_dim1], &c__1, &
		    tau[k]);
#line 266 "zlaqps.f"
	}

#line 268 "zlaqps.f"
	i__1 = rk + k * a_dim1;
#line 268 "zlaqps.f"
	akk.r = a[i__1].r, akk.i = a[i__1].i;
#line 269 "zlaqps.f"
	i__1 = rk + k * a_dim1;
#line 269 "zlaqps.f"
	a[i__1].r = 1., a[i__1].i = 0.;

/*        Compute Kth column of F: */

/*        Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)**H*A(RK:M,K). */

#line 275 "zlaqps.f"
	if (k < *n) {
#line 276 "zlaqps.f"
	    i__1 = *m - rk + 1;
#line 276 "zlaqps.f"
	    i__2 = *n - k;
#line 276 "zlaqps.f"
	    zgemv_("Conjugate transpose", &i__1, &i__2, &tau[k], &a[rk + (k + 
		    1) * a_dim1], lda, &a[rk + k * a_dim1], &c__1, &c_b1, &f[
		    k + 1 + k * f_dim1], &c__1, (ftnlen)19);
#line 279 "zlaqps.f"
	}

/*        Padding F(1:K,K) with zeros. */

#line 283 "zlaqps.f"
	i__1 = k;
#line 283 "zlaqps.f"
	for (j = 1; j <= i__1; ++j) {
#line 284 "zlaqps.f"
	    i__2 = j + k * f_dim1;
#line 284 "zlaqps.f"
	    f[i__2].r = 0., f[i__2].i = 0.;
#line 285 "zlaqps.f"
/* L40: */
#line 285 "zlaqps.f"
	}

/*        Incremental updating of F: */
/*        F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)**H */
/*                    *A(RK:M,K). */

#line 291 "zlaqps.f"
	if (k > 1) {
#line 292 "zlaqps.f"
	    i__1 = *m - rk + 1;
#line 292 "zlaqps.f"
	    i__2 = k - 1;
#line 292 "zlaqps.f"
	    i__3 = k;
#line 292 "zlaqps.f"
	    z__1.r = -tau[i__3].r, z__1.i = -tau[i__3].i;
#line 292 "zlaqps.f"
	    zgemv_("Conjugate transpose", &i__1, &i__2, &z__1, &a[rk + a_dim1]
		    , lda, &a[rk + k * a_dim1], &c__1, &c_b1, &auxv[1], &c__1,
		     (ftnlen)19);

#line 296 "zlaqps.f"
	    i__1 = k - 1;
#line 296 "zlaqps.f"
	    zgemv_("No transpose", n, &i__1, &c_b2, &f[f_dim1 + 1], ldf, &
		    auxv[1], &c__1, &c_b2, &f[k * f_dim1 + 1], &c__1, (ftnlen)
		    12);
#line 298 "zlaqps.f"
	}

/*        Update the current row of A: */
/*        A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)**H. */

#line 303 "zlaqps.f"
	if (k < *n) {
#line 304 "zlaqps.f"
	    i__1 = *n - k;
#line 304 "zlaqps.f"
	    z__1.r = -1., z__1.i = -0.;
#line 304 "zlaqps.f"
	    zgemm_("No transpose", "Conjugate transpose", &c__1, &i__1, &k, &
		    z__1, &a[rk + a_dim1], lda, &f[k + 1 + f_dim1], ldf, &
		    c_b2, &a[rk + (k + 1) * a_dim1], lda, (ftnlen)12, (ftnlen)
		    19);
#line 307 "zlaqps.f"
	}

/*        Update partial column norms. */

#line 311 "zlaqps.f"
	if (rk < lastrk) {
#line 312 "zlaqps.f"
	    i__1 = *n;
#line 312 "zlaqps.f"
	    for (j = k + 1; j <= i__1; ++j) {
#line 313 "zlaqps.f"
		if (vn1[j] != 0.) {

/*                 NOTE: The following 4 lines follow from the analysis in */
/*                 Lapack Working Note 176. */

#line 318 "zlaqps.f"
		    temp = z_abs(&a[rk + j * a_dim1]) / vn1[j];
/* Computing MAX */
#line 319 "zlaqps.f"
		    d__1 = 0., d__2 = (temp + 1.) * (1. - temp);
#line 319 "zlaqps.f"
		    temp = max(d__1,d__2);
/* Computing 2nd power */
#line 320 "zlaqps.f"
		    d__1 = vn1[j] / vn2[j];
#line 320 "zlaqps.f"
		    temp2 = temp * (d__1 * d__1);
#line 321 "zlaqps.f"
		    if (temp2 <= tol3z) {
#line 322 "zlaqps.f"
			vn2[j] = (doublereal) lsticc;
#line 323 "zlaqps.f"
			lsticc = j;
#line 324 "zlaqps.f"
		    } else {
#line 325 "zlaqps.f"
			vn1[j] *= sqrt(temp);
#line 326 "zlaqps.f"
		    }
#line 327 "zlaqps.f"
		}
#line 328 "zlaqps.f"
/* L50: */
#line 328 "zlaqps.f"
	    }
#line 329 "zlaqps.f"
	}

#line 331 "zlaqps.f"
	i__1 = rk + k * a_dim1;
#line 331 "zlaqps.f"
	a[i__1].r = akk.r, a[i__1].i = akk.i;

/*        End of while loop. */

#line 335 "zlaqps.f"
	goto L10;
#line 336 "zlaqps.f"
    }
#line 337 "zlaqps.f"
    *kb = k;
#line 338 "zlaqps.f"
    rk = *offset + *kb;

/*     Apply the block reflector to the rest of the matrix: */
/*     A(OFFSET+KB+1:M,KB+1:N) := A(OFFSET+KB+1:M,KB+1:N) - */
/*                         A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)**H. */

/* Computing MIN */
#line 344 "zlaqps.f"
    i__1 = *n, i__2 = *m - *offset;
#line 344 "zlaqps.f"
    if (*kb < min(i__1,i__2)) {
#line 345 "zlaqps.f"
	i__1 = *m - rk;
#line 345 "zlaqps.f"
	i__2 = *n - *kb;
#line 345 "zlaqps.f"
	z__1.r = -1., z__1.i = -0.;
#line 345 "zlaqps.f"
	zgemm_("No transpose", "Conjugate transpose", &i__1, &i__2, kb, &z__1,
		 &a[rk + 1 + a_dim1], lda, &f[*kb + 1 + f_dim1], ldf, &c_b2, &
		a[rk + 1 + (*kb + 1) * a_dim1], lda, (ftnlen)12, (ftnlen)19);
#line 348 "zlaqps.f"
    }

/*     Recomputation of difficult columns. */

#line 352 "zlaqps.f"
L60:
#line 353 "zlaqps.f"
    if (lsticc > 0) {
#line 354 "zlaqps.f"
	itemp = i_dnnt(&vn2[lsticc]);
#line 355 "zlaqps.f"
	i__1 = *m - rk;
#line 355 "zlaqps.f"
	vn1[lsticc] = dznrm2_(&i__1, &a[rk + 1 + lsticc * a_dim1], &c__1);

/*        NOTE: The computation of VN1( LSTICC ) relies on the fact that */
/*        SNRM2 does not fail on vectors with norm below the value of */
/*        SQRT(DLAMCH('S')) */

#line 361 "zlaqps.f"
	vn2[lsticc] = vn1[lsticc];
#line 362 "zlaqps.f"
	lsticc = itemp;
#line 363 "zlaqps.f"
	goto L60;
#line 364 "zlaqps.f"
    }

#line 366 "zlaqps.f"
    return 0;

/*     End of ZLAQPS */

} /* zlaqps_ */

