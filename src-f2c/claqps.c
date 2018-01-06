#line 1 "claqps.f"
/* claqps.f -- translated by f2c (version 20100827).
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

#line 1 "claqps.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CLAQPS computes a step of QR factorization with column pivoting of a real m-by-n matrix A by us
ing BLAS level 3. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLAQPS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqps.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqps.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqps.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLAQPS( M, N, OFFSET, NB, KB, A, LDA, JPVT, TAU, VN1, */
/*                          VN2, AUXV, F, LDF ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            KB, LDA, LDF, M, N, NB, OFFSET */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            JPVT( * ) */
/*       REAL               VN1( * ), VN2( * ) */
/*       COMPLEX            A( LDA, * ), AUXV( * ), F( LDF, * ), TAU( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAQPS computes a step of QR factorization with column pivoting */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          TAU is COMPLEX array, dimension (KB) */
/* >          The scalar factors of the elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VN1 */
/* > \verbatim */
/* >          VN1 is REAL array, dimension (N) */
/* >          The vector with the partial column norms. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VN2 */
/* > \verbatim */
/* >          VN2 is REAL array, dimension (N) */
/* >          The vector with the exact column norms. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AUXV */
/* > \verbatim */
/* >          AUXV is COMPLEX array, dimension (NB) */
/* >          Auxiliar vector. */
/* > \endverbatim */
/* > */
/* > \param[in,out] F */
/* > \verbatim */
/* >          F is COMPLEX array, dimension (LDF,NB) */
/* >          Matrix  F**H = L * Y**H * A. */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain */
/* >    X. Sun, Computer Science Dept., Duke University, USA */
/* > */
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
/* Subroutine */ int claqps_(integer *m, integer *n, integer *offset, integer 
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
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    cswap_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *);
    static integer itemp;
    extern doublereal scnrm2_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int clarfg_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *);
    extern doublereal slamch_(char *, ftnlen);
    static integer lsticc;
    extern integer isamax_(integer *, doublereal *, integer *);
    static integer lastrk;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
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

#line 222 "claqps.f"
    /* Parameter adjustments */
#line 222 "claqps.f"
    a_dim1 = *lda;
#line 222 "claqps.f"
    a_offset = 1 + a_dim1;
#line 222 "claqps.f"
    a -= a_offset;
#line 222 "claqps.f"
    --jpvt;
#line 222 "claqps.f"
    --tau;
#line 222 "claqps.f"
    --vn1;
#line 222 "claqps.f"
    --vn2;
#line 222 "claqps.f"
    --auxv;
#line 222 "claqps.f"
    f_dim1 = *ldf;
#line 222 "claqps.f"
    f_offset = 1 + f_dim1;
#line 222 "claqps.f"
    f -= f_offset;
#line 222 "claqps.f"

#line 222 "claqps.f"
    /* Function Body */
/* Computing MIN */
#line 222 "claqps.f"
    i__1 = *m, i__2 = *n + *offset;
#line 222 "claqps.f"
    lastrk = min(i__1,i__2);
#line 223 "claqps.f"
    lsticc = 0;
#line 224 "claqps.f"
    k = 0;
#line 225 "claqps.f"
    tol3z = sqrt(slamch_("Epsilon", (ftnlen)7));

/*     Beginning of while loop. */

#line 229 "claqps.f"
L10:
#line 230 "claqps.f"
    if (k < *nb && lsticc == 0) {
#line 231 "claqps.f"
	++k;
#line 232 "claqps.f"
	rk = *offset + k;

/*        Determine ith pivot column and swap if necessary */

#line 236 "claqps.f"
	i__1 = *n - k + 1;
#line 236 "claqps.f"
	pvt = k - 1 + isamax_(&i__1, &vn1[k], &c__1);
#line 237 "claqps.f"
	if (pvt != k) {
#line 238 "claqps.f"
	    cswap_(m, &a[pvt * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
#line 239 "claqps.f"
	    i__1 = k - 1;
#line 239 "claqps.f"
	    cswap_(&i__1, &f[pvt + f_dim1], ldf, &f[k + f_dim1], ldf);
#line 240 "claqps.f"
	    itemp = jpvt[pvt];
#line 241 "claqps.f"
	    jpvt[pvt] = jpvt[k];
#line 242 "claqps.f"
	    jpvt[k] = itemp;
#line 243 "claqps.f"
	    vn1[pvt] = vn1[k];
#line 244 "claqps.f"
	    vn2[pvt] = vn2[k];
#line 245 "claqps.f"
	}

/*        Apply previous Householder reflectors to column K: */
/*        A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)**H. */

#line 250 "claqps.f"
	if (k > 1) {
#line 251 "claqps.f"
	    i__1 = k - 1;
#line 251 "claqps.f"
	    for (j = 1; j <= i__1; ++j) {
#line 252 "claqps.f"
		i__2 = k + j * f_dim1;
#line 252 "claqps.f"
		d_cnjg(&z__1, &f[k + j * f_dim1]);
#line 252 "claqps.f"
		f[i__2].r = z__1.r, f[i__2].i = z__1.i;
#line 253 "claqps.f"
/* L20: */
#line 253 "claqps.f"
	    }
#line 254 "claqps.f"
	    i__1 = *m - rk + 1;
#line 254 "claqps.f"
	    i__2 = k - 1;
#line 254 "claqps.f"
	    z__1.r = -1., z__1.i = -0.;
#line 254 "claqps.f"
	    cgemv_("No transpose", &i__1, &i__2, &z__1, &a[rk + a_dim1], lda, 
		    &f[k + f_dim1], ldf, &c_b2, &a[rk + k * a_dim1], &c__1, (
		    ftnlen)12);
#line 256 "claqps.f"
	    i__1 = k - 1;
#line 256 "claqps.f"
	    for (j = 1; j <= i__1; ++j) {
#line 257 "claqps.f"
		i__2 = k + j * f_dim1;
#line 257 "claqps.f"
		d_cnjg(&z__1, &f[k + j * f_dim1]);
#line 257 "claqps.f"
		f[i__2].r = z__1.r, f[i__2].i = z__1.i;
#line 258 "claqps.f"
/* L30: */
#line 258 "claqps.f"
	    }
#line 259 "claqps.f"
	}

/*        Generate elementary reflector H(k). */

#line 263 "claqps.f"
	if (rk < *m) {
#line 264 "claqps.f"
	    i__1 = *m - rk + 1;
#line 264 "claqps.f"
	    clarfg_(&i__1, &a[rk + k * a_dim1], &a[rk + 1 + k * a_dim1], &
		    c__1, &tau[k]);
#line 265 "claqps.f"
	} else {
#line 266 "claqps.f"
	    clarfg_(&c__1, &a[rk + k * a_dim1], &a[rk + k * a_dim1], &c__1, &
		    tau[k]);
#line 267 "claqps.f"
	}

#line 269 "claqps.f"
	i__1 = rk + k * a_dim1;
#line 269 "claqps.f"
	akk.r = a[i__1].r, akk.i = a[i__1].i;
#line 270 "claqps.f"
	i__1 = rk + k * a_dim1;
#line 270 "claqps.f"
	a[i__1].r = 1., a[i__1].i = 0.;

/*        Compute Kth column of F: */

/*        Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)**H*A(RK:M,K). */

#line 276 "claqps.f"
	if (k < *n) {
#line 277 "claqps.f"
	    i__1 = *m - rk + 1;
#line 277 "claqps.f"
	    i__2 = *n - k;
#line 277 "claqps.f"
	    cgemv_("Conjugate transpose", &i__1, &i__2, &tau[k], &a[rk + (k + 
		    1) * a_dim1], lda, &a[rk + k * a_dim1], &c__1, &c_b1, &f[
		    k + 1 + k * f_dim1], &c__1, (ftnlen)19);
#line 280 "claqps.f"
	}

/*        Padding F(1:K,K) with zeros. */

#line 284 "claqps.f"
	i__1 = k;
#line 284 "claqps.f"
	for (j = 1; j <= i__1; ++j) {
#line 285 "claqps.f"
	    i__2 = j + k * f_dim1;
#line 285 "claqps.f"
	    f[i__2].r = 0., f[i__2].i = 0.;
#line 286 "claqps.f"
/* L40: */
#line 286 "claqps.f"
	}

/*        Incremental updating of F: */
/*        F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)**H */
/*                    *A(RK:M,K). */

#line 292 "claqps.f"
	if (k > 1) {
#line 293 "claqps.f"
	    i__1 = *m - rk + 1;
#line 293 "claqps.f"
	    i__2 = k - 1;
#line 293 "claqps.f"
	    i__3 = k;
#line 293 "claqps.f"
	    z__1.r = -tau[i__3].r, z__1.i = -tau[i__3].i;
#line 293 "claqps.f"
	    cgemv_("Conjugate transpose", &i__1, &i__2, &z__1, &a[rk + a_dim1]
		    , lda, &a[rk + k * a_dim1], &c__1, &c_b1, &auxv[1], &c__1,
		     (ftnlen)19);

#line 297 "claqps.f"
	    i__1 = k - 1;
#line 297 "claqps.f"
	    cgemv_("No transpose", n, &i__1, &c_b2, &f[f_dim1 + 1], ldf, &
		    auxv[1], &c__1, &c_b2, &f[k * f_dim1 + 1], &c__1, (ftnlen)
		    12);
#line 299 "claqps.f"
	}

/*        Update the current row of A: */
/*        A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)**H. */

#line 304 "claqps.f"
	if (k < *n) {
#line 305 "claqps.f"
	    i__1 = *n - k;
#line 305 "claqps.f"
	    z__1.r = -1., z__1.i = -0.;
#line 305 "claqps.f"
	    cgemm_("No transpose", "Conjugate transpose", &c__1, &i__1, &k, &
		    z__1, &a[rk + a_dim1], lda, &f[k + 1 + f_dim1], ldf, &
		    c_b2, &a[rk + (k + 1) * a_dim1], lda, (ftnlen)12, (ftnlen)
		    19);
#line 308 "claqps.f"
	}

/*        Update partial column norms. */

#line 312 "claqps.f"
	if (rk < lastrk) {
#line 313 "claqps.f"
	    i__1 = *n;
#line 313 "claqps.f"
	    for (j = k + 1; j <= i__1; ++j) {
#line 314 "claqps.f"
		if (vn1[j] != 0.) {

/*                 NOTE: The following 4 lines follow from the analysis in */
/*                 Lapack Working Note 176. */

#line 319 "claqps.f"
		    temp = z_abs(&a[rk + j * a_dim1]) / vn1[j];
/* Computing MAX */
#line 320 "claqps.f"
		    d__1 = 0., d__2 = (temp + 1.) * (1. - temp);
#line 320 "claqps.f"
		    temp = max(d__1,d__2);
/* Computing 2nd power */
#line 321 "claqps.f"
		    d__1 = vn1[j] / vn2[j];
#line 321 "claqps.f"
		    temp2 = temp * (d__1 * d__1);
#line 322 "claqps.f"
		    if (temp2 <= tol3z) {
#line 323 "claqps.f"
			vn2[j] = (doublereal) lsticc;
#line 324 "claqps.f"
			lsticc = j;
#line 325 "claqps.f"
		    } else {
#line 326 "claqps.f"
			vn1[j] *= sqrt(temp);
#line 327 "claqps.f"
		    }
#line 328 "claqps.f"
		}
#line 329 "claqps.f"
/* L50: */
#line 329 "claqps.f"
	    }
#line 330 "claqps.f"
	}

#line 332 "claqps.f"
	i__1 = rk + k * a_dim1;
#line 332 "claqps.f"
	a[i__1].r = akk.r, a[i__1].i = akk.i;

/*        End of while loop. */

#line 336 "claqps.f"
	goto L10;
#line 337 "claqps.f"
    }
#line 338 "claqps.f"
    *kb = k;
#line 339 "claqps.f"
    rk = *offset + *kb;

/*     Apply the block reflector to the rest of the matrix: */
/*     A(OFFSET+KB+1:M,KB+1:N) := A(OFFSET+KB+1:M,KB+1:N) - */
/*                         A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)**H. */

/* Computing MIN */
#line 345 "claqps.f"
    i__1 = *n, i__2 = *m - *offset;
#line 345 "claqps.f"
    if (*kb < min(i__1,i__2)) {
#line 346 "claqps.f"
	i__1 = *m - rk;
#line 346 "claqps.f"
	i__2 = *n - *kb;
#line 346 "claqps.f"
	z__1.r = -1., z__1.i = -0.;
#line 346 "claqps.f"
	cgemm_("No transpose", "Conjugate transpose", &i__1, &i__2, kb, &z__1,
		 &a[rk + 1 + a_dim1], lda, &f[*kb + 1 + f_dim1], ldf, &c_b2, &
		a[rk + 1 + (*kb + 1) * a_dim1], lda, (ftnlen)12, (ftnlen)19);
#line 349 "claqps.f"
    }

/*     Recomputation of difficult columns. */

#line 353 "claqps.f"
L60:
#line 354 "claqps.f"
    if (lsticc > 0) {
#line 355 "claqps.f"
	itemp = i_dnnt(&vn2[lsticc]);
#line 356 "claqps.f"
	i__1 = *m - rk;
#line 356 "claqps.f"
	vn1[lsticc] = scnrm2_(&i__1, &a[rk + 1 + lsticc * a_dim1], &c__1);

/*        NOTE: The computation of VN1( LSTICC ) relies on the fact that */
/*        SNRM2 does not fail on vectors with norm below the value of */
/*        SQRT(DLAMCH('S')) */

#line 362 "claqps.f"
	vn2[lsticc] = vn1[lsticc];
#line 363 "claqps.f"
	lsticc = itemp;
#line 364 "claqps.f"
	goto L60;
#line 365 "claqps.f"
    }

#line 367 "claqps.f"
    return 0;

/*     End of CLAQPS */

} /* claqps_ */

