#line 1 "slaqps.f"
/* slaqps.f -- translated by f2c (version 20100827).
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

#line 1 "slaqps.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b8 = -1.;
static doublereal c_b9 = 1.;
static doublereal c_b16 = 0.;

/* > \brief \b SLAQPS computes a step of QR factorization with column pivoting of a real m-by-n matrix A by us
ing BLAS level 3. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAQPS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqps.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqps.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqps.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAQPS( M, N, OFFSET, NB, KB, A, LDA, JPVT, TAU, VN1, */
/*                          VN2, AUXV, F, LDF ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            KB, LDA, LDF, M, N, NB, OFFSET */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            JPVT( * ) */
/*       REAL               A( LDA, * ), AUXV( * ), F( LDF, * ), TAU( * ), */
/*      $                   VN1( * ), VN2( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAQPS computes a step of QR factorization with column pivoting */
/* > of a real M-by-N matrix A by using Blas-3.  It tries to factorize */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          TAU is REAL array, dimension (KB) */
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
/* >          AUXV is REAL array, dimension (NB) */
/* >          Auxiliar vector. */
/* > \endverbatim */
/* > */
/* > \param[in,out] F */
/* > \verbatim */
/* >          F is REAL array, dimension (LDF,NB) */
/* >          Matrix F**T = L*Y**T*A. */
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

/* > \ingroup realOTHERauxiliary */

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
/* Subroutine */ int slaqps_(integer *m, integer *n, integer *offset, integer 
	*nb, integer *kb, doublereal *a, integer *lda, integer *jpvt, 
	doublereal *tau, doublereal *vn1, doublereal *vn2, doublereal *auxv, 
	doublereal *f, integer *ldf)
{
    /* System generated locals */
    integer a_dim1, a_offset, f_dim1, f_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);
    integer i_dnnt(doublereal *);

    /* Local variables */
    static integer j, k, rk;
    static doublereal akk;
    static integer pvt;
    static doublereal temp, temp2;
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    static doublereal tol3z;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer itemp;
    extern /* Subroutine */ int sgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), sswap_(integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int slarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
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

#line 218 "slaqps.f"
    /* Parameter adjustments */
#line 218 "slaqps.f"
    a_dim1 = *lda;
#line 218 "slaqps.f"
    a_offset = 1 + a_dim1;
#line 218 "slaqps.f"
    a -= a_offset;
#line 218 "slaqps.f"
    --jpvt;
#line 218 "slaqps.f"
    --tau;
#line 218 "slaqps.f"
    --vn1;
#line 218 "slaqps.f"
    --vn2;
#line 218 "slaqps.f"
    --auxv;
#line 218 "slaqps.f"
    f_dim1 = *ldf;
#line 218 "slaqps.f"
    f_offset = 1 + f_dim1;
#line 218 "slaqps.f"
    f -= f_offset;
#line 218 "slaqps.f"

#line 218 "slaqps.f"
    /* Function Body */
/* Computing MIN */
#line 218 "slaqps.f"
    i__1 = *m, i__2 = *n + *offset;
#line 218 "slaqps.f"
    lastrk = min(i__1,i__2);
#line 219 "slaqps.f"
    lsticc = 0;
#line 220 "slaqps.f"
    k = 0;
#line 221 "slaqps.f"
    tol3z = sqrt(slamch_("Epsilon", (ftnlen)7));

/*     Beginning of while loop. */

#line 225 "slaqps.f"
L10:
#line 226 "slaqps.f"
    if (k < *nb && lsticc == 0) {
#line 227 "slaqps.f"
	++k;
#line 228 "slaqps.f"
	rk = *offset + k;

/*        Determine ith pivot column and swap if necessary */

#line 232 "slaqps.f"
	i__1 = *n - k + 1;
#line 232 "slaqps.f"
	pvt = k - 1 + isamax_(&i__1, &vn1[k], &c__1);
#line 233 "slaqps.f"
	if (pvt != k) {
#line 234 "slaqps.f"
	    sswap_(m, &a[pvt * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
#line 235 "slaqps.f"
	    i__1 = k - 1;
#line 235 "slaqps.f"
	    sswap_(&i__1, &f[pvt + f_dim1], ldf, &f[k + f_dim1], ldf);
#line 236 "slaqps.f"
	    itemp = jpvt[pvt];
#line 237 "slaqps.f"
	    jpvt[pvt] = jpvt[k];
#line 238 "slaqps.f"
	    jpvt[k] = itemp;
#line 239 "slaqps.f"
	    vn1[pvt] = vn1[k];
#line 240 "slaqps.f"
	    vn2[pvt] = vn2[k];
#line 241 "slaqps.f"
	}

/*        Apply previous Householder reflectors to column K: */
/*        A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)**T. */

#line 246 "slaqps.f"
	if (k > 1) {
#line 247 "slaqps.f"
	    i__1 = *m - rk + 1;
#line 247 "slaqps.f"
	    i__2 = k - 1;
#line 247 "slaqps.f"
	    sgemv_("No transpose", &i__1, &i__2, &c_b8, &a[rk + a_dim1], lda, 
		    &f[k + f_dim1], ldf, &c_b9, &a[rk + k * a_dim1], &c__1, (
		    ftnlen)12);
#line 249 "slaqps.f"
	}

/*        Generate elementary reflector H(k). */

#line 253 "slaqps.f"
	if (rk < *m) {
#line 254 "slaqps.f"
	    i__1 = *m - rk + 1;
#line 254 "slaqps.f"
	    slarfg_(&i__1, &a[rk + k * a_dim1], &a[rk + 1 + k * a_dim1], &
		    c__1, &tau[k]);
#line 255 "slaqps.f"
	} else {
#line 256 "slaqps.f"
	    slarfg_(&c__1, &a[rk + k * a_dim1], &a[rk + k * a_dim1], &c__1, &
		    tau[k]);
#line 257 "slaqps.f"
	}

#line 259 "slaqps.f"
	akk = a[rk + k * a_dim1];
#line 260 "slaqps.f"
	a[rk + k * a_dim1] = 1.;

/*        Compute Kth column of F: */

/*        Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)**T*A(RK:M,K). */

#line 266 "slaqps.f"
	if (k < *n) {
#line 267 "slaqps.f"
	    i__1 = *m - rk + 1;
#line 267 "slaqps.f"
	    i__2 = *n - k;
#line 267 "slaqps.f"
	    sgemv_("Transpose", &i__1, &i__2, &tau[k], &a[rk + (k + 1) * 
		    a_dim1], lda, &a[rk + k * a_dim1], &c__1, &c_b16, &f[k + 
		    1 + k * f_dim1], &c__1, (ftnlen)9);
#line 270 "slaqps.f"
	}

/*        Padding F(1:K,K) with zeros. */

#line 274 "slaqps.f"
	i__1 = k;
#line 274 "slaqps.f"
	for (j = 1; j <= i__1; ++j) {
#line 275 "slaqps.f"
	    f[j + k * f_dim1] = 0.;
#line 276 "slaqps.f"
/* L20: */
#line 276 "slaqps.f"
	}

/*        Incremental updating of F: */
/*        F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)**T */
/*                    *A(RK:M,K). */

#line 282 "slaqps.f"
	if (k > 1) {
#line 283 "slaqps.f"
	    i__1 = *m - rk + 1;
#line 283 "slaqps.f"
	    i__2 = k - 1;
#line 283 "slaqps.f"
	    d__1 = -tau[k];
#line 283 "slaqps.f"
	    sgemv_("Transpose", &i__1, &i__2, &d__1, &a[rk + a_dim1], lda, &a[
		    rk + k * a_dim1], &c__1, &c_b16, &auxv[1], &c__1, (ftnlen)
		    9);

#line 286 "slaqps.f"
	    i__1 = k - 1;
#line 286 "slaqps.f"
	    sgemv_("No transpose", n, &i__1, &c_b9, &f[f_dim1 + 1], ldf, &
		    auxv[1], &c__1, &c_b9, &f[k * f_dim1 + 1], &c__1, (ftnlen)
		    12);
#line 288 "slaqps.f"
	}

/*        Update the current row of A: */
/*        A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)**T. */

#line 293 "slaqps.f"
	if (k < *n) {
#line 294 "slaqps.f"
	    i__1 = *n - k;
#line 294 "slaqps.f"
	    sgemv_("No transpose", &i__1, &k, &c_b8, &f[k + 1 + f_dim1], ldf, 
		    &a[rk + a_dim1], lda, &c_b9, &a[rk + (k + 1) * a_dim1], 
		    lda, (ftnlen)12);
#line 296 "slaqps.f"
	}

/*        Update partial column norms. */

#line 300 "slaqps.f"
	if (rk < lastrk) {
#line 301 "slaqps.f"
	    i__1 = *n;
#line 301 "slaqps.f"
	    for (j = k + 1; j <= i__1; ++j) {
#line 302 "slaqps.f"
		if (vn1[j] != 0.) {

/*                 NOTE: The following 4 lines follow from the analysis in */
/*                 Lapack Working Note 176. */

#line 307 "slaqps.f"
		    temp = (d__1 = a[rk + j * a_dim1], abs(d__1)) / vn1[j];
/* Computing MAX */
#line 308 "slaqps.f"
		    d__1 = 0., d__2 = (temp + 1.) * (1. - temp);
#line 308 "slaqps.f"
		    temp = max(d__1,d__2);
/* Computing 2nd power */
#line 309 "slaqps.f"
		    d__1 = vn1[j] / vn2[j];
#line 309 "slaqps.f"
		    temp2 = temp * (d__1 * d__1);
#line 310 "slaqps.f"
		    if (temp2 <= tol3z) {
#line 311 "slaqps.f"
			vn2[j] = (doublereal) lsticc;
#line 312 "slaqps.f"
			lsticc = j;
#line 313 "slaqps.f"
		    } else {
#line 314 "slaqps.f"
			vn1[j] *= sqrt(temp);
#line 315 "slaqps.f"
		    }
#line 316 "slaqps.f"
		}
#line 317 "slaqps.f"
/* L30: */
#line 317 "slaqps.f"
	    }
#line 318 "slaqps.f"
	}

#line 320 "slaqps.f"
	a[rk + k * a_dim1] = akk;

/*        End of while loop. */

#line 324 "slaqps.f"
	goto L10;
#line 325 "slaqps.f"
    }
#line 326 "slaqps.f"
    *kb = k;
#line 327 "slaqps.f"
    rk = *offset + *kb;

/*     Apply the block reflector to the rest of the matrix: */
/*     A(OFFSET+KB+1:M,KB+1:N) := A(OFFSET+KB+1:M,KB+1:N) - */
/*                         A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)**T. */

/* Computing MIN */
#line 333 "slaqps.f"
    i__1 = *n, i__2 = *m - *offset;
#line 333 "slaqps.f"
    if (*kb < min(i__1,i__2)) {
#line 334 "slaqps.f"
	i__1 = *m - rk;
#line 334 "slaqps.f"
	i__2 = *n - *kb;
#line 334 "slaqps.f"
	sgemm_("No transpose", "Transpose", &i__1, &i__2, kb, &c_b8, &a[rk + 
		1 + a_dim1], lda, &f[*kb + 1 + f_dim1], ldf, &c_b9, &a[rk + 1 
		+ (*kb + 1) * a_dim1], lda, (ftnlen)12, (ftnlen)9);
#line 337 "slaqps.f"
    }

/*     Recomputation of difficult columns. */

#line 341 "slaqps.f"
L40:
#line 342 "slaqps.f"
    if (lsticc > 0) {
#line 343 "slaqps.f"
	itemp = i_dnnt(&vn2[lsticc]);
#line 344 "slaqps.f"
	i__1 = *m - rk;
#line 344 "slaqps.f"
	vn1[lsticc] = snrm2_(&i__1, &a[rk + 1 + lsticc * a_dim1], &c__1);

/*        NOTE: The computation of VN1( LSTICC ) relies on the fact that */
/*        SNRM2 does not fail on vectors with norm below the value of */
/*        SQRT(DLAMCH('S')) */

#line 350 "slaqps.f"
	vn2[lsticc] = vn1[lsticc];
#line 351 "slaqps.f"
	lsticc = itemp;
#line 352 "slaqps.f"
	goto L40;
#line 353 "slaqps.f"
    }

#line 355 "slaqps.f"
    return 0;

/*     End of SLAQPS */

} /* slaqps_ */

