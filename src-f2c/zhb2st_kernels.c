#line 1 "zhb2st_kernels.f"
/* zhb2st_kernels.f -- translated by f2c (version 20100827).
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

#line 1 "zhb2st_kernels.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZHB2ST_KERNELS */

/*  @precisions fortran z -> s d c */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHB2ST_KERNELS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhb2st_
kernels.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhb2st_
kernels.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhb2st_
kernels.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE  ZHB2ST_KERNELS( UPLO, WANTZ, TTYPE, */
/*                                   ST, ED, SWEEP, N, NB, IB, */
/*                                   A, LDA, V, TAU, LDVT, WORK) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       LOGICAL            WANTZ */
/*       INTEGER            TTYPE, ST, ED, SWEEP, N, NB, IB, LDA, LDVT */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), V( * ), */
/*                          TAU( * ), WORK( * ) */

/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHB2ST_KERNELS is an internal routine used by the ZHETRD_HB2ST */
/* > subroutine. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > @param[in] n */
/* >          The order of the matrix A. */
/* > */
/* > @param[in] nb */
/* >          The size of the band. */
/* > */
/* > @param[in, out] A */
/* >          A pointer to the matrix A. */
/* > */
/* > @param[in] lda */
/* >          The leading dimension of the matrix A. */
/* > */
/* > @param[out] V */
/* >          COMPLEX*16 array, dimension 2*n if eigenvalues only are */
/* >          requested or to be queried for vectors. */
/* > */
/* > @param[out] TAU */
/* >          COMPLEX*16 array, dimension (2*n). */
/* >          The scalar factors of the Householder reflectors are stored */
/* >          in this array. */
/* > */
/* > @param[in] st */
/* >          internal parameter for indices. */
/* > */
/* > @param[in] ed */
/* >          internal parameter for indices. */
/* > */
/* > @param[in] sweep */
/* >          internal parameter for indices. */
/* > */
/* > @param[in] Vblksiz */
/* >          internal parameter for indices. */
/* > */
/* > @param[in] wantz */
/* >          logical which indicate if Eigenvalue are requested or both */
/* >          Eigenvalue/Eigenvectors. */
/* > */
/* > @param[in] work */
/* >          Workspace of size nb. */
/* > */
/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Implemented by Azzam Haidar. */
/* > */
/* >  All details are available on technical report, SC11, SC13 papers. */
/* > */
/* >  Azzam Haidar, Hatem Ltaief, and Jack Dongarra. */
/* >  Parallel reduction to condensed forms for symmetric eigenvalue problems */
/* >  using aggregated fine-grained and memory-aware kernels. In Proceedings */
/* >  of 2011 International Conference for High Performance Computing, */
/* >  Networking, Storage and Analysis (SC '11), New York, NY, USA, */
/* >  Article 8 , 11 pages. */
/* >  http://doi.acm.org/10.1145/2063384.2063394 */
/* > */
/* >  A. Haidar, J. Kurzak, P. Luszczek, 2013. */
/* >  An improved parallel singular value algorithm and its implementation */
/* >  for multicore hardware, In Proceedings of 2013 International Conference */
/* >  for High Performance Computing, Networking, Storage and Analysis (SC '13). */
/* >  Denver, Colorado, USA, 2013. */
/* >  Article 90, 12 pages. */
/* >  http://doi.acm.org/10.1145/2503210.2503292 */
/* > */
/* >  A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra. */
/* >  A novel hybrid CPU-GPU generalized eigensolver for electronic structure */
/* >  calculations based on fine-grained memory aware tasks. */
/* >  International Journal of High Performance Computing Applications. */
/* >  Volume 28 Issue 2, Pages 196-209, May 2014. */
/* >  http://hpc.sagepub.com/content/28/2/196 */
/* > */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zhb2st_kernels__(char *uplo, logical *wantz, integer *
	ttype, integer *st, integer *ed, integer *sweep, integer *n, integer *
	nb, integer *ib, doublecomplex *a, integer *lda, doublecomplex *v, 
	doublecomplex *tau, integer *ldvt, doublecomplex *work, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j1, j2, lm, ln;
    static doublecomplex ctmp;
    static integer dpos, vpos;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    static integer ajeter;
    extern /* Subroutine */ int zlarfg_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *);
    static integer ofdpos;
    extern /* Subroutine */ int zlarfx_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, ftnlen), zlarfy_(char *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen);
    static integer taupos;



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
/*     .. External Functions .. */
/*     .. */
/*     .. */
/*     .. Executable Statements .. */

#line 171 "zhb2st_kernels.f"
    /* Parameter adjustments */
#line 171 "zhb2st_kernels.f"
    a_dim1 = *lda;
#line 171 "zhb2st_kernels.f"
    a_offset = 1 + a_dim1;
#line 171 "zhb2st_kernels.f"
    a -= a_offset;
#line 171 "zhb2st_kernels.f"
    --v;
#line 171 "zhb2st_kernels.f"
    --tau;
#line 171 "zhb2st_kernels.f"
    --work;
#line 171 "zhb2st_kernels.f"

#line 171 "zhb2st_kernels.f"
    /* Function Body */
#line 171 "zhb2st_kernels.f"
    ajeter = *ib + *ldvt;
#line 172 "zhb2st_kernels.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 174 "zhb2st_kernels.f"
    if (upper) {
#line 175 "zhb2st_kernels.f"
	dpos = (*nb << 1) + 1;
#line 176 "zhb2st_kernels.f"
	ofdpos = *nb << 1;
#line 177 "zhb2st_kernels.f"
    } else {
#line 178 "zhb2st_kernels.f"
	dpos = 1;
#line 179 "zhb2st_kernels.f"
	ofdpos = 2;
#line 180 "zhb2st_kernels.f"
    }

/*     Upper case */

#line 185 "zhb2st_kernels.f"
    if (upper) {

#line 187 "zhb2st_kernels.f"
	if (*wantz) {
#line 188 "zhb2st_kernels.f"
	    vpos = (*sweep - 1) % 2 * *n + *st;
#line 189 "zhb2st_kernels.f"
	    taupos = (*sweep - 1) % 2 * *n + *st;
#line 190 "zhb2st_kernels.f"
	} else {
#line 191 "zhb2st_kernels.f"
	    vpos = (*sweep - 1) % 2 * *n + *st;
#line 192 "zhb2st_kernels.f"
	    taupos = (*sweep - 1) % 2 * *n + *st;
#line 193 "zhb2st_kernels.f"
	}

#line 195 "zhb2st_kernels.f"
	if (*ttype == 1) {
#line 196 "zhb2st_kernels.f"
	    lm = *ed - *st + 1;

#line 198 "zhb2st_kernels.f"
	    i__1 = vpos;
#line 198 "zhb2st_kernels.f"
	    v[i__1].r = 1., v[i__1].i = 0.;
#line 199 "zhb2st_kernels.f"
	    i__1 = lm - 1;
#line 199 "zhb2st_kernels.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 200 "zhb2st_kernels.f"
		i__2 = vpos + i__;
#line 200 "zhb2st_kernels.f"
		d_cnjg(&z__1, &a[ofdpos - i__ + (*st + i__) * a_dim1]);
#line 200 "zhb2st_kernels.f"
		v[i__2].r = z__1.r, v[i__2].i = z__1.i;
#line 201 "zhb2st_kernels.f"
		i__2 = ofdpos - i__ + (*st + i__) * a_dim1;
#line 201 "zhb2st_kernels.f"
		a[i__2].r = 0., a[i__2].i = 0.;
#line 202 "zhb2st_kernels.f"
/* L10: */
#line 202 "zhb2st_kernels.f"
	    }
#line 203 "zhb2st_kernels.f"
	    d_cnjg(&z__1, &a[ofdpos + *st * a_dim1]);
#line 203 "zhb2st_kernels.f"
	    ctmp.r = z__1.r, ctmp.i = z__1.i;
#line 204 "zhb2st_kernels.f"
	    zlarfg_(&lm, &ctmp, &v[vpos + 1], &c__1, &tau[taupos]);
#line 206 "zhb2st_kernels.f"
	    i__1 = ofdpos + *st * a_dim1;
#line 206 "zhb2st_kernels.f"
	    a[i__1].r = ctmp.r, a[i__1].i = ctmp.i;

#line 208 "zhb2st_kernels.f"
	    lm = *ed - *st + 1;
#line 209 "zhb2st_kernels.f"
	    d_cnjg(&z__1, &tau[taupos]);
#line 209 "zhb2st_kernels.f"
	    i__1 = *lda - 1;
#line 209 "zhb2st_kernels.f"
	    zlarfy_(uplo, &lm, &v[vpos], &c__1, &z__1, &a[dpos + *st * a_dim1]
		    , &i__1, &work[1], (ftnlen)1);
#line 212 "zhb2st_kernels.f"
	}

#line 214 "zhb2st_kernels.f"
	if (*ttype == 3) {

#line 216 "zhb2st_kernels.f"
	    lm = *ed - *st + 1;
#line 217 "zhb2st_kernels.f"
	    d_cnjg(&z__1, &tau[taupos]);
#line 217 "zhb2st_kernels.f"
	    i__1 = *lda - 1;
#line 217 "zhb2st_kernels.f"
	    zlarfy_(uplo, &lm, &v[vpos], &c__1, &z__1, &a[dpos + *st * a_dim1]
		    , &i__1, &work[1], (ftnlen)1);
#line 220 "zhb2st_kernels.f"
	}

#line 222 "zhb2st_kernels.f"
	if (*ttype == 2) {
#line 223 "zhb2st_kernels.f"
	    j1 = *ed + 1;
/* Computing MIN */
#line 224 "zhb2st_kernels.f"
	    i__1 = *ed + *nb;
#line 224 "zhb2st_kernels.f"
	    j2 = min(i__1,*n);
#line 225 "zhb2st_kernels.f"
	    ln = *ed - *st + 1;
#line 226 "zhb2st_kernels.f"
	    lm = j2 - j1 + 1;
#line 227 "zhb2st_kernels.f"
	    if (lm > 0) {
#line 228 "zhb2st_kernels.f"
		d_cnjg(&z__1, &tau[taupos]);
#line 228 "zhb2st_kernels.f"
		i__1 = *lda - 1;
#line 228 "zhb2st_kernels.f"
		zlarfx_("Left", &ln, &lm, &v[vpos], &z__1, &a[dpos - *nb + j1 
			* a_dim1], &i__1, &work[1], (ftnlen)4);

#line 232 "zhb2st_kernels.f"
		if (*wantz) {
#line 233 "zhb2st_kernels.f"
		    vpos = (*sweep - 1) % 2 * *n + j1;
#line 234 "zhb2st_kernels.f"
		    taupos = (*sweep - 1) % 2 * *n + j1;
#line 235 "zhb2st_kernels.f"
		} else {
#line 236 "zhb2st_kernels.f"
		    vpos = (*sweep - 1) % 2 * *n + j1;
#line 237 "zhb2st_kernels.f"
		    taupos = (*sweep - 1) % 2 * *n + j1;
#line 238 "zhb2st_kernels.f"
		}

#line 240 "zhb2st_kernels.f"
		i__1 = vpos;
#line 240 "zhb2st_kernels.f"
		v[i__1].r = 1., v[i__1].i = 0.;
#line 241 "zhb2st_kernels.f"
		i__1 = lm - 1;
#line 241 "zhb2st_kernels.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 242 "zhb2st_kernels.f"
		    i__2 = vpos + i__;
#line 242 "zhb2st_kernels.f"
		    d_cnjg(&z__1, &a[dpos - *nb - i__ + (j1 + i__) * a_dim1]);
#line 242 "zhb2st_kernels.f"
		    v[i__2].r = z__1.r, v[i__2].i = z__1.i;
#line 244 "zhb2st_kernels.f"
		    i__2 = dpos - *nb - i__ + (j1 + i__) * a_dim1;
#line 244 "zhb2st_kernels.f"
		    a[i__2].r = 0., a[i__2].i = 0.;
#line 245 "zhb2st_kernels.f"
/* L30: */
#line 245 "zhb2st_kernels.f"
		}
#line 246 "zhb2st_kernels.f"
		d_cnjg(&z__1, &a[dpos - *nb + j1 * a_dim1]);
#line 246 "zhb2st_kernels.f"
		ctmp.r = z__1.r, ctmp.i = z__1.i;
#line 247 "zhb2st_kernels.f"
		zlarfg_(&lm, &ctmp, &v[vpos + 1], &c__1, &tau[taupos]);
#line 248 "zhb2st_kernels.f"
		i__1 = dpos - *nb + j1 * a_dim1;
#line 248 "zhb2st_kernels.f"
		a[i__1].r = ctmp.r, a[i__1].i = ctmp.i;

#line 250 "zhb2st_kernels.f"
		i__1 = ln - 1;
#line 250 "zhb2st_kernels.f"
		i__2 = *lda - 1;
#line 250 "zhb2st_kernels.f"
		zlarfx_("Right", &i__1, &lm, &v[vpos], &tau[taupos], &a[dpos 
			- *nb + 1 + j1 * a_dim1], &i__2, &work[1], (ftnlen)5);
#line 253 "zhb2st_kernels.f"
	    }
#line 254 "zhb2st_kernels.f"
	}

/*     Lower case */

#line 258 "zhb2st_kernels.f"
    } else {

#line 260 "zhb2st_kernels.f"
	if (*wantz) {
#line 261 "zhb2st_kernels.f"
	    vpos = (*sweep - 1) % 2 * *n + *st;
#line 262 "zhb2st_kernels.f"
	    taupos = (*sweep - 1) % 2 * *n + *st;
#line 263 "zhb2st_kernels.f"
	} else {
#line 264 "zhb2st_kernels.f"
	    vpos = (*sweep - 1) % 2 * *n + *st;
#line 265 "zhb2st_kernels.f"
	    taupos = (*sweep - 1) % 2 * *n + *st;
#line 266 "zhb2st_kernels.f"
	}

#line 268 "zhb2st_kernels.f"
	if (*ttype == 1) {
#line 269 "zhb2st_kernels.f"
	    lm = *ed - *st + 1;

#line 271 "zhb2st_kernels.f"
	    i__1 = vpos;
#line 271 "zhb2st_kernels.f"
	    v[i__1].r = 1., v[i__1].i = 0.;
#line 272 "zhb2st_kernels.f"
	    i__1 = lm - 1;
#line 272 "zhb2st_kernels.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 273 "zhb2st_kernels.f"
		i__2 = vpos + i__;
#line 273 "zhb2st_kernels.f"
		i__3 = ofdpos + i__ + (*st - 1) * a_dim1;
#line 273 "zhb2st_kernels.f"
		v[i__2].r = a[i__3].r, v[i__2].i = a[i__3].i;
#line 274 "zhb2st_kernels.f"
		i__2 = ofdpos + i__ + (*st - 1) * a_dim1;
#line 274 "zhb2st_kernels.f"
		a[i__2].r = 0., a[i__2].i = 0.;
#line 275 "zhb2st_kernels.f"
/* L20: */
#line 275 "zhb2st_kernels.f"
	    }
#line 276 "zhb2st_kernels.f"
	    zlarfg_(&lm, &a[ofdpos + (*st - 1) * a_dim1], &v[vpos + 1], &c__1,
		     &tau[taupos]);

#line 279 "zhb2st_kernels.f"
	    lm = *ed - *st + 1;

#line 281 "zhb2st_kernels.f"
	    d_cnjg(&z__1, &tau[taupos]);
#line 281 "zhb2st_kernels.f"
	    i__1 = *lda - 1;
#line 281 "zhb2st_kernels.f"
	    zlarfy_(uplo, &lm, &v[vpos], &c__1, &z__1, &a[dpos + *st * a_dim1]
		    , &i__1, &work[1], (ftnlen)1);
#line 285 "zhb2st_kernels.f"
	}

#line 287 "zhb2st_kernels.f"
	if (*ttype == 3) {
#line 288 "zhb2st_kernels.f"
	    lm = *ed - *st + 1;

#line 290 "zhb2st_kernels.f"
	    d_cnjg(&z__1, &tau[taupos]);
#line 290 "zhb2st_kernels.f"
	    i__1 = *lda - 1;
#line 290 "zhb2st_kernels.f"
	    zlarfy_(uplo, &lm, &v[vpos], &c__1, &z__1, &a[dpos + *st * a_dim1]
		    , &i__1, &work[1], (ftnlen)1);
#line 294 "zhb2st_kernels.f"
	}

#line 296 "zhb2st_kernels.f"
	if (*ttype == 2) {
#line 297 "zhb2st_kernels.f"
	    j1 = *ed + 1;
/* Computing MIN */
#line 298 "zhb2st_kernels.f"
	    i__1 = *ed + *nb;
#line 298 "zhb2st_kernels.f"
	    j2 = min(i__1,*n);
#line 299 "zhb2st_kernels.f"
	    ln = *ed - *st + 1;
#line 300 "zhb2st_kernels.f"
	    lm = j2 - j1 + 1;

#line 302 "zhb2st_kernels.f"
	    if (lm > 0) {
#line 303 "zhb2st_kernels.f"
		i__1 = *lda - 1;
#line 303 "zhb2st_kernels.f"
		zlarfx_("Right", &lm, &ln, &v[vpos], &tau[taupos], &a[dpos + *
			nb + *st * a_dim1], &i__1, &work[1], (ftnlen)5);

#line 307 "zhb2st_kernels.f"
		if (*wantz) {
#line 308 "zhb2st_kernels.f"
		    vpos = (*sweep - 1) % 2 * *n + j1;
#line 309 "zhb2st_kernels.f"
		    taupos = (*sweep - 1) % 2 * *n + j1;
#line 310 "zhb2st_kernels.f"
		} else {
#line 311 "zhb2st_kernels.f"
		    vpos = (*sweep - 1) % 2 * *n + j1;
#line 312 "zhb2st_kernels.f"
		    taupos = (*sweep - 1) % 2 * *n + j1;
#line 313 "zhb2st_kernels.f"
		}

#line 315 "zhb2st_kernels.f"
		i__1 = vpos;
#line 315 "zhb2st_kernels.f"
		v[i__1].r = 1., v[i__1].i = 0.;
#line 316 "zhb2st_kernels.f"
		i__1 = lm - 1;
#line 316 "zhb2st_kernels.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 317 "zhb2st_kernels.f"
		    i__2 = vpos + i__;
#line 317 "zhb2st_kernels.f"
		    i__3 = dpos + *nb + i__ + *st * a_dim1;
#line 317 "zhb2st_kernels.f"
		    v[i__2].r = a[i__3].r, v[i__2].i = a[i__3].i;
#line 318 "zhb2st_kernels.f"
		    i__2 = dpos + *nb + i__ + *st * a_dim1;
#line 318 "zhb2st_kernels.f"
		    a[i__2].r = 0., a[i__2].i = 0.;
#line 319 "zhb2st_kernels.f"
/* L40: */
#line 319 "zhb2st_kernels.f"
		}
#line 320 "zhb2st_kernels.f"
		zlarfg_(&lm, &a[dpos + *nb + *st * a_dim1], &v[vpos + 1], &
			c__1, &tau[taupos]);

#line 323 "zhb2st_kernels.f"
		i__1 = ln - 1;
#line 323 "zhb2st_kernels.f"
		d_cnjg(&z__1, &tau[taupos]);
#line 323 "zhb2st_kernels.f"
		i__2 = *lda - 1;
#line 323 "zhb2st_kernels.f"
		zlarfx_("Left", &lm, &i__1, &v[vpos], &z__1, &a[dpos + *nb - 
			1 + (*st + 1) * a_dim1], &i__2, &work[1], (ftnlen)4);
#line 327 "zhb2st_kernels.f"
	    }
#line 328 "zhb2st_kernels.f"
	}
#line 329 "zhb2st_kernels.f"
    }

#line 331 "zhb2st_kernels.f"
    return 0;

/*     END OF ZHB2ST_KERNELS */

} /* zhb2st_kernels__ */

