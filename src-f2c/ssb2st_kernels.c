#line 1 "ssb2st_kernels.f"
/* ssb2st_kernels.f -- translated by f2c (version 20100827).
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

#line 1 "ssb2st_kernels.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SSB2ST_KERNELS */

/*  @generated from zhb2st_kernels.f, fortran z -> s, Wed Dec  7 08:22:40 2016 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSB2ST_KERNELS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssb2st_
kernels.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssb2st_
kernels.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssb2st_
kernels.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE  SSB2ST_KERNELS( UPLO, WANTZ, TTYPE, */
/*                                   ST, ED, SWEEP, N, NB, IB, */
/*                                   A, LDA, V, TAU, LDVT, WORK) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       LOGICAL            WANTZ */
/*       INTEGER            TTYPE, ST, ED, SWEEP, N, NB, IB, LDA, LDVT */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), V( * ), */
/*                          TAU( * ), WORK( * ) */

/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSB2ST_KERNELS is an internal routine used by the SSYTRD_SB2ST */
/* > subroutine. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* >          WANTZ is LOGICAL which indicate if Eigenvalue are requested or both */
/* >          Eigenvalue/Eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] TTYPE */
/* > \verbatim */
/* >          TTYPE is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] ST */
/* > \verbatim */
/* >          ST is INTEGER */
/* >          internal parameter for indices. */
/* > \endverbatim */
/* > */
/* > \param[in] ED */
/* > \verbatim */
/* >          ED is INTEGER */
/* >          internal parameter for indices. */
/* > \endverbatim */
/* > */
/* > \param[in] SWEEP */
/* > \verbatim */
/* >          SWEEP is INTEGER */
/* >          internal parameter for indices. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER. The order of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER. The size of the band. */
/* > \endverbatim */
/* > */
/* > \param[in] IB */
/* > \verbatim */
/* >          IB is INTEGER. */
/* > \endverbatim */
/* > */
/* > \param[in, out] A */
/* > \verbatim */
/* >          A is REAL array. A pointer to the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER. The leading dimension of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* >          V is REAL array, dimension 2*n if eigenvalues only are */
/* >          requested or to be queried for vectors. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is REAL array, dimension (2*n). */
/* >          The scalar factors of the Householder reflectors are stored */
/* >          in this array. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT */
/* > \verbatim */
/* >          LDVT is INTEGER. */
/* > \endverbatim */
/* > */
/* > \param[in] WORK */
/* > \verbatim */
/* >          WORK is REAL array. Workspace of size nb. */
/* > \endverbatim */
/* > @param[in] n */
/* >          The order of the matrix A. */
/* > */
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
/* Subroutine */ int ssb2st_kernels__(char *uplo, logical *wantz, integer *
	ttype, integer *st, integer *ed, integer *sweep, integer *n, integer *
	nb, integer *ib, doublereal *a, integer *lda, doublereal *v, 
	doublereal *tau, integer *ldvt, doublereal *work, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j1, j2, lm, ln;
    static doublereal ctmp;
    static integer dpos, vpos;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    static integer ajeter;
    extern /* Subroutine */ int slarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    static integer ofdpos;
    extern /* Subroutine */ int slarfx_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     ftnlen), slarfy_(char *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, ftnlen);
    static integer taupos;



/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

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

#line 216 "ssb2st_kernels.f"
    /* Parameter adjustments */
#line 216 "ssb2st_kernels.f"
    a_dim1 = *lda;
#line 216 "ssb2st_kernels.f"
    a_offset = 1 + a_dim1;
#line 216 "ssb2st_kernels.f"
    a -= a_offset;
#line 216 "ssb2st_kernels.f"
    --v;
#line 216 "ssb2st_kernels.f"
    --tau;
#line 216 "ssb2st_kernels.f"
    --work;
#line 216 "ssb2st_kernels.f"

#line 216 "ssb2st_kernels.f"
    /* Function Body */
#line 216 "ssb2st_kernels.f"
    ajeter = *ib + *ldvt;
#line 217 "ssb2st_kernels.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 219 "ssb2st_kernels.f"
    if (upper) {
#line 220 "ssb2st_kernels.f"
	dpos = (*nb << 1) + 1;
#line 221 "ssb2st_kernels.f"
	ofdpos = *nb << 1;
#line 222 "ssb2st_kernels.f"
    } else {
#line 223 "ssb2st_kernels.f"
	dpos = 1;
#line 224 "ssb2st_kernels.f"
	ofdpos = 2;
#line 225 "ssb2st_kernels.f"
    }

/*     Upper case */

#line 230 "ssb2st_kernels.f"
    if (upper) {

#line 232 "ssb2st_kernels.f"
	if (*wantz) {
#line 233 "ssb2st_kernels.f"
	    vpos = (*sweep - 1) % 2 * *n + *st;
#line 234 "ssb2st_kernels.f"
	    taupos = (*sweep - 1) % 2 * *n + *st;
#line 235 "ssb2st_kernels.f"
	} else {
#line 236 "ssb2st_kernels.f"
	    vpos = (*sweep - 1) % 2 * *n + *st;
#line 237 "ssb2st_kernels.f"
	    taupos = (*sweep - 1) % 2 * *n + *st;
#line 238 "ssb2st_kernels.f"
	}

#line 240 "ssb2st_kernels.f"
	if (*ttype == 1) {
#line 241 "ssb2st_kernels.f"
	    lm = *ed - *st + 1;

#line 243 "ssb2st_kernels.f"
	    v[vpos] = 1.;
#line 244 "ssb2st_kernels.f"
	    i__1 = lm - 1;
#line 244 "ssb2st_kernels.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 245 "ssb2st_kernels.f"
		v[vpos + i__] = a[ofdpos - i__ + (*st + i__) * a_dim1];
#line 246 "ssb2st_kernels.f"
		a[ofdpos - i__ + (*st + i__) * a_dim1] = 0.;
#line 247 "ssb2st_kernels.f"
/* L10: */
#line 247 "ssb2st_kernels.f"
	    }
#line 248 "ssb2st_kernels.f"
	    ctmp = a[ofdpos + *st * a_dim1];
#line 249 "ssb2st_kernels.f"
	    slarfg_(&lm, &ctmp, &v[vpos + 1], &c__1, &tau[taupos]);
#line 251 "ssb2st_kernels.f"
	    a[ofdpos + *st * a_dim1] = ctmp;

#line 253 "ssb2st_kernels.f"
	    lm = *ed - *st + 1;
#line 254 "ssb2st_kernels.f"
	    d__1 = tau[taupos];
#line 254 "ssb2st_kernels.f"
	    i__1 = *lda - 1;
#line 254 "ssb2st_kernels.f"
	    slarfy_(uplo, &lm, &v[vpos], &c__1, &d__1, &a[dpos + *st * a_dim1]
		    , &i__1, &work[1], (ftnlen)1);
#line 257 "ssb2st_kernels.f"
	}

#line 259 "ssb2st_kernels.f"
	if (*ttype == 3) {

#line 261 "ssb2st_kernels.f"
	    lm = *ed - *st + 1;
#line 262 "ssb2st_kernels.f"
	    d__1 = tau[taupos];
#line 262 "ssb2st_kernels.f"
	    i__1 = *lda - 1;
#line 262 "ssb2st_kernels.f"
	    slarfy_(uplo, &lm, &v[vpos], &c__1, &d__1, &a[dpos + *st * a_dim1]
		    , &i__1, &work[1], (ftnlen)1);
#line 265 "ssb2st_kernels.f"
	}

#line 267 "ssb2st_kernels.f"
	if (*ttype == 2) {
#line 268 "ssb2st_kernels.f"
	    j1 = *ed + 1;
/* Computing MIN */
#line 269 "ssb2st_kernels.f"
	    i__1 = *ed + *nb;
#line 269 "ssb2st_kernels.f"
	    j2 = min(i__1,*n);
#line 270 "ssb2st_kernels.f"
	    ln = *ed - *st + 1;
#line 271 "ssb2st_kernels.f"
	    lm = j2 - j1 + 1;
#line 272 "ssb2st_kernels.f"
	    if (lm > 0) {
#line 273 "ssb2st_kernels.f"
		d__1 = tau[taupos];
#line 273 "ssb2st_kernels.f"
		i__1 = *lda - 1;
#line 273 "ssb2st_kernels.f"
		slarfx_("Left", &ln, &lm, &v[vpos], &d__1, &a[dpos - *nb + j1 
			* a_dim1], &i__1, &work[1], (ftnlen)4);

#line 277 "ssb2st_kernels.f"
		if (*wantz) {
#line 278 "ssb2st_kernels.f"
		    vpos = (*sweep - 1) % 2 * *n + j1;
#line 279 "ssb2st_kernels.f"
		    taupos = (*sweep - 1) % 2 * *n + j1;
#line 280 "ssb2st_kernels.f"
		} else {
#line 281 "ssb2st_kernels.f"
		    vpos = (*sweep - 1) % 2 * *n + j1;
#line 282 "ssb2st_kernels.f"
		    taupos = (*sweep - 1) % 2 * *n + j1;
#line 283 "ssb2st_kernels.f"
		}

#line 285 "ssb2st_kernels.f"
		v[vpos] = 1.;
#line 286 "ssb2st_kernels.f"
		i__1 = lm - 1;
#line 286 "ssb2st_kernels.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 287 "ssb2st_kernels.f"
		    v[vpos + i__] = a[dpos - *nb - i__ + (j1 + i__) * a_dim1];
#line 289 "ssb2st_kernels.f"
		    a[dpos - *nb - i__ + (j1 + i__) * a_dim1] = 0.;
#line 290 "ssb2st_kernels.f"
/* L30: */
#line 290 "ssb2st_kernels.f"
		}
#line 291 "ssb2st_kernels.f"
		ctmp = a[dpos - *nb + j1 * a_dim1];
#line 292 "ssb2st_kernels.f"
		slarfg_(&lm, &ctmp, &v[vpos + 1], &c__1, &tau[taupos]);
#line 293 "ssb2st_kernels.f"
		a[dpos - *nb + j1 * a_dim1] = ctmp;

#line 295 "ssb2st_kernels.f"
		i__1 = ln - 1;
#line 295 "ssb2st_kernels.f"
		i__2 = *lda - 1;
#line 295 "ssb2st_kernels.f"
		slarfx_("Right", &i__1, &lm, &v[vpos], &tau[taupos], &a[dpos 
			- *nb + 1 + j1 * a_dim1], &i__2, &work[1], (ftnlen)5);
#line 298 "ssb2st_kernels.f"
	    }
#line 299 "ssb2st_kernels.f"
	}

/*     Lower case */

#line 303 "ssb2st_kernels.f"
    } else {

#line 305 "ssb2st_kernels.f"
	if (*wantz) {
#line 306 "ssb2st_kernels.f"
	    vpos = (*sweep - 1) % 2 * *n + *st;
#line 307 "ssb2st_kernels.f"
	    taupos = (*sweep - 1) % 2 * *n + *st;
#line 308 "ssb2st_kernels.f"
	} else {
#line 309 "ssb2st_kernels.f"
	    vpos = (*sweep - 1) % 2 * *n + *st;
#line 310 "ssb2st_kernels.f"
	    taupos = (*sweep - 1) % 2 * *n + *st;
#line 311 "ssb2st_kernels.f"
	}

#line 313 "ssb2st_kernels.f"
	if (*ttype == 1) {
#line 314 "ssb2st_kernels.f"
	    lm = *ed - *st + 1;

#line 316 "ssb2st_kernels.f"
	    v[vpos] = 1.;
#line 317 "ssb2st_kernels.f"
	    i__1 = lm - 1;
#line 317 "ssb2st_kernels.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 318 "ssb2st_kernels.f"
		v[vpos + i__] = a[ofdpos + i__ + (*st - 1) * a_dim1];
#line 319 "ssb2st_kernels.f"
		a[ofdpos + i__ + (*st - 1) * a_dim1] = 0.;
#line 320 "ssb2st_kernels.f"
/* L20: */
#line 320 "ssb2st_kernels.f"
	    }
#line 321 "ssb2st_kernels.f"
	    slarfg_(&lm, &a[ofdpos + (*st - 1) * a_dim1], &v[vpos + 1], &c__1,
		     &tau[taupos]);

#line 324 "ssb2st_kernels.f"
	    lm = *ed - *st + 1;

#line 326 "ssb2st_kernels.f"
	    d__1 = tau[taupos];
#line 326 "ssb2st_kernels.f"
	    i__1 = *lda - 1;
#line 326 "ssb2st_kernels.f"
	    slarfy_(uplo, &lm, &v[vpos], &c__1, &d__1, &a[dpos + *st * a_dim1]
		    , &i__1, &work[1], (ftnlen)1);
#line 330 "ssb2st_kernels.f"
	}

#line 332 "ssb2st_kernels.f"
	if (*ttype == 3) {
#line 333 "ssb2st_kernels.f"
	    lm = *ed - *st + 1;

#line 335 "ssb2st_kernels.f"
	    d__1 = tau[taupos];
#line 335 "ssb2st_kernels.f"
	    i__1 = *lda - 1;
#line 335 "ssb2st_kernels.f"
	    slarfy_(uplo, &lm, &v[vpos], &c__1, &d__1, &a[dpos + *st * a_dim1]
		    , &i__1, &work[1], (ftnlen)1);
#line 339 "ssb2st_kernels.f"
	}

#line 341 "ssb2st_kernels.f"
	if (*ttype == 2) {
#line 342 "ssb2st_kernels.f"
	    j1 = *ed + 1;
/* Computing MIN */
#line 343 "ssb2st_kernels.f"
	    i__1 = *ed + *nb;
#line 343 "ssb2st_kernels.f"
	    j2 = min(i__1,*n);
#line 344 "ssb2st_kernels.f"
	    ln = *ed - *st + 1;
#line 345 "ssb2st_kernels.f"
	    lm = j2 - j1 + 1;

#line 347 "ssb2st_kernels.f"
	    if (lm > 0) {
#line 348 "ssb2st_kernels.f"
		i__1 = *lda - 1;
#line 348 "ssb2st_kernels.f"
		slarfx_("Right", &lm, &ln, &v[vpos], &tau[taupos], &a[dpos + *
			nb + *st * a_dim1], &i__1, &work[1], (ftnlen)5);

#line 352 "ssb2st_kernels.f"
		if (*wantz) {
#line 353 "ssb2st_kernels.f"
		    vpos = (*sweep - 1) % 2 * *n + j1;
#line 354 "ssb2st_kernels.f"
		    taupos = (*sweep - 1) % 2 * *n + j1;
#line 355 "ssb2st_kernels.f"
		} else {
#line 356 "ssb2st_kernels.f"
		    vpos = (*sweep - 1) % 2 * *n + j1;
#line 357 "ssb2st_kernels.f"
		    taupos = (*sweep - 1) % 2 * *n + j1;
#line 358 "ssb2st_kernels.f"
		}

#line 360 "ssb2st_kernels.f"
		v[vpos] = 1.;
#line 361 "ssb2st_kernels.f"
		i__1 = lm - 1;
#line 361 "ssb2st_kernels.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 362 "ssb2st_kernels.f"
		    v[vpos + i__] = a[dpos + *nb + i__ + *st * a_dim1];
#line 363 "ssb2st_kernels.f"
		    a[dpos + *nb + i__ + *st * a_dim1] = 0.;
#line 364 "ssb2st_kernels.f"
/* L40: */
#line 364 "ssb2st_kernels.f"
		}
#line 365 "ssb2st_kernels.f"
		slarfg_(&lm, &a[dpos + *nb + *st * a_dim1], &v[vpos + 1], &
			c__1, &tau[taupos]);

#line 368 "ssb2st_kernels.f"
		i__1 = ln - 1;
#line 368 "ssb2st_kernels.f"
		d__1 = tau[taupos];
#line 368 "ssb2st_kernels.f"
		i__2 = *lda - 1;
#line 368 "ssb2st_kernels.f"
		slarfx_("Left", &lm, &i__1, &v[vpos], &d__1, &a[dpos + *nb - 
			1 + (*st + 1) * a_dim1], &i__2, &work[1], (ftnlen)4);
#line 372 "ssb2st_kernels.f"
	    }
#line 373 "ssb2st_kernels.f"
	}
#line 374 "ssb2st_kernels.f"
    }

#line 376 "ssb2st_kernels.f"
    return 0;

/*     END OF SSB2ST_KERNELS */

} /* ssb2st_kernels__ */

