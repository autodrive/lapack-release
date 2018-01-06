#line 1 "zggbak.f"
/* zggbak.f -- translated by f2c (version 20100827).
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

#line 1 "zggbak.f"
/* > \brief \b ZGGBAK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGGBAK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zggbak.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zggbak.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zggbak.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V, */
/*                          LDV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOB, SIDE */
/*       INTEGER            IHI, ILO, INFO, LDV, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   LSCALE( * ), RSCALE( * ) */
/*       COMPLEX*16         V( LDV, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGGBAK forms the right or left eigenvectors of a complex generalized */
/* > eigenvalue problem A*x = lambda*B*x, by backward transformation on */
/* > the computed eigenvectors of the balanced pair of matrices output by */
/* > ZGGBAL. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOB */
/* > \verbatim */
/* >          JOB is CHARACTER*1 */
/* >          Specifies the type of backward transformation required: */
/* >          = 'N':  do nothing, return immediately; */
/* >          = 'P':  do backward transformation for permutation only; */
/* >          = 'S':  do backward transformation for scaling only; */
/* >          = 'B':  do backward transformations for both permutation and */
/* >                  scaling. */
/* >          JOB must be the same as the argument JOB supplied to ZGGBAL. */
/* > \endverbatim */
/* > */
/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'R':  V contains right eigenvectors; */
/* >          = 'L':  V contains left eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of rows of the matrix V.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* >          ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* >          IHI is INTEGER */
/* >          The integers ILO and IHI determined by ZGGBAL. */
/* >          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in] LSCALE */
/* > \verbatim */
/* >          LSCALE is DOUBLE PRECISION array, dimension (N) */
/* >          Details of the permutations and/or scaling factors applied */
/* >          to the left side of A and B, as returned by ZGGBAL. */
/* > \endverbatim */
/* > */
/* > \param[in] RSCALE */
/* > \verbatim */
/* >          RSCALE is DOUBLE PRECISION array, dimension (N) */
/* >          Details of the permutations and/or scaling factors applied */
/* >          to the right side of A and B, as returned by ZGGBAL. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of columns of the matrix V.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] V */
/* > \verbatim */
/* >          V is COMPLEX*16 array, dimension (LDV,M) */
/* >          On entry, the matrix of right or left eigenvectors to be */
/* >          transformed, as returned by ZTGEVC. */
/* >          On exit, V is overwritten by the transformed eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the matrix V. LDV >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16GBcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  See R.C. Ward, Balancing the generalized eigenvalue problem, */
/* >                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zggbak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, doublereal *lscale, doublereal *rscale, integer *m, 
	doublecomplex *v, integer *ldv, integer *info, ftnlen job_len, ftnlen 
	side_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, i__1;

    /* Local variables */
    static integer i__, k;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical leftv;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), xerbla_(char *, integer *, ftnlen), 
	    zdscal_(integer *, doublereal *, doublecomplex *, integer *);
    static logical rightv;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 185 "zggbak.f"
    /* Parameter adjustments */
#line 185 "zggbak.f"
    --lscale;
#line 185 "zggbak.f"
    --rscale;
#line 185 "zggbak.f"
    v_dim1 = *ldv;
#line 185 "zggbak.f"
    v_offset = 1 + v_dim1;
#line 185 "zggbak.f"
    v -= v_offset;
#line 185 "zggbak.f"

#line 185 "zggbak.f"
    /* Function Body */
#line 185 "zggbak.f"
    rightv = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 186 "zggbak.f"
    leftv = lsame_(side, "L", (ftnlen)1, (ftnlen)1);

#line 188 "zggbak.f"
    *info = 0;
#line 189 "zggbak.f"
    if (! lsame_(job, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(job, "P", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(job, "S", (ftnlen)1, (ftnlen)1) 
	    && ! lsame_(job, "B", (ftnlen)1, (ftnlen)1)) {
#line 191 "zggbak.f"
	*info = -1;
#line 192 "zggbak.f"
    } else if (! rightv && ! leftv) {
#line 193 "zggbak.f"
	*info = -2;
#line 194 "zggbak.f"
    } else if (*n < 0) {
#line 195 "zggbak.f"
	*info = -3;
#line 196 "zggbak.f"
    } else if (*ilo < 1) {
#line 197 "zggbak.f"
	*info = -4;
#line 198 "zggbak.f"
    } else if (*n == 0 && *ihi == 0 && *ilo != 1) {
#line 199 "zggbak.f"
	*info = -4;
#line 200 "zggbak.f"
    } else if (*n > 0 && (*ihi < *ilo || *ihi > max(1,*n))) {
#line 202 "zggbak.f"
	*info = -5;
#line 203 "zggbak.f"
    } else if (*n == 0 && *ilo == 1 && *ihi != 0) {
#line 204 "zggbak.f"
	*info = -5;
#line 205 "zggbak.f"
    } else if (*m < 0) {
#line 206 "zggbak.f"
	*info = -8;
#line 207 "zggbak.f"
    } else if (*ldv < max(1,*n)) {
#line 208 "zggbak.f"
	*info = -10;
#line 209 "zggbak.f"
    }
#line 210 "zggbak.f"
    if (*info != 0) {
#line 211 "zggbak.f"
	i__1 = -(*info);
#line 211 "zggbak.f"
	xerbla_("ZGGBAK", &i__1, (ftnlen)6);
#line 212 "zggbak.f"
	return 0;
#line 213 "zggbak.f"
    }

/*     Quick return if possible */

#line 217 "zggbak.f"
    if (*n == 0) {
#line 217 "zggbak.f"
	return 0;
#line 217 "zggbak.f"
    }
#line 219 "zggbak.f"
    if (*m == 0) {
#line 219 "zggbak.f"
	return 0;
#line 219 "zggbak.f"
    }
#line 221 "zggbak.f"
    if (lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 221 "zggbak.f"
	return 0;
#line 221 "zggbak.f"
    }

#line 224 "zggbak.f"
    if (*ilo == *ihi) {
#line 224 "zggbak.f"
	goto L30;
#line 224 "zggbak.f"
    }

/*     Backward balance */

#line 229 "zggbak.f"
    if (lsame_(job, "S", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (ftnlen)1, 
	    (ftnlen)1)) {

/*        Backward transformation on right eigenvectors */

#line 233 "zggbak.f"
	if (rightv) {
#line 234 "zggbak.f"
	    i__1 = *ihi;
#line 234 "zggbak.f"
	    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 235 "zggbak.f"
		zdscal_(m, &rscale[i__], &v[i__ + v_dim1], ldv);
#line 236 "zggbak.f"
/* L10: */
#line 236 "zggbak.f"
	    }
#line 237 "zggbak.f"
	}

/*        Backward transformation on left eigenvectors */

#line 241 "zggbak.f"
	if (leftv) {
#line 242 "zggbak.f"
	    i__1 = *ihi;
#line 242 "zggbak.f"
	    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 243 "zggbak.f"
		zdscal_(m, &lscale[i__], &v[i__ + v_dim1], ldv);
#line 244 "zggbak.f"
/* L20: */
#line 244 "zggbak.f"
	    }
#line 245 "zggbak.f"
	}
#line 246 "zggbak.f"
    }

/*     Backward permutation */

#line 250 "zggbak.f"
L30:
#line 251 "zggbak.f"
    if (lsame_(job, "P", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (ftnlen)1, 
	    (ftnlen)1)) {

/*        Backward permutation on right eigenvectors */

#line 255 "zggbak.f"
	if (rightv) {
#line 256 "zggbak.f"
	    if (*ilo == 1) {
#line 256 "zggbak.f"
		goto L50;
#line 256 "zggbak.f"
	    }
#line 258 "zggbak.f"
	    for (i__ = *ilo - 1; i__ >= 1; --i__) {
#line 259 "zggbak.f"
		k = (integer) rscale[i__];
#line 260 "zggbak.f"
		if (k == i__) {
#line 260 "zggbak.f"
		    goto L40;
#line 260 "zggbak.f"
		}
#line 262 "zggbak.f"
		zswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 263 "zggbak.f"
L40:
#line 263 "zggbak.f"
		;
#line 263 "zggbak.f"
	    }

#line 265 "zggbak.f"
L50:
#line 266 "zggbak.f"
	    if (*ihi == *n) {
#line 266 "zggbak.f"
		goto L70;
#line 266 "zggbak.f"
	    }
#line 268 "zggbak.f"
	    i__1 = *n;
#line 268 "zggbak.f"
	    for (i__ = *ihi + 1; i__ <= i__1; ++i__) {
#line 269 "zggbak.f"
		k = (integer) rscale[i__];
#line 270 "zggbak.f"
		if (k == i__) {
#line 270 "zggbak.f"
		    goto L60;
#line 270 "zggbak.f"
		}
#line 272 "zggbak.f"
		zswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 273 "zggbak.f"
L60:
#line 273 "zggbak.f"
		;
#line 273 "zggbak.f"
	    }
#line 274 "zggbak.f"
	}

/*        Backward permutation on left eigenvectors */

#line 278 "zggbak.f"
L70:
#line 279 "zggbak.f"
	if (leftv) {
#line 280 "zggbak.f"
	    if (*ilo == 1) {
#line 280 "zggbak.f"
		goto L90;
#line 280 "zggbak.f"
	    }
#line 282 "zggbak.f"
	    for (i__ = *ilo - 1; i__ >= 1; --i__) {
#line 283 "zggbak.f"
		k = (integer) lscale[i__];
#line 284 "zggbak.f"
		if (k == i__) {
#line 284 "zggbak.f"
		    goto L80;
#line 284 "zggbak.f"
		}
#line 286 "zggbak.f"
		zswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 287 "zggbak.f"
L80:
#line 287 "zggbak.f"
		;
#line 287 "zggbak.f"
	    }

#line 289 "zggbak.f"
L90:
#line 290 "zggbak.f"
	    if (*ihi == *n) {
#line 290 "zggbak.f"
		goto L110;
#line 290 "zggbak.f"
	    }
#line 292 "zggbak.f"
	    i__1 = *n;
#line 292 "zggbak.f"
	    for (i__ = *ihi + 1; i__ <= i__1; ++i__) {
#line 293 "zggbak.f"
		k = (integer) lscale[i__];
#line 294 "zggbak.f"
		if (k == i__) {
#line 294 "zggbak.f"
		    goto L100;
#line 294 "zggbak.f"
		}
#line 296 "zggbak.f"
		zswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 297 "zggbak.f"
L100:
#line 297 "zggbak.f"
		;
#line 297 "zggbak.f"
	    }
#line 298 "zggbak.f"
	}
#line 299 "zggbak.f"
    }

#line 301 "zggbak.f"
L110:

#line 303 "zggbak.f"
    return 0;

/*     End of ZGGBAK */

} /* zggbak_ */

