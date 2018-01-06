#line 1 "cggbak.f"
/* cggbak.f -- translated by f2c (version 20100827).
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

#line 1 "cggbak.f"
/* > \brief \b CGGBAK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGGBAK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cggbak.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cggbak.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cggbak.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V, */
/*                          LDV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOB, SIDE */
/*       INTEGER            IHI, ILO, INFO, LDV, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               LSCALE( * ), RSCALE( * ) */
/*       COMPLEX            V( LDV, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGGBAK forms the right or left eigenvectors of a complex generalized */
/* > eigenvalue problem A*x = lambda*B*x, by backward transformation on */
/* > the computed eigenvectors of the balanced pair of matrices output by */
/* > CGGBAL. */
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
/* >          JOB must be the same as the argument JOB supplied to CGGBAL. */
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
/* >          The integers ILO and IHI determined by CGGBAL. */
/* >          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in] LSCALE */
/* > \verbatim */
/* >          LSCALE is REAL array, dimension (N) */
/* >          Details of the permutations and/or scaling factors applied */
/* >          to the left side of A and B, as returned by CGGBAL. */
/* > \endverbatim */
/* > */
/* > \param[in] RSCALE */
/* > \verbatim */
/* >          RSCALE is REAL array, dimension (N) */
/* >          Details of the permutations and/or scaling factors applied */
/* >          to the right side of A and B, as returned by CGGBAL. */
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
/* >          V is COMPLEX array, dimension (LDV,M) */
/* >          On entry, the matrix of right or left eigenvectors to be */
/* >          transformed, as returned by CTGEVC. */
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

/* > \ingroup complexGBcomputational */

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
/* Subroutine */ int cggbak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, doublereal *lscale, doublereal *rscale, integer *m, 
	doublecomplex *v, integer *ldv, integer *info, ftnlen job_len, ftnlen 
	side_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, i__1;

    /* Local variables */
    static integer i__, k;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical leftv;
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), xerbla_(char *, integer *, ftnlen);
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

#line 185 "cggbak.f"
    /* Parameter adjustments */
#line 185 "cggbak.f"
    --lscale;
#line 185 "cggbak.f"
    --rscale;
#line 185 "cggbak.f"
    v_dim1 = *ldv;
#line 185 "cggbak.f"
    v_offset = 1 + v_dim1;
#line 185 "cggbak.f"
    v -= v_offset;
#line 185 "cggbak.f"

#line 185 "cggbak.f"
    /* Function Body */
#line 185 "cggbak.f"
    rightv = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 186 "cggbak.f"
    leftv = lsame_(side, "L", (ftnlen)1, (ftnlen)1);

#line 188 "cggbak.f"
    *info = 0;
#line 189 "cggbak.f"
    if (! lsame_(job, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(job, "P", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(job, "S", (ftnlen)1, (ftnlen)1) 
	    && ! lsame_(job, "B", (ftnlen)1, (ftnlen)1)) {
#line 191 "cggbak.f"
	*info = -1;
#line 192 "cggbak.f"
    } else if (! rightv && ! leftv) {
#line 193 "cggbak.f"
	*info = -2;
#line 194 "cggbak.f"
    } else if (*n < 0) {
#line 195 "cggbak.f"
	*info = -3;
#line 196 "cggbak.f"
    } else if (*ilo < 1) {
#line 197 "cggbak.f"
	*info = -4;
#line 198 "cggbak.f"
    } else if (*n == 0 && *ihi == 0 && *ilo != 1) {
#line 199 "cggbak.f"
	*info = -4;
#line 200 "cggbak.f"
    } else if (*n > 0 && (*ihi < *ilo || *ihi > max(1,*n))) {
#line 202 "cggbak.f"
	*info = -5;
#line 203 "cggbak.f"
    } else if (*n == 0 && *ilo == 1 && *ihi != 0) {
#line 204 "cggbak.f"
	*info = -5;
#line 205 "cggbak.f"
    } else if (*m < 0) {
#line 206 "cggbak.f"
	*info = -8;
#line 207 "cggbak.f"
    } else if (*ldv < max(1,*n)) {
#line 208 "cggbak.f"
	*info = -10;
#line 209 "cggbak.f"
    }
#line 210 "cggbak.f"
    if (*info != 0) {
#line 211 "cggbak.f"
	i__1 = -(*info);
#line 211 "cggbak.f"
	xerbla_("CGGBAK", &i__1, (ftnlen)6);
#line 212 "cggbak.f"
	return 0;
#line 213 "cggbak.f"
    }

/*     Quick return if possible */

#line 217 "cggbak.f"
    if (*n == 0) {
#line 217 "cggbak.f"
	return 0;
#line 217 "cggbak.f"
    }
#line 219 "cggbak.f"
    if (*m == 0) {
#line 219 "cggbak.f"
	return 0;
#line 219 "cggbak.f"
    }
#line 221 "cggbak.f"
    if (lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 221 "cggbak.f"
	return 0;
#line 221 "cggbak.f"
    }

#line 224 "cggbak.f"
    if (*ilo == *ihi) {
#line 224 "cggbak.f"
	goto L30;
#line 224 "cggbak.f"
    }

/*     Backward balance */

#line 229 "cggbak.f"
    if (lsame_(job, "S", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (ftnlen)1, 
	    (ftnlen)1)) {

/*        Backward transformation on right eigenvectors */

#line 233 "cggbak.f"
	if (rightv) {
#line 234 "cggbak.f"
	    i__1 = *ihi;
#line 234 "cggbak.f"
	    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 235 "cggbak.f"
		csscal_(m, &rscale[i__], &v[i__ + v_dim1], ldv);
#line 236 "cggbak.f"
/* L10: */
#line 236 "cggbak.f"
	    }
#line 237 "cggbak.f"
	}

/*        Backward transformation on left eigenvectors */

#line 241 "cggbak.f"
	if (leftv) {
#line 242 "cggbak.f"
	    i__1 = *ihi;
#line 242 "cggbak.f"
	    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 243 "cggbak.f"
		csscal_(m, &lscale[i__], &v[i__ + v_dim1], ldv);
#line 244 "cggbak.f"
/* L20: */
#line 244 "cggbak.f"
	    }
#line 245 "cggbak.f"
	}
#line 246 "cggbak.f"
    }

/*     Backward permutation */

#line 250 "cggbak.f"
L30:
#line 251 "cggbak.f"
    if (lsame_(job, "P", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (ftnlen)1, 
	    (ftnlen)1)) {

/*        Backward permutation on right eigenvectors */

#line 255 "cggbak.f"
	if (rightv) {
#line 256 "cggbak.f"
	    if (*ilo == 1) {
#line 256 "cggbak.f"
		goto L50;
#line 256 "cggbak.f"
	    }
#line 258 "cggbak.f"
	    for (i__ = *ilo - 1; i__ >= 1; --i__) {
#line 259 "cggbak.f"
		k = (integer) rscale[i__];
#line 260 "cggbak.f"
		if (k == i__) {
#line 260 "cggbak.f"
		    goto L40;
#line 260 "cggbak.f"
		}
#line 262 "cggbak.f"
		cswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 263 "cggbak.f"
L40:
#line 263 "cggbak.f"
		;
#line 263 "cggbak.f"
	    }

#line 265 "cggbak.f"
L50:
#line 266 "cggbak.f"
	    if (*ihi == *n) {
#line 266 "cggbak.f"
		goto L70;
#line 266 "cggbak.f"
	    }
#line 268 "cggbak.f"
	    i__1 = *n;
#line 268 "cggbak.f"
	    for (i__ = *ihi + 1; i__ <= i__1; ++i__) {
#line 269 "cggbak.f"
		k = (integer) rscale[i__];
#line 270 "cggbak.f"
		if (k == i__) {
#line 270 "cggbak.f"
		    goto L60;
#line 270 "cggbak.f"
		}
#line 272 "cggbak.f"
		cswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 273 "cggbak.f"
L60:
#line 273 "cggbak.f"
		;
#line 273 "cggbak.f"
	    }
#line 274 "cggbak.f"
	}

/*        Backward permutation on left eigenvectors */

#line 278 "cggbak.f"
L70:
#line 279 "cggbak.f"
	if (leftv) {
#line 280 "cggbak.f"
	    if (*ilo == 1) {
#line 280 "cggbak.f"
		goto L90;
#line 280 "cggbak.f"
	    }
#line 282 "cggbak.f"
	    for (i__ = *ilo - 1; i__ >= 1; --i__) {
#line 283 "cggbak.f"
		k = (integer) lscale[i__];
#line 284 "cggbak.f"
		if (k == i__) {
#line 284 "cggbak.f"
		    goto L80;
#line 284 "cggbak.f"
		}
#line 286 "cggbak.f"
		cswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 287 "cggbak.f"
L80:
#line 287 "cggbak.f"
		;
#line 287 "cggbak.f"
	    }

#line 289 "cggbak.f"
L90:
#line 290 "cggbak.f"
	    if (*ihi == *n) {
#line 290 "cggbak.f"
		goto L110;
#line 290 "cggbak.f"
	    }
#line 292 "cggbak.f"
	    i__1 = *n;
#line 292 "cggbak.f"
	    for (i__ = *ihi + 1; i__ <= i__1; ++i__) {
#line 293 "cggbak.f"
		k = (integer) lscale[i__];
#line 294 "cggbak.f"
		if (k == i__) {
#line 294 "cggbak.f"
		    goto L100;
#line 294 "cggbak.f"
		}
#line 296 "cggbak.f"
		cswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 297 "cggbak.f"
L100:
#line 297 "cggbak.f"
		;
#line 297 "cggbak.f"
	    }
#line 298 "cggbak.f"
	}
#line 299 "cggbak.f"
    }

#line 301 "cggbak.f"
L110:

#line 303 "cggbak.f"
    return 0;

/*     End of CGGBAK */

} /* cggbak_ */

