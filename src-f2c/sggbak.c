#line 1 "sggbak.f"
/* sggbak.f -- translated by f2c (version 20100827).
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

#line 1 "sggbak.f"
/* > \brief \b SGGBAK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGGBAK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sggbak.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sggbak.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sggbak.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V, */
/*                          LDV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOB, SIDE */
/*       INTEGER            IHI, ILO, INFO, LDV, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               LSCALE( * ), RSCALE( * ), V( LDV, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGGBAK forms the right or left eigenvectors of a real generalized */
/* > eigenvalue problem A*x = lambda*B*x, by backward transformation on */
/* > the computed eigenvectors of the balanced pair of matrices output by */
/* > SGGBAL. */
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
/* >          JOB must be the same as the argument JOB supplied to SGGBAL. */
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
/* >          The integers ILO and IHI determined by SGGBAL. */
/* >          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in] LSCALE */
/* > \verbatim */
/* >          LSCALE is REAL array, dimension (N) */
/* >          Details of the permutations and/or scaling factors applied */
/* >          to the left side of A and B, as returned by SGGBAL. */
/* > \endverbatim */
/* > */
/* > \param[in] RSCALE */
/* > \verbatim */
/* >          RSCALE is REAL array, dimension (N) */
/* >          Details of the permutations and/or scaling factors applied */
/* >          to the right side of A and B, as returned by SGGBAL. */
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
/* >          V is REAL array, dimension (LDV,M) */
/* >          On entry, the matrix of right or left eigenvectors to be */
/* >          transformed, as returned by STGEVC. */
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

/* > \date November 2011 */

/* > \ingroup realGBcomputational */

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
/* Subroutine */ int sggbak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, doublereal *lscale, doublereal *rscale, integer *m, 
	doublereal *v, integer *ldv, integer *info, ftnlen job_len, ftnlen 
	side_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, i__1;

    /* Local variables */
    static integer i__, k;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical leftv;
    extern /* Subroutine */ int sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), xerbla_(char *, integer *, ftnlen);
    static logical rightv;


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 183 "sggbak.f"
    /* Parameter adjustments */
#line 183 "sggbak.f"
    --lscale;
#line 183 "sggbak.f"
    --rscale;
#line 183 "sggbak.f"
    v_dim1 = *ldv;
#line 183 "sggbak.f"
    v_offset = 1 + v_dim1;
#line 183 "sggbak.f"
    v -= v_offset;
#line 183 "sggbak.f"

#line 183 "sggbak.f"
    /* Function Body */
#line 183 "sggbak.f"
    rightv = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 184 "sggbak.f"
    leftv = lsame_(side, "L", (ftnlen)1, (ftnlen)1);

#line 186 "sggbak.f"
    *info = 0;
#line 187 "sggbak.f"
    if (! lsame_(job, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(job, "P", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(job, "S", (ftnlen)1, (ftnlen)1) 
	    && ! lsame_(job, "B", (ftnlen)1, (ftnlen)1)) {
#line 189 "sggbak.f"
	*info = -1;
#line 190 "sggbak.f"
    } else if (! rightv && ! leftv) {
#line 191 "sggbak.f"
	*info = -2;
#line 192 "sggbak.f"
    } else if (*n < 0) {
#line 193 "sggbak.f"
	*info = -3;
#line 194 "sggbak.f"
    } else if (*ilo < 1) {
#line 195 "sggbak.f"
	*info = -4;
#line 196 "sggbak.f"
    } else if (*n == 0 && *ihi == 0 && *ilo != 1) {
#line 197 "sggbak.f"
	*info = -4;
#line 198 "sggbak.f"
    } else if (*n > 0 && (*ihi < *ilo || *ihi > max(1,*n))) {
#line 200 "sggbak.f"
	*info = -5;
#line 201 "sggbak.f"
    } else if (*n == 0 && *ilo == 1 && *ihi != 0) {
#line 202 "sggbak.f"
	*info = -5;
#line 203 "sggbak.f"
    } else if (*m < 0) {
#line 204 "sggbak.f"
	*info = -8;
#line 205 "sggbak.f"
    } else if (*ldv < max(1,*n)) {
#line 206 "sggbak.f"
	*info = -10;
#line 207 "sggbak.f"
    }
#line 208 "sggbak.f"
    if (*info != 0) {
#line 209 "sggbak.f"
	i__1 = -(*info);
#line 209 "sggbak.f"
	xerbla_("SGGBAK", &i__1, (ftnlen)6);
#line 210 "sggbak.f"
	return 0;
#line 211 "sggbak.f"
    }

/*     Quick return if possible */

#line 215 "sggbak.f"
    if (*n == 0) {
#line 215 "sggbak.f"
	return 0;
#line 215 "sggbak.f"
    }
#line 217 "sggbak.f"
    if (*m == 0) {
#line 217 "sggbak.f"
	return 0;
#line 217 "sggbak.f"
    }
#line 219 "sggbak.f"
    if (lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 219 "sggbak.f"
	return 0;
#line 219 "sggbak.f"
    }

#line 222 "sggbak.f"
    if (*ilo == *ihi) {
#line 222 "sggbak.f"
	goto L30;
#line 222 "sggbak.f"
    }

/*     Backward balance */

#line 227 "sggbak.f"
    if (lsame_(job, "S", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (ftnlen)1, 
	    (ftnlen)1)) {

/*        Backward transformation on right eigenvectors */

#line 231 "sggbak.f"
	if (rightv) {
#line 232 "sggbak.f"
	    i__1 = *ihi;
#line 232 "sggbak.f"
	    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 233 "sggbak.f"
		sscal_(m, &rscale[i__], &v[i__ + v_dim1], ldv);
#line 234 "sggbak.f"
/* L10: */
#line 234 "sggbak.f"
	    }
#line 235 "sggbak.f"
	}

/*        Backward transformation on left eigenvectors */

#line 239 "sggbak.f"
	if (leftv) {
#line 240 "sggbak.f"
	    i__1 = *ihi;
#line 240 "sggbak.f"
	    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 241 "sggbak.f"
		sscal_(m, &lscale[i__], &v[i__ + v_dim1], ldv);
#line 242 "sggbak.f"
/* L20: */
#line 242 "sggbak.f"
	    }
#line 243 "sggbak.f"
	}
#line 244 "sggbak.f"
    }

/*     Backward permutation */

#line 248 "sggbak.f"
L30:
#line 249 "sggbak.f"
    if (lsame_(job, "P", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (ftnlen)1, 
	    (ftnlen)1)) {

/*        Backward permutation on right eigenvectors */

#line 253 "sggbak.f"
	if (rightv) {
#line 254 "sggbak.f"
	    if (*ilo == 1) {
#line 254 "sggbak.f"
		goto L50;
#line 254 "sggbak.f"
	    }

#line 257 "sggbak.f"
	    for (i__ = *ilo - 1; i__ >= 1; --i__) {
#line 258 "sggbak.f"
		k = (integer) rscale[i__];
#line 259 "sggbak.f"
		if (k == i__) {
#line 259 "sggbak.f"
		    goto L40;
#line 259 "sggbak.f"
		}
#line 261 "sggbak.f"
		sswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 262 "sggbak.f"
L40:
#line 262 "sggbak.f"
		;
#line 262 "sggbak.f"
	    }

#line 264 "sggbak.f"
L50:
#line 265 "sggbak.f"
	    if (*ihi == *n) {
#line 265 "sggbak.f"
		goto L70;
#line 265 "sggbak.f"
	    }
#line 267 "sggbak.f"
	    i__1 = *n;
#line 267 "sggbak.f"
	    for (i__ = *ihi + 1; i__ <= i__1; ++i__) {
#line 268 "sggbak.f"
		k = (integer) rscale[i__];
#line 269 "sggbak.f"
		if (k == i__) {
#line 269 "sggbak.f"
		    goto L60;
#line 269 "sggbak.f"
		}
#line 271 "sggbak.f"
		sswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 272 "sggbak.f"
L60:
#line 272 "sggbak.f"
		;
#line 272 "sggbak.f"
	    }
#line 273 "sggbak.f"
	}

/*        Backward permutation on left eigenvectors */

#line 277 "sggbak.f"
L70:
#line 278 "sggbak.f"
	if (leftv) {
#line 279 "sggbak.f"
	    if (*ilo == 1) {
#line 279 "sggbak.f"
		goto L90;
#line 279 "sggbak.f"
	    }
#line 281 "sggbak.f"
	    for (i__ = *ilo - 1; i__ >= 1; --i__) {
#line 282 "sggbak.f"
		k = (integer) lscale[i__];
#line 283 "sggbak.f"
		if (k == i__) {
#line 283 "sggbak.f"
		    goto L80;
#line 283 "sggbak.f"
		}
#line 285 "sggbak.f"
		sswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 286 "sggbak.f"
L80:
#line 286 "sggbak.f"
		;
#line 286 "sggbak.f"
	    }

#line 288 "sggbak.f"
L90:
#line 289 "sggbak.f"
	    if (*ihi == *n) {
#line 289 "sggbak.f"
		goto L110;
#line 289 "sggbak.f"
	    }
#line 291 "sggbak.f"
	    i__1 = *n;
#line 291 "sggbak.f"
	    for (i__ = *ihi + 1; i__ <= i__1; ++i__) {
#line 292 "sggbak.f"
		k = (integer) lscale[i__];
#line 293 "sggbak.f"
		if (k == i__) {
#line 293 "sggbak.f"
		    goto L100;
#line 293 "sggbak.f"
		}
#line 295 "sggbak.f"
		sswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 296 "sggbak.f"
L100:
#line 296 "sggbak.f"
		;
#line 296 "sggbak.f"
	    }
#line 297 "sggbak.f"
	}
#line 298 "sggbak.f"
    }

#line 300 "sggbak.f"
L110:

#line 302 "sggbak.f"
    return 0;

/*     End of SGGBAK */

} /* sggbak_ */

