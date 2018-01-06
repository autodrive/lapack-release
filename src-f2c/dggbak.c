#line 1 "dggbak.f"
/* dggbak.f -- translated by f2c (version 20100827).
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

#line 1 "dggbak.f"
/* > \brief \b DGGBAK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGGBAK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggbak.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggbak.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggbak.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V, */
/*                          LDV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOB, SIDE */
/*       INTEGER            IHI, ILO, INFO, LDV, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   LSCALE( * ), RSCALE( * ), V( LDV, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGGBAK forms the right or left eigenvectors of a real generalized */
/* > eigenvalue problem A*x = lambda*B*x, by backward transformation on */
/* > the computed eigenvectors of the balanced pair of matrices output by */
/* > DGGBAL. */
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
/* >          JOB must be the same as the argument JOB supplied to DGGBAL. */
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
/* >          The integers ILO and IHI determined by DGGBAL. */
/* >          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in] LSCALE */
/* > \verbatim */
/* >          LSCALE is DOUBLE PRECISION array, dimension (N) */
/* >          Details of the permutations and/or scaling factors applied */
/* >          to the left side of A and B, as returned by DGGBAL. */
/* > \endverbatim */
/* > */
/* > \param[in] RSCALE */
/* > \verbatim */
/* >          RSCALE is DOUBLE PRECISION array, dimension (N) */
/* >          Details of the permutations and/or scaling factors applied */
/* >          to the right side of A and B, as returned by DGGBAL. */
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
/* >          V is DOUBLE PRECISION array, dimension (LDV,M) */
/* >          On entry, the matrix of right or left eigenvectors to be */
/* >          transformed, as returned by DTGEVC. */
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

/* > \ingroup doubleGBcomputational */

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
/* Subroutine */ int dggbak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, doublereal *lscale, doublereal *rscale, integer *m, 
	doublereal *v, integer *ldv, integer *info, ftnlen job_len, ftnlen 
	side_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, i__1;

    /* Local variables */
    static integer i__, k;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical leftv;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
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

#line 183 "dggbak.f"
    /* Parameter adjustments */
#line 183 "dggbak.f"
    --lscale;
#line 183 "dggbak.f"
    --rscale;
#line 183 "dggbak.f"
    v_dim1 = *ldv;
#line 183 "dggbak.f"
    v_offset = 1 + v_dim1;
#line 183 "dggbak.f"
    v -= v_offset;
#line 183 "dggbak.f"

#line 183 "dggbak.f"
    /* Function Body */
#line 183 "dggbak.f"
    rightv = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 184 "dggbak.f"
    leftv = lsame_(side, "L", (ftnlen)1, (ftnlen)1);

#line 186 "dggbak.f"
    *info = 0;
#line 187 "dggbak.f"
    if (! lsame_(job, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(job, "P", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(job, "S", (ftnlen)1, (ftnlen)1) 
	    && ! lsame_(job, "B", (ftnlen)1, (ftnlen)1)) {
#line 189 "dggbak.f"
	*info = -1;
#line 190 "dggbak.f"
    } else if (! rightv && ! leftv) {
#line 191 "dggbak.f"
	*info = -2;
#line 192 "dggbak.f"
    } else if (*n < 0) {
#line 193 "dggbak.f"
	*info = -3;
#line 194 "dggbak.f"
    } else if (*ilo < 1) {
#line 195 "dggbak.f"
	*info = -4;
#line 196 "dggbak.f"
    } else if (*n == 0 && *ihi == 0 && *ilo != 1) {
#line 197 "dggbak.f"
	*info = -4;
#line 198 "dggbak.f"
    } else if (*n > 0 && (*ihi < *ilo || *ihi > max(1,*n))) {
#line 200 "dggbak.f"
	*info = -5;
#line 201 "dggbak.f"
    } else if (*n == 0 && *ilo == 1 && *ihi != 0) {
#line 202 "dggbak.f"
	*info = -5;
#line 203 "dggbak.f"
    } else if (*m < 0) {
#line 204 "dggbak.f"
	*info = -8;
#line 205 "dggbak.f"
    } else if (*ldv < max(1,*n)) {
#line 206 "dggbak.f"
	*info = -10;
#line 207 "dggbak.f"
    }
#line 208 "dggbak.f"
    if (*info != 0) {
#line 209 "dggbak.f"
	i__1 = -(*info);
#line 209 "dggbak.f"
	xerbla_("DGGBAK", &i__1, (ftnlen)6);
#line 210 "dggbak.f"
	return 0;
#line 211 "dggbak.f"
    }

/*     Quick return if possible */

#line 215 "dggbak.f"
    if (*n == 0) {
#line 215 "dggbak.f"
	return 0;
#line 215 "dggbak.f"
    }
#line 217 "dggbak.f"
    if (*m == 0) {
#line 217 "dggbak.f"
	return 0;
#line 217 "dggbak.f"
    }
#line 219 "dggbak.f"
    if (lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 219 "dggbak.f"
	return 0;
#line 219 "dggbak.f"
    }

#line 222 "dggbak.f"
    if (*ilo == *ihi) {
#line 222 "dggbak.f"
	goto L30;
#line 222 "dggbak.f"
    }

/*     Backward balance */

#line 227 "dggbak.f"
    if (lsame_(job, "S", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (ftnlen)1, 
	    (ftnlen)1)) {

/*        Backward transformation on right eigenvectors */

#line 231 "dggbak.f"
	if (rightv) {
#line 232 "dggbak.f"
	    i__1 = *ihi;
#line 232 "dggbak.f"
	    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 233 "dggbak.f"
		dscal_(m, &rscale[i__], &v[i__ + v_dim1], ldv);
#line 234 "dggbak.f"
/* L10: */
#line 234 "dggbak.f"
	    }
#line 235 "dggbak.f"
	}

/*        Backward transformation on left eigenvectors */

#line 239 "dggbak.f"
	if (leftv) {
#line 240 "dggbak.f"
	    i__1 = *ihi;
#line 240 "dggbak.f"
	    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 241 "dggbak.f"
		dscal_(m, &lscale[i__], &v[i__ + v_dim1], ldv);
#line 242 "dggbak.f"
/* L20: */
#line 242 "dggbak.f"
	    }
#line 243 "dggbak.f"
	}
#line 244 "dggbak.f"
    }

/*     Backward permutation */

#line 248 "dggbak.f"
L30:
#line 249 "dggbak.f"
    if (lsame_(job, "P", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (ftnlen)1, 
	    (ftnlen)1)) {

/*        Backward permutation on right eigenvectors */

#line 253 "dggbak.f"
	if (rightv) {
#line 254 "dggbak.f"
	    if (*ilo == 1) {
#line 254 "dggbak.f"
		goto L50;
#line 254 "dggbak.f"
	    }

#line 257 "dggbak.f"
	    for (i__ = *ilo - 1; i__ >= 1; --i__) {
#line 258 "dggbak.f"
		k = (integer) rscale[i__];
#line 259 "dggbak.f"
		if (k == i__) {
#line 259 "dggbak.f"
		    goto L40;
#line 259 "dggbak.f"
		}
#line 261 "dggbak.f"
		dswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 262 "dggbak.f"
L40:
#line 262 "dggbak.f"
		;
#line 262 "dggbak.f"
	    }

#line 264 "dggbak.f"
L50:
#line 265 "dggbak.f"
	    if (*ihi == *n) {
#line 265 "dggbak.f"
		goto L70;
#line 265 "dggbak.f"
	    }
#line 267 "dggbak.f"
	    i__1 = *n;
#line 267 "dggbak.f"
	    for (i__ = *ihi + 1; i__ <= i__1; ++i__) {
#line 268 "dggbak.f"
		k = (integer) rscale[i__];
#line 269 "dggbak.f"
		if (k == i__) {
#line 269 "dggbak.f"
		    goto L60;
#line 269 "dggbak.f"
		}
#line 271 "dggbak.f"
		dswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 272 "dggbak.f"
L60:
#line 272 "dggbak.f"
		;
#line 272 "dggbak.f"
	    }
#line 273 "dggbak.f"
	}

/*        Backward permutation on left eigenvectors */

#line 277 "dggbak.f"
L70:
#line 278 "dggbak.f"
	if (leftv) {
#line 279 "dggbak.f"
	    if (*ilo == 1) {
#line 279 "dggbak.f"
		goto L90;
#line 279 "dggbak.f"
	    }
#line 281 "dggbak.f"
	    for (i__ = *ilo - 1; i__ >= 1; --i__) {
#line 282 "dggbak.f"
		k = (integer) lscale[i__];
#line 283 "dggbak.f"
		if (k == i__) {
#line 283 "dggbak.f"
		    goto L80;
#line 283 "dggbak.f"
		}
#line 285 "dggbak.f"
		dswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 286 "dggbak.f"
L80:
#line 286 "dggbak.f"
		;
#line 286 "dggbak.f"
	    }

#line 288 "dggbak.f"
L90:
#line 289 "dggbak.f"
	    if (*ihi == *n) {
#line 289 "dggbak.f"
		goto L110;
#line 289 "dggbak.f"
	    }
#line 291 "dggbak.f"
	    i__1 = *n;
#line 291 "dggbak.f"
	    for (i__ = *ihi + 1; i__ <= i__1; ++i__) {
#line 292 "dggbak.f"
		k = (integer) lscale[i__];
#line 293 "dggbak.f"
		if (k == i__) {
#line 293 "dggbak.f"
		    goto L100;
#line 293 "dggbak.f"
		}
#line 295 "dggbak.f"
		dswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 296 "dggbak.f"
L100:
#line 296 "dggbak.f"
		;
#line 296 "dggbak.f"
	    }
#line 297 "dggbak.f"
	}
#line 298 "dggbak.f"
    }

#line 300 "dggbak.f"
L110:

#line 302 "dggbak.f"
    return 0;

/*     End of DGGBAK */

} /* dggbak_ */

