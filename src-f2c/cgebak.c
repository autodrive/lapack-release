#line 1 "cgebak.f"
/* cgebak.f -- translated by f2c (version 20100827).
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

#line 1 "cgebak.f"
/* > \brief \b CGEBAK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGEBAK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgebak.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgebak.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgebak.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOB, SIDE */
/*       INTEGER            IHI, ILO, INFO, LDV, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               SCALE( * ) */
/*       COMPLEX            V( LDV, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEBAK forms the right or left eigenvectors of a complex general */
/* > matrix by backward transformation on the computed eigenvectors of the */
/* > balanced matrix output by CGEBAL. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOB */
/* > \verbatim */
/* >          JOB is CHARACTER*1 */
/* >          Specifies the type of backward transformation required: */
/* >          = 'N', do nothing, return immediately; */
/* >          = 'P', do backward transformation for permutation only; */
/* >          = 'S', do backward transformation for scaling only; */
/* >          = 'B', do backward transformations for both permutation and */
/* >                 scaling. */
/* >          JOB must be the same as the argument JOB supplied to CGEBAL. */
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
/* >          The integers ILO and IHI determined by CGEBAL. */
/* >          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in] SCALE */
/* > \verbatim */
/* >          SCALE is REAL array, dimension (N) */
/* >          Details of the permutation and scaling factors, as returned */
/* >          by CGEBAL. */
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
/* >          transformed, as returned by CHSEIN or CTREVC. */
/* >          On exit, V is overwritten by the transformed eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V. LDV >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int cgebak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, doublereal *scale, integer *m, doublecomplex *v, 
	integer *ldv, integer *info, ftnlen job_len, ftnlen side_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, i__1;

    /* Local variables */
    static integer i__, k;
    static doublereal s;
    static integer ii;
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

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test the input parameters */

#line 173 "cgebak.f"
    /* Parameter adjustments */
#line 173 "cgebak.f"
    --scale;
#line 173 "cgebak.f"
    v_dim1 = *ldv;
#line 173 "cgebak.f"
    v_offset = 1 + v_dim1;
#line 173 "cgebak.f"
    v -= v_offset;
#line 173 "cgebak.f"

#line 173 "cgebak.f"
    /* Function Body */
#line 173 "cgebak.f"
    rightv = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 174 "cgebak.f"
    leftv = lsame_(side, "L", (ftnlen)1, (ftnlen)1);

#line 176 "cgebak.f"
    *info = 0;
#line 177 "cgebak.f"
    if (! lsame_(job, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(job, "P", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(job, "S", (ftnlen)1, (ftnlen)1) 
	    && ! lsame_(job, "B", (ftnlen)1, (ftnlen)1)) {
#line 179 "cgebak.f"
	*info = -1;
#line 180 "cgebak.f"
    } else if (! rightv && ! leftv) {
#line 181 "cgebak.f"
	*info = -2;
#line 182 "cgebak.f"
    } else if (*n < 0) {
#line 183 "cgebak.f"
	*info = -3;
#line 184 "cgebak.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 185 "cgebak.f"
	*info = -4;
#line 186 "cgebak.f"
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
#line 187 "cgebak.f"
	*info = -5;
#line 188 "cgebak.f"
    } else if (*m < 0) {
#line 189 "cgebak.f"
	*info = -7;
#line 190 "cgebak.f"
    } else if (*ldv < max(1,*n)) {
#line 191 "cgebak.f"
	*info = -9;
#line 192 "cgebak.f"
    }
#line 193 "cgebak.f"
    if (*info != 0) {
#line 194 "cgebak.f"
	i__1 = -(*info);
#line 194 "cgebak.f"
	xerbla_("CGEBAK", &i__1, (ftnlen)6);
#line 195 "cgebak.f"
	return 0;
#line 196 "cgebak.f"
    }

/*     Quick return if possible */

#line 200 "cgebak.f"
    if (*n == 0) {
#line 200 "cgebak.f"
	return 0;
#line 200 "cgebak.f"
    }
#line 202 "cgebak.f"
    if (*m == 0) {
#line 202 "cgebak.f"
	return 0;
#line 202 "cgebak.f"
    }
#line 204 "cgebak.f"
    if (lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 204 "cgebak.f"
	return 0;
#line 204 "cgebak.f"
    }

#line 207 "cgebak.f"
    if (*ilo == *ihi) {
#line 207 "cgebak.f"
	goto L30;
#line 207 "cgebak.f"
    }

/*     Backward balance */

#line 212 "cgebak.f"
    if (lsame_(job, "S", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (ftnlen)1, 
	    (ftnlen)1)) {

#line 214 "cgebak.f"
	if (rightv) {
#line 215 "cgebak.f"
	    i__1 = *ihi;
#line 215 "cgebak.f"
	    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 216 "cgebak.f"
		s = scale[i__];
#line 217 "cgebak.f"
		csscal_(m, &s, &v[i__ + v_dim1], ldv);
#line 218 "cgebak.f"
/* L10: */
#line 218 "cgebak.f"
	    }
#line 219 "cgebak.f"
	}

#line 221 "cgebak.f"
	if (leftv) {
#line 222 "cgebak.f"
	    i__1 = *ihi;
#line 222 "cgebak.f"
	    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 223 "cgebak.f"
		s = 1. / scale[i__];
#line 224 "cgebak.f"
		csscal_(m, &s, &v[i__ + v_dim1], ldv);
#line 225 "cgebak.f"
/* L20: */
#line 225 "cgebak.f"
	    }
#line 226 "cgebak.f"
	}

#line 228 "cgebak.f"
    }

/*     Backward permutation */

/*     For  I = ILO-1 step -1 until 1, */
/*              IHI+1 step 1 until N do -- */

#line 235 "cgebak.f"
L30:
#line 236 "cgebak.f"
    if (lsame_(job, "P", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (ftnlen)1, 
	    (ftnlen)1)) {
#line 237 "cgebak.f"
	if (rightv) {
#line 238 "cgebak.f"
	    i__1 = *n;
#line 238 "cgebak.f"
	    for (ii = 1; ii <= i__1; ++ii) {
#line 239 "cgebak.f"
		i__ = ii;
#line 240 "cgebak.f"
		if (i__ >= *ilo && i__ <= *ihi) {
#line 240 "cgebak.f"
		    goto L40;
#line 240 "cgebak.f"
		}
#line 242 "cgebak.f"
		if (i__ < *ilo) {
#line 242 "cgebak.f"
		    i__ = *ilo - ii;
#line 242 "cgebak.f"
		}
#line 244 "cgebak.f"
		k = (integer) scale[i__];
#line 245 "cgebak.f"
		if (k == i__) {
#line 245 "cgebak.f"
		    goto L40;
#line 245 "cgebak.f"
		}
#line 247 "cgebak.f"
		cswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 248 "cgebak.f"
L40:
#line 248 "cgebak.f"
		;
#line 248 "cgebak.f"
	    }
#line 249 "cgebak.f"
	}

#line 251 "cgebak.f"
	if (leftv) {
#line 252 "cgebak.f"
	    i__1 = *n;
#line 252 "cgebak.f"
	    for (ii = 1; ii <= i__1; ++ii) {
#line 253 "cgebak.f"
		i__ = ii;
#line 254 "cgebak.f"
		if (i__ >= *ilo && i__ <= *ihi) {
#line 254 "cgebak.f"
		    goto L50;
#line 254 "cgebak.f"
		}
#line 256 "cgebak.f"
		if (i__ < *ilo) {
#line 256 "cgebak.f"
		    i__ = *ilo - ii;
#line 256 "cgebak.f"
		}
#line 258 "cgebak.f"
		k = (integer) scale[i__];
#line 259 "cgebak.f"
		if (k == i__) {
#line 259 "cgebak.f"
		    goto L50;
#line 259 "cgebak.f"
		}
#line 261 "cgebak.f"
		cswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 262 "cgebak.f"
L50:
#line 262 "cgebak.f"
		;
#line 262 "cgebak.f"
	    }
#line 263 "cgebak.f"
	}
#line 264 "cgebak.f"
    }

#line 266 "cgebak.f"
    return 0;

/*     End of CGEBAK */

} /* cgebak_ */

