#line 1 "sgebak.f"
/* sgebak.f -- translated by f2c (version 20100827).
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

#line 1 "sgebak.f"
/* > \brief \b SGEBAK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGEBAK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgebak.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgebak.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgebak.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOB, SIDE */
/*       INTEGER            IHI, ILO, INFO, LDV, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               V( LDV, * ), SCALE( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEBAK forms the right or left eigenvectors of a real general matrix */
/* > by backward transformation on the computed eigenvectors of the */
/* > balanced matrix output by SGEBAL. */
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
/* >          JOB must be the same as the argument JOB supplied to SGEBAL. */
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
/* >          The integers ILO and IHI determined by SGEBAL. */
/* >          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in] SCALE */
/* > \verbatim */
/* >          SCALE is REAL array, dimension (N) */
/* >          Details of the permutation and scaling factors, as returned */
/* >          by SGEBAL. */
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
/* >          transformed, as returned by SHSEIN or STREVC. */
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

/* > \ingroup realGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int sgebak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, doublereal *scale, integer *m, doublereal *v, integer *
	ldv, integer *info, ftnlen job_len, ftnlen side_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, i__1;

    /* Local variables */
    static integer i__, k;
    static doublereal s;
    static integer ii;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical leftv;
    extern /* Subroutine */ int sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), xerbla_(char *, integer *, ftnlen);
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

#line 171 "sgebak.f"
    /* Parameter adjustments */
#line 171 "sgebak.f"
    --scale;
#line 171 "sgebak.f"
    v_dim1 = *ldv;
#line 171 "sgebak.f"
    v_offset = 1 + v_dim1;
#line 171 "sgebak.f"
    v -= v_offset;
#line 171 "sgebak.f"

#line 171 "sgebak.f"
    /* Function Body */
#line 171 "sgebak.f"
    rightv = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 172 "sgebak.f"
    leftv = lsame_(side, "L", (ftnlen)1, (ftnlen)1);

#line 174 "sgebak.f"
    *info = 0;
#line 175 "sgebak.f"
    if (! lsame_(job, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(job, "P", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(job, "S", (ftnlen)1, (ftnlen)1) 
	    && ! lsame_(job, "B", (ftnlen)1, (ftnlen)1)) {
#line 177 "sgebak.f"
	*info = -1;
#line 178 "sgebak.f"
    } else if (! rightv && ! leftv) {
#line 179 "sgebak.f"
	*info = -2;
#line 180 "sgebak.f"
    } else if (*n < 0) {
#line 181 "sgebak.f"
	*info = -3;
#line 182 "sgebak.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 183 "sgebak.f"
	*info = -4;
#line 184 "sgebak.f"
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
#line 185 "sgebak.f"
	*info = -5;
#line 186 "sgebak.f"
    } else if (*m < 0) {
#line 187 "sgebak.f"
	*info = -7;
#line 188 "sgebak.f"
    } else if (*ldv < max(1,*n)) {
#line 189 "sgebak.f"
	*info = -9;
#line 190 "sgebak.f"
    }
#line 191 "sgebak.f"
    if (*info != 0) {
#line 192 "sgebak.f"
	i__1 = -(*info);
#line 192 "sgebak.f"
	xerbla_("SGEBAK", &i__1, (ftnlen)6);
#line 193 "sgebak.f"
	return 0;
#line 194 "sgebak.f"
    }

/*     Quick return if possible */

#line 198 "sgebak.f"
    if (*n == 0) {
#line 198 "sgebak.f"
	return 0;
#line 198 "sgebak.f"
    }
#line 200 "sgebak.f"
    if (*m == 0) {
#line 200 "sgebak.f"
	return 0;
#line 200 "sgebak.f"
    }
#line 202 "sgebak.f"
    if (lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 202 "sgebak.f"
	return 0;
#line 202 "sgebak.f"
    }

#line 205 "sgebak.f"
    if (*ilo == *ihi) {
#line 205 "sgebak.f"
	goto L30;
#line 205 "sgebak.f"
    }

/*     Backward balance */

#line 210 "sgebak.f"
    if (lsame_(job, "S", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (ftnlen)1, 
	    (ftnlen)1)) {

#line 212 "sgebak.f"
	if (rightv) {
#line 213 "sgebak.f"
	    i__1 = *ihi;
#line 213 "sgebak.f"
	    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 214 "sgebak.f"
		s = scale[i__];
#line 215 "sgebak.f"
		sscal_(m, &s, &v[i__ + v_dim1], ldv);
#line 216 "sgebak.f"
/* L10: */
#line 216 "sgebak.f"
	    }
#line 217 "sgebak.f"
	}

#line 219 "sgebak.f"
	if (leftv) {
#line 220 "sgebak.f"
	    i__1 = *ihi;
#line 220 "sgebak.f"
	    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 221 "sgebak.f"
		s = 1. / scale[i__];
#line 222 "sgebak.f"
		sscal_(m, &s, &v[i__ + v_dim1], ldv);
#line 223 "sgebak.f"
/* L20: */
#line 223 "sgebak.f"
	    }
#line 224 "sgebak.f"
	}

#line 226 "sgebak.f"
    }

/*     Backward permutation */

/*     For  I = ILO-1 step -1 until 1, */
/*              IHI+1 step 1 until N do -- */

#line 233 "sgebak.f"
L30:
#line 234 "sgebak.f"
    if (lsame_(job, "P", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (ftnlen)1, 
	    (ftnlen)1)) {
#line 235 "sgebak.f"
	if (rightv) {
#line 236 "sgebak.f"
	    i__1 = *n;
#line 236 "sgebak.f"
	    for (ii = 1; ii <= i__1; ++ii) {
#line 237 "sgebak.f"
		i__ = ii;
#line 238 "sgebak.f"
		if (i__ >= *ilo && i__ <= *ihi) {
#line 238 "sgebak.f"
		    goto L40;
#line 238 "sgebak.f"
		}
#line 240 "sgebak.f"
		if (i__ < *ilo) {
#line 240 "sgebak.f"
		    i__ = *ilo - ii;
#line 240 "sgebak.f"
		}
#line 242 "sgebak.f"
		k = (integer) scale[i__];
#line 243 "sgebak.f"
		if (k == i__) {
#line 243 "sgebak.f"
		    goto L40;
#line 243 "sgebak.f"
		}
#line 245 "sgebak.f"
		sswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 246 "sgebak.f"
L40:
#line 246 "sgebak.f"
		;
#line 246 "sgebak.f"
	    }
#line 247 "sgebak.f"
	}

#line 249 "sgebak.f"
	if (leftv) {
#line 250 "sgebak.f"
	    i__1 = *n;
#line 250 "sgebak.f"
	    for (ii = 1; ii <= i__1; ++ii) {
#line 251 "sgebak.f"
		i__ = ii;
#line 252 "sgebak.f"
		if (i__ >= *ilo && i__ <= *ihi) {
#line 252 "sgebak.f"
		    goto L50;
#line 252 "sgebak.f"
		}
#line 254 "sgebak.f"
		if (i__ < *ilo) {
#line 254 "sgebak.f"
		    i__ = *ilo - ii;
#line 254 "sgebak.f"
		}
#line 256 "sgebak.f"
		k = (integer) scale[i__];
#line 257 "sgebak.f"
		if (k == i__) {
#line 257 "sgebak.f"
		    goto L50;
#line 257 "sgebak.f"
		}
#line 259 "sgebak.f"
		sswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 260 "sgebak.f"
L50:
#line 260 "sgebak.f"
		;
#line 260 "sgebak.f"
	    }
#line 261 "sgebak.f"
	}
#line 262 "sgebak.f"
    }

#line 264 "sgebak.f"
    return 0;

/*     End of SGEBAK */

} /* sgebak_ */

