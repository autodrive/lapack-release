#line 1 "dgebak.f"
/* dgebak.f -- translated by f2c (version 20100827).
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

#line 1 "dgebak.f"
/* > \brief \b DGEBAK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEBAK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgebak.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgebak.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgebak.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOB, SIDE */
/*       INTEGER            IHI, ILO, INFO, LDV, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   SCALE( * ), V( LDV, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEBAK forms the right or left eigenvectors of a real general matrix */
/* > by backward transformation on the computed eigenvectors of the */
/* > balanced matrix output by DGEBAL. */
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
/* >          JOB must be the same as the argument JOB supplied to DGEBAL. */
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
/* >          The integers ILO and IHI determined by DGEBAL. */
/* >          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in] SCALE */
/* > \verbatim */
/* >          SCALE is DOUBLE PRECISION array, dimension (N) */
/* >          Details of the permutation and scaling factors, as returned */
/* >          by DGEBAL. */
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
/* >          transformed, as returned by DHSEIN or DTREVC. */
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

/* > \ingroup doubleGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int dgebak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, doublereal *scale, integer *m, doublereal *v, integer *
	ldv, integer *info, ftnlen job_len, ftnlen side_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, i__1;

    /* Local variables */
    static integer i__, k;
    static doublereal s;
    static integer ii;
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

#line 171 "dgebak.f"
    /* Parameter adjustments */
#line 171 "dgebak.f"
    --scale;
#line 171 "dgebak.f"
    v_dim1 = *ldv;
#line 171 "dgebak.f"
    v_offset = 1 + v_dim1;
#line 171 "dgebak.f"
    v -= v_offset;
#line 171 "dgebak.f"

#line 171 "dgebak.f"
    /* Function Body */
#line 171 "dgebak.f"
    rightv = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 172 "dgebak.f"
    leftv = lsame_(side, "L", (ftnlen)1, (ftnlen)1);

#line 174 "dgebak.f"
    *info = 0;
#line 175 "dgebak.f"
    if (! lsame_(job, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(job, "P", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(job, "S", (ftnlen)1, (ftnlen)1) 
	    && ! lsame_(job, "B", (ftnlen)1, (ftnlen)1)) {
#line 177 "dgebak.f"
	*info = -1;
#line 178 "dgebak.f"
    } else if (! rightv && ! leftv) {
#line 179 "dgebak.f"
	*info = -2;
#line 180 "dgebak.f"
    } else if (*n < 0) {
#line 181 "dgebak.f"
	*info = -3;
#line 182 "dgebak.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 183 "dgebak.f"
	*info = -4;
#line 184 "dgebak.f"
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
#line 185 "dgebak.f"
	*info = -5;
#line 186 "dgebak.f"
    } else if (*m < 0) {
#line 187 "dgebak.f"
	*info = -7;
#line 188 "dgebak.f"
    } else if (*ldv < max(1,*n)) {
#line 189 "dgebak.f"
	*info = -9;
#line 190 "dgebak.f"
    }
#line 191 "dgebak.f"
    if (*info != 0) {
#line 192 "dgebak.f"
	i__1 = -(*info);
#line 192 "dgebak.f"
	xerbla_("DGEBAK", &i__1, (ftnlen)6);
#line 193 "dgebak.f"
	return 0;
#line 194 "dgebak.f"
    }

/*     Quick return if possible */

#line 198 "dgebak.f"
    if (*n == 0) {
#line 198 "dgebak.f"
	return 0;
#line 198 "dgebak.f"
    }
#line 200 "dgebak.f"
    if (*m == 0) {
#line 200 "dgebak.f"
	return 0;
#line 200 "dgebak.f"
    }
#line 202 "dgebak.f"
    if (lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 202 "dgebak.f"
	return 0;
#line 202 "dgebak.f"
    }

#line 205 "dgebak.f"
    if (*ilo == *ihi) {
#line 205 "dgebak.f"
	goto L30;
#line 205 "dgebak.f"
    }

/*     Backward balance */

#line 210 "dgebak.f"
    if (lsame_(job, "S", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (ftnlen)1, 
	    (ftnlen)1)) {

#line 212 "dgebak.f"
	if (rightv) {
#line 213 "dgebak.f"
	    i__1 = *ihi;
#line 213 "dgebak.f"
	    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 214 "dgebak.f"
		s = scale[i__];
#line 215 "dgebak.f"
		dscal_(m, &s, &v[i__ + v_dim1], ldv);
#line 216 "dgebak.f"
/* L10: */
#line 216 "dgebak.f"
	    }
#line 217 "dgebak.f"
	}

#line 219 "dgebak.f"
	if (leftv) {
#line 220 "dgebak.f"
	    i__1 = *ihi;
#line 220 "dgebak.f"
	    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 221 "dgebak.f"
		s = 1. / scale[i__];
#line 222 "dgebak.f"
		dscal_(m, &s, &v[i__ + v_dim1], ldv);
#line 223 "dgebak.f"
/* L20: */
#line 223 "dgebak.f"
	    }
#line 224 "dgebak.f"
	}

#line 226 "dgebak.f"
    }

/*     Backward permutation */

/*     For  I = ILO-1 step -1 until 1, */
/*              IHI+1 step 1 until N do -- */

#line 233 "dgebak.f"
L30:
#line 234 "dgebak.f"
    if (lsame_(job, "P", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (ftnlen)1, 
	    (ftnlen)1)) {
#line 235 "dgebak.f"
	if (rightv) {
#line 236 "dgebak.f"
	    i__1 = *n;
#line 236 "dgebak.f"
	    for (ii = 1; ii <= i__1; ++ii) {
#line 237 "dgebak.f"
		i__ = ii;
#line 238 "dgebak.f"
		if (i__ >= *ilo && i__ <= *ihi) {
#line 238 "dgebak.f"
		    goto L40;
#line 238 "dgebak.f"
		}
#line 240 "dgebak.f"
		if (i__ < *ilo) {
#line 240 "dgebak.f"
		    i__ = *ilo - ii;
#line 240 "dgebak.f"
		}
#line 242 "dgebak.f"
		k = (integer) scale[i__];
#line 243 "dgebak.f"
		if (k == i__) {
#line 243 "dgebak.f"
		    goto L40;
#line 243 "dgebak.f"
		}
#line 245 "dgebak.f"
		dswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 246 "dgebak.f"
L40:
#line 246 "dgebak.f"
		;
#line 246 "dgebak.f"
	    }
#line 247 "dgebak.f"
	}

#line 249 "dgebak.f"
	if (leftv) {
#line 250 "dgebak.f"
	    i__1 = *n;
#line 250 "dgebak.f"
	    for (ii = 1; ii <= i__1; ++ii) {
#line 251 "dgebak.f"
		i__ = ii;
#line 252 "dgebak.f"
		if (i__ >= *ilo && i__ <= *ihi) {
#line 252 "dgebak.f"
		    goto L50;
#line 252 "dgebak.f"
		}
#line 254 "dgebak.f"
		if (i__ < *ilo) {
#line 254 "dgebak.f"
		    i__ = *ilo - ii;
#line 254 "dgebak.f"
		}
#line 256 "dgebak.f"
		k = (integer) scale[i__];
#line 257 "dgebak.f"
		if (k == i__) {
#line 257 "dgebak.f"
		    goto L50;
#line 257 "dgebak.f"
		}
#line 259 "dgebak.f"
		dswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
#line 260 "dgebak.f"
L50:
#line 260 "dgebak.f"
		;
#line 260 "dgebak.f"
	    }
#line 261 "dgebak.f"
	}
#line 262 "dgebak.f"
    }

#line 264 "dgebak.f"
    return 0;

/*     End of DGEBAK */

} /* dgebak_ */

