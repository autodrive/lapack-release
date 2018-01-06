#line 1 "sorghr.f"
/* sorghr.f -- translated by f2c (version 20100827).
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

#line 1 "sorghr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b SORGHR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORGHR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorghr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorghr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorghr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHI, ILO, INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SORGHR generates a real orthogonal matrix Q which is defined as the */
/* > product of IHI-ILO elementary reflectors of order N, as returned by */
/* > SGEHRD: */
/* > */
/* > Q = H(ilo) H(ilo+1) . . . H(ihi-1). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix Q. N >= 0. */
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
/* > */
/* >          ILO and IHI must have the same values as in the previous call */
/* >          of SGEHRD. Q is equal to the unit matrix except in the */
/* >          submatrix Q(ilo+1:ihi,ilo+1:ihi). */
/* >          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the vectors which define the elementary reflectors, */
/* >          as returned by SGEHRD. */
/* >          On exit, the N-by-N orthogonal matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is REAL array, dimension (N-1) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by SGEHRD. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. LWORK >= IHI-ILO. */
/* >          For optimum performance LWORK >= (IHI-ILO)*NB, where NB is */
/* >          the optimal blocksize. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sorghr_(integer *n, integer *ilo, integer *ihi, 
	doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, nb, nh, iinfo;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int sorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static integer lwkopt;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 165 "sorghr.f"
    /* Parameter adjustments */
#line 165 "sorghr.f"
    a_dim1 = *lda;
#line 165 "sorghr.f"
    a_offset = 1 + a_dim1;
#line 165 "sorghr.f"
    a -= a_offset;
#line 165 "sorghr.f"
    --tau;
#line 165 "sorghr.f"
    --work;
#line 165 "sorghr.f"

#line 165 "sorghr.f"
    /* Function Body */
#line 165 "sorghr.f"
    *info = 0;
#line 166 "sorghr.f"
    nh = *ihi - *ilo;
#line 167 "sorghr.f"
    lquery = *lwork == -1;
#line 168 "sorghr.f"
    if (*n < 0) {
#line 169 "sorghr.f"
	*info = -1;
#line 170 "sorghr.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 171 "sorghr.f"
	*info = -2;
#line 172 "sorghr.f"
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
#line 173 "sorghr.f"
	*info = -3;
#line 174 "sorghr.f"
    } else if (*lda < max(1,*n)) {
#line 175 "sorghr.f"
	*info = -5;
#line 176 "sorghr.f"
    } else if (*lwork < max(1,nh) && ! lquery) {
#line 177 "sorghr.f"
	*info = -8;
#line 178 "sorghr.f"
    }

#line 180 "sorghr.f"
    if (*info == 0) {
#line 181 "sorghr.f"
	nb = ilaenv_(&c__1, "SORGQR", " ", &nh, &nh, &nh, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 182 "sorghr.f"
	lwkopt = max(1,nh) * nb;
#line 183 "sorghr.f"
	work[1] = (doublereal) lwkopt;
#line 184 "sorghr.f"
    }

#line 186 "sorghr.f"
    if (*info != 0) {
#line 187 "sorghr.f"
	i__1 = -(*info);
#line 187 "sorghr.f"
	xerbla_("SORGHR", &i__1, (ftnlen)6);
#line 188 "sorghr.f"
	return 0;
#line 189 "sorghr.f"
    } else if (lquery) {
#line 190 "sorghr.f"
	return 0;
#line 191 "sorghr.f"
    }

/*     Quick return if possible */

#line 195 "sorghr.f"
    if (*n == 0) {
#line 196 "sorghr.f"
	work[1] = 1.;
#line 197 "sorghr.f"
	return 0;
#line 198 "sorghr.f"
    }

/*     Shift the vectors which define the elementary reflectors one */
/*     column to the right, and set the first ilo and the last n-ihi */
/*     rows and columns to those of the unit matrix */

#line 204 "sorghr.f"
    i__1 = *ilo + 1;
#line 204 "sorghr.f"
    for (j = *ihi; j >= i__1; --j) {
#line 205 "sorghr.f"
	i__2 = j - 1;
#line 205 "sorghr.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 206 "sorghr.f"
	    a[i__ + j * a_dim1] = 0.;
#line 207 "sorghr.f"
/* L10: */
#line 207 "sorghr.f"
	}
#line 208 "sorghr.f"
	i__2 = *ihi;
#line 208 "sorghr.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 209 "sorghr.f"
	    a[i__ + j * a_dim1] = a[i__ + (j - 1) * a_dim1];
#line 210 "sorghr.f"
/* L20: */
#line 210 "sorghr.f"
	}
#line 211 "sorghr.f"
	i__2 = *n;
#line 211 "sorghr.f"
	for (i__ = *ihi + 1; i__ <= i__2; ++i__) {
#line 212 "sorghr.f"
	    a[i__ + j * a_dim1] = 0.;
#line 213 "sorghr.f"
/* L30: */
#line 213 "sorghr.f"
	}
#line 214 "sorghr.f"
/* L40: */
#line 214 "sorghr.f"
    }
#line 215 "sorghr.f"
    i__1 = *ilo;
#line 215 "sorghr.f"
    for (j = 1; j <= i__1; ++j) {
#line 216 "sorghr.f"
	i__2 = *n;
#line 216 "sorghr.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 217 "sorghr.f"
	    a[i__ + j * a_dim1] = 0.;
#line 218 "sorghr.f"
/* L50: */
#line 218 "sorghr.f"
	}
#line 219 "sorghr.f"
	a[j + j * a_dim1] = 1.;
#line 220 "sorghr.f"
/* L60: */
#line 220 "sorghr.f"
    }
#line 221 "sorghr.f"
    i__1 = *n;
#line 221 "sorghr.f"
    for (j = *ihi + 1; j <= i__1; ++j) {
#line 222 "sorghr.f"
	i__2 = *n;
#line 222 "sorghr.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 223 "sorghr.f"
	    a[i__ + j * a_dim1] = 0.;
#line 224 "sorghr.f"
/* L70: */
#line 224 "sorghr.f"
	}
#line 225 "sorghr.f"
	a[j + j * a_dim1] = 1.;
#line 226 "sorghr.f"
/* L80: */
#line 226 "sorghr.f"
    }

#line 228 "sorghr.f"
    if (nh > 0) {

/*        Generate Q(ilo+1:ihi,ilo+1:ihi) */

#line 232 "sorghr.f"
	sorgqr_(&nh, &nh, &nh, &a[*ilo + 1 + (*ilo + 1) * a_dim1], lda, &tau[*
		ilo], &work[1], lwork, &iinfo);
#line 234 "sorghr.f"
    }
#line 235 "sorghr.f"
    work[1] = (doublereal) lwkopt;
#line 236 "sorghr.f"
    return 0;

/*     End of SORGHR */

} /* sorghr_ */

