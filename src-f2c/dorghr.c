#line 1 "dorghr.f"
/* dorghr.f -- translated by f2c (version 20100827).
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

#line 1 "dorghr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b DORGHR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORGHR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorghr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorghr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorghr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHI, ILO, INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DORGHR generates a real orthogonal matrix Q which is defined as the */
/* > product of IHI-ILO elementary reflectors of order N, as returned by */
/* > DGEHRD: */
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
/* >          of DGEHRD. Q is equal to the unit matrix except in the */
/* >          submatrix Q(ilo+1:ihi,ilo+1:ihi). */
/* >          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the vectors which define the elementary reflectors, */
/* >          as returned by DGEHRD. */
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
/* >          TAU is DOUBLE PRECISION array, dimension (N-1) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by DGEHRD. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
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

/* > \date December 2016 */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dorghr_(integer *n, integer *ilo, integer *ihi, 
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
    extern /* Subroutine */ int dorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static integer lwkopt;
    static logical lquery;


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
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 165 "dorghr.f"
    /* Parameter adjustments */
#line 165 "dorghr.f"
    a_dim1 = *lda;
#line 165 "dorghr.f"
    a_offset = 1 + a_dim1;
#line 165 "dorghr.f"
    a -= a_offset;
#line 165 "dorghr.f"
    --tau;
#line 165 "dorghr.f"
    --work;
#line 165 "dorghr.f"

#line 165 "dorghr.f"
    /* Function Body */
#line 165 "dorghr.f"
    *info = 0;
#line 166 "dorghr.f"
    nh = *ihi - *ilo;
#line 167 "dorghr.f"
    lquery = *lwork == -1;
#line 168 "dorghr.f"
    if (*n < 0) {
#line 169 "dorghr.f"
	*info = -1;
#line 170 "dorghr.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 171 "dorghr.f"
	*info = -2;
#line 172 "dorghr.f"
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
#line 173 "dorghr.f"
	*info = -3;
#line 174 "dorghr.f"
    } else if (*lda < max(1,*n)) {
#line 175 "dorghr.f"
	*info = -5;
#line 176 "dorghr.f"
    } else if (*lwork < max(1,nh) && ! lquery) {
#line 177 "dorghr.f"
	*info = -8;
#line 178 "dorghr.f"
    }

#line 180 "dorghr.f"
    if (*info == 0) {
#line 181 "dorghr.f"
	nb = ilaenv_(&c__1, "DORGQR", " ", &nh, &nh, &nh, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 182 "dorghr.f"
	lwkopt = max(1,nh) * nb;
#line 183 "dorghr.f"
	work[1] = (doublereal) lwkopt;
#line 184 "dorghr.f"
    }

#line 186 "dorghr.f"
    if (*info != 0) {
#line 187 "dorghr.f"
	i__1 = -(*info);
#line 187 "dorghr.f"
	xerbla_("DORGHR", &i__1, (ftnlen)6);
#line 188 "dorghr.f"
	return 0;
#line 189 "dorghr.f"
    } else if (lquery) {
#line 190 "dorghr.f"
	return 0;
#line 191 "dorghr.f"
    }

/*     Quick return if possible */

#line 195 "dorghr.f"
    if (*n == 0) {
#line 196 "dorghr.f"
	work[1] = 1.;
#line 197 "dorghr.f"
	return 0;
#line 198 "dorghr.f"
    }

/*     Shift the vectors which define the elementary reflectors one */
/*     column to the right, and set the first ilo and the last n-ihi */
/*     rows and columns to those of the unit matrix */

#line 204 "dorghr.f"
    i__1 = *ilo + 1;
#line 204 "dorghr.f"
    for (j = *ihi; j >= i__1; --j) {
#line 205 "dorghr.f"
	i__2 = j - 1;
#line 205 "dorghr.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 206 "dorghr.f"
	    a[i__ + j * a_dim1] = 0.;
#line 207 "dorghr.f"
/* L10: */
#line 207 "dorghr.f"
	}
#line 208 "dorghr.f"
	i__2 = *ihi;
#line 208 "dorghr.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 209 "dorghr.f"
	    a[i__ + j * a_dim1] = a[i__ + (j - 1) * a_dim1];
#line 210 "dorghr.f"
/* L20: */
#line 210 "dorghr.f"
	}
#line 211 "dorghr.f"
	i__2 = *n;
#line 211 "dorghr.f"
	for (i__ = *ihi + 1; i__ <= i__2; ++i__) {
#line 212 "dorghr.f"
	    a[i__ + j * a_dim1] = 0.;
#line 213 "dorghr.f"
/* L30: */
#line 213 "dorghr.f"
	}
#line 214 "dorghr.f"
/* L40: */
#line 214 "dorghr.f"
    }
#line 215 "dorghr.f"
    i__1 = *ilo;
#line 215 "dorghr.f"
    for (j = 1; j <= i__1; ++j) {
#line 216 "dorghr.f"
	i__2 = *n;
#line 216 "dorghr.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 217 "dorghr.f"
	    a[i__ + j * a_dim1] = 0.;
#line 218 "dorghr.f"
/* L50: */
#line 218 "dorghr.f"
	}
#line 219 "dorghr.f"
	a[j + j * a_dim1] = 1.;
#line 220 "dorghr.f"
/* L60: */
#line 220 "dorghr.f"
    }
#line 221 "dorghr.f"
    i__1 = *n;
#line 221 "dorghr.f"
    for (j = *ihi + 1; j <= i__1; ++j) {
#line 222 "dorghr.f"
	i__2 = *n;
#line 222 "dorghr.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 223 "dorghr.f"
	    a[i__ + j * a_dim1] = 0.;
#line 224 "dorghr.f"
/* L70: */
#line 224 "dorghr.f"
	}
#line 225 "dorghr.f"
	a[j + j * a_dim1] = 1.;
#line 226 "dorghr.f"
/* L80: */
#line 226 "dorghr.f"
    }

#line 228 "dorghr.f"
    if (nh > 0) {

/*        Generate Q(ilo+1:ihi,ilo+1:ihi) */

#line 232 "dorghr.f"
	dorgqr_(&nh, &nh, &nh, &a[*ilo + 1 + (*ilo + 1) * a_dim1], lda, &tau[*
		ilo], &work[1], lwork, &iinfo);
#line 234 "dorghr.f"
    }
#line 235 "dorghr.f"
    work[1] = (doublereal) lwkopt;
#line 236 "dorghr.f"
    return 0;

/*     End of DORGHR */

} /* dorghr_ */

