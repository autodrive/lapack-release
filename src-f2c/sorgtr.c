#line 1 "sorgtr.f"
/* sorgtr.f -- translated by f2c (version 20100827).
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

#line 1 "sorgtr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b SORGTR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORGTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgtr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgtr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgtr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SORGTR generates a real orthogonal matrix Q which is defined as the */
/* > product of n-1 elementary reflectors of order N, as returned by */
/* > SSYTRD: */
/* > */
/* > if UPLO = 'U', Q = H(n-1) . . . H(2) H(1), */
/* > */
/* > if UPLO = 'L', Q = H(1) H(2) . . . H(n-1). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U': Upper triangle of A contains elementary reflectors */
/* >                 from SSYTRD; */
/* >          = 'L': Lower triangle of A contains elementary reflectors */
/* >                 from SSYTRD. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix Q. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the vectors which define the elementary reflectors, */
/* >          as returned by SSYTRD. */
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
/* >          reflector H(i), as returned by SSYTRD. */
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
/* >          The dimension of the array WORK. LWORK >= max(1,N-1). */
/* >          For optimum performance LWORK >= (N-1)*NB, where NB is */
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

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sorgtr_(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info,
	 ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, nb;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int sorgql_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), sorgqr_(integer *, integer *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *, integer *, integer *);
    static logical lquery;
    static integer lwkopt;


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

/*     Test the input arguments */

#line 164 "sorgtr.f"
    /* Parameter adjustments */
#line 164 "sorgtr.f"
    a_dim1 = *lda;
#line 164 "sorgtr.f"
    a_offset = 1 + a_dim1;
#line 164 "sorgtr.f"
    a -= a_offset;
#line 164 "sorgtr.f"
    --tau;
#line 164 "sorgtr.f"
    --work;
#line 164 "sorgtr.f"

#line 164 "sorgtr.f"
    /* Function Body */
#line 164 "sorgtr.f"
    *info = 0;
#line 165 "sorgtr.f"
    lquery = *lwork == -1;
#line 166 "sorgtr.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 167 "sorgtr.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 168 "sorgtr.f"
	*info = -1;
#line 169 "sorgtr.f"
    } else if (*n < 0) {
#line 170 "sorgtr.f"
	*info = -2;
#line 171 "sorgtr.f"
    } else if (*lda < max(1,*n)) {
#line 172 "sorgtr.f"
	*info = -4;
#line 173 "sorgtr.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 173 "sorgtr.f"
	i__1 = 1, i__2 = *n - 1;
#line 173 "sorgtr.f"
	if (*lwork < max(i__1,i__2) && ! lquery) {
#line 174 "sorgtr.f"
	    *info = -7;
#line 175 "sorgtr.f"
	}
#line 175 "sorgtr.f"
    }

#line 177 "sorgtr.f"
    if (*info == 0) {
#line 178 "sorgtr.f"
	if (upper) {
#line 179 "sorgtr.f"
	    i__1 = *n - 1;
#line 179 "sorgtr.f"
	    i__2 = *n - 1;
#line 179 "sorgtr.f"
	    i__3 = *n - 1;
#line 179 "sorgtr.f"
	    nb = ilaenv_(&c__1, "SORGQL", " ", &i__1, &i__2, &i__3, &c_n1, (
		    ftnlen)6, (ftnlen)1);
#line 180 "sorgtr.f"
	} else {
#line 181 "sorgtr.f"
	    i__1 = *n - 1;
#line 181 "sorgtr.f"
	    i__2 = *n - 1;
#line 181 "sorgtr.f"
	    i__3 = *n - 1;
#line 181 "sorgtr.f"
	    nb = ilaenv_(&c__1, "SORGQR", " ", &i__1, &i__2, &i__3, &c_n1, (
		    ftnlen)6, (ftnlen)1);
#line 182 "sorgtr.f"
	}
/* Computing MAX */
#line 183 "sorgtr.f"
	i__1 = 1, i__2 = *n - 1;
#line 183 "sorgtr.f"
	lwkopt = max(i__1,i__2) * nb;
#line 184 "sorgtr.f"
	work[1] = (doublereal) lwkopt;
#line 185 "sorgtr.f"
    }

#line 187 "sorgtr.f"
    if (*info != 0) {
#line 188 "sorgtr.f"
	i__1 = -(*info);
#line 188 "sorgtr.f"
	xerbla_("SORGTR", &i__1, (ftnlen)6);
#line 189 "sorgtr.f"
	return 0;
#line 190 "sorgtr.f"
    } else if (lquery) {
#line 191 "sorgtr.f"
	return 0;
#line 192 "sorgtr.f"
    }

/*     Quick return if possible */

#line 196 "sorgtr.f"
    if (*n == 0) {
#line 197 "sorgtr.f"
	work[1] = 1.;
#line 198 "sorgtr.f"
	return 0;
#line 199 "sorgtr.f"
    }

#line 201 "sorgtr.f"
    if (upper) {

/*        Q was determined by a call to SSYTRD with UPLO = 'U' */

/*        Shift the vectors which define the elementary reflectors one */
/*        column to the left, and set the last row and column of Q to */
/*        those of the unit matrix */

#line 209 "sorgtr.f"
	i__1 = *n - 1;
#line 209 "sorgtr.f"
	for (j = 1; j <= i__1; ++j) {
#line 210 "sorgtr.f"
	    i__2 = j - 1;
#line 210 "sorgtr.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 211 "sorgtr.f"
		a[i__ + j * a_dim1] = a[i__ + (j + 1) * a_dim1];
#line 212 "sorgtr.f"
/* L10: */
#line 212 "sorgtr.f"
	    }
#line 213 "sorgtr.f"
	    a[*n + j * a_dim1] = 0.;
#line 214 "sorgtr.f"
/* L20: */
#line 214 "sorgtr.f"
	}
#line 215 "sorgtr.f"
	i__1 = *n - 1;
#line 215 "sorgtr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 216 "sorgtr.f"
	    a[i__ + *n * a_dim1] = 0.;
#line 217 "sorgtr.f"
/* L30: */
#line 217 "sorgtr.f"
	}
#line 218 "sorgtr.f"
	a[*n + *n * a_dim1] = 1.;

/*        Generate Q(1:n-1,1:n-1) */

#line 222 "sorgtr.f"
	i__1 = *n - 1;
#line 222 "sorgtr.f"
	i__2 = *n - 1;
#line 222 "sorgtr.f"
	i__3 = *n - 1;
#line 222 "sorgtr.f"
	sorgql_(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &work[1], 
		lwork, &iinfo);

#line 224 "sorgtr.f"
    } else {

/*        Q was determined by a call to SSYTRD with UPLO = 'L'. */

/*        Shift the vectors which define the elementary reflectors one */
/*        column to the right, and set the first row and column of Q to */
/*        those of the unit matrix */

#line 232 "sorgtr.f"
	for (j = *n; j >= 2; --j) {
#line 233 "sorgtr.f"
	    a[j * a_dim1 + 1] = 0.;
#line 234 "sorgtr.f"
	    i__1 = *n;
#line 234 "sorgtr.f"
	    for (i__ = j + 1; i__ <= i__1; ++i__) {
#line 235 "sorgtr.f"
		a[i__ + j * a_dim1] = a[i__ + (j - 1) * a_dim1];
#line 236 "sorgtr.f"
/* L40: */
#line 236 "sorgtr.f"
	    }
#line 237 "sorgtr.f"
/* L50: */
#line 237 "sorgtr.f"
	}
#line 238 "sorgtr.f"
	a[a_dim1 + 1] = 1.;
#line 239 "sorgtr.f"
	i__1 = *n;
#line 239 "sorgtr.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 240 "sorgtr.f"
	    a[i__ + a_dim1] = 0.;
#line 241 "sorgtr.f"
/* L60: */
#line 241 "sorgtr.f"
	}
#line 242 "sorgtr.f"
	if (*n > 1) {

/*           Generate Q(2:n,2:n) */

#line 246 "sorgtr.f"
	    i__1 = *n - 1;
#line 246 "sorgtr.f"
	    i__2 = *n - 1;
#line 246 "sorgtr.f"
	    i__3 = *n - 1;
#line 246 "sorgtr.f"
	    sorgqr_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[1], 
		    &work[1], lwork, &iinfo);
#line 248 "sorgtr.f"
	}
#line 249 "sorgtr.f"
    }
#line 250 "sorgtr.f"
    work[1] = (doublereal) lwkopt;
#line 251 "sorgtr.f"
    return 0;

/*     End of SORGTR */

} /* sorgtr_ */

