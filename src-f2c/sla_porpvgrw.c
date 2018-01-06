#line 1 "sla_porpvgrw.f"
/* sla_porpvgrw.f -- translated by f2c (version 20100827).
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

#line 1 "sla_porpvgrw.f"
/* > \brief \b SLA_PORPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a symmetric or Her
mitian positive-definite matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLA_PORPVGRW + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_por
pvgrw.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_por
pvgrw.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_por
pvgrw.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION SLA_PORPVGRW( UPLO, NCOLS, A, LDA, AF, LDAF, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER*1        UPLO */
/*       INTEGER            NCOLS, LDA, LDAF */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), AF( LDAF, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > */
/* > SLA_PORPVGRW computes the reciprocal pivot growth factor */
/* > norm(A)/norm(U). The "max absolute element" norm is used. If this is */
/* > much less than 1, the stability of the LU factorization of the */
/* > (equilibrated) matrix A could be poor. This also means that the */
/* > solution X, estimated condition numbers, and error bounds could be */
/* > unreliable. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >       = 'U':  Upper triangle of A is stored; */
/* >       = 'L':  Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] NCOLS */
/* > \verbatim */
/* >          NCOLS is INTEGER */
/* >     The number of columns of the matrix A. NCOLS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >     On entry, the N-by-N matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >     The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] AF */
/* > \verbatim */
/* >          AF is REAL array, dimension (LDAF,N) */
/* >     The triangular factor U or L from the Cholesky factorization */
/* >     A = U**T*U or A = L*L**T, as computed by SPOTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* >          LDAF is INTEGER */
/* >     The leading dimension of the array AF.  LDAF >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (2*N) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realPOcomputational */

/*  ===================================================================== */
doublereal sla_porpvgrw__(char *uplo, integer *ncols, doublereal *a, integer *
	lda, doublereal *af, integer *ldaf, doublereal *work, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j;
    static doublereal amax, umax;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    static doublereal rpvgrw;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 136 "sla_porpvgrw.f"
    /* Parameter adjustments */
#line 136 "sla_porpvgrw.f"
    a_dim1 = *lda;
#line 136 "sla_porpvgrw.f"
    a_offset = 1 + a_dim1;
#line 136 "sla_porpvgrw.f"
    a -= a_offset;
#line 136 "sla_porpvgrw.f"
    af_dim1 = *ldaf;
#line 136 "sla_porpvgrw.f"
    af_offset = 1 + af_dim1;
#line 136 "sla_porpvgrw.f"
    af -= af_offset;
#line 136 "sla_porpvgrw.f"
    --work;
#line 136 "sla_porpvgrw.f"

#line 136 "sla_porpvgrw.f"
    /* Function Body */
#line 136 "sla_porpvgrw.f"
    upper = lsame_("Upper", uplo, (ftnlen)5, (ftnlen)1);

/*     SPOTRF will have factored only the NCOLSxNCOLS leading minor, so */
/*     we restrict the growth search to that minor and use only the first */
/*     2*NCOLS workspace entries. */

#line 142 "sla_porpvgrw.f"
    rpvgrw = 1.;
#line 143 "sla_porpvgrw.f"
    i__1 = *ncols << 1;
#line 143 "sla_porpvgrw.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 144 "sla_porpvgrw.f"
	work[i__] = 0.;
#line 145 "sla_porpvgrw.f"
    }

/*     Find the max magnitude entry of each column. */

#line 149 "sla_porpvgrw.f"
    if (upper) {
#line 150 "sla_porpvgrw.f"
	i__1 = *ncols;
#line 150 "sla_porpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 151 "sla_porpvgrw.f"
	    i__2 = j;
#line 151 "sla_porpvgrw.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 152 "sla_porpvgrw.f"
		d__2 = (d__1 = a[i__ + j * a_dim1], abs(d__1)), d__3 = work[*
			ncols + j];
#line 152 "sla_porpvgrw.f"
		work[*ncols + j] = max(d__2,d__3);
#line 154 "sla_porpvgrw.f"
	    }
#line 155 "sla_porpvgrw.f"
	}
#line 156 "sla_porpvgrw.f"
    } else {
#line 157 "sla_porpvgrw.f"
	i__1 = *ncols;
#line 157 "sla_porpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 158 "sla_porpvgrw.f"
	    i__2 = *ncols;
#line 158 "sla_porpvgrw.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 159 "sla_porpvgrw.f"
		d__2 = (d__1 = a[i__ + j * a_dim1], abs(d__1)), d__3 = work[*
			ncols + j];
#line 159 "sla_porpvgrw.f"
		work[*ncols + j] = max(d__2,d__3);
#line 161 "sla_porpvgrw.f"
	    }
#line 162 "sla_porpvgrw.f"
	}
#line 163 "sla_porpvgrw.f"
    }

/*     Now find the max magnitude entry of each column of the factor in */
/*     AF.  No pivoting, so no permutations. */

#line 168 "sla_porpvgrw.f"
    if (lsame_("Upper", uplo, (ftnlen)5, (ftnlen)1)) {
#line 169 "sla_porpvgrw.f"
	i__1 = *ncols;
#line 169 "sla_porpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 170 "sla_porpvgrw.f"
	    i__2 = j;
#line 170 "sla_porpvgrw.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 171 "sla_porpvgrw.f"
		d__2 = (d__1 = af[i__ + j * af_dim1], abs(d__1)), d__3 = work[
			j];
#line 171 "sla_porpvgrw.f"
		work[j] = max(d__2,d__3);
#line 172 "sla_porpvgrw.f"
	    }
#line 173 "sla_porpvgrw.f"
	}
#line 174 "sla_porpvgrw.f"
    } else {
#line 175 "sla_porpvgrw.f"
	i__1 = *ncols;
#line 175 "sla_porpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 176 "sla_porpvgrw.f"
	    i__2 = *ncols;
#line 176 "sla_porpvgrw.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 177 "sla_porpvgrw.f"
		d__2 = (d__1 = af[i__ + j * af_dim1], abs(d__1)), d__3 = work[
			j];
#line 177 "sla_porpvgrw.f"
		work[j] = max(d__2,d__3);
#line 178 "sla_porpvgrw.f"
	    }
#line 179 "sla_porpvgrw.f"
	}
#line 180 "sla_porpvgrw.f"
    }

/*     Compute the *inverse* of the max element growth factor.  Dividing */
/*     by zero would imply the largest entry of the factor's column is */
/*     zero.  Than can happen when either the column of A is zero or */
/*     massive pivots made the factor underflow to zero.  Neither counts */
/*     as growth in itself, so simply ignore terms with zero */
/*     denominators. */

#line 189 "sla_porpvgrw.f"
    if (lsame_("Upper", uplo, (ftnlen)5, (ftnlen)1)) {
#line 190 "sla_porpvgrw.f"
	i__1 = *ncols;
#line 190 "sla_porpvgrw.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 191 "sla_porpvgrw.f"
	    umax = work[i__];
#line 192 "sla_porpvgrw.f"
	    amax = work[*ncols + i__];
#line 193 "sla_porpvgrw.f"
	    if (umax != 0.) {
/* Computing MIN */
#line 194 "sla_porpvgrw.f"
		d__1 = amax / umax;
#line 194 "sla_porpvgrw.f"
		rpvgrw = min(d__1,rpvgrw);
#line 195 "sla_porpvgrw.f"
	    }
#line 196 "sla_porpvgrw.f"
	}
#line 197 "sla_porpvgrw.f"
    } else {
#line 198 "sla_porpvgrw.f"
	i__1 = *ncols;
#line 198 "sla_porpvgrw.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 199 "sla_porpvgrw.f"
	    umax = work[i__];
#line 200 "sla_porpvgrw.f"
	    amax = work[*ncols + i__];
#line 201 "sla_porpvgrw.f"
	    if (umax != 0.) {
/* Computing MIN */
#line 202 "sla_porpvgrw.f"
		d__1 = amax / umax;
#line 202 "sla_porpvgrw.f"
		rpvgrw = min(d__1,rpvgrw);
#line 203 "sla_porpvgrw.f"
	    }
#line 204 "sla_porpvgrw.f"
	}
#line 205 "sla_porpvgrw.f"
    }
#line 207 "sla_porpvgrw.f"
    ret_val = rpvgrw;
#line 208 "sla_porpvgrw.f"
    return ret_val;
} /* sla_porpvgrw__ */

