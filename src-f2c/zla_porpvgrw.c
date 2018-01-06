#line 1 "zla_porpvgrw.f"
/* zla_porpvgrw.f -- translated by f2c (version 20100827).
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

#line 1 "zla_porpvgrw.f"
/* > \brief \b ZLA_PORPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a symmetric or Her
mitian positive-definite matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLA_PORPVGRW + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_por
pvgrw.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_por
pvgrw.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_por
pvgrw.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLA_PORPVGRW( UPLO, NCOLS, A, LDA, AF, */
/*                                               LDAF, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER*1        UPLO */
/*       INTEGER            NCOLS, LDA, LDAF */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), AF( LDAF, * ) */
/*       DOUBLE PRECISION   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > */
/* > ZLA_PORPVGRW computes the reciprocal pivot growth factor */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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
/* >          AF is COMPLEX*16 array, dimension (LDAF,N) */
/* >     The triangular factor U or L from the Cholesky factorization */
/* >     A = U**T*U or A = L*L**T, as computed by ZPOTRF. */
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
/* >          WORK is DOUBLE PRECISION array, dimension (2*N) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup complex16POcomputational */

/*  ===================================================================== */
doublereal zla_porpvgrw__(char *uplo, integer *ncols, doublecomplex *a, 
	integer *lda, doublecomplex *af, integer *ldaf, doublereal *work, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2, i__3;
    doublereal ret_val, d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublereal amax, umax;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    static doublereal rpvgrw;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function Definitions .. */
/*     .. */
/*     .. Executable Statements .. */
#line 146 "zla_porpvgrw.f"
    /* Parameter adjustments */
#line 146 "zla_porpvgrw.f"
    a_dim1 = *lda;
#line 146 "zla_porpvgrw.f"
    a_offset = 1 + a_dim1;
#line 146 "zla_porpvgrw.f"
    a -= a_offset;
#line 146 "zla_porpvgrw.f"
    af_dim1 = *ldaf;
#line 146 "zla_porpvgrw.f"
    af_offset = 1 + af_dim1;
#line 146 "zla_porpvgrw.f"
    af -= af_offset;
#line 146 "zla_porpvgrw.f"
    --work;
#line 146 "zla_porpvgrw.f"

#line 146 "zla_porpvgrw.f"
    /* Function Body */
#line 146 "zla_porpvgrw.f"
    upper = lsame_("Upper", uplo, (ftnlen)5, (ftnlen)1);

/*     DPOTRF will have factored only the NCOLSxNCOLS leading minor, so */
/*     we restrict the growth search to that minor and use only the first */
/*     2*NCOLS workspace entries. */

#line 152 "zla_porpvgrw.f"
    rpvgrw = 1.;
#line 153 "zla_porpvgrw.f"
    i__1 = *ncols << 1;
#line 153 "zla_porpvgrw.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 154 "zla_porpvgrw.f"
	work[i__] = 0.;
#line 155 "zla_porpvgrw.f"
    }

/*     Find the max magnitude entry of each column. */

#line 159 "zla_porpvgrw.f"
    if (upper) {
#line 160 "zla_porpvgrw.f"
	i__1 = *ncols;
#line 160 "zla_porpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 161 "zla_porpvgrw.f"
	    i__2 = j;
#line 161 "zla_porpvgrw.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 162 "zla_porpvgrw.f"
		i__3 = i__ + j * a_dim1;
#line 162 "zla_porpvgrw.f"
		d__3 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			+ j * a_dim1]), abs(d__2)), d__4 = work[*ncols + j];
#line 162 "zla_porpvgrw.f"
		work[*ncols + j] = max(d__3,d__4);
#line 164 "zla_porpvgrw.f"
	    }
#line 165 "zla_porpvgrw.f"
	}
#line 166 "zla_porpvgrw.f"
    } else {
#line 167 "zla_porpvgrw.f"
	i__1 = *ncols;
#line 167 "zla_porpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 168 "zla_porpvgrw.f"
	    i__2 = *ncols;
#line 168 "zla_porpvgrw.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 169 "zla_porpvgrw.f"
		i__3 = i__ + j * a_dim1;
#line 169 "zla_porpvgrw.f"
		d__3 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			+ j * a_dim1]), abs(d__2)), d__4 = work[*ncols + j];
#line 169 "zla_porpvgrw.f"
		work[*ncols + j] = max(d__3,d__4);
#line 171 "zla_porpvgrw.f"
	    }
#line 172 "zla_porpvgrw.f"
	}
#line 173 "zla_porpvgrw.f"
    }

/*     Now find the max magnitude entry of each column of the factor in */
/*     AF.  No pivoting, so no permutations. */

#line 178 "zla_porpvgrw.f"
    if (lsame_("Upper", uplo, (ftnlen)5, (ftnlen)1)) {
#line 179 "zla_porpvgrw.f"
	i__1 = *ncols;
#line 179 "zla_porpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 180 "zla_porpvgrw.f"
	    i__2 = j;
#line 180 "zla_porpvgrw.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 181 "zla_porpvgrw.f"
		i__3 = i__ + j * af_dim1;
#line 181 "zla_porpvgrw.f"
		d__3 = (d__1 = af[i__3].r, abs(d__1)) + (d__2 = d_imag(&af[
			i__ + j * af_dim1]), abs(d__2)), d__4 = work[j];
#line 181 "zla_porpvgrw.f"
		work[j] = max(d__3,d__4);
#line 182 "zla_porpvgrw.f"
	    }
#line 183 "zla_porpvgrw.f"
	}
#line 184 "zla_porpvgrw.f"
    } else {
#line 185 "zla_porpvgrw.f"
	i__1 = *ncols;
#line 185 "zla_porpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 186 "zla_porpvgrw.f"
	    i__2 = *ncols;
#line 186 "zla_porpvgrw.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 187 "zla_porpvgrw.f"
		i__3 = i__ + j * af_dim1;
#line 187 "zla_porpvgrw.f"
		d__3 = (d__1 = af[i__3].r, abs(d__1)) + (d__2 = d_imag(&af[
			i__ + j * af_dim1]), abs(d__2)), d__4 = work[j];
#line 187 "zla_porpvgrw.f"
		work[j] = max(d__3,d__4);
#line 188 "zla_porpvgrw.f"
	    }
#line 189 "zla_porpvgrw.f"
	}
#line 190 "zla_porpvgrw.f"
    }

/*     Compute the *inverse* of the max element growth factor.  Dividing */
/*     by zero would imply the largest entry of the factor's column is */
/*     zero.  Than can happen when either the column of A is zero or */
/*     massive pivots made the factor underflow to zero.  Neither counts */
/*     as growth in itself, so simply ignore terms with zero */
/*     denominators. */

#line 199 "zla_porpvgrw.f"
    if (lsame_("Upper", uplo, (ftnlen)5, (ftnlen)1)) {
#line 200 "zla_porpvgrw.f"
	i__1 = *ncols;
#line 200 "zla_porpvgrw.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 201 "zla_porpvgrw.f"
	    umax = work[i__];
#line 202 "zla_porpvgrw.f"
	    amax = work[*ncols + i__];
#line 203 "zla_porpvgrw.f"
	    if (umax != 0.) {
/* Computing MIN */
#line 204 "zla_porpvgrw.f"
		d__1 = amax / umax;
#line 204 "zla_porpvgrw.f"
		rpvgrw = min(d__1,rpvgrw);
#line 205 "zla_porpvgrw.f"
	    }
#line 206 "zla_porpvgrw.f"
	}
#line 207 "zla_porpvgrw.f"
    } else {
#line 208 "zla_porpvgrw.f"
	i__1 = *ncols;
#line 208 "zla_porpvgrw.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 209 "zla_porpvgrw.f"
	    umax = work[i__];
#line 210 "zla_porpvgrw.f"
	    amax = work[*ncols + i__];
#line 211 "zla_porpvgrw.f"
	    if (umax != 0.) {
/* Computing MIN */
#line 212 "zla_porpvgrw.f"
		d__1 = amax / umax;
#line 212 "zla_porpvgrw.f"
		rpvgrw = min(d__1,rpvgrw);
#line 213 "zla_porpvgrw.f"
	    }
#line 214 "zla_porpvgrw.f"
	}
#line 215 "zla_porpvgrw.f"
    }
#line 217 "zla_porpvgrw.f"
    ret_val = rpvgrw;
#line 218 "zla_porpvgrw.f"
    return ret_val;
} /* zla_porpvgrw__ */

