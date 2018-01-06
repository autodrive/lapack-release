#line 1 "cla_porpvgrw.f"
/* cla_porpvgrw.f -- translated by f2c (version 20100827).
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

#line 1 "cla_porpvgrw.f"
/* > \brief \b CLA_PORPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a symmetric or Her
mitian positive-definite matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLA_PORPVGRW + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_por
pvgrw.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_por
pvgrw.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_por
pvgrw.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION CLA_PORPVGRW( UPLO, NCOLS, A, LDA, AF, LDAF, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER*1        UPLO */
/*       INTEGER            NCOLS, LDA, LDAF */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), AF( LDAF, * ) */
/*       REAL               WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > */
/* > CLA_PORPVGRW computes the reciprocal pivot growth factor */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          AF is COMPLEX array, dimension (LDAF,N) */
/* >     The triangular factor U or L from the Cholesky factorization */
/* >     A = U**T*U or A = L*L**T, as computed by CPOTRF. */
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

/* > \date June 2016 */

/* > \ingroup complexPOcomputational */

/*  ===================================================================== */
doublereal cla_porpvgrw__(char *uplo, integer *ncols, doublecomplex *a, 
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
#line 144 "cla_porpvgrw.f"
    /* Parameter adjustments */
#line 144 "cla_porpvgrw.f"
    a_dim1 = *lda;
#line 144 "cla_porpvgrw.f"
    a_offset = 1 + a_dim1;
#line 144 "cla_porpvgrw.f"
    a -= a_offset;
#line 144 "cla_porpvgrw.f"
    af_dim1 = *ldaf;
#line 144 "cla_porpvgrw.f"
    af_offset = 1 + af_dim1;
#line 144 "cla_porpvgrw.f"
    af -= af_offset;
#line 144 "cla_porpvgrw.f"
    --work;
#line 144 "cla_porpvgrw.f"

#line 144 "cla_porpvgrw.f"
    /* Function Body */
#line 144 "cla_porpvgrw.f"
    upper = lsame_("Upper", uplo, (ftnlen)5, (ftnlen)1);

/*     SPOTRF will have factored only the NCOLSxNCOLS leading minor, so */
/*     we restrict the growth search to that minor and use only the first */
/*     2*NCOLS workspace entries. */

#line 150 "cla_porpvgrw.f"
    rpvgrw = 1.;
#line 151 "cla_porpvgrw.f"
    i__1 = *ncols << 1;
#line 151 "cla_porpvgrw.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 152 "cla_porpvgrw.f"
	work[i__] = 0.;
#line 153 "cla_porpvgrw.f"
    }

/*     Find the max magnitude entry of each column. */

#line 157 "cla_porpvgrw.f"
    if (upper) {
#line 158 "cla_porpvgrw.f"
	i__1 = *ncols;
#line 158 "cla_porpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 159 "cla_porpvgrw.f"
	    i__2 = j;
#line 159 "cla_porpvgrw.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 160 "cla_porpvgrw.f"
		i__3 = i__ + j * a_dim1;
#line 160 "cla_porpvgrw.f"
		d__3 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			+ j * a_dim1]), abs(d__2)), d__4 = work[*ncols + j];
#line 160 "cla_porpvgrw.f"
		work[*ncols + j] = max(d__3,d__4);
#line 162 "cla_porpvgrw.f"
	    }
#line 163 "cla_porpvgrw.f"
	}
#line 164 "cla_porpvgrw.f"
    } else {
#line 165 "cla_porpvgrw.f"
	i__1 = *ncols;
#line 165 "cla_porpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 166 "cla_porpvgrw.f"
	    i__2 = *ncols;
#line 166 "cla_porpvgrw.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 167 "cla_porpvgrw.f"
		i__3 = i__ + j * a_dim1;
#line 167 "cla_porpvgrw.f"
		d__3 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			+ j * a_dim1]), abs(d__2)), d__4 = work[*ncols + j];
#line 167 "cla_porpvgrw.f"
		work[*ncols + j] = max(d__3,d__4);
#line 169 "cla_porpvgrw.f"
	    }
#line 170 "cla_porpvgrw.f"
	}
#line 171 "cla_porpvgrw.f"
    }

/*     Now find the max magnitude entry of each column of the factor in */
/*     AF.  No pivoting, so no permutations. */

#line 176 "cla_porpvgrw.f"
    if (lsame_("Upper", uplo, (ftnlen)5, (ftnlen)1)) {
#line 177 "cla_porpvgrw.f"
	i__1 = *ncols;
#line 177 "cla_porpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 178 "cla_porpvgrw.f"
	    i__2 = j;
#line 178 "cla_porpvgrw.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 179 "cla_porpvgrw.f"
		i__3 = i__ + j * af_dim1;
#line 179 "cla_porpvgrw.f"
		d__3 = (d__1 = af[i__3].r, abs(d__1)) + (d__2 = d_imag(&af[
			i__ + j * af_dim1]), abs(d__2)), d__4 = work[j];
#line 179 "cla_porpvgrw.f"
		work[j] = max(d__3,d__4);
#line 180 "cla_porpvgrw.f"
	    }
#line 181 "cla_porpvgrw.f"
	}
#line 182 "cla_porpvgrw.f"
    } else {
#line 183 "cla_porpvgrw.f"
	i__1 = *ncols;
#line 183 "cla_porpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 184 "cla_porpvgrw.f"
	    i__2 = *ncols;
#line 184 "cla_porpvgrw.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 185 "cla_porpvgrw.f"
		i__3 = i__ + j * af_dim1;
#line 185 "cla_porpvgrw.f"
		d__3 = (d__1 = af[i__3].r, abs(d__1)) + (d__2 = d_imag(&af[
			i__ + j * af_dim1]), abs(d__2)), d__4 = work[j];
#line 185 "cla_porpvgrw.f"
		work[j] = max(d__3,d__4);
#line 186 "cla_porpvgrw.f"
	    }
#line 187 "cla_porpvgrw.f"
	}
#line 188 "cla_porpvgrw.f"
    }

/*     Compute the *inverse* of the max element growth factor.  Dividing */
/*     by zero would imply the largest entry of the factor's column is */
/*     zero.  Than can happen when either the column of A is zero or */
/*     massive pivots made the factor underflow to zero.  Neither counts */
/*     as growth in itself, so simply ignore terms with zero */
/*     denominators. */

#line 197 "cla_porpvgrw.f"
    if (lsame_("Upper", uplo, (ftnlen)5, (ftnlen)1)) {
#line 198 "cla_porpvgrw.f"
	i__1 = *ncols;
#line 198 "cla_porpvgrw.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 199 "cla_porpvgrw.f"
	    umax = work[i__];
#line 200 "cla_porpvgrw.f"
	    amax = work[*ncols + i__];
#line 201 "cla_porpvgrw.f"
	    if (umax != 0.) {
/* Computing MIN */
#line 202 "cla_porpvgrw.f"
		d__1 = amax / umax;
#line 202 "cla_porpvgrw.f"
		rpvgrw = min(d__1,rpvgrw);
#line 203 "cla_porpvgrw.f"
	    }
#line 204 "cla_porpvgrw.f"
	}
#line 205 "cla_porpvgrw.f"
    } else {
#line 206 "cla_porpvgrw.f"
	i__1 = *ncols;
#line 206 "cla_porpvgrw.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 207 "cla_porpvgrw.f"
	    umax = work[i__];
#line 208 "cla_porpvgrw.f"
	    amax = work[*ncols + i__];
#line 209 "cla_porpvgrw.f"
	    if (umax != 0.) {
/* Computing MIN */
#line 210 "cla_porpvgrw.f"
		d__1 = amax / umax;
#line 210 "cla_porpvgrw.f"
		rpvgrw = min(d__1,rpvgrw);
#line 211 "cla_porpvgrw.f"
	    }
#line 212 "cla_porpvgrw.f"
	}
#line 213 "cla_porpvgrw.f"
    }
#line 215 "cla_porpvgrw.f"
    ret_val = rpvgrw;
#line 216 "cla_porpvgrw.f"
    return ret_val;
} /* cla_porpvgrw__ */

