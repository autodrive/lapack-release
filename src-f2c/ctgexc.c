#line 1 "ctgexc.f"
/* ctgexc.f -- translated by f2c (version 20100827).
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

#line 1 "ctgexc.f"
/* > \brief \b CTGEXC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTGEXC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctgexc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctgexc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctgexc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, */
/*                          LDZ, IFST, ILST, INFO ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            WANTQ, WANTZ */
/*       INTEGER            IFST, ILST, INFO, LDA, LDB, LDQ, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTGEXC reorders the generalized Schur decomposition of a complex */
/* > matrix pair (A,B), using an unitary equivalence transformation */
/* > (A, B) := Q * (A, B) * Z**H, so that the diagonal block of (A, B) with */
/* > row index IFST is moved to row ILST. */
/* > */
/* > (A, B) must be in generalized Schur canonical form, that is, A and */
/* > B are both upper triangular. */
/* > */
/* > Optionally, the matrices Q and Z of generalized Schur vectors are */
/* > updated. */
/* > */
/* >        Q(in) * A(in) * Z(in)**H = Q(out) * A(out) * Z(out)**H */
/* >        Q(in) * B(in) * Z(in)**H = Q(out) * B(out) * Z(out)**H */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] WANTQ */
/* > \verbatim */
/* >          WANTQ is LOGICAL */
/* >          .TRUE. : update the left transformation matrix Q; */
/* >          .FALSE.: do not update Q. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* >          WANTZ is LOGICAL */
/* >          .TRUE. : update the right transformation matrix Z; */
/* >          .FALSE.: do not update Z. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices A and B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the upper triangular matrix A in the pair (A, B). */
/* >          On exit, the updated matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,N) */
/* >          On entry, the upper triangular matrix B in the pair (A, B). */
/* >          On exit, the updated matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is COMPLEX array, dimension (LDQ,N) */
/* >          On entry, if WANTQ = .TRUE., the unitary matrix Q. */
/* >          On exit, the updated matrix Q. */
/* >          If WANTQ = .FALSE., Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q. LDQ >= 1; */
/* >          If WANTQ = .TRUE., LDQ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is COMPLEX array, dimension (LDZ,N) */
/* >          On entry, if WANTZ = .TRUE., the unitary matrix Z. */
/* >          On exit, the updated matrix Z. */
/* >          If WANTZ = .FALSE., Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z. LDZ >= 1; */
/* >          If WANTZ = .TRUE., LDZ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] IFST */
/* > \verbatim */
/* >          IFST is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in,out] ILST */
/* > \verbatim */
/* >          ILST is INTEGER */
/* >          Specify the reordering of the diagonal blocks of (A, B). */
/* >          The block with row index IFST is moved to row ILST, by a */
/* >          sequence of swapping between adjacent blocks. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >           =0:  Successful exit. */
/* >           <0:  if INFO = -i, the i-th argument had an illegal value. */
/* >           =1:  The transformed matrix pair (A, B) would be too far */
/* >                from generalized Schur form; the problem is ill- */
/* >                conditioned. (A, B) may have been partially reordered, */
/* >                and ILST points to the first row of the current */
/* >                position of the block being moved. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2017 */

/* > \ingroup complexGEcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* >     Umea University, S-901 87 Umea, Sweden. */

/* > \par References: */
/*  ================ */
/* > */
/* >  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the */
/* >      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in */
/* >      M.S. Moonen et al (eds), Linear Algebra for Large Scale and */
/* >      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218. */
/* > \n */
/* >  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified */
/* >      Eigenvalues of a Regular Matrix Pair (A, B) and Condition */
/* >      Estimation: Theory, Algorithms and Software, Report */
/* >      UMINF - 94.04, Department of Computing Science, Umea University, */
/* >      S-901 87 Umea, Sweden, 1994. Also as LAPACK Working Note 87. */
/* >      To appear in Numerical Algorithms, 1996. */
/* > \n */
/* >  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software */
/* >      for Solving the Generalized Sylvester Equation and Estimating the */
/* >      Separation between Regular Matrix Pairs, Report UMINF - 93.23, */
/* >      Department of Computing Science, Umea University, S-901 87 Umea, */
/* >      Sweden, December 1993, Revised April 1994, Also as LAPACK working */
/* >      Note 75. To appear in ACM Trans. on Math. Software, Vol 22, No 1, */
/* >      1996. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ctgexc_(logical *wantq, logical *wantz, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *q, integer *ldq, doublecomplex *z__, integer *ldz, 
	integer *ifst, integer *ilst, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1;

    /* Local variables */
    static integer here;
    extern /* Subroutine */ int ctgex2_(logical *, logical *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     integer *), xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and test input arguments. */
#line 231 "ctgexc.f"
    /* Parameter adjustments */
#line 231 "ctgexc.f"
    a_dim1 = *lda;
#line 231 "ctgexc.f"
    a_offset = 1 + a_dim1;
#line 231 "ctgexc.f"
    a -= a_offset;
#line 231 "ctgexc.f"
    b_dim1 = *ldb;
#line 231 "ctgexc.f"
    b_offset = 1 + b_dim1;
#line 231 "ctgexc.f"
    b -= b_offset;
#line 231 "ctgexc.f"
    q_dim1 = *ldq;
#line 231 "ctgexc.f"
    q_offset = 1 + q_dim1;
#line 231 "ctgexc.f"
    q -= q_offset;
#line 231 "ctgexc.f"
    z_dim1 = *ldz;
#line 231 "ctgexc.f"
    z_offset = 1 + z_dim1;
#line 231 "ctgexc.f"
    z__ -= z_offset;
#line 231 "ctgexc.f"

#line 231 "ctgexc.f"
    /* Function Body */
#line 231 "ctgexc.f"
    *info = 0;
#line 232 "ctgexc.f"
    if (*n < 0) {
#line 233 "ctgexc.f"
	*info = -3;
#line 234 "ctgexc.f"
    } else if (*lda < max(1,*n)) {
#line 235 "ctgexc.f"
	*info = -5;
#line 236 "ctgexc.f"
    } else if (*ldb < max(1,*n)) {
#line 237 "ctgexc.f"
	*info = -7;
#line 238 "ctgexc.f"
    } else if (*ldq < 1 || *wantq && *ldq < max(1,*n)) {
#line 239 "ctgexc.f"
	*info = -9;
#line 240 "ctgexc.f"
    } else if (*ldz < 1 || *wantz && *ldz < max(1,*n)) {
#line 241 "ctgexc.f"
	*info = -11;
#line 242 "ctgexc.f"
    } else if (*ifst < 1 || *ifst > *n) {
#line 243 "ctgexc.f"
	*info = -12;
#line 244 "ctgexc.f"
    } else if (*ilst < 1 || *ilst > *n) {
#line 245 "ctgexc.f"
	*info = -13;
#line 246 "ctgexc.f"
    }
#line 247 "ctgexc.f"
    if (*info != 0) {
#line 248 "ctgexc.f"
	i__1 = -(*info);
#line 248 "ctgexc.f"
	xerbla_("CTGEXC", &i__1, (ftnlen)6);
#line 249 "ctgexc.f"
	return 0;
#line 250 "ctgexc.f"
    }

/*     Quick return if possible */

#line 254 "ctgexc.f"
    if (*n <= 1) {
#line 254 "ctgexc.f"
	return 0;
#line 254 "ctgexc.f"
    }
#line 256 "ctgexc.f"
    if (*ifst == *ilst) {
#line 256 "ctgexc.f"
	return 0;
#line 256 "ctgexc.f"
    }

#line 259 "ctgexc.f"
    if (*ifst < *ilst) {

#line 261 "ctgexc.f"
	here = *ifst;

#line 263 "ctgexc.f"
L10:

/*        Swap with next one below */

#line 267 "ctgexc.f"
	ctgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[
		q_offset], ldq, &z__[z_offset], ldz, &here, info);
#line 269 "ctgexc.f"
	if (*info != 0) {
#line 270 "ctgexc.f"
	    *ilst = here;
#line 271 "ctgexc.f"
	    return 0;
#line 272 "ctgexc.f"
	}
#line 273 "ctgexc.f"
	++here;
#line 274 "ctgexc.f"
	if (here < *ilst) {
#line 274 "ctgexc.f"
	    goto L10;
#line 274 "ctgexc.f"
	}
#line 276 "ctgexc.f"
	--here;
#line 277 "ctgexc.f"
    } else {
#line 278 "ctgexc.f"
	here = *ifst - 1;

#line 280 "ctgexc.f"
L20:

/*        Swap with next one above */

#line 284 "ctgexc.f"
	ctgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[
		q_offset], ldq, &z__[z_offset], ldz, &here, info);
#line 286 "ctgexc.f"
	if (*info != 0) {
#line 287 "ctgexc.f"
	    *ilst = here;
#line 288 "ctgexc.f"
	    return 0;
#line 289 "ctgexc.f"
	}
#line 290 "ctgexc.f"
	--here;
#line 291 "ctgexc.f"
	if (here >= *ilst) {
#line 291 "ctgexc.f"
	    goto L20;
#line 291 "ctgexc.f"
	}
#line 293 "ctgexc.f"
	++here;
#line 294 "ctgexc.f"
    }
#line 295 "ctgexc.f"
    *ilst = here;
#line 296 "ctgexc.f"
    return 0;

/*     End of CTGEXC */

} /* ctgexc_ */

