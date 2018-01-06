#line 1 "dlaebz.f"
/* dlaebz.f -- translated by f2c (version 20100827).
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

#line 1 "dlaebz.f"
/* > \brief \b DLAEBZ computes the number of eigenvalues of a real symmetric tridiagonal matrix which are less
 than or equal to a given value, and performs other tasks required by the routine sstebz. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAEBZ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaebz.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaebz.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaebz.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAEBZ( IJOB, NITMAX, N, MMAX, MINP, NBMIN, ABSTOL, */
/*                          RELTOL, PIVMIN, D, E, E2, NVAL, AB, C, MOUT, */
/*                          NAB, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IJOB, INFO, MINP, MMAX, MOUT, N, NBMIN, NITMAX */
/*       DOUBLE PRECISION   ABSTOL, PIVMIN, RELTOL */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ), NAB( MMAX, * ), NVAL( * ) */
/*       DOUBLE PRECISION   AB( MMAX, * ), C( * ), D( * ), E( * ), E2( * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAEBZ contains the iteration loops which compute and use the */
/* > function N(w), which is the count of eigenvalues of a symmetric */
/* > tridiagonal matrix T less than or equal to its argument  w.  It */
/* > performs a choice of two types of loops: */
/* > */
/* > IJOB=1, followed by */
/* > IJOB=2: It takes as input a list of intervals and returns a list of */
/* >         sufficiently small intervals whose union contains the same */
/* >         eigenvalues as the union of the original intervals. */
/* >         The input intervals are (AB(j,1),AB(j,2)], j=1,...,MINP. */
/* >         The output interval (AB(j,1),AB(j,2)] will contain */
/* >         eigenvalues NAB(j,1)+1,...,NAB(j,2), where 1 <= j <= MOUT. */
/* > */
/* > IJOB=3: It performs a binary search in each input interval */
/* >         (AB(j,1),AB(j,2)] for a point  w(j)  such that */
/* >         N(w(j))=NVAL(j), and uses  C(j)  as the starting point of */
/* >         the search.  If such a w(j) is found, then on output */
/* >         AB(j,1)=AB(j,2)=w.  If no such w(j) is found, then on output */
/* >         (AB(j,1),AB(j,2)] will be a small interval containing the */
/* >         point where N(w) jumps through NVAL(j), unless that point */
/* >         lies outside the initial interval. */
/* > */
/* > Note that the intervals are in all cases half-open intervals, */
/* > i.e., of the form  (a,b] , which includes  b  but not  a . */
/* > */
/* > To avoid underflow, the matrix should be scaled so that its largest */
/* > element is no greater than  overflow**(1/2) * underflow**(1/4) */
/* > in absolute value.  To assure the most accurate computation */
/* > of small eigenvalues, the matrix should be scaled to be */
/* > not much smaller than that, either. */
/* > */
/* > See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal */
/* > Matrix", Report CS41, Computer Science Dept., Stanford */
/* > University, July 21, 1966 */
/* > */
/* > Note: the arguments are, in general, *not* checked for unreasonable */
/* > values. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] IJOB */
/* > \verbatim */
/* >          IJOB is INTEGER */
/* >          Specifies what is to be done: */
/* >          = 1:  Compute NAB for the initial intervals. */
/* >          = 2:  Perform bisection iteration to find eigenvalues of T. */
/* >          = 3:  Perform bisection iteration to invert N(w), i.e., */
/* >                to find a point which has a specified number of */
/* >                eigenvalues of T to its left. */
/* >          Other values will cause DLAEBZ to return with INFO=-1. */
/* > \endverbatim */
/* > */
/* > \param[in] NITMAX */
/* > \verbatim */
/* >          NITMAX is INTEGER */
/* >          The maximum number of "levels" of bisection to be */
/* >          performed, i.e., an interval of width W will not be made */
/* >          smaller than 2^(-NITMAX) * W.  If not all intervals */
/* >          have converged after NITMAX iterations, then INFO is set */
/* >          to the number of non-converged intervals. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The dimension n of the tridiagonal matrix T.  It must be at */
/* >          least 1. */
/* > \endverbatim */
/* > */
/* > \param[in] MMAX */
/* > \verbatim */
/* >          MMAX is INTEGER */
/* >          The maximum number of intervals.  If more than MMAX intervals */
/* >          are generated, then DLAEBZ will quit with INFO=MMAX+1. */
/* > \endverbatim */
/* > */
/* > \param[in] MINP */
/* > \verbatim */
/* >          MINP is INTEGER */
/* >          The initial number of intervals.  It may not be greater than */
/* >          MMAX. */
/* > \endverbatim */
/* > */
/* > \param[in] NBMIN */
/* > \verbatim */
/* >          NBMIN is INTEGER */
/* >          The smallest number of intervals that should be processed */
/* >          using a vector loop.  If zero, then only the scalar loop */
/* >          will be used. */
/* > \endverbatim */
/* > */
/* > \param[in] ABSTOL */
/* > \verbatim */
/* >          ABSTOL is DOUBLE PRECISION */
/* >          The minimum (absolute) width of an interval.  When an */
/* >          interval is narrower than ABSTOL, or than RELTOL times the */
/* >          larger (in magnitude) endpoint, then it is considered to be */
/* >          sufficiently small, i.e., converged.  This must be at least */
/* >          zero. */
/* > \endverbatim */
/* > */
/* > \param[in] RELTOL */
/* > \verbatim */
/* >          RELTOL is DOUBLE PRECISION */
/* >          The minimum relative width of an interval.  When an interval */
/* >          is narrower than ABSTOL, or than RELTOL times the larger (in */
/* >          magnitude) endpoint, then it is considered to be */
/* >          sufficiently small, i.e., converged.  Note: this should */
/* >          always be at least radix*machine epsilon. */
/* > \endverbatim */
/* > */
/* > \param[in] PIVMIN */
/* > \verbatim */
/* >          PIVMIN is DOUBLE PRECISION */
/* >          The minimum absolute value of a "pivot" in the Sturm */
/* >          sequence loop. */
/* >          This must be at least  max |e(j)**2|*safe_min  and at */
/* >          least safe_min, where safe_min is at least */
/* >          the smallest number that can divide one without overflow. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N) */
/* >          The offdiagonal elements of the tridiagonal matrix T in */
/* >          positions 1 through N-1.  E(N) is arbitrary. */
/* > \endverbatim */
/* > */
/* > \param[in] E2 */
/* > \verbatim */
/* >          E2 is DOUBLE PRECISION array, dimension (N) */
/* >          The squares of the offdiagonal elements of the tridiagonal */
/* >          matrix T.  E2(N) is ignored. */
/* > \endverbatim */
/* > */
/* > \param[in,out] NVAL */
/* > \verbatim */
/* >          NVAL is INTEGER array, dimension (MINP) */
/* >          If IJOB=1 or 2, not referenced. */
/* >          If IJOB=3, the desired values of N(w).  The elements of NVAL */
/* >          will be reordered to correspond with the intervals in AB. */
/* >          Thus, NVAL(j) on output will not, in general be the same as */
/* >          NVAL(j) on input, but it will correspond with the interval */
/* >          (AB(j,1),AB(j,2)] on output. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is DOUBLE PRECISION array, dimension (MMAX,2) */
/* >          The endpoints of the intervals.  AB(j,1) is  a(j), the left */
/* >          endpoint of the j-th interval, and AB(j,2) is b(j), the */
/* >          right endpoint of the j-th interval.  The input intervals */
/* >          will, in general, be modified, split, and reordered by the */
/* >          calculation. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (MMAX) */
/* >          If IJOB=1, ignored. */
/* >          If IJOB=2, workspace. */
/* >          If IJOB=3, then on input C(j) should be initialized to the */
/* >          first search point in the binary search. */
/* > \endverbatim */
/* > */
/* > \param[out] MOUT */
/* > \verbatim */
/* >          MOUT is INTEGER */
/* >          If IJOB=1, the number of eigenvalues in the intervals. */
/* >          If IJOB=2 or 3, the number of intervals output. */
/* >          If IJOB=3, MOUT will equal MINP. */
/* > \endverbatim */
/* > */
/* > \param[in,out] NAB */
/* > \verbatim */
/* >          NAB is INTEGER array, dimension (MMAX,2) */
/* >          If IJOB=1, then on output NAB(i,j) will be set to N(AB(i,j)). */
/* >          If IJOB=2, then on input, NAB(i,j) should be set.  It must */
/* >             satisfy the condition: */
/* >             N(AB(i,1)) <= NAB(i,1) <= NAB(i,2) <= N(AB(i,2)), */
/* >             which means that in interval i only eigenvalues */
/* >             NAB(i,1)+1,...,NAB(i,2) will be considered.  Usually, */
/* >             NAB(i,j)=N(AB(i,j)), from a previous call to DLAEBZ with */
/* >             IJOB=1. */
/* >             On output, NAB(i,j) will contain */
/* >             max(na(k),min(nb(k),N(AB(i,j)))), where k is the index of */
/* >             the input interval that the output interval */
/* >             (AB(j,1),AB(j,2)] came from, and na(k) and nb(k) are the */
/* >             the input values of NAB(k,1) and NAB(k,2). */
/* >          If IJOB=3, then on output, NAB(i,j) contains N(AB(i,j)), */
/* >             unless N(w) > NVAL(i) for all search points  w , in which */
/* >             case NAB(i,1) will not be modified, i.e., the output */
/* >             value will be the same as the input value (modulo */
/* >             reorderings -- see NVAL and AB), or unless N(w) < NVAL(i) */
/* >             for all search points  w , in which case NAB(i,2) will */
/* >             not be modified.  Normally, NAB should be set to some */
/* >             distinctive value(s) before DLAEBZ is called. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MMAX) */
/* >          Workspace. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (MMAX) */
/* >          Workspace. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:       All intervals converged. */
/* >          = 1--MMAX: The last INFO intervals did not converge. */
/* >          = MMAX+1:  More than MMAX intervals were generated. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup OTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >      This routine is intended to be called only by other LAPACK */
/* >  routines, thus the interface is less user-friendly.  It is intended */
/* >  for two purposes: */
/* > */
/* >  (a) finding eigenvalues.  In this case, DLAEBZ should have one or */
/* >      more initial intervals set up in AB, and DLAEBZ should be called */
/* >      with IJOB=1.  This sets up NAB, and also counts the eigenvalues. */
/* >      Intervals with no eigenvalues would usually be thrown out at */
/* >      this point.  Also, if not all the eigenvalues in an interval i */
/* >      are desired, NAB(i,1) can be increased or NAB(i,2) decreased. */
/* >      For example, set NAB(i,1)=NAB(i,2)-1 to get the largest */
/* >      eigenvalue.  DLAEBZ is then called with IJOB=2 and MMAX */
/* >      no smaller than the value of MOUT returned by the call with */
/* >      IJOB=1.  After this (IJOB=2) call, eigenvalues NAB(i,1)+1 */
/* >      through NAB(i,2) are approximately AB(i,1) (or AB(i,2)) to the */
/* >      tolerance specified by ABSTOL and RELTOL. */
/* > */
/* >  (b) finding an interval (a',b'] containing eigenvalues w(f),...,w(l). */
/* >      In this case, start with a Gershgorin interval  (a,b).  Set up */
/* >      AB to contain 2 search intervals, both initially (a,b).  One */
/* >      NVAL element should contain  f-1  and the other should contain  l */
/* >      , while C should contain a and b, resp.  NAB(i,1) should be -1 */
/* >      and NAB(i,2) should be N+1, to flag an error if the desired */
/* >      interval does not lie in (a,b).  DLAEBZ is then called with */
/* >      IJOB=3.  On exit, if w(f-1) < w(f), then one of the intervals -- */
/* >      j -- will have AB(j,1)=AB(j,2) and NAB(j,1)=NAB(j,2)=f-1, while */
/* >      if, to the specified tolerance, w(f-k)=...=w(f+r), k > 0 and r */
/* >      >= 0, then the interval will have  N(AB(j,1))=NAB(j,1)=f-k and */
/* >      N(AB(j,2))=NAB(j,2)=f+r.  The cases w(l) < w(l+1) and */
/* >      w(l-r)=...=w(l+k) are handled similarly. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlaebz_(integer *ijob, integer *nitmax, integer *n, 
	integer *mmax, integer *minp, integer *nbmin, doublereal *abstol, 
	doublereal *reltol, doublereal *pivmin, doublereal *d__, doublereal *
	e, doublereal *e2, integer *nval, doublereal *ab, doublereal *c__, 
	integer *mout, integer *nab, doublereal *work, integer *iwork, 
	integer *info)
{
    /* System generated locals */
    integer nab_dim1, nab_offset, ab_dim1, ab_offset, i__1, i__2, i__3, i__4, 
	    i__5, i__6;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer j, kf, ji, kl, jp, jit;
    static doublereal tmp1, tmp2;
    static integer itmp1, itmp2, kfnew, klnew;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Check for Errors */

#line 356 "dlaebz.f"
    /* Parameter adjustments */
#line 356 "dlaebz.f"
    nab_dim1 = *mmax;
#line 356 "dlaebz.f"
    nab_offset = 1 + nab_dim1;
#line 356 "dlaebz.f"
    nab -= nab_offset;
#line 356 "dlaebz.f"
    ab_dim1 = *mmax;
#line 356 "dlaebz.f"
    ab_offset = 1 + ab_dim1;
#line 356 "dlaebz.f"
    ab -= ab_offset;
#line 356 "dlaebz.f"
    --d__;
#line 356 "dlaebz.f"
    --e;
#line 356 "dlaebz.f"
    --e2;
#line 356 "dlaebz.f"
    --nval;
#line 356 "dlaebz.f"
    --c__;
#line 356 "dlaebz.f"
    --work;
#line 356 "dlaebz.f"
    --iwork;
#line 356 "dlaebz.f"

#line 356 "dlaebz.f"
    /* Function Body */
#line 356 "dlaebz.f"
    *info = 0;
#line 357 "dlaebz.f"
    if (*ijob < 1 || *ijob > 3) {
#line 358 "dlaebz.f"
	*info = -1;
#line 359 "dlaebz.f"
	return 0;
#line 360 "dlaebz.f"
    }

/*     Initialize NAB */

#line 364 "dlaebz.f"
    if (*ijob == 1) {

/*        Compute the number of eigenvalues in the initial intervals. */

#line 368 "dlaebz.f"
	*mout = 0;
#line 369 "dlaebz.f"
	i__1 = *minp;
#line 369 "dlaebz.f"
	for (ji = 1; ji <= i__1; ++ji) {
#line 370 "dlaebz.f"
	    for (jp = 1; jp <= 2; ++jp) {
#line 371 "dlaebz.f"
		tmp1 = d__[1] - ab[ji + jp * ab_dim1];
#line 372 "dlaebz.f"
		if (abs(tmp1) < *pivmin) {
#line 372 "dlaebz.f"
		    tmp1 = -(*pivmin);
#line 372 "dlaebz.f"
		}
#line 374 "dlaebz.f"
		nab[ji + jp * nab_dim1] = 0;
#line 375 "dlaebz.f"
		if (tmp1 <= 0.) {
#line 375 "dlaebz.f"
		    nab[ji + jp * nab_dim1] = 1;
#line 375 "dlaebz.f"
		}

#line 378 "dlaebz.f"
		i__2 = *n;
#line 378 "dlaebz.f"
		for (j = 2; j <= i__2; ++j) {
#line 379 "dlaebz.f"
		    tmp1 = d__[j] - e2[j - 1] / tmp1 - ab[ji + jp * ab_dim1];
#line 380 "dlaebz.f"
		    if (abs(tmp1) < *pivmin) {
#line 380 "dlaebz.f"
			tmp1 = -(*pivmin);
#line 380 "dlaebz.f"
		    }
#line 382 "dlaebz.f"
		    if (tmp1 <= 0.) {
#line 382 "dlaebz.f"
			++nab[ji + jp * nab_dim1];
#line 382 "dlaebz.f"
		    }
#line 384 "dlaebz.f"
/* L10: */
#line 384 "dlaebz.f"
		}
#line 385 "dlaebz.f"
/* L20: */
#line 385 "dlaebz.f"
	    }
#line 386 "dlaebz.f"
	    *mout = *mout + nab[ji + (nab_dim1 << 1)] - nab[ji + nab_dim1];
#line 387 "dlaebz.f"
/* L30: */
#line 387 "dlaebz.f"
	}
#line 388 "dlaebz.f"
	return 0;
#line 389 "dlaebz.f"
    }

/*     Initialize for loop */

/*     KF and KL have the following meaning: */
/*        Intervals 1,...,KF-1 have converged. */
/*        Intervals KF,...,KL  still need to be refined. */

#line 397 "dlaebz.f"
    kf = 1;
#line 398 "dlaebz.f"
    kl = *minp;

/*     If IJOB=2, initialize C. */
/*     If IJOB=3, use the user-supplied starting point. */

#line 403 "dlaebz.f"
    if (*ijob == 2) {
#line 404 "dlaebz.f"
	i__1 = *minp;
#line 404 "dlaebz.f"
	for (ji = 1; ji <= i__1; ++ji) {
#line 405 "dlaebz.f"
	    c__[ji] = (ab[ji + ab_dim1] + ab[ji + (ab_dim1 << 1)]) * .5;
#line 406 "dlaebz.f"
/* L40: */
#line 406 "dlaebz.f"
	}
#line 407 "dlaebz.f"
    }

/*     Iteration loop */

#line 411 "dlaebz.f"
    i__1 = *nitmax;
#line 411 "dlaebz.f"
    for (jit = 1; jit <= i__1; ++jit) {

/*        Loop over intervals */

#line 415 "dlaebz.f"
	if (kl - kf + 1 >= *nbmin && *nbmin > 0) {

/*           Begin of Parallel Version of the loop */

#line 419 "dlaebz.f"
	    i__2 = kl;
#line 419 "dlaebz.f"
	    for (ji = kf; ji <= i__2; ++ji) {

/*              Compute N(c), the number of eigenvalues less than c */

#line 423 "dlaebz.f"
		work[ji] = d__[1] - c__[ji];
#line 424 "dlaebz.f"
		iwork[ji] = 0;
#line 425 "dlaebz.f"
		if (work[ji] <= *pivmin) {
#line 426 "dlaebz.f"
		    iwork[ji] = 1;
/* Computing MIN */
#line 427 "dlaebz.f"
		    d__1 = work[ji], d__2 = -(*pivmin);
#line 427 "dlaebz.f"
		    work[ji] = min(d__1,d__2);
#line 428 "dlaebz.f"
		}

#line 430 "dlaebz.f"
		i__3 = *n;
#line 430 "dlaebz.f"
		for (j = 2; j <= i__3; ++j) {
#line 431 "dlaebz.f"
		    work[ji] = d__[j] - e2[j - 1] / work[ji] - c__[ji];
#line 432 "dlaebz.f"
		    if (work[ji] <= *pivmin) {
#line 433 "dlaebz.f"
			++iwork[ji];
/* Computing MIN */
#line 434 "dlaebz.f"
			d__1 = work[ji], d__2 = -(*pivmin);
#line 434 "dlaebz.f"
			work[ji] = min(d__1,d__2);
#line 435 "dlaebz.f"
		    }
#line 436 "dlaebz.f"
/* L50: */
#line 436 "dlaebz.f"
		}
#line 437 "dlaebz.f"
/* L60: */
#line 437 "dlaebz.f"
	    }

#line 439 "dlaebz.f"
	    if (*ijob <= 2) {

/*              IJOB=2: Choose all intervals containing eigenvalues. */

#line 443 "dlaebz.f"
		klnew = kl;
#line 444 "dlaebz.f"
		i__2 = kl;
#line 444 "dlaebz.f"
		for (ji = kf; ji <= i__2; ++ji) {

/*                 Insure that N(w) is monotone */

/* Computing MIN */
/* Computing MAX */
#line 448 "dlaebz.f"
		    i__5 = nab[ji + nab_dim1], i__6 = iwork[ji];
#line 448 "dlaebz.f"
		    i__3 = nab[ji + (nab_dim1 << 1)], i__4 = max(i__5,i__6);
#line 448 "dlaebz.f"
		    iwork[ji] = min(i__3,i__4);

/*                 Update the Queue -- add intervals if both halves */
/*                 contain eigenvalues. */

#line 454 "dlaebz.f"
		    if (iwork[ji] == nab[ji + (nab_dim1 << 1)]) {

/*                    No eigenvalue in the upper interval: */
/*                    just use the lower interval. */

#line 459 "dlaebz.f"
			ab[ji + (ab_dim1 << 1)] = c__[ji];

#line 461 "dlaebz.f"
		    } else if (iwork[ji] == nab[ji + nab_dim1]) {

/*                    No eigenvalue in the lower interval: */
/*                    just use the upper interval. */

#line 466 "dlaebz.f"
			ab[ji + ab_dim1] = c__[ji];
#line 467 "dlaebz.f"
		    } else {
#line 468 "dlaebz.f"
			++klnew;
#line 469 "dlaebz.f"
			if (klnew <= *mmax) {

/*                       Eigenvalue in both intervals -- add upper to */
/*                       queue. */

#line 474 "dlaebz.f"
			    ab[klnew + (ab_dim1 << 1)] = ab[ji + (ab_dim1 << 
				    1)];
#line 475 "dlaebz.f"
			    nab[klnew + (nab_dim1 << 1)] = nab[ji + (nab_dim1 
				    << 1)];
#line 476 "dlaebz.f"
			    ab[klnew + ab_dim1] = c__[ji];
#line 477 "dlaebz.f"
			    nab[klnew + nab_dim1] = iwork[ji];
#line 478 "dlaebz.f"
			    ab[ji + (ab_dim1 << 1)] = c__[ji];
#line 479 "dlaebz.f"
			    nab[ji + (nab_dim1 << 1)] = iwork[ji];
#line 480 "dlaebz.f"
			} else {
#line 481 "dlaebz.f"
			    *info = *mmax + 1;
#line 482 "dlaebz.f"
			}
#line 483 "dlaebz.f"
		    }
#line 484 "dlaebz.f"
/* L70: */
#line 484 "dlaebz.f"
		}
#line 485 "dlaebz.f"
		if (*info != 0) {
#line 485 "dlaebz.f"
		    return 0;
#line 485 "dlaebz.f"
		}
#line 487 "dlaebz.f"
		kl = klnew;
#line 488 "dlaebz.f"
	    } else {

/*              IJOB=3: Binary search.  Keep only the interval containing */
/*                      w   s.t. N(w) = NVAL */

#line 493 "dlaebz.f"
		i__2 = kl;
#line 493 "dlaebz.f"
		for (ji = kf; ji <= i__2; ++ji) {
#line 494 "dlaebz.f"
		    if (iwork[ji] <= nval[ji]) {
#line 495 "dlaebz.f"
			ab[ji + ab_dim1] = c__[ji];
#line 496 "dlaebz.f"
			nab[ji + nab_dim1] = iwork[ji];
#line 497 "dlaebz.f"
		    }
#line 498 "dlaebz.f"
		    if (iwork[ji] >= nval[ji]) {
#line 499 "dlaebz.f"
			ab[ji + (ab_dim1 << 1)] = c__[ji];
#line 500 "dlaebz.f"
			nab[ji + (nab_dim1 << 1)] = iwork[ji];
#line 501 "dlaebz.f"
		    }
#line 502 "dlaebz.f"
/* L80: */
#line 502 "dlaebz.f"
		}
#line 503 "dlaebz.f"
	    }

#line 505 "dlaebz.f"
	} else {

/*           End of Parallel Version of the loop */

/*           Begin of Serial Version of the loop */

#line 511 "dlaebz.f"
	    klnew = kl;
#line 512 "dlaebz.f"
	    i__2 = kl;
#line 512 "dlaebz.f"
	    for (ji = kf; ji <= i__2; ++ji) {

/*              Compute N(w), the number of eigenvalues less than w */

#line 516 "dlaebz.f"
		tmp1 = c__[ji];
#line 517 "dlaebz.f"
		tmp2 = d__[1] - tmp1;
#line 518 "dlaebz.f"
		itmp1 = 0;
#line 519 "dlaebz.f"
		if (tmp2 <= *pivmin) {
#line 520 "dlaebz.f"
		    itmp1 = 1;
/* Computing MIN */
#line 521 "dlaebz.f"
		    d__1 = tmp2, d__2 = -(*pivmin);
#line 521 "dlaebz.f"
		    tmp2 = min(d__1,d__2);
#line 522 "dlaebz.f"
		}

#line 524 "dlaebz.f"
		i__3 = *n;
#line 524 "dlaebz.f"
		for (j = 2; j <= i__3; ++j) {
#line 525 "dlaebz.f"
		    tmp2 = d__[j] - e2[j - 1] / tmp2 - tmp1;
#line 526 "dlaebz.f"
		    if (tmp2 <= *pivmin) {
#line 527 "dlaebz.f"
			++itmp1;
/* Computing MIN */
#line 528 "dlaebz.f"
			d__1 = tmp2, d__2 = -(*pivmin);
#line 528 "dlaebz.f"
			tmp2 = min(d__1,d__2);
#line 529 "dlaebz.f"
		    }
#line 530 "dlaebz.f"
/* L90: */
#line 530 "dlaebz.f"
		}

#line 532 "dlaebz.f"
		if (*ijob <= 2) {

/*                 IJOB=2: Choose all intervals containing eigenvalues. */

/*                 Insure that N(w) is monotone */

/* Computing MIN */
/* Computing MAX */
#line 538 "dlaebz.f"
		    i__5 = nab[ji + nab_dim1];
#line 538 "dlaebz.f"
		    i__3 = nab[ji + (nab_dim1 << 1)], i__4 = max(i__5,itmp1);
#line 538 "dlaebz.f"
		    itmp1 = min(i__3,i__4);

/*                 Update the Queue -- add intervals if both halves */
/*                 contain eigenvalues. */

#line 544 "dlaebz.f"
		    if (itmp1 == nab[ji + (nab_dim1 << 1)]) {

/*                    No eigenvalue in the upper interval: */
/*                    just use the lower interval. */

#line 549 "dlaebz.f"
			ab[ji + (ab_dim1 << 1)] = tmp1;

#line 551 "dlaebz.f"
		    } else if (itmp1 == nab[ji + nab_dim1]) {

/*                    No eigenvalue in the lower interval: */
/*                    just use the upper interval. */

#line 556 "dlaebz.f"
			ab[ji + ab_dim1] = tmp1;
#line 557 "dlaebz.f"
		    } else if (klnew < *mmax) {

/*                    Eigenvalue in both intervals -- add upper to queue. */

#line 561 "dlaebz.f"
			++klnew;
#line 562 "dlaebz.f"
			ab[klnew + (ab_dim1 << 1)] = ab[ji + (ab_dim1 << 1)];
#line 563 "dlaebz.f"
			nab[klnew + (nab_dim1 << 1)] = nab[ji + (nab_dim1 << 
				1)];
#line 564 "dlaebz.f"
			ab[klnew + ab_dim1] = tmp1;
#line 565 "dlaebz.f"
			nab[klnew + nab_dim1] = itmp1;
#line 566 "dlaebz.f"
			ab[ji + (ab_dim1 << 1)] = tmp1;
#line 567 "dlaebz.f"
			nab[ji + (nab_dim1 << 1)] = itmp1;
#line 568 "dlaebz.f"
		    } else {
#line 569 "dlaebz.f"
			*info = *mmax + 1;
#line 570 "dlaebz.f"
			return 0;
#line 571 "dlaebz.f"
		    }
#line 572 "dlaebz.f"
		} else {

/*                 IJOB=3: Binary search.  Keep only the interval */
/*                         containing  w  s.t. N(w) = NVAL */

#line 577 "dlaebz.f"
		    if (itmp1 <= nval[ji]) {
#line 578 "dlaebz.f"
			ab[ji + ab_dim1] = tmp1;
#line 579 "dlaebz.f"
			nab[ji + nab_dim1] = itmp1;
#line 580 "dlaebz.f"
		    }
#line 581 "dlaebz.f"
		    if (itmp1 >= nval[ji]) {
#line 582 "dlaebz.f"
			ab[ji + (ab_dim1 << 1)] = tmp1;
#line 583 "dlaebz.f"
			nab[ji + (nab_dim1 << 1)] = itmp1;
#line 584 "dlaebz.f"
		    }
#line 585 "dlaebz.f"
		}
#line 586 "dlaebz.f"
/* L100: */
#line 586 "dlaebz.f"
	    }
#line 587 "dlaebz.f"
	    kl = klnew;

#line 589 "dlaebz.f"
	}

/*        Check for convergence */

#line 593 "dlaebz.f"
	kfnew = kf;
#line 594 "dlaebz.f"
	i__2 = kl;
#line 594 "dlaebz.f"
	for (ji = kf; ji <= i__2; ++ji) {
#line 595 "dlaebz.f"
	    tmp1 = (d__1 = ab[ji + (ab_dim1 << 1)] - ab[ji + ab_dim1], abs(
		    d__1));
/* Computing MAX */
#line 596 "dlaebz.f"
	    d__3 = (d__1 = ab[ji + (ab_dim1 << 1)], abs(d__1)), d__4 = (d__2 =
		     ab[ji + ab_dim1], abs(d__2));
#line 596 "dlaebz.f"
	    tmp2 = max(d__3,d__4);
/* Computing MAX */
#line 597 "dlaebz.f"
	    d__1 = max(*abstol,*pivmin), d__2 = *reltol * tmp2;
#line 597 "dlaebz.f"
	    if (tmp1 < max(d__1,d__2) || nab[ji + nab_dim1] >= nab[ji + (
		    nab_dim1 << 1)]) {

/*              Converged -- Swap with position KFNEW, */
/*                           then increment KFNEW */

#line 603 "dlaebz.f"
		if (ji > kfnew) {
#line 604 "dlaebz.f"
		    tmp1 = ab[ji + ab_dim1];
#line 605 "dlaebz.f"
		    tmp2 = ab[ji + (ab_dim1 << 1)];
#line 606 "dlaebz.f"
		    itmp1 = nab[ji + nab_dim1];
#line 607 "dlaebz.f"
		    itmp2 = nab[ji + (nab_dim1 << 1)];
#line 608 "dlaebz.f"
		    ab[ji + ab_dim1] = ab[kfnew + ab_dim1];
#line 609 "dlaebz.f"
		    ab[ji + (ab_dim1 << 1)] = ab[kfnew + (ab_dim1 << 1)];
#line 610 "dlaebz.f"
		    nab[ji + nab_dim1] = nab[kfnew + nab_dim1];
#line 611 "dlaebz.f"
		    nab[ji + (nab_dim1 << 1)] = nab[kfnew + (nab_dim1 << 1)];
#line 612 "dlaebz.f"
		    ab[kfnew + ab_dim1] = tmp1;
#line 613 "dlaebz.f"
		    ab[kfnew + (ab_dim1 << 1)] = tmp2;
#line 614 "dlaebz.f"
		    nab[kfnew + nab_dim1] = itmp1;
#line 615 "dlaebz.f"
		    nab[kfnew + (nab_dim1 << 1)] = itmp2;
#line 616 "dlaebz.f"
		    if (*ijob == 3) {
#line 617 "dlaebz.f"
			itmp1 = nval[ji];
#line 618 "dlaebz.f"
			nval[ji] = nval[kfnew];
#line 619 "dlaebz.f"
			nval[kfnew] = itmp1;
#line 620 "dlaebz.f"
		    }
#line 621 "dlaebz.f"
		}
#line 622 "dlaebz.f"
		++kfnew;
#line 623 "dlaebz.f"
	    }
#line 624 "dlaebz.f"
/* L110: */
#line 624 "dlaebz.f"
	}
#line 625 "dlaebz.f"
	kf = kfnew;

/*        Choose Midpoints */

#line 629 "dlaebz.f"
	i__2 = kl;
#line 629 "dlaebz.f"
	for (ji = kf; ji <= i__2; ++ji) {
#line 630 "dlaebz.f"
	    c__[ji] = (ab[ji + ab_dim1] + ab[ji + (ab_dim1 << 1)]) * .5;
#line 631 "dlaebz.f"
/* L120: */
#line 631 "dlaebz.f"
	}

/*        If no more intervals to refine, quit. */

#line 635 "dlaebz.f"
	if (kf > kl) {
#line 635 "dlaebz.f"
	    goto L140;
#line 635 "dlaebz.f"
	}
#line 637 "dlaebz.f"
/* L130: */
#line 637 "dlaebz.f"
    }

/*     Converged */

#line 641 "dlaebz.f"
L140:
/* Computing MAX */
#line 642 "dlaebz.f"
    i__1 = kl + 1 - kf;
#line 642 "dlaebz.f"
    *info = max(i__1,0);
#line 643 "dlaebz.f"
    *mout = kl;

#line 645 "dlaebz.f"
    return 0;

/*     End of DLAEBZ */

} /* dlaebz_ */

