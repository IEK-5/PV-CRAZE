/* this is a collection of routines for fitting splines from 
 * the math77/mathc90 library. I am not sure I should do this. On the 
 * one hand I should just put in the math77/mathc90 library as a 
 * rerquirement. However, it seems few systems distribute readily 
 * installable math77 libraries. Also the mathc90 library does not 
 * provide a very portable build system. So if I want to distribute my 
 * code I practically have to distribute the mathc90 library with it. 
 * In that case I might just rip out this one routine I need. So, 
 * sorry... Some minor modifications regarding double precission and 
 * error handling.
 */
/* The following routines are provided under the following License
 * Copyright © 1996 California Institute of Technology, Pasadena, 
 * California. ALL RIGHTS RESERVED. Based on Government Sponsored 
 * Research NAS7-03001. Redistribution and use in source and binary 
 * forms, with or without modification, are permitted provided that the 
 * following conditions are met: Redistributions of source code must 
 * retain this copyright notice, this list of conditions and the 
 * following disclaimer. Redistributions in binary form must reproduce 
 * the above copyright notice, this list of conditions and the 
 * following disclaimer in the documentation and/or other materials 
 * provided with the distribution. Neither the name of the California 
 * Institute of Technology (Caltech) nor the names of its contributors 
 * may be used to endorse or promote products derived from this 
 * software without specific prior written permission. THIS SOFTWARE IS 
 * PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF 
 * SUCH DAMAGE.
 * 
 * For those codes indicated with a Math a la Carte copyright, the same 
 * rules apply, except without the full force of the Caltech legal team.
 * 
 * When citing this software we request that you also mention the names 
 * of the people who wrote the software you are using.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "splinefit.h"

/* some macros */
#define SQ(x)   ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define DOCNT(i,t,n) (_d_l=(n), (_d_m=(t-(i)+_d_l)/_d_l) > 0 ? _d_m : 0L )
#define  max(a,b) ( (a) > (b) ? (a) : (b) )
#define  min(a,b) ( (a) < (b) ? (a) : (b) )
		/* PARAMETER translations */
#define TRUE   1
#define FALSE  0
#define	NBAND	4
#define	NBAND1	(NBAND + 1)
#define	ONE	1.0e0
#define	ZERO	0.0e0
#define	C18	18.0e0
#define	C36	36.0e0
#define	C6	6.0e0
		/* end of PARAMETER translations */

typedef long   LOGICAL32;


double /*FUNCTION*/ snrm2(
long n,
double x[],
long incx)
{
	long int _d_l, _d_m, _do0, _do1, _do2, _do3, i, 
	 id, iemn, iemx, nn;
	double snrm2_v, sum, sumx, xa;
	static double asmall, dbig, dsmall, fbig, fsmall, tbig, tsmall;
	static double abig = 0.e0;
		/* OFFSET Vectors w/subscript range: 1 to dimension */
	double *const X = &x[0] - 1;
		/* end of OFFSET VECTORS */
 
	/* Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
	 * ALL RIGHTS RESERVED.
	 * Based on Government Sponsored Research NAS7-03001.
	 *>> 1998-05-11 SNRM2  Krogh   Minor changes for conversion to C.
	 *>> 1996-08-29 SNRM2  Krogh   Coded an entirely different algorithm
	 *>> ....
	 *>> 1985-08-02 SNRM2  Lawson  Initial code.
	 *
	 *     EUCLIDEAN NORM OF THE N-VECTOR STORED IN X() WITH STORAGE
	 *     INCREMENT INCX .
	 *     IF    N .LE. 0 RETURN WITH RESULT = 0.
	 *     IF N .GE. 1 THEN INCX MUST BE .GE. 1
	 *
	 *           C.L.LAWSON, 1978 JAN 08
	 *
	 * New algorithm avoids underflow as well as overflow, and avoids divides
	 * as much as possible.  F. Krogh, August 29, 1996.
	 *
	 * ************************* Variable Definitions ***********************
	 *
	 * Let b be the base for floating point, b**E be the smallest integer
	 * power of b that overflows, b**e be the smallest integer power of b
	 * that does not underflow, and d be the number of base b digits in a
	 * a floating point number.  The descriptions below use b, e, E, and d.
	 *
	 * ABIG = b ** ((E+5) / 4)).    This is used as a factor to get the
	 *   final answer when combining SUMX from big XA, with SUM.
	 * ASMALL = b ** ((e+d-E+5) / 4)).    This is used as a factor to get the
	 *   final answer when combining SUMX from small XA, with SUM.
	 * DBIG = b ** (2*((E+5)/4)).  This is used as a factor to get the
	 *   final answer when only SUMX formed from big XA's is needed.
	 * DSMALL = b ** (2*((e+d-E+5)/4)).  This is used as a factor to get the
	 *   final answer when only SUMX formed from big XA's is needed.
	 * FBIG = b ** (-2*((E+5)/4)).  This is used as a multiplier when
	 *   accumulating SUMX for big XA.
	 * FSMALL = b ** (-2*((e+d-E+5)/4)).  This is used as a multiplier when
	 *   accumulating SUMX for small XA.
	 * I      Temporary index.
	 * ID     Number of base I1MACH(10) digits in Floating point number.
	 * IEMN   The minimum floating point exponent.
	 * IEMX   The maximum floating point exponent.
	 * INCX   Input, the increment (>0) between elements of X.
	 * N      Input, the number of elements in X.
	 * NN     N * INCX = last index processed.
	 * SUM    Place where sum of XA**2 is accumlated.
	 * SUMX   Place where scaled sum of XA**2 is accumlated.  In first loop
	 *   this is for XA**2 that would likely underflow, and in second loop
	 *   this is for XA**2 that would likely overflow.
	 * TBIG = b ** ((E - d - 1) / 2).  If XA > TBIG, number is "big".  Note
	 *   that for such XA's, (FBIG * XA) ** 2 should not underflow and the
	 *   accumlation overflows only if the final result would.
	 * TSMALL = b ** ((e+1) / 2).  If XA <= TSMALL, number is "small".  Note
	 *   that for such XA's, (FSMALL * XA)**2 should not underflow and the
	 *   accumlation should not overlflow.
	 * X      Input, we are getting the L2 norm of X.
	 * XA     Contains base of floating point numbers when getting saved
	 *   parameters.  Later contains abs(X(I)).
	 *     ------------------------------------------------------------------
	 *--S replaces "?": ?NRM2
	 *     ------------------------------------------------------------------ */
	/* Values of these saved parameters for D.P. IEEE arithmetic are:
	 *   ABIG= .2315841784746324E+78     DBIG= .5363123171977043E+155
	 *   FBIG= .1864585182800050E-154    TBIG= .9989595361011182E+146
	 * ASMALL= .4887898181599363E-149  DSMALL= .2389154863368240E-298
	 * FSMALL= .4185580496821357E+299  TSMALL= .2983336292480080E-153
	 * Values of these saved parameters for S.P. IEEE arithmetic are:
	 *   ABIG= .8589935E+10      DBIG= .7378698E+20
	 *   FBIG= .1355253E-19      TBIG= .2251800E+16
	 * ASMALL= .1387779E-16    DSMALL= .1925930E-33
	 * FSMALL= .5192297E+34    TSMALL= .2168404E-18
	 * */
	/*     ------------------------------------------------------------------
	 *
	 * */
	if (abig == 0.e0)
	{
		/*++ Code for (.N. == 'D') is inactive
		 *         IEMX = I1MACH(16)
		 *         IEMN = I1MACH(15)
		 *         ID = I1MACH(14)
		 *++ Code for (.N. == 'S') is active */
		iemx = DBL_MAX_EXP;
		iemn = DBL_MIN_EXP;
		id = DBL_MANT_DIG;
		/*++ END */
		xa = (double)( FLT_RADIX );
		abig = pow(xa,(iemx + 5)/4);
		dbig = SQ(abig);
		fbig = 1.e0/dbig;
		tbig = pow(xa,(iemx - id - 1)/2);
		asmall = pow(xa,(iemn + id - iemx + 5)/4);
		dsmall = SQ(asmall);
		fsmall = 1.e0/dsmall;
		tsmall = pow(xa,(iemn + 1)/2);
	}
	sum = 0.e0;
	if (n > 0)
	{
		nn = n*incx;
		sumx = 0.e0;
		/*                      Loop when no big number yet encountered. */
		for (i = 1, _do0=DOCNT(i,nn,_do1 = incx); _do0 > 0; i += _do1, _do0--)
		{
			xa = fabs( X[i] );
			if (xa < tsmall)
			{
				sumx += pow(fsmall*xa,2);
			}
			else
			{
				if (xa > tbig)
					goto L_200;
				sum += SQ(xa);
			}
		}
		if (sum != 0.e0)
		{
			if (sumx >= 1.e0)
			{
				if (sum < 1.e0)
				{
					sum = asmall*sqrtf( fsmall*sum + dsmall*sumx );
					goto L_400;
				}
			}
			sum = sqrtf( sum );
		}
		else
		{
			sum = dsmall*sqrtf( sumx );
		}
		goto L_400;
 
L_200:
		sumx = 0.e0;
		/*                      Loop when we have at least one big number. */
		for (i = i, _do2=DOCNT(i,nn,_do3 = incx); _do2 > 0; i += _do3, _do2--)
		{
			xa = fabs( X[i] );
			if (xa > tsmall)
			{
				if (xa > tbig)
				{
					sumx += pow(fbig*xa,2);
				}
				else
				{
					sum += SQ(xa);
				}
			}
		}
		if ((sumx <= 1.e10) && (sum >= 1.e-10))
		{
			sum = abig*sqrtf( fbig*sum + dbig*sumx );
		}
		else
		{
			sum = dbig*sqrtf( sumx );
		}
	}
L_400:
	;
	snrm2_v = sum;
	return( snrm2_v );
} /* end of function */
 


void /*FUNCTION*/ shtcc(
long mode,
long lpivot,
long l1,
long m,
double u[],
double *uparam,
double c[],
long ldc,
long nvc)
{
	long int i, iul0, j, jcbase, jcpiv;
	double b, binv, fac, hold, sum, vnorm;
		/* OFFSET Vectors w/subscript range: 1 to dimension */
	double *const C = &c[0] - 1;
	double *const U = &u[0] - 1;
		/* end of OFFSET VECTORS */
 
	/* Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
	 * ALL RIGHTS RESERVED.
	 * Based on Government Sponsored Research NAS7-03001.
	 *>> 1996-03-30 SHTCC  Krogh   Added external statement.
	 *>> 1994-11-11 SHTCC  Krogh   Declared all vars.
	 *>> 1994-10-20 SHTCC  Krogh  Changes to use M77CON
	 *>> 1987-08-19 SHTCC  Lawson  Initial code.
	 *--S replaces "?": ?HTCC, ?NRM2
	 *
	 *        Construction and/or application of a single Householder
	 *     transformation..     Q = I + U*(U**T)/b
	 *     where I is the MxM identity matrix, b is a scalar, and U is an
	 *     M-dimensional Householder vector.
	 *        All vectors are M-vectors but only the components in positions
	 *     LPIVOT, and L1 through M, will be referenced.
	 *        This version, identified by CC at the end of its name, is
	 *     specialized for the Column-Column case, i.e. U() is a vector or
	 *     a column of a matrix and C() is regarded a containing a set
	 *     of column vectors to which transformations will be applied.
	 *     ------------------------------------------------------------------
	 *                         Subroutine arguments
	 *
	 *     MODE  [in]  = 1 OR 2   When MODE = 1 this subr determines the
	 *           parameters for a Householder transformation and applies
	 *           the transformation to NVC vectors.  When MODE = 2 this
	 *           subr applies a previously determined Householder
	 *           transformation.
	 *     LPIVOT  [in]  The index of the pivot element.
	 *     L1,M  [in]  If L1 .le. M  elements LPIVOT and L1 through M will
	 *           be referenced.  If L1 .gt. M the subroutine returns
	 *           immediately.  This may be regarded
	 *           as performing an identity transformation.
	 *     U()  [inout]  Contains an M-dimensional vector with storage
	 *           spacing of 1 between elements.
	 *           When MODE = 1 this is the vector from which Householder
	 *           parameters are to be determined.
	 *           When MODE = 2 this is the result from previous computation
	 *           with MODE = 1.
	 *     UPARAM  [inout]  Holds a value that supplements the
	 *           contents of U() to complete the definition of a
	 *           Householder transformation.  Computed when MODE = 1 and
	 *           reused when MODE = 2.
	 *           UPARAM is the pivot component of the Householder U-vector.
	 *     C()  [inout]  On entry contains a set of NVC M-vectors to which a
	 *          Householder transformation is to be applied.
	 *          On exit contains the set of transformed vectors.
	 *          These vectors are the columns of an M x NVC matrix in C(,).
	 *     LDC  [in]  Leading dimension of C(,).  Require LDC .ge. M.
	 *     NVC  [in]  Number of vectors in C(,) to be transformed.
	 *           If NVC .le. 0 no reference will be made to the array C(,).
	 *     ------------------------------------------------------------------
	 *     Subprograms referenced: SNRM2
	 *     ------------------------------------------------------------------
	 *          This code was originally developed by Charles L. Lawson and
	 *     Richard J. Hanson at Jet Propulsion Laboratory in 1973.  The
	 *     original code was described and listed in the book,
	 *
	 *                  Solving Least Squares Problems
	 *                  C. L. Lawson and R. J. Hanson
	 *                  Prentice-Hall, 1974
	 *
	 *     Feb, 1985, C. L. Lawson & S. Y. Chan, JPL.  Adapted code from the
	 *     Lawson & Hanson book to Fortran 77 for use in the JPL MATH77
	 *     library.
	 *     Prefixing subprogram names with S or D for s.p. or d.p. versions.
	 *     Using generic names for intrinsic functions.
	 *     Adding calls to BLAS and MATH77 error processing subrs in some
	 *     program units.
	 *     July, 1987. CLL.  Changed user interface so method of specifying
	 *     column/row storage options is more language-independent.
	 *     ------------------------------------------------------------------ */
	/*     ------------------------------------------------------------------ */
	if ((0 >= lpivot || lpivot >= l1) || l1 > m)
		return;
	if (mode == 1)
	{
		/*                            ****** CONSTRUCT THE TRANSFORMATION. ****** */
		iul0 = l1 - 1;
		if (iul0 == lpivot)
		{
			vnorm = snrm2( m - l1 + 2, &U[iul0], 1 );
		}
		else
		{
			hold = U[iul0];
			U[iul0] = U[lpivot];
			vnorm = snrm2( m - l1 + 2, &U[iul0], 1 );
			U[iul0] = hold;
		}
 
		if (U[lpivot] > ZERO)
			vnorm = -vnorm;
		*uparam = U[lpivot] - vnorm;
		U[lpivot] = vnorm;
	}
	/*            ****** Apply the transformation  I + U*(U**T)/B  to C. ****
	 * */
	if (nvc <= 0)
		return;
	b = *uparam*U[lpivot];
	/*                                 Here B .le. 0.  If B  =  0., return. */
	if (b == ZERO)
		return;
	binv = ONE/b;
	jcbase = 0;
 
	for (j = 1; j <= nvc; j++)
	{
		jcpiv = jcbase + lpivot;
		sum = C[jcpiv]**uparam;
		for (i = l1; i <= m; i++)
		{
			sum += C[jcbase + i]*U[i];
		}
		if (sum != ZERO)
		{
			fac = sum*binv;
			C[jcpiv] += fac**uparam;
			for (i = l1; i <= m; i++)
			{
				C[jcbase + i] += fac*U[i];
			}
		}
		jcbase += ldc;
	}
	return;
} /* end of function */
 




void /*FUNCTION*/ sc2bas(
double x,
long nderiv,
long iseg,
LOGICAL32 *newseg,
double b[],
long nb,
double p[])
{
	long int i, i2, ib, ib1, ib2, is, ism4, iu, j, mb, ndp1;
	double cs, flnseg, h, xs;
	static double scale, u[8];
		/* OFFSET Vectors w/subscript range: 1 to dimension */
	double *const B = &b[0] - 1;
	double *const P = &p[0] - 1;
	double *const U = &u[0] - 1;
		/* end of OFFSET VECTORS */
 
	/* Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
	 * ALL RIGHTS RESERVED.
	 * Based on Government Sponsored Research NAS7-03001.
	 *>> 1994-11-11 SC2BAS Krogh  Declared all vars.
	 *>> 1994-10-20 SC2BAS Krogh  Changes to use M77CON
	 *>> 1987-12-09 SC2BAS Lawson Initial code.
	 *--S replaces "?": ?C2BAS
	 *
	 *          Given X in the segment ISEG, this subr evaluates at X the
	 *     four cubic B-spline basis functions that are not identically
	 *     zero throughout segment ISEG.
	 *          This is a reorganization of a set of three subroutines
	 *     SBAS1, SBAS2, & SBAS3 used in LIB*JPL$ from about 1970 on.
	 *     The reorganization improves portability and brings
	 *     conformance to Fortran 77.
	 *     More features are avalable from a similar subr in a package
	 *     due to Carl de Boor, however this subr is being preserved
	 *     for convenience of usage with some other JPL subrs --
	 *     particularly the size and interpretation of the B() array.
	 *          C. L. Lawson, JPL, 1984 July 9, July 1987.
	 *     ------------------------------------------------------------------
	 *               A R G U M E N T S
	 *
	 *     X  [in]  Argument at which basis fcns are to be evaluated.
	 *           Must satisfy B(ISEG) .le. X .le. B(ISEG+1)
	 *
	 *     NDERIV  [in]   = 0, 1, 2, or 3.  Indicates the order of derivative
	 *              to be computed.  0 selects function value,
	 *              1 selects 1st derivative, etc.
	 *
	 *     ISEG  [in] Index of current segment. Require 1 .le. ISEG .le. NB-1
	 *               Segment ISEG is the interval [B(ISEG), B(ISEG+1)].
	 *               The current X must be in this segment.
	 *
	 *     NEWSEG  [inout]  = .true. means ISEG is not the same as on
	 *             previous call.  In this case this subr will set new values
	 *             into U() and reset NEWSEG = .false.  The user prog must
	 *             set NEWSEG = .true. on the first call to this subr, and
	 *             on any subsequent call at which the user prog has set
	 *             ISEG to a new value.
	 *
	 *             = .false. means this subr can assume the values in U()
	 *             are correct for the current ISEG.
	 *
	 *     B()  [in]  Array of knots, including the left and right
	 *               endpoints.  Contains B(I), I = 1, NB.  Set by user.
	 *               Not changed by this subr.  Must be strictly increasing,
	 *               i.e. require B(i) .lt. B(i+1), i = 1, NB-1.
	 *
	 *     NB  [in]  Number of knots, including endpoints. The number
	 *               of degrees of freedom of the set of cubic splines
	 *               defined over this set of knots is NB+2.
	 *
	 *     P()  [out]  Array of length 4 into which this subr will store the
	 *             values of the 4 basis fcns that are not identically
	 *             zero over segment ISEG.  The indices of these basis
	 *             fcns are ISEG, ISEG+1, ISEG+2, & ISEG+3.  The values
	 *             will be stored in (P(i), i=1,4).
	 *     ------------------------------------------------------------------
	 *               IMPORTANT INTERNAL VARIABLES
	 *
	 *     U()       Array of length 8 in which this subr saves a scaled
	 *               record of the knots near the current segment.
	 *
	 *     SCALE  Scale factor for abcissas.  Computed as the reciprocal
	 *            of the average segment length.
	 *     ------------------------------------------------------------------
	 * */
	/*     ------------------------------------------------------------------ */
	mb = nb;
	/*                     The following statement does a type conversion. */
	flnseg = mb - 1;
	if (*newseg)
	{
		*newseg = FALSE;
		scale = flnseg/(B[nb] - B[1]);
		is = iseg;
		ism4 = is - 4;
		ib1 = max( 1, is - 3 );
		ib2 = min( mb, is + 4 );
		/*                               COPY PART OF B() INTO U() */
		for (ib = ib1; ib <= ib2; ib++)
		{
			iu = ib - ism4;
			U[iu] = scale*B[ib];
		}
		/*                               ADJOIN ARTIFICIAL SEGMENTS AT LEFT END */
		iu = ib1 - ism4 - 1;
		if (iu != 0)
		{
			h = U[iu + 2] - U[iu + 1];
L_30:
			U[iu] = U[iu + 1] - h;
			iu -= 1;
			if (iu > 0)
				goto L_30;
		}
		/*                               ADJOIN ARTIFICIAL SEGMENTS AT RIGHT END */
		iu = ib2 - ism4 + 1;
		if (iu <= 8)
		{
			h = U[iu - 1] - U[iu - 2];
L_50:
			U[iu] = U[iu - 1] + h;
			iu += 1;
			if (iu <= 8)
				goto L_50;
		}
	}
	/*     ------------------------------------------------------------------ */
	ndp1 = nderiv + 1;
	xs = x*scale;
	switch (ndp1)
	{
		case 1: goto L_70;
		case 2: goto L_90;
		case 3: goto L_110;
		case 4: goto L_130;
	}
L_70:
	for (i = 1; i <= 4; i++)
	{
		P[i] = pow(xs - U[i],3);
	}
	cs = C6;
	goto L_150;
L_90:
	for (i = 1; i <= 4; i++)
	{
		P[i] = pow(xs - U[i],2);
	}
	cs = C18*scale;
	goto L_150;
L_110:
	for (i = 1; i <= 4; i++)
	{
		P[i] = xs - U[i];
	}
	cs = C36*SQ(scale);
	goto L_150;
L_130:
	for (i = 1; i <= 4; i++)
	{
		P[i] = C36;
	}
	cs = CUBE(scale);
L_150:
	;
 
	/*                            COMPUTE FOURTH DIVIDED DIFFERENCES
	 * */
	for (j = 1; j <= 4; j++)
	{
		for (i = 1; i <= 3; i++)
		{
			i2 = i + j;
			P[i] = (P[i + 1] - P[i])/(U[i2] - U[i]);
		}
		P[4] = -P[4]/(U[j + 4] - U[4]);
	}
 
	for (i = 1; i <= 4; i++)
	{
		P[i] *= cs;
	}
	return;
} /* end of function */
 

void /*FUNCTION*/ sbsol(
long mode,
double *g,
long ldg,
long nb,
long ir,
long jtprev,
double x[],
long n,
double *rnorm,
long *ierr3)
{
#define G(I_,J_)	(*(g+(I_)*(ldg)+(J_)))
	long int i, i1, i2, ie, ii, irm1, ix, j, jg, l, np1;
	double rsq, s;
		/* OFFSET Vectors w/subscript range: 1 to dimension */
	double *const X = &x[0] - 1;
		/* end of OFFSET VECTORS */
 
	/* Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
	 * ALL RIGHTS RESERVED.
	 * Based on Government Sponsored Research NAS7-03001.
	 *>> 1994-11-11 SBSOL  Krogh  Declared all vars.
	 *>> 1994-10-20 SBSOL  Krogh  Changes to use M77CON
	 *>> 1987-11-24 SBSOL  Lawson Initial code.
	 *--S replaces "?": ?BSOL
	 *     Solution of banded least squares problem following use of _BACC to
	 *     accumulate the data.
	 *     ------------------------------------------------------------------
	 *                            Subroutine arguments
	 *
	 *     MODE  [in]  Flag selecting operation to be done.
	 *          = 1     SOLVE R*X=Y   WHERE R AND Y ARE IN THE G( ) ARRAY
	 *                  AND X WILL BE STORED IN THE X( ) ARRAY.
	 *            2     SOLVE (R**T)*X=Y   WHERE R IS IN G( ),
	 *                  Y IS INITIALLY IN X( ), AND X REPLACES Y IN X( ),
	 *            3     SOLVE R*X=Y   WHERE R IS IN G( ).
	 *                  Y IS INITIALLY IN X( ), AND X REPLACES Y IN X( ).
	 *
	 *     G(,)  [in]  Array in which the transformed problem data has
	 *           been left by _BACC.
	 *     LDG   [in]  Leading dimensioning parameter for G(,).
	 *     NB    [in]  Bandwidth of data matrix.
	 *     IR    [in]  Must have value set by previous call to _BACC.
	 *     JTPREV    [in]  Must have value set by previous call to _BACC.
	 *     X()   [inout]  Array of length at least N.
	 *           On entry with MODE = 2 or 3 contains the y vector.
	 *           In all (non-error) cases contains the solution vector on
	 *           return.
	 *     N     [in]  Specifies the dimension of the desired solution
	 *           vector.  Can be set smaller than the max value the data
	 *           would support.
	 *     RNORM  [out]  Set to the euclidean norm of the residual vector
	 *           when MODE = 1.  Set to zero otherwise.
	 *     IERR3  [out]  =0 means no errors detected.
	 *                  =1 means a diagonal element is zero.
	 *     ------------------------------------------------------------------
	 *     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUN 12
	 *     Reference:   'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974
	 *     Comments in this code refer to Algorithm Step numbers on
	 *     pp 213-217 of the book.
	 *     ------------------------------------------------------------------
	 *     1984 July 11, C. L. Lawson, JPL.  Adapted code from book
	 *     to Fortran 77 for JPL MATH 77 library.
	 *     Prefixing subprogram names with S or D for s.p. or d.p. versions.
	 *     Using generic names for any intrinsic functions.
	 *     July 1987.  Using _HTCC in place of _H12.
	 *     ------------------------------------------------------------------
	 *     Subprograms called: IERM1, IERV1
	 *     ------------------------------------------------------------------ */
	/*     ------------------------------------------------------------------
	 * */
	*ierr3 = 0;
	*rnorm = ZERO;
	if (mode == 1)
	{
		/*                                   ********************* MODE = 1
		 *                                   ALG. STEP 26 */
		for (j = 1; j <= n; j++)
		{
			X[j] = G(nb,j - 1);
		}
		rsq = ZERO;
		np1 = n + 1;
		irm1 = ir - 1;
		if (np1 <= irm1)
		{
			for (j = np1; j <= irm1; j++)
			{
				rsq += SQ(G(nb,j - 1));
			}
			*rnorm = sqrtf( rsq );
		}
	}
 
	if (mode == 2)
		goto L_90;
	/*                                   ********************* MODE = 1 or 3
	 *                                   ALG. STEP 27 */
	for (ii = 1; ii <= n; ii++)
	{
		i = n + 1 - ii;
		/*                                   ALG. STEP 28 */
		s = ZERO;
		l = max( 0, i - jtprev );
		/*                                   ALG. STEP 29 */
		if (i != n)
		{
			/*                                   ALG. STEP 30 */
			ie = min( n + 1 - i, nb );
			for (j = 2; j <= ie; j++)
			{
				jg = j + l;
				ix = i - 1 + j;
				s += G(jg - 1,i - 1)*X[ix];
			}
		}
		/*                                   ALG. STEP 31 */
		if (G(l,i - 1) == ZERO)
			goto L_130;
		X[i] = (X[i] - s)/G(l,i - 1);
	}
	/*                                   ALG. STEP 32 */
	return;
	/*     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 *                                   ********************* MODE = 2 */
L_90:
	for (j = 1; j <= n; j++)
	{
		s = ZERO;
		if (j != 1)
		{
			i1 = max( 1, j - nb + 1 );
			i2 = j - 1;
			for (i = i1; i <= i2; i++)
			{
				l = j - i + 1 + max( 0, i - jtprev );
				s += X[i]*G(l - 1,i - 1);
			}
		}
		l = max( 0, j - jtprev );
		if (G(l,j - 1) == ZERO)
		{
			i = j;
			goto L_130;
		}
		X[j] = (X[j] - s)/G(l,j - 1);
	}
	return;
	/*     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
L_130:
	;
	*ierr3 = i; // B. Pieters: error message is raised in sc2fit
	return;
#undef	G
} /* end of function */
 
 
void /*FUNCTION*/ sbacc(
double *g,
long ldg,
long nb,
long *ir,
long mt,
long jt,
long *jtprev,
long *ierr2)
{
#define G(I_,J_)	(*(g+(I_)*(ldg)+(J_)))
	long int i, ie, ig, ig1, ig2, istep, j, jg, k, kh, kh1,
	 l, lp1, mh, nbp1;
	double uparam;
 
	/* Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
	 * ALL RIGHTS RESERVED.
	 * Based on Government Sponsored Research NAS7-03001.
	 *>> 1994-10-20 SBACC Krogh  Changes to use M77CON
	 *>> 1992-06-16 SBACC CLL
	 *>> 1989-10-20 CLL
	 *>> 1987-11-24 SBACC  Lawson  Initial code.
	 *     Sequential accumulation of data for a banded least-squares
	 *     problem.  For solution use _BSOL.
	 *     ------------------------------------------------------------------
	 *                            Subroutine arguments
	 *
	 *     G(,)  [inout]  Array in which user provides data and this subr
	 *              builds a packed representation of the orthogonally
	 *              transformed equivalent problem.
	 *     LDG   [in]  Leading dimensioning parameter for G(,).
	 *     NB    [in]  Bandwidth of data matrix.  Set on first call for a
	 *              new problem and not changed thereafter.
	 *     IR    [inout]  User must set IR = 1 on 1st call for a new problem.
	 *              User must not thereafter alter IR.  This subr will update
	 *              IR as IR := JT + min(NB+1, MT+max(IR-JT,0)).
	 *              IR identifies the row of G(,) in which user's new block
	 *              of data begins.
	 *     MT    [in]  Number of new rows of data being provided on the
	 *              current entry.
	 *     JT    [in]  Column index in the user's problem of the data being
	 *              provided in column 1 of G(,).
	 *     JTPREV    [inout]  Should not be set or altered by the user.
	 *              This subr will update JTPREV as JTPREV := JT.
	 *              JTPREV identifies a row in G(,) at which the mode of
	 *              packed storage changes.  Data in row JTPREV and below
	 *              are subject to possible change as future data are
	 *              introduced, whereas data above row JTPREV are not.
	 *     IERR2  [inout]  Error status indicator.  Must be set to zero on
	 *              initial call to this subr for a new problem.  Values on
	 *              return have meaning as follows:
	 *           0 means no errors detected.
	 *           1 means MT exceeds LDG - IR + 1
	 *           2 means JT .lt. JTPREV
	 *           3 means JT .gt. min(JTPREV + NB, IR) which causes the
	 *             problem to be structurally singular.  The transformations
	 *             by this subr can be completed in this case but the
	 *             resulting transformed problem cannot be solved simply
	 *             by use of _BSOL.
	 *           4 means the condition of Error 3 above applys but computa-
	 *             tion must be abandonded due to storage limitations
	 *             indicated by MT > LDG - JT + 1.
	 *           5 This subr will STOP because the input value of IERR2 was
	 *             nonzero.
	 *     ------------------------------------------------------------------
	 *     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUN 12
	 *     Reference:   'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974
	 *     Comments in this code refer to Algorithm Step numbers on
	 *     pp 213-217 of the book.
	 *     ------------------------------------------------------------------
	 *     1984 July 11, C. L. Lawson, JPL.  Adapted code from book
	 *     to Fortran 77 for JPL MATH 77 library.
	 *     Prefixing subprogram names with S or D for s.p. or d.p. versions.
	 *     Using generic names for any intrinsic functions.
	 *     July 1987.  Using _HTCC in place of _H12.
	 *     1989-10-20 CLL Moved integer declaration earlier to avoid warning
	 *     msg from Cray compiler.
	 *     1992-06-16 CLL  Using "min" fcn in arg list of call to _HTCC
	 *     to avoid "index out of range" when using a bounds checker, even
	 *     though no reference would be made to the "out of range" location.
	 *     ------------------------------------------------------------------
	 *--S replaces "?": ?BACC, ?HTCC
	 *     Both versions use IERM1, IERV1
	 *     ------------------------------------------------------------------ */
	/*     ------------------------------------------------------------------
	 *
	 *              ALG. STEPS 1-4 ARE PERFORMED EXTERNAL TO THIS SUBROUTINE.
	 * */
	if (mt <= 0)
		return;
	if (*ir == 1)
	{
		/*                                     Initialize for new problem. */
		*jtprev = 1;
		*ierr2 = 0;
	}
	else if (*ierr2 != 0)
	{		
		*ierr2 = 5;
		// exit(0); /* B. Pieters: what?! */
		return;
	}
 
	nbp1 = nb + 1;
	if (mt > ldg - *ir + 1)
	{
		*ierr2 = 1;
		return;
	}
	/*                                             ALG. STEP 5 */
	if (jt == *jtprev)
		goto L_70;
	/*                                             ALG. STEP 6 */
	if (jt < *jtprev)
	{
		*ierr2 = 2;
		return;
	}
	/*                                             ALG. STEP 7 */
	if (jt > min( *jtprev + nb, *ir ))
	{
		*ierr2 = 3;
 
		/*        Note:  At this point we know the problem will be structurally
		 *        singular and thus will not be solvable simply by use of _BSOL.
		 *        This does not preclude, however, completion of the transforma-
		 *        tion to the band triangular form and so we continue.  The
		 *        calling program will be informed of this condition by the set-
		 *        ting of of IERR2 = 3, and thus can quit or continue as desired.
		 *        If this is not due to a usage mistake
		 *        then the user should probably introduce more data, modify his
		 *        or her mathematical model, or use a stabilizing method such as
		 *        is discussed on pp. 218-219 of the L & H book.
		 * */
		if (mt > ldg - jt + 1)
		{
			*ierr2 = 4;
			return;
		}
		/*                                             ALG. STEPS 8-9 */
		for (i = 1; i <= mt; i++)
		{
			ig1 = jt + mt - i;
			ig2 = *ir + mt - i;
			for (j = 1; j <= nbp1; j++)
			{
				G(j - 1,ig1 - 1) = G(j - 1,ig2 - 1);
			}
		}
		/*                                             ALG. STEP 10 */
		ie = jt - *ir;
		for (i = 1; i <= ie; i++)
		{
			ig = *ir + i - 1;
			for (j = 1; j <= nbp1; j++)
			{
				G(j - 1,ig - 1) = ZERO;
			}
		}
		/*                                             ALG. STEP 11 */
		*ir = jt;
	}
 
	/*     Here we must shift some rows left by variable amounts.
	 *     min(NB-1,IR-JTPREV-1) is the number of rows to be shifted.
	 *     K = min(L,JT-JTPREV) is the extent of the left shift
	 *     for row JTPREV+L.
	 *
	 *                                             ALG. STEP 12-13 */
	for (l = 1; l <= min( nb - 1, *ir - *jtprev - 1 ); l++)
	{
		/*                                             ALG. STEP 14 */
		k = min( l, jt - *jtprev );
		/*                                             ALG. STEP 15 */
		lp1 = l + 1;
		ig = *jtprev + l;
		for (i = lp1; i <= nb; i++)
		{
			G(i - k - 1,ig - 1) = G(i - 1,ig - 1);
		}
		/*                                             ALG. STEP 16 */
		for (jg = nbp1 - k; jg <= nb; jg++)
		{
			G(jg - 1,ig - 1) = ZERO;
		}
	}
	/*                                             ALG. STEP 17 */
	*jtprev = jt;
	/*                                             ALG. STEPS 18-19 */
L_70:
	;
 
	/*     MH = No. of rows in the block to which the orthogonal
	 *          triangularization will be applied.
	 *     KH1 = Number of pivot columns to be used in the orghogonal
	 *           transformations.
	 * */
	mh = *ir + mt - *jtprev;
	kh1 = min( nbp1, mh - 1 );
	istep = *ir - *jtprev + 1;
	/*                                             ALG. STEP 20 */
	for (i = 1; i <= kh1; i++)
	{
		shtcc( 1, i, max( i + 1, istep ), mh, &G(i - 1,*jtprev - 1),
		 &uparam, &G(min( i + 1, nbp1 ) - 1,*jtprev - 1), ldg, nbp1 -
		 i );
	}
	/*                                             ALG. STEP 21 */
	kh = min( nbp1, mh );
	*ir = *jtprev + kh;
	/*                                             ALG. STEP 22 */
	if (kh >= nbp1)
	{
		/*                                             ALG. STEP 23 */
		for (i = 1; i <= nb; i++)
		{
			G(i - 1,*ir - 2) = ZERO;
		}
		/*                                             ALG. STEP 24 */
	}
	/*                                             ALG. STEP 25 */
	return;
#undef	G
} /* end of function */

void /*FUNCTION*/ strc2c(
double b[],
long nb,
double c1[],
double ytab[],
double yptab[])
{
	LOGICAL32 newseg;
	long int i, ib, ideriv, iseg, isegm1, mb, nseg;
	double p[4], sum;
		/* OFFSET Vectors w/subscript range: 1 to dimension */
	double *const B = &b[0] - 1;
	double *const C1 = &c1[0] - 1;
	double *const P = &p[0] - 1;
	double *const Yptab = &yptab[0] - 1;
	double *const Ytab = &ytab[0] - 1;
		/* end of OFFSET VECTORS */
 
	/* Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
	 * ALL RIGHTS RESERVED.
	 * Based on Government Sponsored Research NAS7-03001.
	 *>> 1994-11-11 STRC2C Krogh   Declared all vars.
	 *>> 1994-10-20 STRC2C Krogh  Changes to use M77CON
	 *>> 1987-12-09 STRC2C Lawson  Initial code.
	 *--S replaces "?": ?TRC2C, ?C2BAS
	 *          Given C1() containing coefficients relative to a
	 *     cubic B-spline basis, this subr computes and stores into
	 *     YTAB() and YPTAB()  respectively, the value and first derivative
	 *     of the spline curve at each nodal point given in B().
	 *     INPUT..
	 *     (B(I),I=1,NB)     PARTITION POINTS,INCLUDING ENDPOINTS
	 *     NB                NO. OF POINTS IN B().
	 *                       NO. OF SEGMENTS IS NSEG = NB-1
	 *     (C1(I),I=1,NB+2)   Coefficients relative to a cubic B-spline basis
	 *
	 *     OUTPUT..
	 *     (YTAB(i),i=1,NB) and (YPTAB(i),i=1,NB)
	 *                    Values of spline function and its first
	 *                    derivative at the points given in B().
	 *                    YTAB(2*j-1) = f(B(j),   YPTAB(2*j) = fprime(B(j)).
	 *     ------------------------------------------------------------------
	 *     Based on code written by Lawson, Hanson, Lang, and Campbell, JPL,
	 *     1968-1974.
	 *     Modified July 1984 and July 1987 by C. L. Lawson, JPL, for
	 *     inclusion in the MATH 77 library.
	 *     ------------------------------------------------------------------
	 *     Subprogram referenced:  SC2BAS
	 *     ------------------------------------------------------------------ */
	/*     ------------------------------------------------------------------ */
	mb = nb;
	nseg = mb - 1;
	for (ib = 1; ib <= mb; ib++)
	{
		iseg = min( nseg, ib );
		isegm1 = iseg - 1;
		newseg = TRUE;
		for (ideriv = 0; ideriv <= 1; ideriv++)
		{
			sc2bas( B[ib], ideriv, iseg, &newseg, b, nb, p );
			sum = ZERO;
			for (i = 1; i <= 4; i++)
			{
				sum += C1[isegm1 + i]*P[i];
			}
			if (ideriv == 0)
			{
				Ytab[ib] = sum;
			}
			else
			{
				Yptab[ib] = sum;
			}
		}
	}
	return;
} /* end of function */
 
void /*FUNCTION*/ sc2fit(
double xi[],
double yi[],
double sdi[],
long nxy,
double b[],
long nb,
double *w,
long nw,
double yknot[],
double ypknot[],
double *sigfac,
long *ierr1)
{
#define W(I_,J_)	(*(w+(I_)*(nw)+(J_)))
	LOGICAL32 newseg, usewt1;
	long int i, ierr2=0, ierr3, irnow, iseg, j, jpoint, jtprev, k, ksize,
	 n1, nparam;
	double dof, p[4], rnorm, sdijp, wt, wt1;
		/* OFFSET Vectors w/subscript range: 1 to dimension */
	double *const B = &b[0] - 1;
	double *const P = &p[0] - 1;
	double *const Sdi = &sdi[0] - 1;
	double *const Xi = &xi[0] - 1;
	double *const Yi = &yi[0] - 1;
	double *const Yknot = &yknot[0] - 1;
	double *const Ypknot = &ypknot[0] - 1;
		/* end of OFFSET VECTORS */
 
	/* Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
	 * ALL RIGHTS RESERVED.
	 * Based on Government Sponsored Research NAS7-03001.
	 *>> 2000-12-01 SC2FIT Krogh  Dim. SDI(*) instead NXY.
	 *>> 1995-11-21 SC2FIT Krogh  Converted from SFTRAN to Fortran 77.
	 *>> 1994-10-19 SC2FIT Krogh  Changes to use M77CON
	 *>> 1994-01-31 SC2FIT CLL Added test for SDI(i) .le. 0 when SDI(1) > 0.
	 *>> 1990-01-23 CLL Deleted ref to unused variable NX in call to IERM1
	 *>> 1989-10-20 CLL
	 *>> 1987-10-22 SC2FIT Lawson  Initial code.
	 *           Least squares fit to discrete data by a C-2 cubic spline.
	 *     ------------------------------------------------------------------
	 *     Algorithm and program designed by C.L.Lawson and R.J.Hanson.
	 *     The general approach but not the complete code is given in
	 *     'SOLVING LEAST SQUARES PROBLEMS', by Lawson and Hanson,
	 *     publ by Prentice-Hall, 1974.
	 *     Programming and later changes and corrections by Lawson,Hanson,
	 *     T.Lang, and D.Campbell, Sept 1968, Nov 1969, and Aug 1970.
	 *     Modified 1968 Sept 17 to provide C-2 continuity.
	 *     1974 5/21, C.L.Lawson, Fixed bug that caused ISEG to get too big.
	 *     Also changed to exit immidiately if B() array is not
	 *     strictly increasing.
	 *     1984 July 10. Modified for improved portability and to conform
	 *          to Fortran 77.  C. L. Lawson, JPL.
	 *          Added calls to the error message subrs.
	 *     7/23/87 CLL.  Added the IERR1 argument.
	 *     1989-10-20 CLL  Changed code so there are no RETURN statements in
	 *     Sftran procedures.  Previous code had such a RETURN that led to
	 *     a warning diagnostic from a Cray compiler due to unreachable
	 *     CONTINUE statements.  Also introduced "c--" lines for CHGTYP.
	 *     ------------------------------------------------------------------
	 *                     SUBROUTINE ARGUMENTS
	 *
	 *  (XI(i),i=1,NXY)  [in]  Abcissas of data to be fitted.  Require
	 *        this data be ordered so X(i) .le. X(i+1).
	 *
	 *  (YI(i),i=1,NXY)  [in]  Ordinates of data to be fitted.
	 *
	 *  (SDI(i),i=1,NXY) [in]  User may use this array to assign an
	 *                     a priori standard deviation of error to each
	 *       YI(i) value.  The weighted fitting algorithm will take
	 *       account of these.  Optionally the user may set SDI(1) to
	 *       a negative value.  Then this subr will use ABS(SDI(1)) as
	 *       the standard deviation for each YI(i) value.  In this case
	 *       the SDI() array can be dimensioned SDI(1).
	 *       If SDI(1) = 0., the subr issues an error message and returns.
	 *
	 *  NXY  [in]  No. of data pairs, (XI(i), YI(i)), and no. of elts
	 *          in SDI() if SDI(1) is positive.  Require NXY .ge. 4.
	 *
	 *  (B(j),j=1,NB)  [in]  Breakpoints for the spline function,
	 *        including endpoints.  These breakpoints must be
	 *        strictly increasing:  B(j) .lt. B(j+1).
	 *        It is required that all abcissas, XI(i), lie in the
	 *        closed interval, [B(1), B(NB)].
	 *
	 *  NB  [in]  No. of breakpoints, including endpoints.
	 *         The no. of parameters in the least squares problem
	 *         will be NB + 2.
	 *         To have a nonsingular problem one must have
	 *         NXY .ge. NB + 2, and the distribution of the breakpoints
	 *         must not be too skewed relative to the data abcissas.
	 *       If singularity is detected, an error message will be
	 *       issued by the subr that solves the band matrix.
	 *
	 *  W()  [scratch]  Work space dimensioned W(NW,5).
	 *
	 *  NW  [in]  First dimension of W().   Must satisfy NW .ge. NB + 4
	 *          Let KMAX denote the max no. of data abcissas, XI(i),
	 *          in any one breakpoint interval, i.e. between B(j) and
	 *          B(j+1) for some j.   The subr will be more efficient
	 *          if NW is at least NB + 3 + KMAX.
	 *
	 *  (YKNOT(k) and YPKNOT(k),k=1,NB)  [out]  The subr will return
	 *       values defining the fitted C2 spline curve in these arrays.
	 *       These values and first derivatives of the fitted curve at the
	 *       knot abcissae.  YKNOT(j) = f(B(j))   and
	 *       YPKNOT(j) = fprime(B(j))  for  j = 1,...,NB.
	 *       The user can then evaluate the fitted curve at any point by
	 *       Hermite interpolation.  See subrs DHINT or SHINT.
	 *
	 *  SIGFAC  [out]  The subr sets SIGFAC to RNORM / sqrt(DOF) where
	 *      RNORM = sqrt( sum over i of [( (yfit(i) - YI(i))/SDI(i))**2])
	 *       and DOF = max(1, NXY - (NB+2))
	 *
	 *  IERR1  [out]  Error status indicator.  Note that IERR2 comes from
	 *         SBACC and IERR3 comes from SBSOL.
	 * 							0 	// ALL OK
	 * SPLNFIT_NBINV 			1	// nb < 2
	 * SPLNFIT_INSUFFDATA 		2	// nxy < nb + 2
	 * SPLNFIT_KNOTSSTRICINC	3	// knot x coordinates not strictly increasing
	 * SPLNFIT_INSUFFWS			4	// insufficient workspace allocated
	 * SPLNFIT_SUNEQZERO		5	// S[0] must not be 0
	 * SPLNFIT_STRICTINC		6	// Xi's not strictly increasing
	 * SPLNFIT_XIOUTRANGE		7	// Xi's not within knot range
	 * SPLNFIT_ICRFWS			8	// increase workspace
	 * SPLNFIT_SBACC			9	// sbacc fail
	 * SPLNFIT_SNOTGTZERO		10	// require S[1+]>0 if S[0]>0
	 * SPLNFIT_SINGULAR			11	// Singular system
	 *     ------------------------------------------------------------------
	 *             Important internal variables.
	 *
	 *     ISEG     Index of current spline segment, starting with 1 for
	 *              the first segment.
	 *              Also tells the band matrix subroutine the column index
	 *              of the least squares matrix with which the first col
	 *              of the new block of data in G() is to be associated.
	 *     KSIZE    Size of current block.
	 *     JPOINT   Current data pointer.
	 *     ------------------------------------------------------------------
	 *--S replaces "?": ?C2FIT, ?BACC, ?C2BAS, ?ERM1, ?ERV1, ?BSOL, ?TRC2C
	 *     Both versions use ERMSG, IERM1, IERV1
	 *     Lower level subrs needed: (D/S)HTCC, (D/S)NRM2, ERFIN
	 *     ------------------------------------------------------------------ */
	/*     ------------------------------------------------------------------ */
	*ierr1 = 0;
 
	/*          EXIT IMMEDIATELY IF NB .lt. 2  OR  NXY .LT NB+2  OR  IF THE
	 *          BREAKPOINTS ARE NOT STRICTLY INCREASING.
	 * */
	nparam = nb + 2;
	if (nb < 2)
	{
		*ierr1=SPLNFIT_NBINV;
		goto L_300;
	}
	if (nxy < nparam)
	{
		*ierr1 = SPLNFIT_INSUFFDATA;
		goto L_300;
	}
 
	n1 = nb - 1;
	for (i = 1; i <= n1; i++)
	{
		if (B[i] >= B[i + 1])
		{
			*ierr1 = SPLNFIT_KNOTSSTRICINC;
			goto L_300;
		}
	}
 
	/*     Require NW .ge. NB+4
	 * */
	if (nw < nb + 4)
	{
		*ierr1 = SPLNFIT_INSUFFWS; // NW >= NB+4
		goto L_300;
	}
 
	/*     ------------------------------------------------------------------
	 *                                       TEST SDI(1) */
	if (Sdi[1] < ZERO)
	{
		wt1 = -ONE/Sdi[1];
		usewt1 = TRUE;
	}
	else if (Sdi[1] > ZERO)
	{
		usewt1 = FALSE;
	}
	else
	{
		*ierr1 = SPLNFIT_SUNEQZERO;
		goto L_300; // B. Pieters: was a return, but I want no chance knots are undefined
	}
 
	/*                             Test ordering of XI() array.
	 * */
	for (i = 2; i <= nxy; i++)
	{
		if (Xi[i - 1] > Xi[i])
		{
			*ierr1 = SPLNFIT_STRICTINC;
			goto L_300;
		}
	}
 
	/*                             TEST THE FIRST AND LAST BREAKPOINT
	 *                             FOR BRACKETING THE DATA ABCISSAS.
	 * */
	if (B[1] > Xi[1] || B[nb] < Xi[nxy])
	{
		*ierr1 = SPLNFIT_XIOUTRANGE;
		goto L_300;
	}
 
	/*     BEGIN LOOP TO FORM EQUATIONS FOR C2 LEAST SQUARES FIT.
	 * */
	irnow = 1;
	k = 1;
	ksize = 0;
	newseg = TRUE;
	iseg = 1;
	for (jpoint = 1; jpoint <= nxy; jpoint++)
	{
		if (k > nw)
		{
			sbacc( w, nw, NBAND, &irnow, ksize, iseg, &jtprev, &ierr2 );
			if (ierr2 != 0)
			{
				*ierr1 = SPLNFIT_SBACC;
				goto L_200;
			}
 
			if (irnow > nw)
			{
				*ierr1 = SPLNFIT_ICRFWS;
				goto L_300;
			}
			k = irnow;
			ksize = 0;
		}
 
		/*        DO WHILE( XI(JPOINT) .gt. B(ISEG+1) ) */
L_80:
		if (Xi[jpoint] > B[iseg + 1])
		{
			sbacc( w, nw, NBAND, &irnow, ksize, iseg, &jtprev, &ierr2 );
			if (ierr2 != 0)
			{
				*ierr1 = SPLNFIT_SBACC;
				goto L_200;
			}
			ksize = 0;
			k = irnow;
			iseg += 1;
			newseg = TRUE;
			goto L_80;
		}
		/*        END WHILE
		 *                          Build one equation */
		sc2bas( Xi[jpoint], 0, iseg, &newseg, b, nb, p );
		if (usewt1)
		{
			wt = wt1;
		}
		else
		{
			sdijp = Sdi[jpoint];
			if (sdijp > ZERO)
			{
				wt = ONE/sdijp;
			}
			else
			{
				*ierr1 = SPLNFIT_SNOTGTZERO;
				goto L_300;
			}
		}
		for (j = 1; j <= 4; j++)
		{
			W(j - 1,k - 1) = P[j]*wt;
		}
		W(4,k - 1) = Yi[jpoint]*wt;
		/*                          End of build one equation */
		k += 1;
		ksize += 1;
 
	}
	sbacc( w, nw, NBAND, &irnow, ksize, iseg, &jtprev, &ierr2 );
	if (ierr2 != 0)
	{
		*ierr1 = 900 + ierr2;
		goto L_200;
	}
 
	/*     ALL DATA POINTS HAVE BEEN PROCESSED.  CALL FOR SOLUTION.
	 * */
	sbsol( 1, w, nw, NBAND, irnow, jtprev, &W(NBAND1 - 1,0), nparam,
	 &rnorm, &ierr3 );
	if (ierr3 != 0)
	{
		*ierr1 = SPLNFIT_SINGULAR;
		goto L_300;
	}
 
	dof = max( 1, nxy - nparam );
	*sigfac = rnorm/sqrtf( dof );
 
	/*                  TRANSFORM PARAMETERS TO Y,YPRIME BASIS
	 * */
	strc2c( b, nb, &W(4,0), yknot, ypknot );
	return;
 
	/*             Error in _BACC */
L_200:
	//fprintf(stderr,"Error: sc2fit error in SBACC\n");
	/*             Set YKNOT() & YPKNOT() to zero */
L_300:
	for (i = 1; i <= nb; i++)
	{
		Yknot[i] = ZERO;
		Ypknot[i] = ZERO;
	}
	return;
#undef	W
} /* end of function */
 
