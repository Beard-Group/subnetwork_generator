#include "../common.h"
using namespace std;
# include "./spline.hpp"
//****************************************************************************80

double pchst ( double arg1, double arg2 )

//****************************************************************************80
//
//  Purpose:
//
//    PCHST: PCHIP sign-testing routine.
//
//  Discussion:
//
//    This routine essentially computes the sign of ARG1 * ARG2.
//
//    The object is to do this without multiplying ARG1 * ARG2, to avoid
//    possible over/underflow problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 August 2005
//
//  Author:
//
//    Original FORTRAN77 version by Fred Fritsch, Lawrence Livermore National Laboratory.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Fred Fritsch, Ralph Carlson, 
//    Monotone Piecewise Cubic Interpolation,
//    SIAM Journal on Numerical Analysis,
//    Volume 17, Number 2, April 1980, pages 238-246.
//
//  Parameters:
//
//    Input, double ARG1, ARG2, two values to check.
//
//    Output, double PCHST,
//    -1.0, if ARG1 and ARG2 are of opposite sign.
//     0.0, if either argument is zero.
//    +1.0, if ARG1 and ARG2 are of the same sign.
//
{
  double value;

  if ( arg1 == 0.0 )
  {
    value = 0.0;
  }
  else if ( arg1 < 0.0 )
  {
    if ( arg2 < 0.0 )
    {
      value = 1.0;
    }
    else if ( arg2 == 0.0 )
    {
      value = 0.0;
    }
    else if ( 0.0 < arg2 )
    {
      value = -1.0;
    }
  }
  else if ( 0.0 < arg1 )
  {
    if ( arg2 < 0.0 )
    {
      value = -1.0;
    }
    else if ( arg2 == 0.0 )
    {
      value = 0.0;
    }
    else if ( 0.0 < arg2 )
    {
      value = 1.0;
    }
  }

  return value;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  if ( y < x )
  {
    return x;
  } 
  else
  {
    return y;
  }
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  if ( y < x )
  {
    return y;
  } 
  else
  {
    return x;
  }
}
//****************************************************************************80

void spline_pchip_set ( const int n, const vector <double> &x, const vector <double> &f, vector <double> &d )

//****************************************************************************80
//
//  Purpose:
//
//    SPLINE_PCHIP_SET sets derivatives for a piecewise cubic Hermite interpolant.
//
//  Discussion:
//
//    This routine computes what would normally be called a Hermite 
//    interpolant.  However, the user is only required to supply function
//    values, not derivative values as well.  This routine computes
//    "suitable" derivative values, so that the resulting Hermite interpolant
//    has desirable shape and monotonicity properties.
//
//    The interpolant will have an extremum at each point where
//    monotonicity switches direction.
//
//    The resulting piecewise cubic Hermite function may be evaluated
//    by SPLINE_PCHIP_VAL..
//
//    This routine was originally called "PCHIM".
//
//    An "abs" was corrected to a "fabs" on the report of Thomas Beutlich,
//    10 October 2012.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 August 2005
//
//  Author:
//
//    FORTRAN77 original version by Fred Fritsch, Lawrence Livermore National Laboratory.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Fred Fritsch, Ralph Carlson,
//    Monotone Piecewise Cubic Interpolation,
//    SIAM Journal on Numerical Analysis,
//    Volume 17, Number 2, April 1980, pages 238-246.
//
//    Fred Fritsch, Judy Butland,
//    A Method for Constructing Local Monotone Piecewise 
//    Cubic Interpolants,
//    SIAM Journal on Scientific and Statistical Computing,
//    Volume 5, Number 2, 1984, pages 300-304.
//
//  Parameters:
//
//    Input, int N, the number of data points.  N must be at least 2.
//
//    Input, double X[N], the strictly increasing independent
//    variable values.
//
//    Input, double F[N], dependent variable values to be interpolated.  This 
//    routine is designed for monotonic data, but it will work for any F-array.
//    It will force extrema at points where monotonicity switches direction.
//
//    Output, double D[N], the derivative values at the
//    data points.  If the data are monotonic, these values will determine
//    a monotone cubic Hermite function.  
//
{
  double del1;
  double del2;
  double dmax;
  double dmin;
  double drat1;
  double drat2;
  double dsave;
  double h1;
  double h2;
  double hsum;
  double hsumt3;
  int i;
  int ierr;
  int nless1;
  double temp;
  double w1;
  double w2;
//
//  Check the arguments.
//
  if ( n < 2 )
  {
    ierr = -1;
    cerr << "\n";
    cerr << "SPLINE_PCHIP_SET - Fatal error!\n";
    cerr << "  Number of data points less than 2.\n";
    exit ( ierr );
  }

  for ( i = 1; i < n; i++ )
  {
    if ( x[i] <= x[i-1] )
    {
      ierr = -3;
      cerr << "\n";
      cerr << "SPLINE_PCHIP_SET - Fatal error!\n";
      cerr << "  X array not strictly increasing.\n";
      exit ( ierr );
    }
  }

  ierr = 0;
  nless1 = n - 1;
  h1 = x[1] - x[0];
  del1 = ( f[1] - f[0] ) / h1;
  dsave = del1;
//
//  Special case N=2, use linear interpolation.
//
  if ( n == 2 )
  {
    d.push_back(del1);
    d.push_back(del1);
    return;
  }
//
//  Normal case, 3 <= N.
//
  h2 = x[2] - x[1];
  del2 = ( f[2] - f[1] ) / h2;
//
//  Set D(1) via non-centered three point formula, adjusted to be
//  shape preserving.
//
  hsum = h1 + h2;
  w1 = ( h1 + hsum ) / hsum;
  w2 = -h1 / hsum;
  d[0] = w1 * del1 + w2 * del2;

  if ( pchst ( d[0], del1 ) <= 0.0 )
  {
    d[0] = 0.0;
  }
//
//  Need do this check only if monotonicity switches.
//
  else if ( pchst ( del1, del2 ) < 0.0 )
  {
     dmax = 3.0 * del1;

     if ( fabs ( dmax ) < fabs ( d[0] ) )
     {
       d[0] = dmax;
     }

  }
//
//  Loop through interior points.
//
  for ( i = 2; i <= nless1; i++ )
  {
    if ( 2 < i )
    {
      h1 = h2;
      h2 = x[i] - x[i-1];
      hsum = h1 + h2;
      del1 = del2;
      del2 = ( f[i] - f[i-1] ) / h2;
    }
//
//  Set D(I)=0 unless data are strictly monotonic.
//
    d[i-1] = 0.0;

    temp = pchst ( del1, del2 );

    if ( temp < 0.0 )
    {
      ierr = ierr + 1;
      dsave = del2;
    }
//
//  Count number of changes in direction of monotonicity.
//
    else if ( temp == 0.0 )
    {
      if ( del2 != 0.0 )
      {
        if ( pchst ( dsave, del2 ) < 0.0 )
        {
          ierr = ierr + 1;
        }
        dsave = del2;
      }
    }
//
//  Use Brodlie modification of Butland formula.
//
    else
    {
      hsumt3 = 3.0 * hsum;
      w1 = ( hsum + h1 ) / hsumt3;
      w2 = ( hsum + h2 ) / hsumt3;
      dmax = r8_max ( fabs ( del1 ), fabs ( del2 ) );
      dmin = r8_min ( fabs ( del1 ), fabs ( del2 ) );
      drat1 = del1 / dmax;
      drat2 = del2 / dmax;
      d[i-1] = dmin / ( w1 * drat1 + w2 * drat2 );
    }
  }
//
//  Set D(N) via non-centered three point formula, adjusted to be
//  shape preserving.
//
  w1 = -h2 / hsum;
  w2 = ( h2 + hsum ) / hsum;
  d[n-1] = w1 * del1 + w2 * del2;

  if ( pchst ( d[n-1], del2 ) <= 0.0 )
  {
    d[n-1] = 0.0;
  }
  else if ( pchst ( del1, del2 ) < 0.0 )
  {
//
//  Need do this check only if monotonicity switches.
//
    dmax = 3.0 * del2;

    if ( fabs ( dmax ) < fabs ( d[n-1] ) )
    {
      d[n-1] = dmax;
    }

  }
  return;
}
//****************************************************************************80
void spline_deriv ( int gene_index, int N_time_s, vector <double> &t_d_s, vector <double> &x_d_s, vector <double> &deriv_s)
{
     const int start = gene_index*N_time_s;
     const int end = (gene_index + 1) * N_time_s;
     vector <double> x_d_temp ( x_d_s.begin () + start, x_d_s.begin () + end );
     //cout << "Size of x_d_temp: " << x_d_temp.size () << endl;
     vector <double> d_temp (N_time_s);
     spline_pchip_set ( N_time_s, t_d_s, x_d_temp, d_temp );
     //spline_pchip_set ( N_time_points, t_d, x_66, d_66 );
     for ( int i = 0; i < N_time_s; i++ )
     {
         deriv_s[start+i] = d_temp[i]; 
     }

     return;

}
//
void spline_coeff ( const int N_gene, const int N_time_s, const vector <double> &t_d, const vector <double> &x_d, const vector <double> &deriv_s, vector <double> &cub_coeff_spline)
{
   for ( int i = 0; i < N_gene; i++)
   {
       cub_coeff_spline[i*N_time_s*4] = x_d[i*N_time_s];
       // 
       for ( int j = 1; j < N_time_s ; j++ )
       {
           // Estimate the coefficient of the piecewise cubic spline
           double df = 0.0;
           double h = 0.0;
           double del1, del2;

           h = t_d[j] - t_d[j-1];
           df = (x_d[i*N_time_s+j] - x_d[i*N_time_s+j-1]) / h;
 
           // p(x) = c0 + (x-x_0)*(c1+(x-x_0)*(c2+(x-x_0)*c3))
           cub_coeff_spline[i*N_time_s*4+j*4+0] = x_d[i*N_time_s+j-1];
           cub_coeff_spline[i*N_time_s*4+j*4+1] = deriv_s[i*N_time_s+j-1];  
           del1 = (deriv_s[i*N_time_s+j-1] - df) / h;
           del2 = (deriv_s[i*N_time_s+j] - df) / h;
           cub_coeff_spline[i*N_time_s*4+j*4+2] = - ( del1 + del1 + del2);
           cub_coeff_spline[i*N_time_s*4+j*4+3] = (del1 + del2) / h;
        }
   }
}
// Calculate mean and standard deviation
void calc_stats( const state_type &x_spline, state_type &mean, state_type &sd )
{
     double mean_temp = 0.0;
     mean_temp = (accumulate(x_spline.begin(),x_spline.end(),0.0)/double(x_spline.size()));
     mean.push_back(mean_temp);
     //Calculate standard deviation
     double sd_temp = 0.0;
     for ( int k = 0; k < x_spline.size(); k++)
     {
         sd_temp += pow((x_spline[k]-mean_temp),2.0);
     }
     sd_temp = sqrt(sd_temp/double(x_spline.size()));
     sd.push_back(sd_temp);
}
