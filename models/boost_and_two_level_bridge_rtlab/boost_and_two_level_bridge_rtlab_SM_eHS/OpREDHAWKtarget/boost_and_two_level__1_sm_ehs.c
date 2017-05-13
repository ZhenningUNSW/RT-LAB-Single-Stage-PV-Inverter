/*
 * boost_and_two_level__1_sm_ehs.c
 *
 * Code generation for model "boost_and_two_level__1_sm_ehs".
 *
 * Model version              : 1.1045
 * Simulink Coder version : 8.7 (R2014b) 08-Sep-2014
 * C source code generated on : Thu May 11 18:35:05 2017
 *
 * Target selection: rtlab_rtmodel.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "boost_and_two_level__1_sm_ehs.h"
#include "boost_and_two_level__1_sm_ehs_private.h"

const real_T boost_and_two_level__1_sm_ehs_RGND = 0.0;/* real_T ground */

/* Block signals (auto storage) */
B_boost_and_two_level__1_sm_ehs_T boost_and_two_level__1_sm_ehs_B;

/* Continuous states */
X_boost_and_two_level__1_sm_ehs_T boost_and_two_level__1_sm_ehs_X;

/* Block states (auto storage) */
DW_boost_and_two_level__1_sm_ehs_T boost_and_two_level__1_sm_ehs_DW;

/* Real-time model */
RT_MODEL_boost_and_two_level__1_sm_ehs_T boost_and_two_level__1_sm_ehs_M_;
RT_MODEL_boost_and_two_level__1_sm_ehs_T *const boost_and_two_level__1_sm_ehs_M =
  &boost_and_two_level__1_sm_ehs_M_;

/* Lookup Binary Search Utility BINARYSEARCH_real_T */
void BINARYSEARCH_real_T(uint32_T *piLeft, uint32_T *piRght, real_T u, const
  real_T *pData, uint32_T iHi)
{
  /* Find the location of current input value in the data table. */
  *piLeft = 0U;
  *piRght = iHi;
  if (u <= pData[0] ) {
    /* Less than or equal to the smallest point in the table. */
    *piRght = 0U;
  } else if (u >= pData[iHi] ) {
    /* Greater than or equal to the largest point in the table. */
    *piLeft = iHi;
  } else {
    uint32_T i;

    /* Do a binary search. */
    while (( *piRght - *piLeft ) > 1U ) {
      /* Get the average of the left and right indices using to Floor rounding. */
      i = (*piLeft + *piRght) >> 1;

      /* Move either the right index or the left index so that */
      /*  LeftDataPoint <= CurrentValue < RightDataPoint */
      if (u < pData[i] ) {
        *piRght = i;
      } else {
        *piLeft = i;
      }
    }
  }
}

/* Lookup Utility LookUp_real_T_real_T */
void LookUp_real_T_real_T(real_T *pY, const real_T *pYData, real_T u, const
  real_T *pUData, uint32_T iHi)
{
  uint32_T iLeft;
  uint32_T iRght;
  BINARYSEARCH_real_T( &(iLeft), &(iRght), u, pUData, iHi);

  {
    real_T lambda;
    if (pUData[iRght] > pUData[iLeft] ) {
      real_T num;
      real_T den;
      den = pUData[iRght];
      den = den - pUData[iLeft];
      num = u;
      num = num - pUData[iLeft];
      lambda = num / den;
    } else {
      lambda = 0.0;
    }

    {
      real_T yLeftCast;
      real_T yRghtCast;
      yLeftCast = pYData[iLeft];
      yRghtCast = pYData[iRght];
      yLeftCast += lambda * ( yRghtCast - yLeftCast );
      (*pY) = yLeftCast;
    }
  }
}

/*
 * Time delay interpolation routine
 *
 * The linear interpolation is performed using the formula:
 *
 *          (t2 - tMinusDelay)         (tMinusDelay - t1)
 * u(t)  =  ----------------- * u1  +  ------------------- * u2
 *              (t2 - t1)                  (t2 - t1)
 */
real_T rt_TDelayInterpolate(
  real_T tMinusDelay,                  /* tMinusDelay = currentSimTime - delay */
  real_T tStart,
  real_T *tBuf,
  real_T *uBuf,
  int_T bufSz,
  int_T *lastIdx,
  int_T oldestIdx,
  int_T newIdx,
  real_T initOutput,
  boolean_T discrete,
  boolean_T minorStepAndTAtLastMajorOutput)
{
  int_T i;
  real_T yout, t1, t2, u1, u2;

  /*
   * If there is only one data point in the buffer, this data point must be
   * the t= 0 and tMinusDelay > t0, it ask for something unknown. The best
   * guess if initial output as well
   */
  if ((newIdx == 0) && (oldestIdx ==0 ) && (tMinusDelay > tStart))
    return initOutput;

  /*
   * If tMinusDelay is less than zero, should output initial value
   */
  if (tMinusDelay <= tStart)
    return initOutput;

  /* For fixed buffer extrapolation:
   * if tMinusDelay is small than the time at oldestIdx, if discrete, output
   * tailptr value,  else use tailptr and tailptr+1 value to extrapolate
   * It is also for fixed buffer. Note: The same condition can happen for transport delay block where
   * use tStart and and t[tail] other than using t[tail] and t[tail+1].
   * See below
   */
  if ((tMinusDelay <= tBuf[oldestIdx] ) ) {
    if (discrete) {
      return(uBuf[oldestIdx]);
    } else {
      int_T tempIdx= oldestIdx + 1;
      if (oldestIdx == bufSz-1)
        tempIdx = 0;
      t1= tBuf[oldestIdx];
      t2= tBuf[tempIdx];
      u1= uBuf[oldestIdx];
      u2= uBuf[tempIdx];
      if (t2 == t1) {
        if (tMinusDelay >= t2) {
          yout = u2;
        } else {
          yout = u1;
        }
      } else {
        real_T f1 = (t2-tMinusDelay) / (t2-t1);
        real_T f2 = 1.0 - f1;

        /*
         * Use Lagrange's interpolation formula.  Exact outputs at t1, t2.
         */
        yout = f1*u1 + f2*u2;
      }

      return yout;
    }
  }

  /*
   * When block does not have direct feedthrough, we use the table of
   * values to extrapolate off the end of the table for delays that are less
   * than 0 (less then step size).  This is not completely accurate.  The
   * chain of events is as follows for a given time t.  Major output - look
   * in table.  Update - add entry to table.  Now, if we call the output at
   * time t again, there is a new entry in the table. For very small delays,
   * this means that we will have a different answer from the previous call
   * to the output fcn at the same time t.  The following code prevents this
   * from happening.
   */
  if (minorStepAndTAtLastMajorOutput) {
    /* pretend that the new entry has not been added to table */
    if (newIdx != 0) {
      if (*lastIdx == newIdx) {
        (*lastIdx)--;
      }

      newIdx--;
    } else {
      if (*lastIdx == newIdx) {
        *lastIdx = bufSz-1;
      }

      newIdx = bufSz - 1;
    }
  }

  i = *lastIdx;
  if (tBuf[i] < tMinusDelay) {
    /* Look forward starting at last index */
    while (tBuf[i] < tMinusDelay) {
      /* May occur if the delay is less than step-size - extrapolate */
      if (i == newIdx)
        break;
      i = ( i < (bufSz-1) ) ? (i+1) : 0;/* move through buffer */
    }
  } else {
    /*
     * Look backwards starting at last index which can happen when the
     * delay time increases.
     */
    while (tBuf[i] >= tMinusDelay) {
      /*
       * Due to the entry condition at top of function, we
       * should never hit the end.
       */
      i = (i > 0) ? i-1 : (bufSz-1);   /* move through buffer */
    }

    i = ( i < (bufSz-1) ) ? (i+1) : 0;
  }

  *lastIdx = i;
  if (discrete) {
    /*
     * tempEps = 128 * eps;
     * localEps = max(tempEps, tempEps*fabs(tBuf[i]))/2;
     */
    double tempEps = (DBL_EPSILON) * 128.0;
    double localEps = tempEps * fabs(tBuf[i]);
    if (tempEps > localEps) {
      localEps = tempEps;
    }

    localEps = localEps / 2.0;
    if (tMinusDelay >= (tBuf[i] - localEps)) {
      yout = uBuf[i];
    } else {
      if (i == 0) {
        yout = uBuf[bufSz-1];
      } else {
        yout = uBuf[i-1];
      }
    }
  } else {
    if (i == 0) {
      t1 = tBuf[bufSz-1];
      u1 = uBuf[bufSz-1];
    } else {
      t1 = tBuf[i-1];
      u1 = uBuf[i-1];
    }

    t2 = tBuf[i];
    u2 = uBuf[i];
    if (t2 == t1) {
      if (tMinusDelay >= t2) {
        yout = u2;
      } else {
        yout = u1;
      }
    } else {
      real_T f1 = (t2-tMinusDelay) / (t2-t1);
      real_T f2 = 1.0 - f1;

      /*
       * Use Lagrange's interpolation formula.  Exact outputs at t1, t2.
       */
      yout = f1*u1 + f2*u2;
    }
  }

  return(yout);
}

void mul_wide_u32(uint32_T in0, uint32_T in1, uint32_T *ptrOutBitsHi, uint32_T
                  *ptrOutBitsLo)
{
  uint32_T outBitsLo;
  uint32_T in0Lo;
  uint32_T in0Hi;
  uint32_T in1Lo;
  uint32_T in1Hi;
  uint32_T productHiLo;
  uint32_T productLoHi;
  in0Hi = in0 >> 16U;
  in0Lo = in0 & 65535U;
  in1Hi = in1 >> 16U;
  in1Lo = in1 & 65535U;
  productHiLo = in0Hi * in1Lo;
  productLoHi = in0Lo * in1Hi;
  in0Lo *= in1Lo;
  in1Lo = 0U;
  outBitsLo = (productLoHi << 16U) + in0Lo;
  if (outBitsLo < in0Lo) {
    in1Lo = 1U;
  }

  in0Lo = outBitsLo;
  outBitsLo += productHiLo << 16U;
  if (outBitsLo < in0Lo) {
    in1Lo++;
  }

  *ptrOutBitsHi = (((productLoHi >> 16U) + (productHiLo >> 16U)) + in0Hi * in1Hi)
    + in1Lo;
  *ptrOutBitsLo = outBitsLo;
}

uint32_T mul_u32_u32_u32_sr29(uint32_T a, uint32_T b)
{
  uint32_T result;
  uint32_T u32_chi;
  mul_wide_u32(a, b, &u32_chi, &result);
  return u32_chi << 3U | result >> 29U;
}

/*
 * This function updates continuous states using the ODE3 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static const real_T rt_ODE3_A[3] = {
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3] = {
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE3_IntgData *id = (ODE3_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T hB[3];
  int_T i;
  int_T nXc = 2;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  boost_and_two_level__1_sm_ehs_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  boost_and_two_level__1_sm_ehs_output();
  boost_and_two_level__1_sm_ehs_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  boost_and_two_level__1_sm_ehs_output();
  boost_and_two_level__1_sm_ehs_derivatives();

  /* tnew = t + hA(3);
     ynew = y + f*hB(:,3); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

real_T rt_roundd_snf(real_T u)
{
  real_T y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

real_T rt_remd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T tr;
  boolean_T y_0;
  boolean_T y_1;
  y_0 = ((!rtIsNaN(u0)) && (!rtIsInf(u0)));
  y_1 = ((!rtIsNaN(u1)) && (!rtIsInf(u1)));
  if (!(y_0 && y_1)) {
    y = (rtNaN);
  } else {
    if (u1 < 0.0) {
      y = ceil(u1);
    } else {
      y = floor(u1);
    }

    if ((u1 != 0.0) && (u1 != y)) {
      tr = u0 / u1;
      if (fabs(tr - rt_roundd_snf(tr)) <= DBL_EPSILON * fabs(tr)) {
        y = 0.0;
      } else {
        y = fmod(u0, u1);
      }
    } else {
      y = fmod(u0, u1);
    }
  }

  return y;
}

real_T rt_modd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T tmp;
  boolean_T y_0;
  boolean_T y_1;
  if (u1 == 0.0) {
    y = u0;
  } else {
    y_0 = ((!rtIsNaN(u0)) && (!rtIsInf(u0)));
    y_1 = ((!rtIsNaN(u1)) && (!rtIsInf(u1)));
    if (!(y_0 && y_1)) {
      y = (rtNaN);
    } else {
      tmp = u0 / u1;
      if (u1 <= floor(u1)) {
        y = u0 - floor(tmp) * u1;
      } else if (fabs(tmp - rt_roundd_snf(tmp)) <= DBL_EPSILON * fabs(tmp)) {
        y = 0.0;
      } else {
        y = (tmp - floor(tmp)) * u1;
      }
    }
  }

  return y;
}

real_T rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T a;
  a = fabs(u0);
  y = fabs(u1);
  if (a < y) {
    a /= y;
    y *= sqrt(a * a + 1.0);
  } else if (a > y) {
    y /= a;
    y = sqrt(y * y + 1.0) * a;
  } else {
    if (!rtIsNaN(y)) {
      y = a * 1.4142135623730951;
    }
  }

  return y;
}

real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  int32_T tmp;
  int32_T tmp_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u1 > 0.0) {
      tmp = 1;
    } else {
      tmp = -1;
    }

    if (u0 > 0.0) {
      tmp_0 = 1;
    } else {
      tmp_0 = -1;
    }

    y = atan2(tmp_0, tmp);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T tmp;
  real_T tmp_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    tmp = fabs(u0);
    tmp_0 = fabs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = (rtNaN);
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = (rtNaN);
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

/* Model output function */
void boost_and_two_level__1_sm_ehs_output(void)
{
  real_T P;
  real_T dV;
  real_T dP;
  real_T B;
  boolean_T tmp;
  int32_T i;
  real_T u2;
  uint32_T u0;
  uint32_T u1;
  uint32_T u2_0;
  if (rtmIsMajorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
    /* set solver stop time */
    if (!(boost_and_two_level__1_sm_ehs_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&boost_and_two_level__1_sm_ehs_M->solverInfo,
                            ((boost_and_two_level__1_sm_ehs_M->Timing.clockTickH0
        + 1) * boost_and_two_level__1_sm_ehs_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&boost_and_two_level__1_sm_ehs_M->solverInfo,
                            ((boost_and_two_level__1_sm_ehs_M->Timing.clockTick0
        + 1) * boost_and_two_level__1_sm_ehs_M->Timing.stepSize0 +
        boost_and_two_level__1_sm_ehs_M->Timing.clockTickH0 *
        boost_and_two_level__1_sm_ehs_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
    boost_and_two_level__1_sm_ehs_M->Timing.t[0] = rtsiGetT
      (&boost_and_two_level__1_sm_ehs_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
    /* Memory: '<S1>/S-Function' */
    boost_and_two_level__1_sm_ehs_B.SFunction =
      boost_and_two_level__1_sm_ehs_DW.SFunction_PreviousInput;

    /* Sum: '<S1>/Sum' incorporates:
     *  Constant: '<S1>/S-Function1'
     */
    boost_and_two_level__1_sm_ehs_B.Sum =
      boost_and_two_level__1_sm_ehs_P.SFunction1_Value +
      boost_and_two_level__1_sm_ehs_B.SFunction;

    /* Stop: '<S1>/Stop Simulation' */
    if (boost_and_two_level__1_sm_ehs_B.Sum != 0.0) {
      rtmSetStopRequested(boost_and_two_level__1_sm_ehs_M, 1);
    }

    /* End of Stop: '<S1>/Stop Simulation' */

    /* Level2 S-Function Block: '<S7>/Outputs_eHS1_Recv' (sfun_fct_op7160ex1_recv) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[2];
      sfcnOutputs(rts, 1);
    }

    /* Level2 S-Function Block: '<S7>/Convert  Single floating-point (FPGA)  to double' (sfun_SFP2DBL) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[3];
      sfcnOutputs(rts, 1);
    }

    /* Saturate: '<S7>/Saturation' */
    u0 = boost_and_two_level__1_sm_ehs_B.Outputs_eHS1_Recv_o1[7];
    u1 = boost_and_two_level__1_sm_ehs_P.Saturation_LowerSat_f;
    u2_0 = boost_and_two_level__1_sm_ehs_P.Saturation_UpperSat_k;
    if (u0 > u2_0) {
      boost_and_two_level__1_sm_ehs_B.Saturation_k = u2_0;
    } else if (u0 < u1) {
      boost_and_two_level__1_sm_ehs_B.Saturation_k = u1;
    } else {
      boost_and_two_level__1_sm_ehs_B.Saturation_k = u0;
    }

    /* End of Saturate: '<S7>/Saturation' */

    /* Product: '<S7>/Divide' */
    for (i = 0; i < 7; i++) {
      boost_and_two_level__1_sm_ehs_B.Divide[i] =
        boost_and_two_level__1_sm_ehs_B.ConvertSinglefloatingpointFPGAt[i] /
        (real_T)boost_and_two_level__1_sm_ehs_B.Saturation_k;
    }

    /* End of Product: '<S7>/Divide' */

    /* UnitDelay: '<S31>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_B.UnitDelay =
      boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE;

    /* DigitalClock: '<S109>/Digital Clock' */
    boost_and_two_level__1_sm_ehs_B.DigitalClock =
      boost_and_two_level__1_sm_ehs_M->Timing.t[1];

    /* Sum: '<S109>/Add1' incorporates:
     *  Constant: '<S109>/Constant3'
     */
    boost_and_two_level__1_sm_ehs_B.Add1 =
      boost_and_two_level__1_sm_ehs_B.DigitalClock +
      boost_and_two_level__1_sm_ehs_P.Constant3_Value;

    /* Math: '<S109>/Math Function' incorporates:
     *  Constant: '<S109>/Constant1'
     */
    boost_and_two_level__1_sm_ehs_B.MathFunction = rt_remd_snf
      (boost_and_two_level__1_sm_ehs_B.Add1,
       boost_and_two_level__1_sm_ehs_P.Constant1_Value_l);

    /* Gain: '<S109>/1\ib1' */
    boost_and_two_level__1_sm_ehs_B.ib1 =
      boost_and_two_level__1_sm_ehs_P.ib1_Gain *
      boost_and_two_level__1_sm_ehs_B.MathFunction;

    /* Lookup: '<S109>/Lookup Table'
     * About '<S109>/Lookup Table':
     * Input0  Data Type:  Floating Point real_T
     * Output0 Data Type:  Floating Point real_T
     * Lookup Method: Linear_Endpoint
     *
     * XData parameter uses the same data type and scaling as Input0
     * YData parameter uses the same data type and scaling as Output0
     */
    LookUp_real_T_real_T( &(boost_and_two_level__1_sm_ehs_B.LookupTable),
                         boost_and_two_level__1_sm_ehs_P.LookupTable_YData,
                         boost_and_two_level__1_sm_ehs_B.ib1,
                         boost_and_two_level__1_sm_ehs_P.LookupTable_XData, 2U);

    /* Sum: '<S109>/Add3' incorporates:
     *  Constant: '<S109>/Constant2'
     */
    boost_and_two_level__1_sm_ehs_B.Add3 =
      boost_and_two_level__1_sm_ehs_B.LookupTable -
      boost_and_two_level__1_sm_ehs_P.Constant2_Value_f;

    /* Sum: '<S87>/Add3' incorporates:
     *  Constant: '<S37>/Constant10'
     */
    boost_and_two_level__1_sm_ehs_B.Add3_g =
      boost_and_two_level__1_sm_ehs_P.PWM_Generator_MinMax[1] -
      boost_and_two_level__1_sm_ehs_P.PWM_Generator_MinMax[0];

    /* Gain: '<S87>/Gain1' */
    boost_and_two_level__1_sm_ehs_B.Gain1 =
      boost_and_two_level__1_sm_ehs_P.Gain1_Gain_a *
      boost_and_two_level__1_sm_ehs_B.Add3_g;

    /* Product: '<S87>/MUL1' */
    boost_and_two_level__1_sm_ehs_B.MUL1 = boost_and_two_level__1_sm_ehs_B.Add3 *
      boost_and_two_level__1_sm_ehs_B.Gain1;

    /* Sum: '<S87>/Add4' incorporates:
     *  Constant: '<S37>/Constant10'
     */
    boost_and_two_level__1_sm_ehs_B.Add4 =
      (boost_and_two_level__1_sm_ehs_P.PWM_Generator_MinMax[0] +
       boost_and_two_level__1_sm_ehs_B.MUL1) +
      boost_and_two_level__1_sm_ehs_B.Gain1;

    /* RelationalOperator: '<S91>/Relational Operator2' */
    boost_and_two_level__1_sm_ehs_B.RelationalOperator2 =
      (boost_and_two_level__1_sm_ehs_B.UnitDelay >=
       boost_and_two_level__1_sm_ehs_B.Add4);

    /* Logic: '<S91>/Logical Operator' */
    boost_and_two_level__1_sm_ehs_B.LogicalOperator =
      !boost_and_two_level__1_sm_ehs_B.RelationalOperator2;

    /* Logic: '<S37>/Logical Operator4' */
    boost_and_two_level__1_sm_ehs_B.LogicalOperator4[0] =
      !boost_and_two_level__1_sm_ehs_B.RelationalOperator2;
    boost_and_two_level__1_sm_ehs_B.LogicalOperator4[1] =
      !boost_and_two_level__1_sm_ehs_B.LogicalOperator;

    /* DataTypeConversion: '<S37>/Data Type Conversion' */
    boost_and_two_level__1_sm_ehs_B.DataTypeConversion[0] =
      boost_and_two_level__1_sm_ehs_B.RelationalOperator2;
    boost_and_two_level__1_sm_ehs_B.DataTypeConversion[1] =
      boost_and_two_level__1_sm_ehs_B.LogicalOperator4[0];
    boost_and_two_level__1_sm_ehs_B.DataTypeConversion[2] =
      boost_and_two_level__1_sm_ehs_B.LogicalOperator;
    boost_and_two_level__1_sm_ehs_B.DataTypeConversion[3] =
      boost_and_two_level__1_sm_ehs_B.LogicalOperator4[1];

    /* Level2 S-Function Block: '<S10>/S-Function' (RECV_Param) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[4];
      sfcnOutputs(rts, 1);
    }

    /* Relay: '<S7>/Relay1' */
    if (boost_and_two_level__1_sm_ehs_B.SFunction_l[3] >=
        boost_and_two_level__1_sm_ehs_P.Relay1_OnVal) {
      boost_and_two_level__1_sm_ehs_DW.Relay1_Mode = true;
    } else {
      if (boost_and_two_level__1_sm_ehs_B.SFunction_l[3] <=
          boost_and_two_level__1_sm_ehs_P.Relay1_OffVal) {
        boost_and_two_level__1_sm_ehs_DW.Relay1_Mode = false;
      }
    }

    if (boost_and_two_level__1_sm_ehs_DW.Relay1_Mode) {
      boost_and_two_level__1_sm_ehs_B.Relay1 =
        boost_and_two_level__1_sm_ehs_P.Relay1_YOn;
    } else {
      boost_and_two_level__1_sm_ehs_B.Relay1 =
        boost_and_two_level__1_sm_ehs_P.Relay1_YOff;
    }

    /* End of Relay: '<S7>/Relay1' */

    /* DataTypeConversion: '<S7>/Data Type Conversion1' */
    B = floor(boost_and_two_level__1_sm_ehs_B.Relay1);
    if (rtIsNaN(B) || rtIsInf(B)) {
      B = 0.0;
    } else {
      B = fmod(B, 4.294967296E+9);
    }

    boost_and_two_level__1_sm_ehs_B.DataTypeConversion1_h = B < 0.0 ? (uint32_T)
      -(int32_T)(uint32_T)-B : (uint32_T)B;

    /* End of DataTypeConversion: '<S7>/Data Type Conversion1' */

    /* Saturate: '<S7>/sat_scn' */
    u1 = boost_and_two_level__1_sm_ehs_P.sat_scn_LowerSat;
    if (0U < u1) {
      boost_and_two_level__1_sm_ehs_B.sat_scn = u1;
    } else {
      boost_and_two_level__1_sm_ehs_B.sat_scn = 0U;
    }

    /* End of Saturate: '<S7>/sat_scn' */

    /* Gain: '<S7>/shift_2bits' */
    boost_and_two_level__1_sm_ehs_B.shift_2bits = mul_u32_u32_u32_sr29
      (boost_and_two_level__1_sm_ehs_P.shift_2bits_Gain,
       boost_and_two_level__1_sm_ehs_B.sat_scn);

    /* Sum: '<S7>/Add' */
    boost_and_two_level__1_sm_ehs_B.Add_l =
      boost_and_two_level__1_sm_ehs_B.DataTypeConversion1_h +
      boost_and_two_level__1_sm_ehs_B.shift_2bits;

    /* UnitDelay: '<S25>/Delay Input1' */
    boost_and_two_level__1_sm_ehs_B.Uk1 =
      boost_and_two_level__1_sm_ehs_DW.DelayInput1_DSTATE;

    /* RelationalOperator: '<S25>/FixPt Relational Operator' */
    boost_and_two_level__1_sm_ehs_B.FixPtRelationalOperator = (uint8_T)
      (boost_and_two_level__1_sm_ehs_B.Add_l !=
       boost_and_two_level__1_sm_ehs_B.Uk1);

    /* Logic: '<S7>/Logical Operator' */
    boost_and_two_level__1_sm_ehs_B.LogicalOperator_c =
      (boost_and_two_level__1_sm_ehs_B.FixPtRelationalOperator != 0);

    /* DataTypeConversion: '<S7>/Data Type Conversion8' */
    boost_and_two_level__1_sm_ehs_B.DataTypeConversion8 =
      boost_and_two_level__1_sm_ehs_B.LogicalOperator_c;

    /* Level2 S-Function Block: '<S7>/eHS_rst_loadin' (sfun_fct_op7160ex1_load_in) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[5];
      sfcnOutputs(rts, 1);
    }

    /* Level2 S-Function Block: '<S7>/Automated_Solver_Mat_Initialisation_1' (sfun_efs_solver_cfg) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[6];
      sfcnOutputs(rts, 1);
    }
  }

  /* Sin: '<S9>/Sine Wave Function' */
  boost_and_two_level__1_sm_ehs_B.SineWaveFunction = sin
    (boost_and_two_level__1_sm_ehs_P.SineWaveFunction_Freq *
     boost_and_two_level__1_sm_ehs_M->Timing.t[0] +
     boost_and_two_level__1_sm_ehs_P.SineWaveFunction_Phase) *
    boost_and_two_level__1_sm_ehs_P.SineWaveFunction_Amp +
    boost_and_two_level__1_sm_ehs_P.SineWaveFunction_Bias;
  if (rtmIsMajorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
    /* Memory: '<S9>/Memory1' */
    boost_and_two_level__1_sm_ehs_B.Memory1 =
      boost_and_two_level__1_sm_ehs_DW.Memory1_PreviousInput;

    /* Memory: '<S9>/Memory' */
    boost_and_two_level__1_sm_ehs_B.Memory =
      boost_and_two_level__1_sm_ehs_DW.Memory_PreviousInput;
  }

  /* Level2 S-Function Block: '<S7>/Convert double to  Single floating-point (FPGA)' (sfun_DBL2SFP) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[7];
    sfcnOutputs(rts, 0);
  }

  if (rtmIsMajorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
    /* Level2 S-Function Block: '<S7>/Inputs_eHS1_Send' (sfun_fct_op7160ex1_send) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[8];
      sfcnOutputs(rts, 1);
    }
  }

  /* Integrator: '<S15>/integrator' */
  boost_and_two_level__1_sm_ehs_B.integrator =
    boost_and_two_level__1_sm_ehs_X.integrator_CSTATE;

  /* TransportDelay: '<S15>/Transport Delay' */
  {
    real_T **uBuffer = (real_T**)
      &boost_and_two_level__1_sm_ehs_DW.TransportDelay_PWORK.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)
      &boost_and_two_level__1_sm_ehs_DW.TransportDelay_PWORK.TUbufferPtrs[1];
    real_T simTime = boost_and_two_level__1_sm_ehs_M->Timing.t[0];
    real_T tMinusDelay = simTime -
      (boost_and_two_level__1_sm_ehs_P.TransportDelay_Delay);
    boost_and_two_level__1_sm_ehs_B.TransportDelay = rt_TDelayInterpolate(
      tMinusDelay,
      0.0,
      *tBuffer,
      *uBuffer,
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.CircularBufSize,
      &boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.Last,
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.Tail,
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.Head,
      boost_and_two_level__1_sm_ehs_P.TransportDelay_InitOutput,
      0,
      0);
  }

  /* Clock: '<S15>/Clock' */
  boost_and_two_level__1_sm_ehs_B.Clock =
    boost_and_two_level__1_sm_ehs_M->Timing.t[0];

  /* RelationalOperator: '<S15>/Relational Operator' incorporates:
   *  Constant: '<S15>/K1'
   */
  boost_and_two_level__1_sm_ehs_B.RelationalOperator =
    (boost_and_two_level__1_sm_ehs_B.Clock >=
     boost_and_two_level__1_sm_ehs_P.K1_Value);
  if (rtmIsMajorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
    /* Memory: '<S15>/Memory' */
    boost_and_two_level__1_sm_ehs_B.Memory_g =
      boost_and_two_level__1_sm_ehs_DW.Memory_PreviousInput_n;
  }

  /* Switch: '<S15>/Switch' */
  if (boost_and_two_level__1_sm_ehs_B.RelationalOperator) {
    /* Sum: '<S15>/Sum' */
    boost_and_two_level__1_sm_ehs_B.Sum_h =
      boost_and_two_level__1_sm_ehs_B.integrator -
      boost_and_two_level__1_sm_ehs_B.TransportDelay;

    /* Gain: '<S15>/Gain' */
    boost_and_two_level__1_sm_ehs_B.Gain_a =
      boost_and_two_level__1_sm_ehs_P.Gain_Gain *
      boost_and_two_level__1_sm_ehs_B.Sum_h;
    boost_and_two_level__1_sm_ehs_B.Switch =
      boost_and_two_level__1_sm_ehs_B.Gain_a;
  } else {
    boost_and_two_level__1_sm_ehs_B.Switch =
      boost_and_two_level__1_sm_ehs_B.Memory_g;
  }

  /* End of Switch: '<S15>/Switch' */

  /* Integrator: '<S16>/integrator' */
  boost_and_two_level__1_sm_ehs_B.integrator_p =
    boost_and_two_level__1_sm_ehs_X.integrator_CSTATE_b;

  /* TransportDelay: '<S16>/Transport Delay' */
  {
    real_T **uBuffer = (real_T**)
      &boost_and_two_level__1_sm_ehs_DW.TransportDelay_PWORK_d.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)
      &boost_and_two_level__1_sm_ehs_DW.TransportDelay_PWORK_d.TUbufferPtrs[1];
    real_T simTime = boost_and_two_level__1_sm_ehs_M->Timing.t[0];
    real_T tMinusDelay = simTime -
      (boost_and_two_level__1_sm_ehs_P.TransportDelay_Delay_b);
    boost_and_two_level__1_sm_ehs_B.TransportDelay_p = rt_TDelayInterpolate(
      tMinusDelay,
      0.0,
      *tBuffer,
      *uBuffer,
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.CircularBufSize,
      &boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.Last,
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.Tail,
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.Head,
      boost_and_two_level__1_sm_ehs_P.TransportDelay_InitOutput_f,
      0,
      0);
  }

  /* Clock: '<S16>/Clock' */
  boost_and_two_level__1_sm_ehs_B.Clock_m =
    boost_and_two_level__1_sm_ehs_M->Timing.t[0];

  /* RelationalOperator: '<S16>/Relational Operator' incorporates:
   *  Constant: '<S16>/K1'
   */
  boost_and_two_level__1_sm_ehs_B.RelationalOperator_b =
    (boost_and_two_level__1_sm_ehs_B.Clock_m >=
     boost_and_two_level__1_sm_ehs_P.K1_Value_k);
  if (rtmIsMajorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
    /* Memory: '<S16>/Memory' */
    boost_and_two_level__1_sm_ehs_B.Memory_e =
      boost_and_two_level__1_sm_ehs_DW.Memory_PreviousInput_e;
  }

  /* Switch: '<S16>/Switch' */
  if (boost_and_two_level__1_sm_ehs_B.RelationalOperator_b) {
    /* Sum: '<S16>/Sum' */
    boost_and_two_level__1_sm_ehs_B.Sum_ns =
      boost_and_two_level__1_sm_ehs_B.integrator_p -
      boost_and_two_level__1_sm_ehs_B.TransportDelay_p;

    /* Gain: '<S16>/Gain' */
    boost_and_two_level__1_sm_ehs_B.Gain_k =
      boost_and_two_level__1_sm_ehs_P.Gain_Gain_l *
      boost_and_two_level__1_sm_ehs_B.Sum_ns;
    boost_and_two_level__1_sm_ehs_B.Switch_e =
      boost_and_two_level__1_sm_ehs_B.Gain_k;
  } else {
    boost_and_two_level__1_sm_ehs_B.Switch_e =
      boost_and_two_level__1_sm_ehs_B.Memory_e;
  }

  /* End of Switch: '<S16>/Switch' */

  /* Product: '<S5>/Divide' */
  boost_and_two_level__1_sm_ehs_B.P_PV =
    boost_and_two_level__1_sm_ehs_B.Switch_e *
    boost_and_two_level__1_sm_ehs_B.Switch;
  if (rtmIsMajorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
    /* Level2 S-Function Block: '<S2>/Outputs_eHS1_Recv' (sfun_fct_op7160ex1_recv) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[9];
      sfcnOutputs(rts, 1);
    }

    /* Level2 S-Function Block: '<S2>/sfp2dbl' (sfun_SFP2DBL) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[10];
      sfcnOutputs(rts, 1);
    }

    /* Level2 S-Function Block: '<S2>/OpMonitor' (opmonitor) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[11];
      sfcnOutputs(rts, 1);
    }

    /* Level2 S-Function Block: '<S30>/S-Function' (OP_SEND) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[12];
      sfcnOutputs(rts, 1);
    }

    /* Level2 S-Function Block: '<S17>/Convert double to  Single floating-point (FPGA)' (sfun_DBL2SFP) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[13];
      sfcnOutputs(rts, 1);
    }

    /* Level2 S-Function Block: '<S17>/Convert double to  Single floating-point (FPGA)1' (sfun_DBL2SFP) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[14];
      sfcnOutputs(rts, 1);
    }

    /* Memory: '<S19>/Memory1' */
    boost_and_two_level__1_sm_ehs_B.Memory1_o =
      boost_and_two_level__1_sm_ehs_DW.Memory1_PreviousInput_h;
    for (i = 0; i < 16; i++) {
      /* Switch: '<S22>/Switch' incorporates:
       *  Constant: '<S22>/Constant'
       */
      if (boost_and_two_level__1_sm_ehs_P.Constant_Value_gg) {
        /* Saturate: '<S22>/Saturation' */
        B = boost_and_two_level__1_sm_ehs_B.SFunction_l[i + 4];
        dP = boost_and_two_level__1_sm_ehs_P.Saturation_LowerSat_g;
        u2 = boost_and_two_level__1_sm_ehs_P.Saturation_UpperSat_d;
        if (B > u2) {
          B = u2;
        } else {
          if (B < dP) {
            B = dP;
          }
        }

        boost_and_two_level__1_sm_ehs_B.Saturation_f[i] = B;

        /* End of Saturate: '<S22>/Saturation' */

        /* Gain: '<S22>/Gain1' */
        boost_and_two_level__1_sm_ehs_B.Gain1_l[i] =
          boost_and_two_level__1_sm_ehs_P.Gain1_Gain_p *
          boost_and_two_level__1_sm_ehs_B.Saturation_f[i];

        /* DataTypeConversion: '<S22>/Data Type Conversion2' */
        B = boost_and_two_level__1_sm_ehs_B.Gain1_l[i];
        if (B < 0.0) {
          B = ceil(B);
        } else {
          B = floor(B);
        }

        if (rtIsNaN(B) || rtIsInf(B)) {
          B = 0.0;
        } else {
          B = fmod(B, 4.294967296E+9);
        }

        boost_and_two_level__1_sm_ehs_B.DataTypeConversion2[i] = B < 0.0 ?
          -(int32_T)(uint32_T)-B : (int32_T)(uint32_T)B;

        /* End of DataTypeConversion: '<S22>/Data Type Conversion2' */
        boost_and_two_level__1_sm_ehs_B.Switch_n[i] = (uint32_T)
          boost_and_two_level__1_sm_ehs_B.DataTypeConversion2[i];
      } else {
        /* Saturate: '<S22>/Saturation1' */
        B = boost_and_two_level__1_sm_ehs_B.SFunction_l[i + 4];
        dP = boost_and_two_level__1_sm_ehs_P.Saturation1_LowerSat_a;
        u2 = boost_and_two_level__1_sm_ehs_P.Saturation1_UpperSat_f;
        if (B > u2) {
          B = u2;
        } else {
          if (B < dP) {
            B = dP;
          }
        }

        boost_and_two_level__1_sm_ehs_B.Saturation1_p[i] = B;

        /* End of Saturate: '<S22>/Saturation1' */

        /* Gain: '<S22>/Gain2' */
        boost_and_two_level__1_sm_ehs_B.Gain2[i] =
          boost_and_two_level__1_sm_ehs_P.Gain2_Gain_p *
          boost_and_two_level__1_sm_ehs_B.Saturation1_p[i];

        /* DataTypeConversion: '<S22>/Data Type Conversion1' */
        B = boost_and_two_level__1_sm_ehs_B.Gain2[i];
        if (B < 0.0) {
          B = ceil(B);
        } else {
          B = floor(B);
        }

        if (rtIsNaN(B) || rtIsInf(B)) {
          B = 0.0;
        } else {
          B = fmod(B, 4.294967296E+9);
        }

        boost_and_two_level__1_sm_ehs_B.DataTypeConversion1_l[i] = B < 0.0 ?
          (uint32_T)-(int32_T)(uint32_T)-B : (uint32_T)B;

        /* End of DataTypeConversion: '<S22>/Data Type Conversion1' */
        boost_and_two_level__1_sm_ehs_B.Switch_n[i] =
          boost_and_two_level__1_sm_ehs_B.DataTypeConversion1_l[i];
      }

      /* End of Switch: '<S22>/Switch' */

      /* S-Function (sfix_bitop): '<S22>/Bitwise Operator' */
      boost_and_two_level__1_sm_ehs_B.BitwiseOperator[i] =
        boost_and_two_level__1_sm_ehs_B.Switch_n[i] &
        boost_and_two_level__1_sm_ehs_P.BitwiseOperator_BitMask;

      /* DataTypeConversion: '<S22>/Data Type Conversion5' */
      boost_and_two_level__1_sm_ehs_B.DataTypeConversion5[i] =
        boost_and_two_level__1_sm_ehs_B.BitwiseOperator[i];

      /* Switch: '<S20>/Switch' incorporates:
       *  Constant: '<S20>/Constant'
       */
      if (boost_and_two_level__1_sm_ehs_P.Constant_Value_n) {
        /* Saturate: '<S20>/Saturation' */
        B = boost_and_two_level__1_sm_ehs_B.SFunction_l[i + 100];
        dP = boost_and_two_level__1_sm_ehs_P.Saturation_LowerSat;
        u2 = boost_and_two_level__1_sm_ehs_P.Saturation_UpperSat;
        if (B > u2) {
          B = u2;
        } else {
          if (B < dP) {
            B = dP;
          }
        }

        boost_and_two_level__1_sm_ehs_B.Saturation_h[i] = B;

        /* End of Saturate: '<S20>/Saturation' */

        /* Gain: '<S20>/Gain1' */
        boost_and_two_level__1_sm_ehs_B.Gain1_g[i] =
          boost_and_two_level__1_sm_ehs_P.Gain1_Gain *
          boost_and_two_level__1_sm_ehs_B.Saturation_h[i];

        /* DataTypeConversion: '<S20>/Data Type Conversion2' */
        B = boost_and_two_level__1_sm_ehs_B.Gain1_g[i];
        if (B < 0.0) {
          B = ceil(B);
        } else {
          B = floor(B);
        }

        if (rtIsNaN(B) || rtIsInf(B)) {
          B = 0.0;
        } else {
          B = fmod(B, 4.294967296E+9);
        }

        boost_and_two_level__1_sm_ehs_B.DataTypeConversion2_n[i] = B < 0.0 ?
          -(int32_T)(uint32_T)-B : (int32_T)(uint32_T)B;

        /* End of DataTypeConversion: '<S20>/Data Type Conversion2' */
        boost_and_two_level__1_sm_ehs_B.Switch_g[i] = (uint32_T)
          boost_and_two_level__1_sm_ehs_B.DataTypeConversion2_n[i];
      } else {
        /* Saturate: '<S20>/Saturation1' */
        B = boost_and_two_level__1_sm_ehs_B.SFunction_l[i + 100];
        dP = boost_and_two_level__1_sm_ehs_P.Saturation1_LowerSat;
        u2 = boost_and_two_level__1_sm_ehs_P.Saturation1_UpperSat;
        if (B > u2) {
          B = u2;
        } else {
          if (B < dP) {
            B = dP;
          }
        }

        boost_and_two_level__1_sm_ehs_B.Saturation1_cu[i] = B;

        /* End of Saturate: '<S20>/Saturation1' */

        /* Gain: '<S20>/Gain2' */
        boost_and_two_level__1_sm_ehs_B.Gain2_a[i] =
          boost_and_two_level__1_sm_ehs_P.Gain2_Gain *
          boost_and_two_level__1_sm_ehs_B.Saturation1_cu[i];

        /* DataTypeConversion: '<S20>/Data Type Conversion1' */
        B = boost_and_two_level__1_sm_ehs_B.Gain2_a[i];
        if (B < 0.0) {
          B = ceil(B);
        } else {
          B = floor(B);
        }

        if (rtIsNaN(B) || rtIsInf(B)) {
          B = 0.0;
        } else {
          B = fmod(B, 4.294967296E+9);
        }

        boost_and_two_level__1_sm_ehs_B.DataTypeConversion1_hm[i] = B < 0.0 ?
          (uint32_T)-(int32_T)(uint32_T)-B : (uint32_T)B;

        /* End of DataTypeConversion: '<S20>/Data Type Conversion1' */
        boost_and_two_level__1_sm_ehs_B.Switch_g[i] =
          boost_and_two_level__1_sm_ehs_B.DataTypeConversion1_hm[i];
      }

      /* End of Switch: '<S20>/Switch' */

      /* S-Function (sfix_bitop): '<S20>/Bitwise Operator' */
      boost_and_two_level__1_sm_ehs_B.BitwiseOperator_p[i] =
        boost_and_two_level__1_sm_ehs_B.Switch_g[i] &
        boost_and_two_level__1_sm_ehs_P.BitwiseOperator_BitMask_j;

      /* DataTypeConversion: '<S20>/Data Type Conversion5' */
      boost_and_two_level__1_sm_ehs_B.DataTypeConversion5_m[i] =
        boost_and_two_level__1_sm_ehs_B.BitwiseOperator_p[i];

      /* Switch: '<S21>/Switch' incorporates:
       *  Constant: '<S21>/Constant'
       */
      if (boost_and_two_level__1_sm_ehs_P.Constant_Value_kj) {
        /* Saturate: '<S21>/Saturation' */
        B = boost_and_two_level__1_sm_ehs_B.SFunction_l[i + 132];
        dP = boost_and_two_level__1_sm_ehs_P.Saturation_LowerSat_h;
        u2 = boost_and_two_level__1_sm_ehs_P.Saturation_UpperSat_i;
        if (B > u2) {
          B = u2;
        } else {
          if (B < dP) {
            B = dP;
          }
        }

        boost_and_two_level__1_sm_ehs_B.Saturation_g[i] = B;

        /* End of Saturate: '<S21>/Saturation' */

        /* Gain: '<S21>/Gain1' */
        boost_and_two_level__1_sm_ehs_B.Gain1_j[i] =
          boost_and_two_level__1_sm_ehs_P.Gain1_Gain_j *
          boost_and_two_level__1_sm_ehs_B.Saturation_g[i];

        /* DataTypeConversion: '<S21>/Data Type Conversion2' */
        B = boost_and_two_level__1_sm_ehs_B.Gain1_j[i];
        if (B < 0.0) {
          B = ceil(B);
        } else {
          B = floor(B);
        }

        if (rtIsNaN(B) || rtIsInf(B)) {
          B = 0.0;
        } else {
          B = fmod(B, 4.294967296E+9);
        }

        boost_and_two_level__1_sm_ehs_B.DataTypeConversion2_c[i] = B < 0.0 ?
          -(int32_T)(uint32_T)-B : (int32_T)(uint32_T)B;

        /* End of DataTypeConversion: '<S21>/Data Type Conversion2' */
        boost_and_two_level__1_sm_ehs_B.Switch_d[i] = (uint32_T)
          boost_and_two_level__1_sm_ehs_B.DataTypeConversion2_c[i];
      } else {
        /* Saturate: '<S21>/Saturation1' */
        B = boost_and_two_level__1_sm_ehs_B.SFunction_l[i + 132];
        dP = boost_and_two_level__1_sm_ehs_P.Saturation1_LowerSat_n;
        u2 = boost_and_two_level__1_sm_ehs_P.Saturation1_UpperSat_k;
        if (B > u2) {
          B = u2;
        } else {
          if (B < dP) {
            B = dP;
          }
        }

        boost_and_two_level__1_sm_ehs_B.Saturation1_c[i] = B;

        /* End of Saturate: '<S21>/Saturation1' */

        /* Gain: '<S21>/Gain2' */
        boost_and_two_level__1_sm_ehs_B.Gain2_d[i] =
          boost_and_two_level__1_sm_ehs_P.Gain2_Gain_d *
          boost_and_two_level__1_sm_ehs_B.Saturation1_c[i];

        /* DataTypeConversion: '<S21>/Data Type Conversion1' */
        B = boost_and_two_level__1_sm_ehs_B.Gain2_d[i];
        if (B < 0.0) {
          B = ceil(B);
        } else {
          B = floor(B);
        }

        if (rtIsNaN(B) || rtIsInf(B)) {
          B = 0.0;
        } else {
          B = fmod(B, 4.294967296E+9);
        }

        boost_and_two_level__1_sm_ehs_B.DataTypeConversion1_d[i] = B < 0.0 ?
          (uint32_T)-(int32_T)(uint32_T)-B : (uint32_T)B;

        /* End of DataTypeConversion: '<S21>/Data Type Conversion1' */
        boost_and_two_level__1_sm_ehs_B.Switch_d[i] =
          boost_and_two_level__1_sm_ehs_B.DataTypeConversion1_d[i];
      }

      /* End of Switch: '<S21>/Switch' */

      /* S-Function (sfix_bitop): '<S21>/Bitwise Operator' */
      boost_and_two_level__1_sm_ehs_B.BitwiseOperator_pi[i] =
        boost_and_two_level__1_sm_ehs_B.Switch_d[i] &
        boost_and_two_level__1_sm_ehs_P.BitwiseOperator_BitMask_h;

      /* DataTypeConversion: '<S21>/Data Type Conversion5' */
      boost_and_two_level__1_sm_ehs_B.DataTypeConversion5_n[i] =
        boost_and_two_level__1_sm_ehs_B.BitwiseOperator_pi[i];
    }

    /* UnitDelay: '<S24>/Delay Input1' */
    memcpy(&boost_and_two_level__1_sm_ehs_B.Uk1_g[0],
           &boost_and_two_level__1_sm_ehs_DW.DelayInput1_DSTATE_m[0], 80U *
           sizeof(uint32_T));

    /* RelationalOperator: '<S24>/FixPt Relational Operator' */
    for (i = 0; i < 16; i++) {
      boost_and_two_level__1_sm_ehs_B.FixPtRelationalOperator_i[i] = (uint8_T)
        (boost_and_two_level__1_sm_ehs_B.DataTypeConversion5[i] !=
         boost_and_two_level__1_sm_ehs_B.Uk1_g[i]);
    }

    for (i = 0; i < 16; i++) {
      boost_and_two_level__1_sm_ehs_B.FixPtRelationalOperator_i[i + 16] =
        (uint8_T)(boost_and_two_level__1_sm_ehs_B.Uk1_g[i + 16] !=
                  boost_and_two_level__1_sm_ehs_B.ConvertdoubletoSinglefloating_o
                  [i]);
    }

    for (i = 0; i < 16; i++) {
      boost_and_two_level__1_sm_ehs_B.FixPtRelationalOperator_i[i + 32] =
        (uint8_T)(boost_and_two_level__1_sm_ehs_B.Uk1_g[i + 32] !=
                  boost_and_two_level__1_sm_ehs_B.ConvertdoubletoSinglefloatin_ou
                  [i]);
    }

    for (i = 0; i < 16; i++) {
      boost_and_two_level__1_sm_ehs_B.FixPtRelationalOperator_i[i + 48] =
        (uint8_T)(boost_and_two_level__1_sm_ehs_B.Uk1_g[i + 48] !=
                  boost_and_two_level__1_sm_ehs_B.DataTypeConversion5_m[i]);
    }

    for (i = 0; i < 16; i++) {
      boost_and_two_level__1_sm_ehs_B.FixPtRelationalOperator_i[i + 64] =
        (uint8_T)(boost_and_two_level__1_sm_ehs_B.Uk1_g[i + 64] !=
                  boost_and_two_level__1_sm_ehs_B.DataTypeConversion5_n[i]);
    }

    /* End of RelationalOperator: '<S24>/FixPt Relational Operator' */

    /* Logic: '<S23>/Logical Operator' */
    tmp = (boost_and_two_level__1_sm_ehs_B.FixPtRelationalOperator_i[0] != 0);
    for (i = 0; i < 79; i++) {
      tmp = (tmp || (boost_and_two_level__1_sm_ehs_B.FixPtRelationalOperator_i[i
                     + 1] != 0));
    }

    boost_and_two_level__1_sm_ehs_B.LogicalOperator_i = tmp;

    /* End of Logic: '<S23>/Logical Operator' */

    /* DataTypeConversion: '<S23>/Data Type Conversion1' */
    boost_and_two_level__1_sm_ehs_B.DataTypeConversion1 =
      boost_and_two_level__1_sm_ehs_B.LogicalOperator_i;

    /* Logic: '<S19>/Logical Operator1' */
    boost_and_two_level__1_sm_ehs_B.LogicalOperator1 =
      ((boost_and_two_level__1_sm_ehs_B.Memory1_o != 0.0) ||
       (boost_and_two_level__1_sm_ehs_B.DataTypeConversion1 != 0.0));

    /* DataTypeConversion: '<S19>/Data Type Conversion' */
    boost_and_two_level__1_sm_ehs_B.DataTypeConversion_d =
      boost_and_two_level__1_sm_ehs_B.LogicalOperator1;

    /* Memory: '<S19>/Memory2' */
    boost_and_two_level__1_sm_ehs_B.Memory2 =
      boost_and_two_level__1_sm_ehs_DW.Memory2_PreviousInput;

    /* Sum: '<S19>/Add' */
    boost_and_two_level__1_sm_ehs_B.Add =
      boost_and_two_level__1_sm_ehs_B.DataTypeConversion_d +
      boost_and_two_level__1_sm_ehs_B.Memory2;

    /* Switch: '<S19>/Switch' incorporates:
     *  Constant: '<S19>/Constant'
     */
    if (boost_and_two_level__1_sm_ehs_B.Add >=
        boost_and_two_level__1_sm_ehs_P.Counter_max_count) {
      boost_and_two_level__1_sm_ehs_B.Switch_i =
        boost_and_two_level__1_sm_ehs_P.Constant_Value;
    } else {
      boost_and_two_level__1_sm_ehs_B.Switch_i =
        boost_and_two_level__1_sm_ehs_B.Add;
    }

    /* End of Switch: '<S19>/Switch' */

    /* Switch: '<S19>/Switch1' incorporates:
     *  Constant: '<S19>/Constant1'
     *  Constant: '<S19>/Constant2'
     */
    if (boost_and_two_level__1_sm_ehs_B.Switch_i != 0.0) {
      boost_and_two_level__1_sm_ehs_B.Switch1 =
        boost_and_two_level__1_sm_ehs_P.Constant2_Value;
    } else {
      boost_and_two_level__1_sm_ehs_B.Switch1 =
        boost_and_two_level__1_sm_ehs_P.Constant1_Value;
    }

    /* End of Switch: '<S19>/Switch1' */

    /* MultiPortSwitch: '<S17>/Multiport Switch' incorporates:
     *  Constant: '<S17>/addr'
     *  Constant: '<S17>/addr1'
     *  Constant: '<S17>/addr2'
     *  Constant: '<S17>/addr3'
     *  Constant: '<S17>/addr4'
     *  Constant: '<S17>/default'
     */
    switch ((int32_T)boost_and_two_level__1_sm_ehs_B.Add) {
     case 1:
      boost_and_two_level__1_sm_ehs_B.MultiportSwitch[0] =
        boost_and_two_level__1_sm_ehs_P.addr_Value;
      memcpy(&boost_and_two_level__1_sm_ehs_B.MultiportSwitch[1],
             &boost_and_two_level__1_sm_ehs_B.DataTypeConversion5[0], sizeof
             (uint32_T) << 4U);
      break;

     case 2:
      boost_and_two_level__1_sm_ehs_B.MultiportSwitch[0] =
        boost_and_two_level__1_sm_ehs_P.addr1_Value;
      memcpy(&boost_and_two_level__1_sm_ehs_B.MultiportSwitch[1],
             &boost_and_two_level__1_sm_ehs_B.ConvertdoubletoSinglefloating_o[0],
             sizeof(uint32_T) << 4U);
      break;

     case 3:
      boost_and_two_level__1_sm_ehs_B.MultiportSwitch[0] =
        boost_and_two_level__1_sm_ehs_P.addr2_Value;
      memcpy(&boost_and_two_level__1_sm_ehs_B.MultiportSwitch[1],
             &boost_and_two_level__1_sm_ehs_B.ConvertdoubletoSinglefloatin_ou[0],
             sizeof(uint32_T) << 4U);
      break;

     case 4:
      boost_and_two_level__1_sm_ehs_B.MultiportSwitch[0] =
        boost_and_two_level__1_sm_ehs_P.addr3_Value;
      memcpy(&boost_and_two_level__1_sm_ehs_B.MultiportSwitch[1],
             &boost_and_two_level__1_sm_ehs_B.DataTypeConversion5_m[0], sizeof
             (uint32_T) << 4U);
      break;

     case 5:
      boost_and_two_level__1_sm_ehs_B.MultiportSwitch[0] =
        boost_and_two_level__1_sm_ehs_P.addr4_Value;
      memcpy(&boost_and_two_level__1_sm_ehs_B.MultiportSwitch[1],
             &boost_and_two_level__1_sm_ehs_B.DataTypeConversion5_n[0], sizeof
             (uint32_T) << 4U);
      break;

     default:
      memcpy(&boost_and_two_level__1_sm_ehs_B.MultiportSwitch[0],
             &boost_and_two_level__1_sm_ehs_P.default_Value[0], 17U * sizeof
             (uint32_T));
      break;
    }

    /* End of MultiPortSwitch: '<S17>/Multiport Switch' */

    /* SignalConversion: '<S17>/TmpSignal ConversionAtLoadInInport2' incorporates:
     *  Constant: '<S17>/blockID'
     */
    boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtLoadInInpo[0] =
      boost_and_two_level__1_sm_ehs_P.blockID_Value;
    memcpy(&boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtLoadInInpo[1],
           &boost_and_two_level__1_sm_ehs_B.MultiportSwitch[0], 17U * sizeof
           (uint32_T));

    /* Level2 S-Function Block: '<S17>/LoadIn' (sfun_fct_op7160ex1_load_in) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[15];
      sfcnOutputs(rts, 1);
    }

    /* Constant: '<S18>/Constant' */
    memcpy(&boost_and_two_level__1_sm_ehs_B.Constant[0],
           &boost_and_two_level__1_sm_ehs_P.Constant_Value_g[0], 9U * sizeof
           (real_T));

    /* Level2 S-Function Block: '<S17>/Convert double to  Single floating-point (FPGA)2' (sfun_DBL2SFP) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[16];
      sfcnOutputs(rts, 1);
    }

    /* Level2 S-Function Block: '<S17>/DataIn Send' (sfun_fct_op7160ex1_send) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[17];
      sfcnOutputs(rts, 1);
    }

    /* Level2 S-Function Block: '<S2>/OpWriteFile' (opwritefile) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[18];
      sfcnOutputs(rts, 1);
    }

    for (i = 0; i < 10; i++) {
      /* Constant: '<S7>/load_config1' */
      boost_and_two_level__1_sm_ehs_B.load_config1[i] =
        boost_and_two_level__1_sm_ehs_P.load_config1_Value[i];

      /* UnitDelay: '<S28>/Delay Input1' */
      boost_and_two_level__1_sm_ehs_B.Uk1_gd[i] =
        boost_and_two_level__1_sm_ehs_DW.DelayInput1_DSTATE_k[i];

      /* RelationalOperator: '<S28>/FixPt Relational Operator' */
      boost_and_two_level__1_sm_ehs_B.FixPtRelationalOperator_h[i] = (uint8_T)
        (boost_and_two_level__1_sm_ehs_B.load_config1[i] !=
         boost_and_two_level__1_sm_ehs_B.Uk1_gd[i]);
    }

    tmp = (boost_and_two_level__1_sm_ehs_B.FixPtRelationalOperator_h[0] != 0);
    for (i = 0; i < 9; i++) {
      tmp = (tmp || (boost_and_two_level__1_sm_ehs_B.FixPtRelationalOperator_h[i
                     + 1] != 0));
    }

    boost_and_two_level__1_sm_ehs_B.LogicalOperator_g = tmp;

    /* DataTypeConversion: '<S27>/Data Type Conversion1' */
    boost_and_two_level__1_sm_ehs_B.DataTypeConversion1_o =
      boost_and_two_level__1_sm_ehs_B.LogicalOperator_g;

    /* Level2 S-Function Block: '<S9>/RTE Conversion1' (rte_conversion) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[19];
      sfcnOutputs(rts, 1);
    }

    /* Step: '<S9>/Step' */
    dP = boost_and_two_level__1_sm_ehs_M->Timing.t[1];
    if (dP < boost_and_two_level__1_sm_ehs_P.Step_Time) {
      boost_and_two_level__1_sm_ehs_B.Step =
        boost_and_two_level__1_sm_ehs_P.Step_Y0;
    } else {
      boost_and_two_level__1_sm_ehs_B.Step =
        boost_and_two_level__1_sm_ehs_P.Step_YFinal;
    }

    /* End of Step: '<S9>/Step' */

    /* Product: '<S9>/Product1' */
    boost_and_two_level__1_sm_ehs_B.Product1 =
      boost_and_two_level__1_sm_ehs_B.SFunction_l[2] *
      boost_and_two_level__1_sm_ehs_B.Step;

    /* SignalConversion: '<S9>/TmpSignal ConversionAtRTE ConversionInport1' */
    boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtRTEConvers[0] =
      boost_and_two_level__1_sm_ehs_B.Product1;
    boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtRTEConvers[1] =
      boost_and_two_level__1_sm_ehs_B.Product1;
    boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtRTEConvers[2] =
      boost_and_two_level__1_sm_ehs_B.Product1;
    boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtRTEConvers[3] =
      boost_and_two_level__1_sm_ehs_B.Product1;

    /* Level2 S-Function Block: '<S9>/RTE Conversion' (rte_conversion) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[20];
      sfcnOutputs(rts, 1);
    }

    /* Level2 S-Function Block: '<S9>/RTE Logical Operator1' (rte_logical_operator) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[21];
      sfcnOutputs(rts, 1);
    }

    /* Level2 S-Function Block: '<S9>/RTE Ground' (rte_ground) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[22];
      sfcnOutputs(rts, 1);
    }

    /* Level2 S-Function Block: '<S26>/RTE_Conversion_1' (rte_conversion_ophsdio) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[23];
      sfcnOutputs(rts, 1);
    }

    /* Level2 S-Function Block: '<S26>/EventGen_eHS_1' (sfun_op7160ex1_event_generator) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[24];
      sfcnOutputs(rts, 1);
    }

    /* RateTransition: '<S9>/Rate Transition5' */
    if (rtmIsMajorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
      boost_and_two_level__1_sm_ehs_B.RateTransition5 =
        boost_and_two_level__1_sm_ehs_B.Divide[4];
    }

    /* End of RateTransition: '<S9>/Rate Transition5' */

    /* Gain: '<S36>/A->pu' */
    P = boost_and_two_level__1_sm_ehs_P.InverterControl_Vnom_prim;
    B = boost_and_two_level__1_sm_ehs_P.InverterControl_Pnom;
    dV = P / B;
    dV /= 1.4142135623730951;
    boost_and_two_level__1_sm_ehs_B.Apu = dV *
      boost_and_two_level__1_sm_ehs_B.RateTransition5;

    /* UnitDelay: '<S56>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_B.UnitDelay_l =
      boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_b;

    /* Saturate: '<S45>/avoid division by zero' */
    B = boost_and_two_level__1_sm_ehs_B.UnitDelay_l;
    dP = boost_and_two_level__1_sm_ehs_P.avoiddivisionbyzero_LowerSat;
    u2 = boost_and_two_level__1_sm_ehs_P.avoiddivisionbyzero_UpperSat;
    if (B > u2) {
      boost_and_two_level__1_sm_ehs_B.avoiddivisionbyzero = u2;
    } else if (B < dP) {
      boost_and_two_level__1_sm_ehs_B.avoiddivisionbyzero = dP;
    } else {
      boost_and_two_level__1_sm_ehs_B.avoiddivisionbyzero = B;
    }

    /* End of Saturate: '<S45>/avoid division by zero' */

    /* Math: '<S45>/Math Function'
     *
     * About '<S45>/Math Function':
     *  Operator: reciprocal
     */
    P = boost_and_two_level__1_sm_ehs_B.avoiddivisionbyzero;
    boost_and_two_level__1_sm_ehs_B.MathFunction_i = 1.0 / P;

    /* Gain: '<S45>/Gain' */
    boost_and_two_level__1_sm_ehs_B.Gain =
      boost_and_two_level__1_sm_ehs_P.Gain_Gain_e *
      boost_and_two_level__1_sm_ehs_B.MathFunction_i;

    /* Level2 S-Function Block: '<S82>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[25];
      sfcnOutputs(rts, 1);
    }

    /* DiscreteIntegrator: '<S56>/Discrete-Time Integrator' */
    boost_and_two_level__1_sm_ehs_B.DiscreteTimeIntegrator =
      boost_and_two_level__1_sm_ehs_DW.DiscreteTimeIntegrator_DSTATE;

    /* Math: '<S56>/Math Function' incorporates:
     *  Constant: '<S56>/Constant4'
     */
    boost_and_two_level__1_sm_ehs_B.MathFunction_m = rt_modd_snf
      (boost_and_two_level__1_sm_ehs_B.DiscreteTimeIntegrator,
       boost_and_two_level__1_sm_ehs_P.Constant4_Value_o);

    /* RelationalOperator: '<S83>/Compare' incorporates:
     *  Constant: '<S81>/Constant'
     *  Constant: '<S83>/Constant'
     */
    boost_and_two_level__1_sm_ehs_B.Compare = (uint8_T)
      (boost_and_two_level__1_sm_ehs_P.AlphaBetaZerotodq0_Alignment ==
       boost_and_two_level__1_sm_ehs_P.CompareToConstant_const);

    /* Outputs for Enabled SubSystem: '<S81>/Subsystem1' incorporates:
     *  EnablePort: '<S86>/Enable'
     */
    if (boost_and_two_level__1_sm_ehs_B.Compare > 0) {
      /* Fcn: '<S86>/Fcn' */
      boost_and_two_level__1_sm_ehs_B.Fcn = boost_and_two_level__1_sm_ehs_B.Apu *
        cos(boost_and_two_level__1_sm_ehs_B.MathFunction_m) +
        boost_and_two_level__1_sm_ehs_B.SFunction_h * sin
        (boost_and_two_level__1_sm_ehs_B.MathFunction_m);

      /* Fcn: '<S86>/Fcn1' */
      boost_and_two_level__1_sm_ehs_B.Fcn1 =
        -boost_and_two_level__1_sm_ehs_B.Apu * sin
        (boost_and_two_level__1_sm_ehs_B.MathFunction_m) +
        boost_and_two_level__1_sm_ehs_B.SFunction_h * cos
        (boost_and_two_level__1_sm_ehs_B.MathFunction_m);
    }

    /* End of Outputs for SubSystem: '<S81>/Subsystem1' */

    /* RelationalOperator: '<S84>/Compare' incorporates:
     *  Constant: '<S81>/Constant'
     *  Constant: '<S84>/Constant'
     */
    boost_and_two_level__1_sm_ehs_B.Compare_m = (uint8_T)
      (boost_and_two_level__1_sm_ehs_P.AlphaBetaZerotodq0_Alignment ==
       boost_and_two_level__1_sm_ehs_P.CompareToConstant1_const);

    /* Outputs for Enabled SubSystem: '<S81>/Subsystem - pi//2 delay' incorporates:
     *  EnablePort: '<S85>/Enable'
     */
    if (boost_and_two_level__1_sm_ehs_B.Compare_m > 0) {
      /* Fcn: '<S85>/Fcn' */
      boost_and_two_level__1_sm_ehs_B.Fcn_i =
        boost_and_two_level__1_sm_ehs_B.Apu * sin
        (boost_and_two_level__1_sm_ehs_B.MathFunction_m) -
        boost_and_two_level__1_sm_ehs_B.SFunction_h * cos
        (boost_and_two_level__1_sm_ehs_B.MathFunction_m);

      /* Fcn: '<S85>/Fcn1' */
      boost_and_two_level__1_sm_ehs_B.Fcn1_k =
        boost_and_two_level__1_sm_ehs_B.Apu * cos
        (boost_and_two_level__1_sm_ehs_B.MathFunction_m) +
        boost_and_two_level__1_sm_ehs_B.SFunction_h * sin
        (boost_and_two_level__1_sm_ehs_B.MathFunction_m);
    }

    /* End of Outputs for SubSystem: '<S81>/Subsystem - pi//2 delay' */

    /* Step: '<S45>/First cycle of simulation Id=0.92, Iq=0' */
    dP = boost_and_two_level__1_sm_ehs_M->Timing.t[1];
    B = boost_and_two_level__1_sm_ehs_P.InverterControl_Fnom;
    dV = 1.0 / B;
    if (dP < dV) {
      boost_and_two_level__1_sm_ehs_B.FirstcycleofsimulationId092Iq0 =
        boost_and_two_level__1_sm_ehs_P.FirstcycleofsimulationId092Iq0_;
    } else {
      boost_and_two_level__1_sm_ehs_B.FirstcycleofsimulationId092Iq0 =
        boost_and_two_level__1_sm_ehs_P.FirstcycleofsimulationId092Iq_k;
    }

    /* End of Step: '<S45>/First cycle of simulation Id=0.92, Iq=0' */

    /* Switch: '<S45>/Switch' incorporates:
     *  Constant: '<S45>/Constant'
     */
    if (boost_and_two_level__1_sm_ehs_B.FirstcycleofsimulationId092Iq0 != 0.0) {
      /* Switch: '<S81>/Switch' */
      if (boost_and_two_level__1_sm_ehs_B.Compare != 0) {
        boost_and_two_level__1_sm_ehs_B.Switch_cq[0] =
          boost_and_two_level__1_sm_ehs_B.Fcn;
      } else {
        boost_and_two_level__1_sm_ehs_B.Switch_cq[0] =
          boost_and_two_level__1_sm_ehs_B.Fcn_i;
      }

      if (boost_and_two_level__1_sm_ehs_B.Compare != 0) {
        boost_and_two_level__1_sm_ehs_B.Switch_cq[1] =
          boost_and_two_level__1_sm_ehs_B.Fcn1;
      } else {
        boost_and_two_level__1_sm_ehs_B.Switch_cq[1] =
          boost_and_two_level__1_sm_ehs_B.Fcn1_k;
      }

      /* End of Switch: '<S81>/Switch' */
      boost_and_two_level__1_sm_ehs_B.Switch_b[0] =
        boost_and_two_level__1_sm_ehs_B.Switch_cq[0];
      boost_and_two_level__1_sm_ehs_B.Switch_b[1] =
        boost_and_two_level__1_sm_ehs_B.Switch_cq[1];
    } else {
      boost_and_two_level__1_sm_ehs_B.Switch_b[0] =
        boost_and_two_level__1_sm_ehs_P.Constant_Value_k[0];
      boost_and_two_level__1_sm_ehs_B.Switch_b[1] =
        boost_and_two_level__1_sm_ehs_P.Constant_Value_k[1];
    }

    /* End of Switch: '<S45>/Switch' */

    /* Gain: '<S77>/D*u(k)' */
    boost_and_two_level__1_sm_ehs_B.Duk[0] =
      boost_and_two_level__1_sm_ehs_P.Duk_Gain *
      boost_and_two_level__1_sm_ehs_B.Switch_b[0];
    boost_and_two_level__1_sm_ehs_B.Duk[1] =
      boost_and_two_level__1_sm_ehs_P.Duk_Gain *
      boost_and_two_level__1_sm_ehs_B.Switch_b[1];

    /* UnitDelay: '<S77>/Delay_x1' */
    boost_and_two_level__1_sm_ehs_B.x1k[0] =
      boost_and_two_level__1_sm_ehs_DW.Delay_x1_DSTATE[0];
    boost_and_two_level__1_sm_ehs_B.x1k[1] =
      boost_and_two_level__1_sm_ehs_DW.Delay_x1_DSTATE[1];

    /* Gain: '<S80>/C11' */
    boost_and_two_level__1_sm_ehs_B.C11[0] =
      boost_and_two_level__1_sm_ehs_P.C11_Gain *
      boost_and_two_level__1_sm_ehs_B.x1k[0];
    boost_and_two_level__1_sm_ehs_B.C11[1] =
      boost_and_two_level__1_sm_ehs_P.C11_Gain *
      boost_and_two_level__1_sm_ehs_B.x1k[1];

    /* UnitDelay: '<S77>/Delay_x2' */
    boost_and_two_level__1_sm_ehs_B.x2k[0] =
      boost_and_two_level__1_sm_ehs_DW.Delay_x2_DSTATE[0];
    boost_and_two_level__1_sm_ehs_B.x2k[1] =
      boost_and_two_level__1_sm_ehs_DW.Delay_x2_DSTATE[1];

    /* Gain: '<S80>/C12' */
    boost_and_two_level__1_sm_ehs_B.C12[0] =
      boost_and_two_level__1_sm_ehs_P.C12_Gain *
      boost_and_two_level__1_sm_ehs_B.x2k[0];
    boost_and_two_level__1_sm_ehs_B.C12[1] =
      boost_and_two_level__1_sm_ehs_P.C12_Gain *
      boost_and_two_level__1_sm_ehs_B.x2k[1];

    /* Sum: '<S80>/sum2' */
    boost_and_two_level__1_sm_ehs_B.sum2[0] =
      boost_and_two_level__1_sm_ehs_B.C11[0] +
      boost_and_two_level__1_sm_ehs_B.C12[0];
    boost_and_two_level__1_sm_ehs_B.sum2[1] =
      boost_and_two_level__1_sm_ehs_B.C11[1] +
      boost_and_two_level__1_sm_ehs_B.C12[1];

    /* Sum: '<S77>/C*X(k)+D*u(k)' */
    boost_and_two_level__1_sm_ehs_B.yk[0] = boost_and_two_level__1_sm_ehs_B.Duk
      [0] + boost_and_two_level__1_sm_ehs_B.sum2[0];
    boost_and_two_level__1_sm_ehs_B.yk[1] = boost_and_two_level__1_sm_ehs_B.Duk
      [1] + boost_and_two_level__1_sm_ehs_B.sum2[1];

    /* UnitDelay: '<S31>/Unit Delay2' */
    boost_and_two_level__1_sm_ehs_B.UnitDelay2 =
      boost_and_two_level__1_sm_ehs_DW.UnitDelay2_DSTATE;

    /* Sum: '<S33>/Sum' incorporates:
     *  Constant: '<S31>/Iq_ref'
     */
    boost_and_two_level__1_sm_ehs_B.Sum_n[0] =
      boost_and_two_level__1_sm_ehs_B.UnitDelay2 -
      boost_and_two_level__1_sm_ehs_B.yk[0];
    boost_and_two_level__1_sm_ehs_B.Sum_n[1] =
      boost_and_two_level__1_sm_ehs_P.Iq_ref_Value -
      boost_and_two_level__1_sm_ehs_B.yk[1];

    /* Gain: '<S40>/Proportional Gain' */
    boost_and_two_level__1_sm_ehs_B.ProportionalGain[0] =
      boost_and_two_level__1_sm_ehs_P.InverterControl_Kp_Ireg *
      boost_and_two_level__1_sm_ehs_B.Sum_n[0];
    boost_and_two_level__1_sm_ehs_B.ProportionalGain[1] =
      boost_and_two_level__1_sm_ehs_P.InverterControl_Kp_Ireg *
      boost_and_two_level__1_sm_ehs_B.Sum_n[1];

    /* DiscreteIntegrator: '<S40>/Integrator' */
    boost_and_two_level__1_sm_ehs_B.Integrator[0] =
      boost_and_two_level__1_sm_ehs_DW.Integrator_DSTATE[0];
    boost_and_two_level__1_sm_ehs_B.Integrator[1] =
      boost_and_two_level__1_sm_ehs_DW.Integrator_DSTATE[1];

    /* Sum: '<S40>/Sum' */
    boost_and_two_level__1_sm_ehs_B.Sum_j[0] =
      boost_and_two_level__1_sm_ehs_B.ProportionalGain[0] +
      boost_and_two_level__1_sm_ehs_B.Integrator[0];
    boost_and_two_level__1_sm_ehs_B.Sum_j[1] =
      boost_and_two_level__1_sm_ehs_B.ProportionalGain[1] +
      boost_and_two_level__1_sm_ehs_B.Integrator[1];

    /* Saturate: '<S40>/Saturate' */
    B = boost_and_two_level__1_sm_ehs_B.Sum_j[0];
    dP = boost_and_two_level__1_sm_ehs_P.PI_LowerSaturationLimit;
    u2 = boost_and_two_level__1_sm_ehs_P.PI_UpperSaturationLimit;
    if (B > u2) {
      B = u2;
    } else {
      if (B < dP) {
        B = dP;
      }
    }

    boost_and_two_level__1_sm_ehs_B.Saturate[0] = B;
    B = boost_and_two_level__1_sm_ehs_B.Sum_j[1];
    dP = boost_and_two_level__1_sm_ehs_P.PI_LowerSaturationLimit;
    u2 = boost_and_two_level__1_sm_ehs_P.PI_UpperSaturationLimit;
    if (B > u2) {
      B = u2;
    } else {
      if (B < dP) {
        B = dP;
      }
    }

    boost_and_two_level__1_sm_ehs_B.Saturate[1] = B;

    /* End of Saturate: '<S40>/Saturate' */

    /* RateTransition: '<S9>/Rate Transition4' */
    if (rtmIsMajorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
      boost_and_two_level__1_sm_ehs_B.RateTransition4 =
        boost_and_two_level__1_sm_ehs_B.Divide[5];
    }

    /* End of RateTransition: '<S9>/Rate Transition4' */

    /* Gain: '<S36>/V->pu' */
    B = boost_and_two_level__1_sm_ehs_P.InverterControl_Vnom_prim *
      1.4142135623730951;
    dV = 1.0 / B;
    boost_and_two_level__1_sm_ehs_B.Vpu = dV *
      boost_and_two_level__1_sm_ehs_B.RateTransition4;

    /* Trigonometry: '<S41>/Trigonometric Function' */
    boost_and_two_level__1_sm_ehs_B.TrigonometricFunction = sin
      (boost_and_two_level__1_sm_ehs_B.MathFunction_m);

    /* Gain: '<S41>/Gain1' */
    boost_and_two_level__1_sm_ehs_B.Gain1_f =
      boost_and_two_level__1_sm_ehs_P.Gain1_Gain_b2 *
      boost_and_two_level__1_sm_ehs_B.TrigonometricFunction;

    /* Product: '<S41>/Product1' */
    boost_and_two_level__1_sm_ehs_B.Product1_k =
      boost_and_two_level__1_sm_ehs_B.Vpu *
      boost_and_two_level__1_sm_ehs_B.Gain1_f;

    /* DiscreteIntegrator: '<S48>/Integ4' */
    if (boost_and_two_level__1_sm_ehs_DW.Integ4_SYSTEM_ENABLE != 0) {
      boost_and_two_level__1_sm_ehs_B.Integ4 =
        boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE;
    } else {
      boost_and_two_level__1_sm_ehs_B.Integ4 =
        boost_and_two_level__1_sm_ehs_P.Integ4_gainval_l *
        boost_and_two_level__1_sm_ehs_B.Product1_k +
        boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE;
    }

    /* End of DiscreteIntegrator: '<S48>/Integ4' */

    /* Saturate: '<S48>/To avoid division  by zero' */
    B = boost_and_two_level__1_sm_ehs_B.UnitDelay_l;
    dP = boost_and_two_level__1_sm_ehs_P.Toavoiddivisionbyzero_LowerSa_o;
    u2 = boost_and_two_level__1_sm_ehs_P.Toavoiddivisionbyzero_UpperSa_d;
    if (B > u2) {
      boost_and_two_level__1_sm_ehs_B.Freq = u2;
    } else if (B < dP) {
      boost_and_two_level__1_sm_ehs_B.Freq = dP;
    } else {
      boost_and_two_level__1_sm_ehs_B.Freq = B;
    }

    /* End of Saturate: '<S48>/To avoid division  by zero' */

    /* Fcn: '<S48>/Number of samples per cycle' */
    boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle = 1.0 /
      boost_and_two_level__1_sm_ehs_B.Freq / 2.0e-5;

    /* Rounding: '<S48>/Rounding Function' */
    boost_and_two_level__1_sm_ehs_B.RoundingFunction = ceil
      (boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle);

    /* Gain: '<S48>/Gain' */
    boost_and_two_level__1_sm_ehs_B.Delay = boost_and_two_level__1_sm_ehs_P.Ts *
      boost_and_two_level__1_sm_ehs_B.RoundingFunction;

    /* Level2 S-Function Block: '<S50>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[26];
      sfcnOutputs(rts, 1);
    }

    /* UnitDelay: '<S49>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_B.UnitDelay_k =
      boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_d;

    /* DigitalClock: '<S48>/Digital  Clock' */
    boost_and_two_level__1_sm_ehs_B.DigitalClock_e =
      boost_and_two_level__1_sm_ehs_M->Timing.t[1];

    /* RelationalOperator: '<S48>/Relational Operator' incorporates:
     *  Constant: '<S48>/Constant'
     */
    boost_and_two_level__1_sm_ehs_B.RelationalOperator_m =
      (boost_and_two_level__1_sm_ehs_B.DigitalClock_e >=
       boost_and_two_level__1_sm_ehs_P.Constant_Value_p);

    /* UnitDelay: '<S48>/Unit Delay1' */
    boost_and_two_level__1_sm_ehs_B.UnitDelay1 =
      boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE;

    /* Switch: '<S48>/Switch' */
    if (boost_and_two_level__1_sm_ehs_B.RelationalOperator_m) {
      /* Sum: '<S49>/Sum1' */
      boost_and_two_level__1_sm_ehs_B.Sum1_g =
        boost_and_two_level__1_sm_ehs_B.Product1_k -
        boost_and_two_level__1_sm_ehs_B.UnitDelay_k;

      /* Sum: '<S49>/Sum5' */
      boost_and_two_level__1_sm_ehs_B.Sum5_f =
        boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle -
        boost_and_two_level__1_sm_ehs_B.RoundingFunction;

      /* Product: '<S49>/Product5' */
      boost_and_two_level__1_sm_ehs_B.Product5_l =
        boost_and_two_level__1_sm_ehs_B.Sum5_f *
        boost_and_two_level__1_sm_ehs_B.Sum1_g;

      /* Gain: '<S49>/Gain1' */
      boost_and_two_level__1_sm_ehs_B.Gain1_pu =
        boost_and_two_level__1_sm_ehs_P.Gain1_Gain_jj *
        boost_and_two_level__1_sm_ehs_B.Product5_l;

      /* Sum: '<S49>/Sum4' */
      boost_and_two_level__1_sm_ehs_B.Sum4_f =
        boost_and_two_level__1_sm_ehs_B.Gain1_pu +
        boost_and_two_level__1_sm_ehs_B.Product1_k;

      /* Product: '<S49>/Product2' */
      boost_and_two_level__1_sm_ehs_B.Product2_f1 =
        boost_and_two_level__1_sm_ehs_B.Sum5_f /
        boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle;

      /* Product: '<S49>/Product4' */
      boost_and_two_level__1_sm_ehs_B.Product4_k =
        boost_and_two_level__1_sm_ehs_B.Product2_f1 *
        boost_and_two_level__1_sm_ehs_B.Sum4_f;

      /* Sum: '<S48>/Sum7' */
      boost_and_two_level__1_sm_ehs_B.Sum7_p =
        boost_and_two_level__1_sm_ehs_B.Integ4 -
        boost_and_two_level__1_sm_ehs_B.SFunction_f;

      /* Product: '<S48>/Product' */
      boost_and_two_level__1_sm_ehs_B.Meanvalue_o =
        boost_and_two_level__1_sm_ehs_B.Sum7_p *
        boost_and_two_level__1_sm_ehs_B.UnitDelay_l;

      /* Sum: '<S48>/Sum5' */
      boost_and_two_level__1_sm_ehs_B.Sum5_n4 =
        boost_and_two_level__1_sm_ehs_B.Meanvalue_o +
        boost_and_two_level__1_sm_ehs_B.Product4_k;
      boost_and_two_level__1_sm_ehs_B.Switch_ie =
        boost_and_two_level__1_sm_ehs_B.Sum5_n4;
    } else {
      boost_and_two_level__1_sm_ehs_B.Switch_ie =
        boost_and_two_level__1_sm_ehs_B.UnitDelay1;
    }

    /* End of Switch: '<S48>/Switch' */

    /* Trigonometry: '<S41>/Trigonometric Function3' */
    boost_and_two_level__1_sm_ehs_B.TrigonometricFunction3 = cos
      (boost_and_two_level__1_sm_ehs_B.MathFunction_m);

    /* Gain: '<S41>/Gain3' */
    boost_and_two_level__1_sm_ehs_B.Gain3 =
      boost_and_two_level__1_sm_ehs_P.Gain3_Gain_k *
      boost_and_two_level__1_sm_ehs_B.TrigonometricFunction3;

    /* Product: '<S41>/Product2' */
    boost_and_two_level__1_sm_ehs_B.Product2 =
      boost_and_two_level__1_sm_ehs_B.Vpu *
      boost_and_two_level__1_sm_ehs_B.Gain3;

    /* DiscreteIntegrator: '<S51>/Integ4' */
    if (boost_and_two_level__1_sm_ehs_DW.Integ4_SYSTEM_ENABLE_n != 0) {
      boost_and_two_level__1_sm_ehs_B.Integ4_b =
        boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_f;
    } else {
      boost_and_two_level__1_sm_ehs_B.Integ4_b =
        boost_and_two_level__1_sm_ehs_P.Integ4_gainval_g *
        boost_and_two_level__1_sm_ehs_B.Product2 +
        boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_f;
    }

    /* End of DiscreteIntegrator: '<S51>/Integ4' */

    /* Saturate: '<S51>/To avoid division  by zero' */
    B = boost_and_two_level__1_sm_ehs_B.UnitDelay_l;
    dP = boost_and_two_level__1_sm_ehs_P.Toavoiddivisionbyzero_LowerS_h5;
    u2 = boost_and_two_level__1_sm_ehs_P.Toavoiddivisionbyzero_UpperSa_e;
    if (B > u2) {
      boost_and_two_level__1_sm_ehs_B.Freq_g = u2;
    } else if (B < dP) {
      boost_and_two_level__1_sm_ehs_B.Freq_g = dP;
    } else {
      boost_and_two_level__1_sm_ehs_B.Freq_g = B;
    }

    /* End of Saturate: '<S51>/To avoid division  by zero' */

    /* Fcn: '<S51>/Number of samples per cycle' */
    boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_h = 1.0 /
      boost_and_two_level__1_sm_ehs_B.Freq_g / 2.0e-5;

    /* Rounding: '<S51>/Rounding Function' */
    boost_and_two_level__1_sm_ehs_B.RoundingFunction_e = ceil
      (boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_h);

    /* Gain: '<S51>/Gain' */
    boost_and_two_level__1_sm_ehs_B.Delay_p = boost_and_two_level__1_sm_ehs_P.Ts
      * boost_and_two_level__1_sm_ehs_B.RoundingFunction_e;

    /* Level2 S-Function Block: '<S53>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[27];
      sfcnOutputs(rts, 1);
    }

    /* UnitDelay: '<S52>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_B.UnitDelay_h =
      boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_i;

    /* DigitalClock: '<S51>/Digital  Clock' */
    boost_and_two_level__1_sm_ehs_B.DigitalClock_i =
      boost_and_two_level__1_sm_ehs_M->Timing.t[1];

    /* RelationalOperator: '<S51>/Relational Operator' incorporates:
     *  Constant: '<S51>/Constant'
     */
    boost_and_two_level__1_sm_ehs_B.RelationalOperator_a =
      (boost_and_two_level__1_sm_ehs_B.DigitalClock_i >=
       boost_and_two_level__1_sm_ehs_P.Constant_Value_d);

    /* UnitDelay: '<S51>/Unit Delay1' */
    boost_and_two_level__1_sm_ehs_B.UnitDelay1_c =
      boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_i;

    /* Switch: '<S51>/Switch' */
    if (boost_and_two_level__1_sm_ehs_B.RelationalOperator_a) {
      /* Sum: '<S52>/Sum1' */
      boost_and_two_level__1_sm_ehs_B.Sum1_e =
        boost_and_two_level__1_sm_ehs_B.Product2 -
        boost_and_two_level__1_sm_ehs_B.UnitDelay_h;

      /* Sum: '<S52>/Sum5' */
      boost_and_two_level__1_sm_ehs_B.Sum5_p =
        boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_h -
        boost_and_two_level__1_sm_ehs_B.RoundingFunction_e;

      /* Product: '<S52>/Product5' */
      boost_and_two_level__1_sm_ehs_B.Product5_i4 =
        boost_and_two_level__1_sm_ehs_B.Sum5_p *
        boost_and_two_level__1_sm_ehs_B.Sum1_e;

      /* Gain: '<S52>/Gain1' */
      boost_and_two_level__1_sm_ehs_B.Gain1_i =
        boost_and_two_level__1_sm_ehs_P.Gain1_Gain_jb *
        boost_and_two_level__1_sm_ehs_B.Product5_i4;

      /* Sum: '<S52>/Sum4' */
      boost_and_two_level__1_sm_ehs_B.Sum4_l =
        boost_and_two_level__1_sm_ehs_B.Gain1_i +
        boost_and_two_level__1_sm_ehs_B.Product2;

      /* Product: '<S52>/Product2' */
      boost_and_two_level__1_sm_ehs_B.Product2_gq =
        boost_and_two_level__1_sm_ehs_B.Sum5_p /
        boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_h;

      /* Product: '<S52>/Product4' */
      boost_and_two_level__1_sm_ehs_B.Product4_n =
        boost_and_two_level__1_sm_ehs_B.Product2_gq *
        boost_and_two_level__1_sm_ehs_B.Sum4_l;

      /* Sum: '<S51>/Sum7' */
      boost_and_two_level__1_sm_ehs_B.Sum7_m =
        boost_and_two_level__1_sm_ehs_B.Integ4_b -
        boost_and_two_level__1_sm_ehs_B.SFunction_fz;

      /* Product: '<S51>/Product' */
      boost_and_two_level__1_sm_ehs_B.Meanvalue_a =
        boost_and_two_level__1_sm_ehs_B.Sum7_m *
        boost_and_two_level__1_sm_ehs_B.UnitDelay_l;

      /* Sum: '<S51>/Sum5' */
      boost_and_two_level__1_sm_ehs_B.Sum5_an =
        boost_and_two_level__1_sm_ehs_B.Meanvalue_a +
        boost_and_two_level__1_sm_ehs_B.Product4_n;
      boost_and_two_level__1_sm_ehs_B.Switch_c =
        boost_and_two_level__1_sm_ehs_B.Sum5_an;
    } else {
      boost_and_two_level__1_sm_ehs_B.Switch_c =
        boost_and_two_level__1_sm_ehs_B.UnitDelay1_c;
    }

    /* End of Switch: '<S51>/Switch' */

    /* RealImagToComplex: '<S41>/Real-Imag to Complex' */
    boost_and_two_level__1_sm_ehs_B.RealImagtoComplex.re =
      boost_and_two_level__1_sm_ehs_B.Switch_ie;
    boost_and_two_level__1_sm_ehs_B.RealImagtoComplex.im =
      boost_and_two_level__1_sm_ehs_B.Switch_c;

    /* ComplexToMagnitudeAngle: '<S41>/Complex to Magnitude-Angle' */
    boost_and_two_level__1_sm_ehs_B.ComplextoMagnitudeAngle_o1 = rt_hypotd_snf
      (boost_and_two_level__1_sm_ehs_B.RealImagtoComplex.re,
       boost_and_two_level__1_sm_ehs_B.RealImagtoComplex.im);
    boost_and_two_level__1_sm_ehs_B.ComplextoMagnitudeAngle_o2 = rt_atan2d_snf
      (boost_and_two_level__1_sm_ehs_B.RealImagtoComplex.im,
       boost_and_two_level__1_sm_ehs_B.RealImagtoComplex.re);

    /* Gain: '<S41>/Rad->Deg.' */
    boost_and_two_level__1_sm_ehs_B.RadDeg =
      boost_and_two_level__1_sm_ehs_P.RadDeg_Gain_n *
      boost_and_two_level__1_sm_ehs_B.ComplextoMagnitudeAngle_o2;

    /* Gain: '<S36>/to-rad' */
    boost_and_two_level__1_sm_ehs_B.torad =
      boost_and_two_level__1_sm_ehs_P.torad_Gain *
      boost_and_two_level__1_sm_ehs_B.RadDeg;

    /* MagnitudeAngleToComplex: '<S36>/Magnitude-Angle to Complex' */
    P = boost_and_two_level__1_sm_ehs_B.torad;
    dV = boost_and_two_level__1_sm_ehs_B.ComplextoMagnitudeAngle_o1;
    u2 = cos(P);
    B = sin(P);
    boost_and_two_level__1_sm_ehs_B.MagnitudeAngletoComplex.re = dV * u2;
    boost_and_two_level__1_sm_ehs_B.MagnitudeAngletoComplex.im = dV * B;

    /* ComplexToRealImag: '<S36>/Complex to Real-Imag' */
    boost_and_two_level__1_sm_ehs_B.ComplextoRealImag_o1 =
      boost_and_two_level__1_sm_ehs_B.MagnitudeAngletoComplex.re;
    boost_and_two_level__1_sm_ehs_B.ComplextoRealImag_o2 =
      boost_and_two_level__1_sm_ehs_B.MagnitudeAngletoComplex.im;

    /* Gain: '<S33>/Rff ' */
    boost_and_two_level__1_sm_ehs_B.Rff =
      boost_and_two_level__1_sm_ehs_P.Rff_Gain *
      boost_and_two_level__1_sm_ehs_B.UnitDelay2;

    /* Gain: '<S33>/Lff  ' incorporates:
     *  Constant: '<S31>/Iq_ref'
     */
    boost_and_two_level__1_sm_ehs_B.Lff =
      boost_and_two_level__1_sm_ehs_P.Lff_Gain *
      boost_and_two_level__1_sm_ehs_P.Iq_ref_Value;

    /* Sum: '<S33>/Add1' */
    boost_and_two_level__1_sm_ehs_B.Feedforward =
      (boost_and_two_level__1_sm_ehs_B.ComplextoRealImag_o1 +
       boost_and_two_level__1_sm_ehs_B.Rff) -
      boost_and_two_level__1_sm_ehs_B.Lff;

    /* Gain: '<S33>/Rff' incorporates:
     *  Constant: '<S31>/Iq_ref'
     */
    boost_and_two_level__1_sm_ehs_B.Rff_d =
      boost_and_two_level__1_sm_ehs_P.Rff_Gain_o *
      boost_and_two_level__1_sm_ehs_P.Iq_ref_Value;

    /* Gain: '<S33>/Lff' */
    boost_and_two_level__1_sm_ehs_B.Lff_o =
      boost_and_two_level__1_sm_ehs_P.Lff_Gain_b *
      boost_and_two_level__1_sm_ehs_B.UnitDelay2;

    /* Sum: '<S33>/Add3' */
    boost_and_two_level__1_sm_ehs_B.Add3_i =
      (boost_and_two_level__1_sm_ehs_B.ComplextoRealImag_o2 +
       boost_and_two_level__1_sm_ehs_B.Rff_d) +
      boost_and_two_level__1_sm_ehs_B.Lff_o;

    /* Sum: '<S33>/Add2' */
    boost_and_two_level__1_sm_ehs_B.Add2[0] =
      boost_and_two_level__1_sm_ehs_B.Saturate[0] +
      boost_and_two_level__1_sm_ehs_B.Feedforward;
    boost_and_two_level__1_sm_ehs_B.Add2[1] =
      boost_and_two_level__1_sm_ehs_B.Saturate[1] +
      boost_and_two_level__1_sm_ehs_B.Add3_i;

    /* Gain: '<S40>/Integral Gain' */
    boost_and_two_level__1_sm_ehs_B.IntegralGain[0] =
      boost_and_two_level__1_sm_ehs_P.InverterControl_Ki_Ireg *
      boost_and_two_level__1_sm_ehs_B.Sum_n[0];
    boost_and_two_level__1_sm_ehs_B.IntegralGain[1] =
      boost_and_two_level__1_sm_ehs_P.InverterControl_Ki_Ireg *
      boost_and_two_level__1_sm_ehs_B.Sum_n[1];

    /* Saturate: '<S33>/Saturation' */
    B = boost_and_two_level__1_sm_ehs_B.Add2[0];
    dP = boost_and_two_level__1_sm_ehs_P.Saturation_LowerSat_e;
    u2 = boost_and_two_level__1_sm_ehs_P.Saturation_UpperSat_c;
    if (B > u2) {
      B = u2;
    } else {
      if (B < dP) {
        B = dP;
      }
    }

    boost_and_two_level__1_sm_ehs_B.Saturation[0] = B;
    B = boost_and_two_level__1_sm_ehs_B.Add2[1];
    dP = boost_and_two_level__1_sm_ehs_P.Saturation_LowerSat_e;
    u2 = boost_and_two_level__1_sm_ehs_P.Saturation_UpperSat_c;
    if (B > u2) {
      B = u2;
    } else {
      if (B < dP) {
        B = dP;
      }
    }

    boost_and_two_level__1_sm_ehs_B.Saturation[1] = B;

    /* End of Saturate: '<S33>/Saturation' */

    /* RateTransition: '<S9>/Rate Transition7' incorporates:
     *  RateTransition: '<S9>/Rate Transition8'
     */
    if (rtmIsMajorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
      boost_and_two_level__1_sm_ehs_B.RateTransition7 =
        boost_and_two_level__1_sm_ehs_B.Divide[1];
      boost_and_two_level__1_sm_ehs_B.RateTransition8 =
        boost_and_two_level__1_sm_ehs_B.Divide[0];
    }

    /* End of RateTransition: '<S9>/Rate Transition7' */

    /* SignalConversion: '<S35>/TmpSignal ConversionAt SFunction Inport1' incorporates:
     *  Constant: '<S34>/Iph_'
     *  Constant: '<S34>/Iph_1'
     *  Constant: '<S34>/Iph_2'
     *  Constant: '<S34>/Iph_3'
     *  MATLAB Function: '<S31>/MPPT Controller using Perturbe  & Observe technique  '
     */
    boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtSFunctionI[0] =
      boost_and_two_level__1_sm_ehs_P.InverterControl_Vdc_ref_Init;
    boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtSFunctionI[1] =
      boost_and_two_level__1_sm_ehs_P.Iph_1_Value;
    boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtSFunctionI[2] =
      boost_and_two_level__1_sm_ehs_P.Iph_2_Value;
    boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtSFunctionI[3] =
      boost_and_two_level__1_sm_ehs_P.InverterControl_Increment_MPPT;

    /* MATLAB Function: '<S31>/MPPT Controller using Perturbe  & Observe technique  ' incorporates:
     *  Constant: '<S9>/MPPT_On'
     */
    /* MATLAB Function 'SM_eHS/source_and_gates/Inverter Control/MPPT Controller using Perturbe  & Observe technique  ': '<S35>:1' */
    /*  MPPT controller based on the Perturb & Observe algorithm. */
    /*  D output = Reference for DC link voltage (Vdc_ref)  */
    /*  */
    /*  Enabled input = 1 to enable the MPPT controller */
    /*  V input = PV array terminal voltage (V) */
    /*  I input = PV array current (A) */
    /*  */
    /*  Param input: */
    /* Initial value for Vdc_ref */
    /* '<S35>:1:13' */
    /* Maximum value for Vdc_ref */
    /* '<S35>:1:14' */
    /* Minimum value for Vdc_ref */
    /* '<S35>:1:15' */
    /* Increment value used to increase/decrease Vdc_ref */
    /*    */
    if (!boost_and_two_level__1_sm_ehs_DW.Vold_not_empty) {
      /* '<S35>:1:22' */
      boost_and_two_level__1_sm_ehs_DW.Vold_not_empty = true;

      /* '<S35>:1:25' */
      boost_and_two_level__1_sm_ehs_DW.Dold =
        boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtSFunctionI[0];
    }

    /* '<S35>:1:27' */
    P = boost_and_two_level__1_sm_ehs_B.RateTransition7 *
      boost_and_two_level__1_sm_ehs_B.RateTransition8;

    /* '<S35>:1:28' */
    dV = boost_and_two_level__1_sm_ehs_B.RateTransition7 -
      boost_and_two_level__1_sm_ehs_DW.Vold;

    /* '<S35>:1:29' */
    dP = P - boost_and_two_level__1_sm_ehs_DW.Pold;
    if ((dP != 0.0) && (boost_and_two_level__1_sm_ehs_P.MPPT_On_Value != 0.0)) {
      /* '<S35>:1:31' */
      if (dP < 0.0) {
        /* '<S35>:1:32' */
        if (dV < 0.0) {
          /* '<S35>:1:33' */
          /* '<S35>:1:34' */
          boost_and_two_level__1_sm_ehs_B.D =
            boost_and_two_level__1_sm_ehs_DW.Dold +
            boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtSFunctionI[3];
        } else {
          /* '<S35>:1:36' */
          boost_and_two_level__1_sm_ehs_B.D =
            boost_and_two_level__1_sm_ehs_DW.Dold -
            boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtSFunctionI[3];
        }
      } else if (dV < 0.0) {
        /* '<S35>:1:39' */
        /* '<S35>:1:40' */
        boost_and_two_level__1_sm_ehs_B.D =
          boost_and_two_level__1_sm_ehs_DW.Dold -
          boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtSFunctionI[3];
      } else {
        /* '<S35>:1:42' */
        boost_and_two_level__1_sm_ehs_B.D =
          boost_and_two_level__1_sm_ehs_DW.Dold +
          boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtSFunctionI[3];
      }
    } else {
      /* '<S35>:1:45' */
      boost_and_two_level__1_sm_ehs_B.D = boost_and_two_level__1_sm_ehs_DW.Dold;
    }

    if ((boost_and_two_level__1_sm_ehs_B.D >=
         boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtSFunctionI[1]) ||
        (boost_and_two_level__1_sm_ehs_B.D <=
         boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtSFunctionI[2])) {
      /* '<S35>:1:48' */
      /* '<S35>:1:49' */
      boost_and_two_level__1_sm_ehs_B.D = boost_and_two_level__1_sm_ehs_DW.Dold;
    }

    /* '<S35>:1:52' */
    boost_and_two_level__1_sm_ehs_DW.Dold = boost_and_two_level__1_sm_ehs_B.D;

    /* '<S35>:1:53' */
    boost_and_two_level__1_sm_ehs_DW.Vold =
      boost_and_two_level__1_sm_ehs_B.RateTransition7;

    /* '<S35>:1:54' */
    boost_and_two_level__1_sm_ehs_DW.Pold = P;

    /* DigitalClock: '<S54>/Digital  Clock' */
    boost_and_two_level__1_sm_ehs_B.DigitalClock_j =
      boost_and_two_level__1_sm_ehs_M->Timing.t[1];

    /* RateTransition: '<S9>/Rate Transition6' */
    if (rtmIsMajorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
      boost_and_two_level__1_sm_ehs_B.RateTransition6 =
        boost_and_two_level__1_sm_ehs_B.Divide[1];
    }

    /* End of RateTransition: '<S9>/Rate Transition6' */

    /* DiscreteIntegrator: '<S54>/Integ4' */
    if (boost_and_two_level__1_sm_ehs_DW.Integ4_SYSTEM_ENABLE_i != 0) {
      boost_and_two_level__1_sm_ehs_B.Integ4_m =
        boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_g;
    } else {
      boost_and_two_level__1_sm_ehs_B.Integ4_m =
        boost_and_two_level__1_sm_ehs_P.Integ4_gainval_o *
        boost_and_two_level__1_sm_ehs_B.RateTransition6 +
        boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_g;
    }

    /* End of DiscreteIntegrator: '<S54>/Integ4' */

    /* Constant: '<S54>/K1' */
    boost_and_two_level__1_sm_ehs_B.K1 =
      boost_and_two_level__1_sm_ehs_P.K1_Value_f;

    /* Level2 S-Function Block: '<S55>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[28];
      sfcnOutputs(rts, 1);
    }

    /* UnitDelay: '<S54>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_B.UnitDelay_b =
      boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_h;

    /* RelationalOperator: '<S54>/Relational Operator' */
    boost_and_two_level__1_sm_ehs_B.RelationalOperator_h =
      (boost_and_two_level__1_sm_ehs_B.DigitalClock_j >=
       boost_and_two_level__1_sm_ehs_B.K1);

    /* UnitDelay: '<S54>/Unit Delay1' */
    boost_and_two_level__1_sm_ehs_B.UnitDelay1_cm =
      boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_im;

    /* Switch: '<S54>/Switch' */
    if (boost_and_two_level__1_sm_ehs_B.RelationalOperator_h) {
      /* Gain: '<S54>/Gain1' */
      boost_and_two_level__1_sm_ehs_B.Gain1_a =
        boost_and_two_level__1_sm_ehs_P.Gain1_Gain_e *
        boost_and_two_level__1_sm_ehs_B.RateTransition6;

      /* Gain: '<S54>/Gain' */
      boost_and_two_level__1_sm_ehs_B.Gain_o =
        boost_and_two_level__1_sm_ehs_P.Gain_Gain_j *
        boost_and_two_level__1_sm_ehs_B.UnitDelay_b;

      /* Sum: '<S54>/Sum1' */
      boost_and_two_level__1_sm_ehs_B.Correction =
        boost_and_two_level__1_sm_ehs_B.Gain1_a -
        boost_and_two_level__1_sm_ehs_B.Gain_o;

      /* Sum: '<S54>/Sum7' */
      boost_and_two_level__1_sm_ehs_B.Sum7_js =
        boost_and_two_level__1_sm_ehs_B.Integ4_m -
        boost_and_two_level__1_sm_ehs_B.SFunction_c;

      /* Product: '<S54>/Product' incorporates:
       *  Constant: '<S54>/K2'
       */
      boost_and_two_level__1_sm_ehs_B.Mean =
        boost_and_two_level__1_sm_ehs_B.Sum7_js *
        boost_and_two_level__1_sm_ehs_P.K2_Value;

      /* Sum: '<S54>/Sum5' */
      boost_and_two_level__1_sm_ehs_B.Sum5_c =
        boost_and_two_level__1_sm_ehs_B.Mean +
        boost_and_two_level__1_sm_ehs_B.Correction;
      boost_and_two_level__1_sm_ehs_B.Switch_ej =
        boost_and_two_level__1_sm_ehs_B.Sum5_c;
    } else {
      boost_and_two_level__1_sm_ehs_B.Switch_ej =
        boost_and_two_level__1_sm_ehs_B.UnitDelay1_cm;
    }

    /* End of Switch: '<S54>/Switch' */

    /* Outputs for Enabled SubSystem: '<S56>/Automatic Gain Control' incorporates:
     *  EnablePort: '<S57>/Enable'
     */
    if (rtmIsMajorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
      /* Constant: '<S56>/Constant1' */
      if (boost_and_two_level__1_sm_ehs_P.PLL_AGC > 0.0) {
        if (!boost_and_two_level__1_sm_ehs_DW.AutomaticGainControl_MODE) {
          /* Enable for DiscreteIntegrator: '<S64>/Integ4' */
          boost_and_two_level__1_sm_ehs_DW.Integ4_SYSTEM_ENABLE_nz = 1U;

          /* Enable for DiscreteIntegrator: '<S67>/Integ4' */
          boost_and_two_level__1_sm_ehs_DW.Integ4_SYSTEM_ENABLE_d = 1U;
          boost_and_two_level__1_sm_ehs_DW.AutomaticGainControl_MODE = true;
        }
      } else {
        if (boost_and_two_level__1_sm_ehs_DW.AutomaticGainControl_MODE) {
          boost_and_two_level__1_sm_ehs_DW.AutomaticGainControl_MODE = false;
        }
      }

      /* End of Constant: '<S56>/Constant1' */
    }

    if (boost_and_two_level__1_sm_ehs_DW.AutomaticGainControl_MODE) {
      /* Trigonometry: '<S61>/Trigonometric Function' */
      boost_and_two_level__1_sm_ehs_B.TrigonometricFunction_e = sin
        (boost_and_two_level__1_sm_ehs_B.MathFunction_m);

      /* Gain: '<S61>/Gain1' */
      boost_and_two_level__1_sm_ehs_B.Gain1_fj =
        boost_and_two_level__1_sm_ehs_P.Gain1_Gain_o *
        boost_and_two_level__1_sm_ehs_B.TrigonometricFunction_e;

      /* Product: '<S61>/Product1' */
      boost_and_two_level__1_sm_ehs_B.Product1_a =
        boost_and_two_level__1_sm_ehs_B.Vpu *
        boost_and_two_level__1_sm_ehs_B.Gain1_fj;

      /* DiscreteIntegrator: '<S64>/Integ4' */
      if (boost_and_two_level__1_sm_ehs_DW.Integ4_SYSTEM_ENABLE_nz != 0) {
        boost_and_two_level__1_sm_ehs_B.Integ4_ee =
          boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_i;
      } else {
        boost_and_two_level__1_sm_ehs_B.Integ4_ee =
          boost_and_two_level__1_sm_ehs_P.Integ4_gainval *
          boost_and_two_level__1_sm_ehs_B.Product1_a +
          boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_i;
      }

      /* End of DiscreteIntegrator: '<S64>/Integ4' */

      /* Saturate: '<S64>/To avoid division  by zero' */
      B = boost_and_two_level__1_sm_ehs_B.UnitDelay_l;
      dP = boost_and_two_level__1_sm_ehs_P.Toavoiddivisionbyzero_LowerSat;
      u2 = boost_and_two_level__1_sm_ehs_P.Toavoiddivisionbyzero_UpperSat;
      if (B > u2) {
        boost_and_two_level__1_sm_ehs_B.Freq_h = u2;
      } else if (B < dP) {
        boost_and_two_level__1_sm_ehs_B.Freq_h = dP;
      } else {
        boost_and_two_level__1_sm_ehs_B.Freq_h = B;
      }

      /* End of Saturate: '<S64>/To avoid division  by zero' */

      /* Fcn: '<S64>/Number of samples per cycle' */
      boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_o = 1.0 /
        boost_and_two_level__1_sm_ehs_B.Freq_h / 2.0e-5;

      /* Rounding: '<S64>/Rounding Function' */
      boost_and_two_level__1_sm_ehs_B.RoundingFunction_j = ceil
        (boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_o);

      /* Gain: '<S64>/Gain' */
      boost_and_two_level__1_sm_ehs_B.Delay_o =
        boost_and_two_level__1_sm_ehs_P.Ts *
        boost_and_two_level__1_sm_ehs_B.RoundingFunction_j;

      /* Level2 S-Function Block: '<S66>/S-Function' (sfun_discreteVariableDelay) */
      {
        SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[0];
        sfcnOutputs(rts, 1);
      }

      /* UnitDelay: '<S65>/Unit Delay' */
      boost_and_two_level__1_sm_ehs_B.UnitDelay_e =
        boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_c;

      /* DigitalClock: '<S64>/Digital  Clock' */
      boost_and_two_level__1_sm_ehs_B.DigitalClock_g =
        boost_and_two_level__1_sm_ehs_M->Timing.t[1];

      /* RelationalOperator: '<S64>/Relational Operator' incorporates:
       *  Constant: '<S64>/Constant'
       */
      boost_and_two_level__1_sm_ehs_B.RelationalOperator_mw =
        (boost_and_two_level__1_sm_ehs_B.DigitalClock_g >=
         boost_and_two_level__1_sm_ehs_P.Constant_Value_l);

      /* UnitDelay: '<S64>/Unit Delay1' */
      boost_and_two_level__1_sm_ehs_B.UnitDelay1_a =
        boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_eb;

      /* Switch: '<S64>/Switch' */
      if (boost_and_two_level__1_sm_ehs_B.RelationalOperator_mw) {
        /* Sum: '<S65>/Sum1' */
        boost_and_two_level__1_sm_ehs_B.Sum1_i =
          boost_and_two_level__1_sm_ehs_B.Product1_a -
          boost_and_two_level__1_sm_ehs_B.UnitDelay_e;

        /* Sum: '<S65>/Sum5' */
        boost_and_two_level__1_sm_ehs_B.Sum5_a =
          boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_o -
          boost_and_two_level__1_sm_ehs_B.RoundingFunction_j;

        /* Product: '<S65>/Product5' */
        boost_and_two_level__1_sm_ehs_B.Product5_i =
          boost_and_two_level__1_sm_ehs_B.Sum5_a *
          boost_and_two_level__1_sm_ehs_B.Sum1_i;

        /* Gain: '<S65>/Gain1' */
        boost_and_two_level__1_sm_ehs_B.Gain1_p5 =
          boost_and_two_level__1_sm_ehs_P.Gain1_Gain_b *
          boost_and_two_level__1_sm_ehs_B.Product5_i;

        /* Sum: '<S65>/Sum4' */
        boost_and_two_level__1_sm_ehs_B.Sum4_n =
          boost_and_two_level__1_sm_ehs_B.Gain1_p5 +
          boost_and_two_level__1_sm_ehs_B.Product1_a;

        /* Product: '<S65>/Product2' */
        boost_and_two_level__1_sm_ehs_B.Product2_c =
          boost_and_two_level__1_sm_ehs_B.Sum5_a /
          boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_o;

        /* Product: '<S65>/Product4' */
        boost_and_two_level__1_sm_ehs_B.Product4_p =
          boost_and_two_level__1_sm_ehs_B.Product2_c *
          boost_and_two_level__1_sm_ehs_B.Sum4_n;

        /* Sum: '<S64>/Sum7' */
        boost_and_two_level__1_sm_ehs_B.Sum7_j =
          boost_and_two_level__1_sm_ehs_B.Integ4_ee -
          boost_and_two_level__1_sm_ehs_B.SFunction_b;

        /* Product: '<S64>/Product' */
        boost_and_two_level__1_sm_ehs_B.Meanvalue_h =
          boost_and_two_level__1_sm_ehs_B.Sum7_j *
          boost_and_two_level__1_sm_ehs_B.UnitDelay_l;

        /* Sum: '<S64>/Sum5' */
        boost_and_two_level__1_sm_ehs_B.Sum5_d =
          boost_and_two_level__1_sm_ehs_B.Meanvalue_h +
          boost_and_two_level__1_sm_ehs_B.Product4_p;
        boost_and_two_level__1_sm_ehs_B.Switch_h =
          boost_and_two_level__1_sm_ehs_B.Sum5_d;
      } else {
        boost_and_two_level__1_sm_ehs_B.Switch_h =
          boost_and_two_level__1_sm_ehs_B.UnitDelay1_a;
      }

      /* End of Switch: '<S64>/Switch' */

      /* Trigonometry: '<S61>/Trigonometric Function3' */
      boost_and_two_level__1_sm_ehs_B.TrigonometricFunction3_k = cos
        (boost_and_two_level__1_sm_ehs_B.MathFunction_m);

      /* Gain: '<S61>/Gain3' */
      boost_and_two_level__1_sm_ehs_B.Gain3_b =
        boost_and_two_level__1_sm_ehs_P.Gain3_Gain *
        boost_and_two_level__1_sm_ehs_B.TrigonometricFunction3_k;

      /* Product: '<S61>/Product2' */
      boost_and_two_level__1_sm_ehs_B.Product2_b =
        boost_and_two_level__1_sm_ehs_B.Vpu *
        boost_and_two_level__1_sm_ehs_B.Gain3_b;

      /* DiscreteIntegrator: '<S67>/Integ4' */
      if (boost_and_two_level__1_sm_ehs_DW.Integ4_SYSTEM_ENABLE_d != 0) {
        boost_and_two_level__1_sm_ehs_B.Integ4_i =
          boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_c;
      } else {
        boost_and_two_level__1_sm_ehs_B.Integ4_i =
          boost_and_two_level__1_sm_ehs_P.Integ4_gainval_c *
          boost_and_two_level__1_sm_ehs_B.Product2_b +
          boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_c;
      }

      /* End of DiscreteIntegrator: '<S67>/Integ4' */

      /* Saturate: '<S67>/To avoid division  by zero' */
      B = boost_and_two_level__1_sm_ehs_B.UnitDelay_l;
      dP = boost_and_two_level__1_sm_ehs_P.Toavoiddivisionbyzero_LowerSa_h;
      u2 = boost_and_two_level__1_sm_ehs_P.Toavoiddivisionbyzero_UpperSa_i;
      if (B > u2) {
        boost_and_two_level__1_sm_ehs_B.Freq_c = u2;
      } else if (B < dP) {
        boost_and_two_level__1_sm_ehs_B.Freq_c = dP;
      } else {
        boost_and_two_level__1_sm_ehs_B.Freq_c = B;
      }

      /* End of Saturate: '<S67>/To avoid division  by zero' */

      /* Fcn: '<S67>/Number of samples per cycle' */
      boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_e = 1.0 /
        boost_and_two_level__1_sm_ehs_B.Freq_c / 2.0e-5;

      /* Rounding: '<S67>/Rounding Function' */
      boost_and_two_level__1_sm_ehs_B.RoundingFunction_n = ceil
        (boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_e);

      /* Gain: '<S67>/Gain' */
      boost_and_two_level__1_sm_ehs_B.Delay_e =
        boost_and_two_level__1_sm_ehs_P.Ts *
        boost_and_two_level__1_sm_ehs_B.RoundingFunction_n;

      /* Level2 S-Function Block: '<S69>/S-Function' (sfun_discreteVariableDelay) */
      {
        SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[1];
        sfcnOutputs(rts, 1);
      }

      /* UnitDelay: '<S68>/Unit Delay' */
      boost_and_two_level__1_sm_ehs_B.UnitDelay_m =
        boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_fv;

      /* DigitalClock: '<S67>/Digital  Clock' */
      boost_and_two_level__1_sm_ehs_B.DigitalClock_d =
        boost_and_two_level__1_sm_ehs_M->Timing.t[1];

      /* RelationalOperator: '<S67>/Relational Operator' incorporates:
       *  Constant: '<S67>/Constant'
       */
      boost_and_two_level__1_sm_ehs_B.RelationalOperator_am =
        (boost_and_two_level__1_sm_ehs_B.DigitalClock_d >=
         boost_and_two_level__1_sm_ehs_P.Constant_Value_o);

      /* UnitDelay: '<S67>/Unit Delay1' */
      boost_and_two_level__1_sm_ehs_B.UnitDelay1_n =
        boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_k;

      /* Switch: '<S67>/Switch' */
      if (boost_and_two_level__1_sm_ehs_B.RelationalOperator_am) {
        /* Sum: '<S68>/Sum1' */
        boost_and_two_level__1_sm_ehs_B.Sum1_h =
          boost_and_two_level__1_sm_ehs_B.Product2_b -
          boost_and_two_level__1_sm_ehs_B.UnitDelay_m;

        /* Sum: '<S68>/Sum5' */
        boost_and_two_level__1_sm_ehs_B.Sum5_n =
          boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_e -
          boost_and_two_level__1_sm_ehs_B.RoundingFunction_n;

        /* Product: '<S68>/Product5' */
        boost_and_two_level__1_sm_ehs_B.Product5_j =
          boost_and_two_level__1_sm_ehs_B.Sum5_n *
          boost_and_two_level__1_sm_ehs_B.Sum1_h;

        /* Gain: '<S68>/Gain1' */
        boost_and_two_level__1_sm_ehs_B.Gain1_b =
          boost_and_two_level__1_sm_ehs_P.Gain1_Gain_e3 *
          boost_and_two_level__1_sm_ehs_B.Product5_j;

        /* Sum: '<S68>/Sum4' */
        boost_and_two_level__1_sm_ehs_B.Sum4_j =
          boost_and_two_level__1_sm_ehs_B.Gain1_b +
          boost_and_two_level__1_sm_ehs_B.Product2_b;

        /* Product: '<S68>/Product2' */
        boost_and_two_level__1_sm_ehs_B.Product2_m =
          boost_and_two_level__1_sm_ehs_B.Sum5_n /
          boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_e;

        /* Product: '<S68>/Product4' */
        boost_and_two_level__1_sm_ehs_B.Product4_h =
          boost_and_two_level__1_sm_ehs_B.Product2_m *
          boost_and_two_level__1_sm_ehs_B.Sum4_j;

        /* Sum: '<S67>/Sum7' */
        boost_and_two_level__1_sm_ehs_B.Sum7_g =
          boost_and_two_level__1_sm_ehs_B.Integ4_i -
          boost_and_two_level__1_sm_ehs_B.SFunction_j;

        /* Product: '<S67>/Product' */
        boost_and_two_level__1_sm_ehs_B.Meanvalue_b =
          boost_and_two_level__1_sm_ehs_B.Sum7_g *
          boost_and_two_level__1_sm_ehs_B.UnitDelay_l;

        /* Sum: '<S67>/Sum5' */
        boost_and_two_level__1_sm_ehs_B.Sum5_nt =
          boost_and_two_level__1_sm_ehs_B.Meanvalue_b +
          boost_and_two_level__1_sm_ehs_B.Product4_h;
        boost_and_two_level__1_sm_ehs_B.Switch_cf =
          boost_and_two_level__1_sm_ehs_B.Sum5_nt;
      } else {
        boost_and_two_level__1_sm_ehs_B.Switch_cf =
          boost_and_two_level__1_sm_ehs_B.UnitDelay1_n;
      }

      /* End of Switch: '<S67>/Switch' */

      /* RealImagToComplex: '<S61>/Real-Imag to Complex' */
      boost_and_two_level__1_sm_ehs_B.RealImagtoComplex_f.re =
        boost_and_two_level__1_sm_ehs_B.Switch_h;
      boost_and_two_level__1_sm_ehs_B.RealImagtoComplex_f.im =
        boost_and_two_level__1_sm_ehs_B.Switch_cf;

      /* ComplexToMagnitudeAngle: '<S61>/Complex to Magnitude-Angle' */
      boost_and_two_level__1_sm_ehs_B.ComplextoMagnitudeAngle_o1_c =
        rt_hypotd_snf(boost_and_two_level__1_sm_ehs_B.RealImagtoComplex_f.re,
                      boost_and_two_level__1_sm_ehs_B.RealImagtoComplex_f.im);
      boost_and_two_level__1_sm_ehs_B.ComplextoMagnitudeAngle_o2_k =
        rt_atan2d_snf(boost_and_two_level__1_sm_ehs_B.RealImagtoComplex_f.im,
                      boost_and_two_level__1_sm_ehs_B.RealImagtoComplex_f.re);

      /* Gain: '<S61>/Rad->Deg.' */
      boost_and_two_level__1_sm_ehs_B.RadDeg_l =
        boost_and_two_level__1_sm_ehs_P.RadDeg_Gain *
        boost_and_two_level__1_sm_ehs_B.ComplextoMagnitudeAngle_o2_k;

      /* Saturate: '<S57>/Saturation' */
      B = boost_and_two_level__1_sm_ehs_B.ComplextoMagnitudeAngle_o1_c;
      dP = boost_and_two_level__1_sm_ehs_P.Saturation_LowerSat_n;
      u2 = boost_and_two_level__1_sm_ehs_P.Saturation_UpperSat_b;
      if (B > u2) {
        boost_and_two_level__1_sm_ehs_B.Saturation_o = u2;
      } else if (B < dP) {
        boost_and_two_level__1_sm_ehs_B.Saturation_o = dP;
      } else {
        boost_and_two_level__1_sm_ehs_B.Saturation_o = B;
      }

      /* End of Saturate: '<S57>/Saturation' */

      /* Math: '<S57>/Math Function'
       *
       * About '<S57>/Math Function':
       *  Operator: reciprocal
       */
      P = boost_and_two_level__1_sm_ehs_B.Saturation_o;
      boost_and_two_level__1_sm_ehs_B.MathFunction_c = 1.0 / P;
    }

    /* End of Outputs for SubSystem: '<S56>/Automatic Gain Control' */

    /* Trigonometry: '<S56>/Trigonometric Function2' */
    boost_and_two_level__1_sm_ehs_B.TrigonometricFunction2 = cos
      (boost_and_two_level__1_sm_ehs_B.MathFunction_m);

    /* Product: '<S56>/Product1' */
    boost_and_two_level__1_sm_ehs_B.Product1_g =
      boost_and_two_level__1_sm_ehs_B.Vpu *
      boost_and_two_level__1_sm_ehs_B.TrigonometricFunction2;

    /* DiscreteIntegrator: '<S70>/Integ4' */
    if (boost_and_two_level__1_sm_ehs_DW.Integ4_SYSTEM_ENABLE_e != 0) {
      boost_and_two_level__1_sm_ehs_B.Integ4_e =
        boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_k;
    } else {
      boost_and_two_level__1_sm_ehs_B.Integ4_e =
        boost_and_two_level__1_sm_ehs_P.Integ4_gainval_b *
        boost_and_two_level__1_sm_ehs_B.Product1_g +
        boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_k;
    }

    /* End of DiscreteIntegrator: '<S70>/Integ4' */

    /* Saturate: '<S70>/To avoid division  by zero' */
    B = boost_and_two_level__1_sm_ehs_B.UnitDelay_l;
    dP = boost_and_two_level__1_sm_ehs_P.Toavoiddivisionbyzero_LowerSa_e;
    u2 = boost_and_two_level__1_sm_ehs_P.Toavoiddivisionbyzero_UpperS_ix;
    if (B > u2) {
      boost_and_two_level__1_sm_ehs_B.Freq_a = u2;
    } else if (B < dP) {
      boost_and_two_level__1_sm_ehs_B.Freq_a = dP;
    } else {
      boost_and_two_level__1_sm_ehs_B.Freq_a = B;
    }

    /* End of Saturate: '<S70>/To avoid division  by zero' */

    /* Fcn: '<S70>/Number of samples per cycle' */
    boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_g = 1.0 /
      boost_and_two_level__1_sm_ehs_B.Freq_a / 2.0e-5;

    /* Rounding: '<S70>/Rounding Function' */
    boost_and_two_level__1_sm_ehs_B.RoundingFunction_l = ceil
      (boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_g);

    /* Gain: '<S70>/Gain' */
    boost_and_two_level__1_sm_ehs_B.Delay_pz =
      boost_and_two_level__1_sm_ehs_P.Ts *
      boost_and_two_level__1_sm_ehs_B.RoundingFunction_l;

    /* Level2 S-Function Block: '<S72>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[29];
      sfcnOutputs(rts, 1);
    }

    /* UnitDelay: '<S71>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_B.UnitDelay_lg =
      boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_f;

    /* DigitalClock: '<S70>/Digital  Clock' */
    boost_and_two_level__1_sm_ehs_B.DigitalClock_a =
      boost_and_two_level__1_sm_ehs_M->Timing.t[1];

    /* RelationalOperator: '<S70>/Relational Operator' incorporates:
     *  Constant: '<S70>/Constant'
     */
    boost_and_two_level__1_sm_ehs_B.RelationalOperator_k =
      (boost_and_two_level__1_sm_ehs_B.DigitalClock_a >=
       boost_and_two_level__1_sm_ehs_P.Constant_Value_c);

    /* UnitDelay: '<S70>/Unit Delay1' */
    boost_and_two_level__1_sm_ehs_B.UnitDelay1_i =
      boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_e;

    /* Switch: '<S70>/Switch' */
    if (boost_and_two_level__1_sm_ehs_B.RelationalOperator_k) {
      /* Sum: '<S71>/Sum1' */
      boost_and_two_level__1_sm_ehs_B.Sum1 =
        boost_and_two_level__1_sm_ehs_B.Product1_g -
        boost_and_two_level__1_sm_ehs_B.UnitDelay_lg;

      /* Sum: '<S71>/Sum5' */
      boost_and_two_level__1_sm_ehs_B.Sum5 =
        boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_g -
        boost_and_two_level__1_sm_ehs_B.RoundingFunction_l;

      /* Product: '<S71>/Product5' */
      boost_and_two_level__1_sm_ehs_B.Product5 =
        boost_and_two_level__1_sm_ehs_B.Sum5 *
        boost_and_two_level__1_sm_ehs_B.Sum1;

      /* Gain: '<S71>/Gain1' */
      boost_and_two_level__1_sm_ehs_B.Gain1_e =
        boost_and_two_level__1_sm_ehs_P.Gain1_Gain_n *
        boost_and_two_level__1_sm_ehs_B.Product5;

      /* Sum: '<S71>/Sum4' */
      boost_and_two_level__1_sm_ehs_B.Sum4 =
        boost_and_two_level__1_sm_ehs_B.Gain1_e +
        boost_and_two_level__1_sm_ehs_B.Product1_g;

      /* Product: '<S71>/Product2' */
      boost_and_two_level__1_sm_ehs_B.Product2_f =
        boost_and_two_level__1_sm_ehs_B.Sum5 /
        boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_g;

      /* Product: '<S71>/Product4' */
      boost_and_two_level__1_sm_ehs_B.Product4 =
        boost_and_two_level__1_sm_ehs_B.Product2_f *
        boost_and_two_level__1_sm_ehs_B.Sum4;

      /* Sum: '<S70>/Sum7' */
      boost_and_two_level__1_sm_ehs_B.Sum7 =
        boost_and_two_level__1_sm_ehs_B.Integ4_e -
        boost_and_two_level__1_sm_ehs_B.SFunction_i;

      /* Product: '<S70>/Product' */
      boost_and_two_level__1_sm_ehs_B.Meanvalue =
        boost_and_two_level__1_sm_ehs_B.Sum7 *
        boost_and_two_level__1_sm_ehs_B.UnitDelay_l;

      /* Sum: '<S70>/Sum5' */
      boost_and_two_level__1_sm_ehs_B.Sum5_h =
        boost_and_two_level__1_sm_ehs_B.Meanvalue +
        boost_and_two_level__1_sm_ehs_B.Product4;
      boost_and_two_level__1_sm_ehs_B.Switch_cd =
        boost_and_two_level__1_sm_ehs_B.Sum5_h;
    } else {
      boost_and_two_level__1_sm_ehs_B.Switch_cd =
        boost_and_two_level__1_sm_ehs_B.UnitDelay1_i;
    }

    /* End of Switch: '<S70>/Switch' */

    /* Product: '<S56>/Divide' */
    boost_and_two_level__1_sm_ehs_B.Divide_f =
      boost_and_two_level__1_sm_ehs_B.Switch_cd *
      boost_and_two_level__1_sm_ehs_B.MathFunction_c;

    /* DiscreteTransferFcn: '<S58>/Discrete Derivative ' */
    P = boost_and_two_level__1_sm_ehs_B.Divide_f;
    P -= boost_and_two_level__1_sm_ehs_P.DiscreteDerivative_DenCoef[1] *
      boost_and_two_level__1_sm_ehs_DW.DiscreteDerivative_states;
    P /= boost_and_two_level__1_sm_ehs_P.DiscreteDerivative_DenCoef[0];
    boost_and_two_level__1_sm_ehs_DW.DiscreteDerivative_tmp = P;
    B = 1.0;
    P = boost_and_two_level__1_sm_ehs_P.Discrete_Kd;
    B *= P;
    dV = B * boost_and_two_level__1_sm_ehs_DW.DiscreteDerivative_tmp;
    B = (-1.0);
    P = boost_and_two_level__1_sm_ehs_P.Discrete_Kd;
    B *= P;
    dV += B * boost_and_two_level__1_sm_ehs_DW.DiscreteDerivative_states;
    boost_and_two_level__1_sm_ehs_B.DiscreteDerivative = dV;

    /* DiscreteIntegrator: '<S58>/Discrete-Time Integrator' */
    boost_and_two_level__1_sm_ehs_B.DiscreteTimeIntegrator_i =
      boost_and_two_level__1_sm_ehs_DW.DiscreteTimeIntegrator_DSTATE_o;

    /* Gain: '<S58>/Kp4' */
    boost_and_two_level__1_sm_ehs_B.Kp4 =
      boost_and_two_level__1_sm_ehs_P.Discrete_Kp *
      boost_and_two_level__1_sm_ehs_B.Divide_f;

    /* Sum: '<S58>/Sum6' */
    boost_and_two_level__1_sm_ehs_B.Sum6 = (boost_and_two_level__1_sm_ehs_B.Kp4
      + boost_and_two_level__1_sm_ehs_B.DiscreteTimeIntegrator_i) +
      boost_and_two_level__1_sm_ehs_B.DiscreteDerivative;

    /* Saturate: '<S58>/Saturation1' */
    B = boost_and_two_level__1_sm_ehs_B.Sum6;
    dP = boost_and_two_level__1_sm_ehs_P.Saturation1_LowerSat_j;
    u2 = boost_and_two_level__1_sm_ehs_P.Saturation1_UpperSat_c;
    if (B > u2) {
      boost_and_two_level__1_sm_ehs_B.Saturation1 = u2;
    } else if (B < dP) {
      boost_and_two_level__1_sm_ehs_B.Saturation1 = dP;
    } else {
      boost_and_two_level__1_sm_ehs_B.Saturation1 = B;
    }

    /* End of Saturate: '<S58>/Saturation1' */

    /* Gain: '<S56>/Gain10' */
    boost_and_two_level__1_sm_ehs_B.Gain10 =
      boost_and_two_level__1_sm_ehs_P.Gain10_Gain *
      boost_and_two_level__1_sm_ehs_B.Saturation1;

    /* RateLimiter: '<S56>/Rate Limiter' */
    P = boost_and_two_level__1_sm_ehs_B.Gain10 -
      boost_and_two_level__1_sm_ehs_DW.PrevY;
    if (P > boost_and_two_level__1_sm_ehs_P.RateLimiter_RisingLim) {
      boost_and_two_level__1_sm_ehs_B.RateLimiter =
        boost_and_two_level__1_sm_ehs_DW.PrevY +
        boost_and_two_level__1_sm_ehs_P.RateLimiter_RisingLim;
    } else if (P < boost_and_two_level__1_sm_ehs_P.RateLimiter_FallingLim) {
      boost_and_two_level__1_sm_ehs_B.RateLimiter =
        boost_and_two_level__1_sm_ehs_DW.PrevY +
        boost_and_two_level__1_sm_ehs_P.RateLimiter_FallingLim;
    } else {
      boost_and_two_level__1_sm_ehs_B.RateLimiter =
        boost_and_two_level__1_sm_ehs_B.Gain10;
    }

    boost_and_two_level__1_sm_ehs_DW.PrevY =
      boost_and_two_level__1_sm_ehs_B.RateLimiter;

    /* End of RateLimiter: '<S56>/Rate Limiter' */

    /* UnitDelay: '<S73>/Delay_x1' */
    boost_and_two_level__1_sm_ehs_B.x1k_e =
      boost_and_two_level__1_sm_ehs_DW.Delay_x1_DSTATE_p;

    /* Gain: '<S74>/A11' */
    boost_and_two_level__1_sm_ehs_B.A11 =
      boost_and_two_level__1_sm_ehs_P.A11_Gain *
      boost_and_two_level__1_sm_ehs_B.x1k_e;

    /* UnitDelay: '<S73>/Delay_x2' */
    boost_and_two_level__1_sm_ehs_B.x2k_h =
      boost_and_two_level__1_sm_ehs_DW.Delay_x2_DSTATE_l;

    /* Gain: '<S74>/A12' */
    boost_and_two_level__1_sm_ehs_B.A12 =
      boost_and_two_level__1_sm_ehs_P.A12_Gain *
      boost_and_two_level__1_sm_ehs_B.x2k_h;

    /* Gain: '<S74>/A21' */
    boost_and_two_level__1_sm_ehs_B.A21 =
      boost_and_two_level__1_sm_ehs_P.A21_Gain *
      boost_and_two_level__1_sm_ehs_B.x1k_e;

    /* Gain: '<S74>/A22' */
    boost_and_two_level__1_sm_ehs_B.A22 =
      boost_and_two_level__1_sm_ehs_P.A22_Gain *
      boost_and_two_level__1_sm_ehs_B.x2k_h;

    /* Sum: '<S74>/sum2' */
    boost_and_two_level__1_sm_ehs_B.sum2_c = boost_and_two_level__1_sm_ehs_B.A11
      + boost_and_two_level__1_sm_ehs_B.A12;

    /* Sum: '<S74>/sum3' */
    boost_and_two_level__1_sm_ehs_B.sum3 = boost_and_two_level__1_sm_ehs_B.A21 +
      boost_and_two_level__1_sm_ehs_B.A22;

    /* Gain: '<S75>/B11' */
    boost_and_two_level__1_sm_ehs_B.B11 =
      boost_and_two_level__1_sm_ehs_P.B11_Gain *
      boost_and_two_level__1_sm_ehs_B.RateLimiter;

    /* Sum: '<S73>/A*x1(k) + B*u1(k) ' */
    boost_and_two_level__1_sm_ehs_B.x1k1 =
      boost_and_two_level__1_sm_ehs_B.sum2_c +
      boost_and_two_level__1_sm_ehs_B.B11;

    /* Gain: '<S75>/B21' */
    boost_and_two_level__1_sm_ehs_B.B21 =
      boost_and_two_level__1_sm_ehs_P.B21_Gain *
      boost_and_two_level__1_sm_ehs_B.RateLimiter;

    /* Sum: '<S73>/A*x2(k) + B*u2(k)' */
    boost_and_two_level__1_sm_ehs_B.x2k1 = boost_and_two_level__1_sm_ehs_B.sum3
      + boost_and_two_level__1_sm_ehs_B.B21;

    /* Gain: '<S73>/D*u(k)' */
    boost_and_two_level__1_sm_ehs_B.Duk_c =
      boost_and_two_level__1_sm_ehs_P.Duk_Gain_b *
      boost_and_two_level__1_sm_ehs_B.RateLimiter;

    /* Gain: '<S76>/C11' */
    boost_and_two_level__1_sm_ehs_B.C11_f =
      boost_and_two_level__1_sm_ehs_P.C11_Gain_b *
      boost_and_two_level__1_sm_ehs_B.x1k_e;

    /* Gain: '<S76>/C12' */
    boost_and_two_level__1_sm_ehs_B.C12_o =
      boost_and_two_level__1_sm_ehs_P.C12_Gain_f *
      boost_and_two_level__1_sm_ehs_B.x2k_h;

    /* Sum: '<S76>/sum2' */
    boost_and_two_level__1_sm_ehs_B.sum2_n =
      boost_and_two_level__1_sm_ehs_B.C11_f +
      boost_and_two_level__1_sm_ehs_B.C12_o;

    /* Sum: '<S73>/C*X(k)+D*u(k)' */
    boost_and_two_level__1_sm_ehs_B.yk_c = boost_and_two_level__1_sm_ehs_B.Duk_c
      + boost_and_two_level__1_sm_ehs_B.sum2_n;

    /* Gain: '<S78>/A11' */
    boost_and_two_level__1_sm_ehs_B.A11_n[0] =
      boost_and_two_level__1_sm_ehs_P.A11_Gain_k *
      boost_and_two_level__1_sm_ehs_B.x1k[0];
    boost_and_two_level__1_sm_ehs_B.A11_n[1] =
      boost_and_two_level__1_sm_ehs_P.A11_Gain_k *
      boost_and_two_level__1_sm_ehs_B.x1k[1];

    /* Gain: '<S78>/A12' */
    boost_and_two_level__1_sm_ehs_B.A12_e[0] =
      boost_and_two_level__1_sm_ehs_P.A12_Gain_k *
      boost_and_two_level__1_sm_ehs_B.x2k[0];
    boost_and_two_level__1_sm_ehs_B.A12_e[1] =
      boost_and_two_level__1_sm_ehs_P.A12_Gain_k *
      boost_and_two_level__1_sm_ehs_B.x2k[1];

    /* Gain: '<S78>/A21' */
    boost_and_two_level__1_sm_ehs_B.A21_b[0] =
      boost_and_two_level__1_sm_ehs_P.A21_Gain_d *
      boost_and_two_level__1_sm_ehs_B.x1k[0];
    boost_and_two_level__1_sm_ehs_B.A21_b[1] =
      boost_and_two_level__1_sm_ehs_P.A21_Gain_d *
      boost_and_two_level__1_sm_ehs_B.x1k[1];

    /* Gain: '<S78>/A22' */
    boost_and_two_level__1_sm_ehs_B.A22_o[0] =
      boost_and_two_level__1_sm_ehs_P.A22_Gain_k *
      boost_and_two_level__1_sm_ehs_B.x2k[0];
    boost_and_two_level__1_sm_ehs_B.A22_o[1] =
      boost_and_two_level__1_sm_ehs_P.A22_Gain_k *
      boost_and_two_level__1_sm_ehs_B.x2k[1];

    /* Sum: '<S78>/sum2' */
    boost_and_two_level__1_sm_ehs_B.sum2_a[0] =
      boost_and_two_level__1_sm_ehs_B.A11_n[0] +
      boost_and_two_level__1_sm_ehs_B.A12_e[0];
    boost_and_two_level__1_sm_ehs_B.sum2_a[1] =
      boost_and_two_level__1_sm_ehs_B.A11_n[1] +
      boost_and_two_level__1_sm_ehs_B.A12_e[1];

    /* Sum: '<S78>/sum3' */
    boost_and_two_level__1_sm_ehs_B.sum3_f[0] =
      boost_and_two_level__1_sm_ehs_B.A21_b[0] +
      boost_and_two_level__1_sm_ehs_B.A22_o[0];
    boost_and_two_level__1_sm_ehs_B.sum3_f[1] =
      boost_and_two_level__1_sm_ehs_B.A21_b[1] +
      boost_and_two_level__1_sm_ehs_B.A22_o[1];

    /* Gain: '<S79>/B11' */
    boost_and_two_level__1_sm_ehs_B.B11_a[0] =
      boost_and_two_level__1_sm_ehs_P.B11_Gain_f *
      boost_and_two_level__1_sm_ehs_B.Switch_b[0];
    boost_and_two_level__1_sm_ehs_B.B11_a[1] =
      boost_and_two_level__1_sm_ehs_P.B11_Gain_f *
      boost_and_two_level__1_sm_ehs_B.Switch_b[1];

    /* Sum: '<S77>/A*x1(k) + B*u1(k) ' */
    boost_and_two_level__1_sm_ehs_B.x1k1_l[0] =
      boost_and_two_level__1_sm_ehs_B.sum2_a[0] +
      boost_and_two_level__1_sm_ehs_B.B11_a[0];
    boost_and_two_level__1_sm_ehs_B.x1k1_l[1] =
      boost_and_two_level__1_sm_ehs_B.sum2_a[1] +
      boost_and_two_level__1_sm_ehs_B.B11_a[1];

    /* Gain: '<S79>/B21' */
    boost_and_two_level__1_sm_ehs_B.B21_a[0] =
      boost_and_two_level__1_sm_ehs_P.B21_Gain_j *
      boost_and_two_level__1_sm_ehs_B.Switch_b[0];
    boost_and_two_level__1_sm_ehs_B.B21_a[1] =
      boost_and_two_level__1_sm_ehs_P.B21_Gain_j *
      boost_and_two_level__1_sm_ehs_B.Switch_b[1];

    /* Sum: '<S77>/A*x2(k) + B*u2(k)' */
    boost_and_two_level__1_sm_ehs_B.x2k1_j[0] =
      boost_and_two_level__1_sm_ehs_B.sum3_f[0] +
      boost_and_two_level__1_sm_ehs_B.B21_a[0];
    boost_and_two_level__1_sm_ehs_B.x2k1_j[1] =
      boost_and_two_level__1_sm_ehs_B.sum3_f[1] +
      boost_and_two_level__1_sm_ehs_B.B21_a[1];

    /* Constant: '<S45>/Constant1' */
    boost_and_two_level__1_sm_ehs_B.Constant1 =
      boost_and_two_level__1_sm_ehs_P.Constant1_Value_a;

    /* Switch: '<S31>/Switch' incorporates:
     *  Constant: '<S31>/Vnom_dc1'
     *  Constant: '<S9>/MPPT_On'
     */
    if (boost_and_two_level__1_sm_ehs_P.MPPT_On_Value != 0.0) {
      boost_and_two_level__1_sm_ehs_B.Switch_ex =
        boost_and_two_level__1_sm_ehs_B.D;
    } else {
      boost_and_two_level__1_sm_ehs_B.Switch_ex =
        boost_and_two_level__1_sm_ehs_P.InverterControl_Vdc_ref_Init;
    }

    /* End of Switch: '<S31>/Switch' */

    /* Sum: '<S38>/Add1' incorporates:
     *  Constant: '<S38>/Constant2'
     *  Constant: '<S38>/Constant4'
     */
    dV = boost_and_two_level__1_sm_ehs_P.Ts *
      boost_and_two_level__1_sm_ehs_P.InverterControl_Fnom * 6.2831853071795862;
    boost_and_two_level__1_sm_ehs_B.Add1_j =
      (boost_and_two_level__1_sm_ehs_B.MathFunction_m +
       boost_and_two_level__1_sm_ehs_P.Constant2_Value_i) + dV;

    /* UnitDelay: '<S31>/Unit Delay3' */
    boost_and_two_level__1_sm_ehs_B.UnitDelay3[0] =
      boost_and_two_level__1_sm_ehs_DW.UnitDelay3_DSTATE[0];
    boost_and_two_level__1_sm_ehs_B.UnitDelay3[1] =
      boost_and_two_level__1_sm_ehs_DW.UnitDelay3_DSTATE[1];

    /* Gain: '<S38>/Gain1' */
    boost_and_two_level__1_sm_ehs_B.Gain1_p =
      boost_and_two_level__1_sm_ehs_P.Gain1_Gain_ne *
      boost_and_two_level__1_sm_ehs_B.Switch_ej;

    /* Product: '<S38>/Product' incorporates:
     *  Constant: '<S38>/Constant3'
     */
    dV = boost_and_two_level__1_sm_ehs_P.InverterControl_Vnom_prim *
      1.4142135623730951;
    boost_and_two_level__1_sm_ehs_B.Product =
      boost_and_two_level__1_sm_ehs_B.Gain1_p / dV;

    /* Product: '<S38>/Product1' */
    boost_and_two_level__1_sm_ehs_B.Product1_m[0] =
      boost_and_two_level__1_sm_ehs_B.UnitDelay3[0] /
      boost_and_two_level__1_sm_ehs_B.Product;
    boost_and_two_level__1_sm_ehs_B.Product1_m[1] =
      boost_and_two_level__1_sm_ehs_B.UnitDelay3[1] /
      boost_and_two_level__1_sm_ehs_B.Product;

    /* RealImagToComplex: '<S38>/Real-Imag to Complex' */
    boost_and_two_level__1_sm_ehs_B.RealImagtoComplex_j.re =
      boost_and_two_level__1_sm_ehs_B.Product1_m[0];
    boost_and_two_level__1_sm_ehs_B.RealImagtoComplex_j.im =
      boost_and_two_level__1_sm_ehs_B.Product1_m[1];

    /* ComplexToMagnitudeAngle: '<S38>/Complex to Magnitude-Angle' */
    boost_and_two_level__1_sm_ehs_B.ComplextoMagnitudeAngle_o1_k = rt_hypotd_snf
      (boost_and_two_level__1_sm_ehs_B.RealImagtoComplex_j.re,
       boost_and_two_level__1_sm_ehs_B.RealImagtoComplex_j.im);
    boost_and_two_level__1_sm_ehs_B.ComplextoMagnitudeAngle_o2_g = rt_atan2d_snf
      (boost_and_two_level__1_sm_ehs_B.RealImagtoComplex_j.im,
       boost_and_two_level__1_sm_ehs_B.RealImagtoComplex_j.re);

    /* Sum: '<S38>/Add2' */
    boost_and_two_level__1_sm_ehs_B.Add2_l =
      boost_and_two_level__1_sm_ehs_B.ComplextoMagnitudeAngle_o2_g +
      boost_and_two_level__1_sm_ehs_B.Add1_j;

    /* Trigonometry: '<S38>/Trigonometric Function' */
    boost_and_two_level__1_sm_ehs_B.TrigonometricFunction_c = sin
      (boost_and_two_level__1_sm_ehs_B.Add2_l);

    /* Product: '<S38>/Product2' */
    boost_and_two_level__1_sm_ehs_B.Product2_g =
      boost_and_two_level__1_sm_ehs_B.ComplextoMagnitudeAngle_o1_k *
      boost_and_two_level__1_sm_ehs_B.TrigonometricFunction_c;

    /* UnitDelay: '<S31>/Unit Delay1' */
    boost_and_two_level__1_sm_ehs_B.UnitDelay1_k =
      boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_n;

    /* Sum: '<S39>/Sum' */
    boost_and_two_level__1_sm_ehs_B.Sum_p =
      boost_and_two_level__1_sm_ehs_B.Switch_ej -
      boost_and_two_level__1_sm_ehs_B.UnitDelay1_k;

    /* Gain: '<S39>/Rtot_pu2' */
    B = boost_and_two_level__1_sm_ehs_P.InverterControl_Vnom_dc;
    dV = 1.0 / B;
    boost_and_two_level__1_sm_ehs_B.Rtot_pu2 = dV *
      boost_and_two_level__1_sm_ehs_B.Sum_p;

    /* Gain: '<S111>/Integral Gain' */
    boost_and_two_level__1_sm_ehs_B.IntegralGain_m =
      boost_and_two_level__1_sm_ehs_P.InverterControl_Ki_VDCreg *
      boost_and_two_level__1_sm_ehs_B.Rtot_pu2;

    /* DiscreteIntegrator: '<S111>/Integrator' */
    boost_and_two_level__1_sm_ehs_B.Integrator_j =
      boost_and_two_level__1_sm_ehs_P.Integrator_gainval_e *
      boost_and_two_level__1_sm_ehs_B.IntegralGain_m +
      boost_and_two_level__1_sm_ehs_DW.Integrator_DSTATE_h;

    /* Gain: '<S111>/Proportional Gain' */
    boost_and_two_level__1_sm_ehs_B.ProportionalGain_e =
      boost_and_two_level__1_sm_ehs_P.InverterControl_Kp_VDCreg *
      boost_and_two_level__1_sm_ehs_B.Rtot_pu2;

    /* Sum: '<S111>/Sum' */
    boost_and_two_level__1_sm_ehs_B.Sum_e =
      boost_and_two_level__1_sm_ehs_B.ProportionalGain_e +
      boost_and_two_level__1_sm_ehs_B.Integrator_j;

    /* Saturate: '<S111>/Saturate' */
    B = boost_and_two_level__1_sm_ehs_B.Sum_e;
    dP = boost_and_two_level__1_sm_ehs_P.PI_LowerSaturationLimit_j;
    u2 = boost_and_two_level__1_sm_ehs_P.PI_UpperSaturationLimit_a;
    if (B > u2) {
      boost_and_two_level__1_sm_ehs_B.Saturate_d = u2;
    } else if (B < dP) {
      boost_and_two_level__1_sm_ehs_B.Saturate_d = dP;
    } else {
      boost_and_two_level__1_sm_ehs_B.Saturate_d = B;
    }

    /* End of Saturate: '<S111>/Saturate' */

    /* RateTransition: '<S9>/Rate Transition3' incorporates:
     *  RateTransition: '<S9>/Rate Transition1'
     *  RateTransition: '<S9>/Rate Transition2'
     */
    if (rtmIsMajorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
      boost_and_two_level__1_sm_ehs_B.RateTransition3 =
        boost_and_two_level__1_sm_ehs_B.Divide[2];
      boost_and_two_level__1_sm_ehs_B.RateTransition1 =
        boost_and_two_level__1_sm_ehs_B.SFunction_l[0];
      boost_and_two_level__1_sm_ehs_B.RateTransition2 =
        boost_and_two_level__1_sm_ehs_B.SFunction_l[1];
    }

    /* End of RateTransition: '<S9>/Rate Transition3' */

    /* MATLAB Function: '<S9>/MATLAB Function' incorporates:
     *  Constant: '<S9>/Constant1'
     *  Constant: '<S9>/Constant2'
     *  Constant: '<S9>/Constant3'
     *  Constant: '<S9>/Constant4'
     */
    /* MATLAB Function 'SM_eHS/source_and_gates/MATLAB Function': '<S32>:1' */
    /* '<S32>:1:3' */
    /* '<S32>:1:4' */
    /* '<S32>:1:6' */
    /* '<S32>:1:7' */
    /* '<S32>:1:8' */
    /*  bandgap energy(1.12 eV) */
    /*  ideal factor */
    /* '<S32>:1:13' */
    /* '<S32>:1:14' */
    /* '<S32>:1:17' */
    boost_and_two_level__1_sm_ehs_B.Iph =
      (((boost_and_two_level__1_sm_ehs_B.RateTransition2 + 273.0) - 298.0) *
       (boost_and_two_level__1_sm_ehs_P.Constant3_Value_b / 100.0) + 1.0) *
      (boost_and_two_level__1_sm_ehs_B.RateTransition1 / 1000.0 *
       boost_and_two_level__1_sm_ehs_P.Constant1_Value_p);

    /* '<S32>:1:18' */
    /* '<S32>:1:19' */
    /* '<S32>:1:21' */
    boost_and_two_level__1_sm_ehs_B.Io =
      boost_and_two_level__1_sm_ehs_P.Constant1_Value_p / (exp
      (boost_and_two_level__1_sm_ehs_P.Constant2_Value_c /
       (boost_and_two_level__1_sm_ehs_P.Constant4_Value * 1.381E-23 * 1.3 *
        298.0 / 1.6E-19)) - 1.0) * rt_powd_snf
      ((boost_and_two_level__1_sm_ehs_B.RateTransition2 + 273.0) / 298.0, 3.0) *
      exp((0.0033557046979865771 - 1.0 /
           (boost_and_two_level__1_sm_ehs_B.RateTransition2 + 273.0)) *
          9981.61867097421) * (exp
      (boost_and_two_level__1_sm_ehs_B.RateTransition3 /
       (boost_and_two_level__1_sm_ehs_P.Constant4_Value * 1.381E-23 * 1.3 *
        (boost_and_two_level__1_sm_ehs_B.RateTransition2 + 273.0) / 1.6E-19)) -
      1.0);

    /* Saturate: '<S9>/Saturation' */
    B = boost_and_two_level__1_sm_ehs_B.Iph;
    dP = boost_and_two_level__1_sm_ehs_P.Saturation_LowerSat_c;
    u2 = boost_and_two_level__1_sm_ehs_P.Saturation_UpperSat_b1;
    if (B > u2) {
      boost_and_two_level__1_sm_ehs_B.Saturation_p = u2;
    } else if (B < dP) {
      boost_and_two_level__1_sm_ehs_B.Saturation_p = dP;
    } else {
      boost_and_two_level__1_sm_ehs_B.Saturation_p = B;
    }

    /* End of Saturate: '<S9>/Saturation' */
  }

  /* Clock: '<S2>/Clock' */
  boost_and_two_level__1_sm_ehs_B.Clock_o =
    boost_and_two_level__1_sm_ehs_M->Timing.t[0];

  /* RelationalOperator: '<S3>/Compare' incorporates:
   *  Constant: '<S3>/Constant'
   */
  boost_and_two_level__1_sm_ehs_B.Compare_p =
    (boost_and_two_level__1_sm_ehs_B.Clock_o >=
     boost_and_two_level__1_sm_ehs_P.CompareToConstant_const_d);
  if (rtmIsMajorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
    /* Constant: '<S2>/Constant1' */
    boost_and_two_level__1_sm_ehs_B.Constant1_f =
      boost_and_two_level__1_sm_ehs_P.Constant1_Value_h;
  }

  /* DataTypeConversion: '<S2>/Data Type Conversion' */
  boost_and_two_level__1_sm_ehs_B.DataTypeConversion_c =
    boost_and_two_level__1_sm_ehs_B.Compare_p;
  if (rtmIsMajorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
    /* Level2 S-Function Block: '<S2>/OpTrigger' (optrigger) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[30];
      sfcnOutputs(rts, 1);
    }

    /* Level2 S-Function Block: '<S2>/OpCtrl' (sfun_ctrl_op7160ex1) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[31];
      sfcnOutputs(rts, 1);
    }
  }
}

/* Model update function */
void boost_and_two_level__1_sm_ehs_update(void)
{
  int32_T i;
  if (rtmIsMajorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
    /* Update for Memory: '<S1>/S-Function' */
    boost_and_two_level__1_sm_ehs_DW.SFunction_PreviousInput =
      boost_and_two_level__1_sm_ehs_B.Sum;

    /* Update for UnitDelay: '<S31>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE =
      boost_and_two_level__1_sm_ehs_B.Product2_g;

    /* Update for UnitDelay: '<S25>/Delay Input1' */
    boost_and_two_level__1_sm_ehs_DW.DelayInput1_DSTATE =
      boost_and_two_level__1_sm_ehs_B.Add_l;

    /* Update for Memory: '<S9>/Memory1' */
    boost_and_two_level__1_sm_ehs_DW.Memory1_PreviousInput =
      boost_and_two_level__1_sm_ehs_B.Saturation_p;

    /* Update for Memory: '<S9>/Memory' */
    boost_and_two_level__1_sm_ehs_DW.Memory_PreviousInput =
      boost_and_two_level__1_sm_ehs_B.Io;
  }

  /* Update for TransportDelay: '<S15>/Transport Delay' */
  {
    real_T **uBuffer = (real_T**)
      &boost_and_two_level__1_sm_ehs_DW.TransportDelay_PWORK.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)
      &boost_and_two_level__1_sm_ehs_DW.TransportDelay_PWORK.TUbufferPtrs[1];
    real_T simTime = boost_and_two_level__1_sm_ehs_M->Timing.t[0];
    boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.Head =
      ((boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.Head <
        (boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.CircularBufSize-1))
       ? (boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.Head+1) : 0);
    if (boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.Head ==
        boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.Tail) {
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.Tail =
        ((boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.Tail <
          (boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.CircularBufSize
           -1)) ? (boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.Tail+1)
         : 0);
    }

    (*tBuffer)[boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.Head] =
      simTime;
    (*uBuffer)[boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.Head] =
      boost_and_two_level__1_sm_ehs_B.integrator;
  }

  if (rtmIsMajorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
    /* Update for Memory: '<S15>/Memory' */
    boost_and_two_level__1_sm_ehs_DW.Memory_PreviousInput_n =
      boost_and_two_level__1_sm_ehs_B.Switch;
  }

  /* Update for TransportDelay: '<S16>/Transport Delay' */
  {
    real_T **uBuffer = (real_T**)
      &boost_and_two_level__1_sm_ehs_DW.TransportDelay_PWORK_d.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)
      &boost_and_two_level__1_sm_ehs_DW.TransportDelay_PWORK_d.TUbufferPtrs[1];
    real_T simTime = boost_and_two_level__1_sm_ehs_M->Timing.t[0];
    boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.Head =
      ((boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.Head <
        (boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.CircularBufSize
         -1)) ? (boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.Head+1)
       : 0);
    if (boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.Head ==
        boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.Tail) {
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.Tail =
        ((boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.Tail <
          (boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.CircularBufSize
           -1)) ? (boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.Tail+
                   1) : 0);
    }

    (*tBuffer)[boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.Head] =
      simTime;
    (*uBuffer)[boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.Head] =
      boost_and_two_level__1_sm_ehs_B.integrator_p;
  }

  if (rtmIsMajorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
    /* Update for Memory: '<S16>/Memory' */
    boost_and_two_level__1_sm_ehs_DW.Memory_PreviousInput_e =
      boost_and_two_level__1_sm_ehs_B.Switch_e;

    /* Update for Memory: '<S19>/Memory1' */
    boost_and_two_level__1_sm_ehs_DW.Memory1_PreviousInput_h =
      boost_and_two_level__1_sm_ehs_B.Switch1;

    /* Update for UnitDelay: '<S24>/Delay Input1' */
    memcpy(&boost_and_two_level__1_sm_ehs_DW.DelayInput1_DSTATE_m[0],
           &boost_and_two_level__1_sm_ehs_B.DataTypeConversion5[0], sizeof
           (uint32_T) << 4U);
    memcpy(&boost_and_two_level__1_sm_ehs_DW.DelayInput1_DSTATE_m[16],
           &boost_and_two_level__1_sm_ehs_B.ConvertdoubletoSinglefloating_o[0],
           sizeof(uint32_T) << 4U);
    memcpy(&boost_and_two_level__1_sm_ehs_DW.DelayInput1_DSTATE_m[32],
           &boost_and_two_level__1_sm_ehs_B.ConvertdoubletoSinglefloatin_ou[0],
           sizeof(uint32_T) << 4U);
    memcpy(&boost_and_two_level__1_sm_ehs_DW.DelayInput1_DSTATE_m[48],
           &boost_and_two_level__1_sm_ehs_B.DataTypeConversion5_m[0], sizeof
           (uint32_T) << 4U);
    memcpy(&boost_and_two_level__1_sm_ehs_DW.DelayInput1_DSTATE_m[64],
           &boost_and_two_level__1_sm_ehs_B.DataTypeConversion5_n[0], sizeof
           (uint32_T) << 4U);

    /* Update for Memory: '<S19>/Memory2' */
    boost_and_two_level__1_sm_ehs_DW.Memory2_PreviousInput =
      boost_and_two_level__1_sm_ehs_B.Switch_i;

    /* Update for UnitDelay: '<S28>/Delay Input1' */
    for (i = 0; i < 10; i++) {
      boost_and_two_level__1_sm_ehs_DW.DelayInput1_DSTATE_k[i] =
        boost_and_two_level__1_sm_ehs_B.load_config1[i];
    }

    /* End of Update for UnitDelay: '<S28>/Delay Input1' */

    /* Update for UnitDelay: '<S56>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_b =
      boost_and_two_level__1_sm_ehs_B.yk_c;

    /* Level2 S-Function Block: '<S82>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[25];
      sfcnUpdate(rts, 1);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Update for DiscreteIntegrator: '<S56>/Discrete-Time Integrator' */
    boost_and_two_level__1_sm_ehs_DW.DiscreteTimeIntegrator_DSTATE +=
      boost_and_two_level__1_sm_ehs_P.DiscreteTimeIntegrator_gainval *
      boost_and_two_level__1_sm_ehs_B.Saturation1;

    /* Update for UnitDelay: '<S77>/Delay_x1' */
    boost_and_two_level__1_sm_ehs_DW.Delay_x1_DSTATE[0] =
      boost_and_two_level__1_sm_ehs_B.x1k1_l[0];
    boost_and_two_level__1_sm_ehs_DW.Delay_x1_DSTATE[1] =
      boost_and_two_level__1_sm_ehs_B.x1k1_l[1];

    /* Update for UnitDelay: '<S77>/Delay_x2' */
    boost_and_two_level__1_sm_ehs_DW.Delay_x2_DSTATE[0] =
      boost_and_two_level__1_sm_ehs_B.x2k1_j[0];
    boost_and_two_level__1_sm_ehs_DW.Delay_x2_DSTATE[1] =
      boost_and_two_level__1_sm_ehs_B.x2k1_j[1];

    /* Update for UnitDelay: '<S31>/Unit Delay2' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay2_DSTATE =
      boost_and_two_level__1_sm_ehs_B.Saturate_d;

    /* Update for DiscreteIntegrator: '<S40>/Integrator' */
    boost_and_two_level__1_sm_ehs_DW.Integrator_DSTATE[0] +=
      boost_and_two_level__1_sm_ehs_P.Integrator_gainval *
      boost_and_two_level__1_sm_ehs_B.IntegralGain[0];
    boost_and_two_level__1_sm_ehs_DW.Integrator_DSTATE[1] +=
      boost_and_two_level__1_sm_ehs_P.Integrator_gainval *
      boost_and_two_level__1_sm_ehs_B.IntegralGain[1];

    /* Update for DiscreteIntegrator: '<S48>/Integ4' */
    boost_and_two_level__1_sm_ehs_DW.Integ4_SYSTEM_ENABLE = 0U;
    boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE =
      boost_and_two_level__1_sm_ehs_P.Integ4_gainval_l *
      boost_and_two_level__1_sm_ehs_B.Product1_k +
      boost_and_two_level__1_sm_ehs_B.Integ4;

    /* Level2 S-Function Block: '<S50>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[26];
      sfcnUpdate(rts, 1);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Update for UnitDelay: '<S49>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_d =
      boost_and_two_level__1_sm_ehs_B.Product1_k;

    /* Update for UnitDelay: '<S48>/Unit Delay1' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE =
      boost_and_two_level__1_sm_ehs_B.Switch_ie;

    /* Update for DiscreteIntegrator: '<S51>/Integ4' */
    boost_and_two_level__1_sm_ehs_DW.Integ4_SYSTEM_ENABLE_n = 0U;
    boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_f =
      boost_and_two_level__1_sm_ehs_P.Integ4_gainval_g *
      boost_and_two_level__1_sm_ehs_B.Product2 +
      boost_and_two_level__1_sm_ehs_B.Integ4_b;

    /* Level2 S-Function Block: '<S53>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[27];
      sfcnUpdate(rts, 1);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Update for UnitDelay: '<S52>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_i =
      boost_and_two_level__1_sm_ehs_B.Product2;

    /* Update for UnitDelay: '<S51>/Unit Delay1' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_i =
      boost_and_two_level__1_sm_ehs_B.Switch_c;

    /* Update for DiscreteIntegrator: '<S54>/Integ4' */
    boost_and_two_level__1_sm_ehs_DW.Integ4_SYSTEM_ENABLE_i = 0U;
    boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_g =
      boost_and_two_level__1_sm_ehs_P.Integ4_gainval_o *
      boost_and_two_level__1_sm_ehs_B.RateTransition6 +
      boost_and_two_level__1_sm_ehs_B.Integ4_m;

    /* Level2 S-Function Block: '<S55>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[28];
      sfcnUpdate(rts, 1);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Update for UnitDelay: '<S54>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_h =
      boost_and_two_level__1_sm_ehs_B.RateTransition6;

    /* Update for UnitDelay: '<S54>/Unit Delay1' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_im =
      boost_and_two_level__1_sm_ehs_B.Switch_ej;

    /* Update for Enabled SubSystem: '<S56>/Automatic Gain Control' incorporates:
     *  Update for EnablePort: '<S57>/Enable'
     */
    if (boost_and_two_level__1_sm_ehs_DW.AutomaticGainControl_MODE) {
      /* Update for DiscreteIntegrator: '<S64>/Integ4' */
      boost_and_two_level__1_sm_ehs_DW.Integ4_SYSTEM_ENABLE_nz = 0U;
      boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_i =
        boost_and_two_level__1_sm_ehs_P.Integ4_gainval *
        boost_and_two_level__1_sm_ehs_B.Product1_a +
        boost_and_two_level__1_sm_ehs_B.Integ4_ee;

      /* Level2 S-Function Block: '<S66>/S-Function' (sfun_discreteVariableDelay) */
      {
        SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[0];
        sfcnUpdate(rts, 1);
        if (ssGetErrorStatus(rts) != (NULL))
          return;
      }

      /* Update for UnitDelay: '<S65>/Unit Delay' */
      boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_c =
        boost_and_two_level__1_sm_ehs_B.Product1_a;

      /* Update for UnitDelay: '<S64>/Unit Delay1' */
      boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_eb =
        boost_and_two_level__1_sm_ehs_B.Switch_h;

      /* Update for DiscreteIntegrator: '<S67>/Integ4' */
      boost_and_two_level__1_sm_ehs_DW.Integ4_SYSTEM_ENABLE_d = 0U;
      boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_c =
        boost_and_two_level__1_sm_ehs_P.Integ4_gainval_c *
        boost_and_two_level__1_sm_ehs_B.Product2_b +
        boost_and_two_level__1_sm_ehs_B.Integ4_i;

      /* Level2 S-Function Block: '<S69>/S-Function' (sfun_discreteVariableDelay) */
      {
        SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[1];
        sfcnUpdate(rts, 1);
        if (ssGetErrorStatus(rts) != (NULL))
          return;
      }

      /* Update for UnitDelay: '<S68>/Unit Delay' */
      boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_fv =
        boost_and_two_level__1_sm_ehs_B.Product2_b;

      /* Update for UnitDelay: '<S67>/Unit Delay1' */
      boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_k =
        boost_and_two_level__1_sm_ehs_B.Switch_cf;
    }

    /* End of Update for SubSystem: '<S56>/Automatic Gain Control' */

    /* Update for DiscreteIntegrator: '<S70>/Integ4' */
    boost_and_two_level__1_sm_ehs_DW.Integ4_SYSTEM_ENABLE_e = 0U;
    boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_k =
      boost_and_two_level__1_sm_ehs_P.Integ4_gainval_b *
      boost_and_two_level__1_sm_ehs_B.Product1_g +
      boost_and_two_level__1_sm_ehs_B.Integ4_e;

    /* Level2 S-Function Block: '<S72>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[29];
      sfcnUpdate(rts, 1);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Update for UnitDelay: '<S71>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_f =
      boost_and_two_level__1_sm_ehs_B.Product1_g;

    /* Update for UnitDelay: '<S70>/Unit Delay1' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_e =
      boost_and_two_level__1_sm_ehs_B.Switch_cd;

    /* Update for DiscreteTransferFcn: '<S58>/Discrete Derivative ' */
    boost_and_two_level__1_sm_ehs_DW.DiscreteDerivative_states =
      boost_and_two_level__1_sm_ehs_DW.DiscreteDerivative_tmp;

    /* Update for DiscreteIntegrator: '<S58>/Discrete-Time Integrator' */
    boost_and_two_level__1_sm_ehs_DW.DiscreteTimeIntegrator_DSTATE_o +=
      boost_and_two_level__1_sm_ehs_P.DiscreteTimeIntegrator_gainva_a *
      boost_and_two_level__1_sm_ehs_B.Divide_f;
    if (boost_and_two_level__1_sm_ehs_DW.DiscreteTimeIntegrator_DSTATE_o >=
        boost_and_two_level__1_sm_ehs_P.DiscreteTimeIntegrator_UpperSat) {
      boost_and_two_level__1_sm_ehs_DW.DiscreteTimeIntegrator_DSTATE_o =
        boost_and_two_level__1_sm_ehs_P.DiscreteTimeIntegrator_UpperSat;
    } else {
      if (boost_and_two_level__1_sm_ehs_DW.DiscreteTimeIntegrator_DSTATE_o <=
          boost_and_two_level__1_sm_ehs_P.DiscreteTimeIntegrator_LowerSat) {
        boost_and_two_level__1_sm_ehs_DW.DiscreteTimeIntegrator_DSTATE_o =
          boost_and_two_level__1_sm_ehs_P.DiscreteTimeIntegrator_LowerSat;
      }
    }

    /* End of Update for DiscreteIntegrator: '<S58>/Discrete-Time Integrator' */

    /* Update for UnitDelay: '<S73>/Delay_x1' */
    boost_and_two_level__1_sm_ehs_DW.Delay_x1_DSTATE_p =
      boost_and_two_level__1_sm_ehs_B.x1k1;

    /* Update for UnitDelay: '<S73>/Delay_x2' */
    boost_and_two_level__1_sm_ehs_DW.Delay_x2_DSTATE_l =
      boost_and_two_level__1_sm_ehs_B.x2k1;

    /* Update for UnitDelay: '<S31>/Unit Delay3' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay3_DSTATE[0] =
      boost_and_two_level__1_sm_ehs_B.Saturation[0];
    boost_and_two_level__1_sm_ehs_DW.UnitDelay3_DSTATE[1] =
      boost_and_two_level__1_sm_ehs_B.Saturation[1];

    /* Update for UnitDelay: '<S31>/Unit Delay1' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_n =
      boost_and_two_level__1_sm_ehs_B.Switch_ex;

    /* Update for DiscreteIntegrator: '<S111>/Integrator' */
    boost_and_two_level__1_sm_ehs_DW.Integrator_DSTATE_h =
      boost_and_two_level__1_sm_ehs_P.Integrator_gainval_e *
      boost_and_two_level__1_sm_ehs_B.IntegralGain_m +
      boost_and_two_level__1_sm_ehs_B.Integrator_j;
  }

  if (rtmIsMajorTimeStep(boost_and_two_level__1_sm_ehs_M)) {
    rt_ertODEUpdateContinuousStates(&boost_and_two_level__1_sm_ehs_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++boost_and_two_level__1_sm_ehs_M->Timing.clockTick0)) {
    ++boost_and_two_level__1_sm_ehs_M->Timing.clockTickH0;
  }

  boost_and_two_level__1_sm_ehs_M->Timing.t[0] = rtsiGetSolverStopTime
    (&boost_and_two_level__1_sm_ehs_M->solverInfo);

  {
    /* Update absolute timer for sample time: [2.0E-5s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++boost_and_two_level__1_sm_ehs_M->Timing.clockTick1)) {
      ++boost_and_two_level__1_sm_ehs_M->Timing.clockTickH1;
    }

    boost_and_two_level__1_sm_ehs_M->Timing.t[1] =
      boost_and_two_level__1_sm_ehs_M->Timing.clockTick1 *
      boost_and_two_level__1_sm_ehs_M->Timing.stepSize1 +
      boost_and_two_level__1_sm_ehs_M->Timing.clockTickH1 *
      boost_and_two_level__1_sm_ehs_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void boost_and_two_level__1_sm_ehs_derivatives(void)
{
  XDot_boost_and_two_level__1_sm_ehs_T *_rtXdot;
  _rtXdot = ((XDot_boost_and_two_level__1_sm_ehs_T *)
             boost_and_two_level__1_sm_ehs_M->ModelData.derivs);

  /* Derivatives for Integrator: '<S15>/integrator' */
  _rtXdot->integrator_CSTATE = boost_and_two_level__1_sm_ehs_B.Divide[0];

  /* Derivatives for Integrator: '<S16>/integrator' */
  _rtXdot->integrator_CSTATE_b = boost_and_two_level__1_sm_ehs_B.Divide[1];
}

/* Model initialize function */
void boost_and_two_level__1_sm_ehs_initialize(void)
{
  {
    int32_T i;

    /* Start for TransportDelay: '<S15>/Transport Delay' */
    {
      real_T *pBuffer =
        &boost_and_two_level__1_sm_ehs_DW.TransportDelay_RWORK.TUbufferArea[0];
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.Tail = 0;
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.Head = 0;
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.Last = 0;
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK.CircularBufSize =
        8192;
      pBuffer[0] = boost_and_two_level__1_sm_ehs_P.TransportDelay_InitOutput;
      pBuffer[8192] = boost_and_two_level__1_sm_ehs_M->Timing.t[0];
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_PWORK.TUbufferPtrs[0] =
        (void *) &pBuffer[0];
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_PWORK.TUbufferPtrs[1] =
        (void *) &pBuffer[8192];
    }

    /* Start for TransportDelay: '<S16>/Transport Delay' */
    {
      real_T *pBuffer =
        &boost_and_two_level__1_sm_ehs_DW.TransportDelay_RWORK_o.TUbufferArea[0];
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.Tail = 0;
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.Head = 0;
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.Last = 0;
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_IWORK_g.CircularBufSize =
        8192;
      pBuffer[0] = boost_and_two_level__1_sm_ehs_P.TransportDelay_InitOutput_f;
      pBuffer[8192] = boost_and_two_level__1_sm_ehs_M->Timing.t[0];
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_PWORK_d.TUbufferPtrs[0] =
        (void *) &pBuffer[0];
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_PWORK_d.TUbufferPtrs[1] =
        (void *) &pBuffer[8192];
    }

    /* Start for Constant: '<S7>/load_config1' */
    for (i = 0; i < 10; i++) {
      boost_and_two_level__1_sm_ehs_B.load_config1[i] =
        boost_and_two_level__1_sm_ehs_P.load_config1_Value[i];
    }

    /* End of Start for Constant: '<S7>/load_config1' */
    /* Level2 S-Function Block: '<S9>/RTE Conversion1' (rte_conversion) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[19];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Level2 S-Function Block: '<S9>/RTE Conversion' (rte_conversion) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[20];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Level2 S-Function Block: '<S9>/RTE Logical Operator1' (rte_logical_operator) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[21];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Level2 S-Function Block: '<S9>/RTE Ground' (rte_ground) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[22];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Level2 S-Function Block: '<S26>/RTE_Conversion_1' (rte_conversion_ophsdio) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[23];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Level2 S-Function Block: '<S82>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[25];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Start for Enabled SubSystem: '<S81>/Subsystem1' */
    /* VirtualOutportStart for Outport: '<S86>/dq' */
    boost_and_two_level__1_sm_ehs_B.Fcn =
      boost_and_two_level__1_sm_ehs_P.dq_Y0_b[0];
    boost_and_two_level__1_sm_ehs_B.Fcn1 =
      boost_and_two_level__1_sm_ehs_P.dq_Y0_b[1];

    /* End of Start for SubSystem: '<S81>/Subsystem1' */

    /* Start for Enabled SubSystem: '<S81>/Subsystem - pi//2 delay' */
    /* VirtualOutportStart for Outport: '<S85>/dq' */
    boost_and_two_level__1_sm_ehs_B.Fcn_i =
      boost_and_two_level__1_sm_ehs_P.dq_Y0[0];
    boost_and_two_level__1_sm_ehs_B.Fcn1_k =
      boost_and_two_level__1_sm_ehs_P.dq_Y0[1];

    /* End of Start for SubSystem: '<S81>/Subsystem - pi//2 delay' */
    /* Level2 S-Function Block: '<S50>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[26];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Level2 S-Function Block: '<S53>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[27];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Start for Constant: '<S54>/K1' */
    boost_and_two_level__1_sm_ehs_B.K1 =
      boost_and_two_level__1_sm_ehs_P.K1_Value_f;

    /* Level2 S-Function Block: '<S55>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[28];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Start for Enabled SubSystem: '<S56>/Automatic Gain Control' */

    /* Level2 S-Function Block: '<S66>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[0];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Level2 S-Function Block: '<S69>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[1];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* End of Start for SubSystem: '<S56>/Automatic Gain Control' */

    /* InitializeConditions for Enabled SubSystem: '<S56>/Automatic Gain Control' */
    /* InitializeConditions for DiscreteIntegrator: '<S64>/Integ4' */
    boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_i =
      boost_and_two_level__1_sm_ehs_P.Integ4_IC;

    /* Level2 S-Function Block: '<S66>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[0];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S65>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_c =
      boost_and_two_level__1_sm_ehs_P.UnitDelay_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S64>/Unit Delay1' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_eb =
      boost_and_two_level__1_sm_ehs_P.UnitDelay1_InitialCondition;

    /* InitializeConditions for DiscreteIntegrator: '<S67>/Integ4' */
    boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_c =
      boost_and_two_level__1_sm_ehs_P.Integ4_IC_b;

    /* Level2 S-Function Block: '<S69>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[1];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S68>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_fv =
      boost_and_two_level__1_sm_ehs_P.UnitDelay_InitialCondition_k;

    /* InitializeConditions for UnitDelay: '<S67>/Unit Delay1' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_k =
      boost_and_two_level__1_sm_ehs_P.UnitDelay1_InitialCondition_n;

    /* End of InitializeConditions for SubSystem: '<S56>/Automatic Gain Control' */

    /* Start for Enabled SubSystem: '<S56>/Automatic Gain Control' */
    /* VirtualOutportStart for Outport: '<S57>/Gain' */
    boost_and_two_level__1_sm_ehs_B.MathFunction_c =
      boost_and_two_level__1_sm_ehs_P.Gain_Y0;

    /* End of Start for SubSystem: '<S56>/Automatic Gain Control' */
    /* Level2 S-Function Block: '<S72>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[29];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Start for Constant: '<S2>/Constant1' */
    boost_and_two_level__1_sm_ehs_B.Constant1_f =
      boost_and_two_level__1_sm_ehs_P.Constant1_Value_h;
  }

  {
    int32_T i;

    /* user code (Initialize function Body) */

    /* System '<Root>' */

    /* Opal-RT Technologies */
    opalSizeDwork = sizeof(rtDWork);

#ifdef USE_RTMODEL

    if (Opal_rtmGetNumBlockIO(pRtModel) != 0)
      opalSizeBlockIO = sizeof(rtB);
    else
      opalSizeBlockIO = 0;
    if (Opal_rtmGetNumBlockParams(pRtModel) != 0)
      opalSizeRTP = sizeof(rtP);
    else
      opalSizeRTP = 0;

#else

    if (ssGetNumBlockIO(rtS) != 0)
      opalSizeBlockIO = sizeof(rtB);
    else
      opalSizeBlockIO = 0;
    if (ssGetNumBlockParams(rtS) != 0)
      opalSizeRTP = sizeof(rtP);
    else
      opalSizeRTP = 0;

#endif

    /* InitializeConditions for Memory: '<S1>/S-Function' */
    boost_and_two_level__1_sm_ehs_DW.SFunction_PreviousInput =
      boost_and_two_level__1_sm_ehs_P.SFunction_X0;

    /* Level2 S-Function Block: '<S7>/Outputs_eHS1_Recv' (sfun_fct_op7160ex1_recv) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[2];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S31>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE =
      boost_and_two_level__1_sm_ehs_P.UnitDelay_InitialCondition_i;

    /* Level2 S-Function Block: '<S10>/S-Function' (RECV_Param) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[4];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S25>/Delay Input1' */
    boost_and_two_level__1_sm_ehs_DW.DelayInput1_DSTATE =
      boost_and_two_level__1_sm_ehs_P.DetectChange_vinit;

    /* Level2 S-Function Block: '<S7>/eHS_rst_loadin' (sfun_fct_op7160ex1_load_in) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[5];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Level2 S-Function Block: '<S7>/Automated_Solver_Mat_Initialisation_1' (sfun_efs_solver_cfg) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[6];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for Memory: '<S9>/Memory1' */
    boost_and_two_level__1_sm_ehs_DW.Memory1_PreviousInput =
      boost_and_two_level__1_sm_ehs_P.Memory1_X0;

    /* InitializeConditions for Memory: '<S9>/Memory' */
    boost_and_two_level__1_sm_ehs_DW.Memory_PreviousInput =
      boost_and_two_level__1_sm_ehs_P.Memory_X0;

    /* Level2 S-Function Block: '<S7>/Inputs_eHS1_Send' (sfun_fct_op7160ex1_send) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[8];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for Integrator: '<S15>/integrator' */
    boost_and_two_level__1_sm_ehs_X.integrator_CSTATE =
      boost_and_two_level__1_sm_ehs_P.integrator_IC;

    /* InitializeConditions for Memory: '<S15>/Memory' */
    boost_and_two_level__1_sm_ehs_DW.Memory_PreviousInput_n =
      boost_and_two_level__1_sm_ehs_P.Memory_X0_a;

    /* InitializeConditions for Integrator: '<S16>/integrator' */
    boost_and_two_level__1_sm_ehs_X.integrator_CSTATE_b =
      boost_and_two_level__1_sm_ehs_P.integrator_IC_g;

    /* InitializeConditions for Memory: '<S16>/Memory' */
    boost_and_two_level__1_sm_ehs_DW.Memory_PreviousInput_e =
      boost_and_two_level__1_sm_ehs_P.Memory_X0_c;

    /* Level2 S-Function Block: '<S2>/Outputs_eHS1_Recv' (sfun_fct_op7160ex1_recv) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[9];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Level2 S-Function Block: '<S2>/OpMonitor' (opmonitor) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[11];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Level2 S-Function Block: '<S30>/S-Function' (OP_SEND) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[12];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for Memory: '<S19>/Memory1' */
    boost_and_two_level__1_sm_ehs_DW.Memory1_PreviousInput_h =
      boost_and_two_level__1_sm_ehs_P.Memory1_X0_c;

    /* InitializeConditions for UnitDelay: '<S24>/Delay Input1' */
    for (i = 0; i < 80; i++) {
      boost_and_two_level__1_sm_ehs_DW.DelayInput1_DSTATE_m[i] =
        boost_and_two_level__1_sm_ehs_P.DetectChange_vinit_i;
    }

    /* End of InitializeConditions for UnitDelay: '<S24>/Delay Input1' */

    /* InitializeConditions for Memory: '<S19>/Memory2' */
    boost_and_two_level__1_sm_ehs_DW.Memory2_PreviousInput =
      boost_and_two_level__1_sm_ehs_P.Memory2_X0;

    /* Level2 S-Function Block: '<S17>/LoadIn' (sfun_fct_op7160ex1_load_in) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[15];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Level2 S-Function Block: '<S17>/DataIn Send' (sfun_fct_op7160ex1_send) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[17];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Level2 S-Function Block: '<S2>/OpWriteFile' (opwritefile) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[18];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S28>/Delay Input1' */
    for (i = 0; i < 10; i++) {
      boost_and_two_level__1_sm_ehs_DW.DelayInput1_DSTATE_k[i] =
        boost_and_two_level__1_sm_ehs_P.DetectChange_vinit_g;
    }

    /* End of InitializeConditions for UnitDelay: '<S28>/Delay Input1' */
    /* Level2 S-Function Block: '<S9>/RTE Conversion1' (rte_conversion) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[19];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Level2 S-Function Block: '<S9>/RTE Conversion' (rte_conversion) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[20];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Level2 S-Function Block: '<S9>/RTE Logical Operator1' (rte_logical_operator) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[21];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Level2 S-Function Block: '<S9>/RTE Ground' (rte_ground) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[22];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Level2 S-Function Block: '<S26>/RTE_Conversion_1' (rte_conversion_ophsdio) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[23];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Level2 S-Function Block: '<S26>/EventGen_eHS_1' (sfun_op7160ex1_event_generator) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[24];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S56>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_b =
      boost_and_two_level__1_sm_ehs_P.UnitDelay_InitialCondition_a;

    /* Level2 S-Function Block: '<S82>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[25];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for DiscreteIntegrator: '<S56>/Discrete-Time Integrator' */
    boost_and_two_level__1_sm_ehs_DW.DiscreteTimeIntegrator_DSTATE =
      boost_and_two_level__1_sm_ehs_P.DiscreteTimeIntegrator_IC;

    /* InitializeConditions for UnitDelay: '<S77>/Delay_x1' */
    boost_and_two_level__1_sm_ehs_DW.Delay_x1_DSTATE[0] =
      boost_and_two_level__1_sm_ehs_P.Delay_x1_InitialCondition[0];
    boost_and_two_level__1_sm_ehs_DW.Delay_x1_DSTATE[1] =
      boost_and_two_level__1_sm_ehs_P.Delay_x1_InitialCondition[1];

    /* InitializeConditions for UnitDelay: '<S77>/Delay_x2' */
    boost_and_two_level__1_sm_ehs_DW.Delay_x2_DSTATE[0] =
      boost_and_two_level__1_sm_ehs_P.Delay_x2_InitialCondition[0];
    boost_and_two_level__1_sm_ehs_DW.Delay_x2_DSTATE[1] =
      boost_and_two_level__1_sm_ehs_P.Delay_x2_InitialCondition[1];

    /* InitializeConditions for UnitDelay: '<S31>/Unit Delay2' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay2_DSTATE =
      boost_and_two_level__1_sm_ehs_P.UnitDelay2_InitialCondition;

    /* InitializeConditions for DiscreteIntegrator: '<S40>/Integrator' */
    boost_and_two_level__1_sm_ehs_DW.Integrator_DSTATE[0] =
      boost_and_two_level__1_sm_ehs_P.Integrator_IC;
    boost_and_two_level__1_sm_ehs_DW.Integrator_DSTATE[1] =
      boost_and_two_level__1_sm_ehs_P.Integrator_IC;

    /* InitializeConditions for DiscreteIntegrator: '<S48>/Integ4' */
    boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE =
      boost_and_two_level__1_sm_ehs_P.Integ4_IC_e;

    /* Level2 S-Function Block: '<S50>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[26];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S49>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_d =
      boost_and_two_level__1_sm_ehs_P.UnitDelay_InitialCondition_f;

    /* InitializeConditions for UnitDelay: '<S48>/Unit Delay1' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE =
      boost_and_two_level__1_sm_ehs_P.UnitDelay1_InitialCondition_j;

    /* InitializeConditions for DiscreteIntegrator: '<S51>/Integ4' */
    boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_f =
      boost_and_two_level__1_sm_ehs_P.Integ4_IC_l;

    /* Level2 S-Function Block: '<S53>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[27];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S52>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_i =
      boost_and_two_level__1_sm_ehs_P.UnitDelay_InitialCondition_d;

    /* InitializeConditions for UnitDelay: '<S51>/Unit Delay1' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_i =
      boost_and_two_level__1_sm_ehs_P.UnitDelay1_InitialCondition_e;

    /* InitializeConditions for MATLAB Function: '<S31>/MPPT Controller using Perturbe  & Observe technique  ' */
    boost_and_two_level__1_sm_ehs_DW.Vold_not_empty = false;
    boost_and_two_level__1_sm_ehs_DW.Vold = 0.0;
    boost_and_two_level__1_sm_ehs_DW.Pold = 0.0;

    /* InitializeConditions for DiscreteIntegrator: '<S54>/Integ4' */
    boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_g =
      boost_and_two_level__1_sm_ehs_P.Integ4_IC_o;

    /* Level2 S-Function Block: '<S55>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[28];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S54>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_h =
      boost_and_two_level__1_sm_ehs_P.UnitDelay_InitialCondition_m;

    /* InitializeConditions for UnitDelay: '<S54>/Unit Delay1' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_im =
      boost_and_two_level__1_sm_ehs_P.UnitDelay1_InitialCondition_m;

    /* InitializeConditions for DiscreteIntegrator: '<S70>/Integ4' */
    boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_k =
      boost_and_two_level__1_sm_ehs_P.Integ4_IC_m;

    /* Level2 S-Function Block: '<S72>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[29];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S71>/Unit Delay' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_f =
      boost_and_two_level__1_sm_ehs_P.UnitDelay_InitialCondition_e;

    /* InitializeConditions for UnitDelay: '<S70>/Unit Delay1' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_e =
      boost_and_two_level__1_sm_ehs_P.UnitDelay1_InitialCondition_m1;

    /* InitializeConditions for DiscreteTransferFcn: '<S58>/Discrete Derivative ' */
    boost_and_two_level__1_sm_ehs_DW.DiscreteDerivative_states =
      boost_and_two_level__1_sm_ehs_P.DiscreteDerivative_InitialState;

    /* InitializeConditions for DiscreteIntegrator: '<S58>/Discrete-Time Integrator' */
    boost_and_two_level__1_sm_ehs_DW.DiscreteTimeIntegrator_DSTATE_o =
      boost_and_two_level__1_sm_ehs_P.Discrete_Init;

    /* InitializeConditions for RateLimiter: '<S56>/Rate Limiter' */
    boost_and_two_level__1_sm_ehs_DW.PrevY =
      boost_and_two_level__1_sm_ehs_P.RateLimiter_IC;

    /* InitializeConditions for UnitDelay: '<S73>/Delay_x1' */
    boost_and_two_level__1_sm_ehs_DW.Delay_x1_DSTATE_p =
      boost_and_two_level__1_sm_ehs_P.Delay_x1_InitialCondition_b;

    /* InitializeConditions for UnitDelay: '<S73>/Delay_x2' */
    boost_and_two_level__1_sm_ehs_DW.Delay_x2_DSTATE_l =
      boost_and_two_level__1_sm_ehs_P.Delay_x2_InitialCondition_h;

    /* InitializeConditions for UnitDelay: '<S31>/Unit Delay3' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay3_DSTATE[0] =
      boost_and_two_level__1_sm_ehs_P.UnitDelay3_InitialCondition;
    boost_and_two_level__1_sm_ehs_DW.UnitDelay3_DSTATE[1] =
      boost_and_two_level__1_sm_ehs_P.UnitDelay3_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S31>/Unit Delay1' */
    boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_n =
      boost_and_two_level__1_sm_ehs_P.UnitDelay1_InitialCondition_eb;

    /* InitializeConditions for DiscreteIntegrator: '<S111>/Integrator' */
    boost_and_two_level__1_sm_ehs_DW.Integrator_DSTATE_h =
      boost_and_two_level__1_sm_ehs_P.Integrator_IC_p;

    /* Level2 S-Function Block: '<S2>/OpTrigger' (optrigger) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[30];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Level2 S-Function Block: '<S2>/OpCtrl' (sfun_ctrl_op7160ex1) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[31];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }
  }

  /* Enable for DiscreteIntegrator: '<S48>/Integ4' */
  boost_and_two_level__1_sm_ehs_DW.Integ4_SYSTEM_ENABLE = 1U;

  /* Enable for DiscreteIntegrator: '<S51>/Integ4' */
  boost_and_two_level__1_sm_ehs_DW.Integ4_SYSTEM_ENABLE_n = 1U;

  /* Enable for DiscreteIntegrator: '<S54>/Integ4' */
  boost_and_two_level__1_sm_ehs_DW.Integ4_SYSTEM_ENABLE_i = 1U;

  /* Enable for DiscreteIntegrator: '<S70>/Integ4' */
  boost_and_two_level__1_sm_ehs_DW.Integ4_SYSTEM_ENABLE_e = 1U;
}

/* Model terminate function */
void boost_and_two_level__1_sm_ehs_terminate(void)
{
  /* Level2 S-Function Block: '<S7>/Outputs_eHS1_Recv' (sfun_fct_op7160ex1_recv) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[2];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S7>/Convert  Single floating-point (FPGA)  to double' (sfun_SFP2DBL) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[3];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S10>/S-Function' (RECV_Param) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[4];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S7>/eHS_rst_loadin' (sfun_fct_op7160ex1_load_in) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[5];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S7>/Automated_Solver_Mat_Initialisation_1' (sfun_efs_solver_cfg) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[6];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S7>/Convert double to  Single floating-point (FPGA)' (sfun_DBL2SFP) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[7];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S7>/Inputs_eHS1_Send' (sfun_fct_op7160ex1_send) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[8];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S2>/Outputs_eHS1_Recv' (sfun_fct_op7160ex1_recv) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[9];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S2>/sfp2dbl' (sfun_SFP2DBL) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[10];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S2>/OpMonitor' (opmonitor) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[11];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S30>/S-Function' (OP_SEND) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[12];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S17>/Convert double to  Single floating-point (FPGA)' (sfun_DBL2SFP) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[13];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S17>/Convert double to  Single floating-point (FPGA)1' (sfun_DBL2SFP) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[14];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S17>/LoadIn' (sfun_fct_op7160ex1_load_in) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[15];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S17>/Convert double to  Single floating-point (FPGA)2' (sfun_DBL2SFP) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[16];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S17>/DataIn Send' (sfun_fct_op7160ex1_send) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[17];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S2>/OpWriteFile' (opwritefile) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[18];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S9>/RTE Conversion1' (rte_conversion) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[19];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S9>/RTE Conversion' (rte_conversion) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[20];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S9>/RTE Logical Operator1' (rte_logical_operator) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[21];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S9>/RTE Ground' (rte_ground) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[22];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S26>/RTE_Conversion_1' (rte_conversion_ophsdio) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[23];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S26>/EventGen_eHS_1' (sfun_op7160ex1_event_generator) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[24];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S82>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[25];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S50>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[26];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S53>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[27];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S55>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[28];
    sfcnTerminate(rts);
  }

  /* Terminate for Enabled SubSystem: '<S56>/Automatic Gain Control' */

  /* Level2 S-Function Block: '<S66>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[0];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S69>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[1];
    sfcnTerminate(rts);
  }

  /* End of Terminate for SubSystem: '<S56>/Automatic Gain Control' */

  /* Level2 S-Function Block: '<S72>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[29];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S2>/OpTrigger' (optrigger) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[30];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S2>/OpCtrl' (sfun_ctrl_op7160ex1) */
  {
    SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[31];
    sfcnTerminate(rts);
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  boost_and_two_level__1_sm_ehs_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  boost_and_two_level__1_sm_ehs_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  boost_and_two_level__1_sm_ehs_initialize();
}

void MdlTerminate(void)
{
  boost_and_two_level__1_sm_ehs_terminate();
}

/* Registration function */
RT_MODEL_boost_and_two_level__1_sm_ehs_T *boost_and_two_level__1_sm_ehs(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  boost_and_two_level__1_sm_ehs_P.Saturation_UpperSat_b = rtInf;
  boost_and_two_level__1_sm_ehs_P.DiscreteTimeIntegrator_UpperSat = rtInf;
  boost_and_two_level__1_sm_ehs_P.Saturation1_UpperSat_c = rtInf;

  /* initialize real-time model */
  (void) memset((void *)boost_and_two_level__1_sm_ehs_M, 0,
                sizeof(RT_MODEL_boost_and_two_level__1_sm_ehs_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&boost_and_two_level__1_sm_ehs_M->solverInfo,
                          &boost_and_two_level__1_sm_ehs_M->Timing.simTimeStep);
    rtsiSetTPtr(&boost_and_two_level__1_sm_ehs_M->solverInfo, &rtmGetTPtr
                (boost_and_two_level__1_sm_ehs_M));
    rtsiSetStepSizePtr(&boost_and_two_level__1_sm_ehs_M->solverInfo,
                       &boost_and_two_level__1_sm_ehs_M->Timing.stepSize0);
    rtsiSetdXPtr(&boost_and_two_level__1_sm_ehs_M->solverInfo,
                 &boost_and_two_level__1_sm_ehs_M->ModelData.derivs);
    rtsiSetContStatesPtr(&boost_and_two_level__1_sm_ehs_M->solverInfo, (real_T **)
                         &boost_and_two_level__1_sm_ehs_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&boost_and_two_level__1_sm_ehs_M->solverInfo,
      &boost_and_two_level__1_sm_ehs_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&boost_and_two_level__1_sm_ehs_M->solverInfo,
                          (&rtmGetErrorStatus(boost_and_two_level__1_sm_ehs_M)));
    rtsiSetRTModelPtr(&boost_and_two_level__1_sm_ehs_M->solverInfo,
                      boost_and_two_level__1_sm_ehs_M);
  }

  rtsiSetSimTimeStep(&boost_and_two_level__1_sm_ehs_M->solverInfo,
                     MAJOR_TIME_STEP);
  boost_and_two_level__1_sm_ehs_M->ModelData.intgData.y =
    boost_and_two_level__1_sm_ehs_M->ModelData.odeY;
  boost_and_two_level__1_sm_ehs_M->ModelData.intgData.f[0] =
    boost_and_two_level__1_sm_ehs_M->ModelData.odeF[0];
  boost_and_two_level__1_sm_ehs_M->ModelData.intgData.f[1] =
    boost_and_two_level__1_sm_ehs_M->ModelData.odeF[1];
  boost_and_two_level__1_sm_ehs_M->ModelData.intgData.f[2] =
    boost_and_two_level__1_sm_ehs_M->ModelData.odeF[2];
  boost_and_two_level__1_sm_ehs_M->ModelData.contStates = ((real_T *)
    &boost_and_two_level__1_sm_ehs_X);
  rtsiSetSolverData(&boost_and_two_level__1_sm_ehs_M->solverInfo, (void *)
                    &boost_and_two_level__1_sm_ehs_M->ModelData.intgData);
  rtsiSetSolverName(&boost_and_two_level__1_sm_ehs_M->solverInfo,"ode3");
  boost_and_two_level__1_sm_ehs_M->solverInfoPtr =
    (&boost_and_two_level__1_sm_ehs_M->solverInfo);

  /* Initialize timing info */
  {
    int_T *mdlTsMap =
      boost_and_two_level__1_sm_ehs_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    boost_and_two_level__1_sm_ehs_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    boost_and_two_level__1_sm_ehs_M->Timing.sampleTimes =
      (&boost_and_two_level__1_sm_ehs_M->Timing.sampleTimesArray[0]);
    boost_and_two_level__1_sm_ehs_M->Timing.offsetTimes =
      (&boost_and_two_level__1_sm_ehs_M->Timing.offsetTimesArray[0]);

    /* task periods */
    boost_and_two_level__1_sm_ehs_M->Timing.sampleTimes[0] = (0.0);
    boost_and_two_level__1_sm_ehs_M->Timing.sampleTimes[1] = (2.0E-5);

    /* task offsets */
    boost_and_two_level__1_sm_ehs_M->Timing.offsetTimes[0] = (0.0);
    boost_and_two_level__1_sm_ehs_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(boost_and_two_level__1_sm_ehs_M,
             &boost_and_two_level__1_sm_ehs_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits =
      boost_and_two_level__1_sm_ehs_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    boost_and_two_level__1_sm_ehs_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(boost_and_two_level__1_sm_ehs_M, -1);
  boost_and_two_level__1_sm_ehs_M->Timing.stepSize0 = 2.0E-5;
  boost_and_two_level__1_sm_ehs_M->Timing.stepSize1 = 2.0E-5;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    boost_and_two_level__1_sm_ehs_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(boost_and_two_level__1_sm_ehs_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(boost_and_two_level__1_sm_ehs_M->rtwLogInfo, (NULL));
    rtliSetLogT(boost_and_two_level__1_sm_ehs_M->rtwLogInfo, "");
    rtliSetLogX(boost_and_two_level__1_sm_ehs_M->rtwLogInfo, "");
    rtliSetLogXFinal(boost_and_two_level__1_sm_ehs_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(boost_and_two_level__1_sm_ehs_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(boost_and_two_level__1_sm_ehs_M->rtwLogInfo, 0);
    rtliSetLogMaxRows(boost_and_two_level__1_sm_ehs_M->rtwLogInfo, 1000);
    rtliSetLogDecimation(boost_and_two_level__1_sm_ehs_M->rtwLogInfo, 1);
    rtliSetLogY(boost_and_two_level__1_sm_ehs_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(boost_and_two_level__1_sm_ehs_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(boost_and_two_level__1_sm_ehs_M->rtwLogInfo, (NULL));
  }

  boost_and_two_level__1_sm_ehs_M->solverInfoPtr =
    (&boost_and_two_level__1_sm_ehs_M->solverInfo);
  boost_and_two_level__1_sm_ehs_M->Timing.stepSize = (2.0E-5);
  rtsiSetFixedStepSize(&boost_and_two_level__1_sm_ehs_M->solverInfo, 2.0E-5);
  rtsiSetSolverMode(&boost_and_two_level__1_sm_ehs_M->solverInfo,
                    SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  boost_and_two_level__1_sm_ehs_M->ModelData.blockIO = ((void *)
    &boost_and_two_level__1_sm_ehs_B);
  (void) memset(((void *) &boost_and_two_level__1_sm_ehs_B), 0,
                sizeof(B_boost_and_two_level__1_sm_ehs_T));

  {
    int_T i;
    for (i = 0; i < 7; i++) {
      boost_and_two_level__1_sm_ehs_B.ConvertSinglefloatingpointFPGAt[i] = 0.0;
    }

    for (i = 0; i < 7; i++) {
      boost_and_two_level__1_sm_ehs_B.Divide[i] = 0.0;
    }

    for (i = 0; i < 164; i++) {
      boost_and_two_level__1_sm_ehs_B.SFunction_l[i] = 0.0;
    }

    for (i = 0; i < 7; i++) {
      boost_and_two_level__1_sm_ehs_B.fpga_raw[i] = 0.0;
    }

    for (i = 0; i < 9; i++) {
      boost_and_two_level__1_sm_ehs_B.Constant[i] = 0.0;
    }

    for (i = 0; i < 16; i++) {
      boost_and_two_level__1_sm_ehs_B.Saturation_f[i] = 0.0;
    }

    for (i = 0; i < 16; i++) {
      boost_and_two_level__1_sm_ehs_B.Gain1_l[i] = 0.0;
    }

    for (i = 0; i < 16; i++) {
      boost_and_two_level__1_sm_ehs_B.Saturation1_p[i] = 0.0;
    }

    for (i = 0; i < 16; i++) {
      boost_and_two_level__1_sm_ehs_B.Gain2[i] = 0.0;
    }

    for (i = 0; i < 16; i++) {
      boost_and_two_level__1_sm_ehs_B.Saturation_g[i] = 0.0;
    }

    for (i = 0; i < 16; i++) {
      boost_and_two_level__1_sm_ehs_B.Gain1_j[i] = 0.0;
    }

    for (i = 0; i < 16; i++) {
      boost_and_two_level__1_sm_ehs_B.Saturation1_c[i] = 0.0;
    }

    for (i = 0; i < 16; i++) {
      boost_and_two_level__1_sm_ehs_B.Gain2_d[i] = 0.0;
    }

    for (i = 0; i < 16; i++) {
      boost_and_two_level__1_sm_ehs_B.Saturation_h[i] = 0.0;
    }

    for (i = 0; i < 16; i++) {
      boost_and_two_level__1_sm_ehs_B.Gain1_g[i] = 0.0;
    }

    for (i = 0; i < 16; i++) {
      boost_and_two_level__1_sm_ehs_B.Saturation1_cu[i] = 0.0;
    }

    for (i = 0; i < 16; i++) {
      boost_and_two_level__1_sm_ehs_B.Gain2_a[i] = 0.0;
    }

    boost_and_two_level__1_sm_ehs_B.RealImagtoComplex.re = 0.0;
    boost_and_two_level__1_sm_ehs_B.RealImagtoComplex.im = 0.0;
    boost_and_two_level__1_sm_ehs_B.MagnitudeAngletoComplex.re = 0.0;
    boost_and_two_level__1_sm_ehs_B.MagnitudeAngletoComplex.im = 0.0;
    boost_and_two_level__1_sm_ehs_B.RealImagtoComplex_j.re = 0.0;
    boost_and_two_level__1_sm_ehs_B.RealImagtoComplex_j.im = 0.0;
    boost_and_two_level__1_sm_ehs_B.RealImagtoComplex_f.re = 0.0;
    boost_and_two_level__1_sm_ehs_B.RealImagtoComplex_f.im = 0.0;
    boost_and_two_level__1_sm_ehs_B.SFunction = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum = 0.0;
    boost_and_two_level__1_sm_ehs_B.Outputs_eHS1_Recv_o2 = 0.0;
    boost_and_two_level__1_sm_ehs_B.UnitDelay = 0.0;
    boost_and_two_level__1_sm_ehs_B.DigitalClock = 0.0;
    boost_and_two_level__1_sm_ehs_B.Add1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.MathFunction = 0.0;
    boost_and_two_level__1_sm_ehs_B.ib1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.LookupTable = 0.0;
    boost_and_two_level__1_sm_ehs_B.Add3 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Add3_g = 0.0;
    boost_and_two_level__1_sm_ehs_B.Gain1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.MUL1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Add4 = 0.0;
    boost_and_two_level__1_sm_ehs_B.DataTypeConversion[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.DataTypeConversion[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.DataTypeConversion[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.DataTypeConversion[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Relay1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.DataTypeConversion8 = 0.0;
    boost_and_two_level__1_sm_ehs_B.eHS_rst_loadin = 0.0;
    boost_and_two_level__1_sm_ehs_B.Automated_Solver_Mat_Initialisa = 0.0;
    boost_and_two_level__1_sm_ehs_B.SineWaveFunction = 0.0;
    boost_and_two_level__1_sm_ehs_B.Memory1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Memory = 0.0;
    boost_and_two_level__1_sm_ehs_B.Inputs_eHS1_Send = 0.0;
    boost_and_two_level__1_sm_ehs_B.integrator = 0.0;
    boost_and_two_level__1_sm_ehs_B.TransportDelay = 0.0;
    boost_and_two_level__1_sm_ehs_B.Clock = 0.0;
    boost_and_two_level__1_sm_ehs_B.Memory_g = 0.0;
    boost_and_two_level__1_sm_ehs_B.Switch = 0.0;
    boost_and_two_level__1_sm_ehs_B.integrator_p = 0.0;
    boost_and_two_level__1_sm_ehs_B.TransportDelay_p = 0.0;
    boost_and_two_level__1_sm_ehs_B.Clock_m = 0.0;
    boost_and_two_level__1_sm_ehs_B.Memory_e = 0.0;
    boost_and_two_level__1_sm_ehs_B.Switch_e = 0.0;
    boost_and_two_level__1_sm_ehs_B.P_PV = 0.0;
    boost_and_two_level__1_sm_ehs_B.Outputs_eHS1_Recv_o2_j = 0.0;
    boost_and_two_level__1_sm_ehs_B.Computation_time = 0.0;
    boost_and_two_level__1_sm_ehs_B.Real_step_size = 0.0;
    boost_and_two_level__1_sm_ehs_B.Idle_time = 0.0;
    boost_and_two_level__1_sm_ehs_B.Num_overruns = 0.0;
    boost_and_two_level__1_sm_ehs_B.Memory1_o = 0.0;
    boost_and_two_level__1_sm_ehs_B.DataTypeConversion1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.DataTypeConversion_d = 0.0;
    boost_and_two_level__1_sm_ehs_B.Memory2 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Add = 0.0;
    boost_and_two_level__1_sm_ehs_B.Switch_i = 0.0;
    boost_and_two_level__1_sm_ehs_B.Switch1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.LoadIn = 0.0;
    boost_and_two_level__1_sm_ehs_B.DataInSend = 0.0;
    boost_and_two_level__1_sm_ehs_B.OpWriteFile_o1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.OpWriteFile_o2 = 0.0;
    boost_and_two_level__1_sm_ehs_B.DataTypeConversion1_o = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTEConversion1[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTEConversion1[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTEConversion1[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTEConversion1[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Step = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtRTEConvers[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtRTEConvers[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtRTEConvers[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtRTEConvers[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTEConversion[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTEConversion[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTEConversion[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTEConversion[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTELogicalOperator1[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTELogicalOperator1[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTELogicalOperator1[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTELogicalOperator1[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTEGround[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTEGround[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTEGround[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTEGround[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o1[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o1[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o1[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o1[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o2[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o2[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o2[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o2[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o3[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o3[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o3[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o3[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o4[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o4[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o4[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o4[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o5[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o5[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o5[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o5[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o6[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o6[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o6[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o6[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o7[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o7[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o7[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o7[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o8[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o8[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o8[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o8[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o9[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o9[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o9[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o9[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o10[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o10[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o10[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o10[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o11[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o11[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o11[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o11[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o12[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o12[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o12[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o12[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o13[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o13[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o13[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o13[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o14[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o14[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o14[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o14[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o15[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o15[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o15[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o15[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o16[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o16[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o16[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o16[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.EventGen_eHS_1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.RateTransition5 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Apu = 0.0;
    boost_and_two_level__1_sm_ehs_B.UnitDelay_l = 0.0;
    boost_and_two_level__1_sm_ehs_B.avoiddivisionbyzero = 0.0;
    boost_and_two_level__1_sm_ehs_B.MathFunction_i = 0.0;
    boost_and_two_level__1_sm_ehs_B.Gain = 0.0;
    boost_and_two_level__1_sm_ehs_B.SFunction_h = 0.0;
    boost_and_two_level__1_sm_ehs_B.DiscreteTimeIntegrator = 0.0;
    boost_and_two_level__1_sm_ehs_B.MathFunction_m = 0.0;
    boost_and_two_level__1_sm_ehs_B.FirstcycleofsimulationId092Iq0 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Switch_b[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Switch_b[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Duk[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Duk[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.x1k[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.x1k[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.C11[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.C11[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.x2k[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.x2k[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.C12[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.C12[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.sum2[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.sum2[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.yk[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.yk[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.UnitDelay2 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum_n[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum_n[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.ProportionalGain[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.ProportionalGain[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Integrator[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Integrator[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum_j[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum_j[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Saturate[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Saturate[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RateTransition4 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Vpu = 0.0;
    boost_and_two_level__1_sm_ehs_B.TrigonometricFunction = 0.0;
    boost_and_two_level__1_sm_ehs_B.Gain1_f = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product1_k = 0.0;
    boost_and_two_level__1_sm_ehs_B.Integ4 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Freq = 0.0;
    boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle = 0.0;
    boost_and_two_level__1_sm_ehs_B.RoundingFunction = 0.0;
    boost_and_two_level__1_sm_ehs_B.Delay = 0.0;
    boost_and_two_level__1_sm_ehs_B.SFunction_f = 0.0;
    boost_and_two_level__1_sm_ehs_B.UnitDelay_k = 0.0;
    boost_and_two_level__1_sm_ehs_B.DigitalClock_e = 0.0;
    boost_and_two_level__1_sm_ehs_B.UnitDelay1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Switch_ie = 0.0;
    boost_and_two_level__1_sm_ehs_B.TrigonometricFunction3 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Gain3 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product2 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Integ4_b = 0.0;
    boost_and_two_level__1_sm_ehs_B.Freq_g = 0.0;
    boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_h = 0.0;
    boost_and_two_level__1_sm_ehs_B.RoundingFunction_e = 0.0;
    boost_and_two_level__1_sm_ehs_B.Delay_p = 0.0;
    boost_and_two_level__1_sm_ehs_B.SFunction_fz = 0.0;
    boost_and_two_level__1_sm_ehs_B.UnitDelay_h = 0.0;
    boost_and_two_level__1_sm_ehs_B.DigitalClock_i = 0.0;
    boost_and_two_level__1_sm_ehs_B.UnitDelay1_c = 0.0;
    boost_and_two_level__1_sm_ehs_B.Switch_c = 0.0;
    boost_and_two_level__1_sm_ehs_B.ComplextoMagnitudeAngle_o1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.ComplextoMagnitudeAngle_o2 = 0.0;
    boost_and_two_level__1_sm_ehs_B.RadDeg = 0.0;
    boost_and_two_level__1_sm_ehs_B.torad = 0.0;
    boost_and_two_level__1_sm_ehs_B.ComplextoRealImag_o1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.ComplextoRealImag_o2 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Rff = 0.0;
    boost_and_two_level__1_sm_ehs_B.Lff = 0.0;
    boost_and_two_level__1_sm_ehs_B.Feedforward = 0.0;
    boost_and_two_level__1_sm_ehs_B.Rff_d = 0.0;
    boost_and_two_level__1_sm_ehs_B.Lff_o = 0.0;
    boost_and_two_level__1_sm_ehs_B.Add3_i = 0.0;
    boost_and_two_level__1_sm_ehs_B.Add2[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Add2[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.IntegralGain[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.IntegralGain[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Saturation[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Saturation[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.RateTransition7 = 0.0;
    boost_and_two_level__1_sm_ehs_B.RateTransition8 = 0.0;
    boost_and_two_level__1_sm_ehs_B.DigitalClock_j = 0.0;
    boost_and_two_level__1_sm_ehs_B.RateTransition6 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Integ4_m = 0.0;
    boost_and_two_level__1_sm_ehs_B.K1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.SFunction_c = 0.0;
    boost_and_two_level__1_sm_ehs_B.UnitDelay_b = 0.0;
    boost_and_two_level__1_sm_ehs_B.UnitDelay1_cm = 0.0;
    boost_and_two_level__1_sm_ehs_B.Switch_ej = 0.0;
    boost_and_two_level__1_sm_ehs_B.TrigonometricFunction2 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product1_g = 0.0;
    boost_and_two_level__1_sm_ehs_B.Integ4_e = 0.0;
    boost_and_two_level__1_sm_ehs_B.Freq_a = 0.0;
    boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_g = 0.0;
    boost_and_two_level__1_sm_ehs_B.RoundingFunction_l = 0.0;
    boost_and_two_level__1_sm_ehs_B.Delay_pz = 0.0;
    boost_and_two_level__1_sm_ehs_B.SFunction_i = 0.0;
    boost_and_two_level__1_sm_ehs_B.UnitDelay_lg = 0.0;
    boost_and_two_level__1_sm_ehs_B.DigitalClock_a = 0.0;
    boost_and_two_level__1_sm_ehs_B.UnitDelay1_i = 0.0;
    boost_and_two_level__1_sm_ehs_B.Switch_cd = 0.0;
    boost_and_two_level__1_sm_ehs_B.Divide_f = 0.0;
    boost_and_two_level__1_sm_ehs_B.DiscreteDerivative = 0.0;
    boost_and_two_level__1_sm_ehs_B.DiscreteTimeIntegrator_i = 0.0;
    boost_and_two_level__1_sm_ehs_B.Kp4 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum6 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Saturation1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Gain10 = 0.0;
    boost_and_two_level__1_sm_ehs_B.RateLimiter = 0.0;
    boost_and_two_level__1_sm_ehs_B.x1k_e = 0.0;
    boost_and_two_level__1_sm_ehs_B.A11 = 0.0;
    boost_and_two_level__1_sm_ehs_B.x2k_h = 0.0;
    boost_and_two_level__1_sm_ehs_B.A12 = 0.0;
    boost_and_two_level__1_sm_ehs_B.A21 = 0.0;
    boost_and_two_level__1_sm_ehs_B.A22 = 0.0;
    boost_and_two_level__1_sm_ehs_B.sum2_c = 0.0;
    boost_and_two_level__1_sm_ehs_B.sum3 = 0.0;
    boost_and_two_level__1_sm_ehs_B.B11 = 0.0;
    boost_and_two_level__1_sm_ehs_B.x1k1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.B21 = 0.0;
    boost_and_two_level__1_sm_ehs_B.x2k1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Duk_c = 0.0;
    boost_and_two_level__1_sm_ehs_B.C11_f = 0.0;
    boost_and_two_level__1_sm_ehs_B.C12_o = 0.0;
    boost_and_two_level__1_sm_ehs_B.sum2_n = 0.0;
    boost_and_two_level__1_sm_ehs_B.yk_c = 0.0;
    boost_and_two_level__1_sm_ehs_B.A11_n[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.A11_n[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.A12_e[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.A12_e[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.A21_b[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.A21_b[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.A22_o[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.A22_o[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.sum2_a[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.sum2_a[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.sum3_f[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.sum3_f[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.B11_a[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.B11_a[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.x1k1_l[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.x1k1_l[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.B21_a[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.B21_a[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.x2k1_j[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.x2k1_j[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Constant1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Switch_ex = 0.0;
    boost_and_two_level__1_sm_ehs_B.Add1_j = 0.0;
    boost_and_two_level__1_sm_ehs_B.UnitDelay3[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.UnitDelay3[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Gain1_p = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product1_m[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product1_m[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.ComplextoMagnitudeAngle_o1_k = 0.0;
    boost_and_two_level__1_sm_ehs_B.ComplextoMagnitudeAngle_o2_g = 0.0;
    boost_and_two_level__1_sm_ehs_B.Add2_l = 0.0;
    boost_and_two_level__1_sm_ehs_B.TrigonometricFunction_c = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product2_g = 0.0;
    boost_and_two_level__1_sm_ehs_B.UnitDelay1_k = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum_p = 0.0;
    boost_and_two_level__1_sm_ehs_B.Rtot_pu2 = 0.0;
    boost_and_two_level__1_sm_ehs_B.IntegralGain_m = 0.0;
    boost_and_two_level__1_sm_ehs_B.Integrator_j = 0.0;
    boost_and_two_level__1_sm_ehs_B.ProportionalGain_e = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum_e = 0.0;
    boost_and_two_level__1_sm_ehs_B.Saturate_d = 0.0;
    boost_and_two_level__1_sm_ehs_B.RateTransition3 = 0.0;
    boost_and_two_level__1_sm_ehs_B.RateTransition1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.RateTransition2 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Saturation_p = 0.0;
    boost_and_two_level__1_sm_ehs_B.Clock_o = 0.0;
    boost_and_two_level__1_sm_ehs_B.Constant1_f = 0.0;
    boost_and_two_level__1_sm_ehs_B.DataTypeConversion_c = 0.0;
    boost_and_two_level__1_sm_ehs_B.OpTrigger = 0.0;
    boost_and_two_level__1_sm_ehs_B.OpCtrl_o1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.OpCtrl_o2[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.OpCtrl_o2[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.OpCtrl_o2[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.OpCtrl_o2[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Iph = 0.0;
    boost_and_two_level__1_sm_ehs_B.Io = 0.0;
    boost_and_two_level__1_sm_ehs_B.Fcn = 0.0;
    boost_and_two_level__1_sm_ehs_B.Fcn1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Fcn_i = 0.0;
    boost_and_two_level__1_sm_ehs_B.Fcn1_k = 0.0;
    boost_and_two_level__1_sm_ehs_B.Switch_cq[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Switch_cq[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum5 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product5 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Gain1_e = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum4 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product2_f = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product4 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum7 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Meanvalue = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum5_h = 0.0;
    boost_and_two_level__1_sm_ehs_B.TrigonometricFunction_e = 0.0;
    boost_and_two_level__1_sm_ehs_B.Gain1_fj = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product1_a = 0.0;
    boost_and_two_level__1_sm_ehs_B.Integ4_ee = 0.0;
    boost_and_two_level__1_sm_ehs_B.Freq_h = 0.0;
    boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_o = 0.0;
    boost_and_two_level__1_sm_ehs_B.RoundingFunction_j = 0.0;
    boost_and_two_level__1_sm_ehs_B.Delay_o = 0.0;
    boost_and_two_level__1_sm_ehs_B.SFunction_b = 0.0;
    boost_and_two_level__1_sm_ehs_B.UnitDelay_e = 0.0;
    boost_and_two_level__1_sm_ehs_B.DigitalClock_g = 0.0;
    boost_and_two_level__1_sm_ehs_B.UnitDelay1_a = 0.0;
    boost_and_two_level__1_sm_ehs_B.Switch_h = 0.0;
    boost_and_two_level__1_sm_ehs_B.TrigonometricFunction3_k = 0.0;
    boost_and_two_level__1_sm_ehs_B.Gain3_b = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product2_b = 0.0;
    boost_and_two_level__1_sm_ehs_B.Integ4_i = 0.0;
    boost_and_two_level__1_sm_ehs_B.Freq_c = 0.0;
    boost_and_two_level__1_sm_ehs_B.Numberofsamplespercycle_e = 0.0;
    boost_and_two_level__1_sm_ehs_B.RoundingFunction_n = 0.0;
    boost_and_two_level__1_sm_ehs_B.Delay_e = 0.0;
    boost_and_two_level__1_sm_ehs_B.SFunction_j = 0.0;
    boost_and_two_level__1_sm_ehs_B.UnitDelay_m = 0.0;
    boost_and_two_level__1_sm_ehs_B.DigitalClock_d = 0.0;
    boost_and_two_level__1_sm_ehs_B.UnitDelay1_n = 0.0;
    boost_and_two_level__1_sm_ehs_B.Switch_cf = 0.0;
    boost_and_two_level__1_sm_ehs_B.ComplextoMagnitudeAngle_o1_c = 0.0;
    boost_and_two_level__1_sm_ehs_B.ComplextoMagnitudeAngle_o2_k = 0.0;
    boost_and_two_level__1_sm_ehs_B.RadDeg_l = 0.0;
    boost_and_two_level__1_sm_ehs_B.Saturation_o = 0.0;
    boost_and_two_level__1_sm_ehs_B.MathFunction_c = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum1_h = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum5_n = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product5_j = 0.0;
    boost_and_two_level__1_sm_ehs_B.Gain1_b = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum4_j = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product2_m = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product4_h = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum7_g = 0.0;
    boost_and_two_level__1_sm_ehs_B.Meanvalue_b = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum5_nt = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum1_i = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum5_a = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product5_i = 0.0;
    boost_and_two_level__1_sm_ehs_B.Gain1_p5 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum4_n = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product2_c = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product4_p = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum7_j = 0.0;
    boost_and_two_level__1_sm_ehs_B.Meanvalue_h = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum5_d = 0.0;
    boost_and_two_level__1_sm_ehs_B.Gain1_a = 0.0;
    boost_and_two_level__1_sm_ehs_B.Gain_o = 0.0;
    boost_and_two_level__1_sm_ehs_B.Correction = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum7_js = 0.0;
    boost_and_two_level__1_sm_ehs_B.Mean = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum5_c = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum1_e = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum5_p = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product5_i4 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Gain1_i = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum4_l = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product2_gq = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product4_n = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum7_m = 0.0;
    boost_and_two_level__1_sm_ehs_B.Meanvalue_a = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum5_an = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum1_g = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum5_f = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product5_l = 0.0;
    boost_and_two_level__1_sm_ehs_B.Gain1_pu = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum4_f = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product2_f1 = 0.0;
    boost_and_two_level__1_sm_ehs_B.Product4_k = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum7_p = 0.0;
    boost_and_two_level__1_sm_ehs_B.Meanvalue_o = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum5_n4 = 0.0;
    boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtSFunctionI[0] = 0.0;
    boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtSFunctionI[1] = 0.0;
    boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtSFunctionI[2] = 0.0;
    boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtSFunctionI[3] = 0.0;
    boost_and_two_level__1_sm_ehs_B.D = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum_ns = 0.0;
    boost_and_two_level__1_sm_ehs_B.Gain_k = 0.0;
    boost_and_two_level__1_sm_ehs_B.Sum_h = 0.0;
    boost_and_two_level__1_sm_ehs_B.Gain_a = 0.0;
  }

  /* parameters */
  boost_and_two_level__1_sm_ehs_M->ModelData.defaultParam = ((real_T *)
    &boost_and_two_level__1_sm_ehs_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &boost_and_two_level__1_sm_ehs_X;
    boost_and_two_level__1_sm_ehs_M->ModelData.contStates = (x);
    (void) memset((void *)&boost_and_two_level__1_sm_ehs_X, 0,
                  sizeof(X_boost_and_two_level__1_sm_ehs_T));
  }

  /* states (dwork) */
  boost_and_two_level__1_sm_ehs_M->ModelData.dwork = ((void *)
    &boost_and_two_level__1_sm_ehs_DW);
  (void) memset((void *)&boost_and_two_level__1_sm_ehs_DW, 0,
                sizeof(DW_boost_and_two_level__1_sm_ehs_T));
  boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE = 0.0;
  boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_b = 0.0;
  boost_and_two_level__1_sm_ehs_DW.DiscreteTimeIntegrator_DSTATE = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Delay_x1_DSTATE[0] = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Delay_x1_DSTATE[1] = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Delay_x2_DSTATE[0] = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Delay_x2_DSTATE[1] = 0.0;
  boost_and_two_level__1_sm_ehs_DW.UnitDelay2_DSTATE = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Integrator_DSTATE[0] = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Integrator_DSTATE[1] = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE = 0.0;
  boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_d = 0.0;
  boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_f = 0.0;
  boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_i = 0.0;
  boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_i = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_g = 0.0;
  boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_h = 0.0;
  boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_im = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_k = 0.0;
  boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_f = 0.0;
  boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_e = 0.0;
  boost_and_two_level__1_sm_ehs_DW.DiscreteDerivative_states = 0.0;
  boost_and_two_level__1_sm_ehs_DW.DiscreteTimeIntegrator_DSTATE_o = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Delay_x1_DSTATE_p = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Delay_x2_DSTATE_l = 0.0;
  boost_and_two_level__1_sm_ehs_DW.UnitDelay3_DSTATE[0] = 0.0;
  boost_and_two_level__1_sm_ehs_DW.UnitDelay3_DSTATE[1] = 0.0;
  boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_n = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Integrator_DSTATE_h = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_i = 0.0;
  boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_c = 0.0;
  boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_eb = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Integ4_DSTATE_c = 0.0;
  boost_and_two_level__1_sm_ehs_DW.UnitDelay_DSTATE_fv = 0.0;
  boost_and_two_level__1_sm_ehs_DW.UnitDelay1_DSTATE_k = 0.0;
  boost_and_two_level__1_sm_ehs_DW.SFunction_PreviousInput = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Memory1_PreviousInput = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Memory_PreviousInput = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Memory_PreviousInput_n = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Memory_PreviousInput_e = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Memory1_PreviousInput_h = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Memory2_PreviousInput = 0.0;
  boost_and_two_level__1_sm_ehs_DW.DiscreteDerivative_tmp = 0.0;
  boost_and_two_level__1_sm_ehs_DW.PrevY = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Vold = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Pold = 0.0;
  boost_and_two_level__1_sm_ehs_DW.Dold = 0.0;
  boost_and_two_level__1_sm_ehs_DW.TransportDelay_RWORK.modelTStart = 0.0;

  {
    int_T i;
    for (i = 0; i < 16384; i++) {
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_RWORK.TUbufferArea[i] =
        0.0;
    }
  }

  boost_and_two_level__1_sm_ehs_DW.TransportDelay_RWORK_o.modelTStart = 0.0;

  {
    int_T i;
    for (i = 0; i < 16384; i++) {
      boost_and_two_level__1_sm_ehs_DW.TransportDelay_RWORK_o.TUbufferArea[i] =
        0.0;
    }
  }

  boost_and_two_level__1_sm_ehs_DW.RTEConversion1_RWORK[0] = 0.0;
  boost_and_two_level__1_sm_ehs_DW.RTEConversion1_RWORK[1] = 0.0;
  boost_and_two_level__1_sm_ehs_DW.RTEConversion1_RWORK[2] = 0.0;
  boost_and_two_level__1_sm_ehs_DW.RTEConversion1_RWORK[3] = 0.0;
  boost_and_two_level__1_sm_ehs_DW.RTEConversion_RWORK[0] = 0.0;
  boost_and_two_level__1_sm_ehs_DW.RTEConversion_RWORK[1] = 0.0;
  boost_and_two_level__1_sm_ehs_DW.RTEConversion_RWORK[2] = 0.0;
  boost_and_two_level__1_sm_ehs_DW.RTEConversion_RWORK[3] = 0.0;
  boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK = 0.0;
  boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK_d = 0.0;
  boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK_g = 0.0;
  boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK_m = 0.0;
  boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK_m4 = 0.0;
  boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK_e = 0.0;
  boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK_p = 0.0;

  /* child S-Function registration */
  {
    RTWSfcnInfo *sfcnInfo =
      &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.sfcnInfo;
    boost_and_two_level__1_sm_ehs_M->sfcnInfo = (sfcnInfo);
    rtssSetErrorStatusPtr(sfcnInfo, (&rtmGetErrorStatus
      (boost_and_two_level__1_sm_ehs_M)));
    rtssSetNumRootSampTimesPtr(sfcnInfo,
      &boost_and_two_level__1_sm_ehs_M->Sizes.numSampTimes);
    boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.taskTimePtrs[0] =
      &(rtmGetTPtr(boost_and_two_level__1_sm_ehs_M)[0]);
    boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.taskTimePtrs[1] =
      &(rtmGetTPtr(boost_and_two_level__1_sm_ehs_M)[1]);
    rtssSetTPtrPtr(sfcnInfo,
                   boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.taskTimePtrs);
    rtssSetTStartPtr(sfcnInfo, &rtmGetTStart(boost_and_two_level__1_sm_ehs_M));
    rtssSetTFinalPtr(sfcnInfo, &rtmGetTFinal(boost_and_two_level__1_sm_ehs_M));
    rtssSetTimeOfLastOutputPtr(sfcnInfo, &rtmGetTimeOfLastOutput
      (boost_and_two_level__1_sm_ehs_M));
    rtssSetStepSizePtr(sfcnInfo,
                       &boost_and_two_level__1_sm_ehs_M->Timing.stepSize);
    rtssSetStopRequestedPtr(sfcnInfo, &rtmGetStopRequested
      (boost_and_two_level__1_sm_ehs_M));
    rtssSetDerivCacheNeedsResetPtr(sfcnInfo,
      &boost_and_two_level__1_sm_ehs_M->ModelData.derivCacheNeedsReset);
    rtssSetZCCacheNeedsResetPtr(sfcnInfo,
      &boost_and_two_level__1_sm_ehs_M->ModelData.zCCacheNeedsReset);
    rtssSetBlkStateChangePtr(sfcnInfo,
      &boost_and_two_level__1_sm_ehs_M->ModelData.blkStateChange);
    rtssSetSampleHitsPtr(sfcnInfo,
                         &boost_and_two_level__1_sm_ehs_M->Timing.sampleHits);
    rtssSetPerTaskSampleHitsPtr(sfcnInfo,
      &boost_and_two_level__1_sm_ehs_M->Timing.perTaskSampleHits);
    rtssSetSimModePtr(sfcnInfo, &boost_and_two_level__1_sm_ehs_M->simMode);
    rtssSetSolverInfoPtr(sfcnInfo,
                         &boost_and_two_level__1_sm_ehs_M->solverInfoPtr);
  }

  boost_and_two_level__1_sm_ehs_M->Sizes.numSFcns = (32);

  /* register each child */
  {
    (void) memset((void *)
                  &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.childSFunctions
                  [0], 0,
                  32*sizeof(SimStruct));
    boost_and_two_level__1_sm_ehs_M->childSfunctions =
      (&boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.childSFunctionPtrs[0]);

    {
      int_T i;
      for (i = 0; i < 32; i++) {
        boost_and_two_level__1_sm_ehs_M->childSfunctions[i] =
          (&boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.childSFunctions[i]);
      }
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S66>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[0];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn0.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn0.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn0.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [0]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [0]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [0]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [0]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn0.inputPortInfo
          [0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn0.UPtrs0;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.Integ4_ee;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn0.UPtrs1;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.Delay_o;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn0.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.SFunction_b));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn0.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P4_Size);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK_e);
      ssSetIWork(rts, (int_T *)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_IWORK_nb);
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_PWORK_lr);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn0.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn0.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK_e);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &boost_and_two_level__1_sm_ehs_DW.SFunction_IWORK_nb);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &boost_and_two_level__1_sm_ehs_DW.SFunction_PWORK_lr);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 1);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetInputPortWidth(rts, 1, 1);
      ssSetInputPortDataType(rts, 1, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 1, 0);
      ssSetInputPortFrameData(rts, 1, 0);
      ssSetOutputPortWidth(rts, 0, 1);
      ssSetOutputPortDataType(rts, 0, SS_DOUBLE);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumIWork(rts, 1);
      ssSetNumPWork(rts, 1);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S69>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[1];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn1.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn1.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn1.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [1]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [1]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [1]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [1]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn1.inputPortInfo
          [0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn1.UPtrs0;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.Integ4_i;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn1.UPtrs1;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.Delay_e;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn1.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.SFunction_j));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn1.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P1_Size_h);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P2_Size_p);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P3_Size_p);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P4_Size_p);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK_p);
      ssSetIWork(rts, (int_T *)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_IWORK_l);
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_PWORK_j);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn1.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn1.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK_p);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &boost_and_two_level__1_sm_ehs_DW.SFunction_IWORK_l);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &boost_and_two_level__1_sm_ehs_DW.SFunction_PWORK_j);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 1);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetInputPortWidth(rts, 1, 1);
      ssSetInputPortDataType(rts, 1, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 1, 0);
      ssSetInputPortFrameData(rts, 1, 0);
      ssSetOutputPortWidth(rts, 0, 1);
      ssSetOutputPortDataType(rts, 0, SS_DOUBLE);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumIWork(rts, 1);
      ssSetNumPWork(rts, 1);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S7>/Outputs_eHS1_Recv (sfun_fct_op7160ex1_recv) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[2];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn2.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn2.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn2.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [2]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [2]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [2]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [2]);
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn2.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 2);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 8);
          ssSetOutputPortSignal(rts, 0, ((uint32_T *)
            boost_and_two_level__1_sm_ehs_B.Outputs_eHS1_Recv_o1));
        }

        /* port 1 */
        {
          _ssSetOutputPortNumDimensions(rts, 1, 1);
          ssSetOutputPortWidth(rts, 1, 1);
          ssSetOutputPortSignal(rts, 1, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.Outputs_eHS1_Recv_o2));
        }
      }

      /* path info */
      ssSetModelName(rts, "Outputs_eHS1_Recv");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/eHSx64 Gen3 CommBlk/Outputs_eHS1_Recv");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn2.params;
        ssSetSFcnParamsCount(rts, 9);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Outputs_eHS1_Recv_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Outputs_eHS1_Recv_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Outputs_eHS1_Recv_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Outputs_eHS1_Recv_P4_Size);
        ssSetSFcnParam(rts, 4, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Outputs_eHS1_Recv_P5_Size);
        ssSetSFcnParam(rts, 5, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Outputs_eHS1_Recv_P6_Size);
        ssSetSFcnParam(rts, 6, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Outputs_eHS1_Recv_P7_Size);
        ssSetSFcnParam(rts, 7, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Outputs_eHS1_Recv_P8_Size);
        ssSetSFcnParam(rts, 8, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Outputs_eHS1_Recv_P9_Size);
      }

      /* work vectors */
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.Outputs_eHS1_Recv_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn2.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn2.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0,
                   &boost_and_two_level__1_sm_ehs_DW.Outputs_eHS1_Recv_PWORK);
      }

      /* registration */
      sfun_fct_op7160ex1_recv(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 1, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);
      _ssSetOutputPortBeingMerged(rts, 1, 0);

      /* Update the BufferDstPort flags for each input port */
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S7>/Convert  Single floating-point (FPGA)  to double (sfun_SFP2DBL) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[3];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn3.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn3.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn3.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [3]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [3]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [3]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [3]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn3.inputPortInfo
          [0]);

        /* port 0 */
        {
          uint32_T const **sfcnUPtrs = (uint32_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn3.UPtrs0;

          {
            int_T i1;
            const uint32_T *u0 =
              &boost_and_two_level__1_sm_ehs_B.Outputs_eHS1_Recv_o1[0];
            for (i1=0; i1 < 7; i1++) {
              sfcnUPtrs[i1] = &u0[i1];
            }
          }

          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 7);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn3.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 7);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.ConvertSinglefloatingpointFPGAt));
        }
      }

      /* path info */
      ssSetModelName(rts, "Convert \nSingle floating-point (FPGA) \nto double");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/eHSx64 Gen3 CommBlk/Convert  Single floating-point (FPGA)  to double");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* registration */
      sfun_SFP2DBL(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 7);
      ssSetInputPortDataType(rts, 0, SS_UINT32);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetOutputPortWidth(rts, 0, 7);
      ssSetOutputPortDataType(rts, 0, SS_DOUBLE);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S10>/S-Function (RECV_Param) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[4];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn4.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn4.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn4.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [4]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [4]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [4]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [4]);
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn4.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 164);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.SFunction_l));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/OpComm/Receive/S-Function");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn4.params;
        ssSetSFcnParamsCount(rts, 2);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P1_Size_f);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P2_Size_k);
      }

      /* registration */
      RECV_Param(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S7>/eHS_rst_loadin (sfun_fct_op7160ex1_load_in) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[5];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn5.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn5.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn5.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [5]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [5]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [5]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [5]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn5.inputPortInfo
          [0]);

        /* port 0 */
        {
          ssSetInputPortRequiredContiguous(rts, 0, 1);
          ssSetInputPortSignal(rts, 0,
                               &boost_and_two_level__1_sm_ehs_B.DataTypeConversion8);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          ssSetInputPortRequiredContiguous(rts, 1, 1);
          ssSetInputPortSignal(rts, 1, &boost_and_two_level__1_sm_ehs_B.Add_l);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn5.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.eHS_rst_loadin));
        }
      }

      /* path info */
      ssSetModelName(rts, "eHS_rst_loadin");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/eHSx64 Gen3 CommBlk/eHS_rst_loadin");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn5.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.eHS_rst_loadin_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.eHS_rst_loadin_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.eHS_rst_loadin_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.eHS_rst_loadin_P4_Size);
      }

      /* work vectors */
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.eHS_rst_loadin_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn5.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn5.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0,
                   &boost_and_two_level__1_sm_ehs_DW.eHS_rst_loadin_PWORK);
      }

      /* registration */
      sfun_fct_op7160ex1_load_in(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S7>/Automated_Solver_Mat_Initialisation_1 (sfun_efs_solver_cfg) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[6];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn6.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn6.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn6.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [6]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [6]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [6]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [6]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn6.inputPortInfo
          [0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn6.UPtrs0;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.DataTypeConversion1_o;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          uint32_T const **sfcnUPtrs = (uint32_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn6.UPtrs1;

          {
            int_T i1;
            const uint32_T *u1 = boost_and_two_level__1_sm_ehs_B.load_config1;
            for (i1=0; i1 < 10; i1++) {
              sfcnUPtrs[i1] = &u1[i1];
            }
          }

          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 10);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn6.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.Automated_Solver_Mat_Initialisa));
        }
      }

      /* path info */
      ssSetModelName(rts, "Automated_Solver_Mat_Initialisation_1");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/eHSx64 Gen3 CommBlk/Automated_Solver_Mat_Initialisation_1");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn6.params;
        ssSetSFcnParamsCount(rts, 3);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Automated_Solver_Mat_Initialisa);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Automated_Solver_Mat_Initiali_o);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Automated_Solver_Mat_Initiali_e);
      }

      /* work vectors */
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.Automated_Solver_Mat_Initialisa
                 [0]);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn6.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn6.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 2);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0,
                   &boost_and_two_level__1_sm_ehs_DW.Automated_Solver_Mat_Initialisa
                   [0]);
      }

      /* registration */
      sfun_efs_solver_cfg(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 1, 10);
      ssSetInputPortDataType(rts, 1, SS_UINT32);
      ssSetInputPortComplexSignal(rts, 1, 0);
      ssSetInputPortFrameData(rts, 1, 0);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S7>/Convert double to  Single floating-point (FPGA) (sfun_DBL2SFP) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[7];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn7.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn7.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn7.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [7]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [7]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [7]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [7]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn7.inputPortInfo
          [0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn7.UPtrs0;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.SineWaveFunction;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.Memory1;
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.Memory;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 3);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn7.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 3);
          ssSetOutputPortSignal(rts, 0, ((uint32_T *)
            boost_and_two_level__1_sm_ehs_B.ConvertdoubletoSinglefloatingpo));
        }
      }

      /* path info */
      ssSetModelName(rts, "Convert double to \nSingle floating-point (FPGA)");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/eHSx64 Gen3 CommBlk/Convert double to  Single floating-point (FPGA)");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* registration */
      sfun_DBL2SFP(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.0);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 3);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetOutputPortWidth(rts, 0, 3);
      ssSetOutputPortDataType(rts, 0, SS_UINT32);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S7>/Inputs_eHS1_Send (sfun_fct_op7160ex1_send) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[8];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn8.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn8.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn8.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [8]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [8]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [8]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [8]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn8.inputPortInfo
          [0]);

        /* port 0 */
        {
          ssSetInputPortRequiredContiguous(rts, 0, 1);
          ssSetInputPortSignal(rts, 0,
                               boost_and_two_level__1_sm_ehs_B.ConvertdoubletoSinglefloatingpo);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 3);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn8.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.Inputs_eHS1_Send));
        }
      }

      /* path info */
      ssSetModelName(rts, "Inputs_eHS1_Send");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/eHSx64 Gen3 CommBlk/Inputs_eHS1_Send");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn8.params;
        ssSetSFcnParamsCount(rts, 8);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Inputs_eHS1_Send_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Inputs_eHS1_Send_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Inputs_eHS1_Send_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Inputs_eHS1_Send_P4_Size);
        ssSetSFcnParam(rts, 4, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Inputs_eHS1_Send_P5_Size);
        ssSetSFcnParam(rts, 5, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Inputs_eHS1_Send_P6_Size);
        ssSetSFcnParam(rts, 6, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Inputs_eHS1_Send_P7_Size);
        ssSetSFcnParam(rts, 7, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Inputs_eHS1_Send_P8_Size);
      }

      /* work vectors */
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.Inputs_eHS1_Send_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn8.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn8.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0,
                   &boost_and_two_level__1_sm_ehs_DW.Inputs_eHS1_Send_PWORK);
      }

      /* registration */
      sfun_fct_op7160ex1_send(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S2>/Outputs_eHS1_Recv (sfun_fct_op7160ex1_recv) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[9];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn9.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn9.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn9.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [9]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [9]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [9]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [9]);
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn9.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 2);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 7);
          ssSetOutputPortSignal(rts, 0, ((uint32_T *)
            boost_and_two_level__1_sm_ehs_B.Outputs_eHS1_Recv_o1_h));
        }

        /* port 1 */
        {
          _ssSetOutputPortNumDimensions(rts, 1, 1);
          ssSetOutputPortWidth(rts, 1, 1);
          ssSetOutputPortSignal(rts, 1, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.Outputs_eHS1_Recv_o2_j));
        }
      }

      /* path info */
      ssSetModelName(rts, "Outputs_eHS1_Recv");
      ssSetPath(rts, "boost_and_two_level__1_sm_ehs/SM_eHS/Outputs_eHS1_Recv");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn9.params;
        ssSetSFcnParamsCount(rts, 9);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Outputs_eHS1_Recv_P1_Size_e);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Outputs_eHS1_Recv_P2_Size_g);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Outputs_eHS1_Recv_P3_Size_a);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Outputs_eHS1_Recv_P4_Size_b);
        ssSetSFcnParam(rts, 4, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Outputs_eHS1_Recv_P5_Size_m);
        ssSetSFcnParam(rts, 5, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Outputs_eHS1_Recv_P6_Size_b);
        ssSetSFcnParam(rts, 6, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Outputs_eHS1_Recv_P7_Size_j);
        ssSetSFcnParam(rts, 7, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Outputs_eHS1_Recv_P8_Size_p);
        ssSetSFcnParam(rts, 8, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.Outputs_eHS1_Recv_P9_Size_m);
      }

      /* work vectors */
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.Outputs_eHS1_Recv_PWORK_d);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn9.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn9.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0,
                   &boost_and_two_level__1_sm_ehs_DW.Outputs_eHS1_Recv_PWORK_d);
      }

      /* registration */
      sfun_fct_op7160ex1_recv(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 1, 0);
      _ssSetOutputPortBeingMerged(rts, 0, 0);
      _ssSetOutputPortBeingMerged(rts, 1, 0);

      /* Update the BufferDstPort flags for each input port */
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S2>/sfp2dbl (sfun_SFP2DBL) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[10];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn10.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn10.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn10.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [10]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [10]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [10]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [10]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn10.inputPortInfo[0]);

        /* port 0 */
        {
          uint32_T const **sfcnUPtrs = (uint32_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn10.UPtrs0;

          {
            int_T i1;
            const uint32_T *u0 =
              boost_and_two_level__1_sm_ehs_B.Outputs_eHS1_Recv_o1_h;
            for (i1=0; i1 < 7; i1++) {
              sfcnUPtrs[i1] = &u0[i1];
            }
          }

          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 7);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn10.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 7);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.fpga_raw));
        }
      }

      /* path info */
      ssSetModelName(rts, "sfp2dbl");
      ssSetPath(rts, "boost_and_two_level__1_sm_ehs/SM_eHS/sfp2dbl");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* registration */
      sfun_SFP2DBL(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 7);
      ssSetInputPortDataType(rts, 0, SS_UINT32);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetOutputPortWidth(rts, 0, 7);
      ssSetOutputPortDataType(rts, 0, SS_DOUBLE);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S2>/OpMonitor (opmonitor) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[11];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn11.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn11.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn11.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [11]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [11]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [11]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [11]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn11.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn11.UPtrs0;
          sfcnUPtrs[0] = (real_T*)&boost_and_two_level__1_sm_ehs_RGND;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn11.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 4);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.Computation_time));
        }

        /* port 1 */
        {
          _ssSetOutputPortNumDimensions(rts, 1, 1);
          ssSetOutputPortWidth(rts, 1, 1);
          ssSetOutputPortSignal(rts, 1, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.Real_step_size));
        }

        /* port 2 */
        {
          _ssSetOutputPortNumDimensions(rts, 2, 1);
          ssSetOutputPortWidth(rts, 2, 1);
          ssSetOutputPortSignal(rts, 2, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.Idle_time));
        }

        /* port 3 */
        {
          _ssSetOutputPortNumDimensions(rts, 3, 1);
          ssSetOutputPortWidth(rts, 3, 1);
          ssSetOutputPortSignal(rts, 3, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.Num_overruns));
        }
      }

      /* path info */
      ssSetModelName(rts, "OpMonitor");
      ssSetPath(rts, "boost_and_two_level__1_sm_ehs/SM_eHS/OpMonitor");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn11.params;
        ssSetSFcnParamsCount(rts, 6);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpMonitor_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpMonitor_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpMonitor_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpMonitor_P4_Size);
        ssSetSFcnParam(rts, 4, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpMonitor_P5_Size);
        ssSetSFcnParam(rts, 5, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpMonitor_P6_Size);
      }

      /* work vectors */
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.OpMonitor_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn11.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn11.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &boost_and_two_level__1_sm_ehs_DW.OpMonitor_PWORK);
      }

      /* registration */
      opmonitor(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 0);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 2, 1);
      _ssSetOutputPortConnected(rts, 3, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);
      _ssSetOutputPortBeingMerged(rts, 1, 0);
      _ssSetOutputPortBeingMerged(rts, 2, 0);
      _ssSetOutputPortBeingMerged(rts, 3, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S30>/S-Function (OP_SEND) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[12];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn12.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn12.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn12.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [12]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [12]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [12]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [12]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn12.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn12.UPtrs0;

          {
            int_T i1;
            const real_T *u0 = &boost_and_two_level__1_sm_ehs_B.Divide[0];
            for (i1=0; i1 < 7; i1++) {
              sfcnUPtrs[i1] = &u0[i1];
            }

            sfcnUPtrs[7] = &boost_and_two_level__1_sm_ehs_B.DataTypeConversion[0];
            sfcnUPtrs[8] = &boost_and_two_level__1_sm_ehs_B.SFunction_l[0];
            sfcnUPtrs[9] = &boost_and_two_level__1_sm_ehs_B.SFunction_l[1];
            sfcnUPtrs[10] = &boost_and_two_level__1_sm_ehs_B.eHS_rst_loadin;
            sfcnUPtrs[11] =
              &boost_and_two_level__1_sm_ehs_B.Automated_Solver_Mat_Initialisa;
            sfcnUPtrs[12] = &boost_and_two_level__1_sm_ehs_B.Inputs_eHS1_Send;
            sfcnUPtrs[13] =
              &boost_and_two_level__1_sm_ehs_B.Outputs_eHS1_Recv_o2;
            sfcnUPtrs[14] = &boost_and_two_level__1_sm_ehs_B.Switch;
            sfcnUPtrs[15] = &boost_and_two_level__1_sm_ehs_B.Switch_e;
            sfcnUPtrs[16] = &boost_and_two_level__1_sm_ehs_B.P_PV;
            sfcnUPtrs[17] = &boost_and_two_level__1_sm_ehs_B.SFunction_l[0];
            u0 = &boost_and_two_level__1_sm_ehs_B.fpga_raw[0];
            for (i1=0; i1 < 7; i1++) {
              sfcnUPtrs[i1+ 18] = &u0[i1];
            }

            sfcnUPtrs[25] = &boost_and_two_level__1_sm_ehs_B.Computation_time;
            sfcnUPtrs[26] = &boost_and_two_level__1_sm_ehs_B.Real_step_size;
            sfcnUPtrs[27] = &boost_and_two_level__1_sm_ehs_B.Idle_time;
            sfcnUPtrs[28] = &boost_and_two_level__1_sm_ehs_B.Num_overruns;
          }

          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 29);
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/rtlab_send_subsystem/Subsystem1/Send1/S-Function");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn12.params;
        ssSetSFcnParamsCount(rts, 1);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P1_Size_l);
      }

      /* work vectors */
      ssSetIWork(rts, (int_T *)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_IWORK[0]);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn12.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn12.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* IWORK */
        ssSetDWorkWidth(rts, 0, 5);
        ssSetDWorkDataType(rts, 0,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &boost_and_two_level__1_sm_ehs_DW.SFunction_IWORK[0]);
      }

      /* registration */
      OP_SEND(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 29);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S17>/Convert double to  Single floating-point (FPGA) (sfun_DBL2SFP) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[13];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn13.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn13.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn13.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [13]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [13]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [13]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [13]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn13.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn13.UPtrs0;

          {
            int_T i1;
            const real_T *u0 = &boost_and_two_level__1_sm_ehs_B.SFunction_l[36];
            for (i1=0; i1 < 16; i1++) {
              sfcnUPtrs[i1] = &u0[i1];
            }
          }

          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 16);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn13.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 16);
          ssSetOutputPortSignal(rts, 0, ((uint32_T *)
            boost_and_two_level__1_sm_ehs_B.ConvertdoubletoSinglefloating_o));
        }
      }

      /* path info */
      ssSetModelName(rts, "Convert double to \nSingle floating-point (FPGA)");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/To Analog Outputs/Aout_Mapping&Adjustments/Convert double to  Single floating-point (FPGA)");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* registration */
      sfun_DBL2SFP(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 16);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetOutputPortWidth(rts, 0, 16);
      ssSetOutputPortDataType(rts, 0, SS_UINT32);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S17>/Convert double to  Single floating-point (FPGA)1 (sfun_DBL2SFP) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[14];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn14.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn14.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn14.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [14]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [14]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [14]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [14]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn14.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn14.UPtrs0;

          {
            int_T i1;
            const real_T *u0 = &boost_and_two_level__1_sm_ehs_B.SFunction_l[68];
            for (i1=0; i1 < 16; i1++) {
              sfcnUPtrs[i1] = &u0[i1];
            }
          }

          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 16);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn14.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 16);
          ssSetOutputPortSignal(rts, 0, ((uint32_T *)
            boost_and_two_level__1_sm_ehs_B.ConvertdoubletoSinglefloatin_ou));
        }
      }

      /* path info */
      ssSetModelName(rts, "Convert double to \nSingle floating-point (FPGA)1");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/To Analog Outputs/Aout_Mapping&Adjustments/Convert double to  Single floating-point (FPGA)1");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* registration */
      sfun_DBL2SFP(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 16);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetOutputPortWidth(rts, 0, 16);
      ssSetOutputPortDataType(rts, 0, SS_UINT32);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S17>/LoadIn (sfun_fct_op7160ex1_load_in) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[15];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn15.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn15.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn15.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [15]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [15]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [15]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [15]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn15.inputPortInfo[0]);

        /* port 0 */
        {
          ssSetInputPortRequiredContiguous(rts, 0, 1);
          ssSetInputPortSignal(rts, 0,
                               &boost_and_two_level__1_sm_ehs_B.DataTypeConversion_d);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          ssSetInputPortRequiredContiguous(rts, 1, 1);
          ssSetInputPortSignal(rts, 1,
                               boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtLoadInInpo);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 18);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn15.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.LoadIn));
        }
      }

      /* path info */
      ssSetModelName(rts, "LoadIn");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/To Analog Outputs/Aout_Mapping&Adjustments/LoadIn");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn15.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.LoadIn_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.LoadIn_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.LoadIn_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.LoadIn_P4_Size);
      }

      /* work vectors */
      ssSetPWork(rts, (void **) &boost_and_two_level__1_sm_ehs_DW.LoadIn_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn15.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn15.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &boost_and_two_level__1_sm_ehs_DW.LoadIn_PWORK);
      }

      /* registration */
      sfun_fct_op7160ex1_load_in(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 0, 0);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S17>/Convert double to  Single floating-point (FPGA)2 (sfun_DBL2SFP) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[16];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn16.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn16.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn16.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [16]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [16]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [16]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [16]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn16.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn16.UPtrs0;

          {
            int_T i1;
            const real_T *u0 = &boost_and_two_level__1_sm_ehs_B.fpga_raw[0];
            for (i1=0; i1 < 7; i1++) {
              sfcnUPtrs[i1] = &u0[i1];
            }

            u0 = &boost_and_two_level__1_sm_ehs_B.Constant[0];
            for (i1=0; i1 < 9; i1++) {
              sfcnUPtrs[i1+ 7] = &u0[i1];
            }
          }

          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 16);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn16.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 16);
          ssSetOutputPortSignal(rts, 0, ((uint32_T *)
            boost_and_two_level__1_sm_ehs_B.ConvertdoubletoSinglefloating_e));
        }
      }

      /* path info */
      ssSetModelName(rts, "Convert double to \nSingle floating-point (FPGA)2");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/To Analog Outputs/Aout_Mapping&Adjustments/Convert double to  Single floating-point (FPGA)2");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* registration */
      sfun_DBL2SFP(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 16);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetOutputPortWidth(rts, 0, 16);
      ssSetOutputPortDataType(rts, 0, SS_UINT32);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S17>/DataIn Send (sfun_fct_op7160ex1_send) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[17];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn17.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn17.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn17.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [17]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [17]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [17]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [17]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn17.inputPortInfo[0]);

        /* port 0 */
        {
          ssSetInputPortRequiredContiguous(rts, 0, 1);
          ssSetInputPortSignal(rts, 0,
                               boost_and_two_level__1_sm_ehs_B.ConvertdoubletoSinglefloating_e);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 16);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn17.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.DataInSend));
        }
      }

      /* path info */
      ssSetModelName(rts, "DataIn Send");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/To Analog Outputs/Aout_Mapping&Adjustments/DataIn Send");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn17.params;
        ssSetSFcnParamsCount(rts, 8);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.DataInSend_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.DataInSend_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.DataInSend_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.DataInSend_P4_Size);
        ssSetSFcnParam(rts, 4, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.DataInSend_P5_Size);
        ssSetSFcnParam(rts, 5, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.DataInSend_P6_Size);
        ssSetSFcnParam(rts, 6, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.DataInSend_P7_Size);
        ssSetSFcnParam(rts, 7, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.DataInSend_P8_Size);
      }

      /* work vectors */
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.DataInSend_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn17.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn17.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &boost_and_two_level__1_sm_ehs_DW.DataInSend_PWORK);
      }

      /* registration */
      sfun_fct_op7160ex1_send(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 0, 0);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S2>/OpWriteFile (opwritefile) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[18];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn18.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn18.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn18.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [18]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [18]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [18]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [18]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn18.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn18.UPtrs0;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.Switch;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.Switch_e;
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.P_PV;
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.SFunction_l[0];
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 4);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn18.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 2);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.OpWriteFile_o1));
        }

        /* port 1 */
        {
          _ssSetOutputPortNumDimensions(rts, 1, 1);
          ssSetOutputPortWidth(rts, 1, 1);
          ssSetOutputPortSignal(rts, 1, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.OpWriteFile_o2));
        }
      }

      /* path info */
      ssSetModelName(rts, "OpWriteFile");
      ssSetPath(rts, "boost_and_two_level__1_sm_ehs/SM_eHS/OpWriteFile");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn18.params;
        ssSetSFcnParamsCount(rts, 12);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpWriteFile_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpWriteFile_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpWriteFile_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpWriteFile_P4_Size);
        ssSetSFcnParam(rts, 4, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpWriteFile_P5_Size);
        ssSetSFcnParam(rts, 5, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpWriteFile_P6_Size);
        ssSetSFcnParam(rts, 6, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpWriteFile_P7_Size);
        ssSetSFcnParam(rts, 7, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpWriteFile_P8_Size);
        ssSetSFcnParam(rts, 8, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpWriteFile_P9_Size);
        ssSetSFcnParam(rts, 9, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpWriteFile_P10_Size);
        ssSetSFcnParam(rts, 10, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpWriteFile_P11_Size);
        ssSetSFcnParam(rts, 11, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpWriteFile_P12_Size);
      }

      /* work vectors */
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.OpWriteFile_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn18.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn18.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &boost_and_two_level__1_sm_ehs_DW.OpWriteFile_PWORK);
      }

      /* registration */
      opwritefile(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 4);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 0, 0);
      _ssSetOutputPortConnected(rts, 1, 0);
      _ssSetOutputPortBeingMerged(rts, 0, 0);
      _ssSetOutputPortBeingMerged(rts, 1, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S9>/RTE Conversion1 (rte_conversion) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[19];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn19.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn19.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn19.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [19]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [19]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [19]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [19]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn19.inputPortInfo[0]);

        /* port 0 */
        {
          ssSetInputPortRequiredContiguous(rts, 0, 1);
          ssSetInputPortSignal(rts, 0,
                               boost_and_two_level__1_sm_ehs_B.DataTypeConversion);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 4);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn19.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 4);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTEConversion1));
        }
      }

      /* path info */
      ssSetModelName(rts, "RTE Conversion1");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/RTE Conversion1");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn19.params;
        ssSetSFcnParamsCount(rts, 5);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTEConversion1_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTEConversion1_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTEConversion1_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTEConversion1_P4_Size);
        ssSetSFcnParam(rts, 4, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTEConversion1_P5_Size);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *)
                 &boost_and_two_level__1_sm_ehs_DW.RTEConversion1_RWORK[0]);
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.RTEConversion1_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn19.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn19.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 2);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 4);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0,
                   &boost_and_two_level__1_sm_ehs_DW.RTEConversion1_RWORK[0]);

        /* PWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1,
                   &boost_and_two_level__1_sm_ehs_DW.RTEConversion1_PWORK);
      }

      /* registration */
      rte_conversion(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 4);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetOutputPortWidth(rts, 0, 4);
      ssSetOutputPortDataType(rts, 0, SS_DOUBLE);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumRWork(rts, 4);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S9>/RTE Conversion (rte_conversion) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[20];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn20.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn20.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn20.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [20]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [20]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [20]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [20]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn20.inputPortInfo[0]);

        /* port 0 */
        {
          ssSetInputPortRequiredContiguous(rts, 0, 1);
          ssSetInputPortSignal(rts, 0,
                               boost_and_two_level__1_sm_ehs_B.TmpSignalConversionAtRTEConvers);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 4);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn20.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 4);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTEConversion));
        }
      }

      /* path info */
      ssSetModelName(rts, "RTE Conversion");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/RTE Conversion");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn20.params;
        ssSetSFcnParamsCount(rts, 5);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTEConversion_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTEConversion_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTEConversion_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTEConversion_P4_Size);
        ssSetSFcnParam(rts, 4, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTEConversion_P5_Size);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *)
                 &boost_and_two_level__1_sm_ehs_DW.RTEConversion_RWORK[0]);
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.RTEConversion_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn20.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn20.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 2);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 4);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0,
                   &boost_and_two_level__1_sm_ehs_DW.RTEConversion_RWORK[0]);

        /* PWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &boost_and_two_level__1_sm_ehs_DW.RTEConversion_PWORK);
      }

      /* registration */
      rte_conversion(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 4);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetOutputPortWidth(rts, 0, 4);
      ssSetOutputPortDataType(rts, 0, SS_DOUBLE);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumRWork(rts, 4);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S9>/RTE Logical Operator1 (rte_logical_operator) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[21];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn21.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn21.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn21.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [21]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [21]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [21]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [21]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn21.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn21.UPtrs0;
          sfcnUPtrs[0] = boost_and_two_level__1_sm_ehs_B.RTEConversion1;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.RTEConversion1[1];
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.RTEConversion1[2];
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.RTEConversion1[3];
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 4);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn21.UPtrs1;
          sfcnUPtrs[0] = boost_and_two_level__1_sm_ehs_B.RTEConversion;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.RTEConversion[1];
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.RTEConversion[2];
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.RTEConversion[3];
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 4);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn21.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 4);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTELogicalOperator1));
        }
      }

      /* path info */
      ssSetModelName(rts, "RTE Logical Operator1");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/RTE Logical Operator1");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn21.params;
        ssSetSFcnParamsCount(rts, 3);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTELogicalOperator1_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTELogicalOperator1_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTELogicalOperator1_P3_Size);
      }

      /* work vectors */
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.RTELogicalOperator1_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn21.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn21.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0,
                   &boost_and_two_level__1_sm_ehs_DW.RTELogicalOperator1_PWORK);
      }

      /* registration */
      rte_logical_operator(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 4);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetInputPortWidth(rts, 1, 4);
      ssSetInputPortDataType(rts, 1, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 1, 0);
      ssSetInputPortFrameData(rts, 1, 0);
      ssSetOutputPortWidth(rts, 0, 4);
      ssSetOutputPortDataType(rts, 0, SS_DOUBLE);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S9>/RTE Ground (rte_ground) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[22];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn22.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn22.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn22.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [22]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [22]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [22]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [22]);
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn22.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 4);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTEGround));
        }
      }

      /* path info */
      ssSetModelName(rts, "RTE Ground");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/RTE Ground");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn22.params;
        ssSetSFcnParamsCount(rts, 1);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTEGround_P1_Size);
      }

      /* work vectors */
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.RTEGround_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn22.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn22.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &boost_and_two_level__1_sm_ehs_DW.RTEGround_PWORK);
      }

      /* registration */
      rte_ground(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetOutputPortWidth(rts, 0, 4);
      ssSetOutputPortDataType(rts, 0, SS_DOUBLE);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S26>/RTE_Conversion_1 (rte_conversion_ophsdio) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[23];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn23.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn23.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn23.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [23]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [23]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [23]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [23]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn23.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn23.UPtrs0;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.RTELogicalOperator1[0];
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.RTELogicalOperator1[1];
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.RTELogicalOperator1[2];
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.RTELogicalOperator1[3];
          sfcnUPtrs[4] = &boost_and_two_level__1_sm_ehs_B.RTEGround[0];
          sfcnUPtrs[5] = &boost_and_two_level__1_sm_ehs_B.RTEGround[1];
          sfcnUPtrs[6] = &boost_and_two_level__1_sm_ehs_B.RTEGround[2];
          sfcnUPtrs[7] = &boost_and_two_level__1_sm_ehs_B.RTEGround[3];
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 8);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn23.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 16);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 4);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o1));
        }

        /* port 1 */
        {
          _ssSetOutputPortNumDimensions(rts, 1, 1);
          ssSetOutputPortWidth(rts, 1, 4);
          ssSetOutputPortSignal(rts, 1, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o2));
        }

        /* port 2 */
        {
          _ssSetOutputPortNumDimensions(rts, 2, 1);
          ssSetOutputPortWidth(rts, 2, 4);
          ssSetOutputPortSignal(rts, 2, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o3));
        }

        /* port 3 */
        {
          _ssSetOutputPortNumDimensions(rts, 3, 1);
          ssSetOutputPortWidth(rts, 3, 4);
          ssSetOutputPortSignal(rts, 3, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o4));
        }

        /* port 4 */
        {
          _ssSetOutputPortNumDimensions(rts, 4, 1);
          ssSetOutputPortWidth(rts, 4, 4);
          ssSetOutputPortSignal(rts, 4, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o5));
        }

        /* port 5 */
        {
          _ssSetOutputPortNumDimensions(rts, 5, 1);
          ssSetOutputPortWidth(rts, 5, 4);
          ssSetOutputPortSignal(rts, 5, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o6));
        }

        /* port 6 */
        {
          _ssSetOutputPortNumDimensions(rts, 6, 1);
          ssSetOutputPortWidth(rts, 6, 4);
          ssSetOutputPortSignal(rts, 6, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o7));
        }

        /* port 7 */
        {
          _ssSetOutputPortNumDimensions(rts, 7, 1);
          ssSetOutputPortWidth(rts, 7, 4);
          ssSetOutputPortSignal(rts, 7, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o8));
        }

        /* port 8 */
        {
          _ssSetOutputPortNumDimensions(rts, 8, 1);
          ssSetOutputPortWidth(rts, 8, 4);
          ssSetOutputPortSignal(rts, 8, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o9));
        }

        /* port 9 */
        {
          _ssSetOutputPortNumDimensions(rts, 9, 1);
          ssSetOutputPortWidth(rts, 9, 4);
          ssSetOutputPortSignal(rts, 9, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o10));
        }

        /* port 10 */
        {
          _ssSetOutputPortNumDimensions(rts, 10, 1);
          ssSetOutputPortWidth(rts, 10, 4);
          ssSetOutputPortSignal(rts, 10, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o11));
        }

        /* port 11 */
        {
          _ssSetOutputPortNumDimensions(rts, 11, 1);
          ssSetOutputPortWidth(rts, 11, 4);
          ssSetOutputPortSignal(rts, 11, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o12));
        }

        /* port 12 */
        {
          _ssSetOutputPortNumDimensions(rts, 12, 1);
          ssSetOutputPortWidth(rts, 12, 4);
          ssSetOutputPortSignal(rts, 12, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o13));
        }

        /* port 13 */
        {
          _ssSetOutputPortNumDimensions(rts, 13, 1);
          ssSetOutputPortWidth(rts, 13, 4);
          ssSetOutputPortSignal(rts, 13, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o14));
        }

        /* port 14 */
        {
          _ssSetOutputPortNumDimensions(rts, 14, 1);
          ssSetOutputPortWidth(rts, 14, 4);
          ssSetOutputPortSignal(rts, 14, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o15));
        }

        /* port 15 */
        {
          _ssSetOutputPortNumDimensions(rts, 15, 1);
          ssSetOutputPortWidth(rts, 15, 4);
          ssSetOutputPortSignal(rts, 15, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o16));
        }
      }

      /* path info */
      ssSetModelName(rts, "RTE_Conversion_1");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/eHSx64 Gen3 CommBlk/c_solver_rte1/RTE_Conversion_1");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn23.params;
        ssSetSFcnParamsCount(rts, 7);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTE_Conversion_1_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTE_Conversion_1_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTE_Conversion_1_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTE_Conversion_1_P4_Size);
        ssSetSFcnParam(rts, 4, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTE_Conversion_1_P5_Size);
        ssSetSFcnParam(rts, 5, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTE_Conversion_1_P6_Size);
        ssSetSFcnParam(rts, 6, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.RTE_Conversion_1_P7_Size);
      }

      /* work vectors */
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.RTE_Conversion_1_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn23.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn23.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0,
                   &boost_and_two_level__1_sm_ehs_DW.RTE_Conversion_1_PWORK);
      }

      /* registration */
      rte_conversion_ophsdio(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 2, 1);
      _ssSetOutputPortConnected(rts, 3, 1);
      _ssSetOutputPortConnected(rts, 4, 1);
      _ssSetOutputPortConnected(rts, 5, 1);
      _ssSetOutputPortConnected(rts, 6, 1);
      _ssSetOutputPortConnected(rts, 7, 1);
      _ssSetOutputPortConnected(rts, 8, 1);
      _ssSetOutputPortConnected(rts, 9, 1);
      _ssSetOutputPortConnected(rts, 10, 1);
      _ssSetOutputPortConnected(rts, 11, 1);
      _ssSetOutputPortConnected(rts, 12, 1);
      _ssSetOutputPortConnected(rts, 13, 1);
      _ssSetOutputPortConnected(rts, 14, 1);
      _ssSetOutputPortConnected(rts, 15, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);
      _ssSetOutputPortBeingMerged(rts, 1, 0);
      _ssSetOutputPortBeingMerged(rts, 2, 0);
      _ssSetOutputPortBeingMerged(rts, 3, 0);
      _ssSetOutputPortBeingMerged(rts, 4, 0);
      _ssSetOutputPortBeingMerged(rts, 5, 0);
      _ssSetOutputPortBeingMerged(rts, 6, 0);
      _ssSetOutputPortBeingMerged(rts, 7, 0);
      _ssSetOutputPortBeingMerged(rts, 8, 0);
      _ssSetOutputPortBeingMerged(rts, 9, 0);
      _ssSetOutputPortBeingMerged(rts, 10, 0);
      _ssSetOutputPortBeingMerged(rts, 11, 0);
      _ssSetOutputPortBeingMerged(rts, 12, 0);
      _ssSetOutputPortBeingMerged(rts, 13, 0);
      _ssSetOutputPortBeingMerged(rts, 14, 0);
      _ssSetOutputPortBeingMerged(rts, 15, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S26>/EventGen_eHS_1 (sfun_op7160ex1_event_generator) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[24];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [24]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [24]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [24]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [24]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 16);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn24.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.UPtrs0;
          sfcnUPtrs[0] = boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o1;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o1[1];
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o1[2];
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o1[3];
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 4);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.UPtrs1;
          sfcnUPtrs[0] = boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o2;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o2[1];
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o2[2];
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o2[3];
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 4);
        }

        /* port 2 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.UPtrs2;
          sfcnUPtrs[0] = boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o3;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o3[1];
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o3[2];
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o3[3];
          ssSetInputPortSignalPtrs(rts, 2, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 2, 1);
          ssSetInputPortWidth(rts, 2, 4);
        }

        /* port 3 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.UPtrs3;
          sfcnUPtrs[0] = boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o4;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o4[1];
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o4[2];
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o4[3];
          ssSetInputPortSignalPtrs(rts, 3, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 3, 1);
          ssSetInputPortWidth(rts, 3, 4);
        }

        /* port 4 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.UPtrs4;
          sfcnUPtrs[0] = boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o5;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o5[1];
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o5[2];
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o5[3];
          ssSetInputPortSignalPtrs(rts, 4, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 4, 1);
          ssSetInputPortWidth(rts, 4, 4);
        }

        /* port 5 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.UPtrs5;
          sfcnUPtrs[0] = boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o6;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o6[1];
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o6[2];
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o6[3];
          ssSetInputPortSignalPtrs(rts, 5, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 5, 1);
          ssSetInputPortWidth(rts, 5, 4);
        }

        /* port 6 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.UPtrs6;
          sfcnUPtrs[0] = boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o7;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o7[1];
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o7[2];
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o7[3];
          ssSetInputPortSignalPtrs(rts, 6, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 6, 1);
          ssSetInputPortWidth(rts, 6, 4);
        }

        /* port 7 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.UPtrs7;
          sfcnUPtrs[0] = boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o8;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o8[1];
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o8[2];
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o8[3];
          ssSetInputPortSignalPtrs(rts, 7, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 7, 1);
          ssSetInputPortWidth(rts, 7, 4);
        }

        /* port 8 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.UPtrs8;
          sfcnUPtrs[0] = boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o9;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o9[1];
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o9[2];
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o9[3];
          ssSetInputPortSignalPtrs(rts, 8, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 8, 1);
          ssSetInputPortWidth(rts, 8, 4);
        }

        /* port 9 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.UPtrs9;
          sfcnUPtrs[0] = boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o10;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o10[1];
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o10[2];
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o10[3];
          ssSetInputPortSignalPtrs(rts, 9, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 9, 1);
          ssSetInputPortWidth(rts, 9, 4);
        }

        /* port 10 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.UPtrs10;
          sfcnUPtrs[0] = boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o11;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o11[1];
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o11[2];
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o11[3];
          ssSetInputPortSignalPtrs(rts, 10, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 10, 1);
          ssSetInputPortWidth(rts, 10, 4);
        }

        /* port 11 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.UPtrs11;
          sfcnUPtrs[0] = boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o12;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o12[1];
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o12[2];
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o12[3];
          ssSetInputPortSignalPtrs(rts, 11, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 11, 1);
          ssSetInputPortWidth(rts, 11, 4);
        }

        /* port 12 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.UPtrs12;
          sfcnUPtrs[0] = boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o13;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o13[1];
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o13[2];
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o13[3];
          ssSetInputPortSignalPtrs(rts, 12, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 12, 1);
          ssSetInputPortWidth(rts, 12, 4);
        }

        /* port 13 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.UPtrs13;
          sfcnUPtrs[0] = boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o14;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o14[1];
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o14[2];
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o14[3];
          ssSetInputPortSignalPtrs(rts, 13, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 13, 1);
          ssSetInputPortWidth(rts, 13, 4);
        }

        /* port 14 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.UPtrs14;
          sfcnUPtrs[0] = boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o15;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o15[1];
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o15[2];
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o15[3];
          ssSetInputPortSignalPtrs(rts, 14, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 14, 1);
          ssSetInputPortWidth(rts, 14, 4);
        }

        /* port 15 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.UPtrs15;
          sfcnUPtrs[0] = boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o16;
          sfcnUPtrs[1] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o16[1];
          sfcnUPtrs[2] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o16[2];
          sfcnUPtrs[3] = &boost_and_two_level__1_sm_ehs_B.RTE_Conversion_1_o16[3];
          ssSetInputPortSignalPtrs(rts, 15, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 15, 1);
          ssSetInputPortWidth(rts, 15, 4);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.EventGen_eHS_1));
        }
      }

      /* path info */
      ssSetModelName(rts, "EventGen_eHS_1");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/eHSx64 Gen3 CommBlk/c_solver_rte1/EventGen_eHS_1");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.params;
        ssSetSFcnParamsCount(rts, 7);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.EventGen_eHS_1_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.EventGen_eHS_1_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.EventGen_eHS_1_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.EventGen_eHS_1_P4_Size);
        ssSetSFcnParam(rts, 4, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.EventGen_eHS_1_P5_Size);
        ssSetSFcnParam(rts, 5, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.EventGen_eHS_1_P6_Size);
        ssSetSFcnParam(rts, 6, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.EventGen_eHS_1_P7_Size);
      }

      /* work vectors */
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.EventGen_eHS_1_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn24.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0,
                   &boost_and_two_level__1_sm_ehs_DW.EventGen_eHS_1_PWORK);
      }

      /* registration */
      sfun_op7160ex1_event_generator(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetInputPortConnected(rts, 2, 1);
      _ssSetInputPortConnected(rts, 3, 1);
      _ssSetInputPortConnected(rts, 4, 1);
      _ssSetInputPortConnected(rts, 5, 1);
      _ssSetInputPortConnected(rts, 6, 1);
      _ssSetInputPortConnected(rts, 7, 1);
      _ssSetInputPortConnected(rts, 8, 1);
      _ssSetInputPortConnected(rts, 9, 1);
      _ssSetInputPortConnected(rts, 10, 1);
      _ssSetInputPortConnected(rts, 11, 1);
      _ssSetInputPortConnected(rts, 12, 1);
      _ssSetInputPortConnected(rts, 13, 1);
      _ssSetInputPortConnected(rts, 14, 1);
      _ssSetInputPortConnected(rts, 15, 1);
      _ssSetOutputPortConnected(rts, 0, 0);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
      ssSetInputPortBufferDstPort(rts, 2, -1);
      ssSetInputPortBufferDstPort(rts, 3, -1);
      ssSetInputPortBufferDstPort(rts, 4, -1);
      ssSetInputPortBufferDstPort(rts, 5, -1);
      ssSetInputPortBufferDstPort(rts, 6, -1);
      ssSetInputPortBufferDstPort(rts, 7, -1);
      ssSetInputPortBufferDstPort(rts, 8, -1);
      ssSetInputPortBufferDstPort(rts, 9, -1);
      ssSetInputPortBufferDstPort(rts, 10, -1);
      ssSetInputPortBufferDstPort(rts, 11, -1);
      ssSetInputPortBufferDstPort(rts, 12, -1);
      ssSetInputPortBufferDstPort(rts, 13, -1);
      ssSetInputPortBufferDstPort(rts, 14, -1);
      ssSetInputPortBufferDstPort(rts, 15, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S82>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[25];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn25.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn25.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn25.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [25]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [25]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [25]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [25]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn25.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn25.UPtrs0;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.Apu;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn25.UPtrs1;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.Gain;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn25.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.SFunction_h));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn25.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P1_Size_i);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P2_Size_i);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P3_Size_i);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P4_Size_n);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK);
      ssSetIWork(rts, (int_T *)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_IWORK_n);
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn25.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn25.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &boost_and_two_level__1_sm_ehs_DW.SFunction_IWORK_n);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &boost_and_two_level__1_sm_ehs_DW.SFunction_PWORK);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 1);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetInputPortWidth(rts, 1, 1);
      ssSetInputPortDataType(rts, 1, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 1, 0);
      ssSetInputPortFrameData(rts, 1, 0);
      ssSetOutputPortWidth(rts, 0, 1);
      ssSetOutputPortDataType(rts, 0, SS_DOUBLE);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumIWork(rts, 1);
      ssSetNumPWork(rts, 1);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S50>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[26];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn26.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn26.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn26.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [26]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [26]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [26]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [26]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn26.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn26.UPtrs0;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.Integ4;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn26.UPtrs1;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.Delay;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn26.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.SFunction_f));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn26.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P1_Size_ll);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P2_Size_c);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P3_Size_b);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P4_Size_f);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK_d);
      ssSetIWork(rts, (int_T *)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_IWORK_b);
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_PWORK_o);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn26.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn26.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK_d);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &boost_and_two_level__1_sm_ehs_DW.SFunction_IWORK_b);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &boost_and_two_level__1_sm_ehs_DW.SFunction_PWORK_o);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 1);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetInputPortWidth(rts, 1, 1);
      ssSetInputPortDataType(rts, 1, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 1, 0);
      ssSetInputPortFrameData(rts, 1, 0);
      ssSetOutputPortWidth(rts, 0, 1);
      ssSetOutputPortDataType(rts, 0, SS_DOUBLE);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumIWork(rts, 1);
      ssSetNumPWork(rts, 1);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S53>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[27];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn27.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn27.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn27.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [27]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [27]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [27]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [27]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn27.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn27.UPtrs0;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.Integ4_b;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn27.UPtrs1;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.Delay_p;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn27.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.SFunction_fz));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn27.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P1_Size_hi);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P2_Size_j);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P3_Size_c);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P4_Size_m);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK_g);
      ssSetIWork(rts, (int_T *)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_IWORK_j);
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_PWORK_e);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn27.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn27.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK_g);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &boost_and_two_level__1_sm_ehs_DW.SFunction_IWORK_j);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &boost_and_two_level__1_sm_ehs_DW.SFunction_PWORK_e);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 1);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetInputPortWidth(rts, 1, 1);
      ssSetInputPortDataType(rts, 1, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 1, 0);
      ssSetInputPortFrameData(rts, 1, 0);
      ssSetOutputPortWidth(rts, 0, 1);
      ssSetOutputPortDataType(rts, 0, SS_DOUBLE);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumIWork(rts, 1);
      ssSetNumPWork(rts, 1);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S55>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[28];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn28.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn28.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn28.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [28]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [28]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [28]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [28]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn28.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn28.UPtrs0;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.Integ4_m;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn28.UPtrs1;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.K1;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn28.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.SFunction_c));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Mean/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn28.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P1_Size_o);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P2_Size_id);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P3_Size_l);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P4_Size_i);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK_m);
      ssSetIWork(rts, (int_T *)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_IWORK_o);
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_PWORK_b);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn28.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn28.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK_m);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &boost_and_two_level__1_sm_ehs_DW.SFunction_IWORK_o);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &boost_and_two_level__1_sm_ehs_DW.SFunction_PWORK_b);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 1);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetInputPortWidth(rts, 1, 1);
      ssSetInputPortDataType(rts, 1, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 1, 0);
      ssSetInputPortFrameData(rts, 1, 0);
      ssSetOutputPortWidth(rts, 0, 1);
      ssSetOutputPortDataType(rts, 0, SS_DOUBLE);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumIWork(rts, 1);
      ssSetNumPWork(rts, 1);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S72>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[29];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn29.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn29.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn29.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [29]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [29]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [29]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [29]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn29.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn29.UPtrs0;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.Integ4_e;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn29.UPtrs1;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.Delay_pz;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn29.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.SFunction_i));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Mean (Variable Frequency)/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn29.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P1_Size_lc);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P2_Size_g);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P3_Size_f);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.SFunction_P4_Size_a);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK_m4);
      ssSetIWork(rts, (int_T *)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_IWORK_k);
      ssSetPWork(rts, (void **)
                 &boost_and_two_level__1_sm_ehs_DW.SFunction_PWORK_l);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn29.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn29.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &boost_and_two_level__1_sm_ehs_DW.SFunction_RWORK_m4);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &boost_and_two_level__1_sm_ehs_DW.SFunction_IWORK_k);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &boost_and_two_level__1_sm_ehs_DW.SFunction_PWORK_l);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 1);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetInputPortWidth(rts, 1, 1);
      ssSetInputPortDataType(rts, 1, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 1, 0);
      ssSetInputPortFrameData(rts, 1, 0);
      ssSetOutputPortWidth(rts, 0, 1);
      ssSetOutputPortDataType(rts, 0, SS_DOUBLE);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumIWork(rts, 1);
      ssSetNumPWork(rts, 1);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S2>/OpTrigger (optrigger) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[30];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn30.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn30.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn30.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [30]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [30]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [30]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [30]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &boost_and_two_level__1_sm_ehs_M->
          NonInlinedSFcns.Sfcn30.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn30.UPtrs0;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.DataTypeConversion_c;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn30.UPtrs1;
          sfcnUPtrs[0] = &boost_and_two_level__1_sm_ehs_B.Constant1_f;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn30.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.OpTrigger));
        }
      }

      /* path info */
      ssSetModelName(rts, "OpTrigger");
      ssSetPath(rts, "boost_and_two_level__1_sm_ehs/SM_eHS/OpTrigger");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn30.params;
        ssSetSFcnParamsCount(rts, 3);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpTrigger_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpTrigger_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpTrigger_P3_Size);
      }

      /* work vectors */
      ssSetIWork(rts, (int_T *)
                 &boost_and_two_level__1_sm_ehs_DW.OpTrigger_IWORK[0]);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn30.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn30.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* IWORK */
        ssSetDWorkWidth(rts, 0, 5);
        ssSetDWorkDataType(rts, 0,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &boost_and_two_level__1_sm_ehs_DW.OpTrigger_IWORK[0]);
      }

      /* registration */
      optrigger(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 0, 0);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
    }

    /* Level2 S-Function Block: boost_and_two_level__1_sm_ehs/<S2>/OpCtrl (sfun_ctrl_op7160ex1) */
    {
      SimStruct *rts = boost_and_two_level__1_sm_ehs_M->childSfunctions[31];

      /* timing info */
      time_T *sfcnPeriod =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn31.sfcnPeriod;
      time_T *sfcnOffset =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn31.sfcnOffset;
      int_T *sfcnTsMap =
        boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn31.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.blkInfo2
                         [31]);
      }

      ssSetRTWSfcnInfo(rts, boost_and_two_level__1_sm_ehs_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods2
                           [31]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.methods3
                           [31]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.statesInfo2
                         [31]);
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn31.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 2);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &boost_and_two_level__1_sm_ehs_B.OpCtrl_o1));
        }

        /* port 1 */
        {
          _ssSetOutputPortNumDimensions(rts, 1, 1);
          ssSetOutputPortWidth(rts, 1, 4);
          ssSetOutputPortSignal(rts, 1, ((real_T *)
            boost_and_two_level__1_sm_ehs_B.OpCtrl_o2));
        }
      }

      /* path info */
      ssSetModelName(rts, "OpCtrl");
      ssSetPath(rts, "boost_and_two_level__1_sm_ehs/SM_eHS/OpCtrl");
      ssSetRTModel(rts,boost_and_two_level__1_sm_ehs_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn31.params;
        ssSetSFcnParamsCount(rts, 12);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpCtrl_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpCtrl_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpCtrl_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpCtrl_P4_Size);
        ssSetSFcnParam(rts, 4, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpCtrl_P5_Size);
        ssSetSFcnParam(rts, 5, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpCtrl_P6_Size);
        ssSetSFcnParam(rts, 6, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpCtrl_P7_Size);
        ssSetSFcnParam(rts, 7, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpCtrl_P8_Size);
        ssSetSFcnParam(rts, 8, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpCtrl_P9_Size);
        ssSetSFcnParam(rts, 9, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpCtrl_P10_Size);
        ssSetSFcnParam(rts, 10, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpCtrl_P11_Size);
        ssSetSFcnParam(rts, 11, (mxArray*)
                       boost_and_two_level__1_sm_ehs_P.OpCtrl_P12_Size);
      }

      /* work vectors */
      ssSetIWork(rts, (int_T *) &boost_and_two_level__1_sm_ehs_DW.OpCtrl_IWORK);
      ssSetPWork(rts, (void **) &boost_and_two_level__1_sm_ehs_DW.OpCtrl_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn31.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &boost_and_two_level__1_sm_ehs_M->NonInlinedSFcns.Sfcn31.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 2);

        /* IWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &boost_and_two_level__1_sm_ehs_DW.OpCtrl_IWORK);

        /* PWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &boost_and_two_level__1_sm_ehs_DW.OpCtrl_PWORK);
      }

      /* registration */
      sfun_ctrl_op7160ex1(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetOutputPortConnected(rts, 0, 0);
      _ssSetOutputPortConnected(rts, 1, 0);
      _ssSetOutputPortBeingMerged(rts, 0, 0);
      _ssSetOutputPortBeingMerged(rts, 1, 0);

      /* Update the BufferDstPort flags for each input port */
    }
  }

  /* Initialize Sizes */
  boost_and_two_level__1_sm_ehs_M->Sizes.numContStates = (2);/* Number of continuous states */
  boost_and_two_level__1_sm_ehs_M->Sizes.numY = (0);/* Number of model outputs */
  boost_and_two_level__1_sm_ehs_M->Sizes.numU = (0);/* Number of model inputs */
  boost_and_two_level__1_sm_ehs_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  boost_and_two_level__1_sm_ehs_M->Sizes.numSampTimes = (2);/* Number of sample times */
  boost_and_two_level__1_sm_ehs_M->Sizes.numBlocks = (457);/* Number of blocks */
  boost_and_two_level__1_sm_ehs_M->Sizes.numBlockIO = (403);/* Number of block outputs */
  boost_and_two_level__1_sm_ehs_M->Sizes.numBlockPrms = (981);/* Sum of parameter "widths" */
  return boost_and_two_level__1_sm_ehs_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
