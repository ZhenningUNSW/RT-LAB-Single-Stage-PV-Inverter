/*
 * my_pv_system_1_sm_master.c
 *
 * Code generation for model "my_pv_system_1_sm_master".
 *
 * Model version              : 1.231
 * Simulink Coder version : 8.7 (R2014b) 08-Sep-2014
 * C source code generated on : Fri May 26 22:43:51 2017
 *
 * Target selection: rtlab_rtmodel.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "my_pv_system_1_sm_master.h"
#include "my_pv_system_1_sm_master_private.h"

const real_T my_pv_system_1_sm_master_RGND = 0.0;/* real_T ground */

/* Block signals (auto storage) */
B_my_pv_system_1_sm_master_T my_pv_system_1_sm_master_B;

/* Continuous states */
X_my_pv_system_1_sm_master_T my_pv_system_1_sm_master_X;

/* Block states (auto storage) */
DW_my_pv_system_1_sm_master_T my_pv_system_1_sm_master_DW;

/* Real-time model */
RT_MODEL_my_pv_system_1_sm_master_T my_pv_system_1_sm_master_M_;
RT_MODEL_my_pv_system_1_sm_master_T *const my_pv_system_1_sm_master_M =
  &my_pv_system_1_sm_master_M_;

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
  my_pv_system_1_sm_master_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  my_pv_system_1_sm_master_output();
  my_pv_system_1_sm_master_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  my_pv_system_1_sm_master_output();
  my_pv_system_1_sm_master_derivatives();

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

/* Model output function */
void my_pv_system_1_sm_master_output(void)
{
  real_T y;
  real_T u1;
  real_T u2;
  if (rtmIsMajorTimeStep(my_pv_system_1_sm_master_M)) {
    /* set solver stop time */
    if (!(my_pv_system_1_sm_master_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&my_pv_system_1_sm_master_M->solverInfo,
                            ((my_pv_system_1_sm_master_M->Timing.clockTickH0 + 1)
        * my_pv_system_1_sm_master_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&my_pv_system_1_sm_master_M->solverInfo,
                            ((my_pv_system_1_sm_master_M->Timing.clockTick0 + 1)
        * my_pv_system_1_sm_master_M->Timing.stepSize0 +
        my_pv_system_1_sm_master_M->Timing.clockTickH0 *
        my_pv_system_1_sm_master_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(my_pv_system_1_sm_master_M)) {
    my_pv_system_1_sm_master_M->Timing.t[0] = rtsiGetT
      (&my_pv_system_1_sm_master_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(my_pv_system_1_sm_master_M)) {
    /* Memory: '<S2>/S-Function' */
    my_pv_system_1_sm_master_B.SFunction =
      my_pv_system_1_sm_master_DW.SFunction_PreviousInput;

    /* Sum: '<S2>/Sum' incorporates:
     *  Constant: '<S2>/S-Function1'
     */
    my_pv_system_1_sm_master_B.Sum = my_pv_system_1_sm_master_P.SFunction1_Value
      + my_pv_system_1_sm_master_B.SFunction;

    /* Stop: '<S2>/Stop Simulation' */
    if (my_pv_system_1_sm_master_B.Sum != 0.0) {
      rtmSetStopRequested(my_pv_system_1_sm_master_M, 1);
    }

    /* End of Stop: '<S2>/Stop Simulation' */

    /* Memory: '<S3>/Memory' */
    my_pv_system_1_sm_master_B.IrradianceWm2 =
      my_pv_system_1_sm_master_DW.Memory_1_PreviousInput;

    /* Memory: '<S3>/Memory' */
    my_pv_system_1_sm_master_B.I_PV =
      my_pv_system_1_sm_master_DW.Memory_2_PreviousInput;

    /* Memory: '<S3>/Memory' */
    my_pv_system_1_sm_master_B.V_PV =
      my_pv_system_1_sm_master_DW.Memory_3_PreviousInput;

    /* Memory: '<S3>/Memory' */
    my_pv_system_1_sm_master_B.I_Diode =
      my_pv_system_1_sm_master_DW.Memory_4_PreviousInput;

    /* Memory: '<S3>/Memory' */
    my_pv_system_1_sm_master_B.Temperature =
      my_pv_system_1_sm_master_DW.Memory_5_PreviousInput;

    /* Memory: '<S3>/Memory' */
    my_pv_system_1_sm_master_B.Vph =
      my_pv_system_1_sm_master_DW.Memory_6_PreviousInput;

    /* Memory: '<S3>/Memory1' */
    my_pv_system_1_sm_master_B.Vinv =
      my_pv_system_1_sm_master_DW.Memory1_1_PreviousInput;

    /* Memory: '<S3>/Memory1' */
    my_pv_system_1_sm_master_B.Vdc =
      my_pv_system_1_sm_master_DW.Memory1_2_PreviousInput;

    /* Memory: '<S3>/Memory1' */
    my_pv_system_1_sm_master_B.Igrid =
      my_pv_system_1_sm_master_DW.Memory1_3_PreviousInput;

    /* Memory: '<S3>/Memory1' */
    my_pv_system_1_sm_master_B.Vgrid =
      my_pv_system_1_sm_master_DW.Memory1_4_PreviousInput;

    /* Memory: '<S3>/Memory1' */
    my_pv_system_1_sm_master_B.Ig =
      my_pv_system_1_sm_master_DW.Memory1_5_PreviousInput;

    /* Memory: '<S3>/Memory1' */
    my_pv_system_1_sm_master_B.P_PV =
      my_pv_system_1_sm_master_DW.Memory1_6_PreviousInput;

    /* Memory: '<S3>/Memory1' */
    my_pv_system_1_sm_master_B.P_PV_MEAN =
      my_pv_system_1_sm_master_DW.Memory1_7_PreviousInput;

    /* Memory: '<S3>/Memory1' */
    my_pv_system_1_sm_master_B.Gate[0] =
      my_pv_system_1_sm_master_DW.Memory1_8_PreviousInput[0];
    my_pv_system_1_sm_master_B.Gate[1] =
      my_pv_system_1_sm_master_DW.Memory1_8_PreviousInput[1];
    my_pv_system_1_sm_master_B.Gate[2] =
      my_pv_system_1_sm_master_DW.Memory1_8_PreviousInput[2];
    my_pv_system_1_sm_master_B.Gate[3] =
      my_pv_system_1_sm_master_DW.Memory1_8_PreviousInput[3];

    /* Outputs for Atomic SubSystem: '<S16>/Subsystem5' */

    /* Level2 S-Function Block: '<S50>/S-Function' (send_rt) */
    {
      SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[0];
      sfcnOutputs(rts, 1);
    }

    /* End of Outputs for SubSystem: '<S16>/Subsystem5' */

    /* Level2 S-Function Block: '<S48>/S-Function' (OP_SEND) */
    {
      SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[1];
      sfcnOutputs(rts, 1);
    }

    /* Level2 S-Function Block: '<S3>/OpMonitor' (opmonitor) */
    {
      SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[2];
      sfcnOutputs(rts, 1);
    }

    /* Level2 S-Function Block: '<S49>/S-Function' (OP_SEND) */
    {
      SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[3];
      sfcnOutputs(rts, 1);
    }

    /* UnitDelay: '<S24>/Unit Delay' */
    my_pv_system_1_sm_master_B.UnitDelay[0] =
      my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[0];
    my_pv_system_1_sm_master_B.UnitDelay[1] =
      my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[1];
    my_pv_system_1_sm_master_B.UnitDelay[2] =
      my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[2];
    my_pv_system_1_sm_master_B.UnitDelay[3] =
      my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[3];

    /* Level2 S-Function Block: '<S51>/S-Function' (recv_rt) */
    {
      SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[4];
      sfcnOutputs(rts, 1);
    }

    /* DataTypeConversion: '<S27>/Data Type Conversion' */
    my_pv_system_1_sm_master_B.DataTypeConversion[0] =
      my_pv_system_1_sm_master_B.SFunction_c[0];
    my_pv_system_1_sm_master_B.DataTypeConversion[1] =
      my_pv_system_1_sm_master_B.SFunction_c[1];
    my_pv_system_1_sm_master_B.DataTypeConversion[2] =
      my_pv_system_1_sm_master_B.SFunction_c[2];
    my_pv_system_1_sm_master_B.DataTypeConversion[3] =
      my_pv_system_1_sm_master_B.SFunction_c[3];

    /* Outputs for Enabled SubSystem: '<S24>/Tail' incorporates:
     *  EnablePort: '<S25>/Enable'
     */
    if (rtmIsMajorTimeStep(my_pv_system_1_sm_master_M)) {
      /* Constant: '<S24>/Not in ARTEMIS' */
      if (my_pv_system_1_sm_master_P.NotinARTEMIS_Value > 0.0) {
        if (!my_pv_system_1_sm_master_DW.Tail_MODE) {
          /* InitializeConditions for DiscreteIntegrator: '<S25>/Discrete-Time Integrator' */
          my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[0] =
            my_pv_system_1_sm_master_P.DiscreteTimeIntegrator_IC;
          my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[1] =
            my_pv_system_1_sm_master_P.DiscreteTimeIntegrator_IC;
          my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[2] =
            my_pv_system_1_sm_master_P.DiscreteTimeIntegrator_IC;
          my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[3] =
            my_pv_system_1_sm_master_P.DiscreteTimeIntegrator_IC;
          my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[0] = 2;
          my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[1] = 2;
          my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[2] = 2;
          my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[3] = 2;

          /* InitializeConditions for UnitDelay: '<S25>/Unit Delay' */
          my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[0] =
            my_pv_system_1_sm_master_P.UnitDelay_InitialCondition;
          my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[1] =
            my_pv_system_1_sm_master_P.UnitDelay_InitialCondition;
          my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[2] =
            my_pv_system_1_sm_master_P.UnitDelay_InitialCondition;
          my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[3] =
            my_pv_system_1_sm_master_P.UnitDelay_InitialCondition;
          my_pv_system_1_sm_master_DW.Tail_MODE = true;
        }
      } else {
        if (my_pv_system_1_sm_master_DW.Tail_MODE) {
          my_pv_system_1_sm_master_DW.Tail_MODE = false;
        }
      }

      /* End of Constant: '<S24>/Not in ARTEMIS' */
    }

    if (my_pv_system_1_sm_master_DW.Tail_MODE) {
      /* DiscreteIntegrator: '<S25>/Discrete-Time Integrator' */
      if ((my_pv_system_1_sm_master_B.DataTypeConversion[0] <= 0.0) &&
          (my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[0] == 1))
      {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[0] =
          my_pv_system_1_sm_master_P.DiscreteTimeIntegrator_IC;
      }

      if ((my_pv_system_1_sm_master_B.DataTypeConversion[1] <= 0.0) &&
          (my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[1] == 1))
      {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[1] =
          my_pv_system_1_sm_master_P.DiscreteTimeIntegrator_IC;
      }

      if ((my_pv_system_1_sm_master_B.DataTypeConversion[2] <= 0.0) &&
          (my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[2] == 1))
      {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[2] =
          my_pv_system_1_sm_master_P.DiscreteTimeIntegrator_IC;
      }

      if ((my_pv_system_1_sm_master_B.DataTypeConversion[3] <= 0.0) &&
          (my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[3] == 1))
      {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[3] =
          my_pv_system_1_sm_master_P.DiscreteTimeIntegrator_IC;
      }

      my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[0] =
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[0];
      my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[1] =
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[1];
      my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[2] =
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[2];
      my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[3] =
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[3];

      /* End of DiscreteIntegrator: '<S25>/Discrete-Time Integrator' */

      /* Gain: '<S28>/-0.9//Tf' */
      y = my_pv_system_1_sm_master_P.Tail_Tf + 2.2204460492503131E-16;
      y = -0.9 / y;
      my_pv_system_1_sm_master_B.u9Tf[0] = y *
        my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[0];
      my_pv_system_1_sm_master_B.u9Tf[1] = y *
        my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[1];
      my_pv_system_1_sm_master_B.u9Tf[2] = y *
        my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[2];
      my_pv_system_1_sm_master_B.u9Tf[3] = y *
        my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[3];

      /* Sum: '<S28>/Add' incorporates:
       *  Constant: '<S28>/Constant'
       */
      my_pv_system_1_sm_master_B.Add[0] =
        my_pv_system_1_sm_master_P.Constant_Value +
        my_pv_system_1_sm_master_B.u9Tf[0];
      my_pv_system_1_sm_master_B.Add[1] =
        my_pv_system_1_sm_master_P.Constant_Value +
        my_pv_system_1_sm_master_B.u9Tf[1];
      my_pv_system_1_sm_master_B.Add[2] =
        my_pv_system_1_sm_master_P.Constant_Value +
        my_pv_system_1_sm_master_B.u9Tf[2];
      my_pv_system_1_sm_master_B.Add[3] =
        my_pv_system_1_sm_master_P.Constant_Value +
        my_pv_system_1_sm_master_B.u9Tf[3];

      /* Saturate: '<S28>/Saturation1' */
      y = my_pv_system_1_sm_master_B.Add[0];
      u1 = my_pv_system_1_sm_master_P.Saturation1_LowerSat;
      u2 = my_pv_system_1_sm_master_P.Saturation1_UpperSat;
      if (y > u2) {
        y = u2;
      } else {
        if (y < u1) {
          y = u1;
        }
      }

      my_pv_system_1_sm_master_B.Saturation1[0] = y;
      y = my_pv_system_1_sm_master_B.Add[1];
      u1 = my_pv_system_1_sm_master_P.Saturation1_LowerSat;
      u2 = my_pv_system_1_sm_master_P.Saturation1_UpperSat;
      if (y > u2) {
        y = u2;
      } else {
        if (y < u1) {
          y = u1;
        }
      }

      my_pv_system_1_sm_master_B.Saturation1[1] = y;
      y = my_pv_system_1_sm_master_B.Add[2];
      u1 = my_pv_system_1_sm_master_P.Saturation1_LowerSat;
      u2 = my_pv_system_1_sm_master_P.Saturation1_UpperSat;
      if (y > u2) {
        y = u2;
      } else {
        if (y < u1) {
          y = u1;
        }
      }

      my_pv_system_1_sm_master_B.Saturation1[2] = y;
      y = my_pv_system_1_sm_master_B.Add[3];
      u1 = my_pv_system_1_sm_master_P.Saturation1_LowerSat;
      u2 = my_pv_system_1_sm_master_P.Saturation1_UpperSat;
      if (y > u2) {
        y = u2;
      } else {
        if (y < u1) {
          y = u1;
        }
      }

      my_pv_system_1_sm_master_B.Saturation1[3] = y;

      /* End of Saturate: '<S28>/Saturation1' */

      /* Sum: '<S28>/Add1' incorporates:
       *  Constant: '<S28>/Constant2'
       */
      y = my_pv_system_1_sm_master_P.Tail_Tf +
        my_pv_system_1_sm_master_P.Tail_Tt;
      my_pv_system_1_sm_master_B.Add1[0] = y -
        my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[0];
      my_pv_system_1_sm_master_B.Add1[1] = y -
        my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[1];
      my_pv_system_1_sm_master_B.Add1[2] = y -
        my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[2];
      my_pv_system_1_sm_master_B.Add1[3] = y -
        my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[3];

      /* Gain: '<S28>/0.1//Tt' */
      y = my_pv_system_1_sm_master_P.Tail_Tt + 2.2204460492503131E-16;
      y = 0.1 / y;
      my_pv_system_1_sm_master_B.uTt[0] = y * my_pv_system_1_sm_master_B.Add1[0];
      my_pv_system_1_sm_master_B.uTt[1] = y * my_pv_system_1_sm_master_B.Add1[1];
      my_pv_system_1_sm_master_B.uTt[2] = y * my_pv_system_1_sm_master_B.Add1[2];
      my_pv_system_1_sm_master_B.uTt[3] = y * my_pv_system_1_sm_master_B.Add1[3];

      /* Saturate: '<S28>/Saturation2' */
      y = my_pv_system_1_sm_master_B.uTt[0];
      u1 = my_pv_system_1_sm_master_P.Saturation2_LowerSat;
      u2 = my_pv_system_1_sm_master_P.Saturation2_UpperSat;
      if (y > u2) {
        y = u2;
      } else {
        if (y < u1) {
          y = u1;
        }
      }

      my_pv_system_1_sm_master_B.Saturation2[0] = y;
      y = my_pv_system_1_sm_master_B.uTt[1];
      u1 = my_pv_system_1_sm_master_P.Saturation2_LowerSat;
      u2 = my_pv_system_1_sm_master_P.Saturation2_UpperSat;
      if (y > u2) {
        y = u2;
      } else {
        if (y < u1) {
          y = u1;
        }
      }

      my_pv_system_1_sm_master_B.Saturation2[1] = y;
      y = my_pv_system_1_sm_master_B.uTt[2];
      u1 = my_pv_system_1_sm_master_P.Saturation2_LowerSat;
      u2 = my_pv_system_1_sm_master_P.Saturation2_UpperSat;
      if (y > u2) {
        y = u2;
      } else {
        if (y < u1) {
          y = u1;
        }
      }

      my_pv_system_1_sm_master_B.Saturation2[2] = y;
      y = my_pv_system_1_sm_master_B.uTt[3];
      u1 = my_pv_system_1_sm_master_P.Saturation2_LowerSat;
      u2 = my_pv_system_1_sm_master_P.Saturation2_UpperSat;
      if (y > u2) {
        y = u2;
      } else {
        if (y < u1) {
          y = u1;
        }
      }

      my_pv_system_1_sm_master_B.Saturation2[3] = y;

      /* End of Saturate: '<S28>/Saturation2' */

      /* Sum: '<S28>/Add2' */
      my_pv_system_1_sm_master_B.Add2[0] =
        my_pv_system_1_sm_master_B.Saturation1[0] +
        my_pv_system_1_sm_master_B.Saturation2[0];
      my_pv_system_1_sm_master_B.Add2[1] =
        my_pv_system_1_sm_master_B.Saturation1[1] +
        my_pv_system_1_sm_master_B.Saturation2[1];
      my_pv_system_1_sm_master_B.Add2[2] =
        my_pv_system_1_sm_master_B.Saturation1[2] +
        my_pv_system_1_sm_master_B.Saturation2[2];
      my_pv_system_1_sm_master_B.Add2[3] =
        my_pv_system_1_sm_master_B.Saturation1[3] +
        my_pv_system_1_sm_master_B.Saturation2[3];

      /* UnitDelay: '<S25>/Unit Delay' */
      my_pv_system_1_sm_master_B.UnitDelay_l[0] =
        my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[0];
      my_pv_system_1_sm_master_B.UnitDelay_l[1] =
        my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[1];
      my_pv_system_1_sm_master_B.UnitDelay_l[2] =
        my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[2];
      my_pv_system_1_sm_master_B.UnitDelay_l[3] =
        my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[3];

      /* Switch: '<S25>/Switch' */
      if (my_pv_system_1_sm_master_B.DataTypeConversion[0] >=
          my_pv_system_1_sm_master_P.Switch_Threshold) {
        my_pv_system_1_sm_master_B.Switch_j[0] =
          my_pv_system_1_sm_master_B.UnitDelay[0];
      } else {
        my_pv_system_1_sm_master_B.Switch_j[0] =
          my_pv_system_1_sm_master_B.UnitDelay_l[0];
      }

      if (my_pv_system_1_sm_master_B.DataTypeConversion[1] >=
          my_pv_system_1_sm_master_P.Switch_Threshold) {
        my_pv_system_1_sm_master_B.Switch_j[1] =
          my_pv_system_1_sm_master_B.UnitDelay[1];
      } else {
        my_pv_system_1_sm_master_B.Switch_j[1] =
          my_pv_system_1_sm_master_B.UnitDelay_l[1];
      }

      if (my_pv_system_1_sm_master_B.DataTypeConversion[2] >=
          my_pv_system_1_sm_master_P.Switch_Threshold) {
        my_pv_system_1_sm_master_B.Switch_j[2] =
          my_pv_system_1_sm_master_B.UnitDelay[2];
      } else {
        my_pv_system_1_sm_master_B.Switch_j[2] =
          my_pv_system_1_sm_master_B.UnitDelay_l[2];
      }

      if (my_pv_system_1_sm_master_B.DataTypeConversion[3] >=
          my_pv_system_1_sm_master_P.Switch_Threshold) {
        my_pv_system_1_sm_master_B.Switch_j[3] =
          my_pv_system_1_sm_master_B.UnitDelay[3];
      } else {
        my_pv_system_1_sm_master_B.Switch_j[3] =
          my_pv_system_1_sm_master_B.UnitDelay_l[3];
      }

      /* End of Switch: '<S25>/Switch' */

      /* Product: '<S25>/Product' incorporates:
       *  Constant: '<S25>/2'
       */
      my_pv_system_1_sm_master_B.Product_d[0] = my_pv_system_1_sm_master_B.Add2
        [0] * my_pv_system_1_sm_master_B.Switch_j[0] *
        my_pv_system_1_sm_master_P._Value_a;
      my_pv_system_1_sm_master_B.Product_d[1] = my_pv_system_1_sm_master_B.Add2
        [1] * my_pv_system_1_sm_master_B.Switch_j[1] *
        my_pv_system_1_sm_master_P._Value_a;
      my_pv_system_1_sm_master_B.Product_d[2] = my_pv_system_1_sm_master_B.Add2
        [2] * my_pv_system_1_sm_master_B.Switch_j[2] *
        my_pv_system_1_sm_master_P._Value_a;
      my_pv_system_1_sm_master_B.Product_d[3] = my_pv_system_1_sm_master_B.Add2
        [3] * my_pv_system_1_sm_master_B.Switch_j[3] *
        my_pv_system_1_sm_master_P._Value_a;
    }

    /* End of Outputs for SubSystem: '<S24>/Tail' */
  }

  /* Sin: '<S44>/AC' */
  u1 = my_pv_system_1_sm_master_P.Vs14400V_Phase * 3.1415926535897931;
  y = u1 / 180.0;
  u1 = 6.2831853071795862 * my_pv_system_1_sm_master_P.Vs14400V_Frequency;
  my_pv_system_1_sm_master_B.AC = sin(u1 * my_pv_system_1_sm_master_M->Timing.t
    [0] + y) * my_pv_system_1_sm_master_P.Vs14400V_Amplitude +
    my_pv_system_1_sm_master_P.AC_Bias;
  if (rtmIsMajorTimeStep(my_pv_system_1_sm_master_M)) {
    /* Level2 S-Function Block: '<S54>/State-Space' (fts5abcd_noncomp) */
    {
      SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[5];
      sfcnOutputs(rts, 1);
    }

    /* Gain: '<S24>/1//Ron' */
    my_pv_system_1_sm_master_B.Ron[0] = my_pv_system_1_sm_master_P.Ron_Gain *
      my_pv_system_1_sm_master_B.StateSpace_o1[0];
    my_pv_system_1_sm_master_B.Ron[1] = my_pv_system_1_sm_master_P.Ron_Gain *
      my_pv_system_1_sm_master_B.StateSpace_o1[1];
    my_pv_system_1_sm_master_B.Ron[2] = my_pv_system_1_sm_master_P.Ron_Gain *
      my_pv_system_1_sm_master_B.StateSpace_o1[2];
    my_pv_system_1_sm_master_B.Ron[3] = my_pv_system_1_sm_master_P.Ron_Gain *
      my_pv_system_1_sm_master_B.StateSpace_o1[3];

    /* Switch: '<S24>/Switch' incorporates:
     *  Constant: '<S24>/0 4'
     */
    if (my_pv_system_1_sm_master_B.StateSpace_o2[0] >=
        my_pv_system_1_sm_master_P.Switch_Threshold_l) {
      my_pv_system_1_sm_master_B.Switch[0] = my_pv_system_1_sm_master_B.Ron[0];
    } else {
      my_pv_system_1_sm_master_B.Switch[0] = my_pv_system_1_sm_master_P.u_Value;
    }

    if (my_pv_system_1_sm_master_B.StateSpace_o2[1] >=
        my_pv_system_1_sm_master_P.Switch_Threshold_l) {
      my_pv_system_1_sm_master_B.Switch[1] = my_pv_system_1_sm_master_B.Ron[1];
    } else {
      my_pv_system_1_sm_master_B.Switch[1] = my_pv_system_1_sm_master_P.u_Value;
    }

    if (my_pv_system_1_sm_master_B.StateSpace_o2[2] >=
        my_pv_system_1_sm_master_P.Switch_Threshold_l) {
      my_pv_system_1_sm_master_B.Switch[2] = my_pv_system_1_sm_master_B.Ron[2];
    } else {
      my_pv_system_1_sm_master_B.Switch[2] = my_pv_system_1_sm_master_P.u_Value;
    }

    if (my_pv_system_1_sm_master_B.StateSpace_o2[3] >=
        my_pv_system_1_sm_master_P.Switch_Threshold_l) {
      my_pv_system_1_sm_master_B.Switch[3] = my_pv_system_1_sm_master_B.Ron[3];
    } else {
      my_pv_system_1_sm_master_B.Switch[3] = my_pv_system_1_sm_master_P.u_Value;
    }

    /* End of Switch: '<S24>/Switch' */

    /* Saturate: '<S24>/Saturation' */
    y = my_pv_system_1_sm_master_B.Switch[0];
    u1 = my_pv_system_1_sm_master_P.Saturation_LowerSat;
    u2 = my_pv_system_1_sm_master_P.Saturation_UpperSat;
    if (y > u2) {
      y = u2;
    } else {
      if (y < u1) {
        y = u1;
      }
    }

    my_pv_system_1_sm_master_B.Saturation[0] = y;
    y = my_pv_system_1_sm_master_B.Switch[1];
    u1 = my_pv_system_1_sm_master_P.Saturation_LowerSat;
    u2 = my_pv_system_1_sm_master_P.Saturation_UpperSat;
    if (y > u2) {
      y = u2;
    } else {
      if (y < u1) {
        y = u1;
      }
    }

    my_pv_system_1_sm_master_B.Saturation[1] = y;
    y = my_pv_system_1_sm_master_B.Switch[2];
    u1 = my_pv_system_1_sm_master_P.Saturation_LowerSat;
    u2 = my_pv_system_1_sm_master_P.Saturation_UpperSat;
    if (y > u2) {
      y = u2;
    } else {
      if (y < u1) {
        y = u1;
      }
    }

    my_pv_system_1_sm_master_B.Saturation[2] = y;
    y = my_pv_system_1_sm_master_B.Switch[3];
    u1 = my_pv_system_1_sm_master_P.Saturation_LowerSat;
    u2 = my_pv_system_1_sm_master_P.Saturation_UpperSat;
    if (y > u2) {
      y = u2;
    } else {
      if (y < u1) {
        y = u1;
      }
    }

    my_pv_system_1_sm_master_B.Saturation[3] = y;

    /* End of Saturate: '<S24>/Saturation' */

    /* Switch: '<S26>/Switch' incorporates:
     *  Constant: '<S26>/Vf Devices & Clamping Diodes'
     *  Constant: '<S26>/Vf Diodes'
     */
    if (my_pv_system_1_sm_master_B.UnitDelay[0] != 0.0) {
      my_pv_system_1_sm_master_B.Switch_e[0] =
        my_pv_system_1_sm_master_P.VfDevicesClampingDiodes_Value[0];
    } else {
      my_pv_system_1_sm_master_B.Switch_e[0] =
        my_pv_system_1_sm_master_P.VfDiodes_Value[0];
    }

    if (my_pv_system_1_sm_master_B.UnitDelay[1] != 0.0) {
      my_pv_system_1_sm_master_B.Switch_e[1] =
        my_pv_system_1_sm_master_P.VfDevicesClampingDiodes_Value[1];
    } else {
      my_pv_system_1_sm_master_B.Switch_e[1] =
        my_pv_system_1_sm_master_P.VfDiodes_Value[1];
    }

    if (my_pv_system_1_sm_master_B.UnitDelay[2] != 0.0) {
      my_pv_system_1_sm_master_B.Switch_e[2] =
        my_pv_system_1_sm_master_P.VfDevicesClampingDiodes_Value[2];
    } else {
      my_pv_system_1_sm_master_B.Switch_e[2] =
        my_pv_system_1_sm_master_P.VfDiodes_Value[2];
    }

    if (my_pv_system_1_sm_master_B.UnitDelay[3] != 0.0) {
      my_pv_system_1_sm_master_B.Switch_e[3] =
        my_pv_system_1_sm_master_P.VfDevicesClampingDiodes_Value[3];
    } else {
      my_pv_system_1_sm_master_B.Switch_e[3] =
        my_pv_system_1_sm_master_P.VfDiodes_Value[3];
    }

    /* End of Switch: '<S26>/Switch' */

    /* Gain: '<S6>/do not delete this gain' */
    my_pv_system_1_sm_master_B.donotdeletethisgain =
      my_pv_system_1_sm_master_P.donotdeletethisgain_Gain *
      my_pv_system_1_sm_master_B.StateSpace_o1[9];
  }

  /* Clock: '<S29>/Clock' */
  my_pv_system_1_sm_master_B.Clock = my_pv_system_1_sm_master_M->Timing.t[0];

  /* Integrator: '<S29>/integrator' */
  my_pv_system_1_sm_master_B.integrator =
    my_pv_system_1_sm_master_X.integrator_CSTATE;

  /* TransportDelay: '<S29>/Transport Delay' */
  {
    real_T **uBuffer = (real_T**)
      &my_pv_system_1_sm_master_DW.TransportDelay_PWORK.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)
      &my_pv_system_1_sm_master_DW.TransportDelay_PWORK.TUbufferPtrs[1];
    real_T simTime = my_pv_system_1_sm_master_M->Timing.t[0];
    real_T tMinusDelay = simTime -
      (my_pv_system_1_sm_master_P.TransportDelay_Delay);
    my_pv_system_1_sm_master_B.TransportDelay = rt_TDelayInterpolate(
      tMinusDelay,
      0.0,
      *tBuffer,
      *uBuffer,
      my_pv_system_1_sm_master_DW.TransportDelay_IWORK.CircularBufSize,
      &my_pv_system_1_sm_master_DW.TransportDelay_IWORK.Last,
      my_pv_system_1_sm_master_DW.TransportDelay_IWORK.Tail,
      my_pv_system_1_sm_master_DW.TransportDelay_IWORK.Head,
      my_pv_system_1_sm_master_P.TransportDelay_InitOutput,
      0,
      0);
  }

  if (rtmIsMajorTimeStep(my_pv_system_1_sm_master_M)) {
    /* Memory: '<S29>/Memory' */
    my_pv_system_1_sm_master_B.Memory =
      my_pv_system_1_sm_master_DW.Memory_PreviousInput;
  }

  /* RelationalOperator: '<S29>/Relational Operator' incorporates:
   *  Constant: '<S29>/K1'
   */
  my_pv_system_1_sm_master_B.RelationalOperator =
    (my_pv_system_1_sm_master_B.Clock >= my_pv_system_1_sm_master_P.K1_Value);

  /* Switch: '<S29>/Switch' */
  if (my_pv_system_1_sm_master_B.RelationalOperator) {
    /* Sum: '<S29>/Sum' */
    my_pv_system_1_sm_master_B.Sum_g = my_pv_system_1_sm_master_B.integrator -
      my_pv_system_1_sm_master_B.TransportDelay;

    /* Gain: '<S29>/Gain' */
    my_pv_system_1_sm_master_B.Gain_f = my_pv_system_1_sm_master_P.Gain_Gain *
      my_pv_system_1_sm_master_B.Sum_g;
    my_pv_system_1_sm_master_B.Switch_d = my_pv_system_1_sm_master_B.Gain_f;
  } else {
    my_pv_system_1_sm_master_B.Switch_d = my_pv_system_1_sm_master_B.Memory;
  }

  /* End of Switch: '<S29>/Switch' */

  /* Clock: '<S30>/Clock' */
  my_pv_system_1_sm_master_B.Clock_c = my_pv_system_1_sm_master_M->Timing.t[0];

  /* Integrator: '<S30>/integrator' */
  my_pv_system_1_sm_master_B.integrator_i =
    my_pv_system_1_sm_master_X.integrator_CSTATE_a;

  /* TransportDelay: '<S30>/Transport Delay' */
  {
    real_T **uBuffer = (real_T**)
      &my_pv_system_1_sm_master_DW.TransportDelay_PWORK_a.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)
      &my_pv_system_1_sm_master_DW.TransportDelay_PWORK_a.TUbufferPtrs[1];
    real_T simTime = my_pv_system_1_sm_master_M->Timing.t[0];
    real_T tMinusDelay = simTime -
      (my_pv_system_1_sm_master_P.TransportDelay_Delay_j);
    my_pv_system_1_sm_master_B.TransportDelay_j = rt_TDelayInterpolate(
      tMinusDelay,
      0.0,
      *tBuffer,
      *uBuffer,
      my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.CircularBufSize,
      &my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.Last,
      my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.Tail,
      my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.Head,
      my_pv_system_1_sm_master_P.TransportDelay_InitOutput_i,
      0,
      0);
  }

  if (rtmIsMajorTimeStep(my_pv_system_1_sm_master_M)) {
    /* Memory: '<S30>/Memory' */
    my_pv_system_1_sm_master_B.Memory_j =
      my_pv_system_1_sm_master_DW.Memory_PreviousInput_h;
  }

  /* RelationalOperator: '<S30>/Relational Operator' incorporates:
   *  Constant: '<S30>/K1'
   */
  my_pv_system_1_sm_master_B.RelationalOperator_h =
    (my_pv_system_1_sm_master_B.Clock_c >= my_pv_system_1_sm_master_P.K1_Value_o);

  /* Switch: '<S30>/Switch' */
  if (my_pv_system_1_sm_master_B.RelationalOperator_h) {
    /* Sum: '<S30>/Sum' */
    my_pv_system_1_sm_master_B.Sum_p = my_pv_system_1_sm_master_B.integrator_i -
      my_pv_system_1_sm_master_B.TransportDelay_j;

    /* Gain: '<S30>/Gain' */
    my_pv_system_1_sm_master_B.Gain = my_pv_system_1_sm_master_P.Gain_Gain_k *
      my_pv_system_1_sm_master_B.Sum_p;
    my_pv_system_1_sm_master_B.Switch_p = my_pv_system_1_sm_master_B.Gain;
  } else {
    my_pv_system_1_sm_master_B.Switch_p = my_pv_system_1_sm_master_B.Memory_j;
  }

  /* End of Switch: '<S30>/Switch' */
  if (rtmIsMajorTimeStep(my_pv_system_1_sm_master_M)) {
    /* Level2 S-Function Block: '<S31>/S-Function' (RECV_Param) */
    {
      SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[6];
      sfcnOutputs(rts, 1);
    }
  }

  /* Product: '<S3>/Product' */
  my_pv_system_1_sm_master_B.Product = my_pv_system_1_sm_master_B.Switch_d *
    my_pv_system_1_sm_master_B.Switch_p;
  if (rtmIsMajorTimeStep(my_pv_system_1_sm_master_M)) {
    /* Gain: '<S35>/do not delete this gain' */
    my_pv_system_1_sm_master_B.donotdeletethisgain_d =
      my_pv_system_1_sm_master_P.donotdeletethisgain_Gain_n *
      my_pv_system_1_sm_master_B.StateSpace_o1[5];

    /* Gain: '<S34>/do not delete this gain' */
    my_pv_system_1_sm_master_B.donotdeletethisgain_a =
      my_pv_system_1_sm_master_P.donotdeletethisgain_Gain_j *
      my_pv_system_1_sm_master_B.StateSpace_o1[11];

    /* Product: '<S3>/Product1' */
    my_pv_system_1_sm_master_B.Product1 =
      my_pv_system_1_sm_master_B.donotdeletethisgain_d *
      my_pv_system_1_sm_master_B.donotdeletethisgain_a;

    /* RateLimiter: '<S12>/Rate Limiter' */
    u1 = my_pv_system_1_sm_master_B.SFunction_f[0] -
      my_pv_system_1_sm_master_DW.PrevY;
    if (u1 > my_pv_system_1_sm_master_P.RateLimiter_RisingLim) {
      my_pv_system_1_sm_master_B.IrradianceWm2_j =
        my_pv_system_1_sm_master_DW.PrevY +
        my_pv_system_1_sm_master_P.RateLimiter_RisingLim;
    } else if (u1 < my_pv_system_1_sm_master_P.RateLimiter_FallingLim) {
      my_pv_system_1_sm_master_B.IrradianceWm2_j =
        my_pv_system_1_sm_master_DW.PrevY +
        my_pv_system_1_sm_master_P.RateLimiter_FallingLim;
    } else {
      my_pv_system_1_sm_master_B.IrradianceWm2_j =
        my_pv_system_1_sm_master_B.SFunction_f[0];
    }

    my_pv_system_1_sm_master_DW.PrevY =
      my_pv_system_1_sm_master_B.IrradianceWm2_j;

    /* End of RateLimiter: '<S12>/Rate Limiter' */

    /* Saturate: '<S12>/Saturation' */
    y = my_pv_system_1_sm_master_B.SFunction_f[1];
    u1 = my_pv_system_1_sm_master_P.Saturation_LowerSat_o;
    u2 = my_pv_system_1_sm_master_P.Saturation_UpperSat_e;
    if (y > u2) {
      my_pv_system_1_sm_master_B.Temperature_l = u2;
    } else if (y < u1) {
      my_pv_system_1_sm_master_B.Temperature_l = u1;
    } else {
      my_pv_system_1_sm_master_B.Temperature_l = y;
    }

    /* End of Saturate: '<S12>/Saturation' */

    /* Gain: '<S36>/do not delete this gain' */
    my_pv_system_1_sm_master_B.donotdeletethisgain_db =
      my_pv_system_1_sm_master_P.donotdeletethisgain_Gain_jy *
      my_pv_system_1_sm_master_B.StateSpace_o1[6];

    /* Gain: '<S13>/do not delete this gain' */
    my_pv_system_1_sm_master_B.donotdeletethisgain_e =
      my_pv_system_1_sm_master_P.donotdeletethisgain_Gain_g *
      my_pv_system_1_sm_master_B.StateSpace_o1[7];

    /* Gain: '<S14>/do not delete this gain' */
    my_pv_system_1_sm_master_B.donotdeletethisgain_g =
      my_pv_system_1_sm_master_P.donotdeletethisgain_Gain_h *
      my_pv_system_1_sm_master_B.StateSpace_o1[8];

    /* Gain: '<S20>/do not delete this gain' */
    my_pv_system_1_sm_master_B.donotdeletethisgain_i =
      my_pv_system_1_sm_master_P.donotdeletethisgain_Gain_c *
      my_pv_system_1_sm_master_B.StateSpace_o1[10];

    /* Gain: '<S21>/do not delete this gain' */
    my_pv_system_1_sm_master_B.donotdeletethisgain_de =
      my_pv_system_1_sm_master_P.donotdeletethisgain_Gain_k *
      my_pv_system_1_sm_master_B.StateSpace_o1[4];
  }
}

/* Model update function */
void my_pv_system_1_sm_master_update(void)
{
  real_T tmp;
  if (rtmIsMajorTimeStep(my_pv_system_1_sm_master_M)) {
    /* Update for Memory: '<S2>/S-Function' */
    my_pv_system_1_sm_master_DW.SFunction_PreviousInput =
      my_pv_system_1_sm_master_B.Sum;

    /* Update for Memory: '<S3>/Memory' */
    my_pv_system_1_sm_master_DW.Memory_1_PreviousInput =
      my_pv_system_1_sm_master_B.IrradianceWm2_j;

    /* Update for Memory: '<S3>/Memory' */
    my_pv_system_1_sm_master_DW.Memory_2_PreviousInput =
      my_pv_system_1_sm_master_B.donotdeletethisgain_a;

    /* Update for Memory: '<S3>/Memory' */
    my_pv_system_1_sm_master_DW.Memory_3_PreviousInput =
      my_pv_system_1_sm_master_B.donotdeletethisgain_d;

    /* Update for Memory: '<S3>/Memory' */
    my_pv_system_1_sm_master_DW.Memory_4_PreviousInput =
      my_pv_system_1_sm_master_B.SFunction_c[5];

    /* Update for Memory: '<S3>/Memory' */
    my_pv_system_1_sm_master_DW.Memory_5_PreviousInput =
      my_pv_system_1_sm_master_B.Temperature_l;

    /* Update for Memory: '<S3>/Memory' */
    my_pv_system_1_sm_master_DW.Memory_6_PreviousInput =
      my_pv_system_1_sm_master_B.donotdeletethisgain_db;

    /* Update for Memory: '<S3>/Memory1' */
    my_pv_system_1_sm_master_DW.Memory1_1_PreviousInput =
      my_pv_system_1_sm_master_B.donotdeletethisgain_g;

    /* Update for Memory: '<S3>/Memory1' */
    my_pv_system_1_sm_master_DW.Memory1_2_PreviousInput =
      my_pv_system_1_sm_master_B.donotdeletethisgain_e;

    /* Update for Memory: '<S3>/Memory1' */
    my_pv_system_1_sm_master_DW.Memory1_3_PreviousInput =
      my_pv_system_1_sm_master_B.donotdeletethisgain_i;

    /* Update for Memory: '<S3>/Memory1' */
    my_pv_system_1_sm_master_DW.Memory1_4_PreviousInput =
      my_pv_system_1_sm_master_B.donotdeletethisgain_de;

    /* Update for Memory: '<S3>/Memory1' */
    my_pv_system_1_sm_master_DW.Memory1_5_PreviousInput =
      my_pv_system_1_sm_master_B.donotdeletethisgain;

    /* Update for Memory: '<S3>/Memory1' */
    my_pv_system_1_sm_master_DW.Memory1_6_PreviousInput =
      my_pv_system_1_sm_master_B.Product1;

    /* Update for Memory: '<S3>/Memory1' */
    my_pv_system_1_sm_master_DW.Memory1_7_PreviousInput =
      my_pv_system_1_sm_master_B.Product;

    /* Update for Memory: '<S3>/Memory1' */
    my_pv_system_1_sm_master_DW.Memory1_8_PreviousInput[0] =
      my_pv_system_1_sm_master_B.SFunction_c[0];
    my_pv_system_1_sm_master_DW.Memory1_8_PreviousInput[1] =
      my_pv_system_1_sm_master_B.SFunction_c[1];
    my_pv_system_1_sm_master_DW.Memory1_8_PreviousInput[2] =
      my_pv_system_1_sm_master_B.SFunction_c[2];
    my_pv_system_1_sm_master_DW.Memory1_8_PreviousInput[3] =
      my_pv_system_1_sm_master_B.SFunction_c[3];

    /* Update for UnitDelay: '<S24>/Unit Delay' */
    my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[0] =
      my_pv_system_1_sm_master_B.Saturation[0];
    my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[1] =
      my_pv_system_1_sm_master_B.Saturation[1];
    my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[2] =
      my_pv_system_1_sm_master_B.Saturation[2];
    my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[3] =
      my_pv_system_1_sm_master_B.Saturation[3];

    /* Update for Enabled SubSystem: '<S24>/Tail' incorporates:
     *  Update for EnablePort: '<S25>/Enable'
     */
    if (my_pv_system_1_sm_master_DW.Tail_MODE) {
      /* Update for DiscreteIntegrator: '<S25>/Discrete-Time Integrator' incorporates:
       *  Constant: '<S25>/1'
       */
      tmp = my_pv_system_1_sm_master_P.DiscreteTimeIntegrator_gainval *
        my_pv_system_1_sm_master_P._Value;
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[0] += tmp;
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[1] += tmp;
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[2] += tmp;
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[3] += tmp;
      if (my_pv_system_1_sm_master_B.DataTypeConversion[0] > 0.0) {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[0] = 1;
      } else if (my_pv_system_1_sm_master_B.DataTypeConversion[0] < 0.0) {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[0] = -1;
      } else if (my_pv_system_1_sm_master_B.DataTypeConversion[0] == 0.0) {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[0] = 0;
      } else {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[0] = 2;
      }

      if (my_pv_system_1_sm_master_B.DataTypeConversion[1] > 0.0) {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[1] = 1;
      } else if (my_pv_system_1_sm_master_B.DataTypeConversion[1] < 0.0) {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[1] = -1;
      } else if (my_pv_system_1_sm_master_B.DataTypeConversion[1] == 0.0) {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[1] = 0;
      } else {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[1] = 2;
      }

      if (my_pv_system_1_sm_master_B.DataTypeConversion[2] > 0.0) {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[2] = 1;
      } else if (my_pv_system_1_sm_master_B.DataTypeConversion[2] < 0.0) {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[2] = -1;
      } else if (my_pv_system_1_sm_master_B.DataTypeConversion[2] == 0.0) {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[2] = 0;
      } else {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[2] = 2;
      }

      if (my_pv_system_1_sm_master_B.DataTypeConversion[3] > 0.0) {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[3] = 1;
      } else if (my_pv_system_1_sm_master_B.DataTypeConversion[3] < 0.0) {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[3] = -1;
      } else if (my_pv_system_1_sm_master_B.DataTypeConversion[3] == 0.0) {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[3] = 0;
      } else {
        my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[3] = 2;
      }

      /* End of Update for DiscreteIntegrator: '<S25>/Discrete-Time Integrator' */

      /* Update for UnitDelay: '<S25>/Unit Delay' */
      my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[0] =
        my_pv_system_1_sm_master_B.Switch_j[0];
      my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[1] =
        my_pv_system_1_sm_master_B.Switch_j[1];
      my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[2] =
        my_pv_system_1_sm_master_B.Switch_j[2];
      my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[3] =
        my_pv_system_1_sm_master_B.Switch_j[3];
    }

    /* End of Update for SubSystem: '<S24>/Tail' */
  }

  /* Update for TransportDelay: '<S29>/Transport Delay' */
  {
    real_T **uBuffer = (real_T**)
      &my_pv_system_1_sm_master_DW.TransportDelay_PWORK.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)
      &my_pv_system_1_sm_master_DW.TransportDelay_PWORK.TUbufferPtrs[1];
    real_T simTime = my_pv_system_1_sm_master_M->Timing.t[0];
    my_pv_system_1_sm_master_DW.TransportDelay_IWORK.Head =
      ((my_pv_system_1_sm_master_DW.TransportDelay_IWORK.Head <
        (my_pv_system_1_sm_master_DW.TransportDelay_IWORK.CircularBufSize-1)) ?
       (my_pv_system_1_sm_master_DW.TransportDelay_IWORK.Head+1) : 0);
    if (my_pv_system_1_sm_master_DW.TransportDelay_IWORK.Head ==
        my_pv_system_1_sm_master_DW.TransportDelay_IWORK.Tail) {
      my_pv_system_1_sm_master_DW.TransportDelay_IWORK.Tail =
        ((my_pv_system_1_sm_master_DW.TransportDelay_IWORK.Tail <
          (my_pv_system_1_sm_master_DW.TransportDelay_IWORK.CircularBufSize-1)) ?
         (my_pv_system_1_sm_master_DW.TransportDelay_IWORK.Tail+1) : 0);
    }

    (*tBuffer)[my_pv_system_1_sm_master_DW.TransportDelay_IWORK.Head] = simTime;
    (*uBuffer)[my_pv_system_1_sm_master_DW.TransportDelay_IWORK.Head] =
      my_pv_system_1_sm_master_B.integrator;
  }

  if (rtmIsMajorTimeStep(my_pv_system_1_sm_master_M)) {
    /* Update for Memory: '<S29>/Memory' */
    my_pv_system_1_sm_master_DW.Memory_PreviousInput =
      my_pv_system_1_sm_master_B.Switch_d;
  }

  /* Update for TransportDelay: '<S30>/Transport Delay' */
  {
    real_T **uBuffer = (real_T**)
      &my_pv_system_1_sm_master_DW.TransportDelay_PWORK_a.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)
      &my_pv_system_1_sm_master_DW.TransportDelay_PWORK_a.TUbufferPtrs[1];
    real_T simTime = my_pv_system_1_sm_master_M->Timing.t[0];
    my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.Head =
      ((my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.Head <
        (my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.CircularBufSize-1)) ?
       (my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.Head+1) : 0);
    if (my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.Head ==
        my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.Tail) {
      my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.Tail =
        ((my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.Tail <
          (my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.CircularBufSize-1))
         ? (my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.Tail+1) : 0);
    }

    (*tBuffer)[my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.Head] =
      simTime;
    (*uBuffer)[my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.Head] =
      my_pv_system_1_sm_master_B.integrator_i;
  }

  if (rtmIsMajorTimeStep(my_pv_system_1_sm_master_M)) {
    /* Update for Memory: '<S30>/Memory' */
    my_pv_system_1_sm_master_DW.Memory_PreviousInput_h =
      my_pv_system_1_sm_master_B.Switch_p;
  }

  if (rtmIsMajorTimeStep(my_pv_system_1_sm_master_M)) {
    rt_ertODEUpdateContinuousStates(&my_pv_system_1_sm_master_M->solverInfo);
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
  if (!(++my_pv_system_1_sm_master_M->Timing.clockTick0)) {
    ++my_pv_system_1_sm_master_M->Timing.clockTickH0;
  }

  my_pv_system_1_sm_master_M->Timing.t[0] = rtsiGetSolverStopTime
    (&my_pv_system_1_sm_master_M->solverInfo);

  {
    /* Update absolute timer for sample time: [5.0E-5s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++my_pv_system_1_sm_master_M->Timing.clockTick1)) {
      ++my_pv_system_1_sm_master_M->Timing.clockTickH1;
    }

    my_pv_system_1_sm_master_M->Timing.t[1] =
      my_pv_system_1_sm_master_M->Timing.clockTick1 *
      my_pv_system_1_sm_master_M->Timing.stepSize1 +
      my_pv_system_1_sm_master_M->Timing.clockTickH1 *
      my_pv_system_1_sm_master_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void my_pv_system_1_sm_master_derivatives(void)
{
  XDot_my_pv_system_1_sm_master_T *_rtXdot;
  _rtXdot = ((XDot_my_pv_system_1_sm_master_T *)
             my_pv_system_1_sm_master_M->ModelData.derivs);

  /* Derivatives for Integrator: '<S29>/integrator' */
  _rtXdot->integrator_CSTATE = my_pv_system_1_sm_master_B.donotdeletethisgain_d;

  /* Derivatives for Integrator: '<S30>/integrator' */
  _rtXdot->integrator_CSTATE_a =
    my_pv_system_1_sm_master_B.donotdeletethisgain_a;
}

/* Model initialize function */
void my_pv_system_1_sm_master_initialize(void)
{
  /* InitializeConditions for Enabled SubSystem: '<S24>/Tail' */
  /* InitializeConditions for DiscreteIntegrator: '<S25>/Discrete-Time Integrator' */
  my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[0] =
    my_pv_system_1_sm_master_P.DiscreteTimeIntegrator_IC;
  my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[1] =
    my_pv_system_1_sm_master_P.DiscreteTimeIntegrator_IC;
  my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[2] =
    my_pv_system_1_sm_master_P.DiscreteTimeIntegrator_IC;
  my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[3] =
    my_pv_system_1_sm_master_P.DiscreteTimeIntegrator_IC;
  my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[0] = 2;
  my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[1] = 2;
  my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[2] = 2;
  my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[3] = 2;

  /* InitializeConditions for UnitDelay: '<S25>/Unit Delay' */
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[0] =
    my_pv_system_1_sm_master_P.UnitDelay_InitialCondition;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[1] =
    my_pv_system_1_sm_master_P.UnitDelay_InitialCondition;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[2] =
    my_pv_system_1_sm_master_P.UnitDelay_InitialCondition;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[3] =
    my_pv_system_1_sm_master_P.UnitDelay_InitialCondition;

  /* End of InitializeConditions for SubSystem: '<S24>/Tail' */

  /* Start for Enabled SubSystem: '<S24>/Tail' */
  /* VirtualOutportStart for Outport: '<S25>/itail' */
  my_pv_system_1_sm_master_B.Product_d[0] = my_pv_system_1_sm_master_P.itail_Y0;
  my_pv_system_1_sm_master_B.Product_d[1] = my_pv_system_1_sm_master_P.itail_Y0;
  my_pv_system_1_sm_master_B.Product_d[2] = my_pv_system_1_sm_master_P.itail_Y0;
  my_pv_system_1_sm_master_B.Product_d[3] = my_pv_system_1_sm_master_P.itail_Y0;

  /* End of Start for SubSystem: '<S24>/Tail' */
  /* Level2 S-Function Block: '<S54>/State-Space' (fts5abcd_noncomp) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[5];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Start for TransportDelay: '<S29>/Transport Delay' */
  {
    real_T *pBuffer =
      &my_pv_system_1_sm_master_DW.TransportDelay_RWORK.TUbufferArea[0];
    my_pv_system_1_sm_master_DW.TransportDelay_IWORK.Tail = 0;
    my_pv_system_1_sm_master_DW.TransportDelay_IWORK.Head = 0;
    my_pv_system_1_sm_master_DW.TransportDelay_IWORK.Last = 0;
    my_pv_system_1_sm_master_DW.TransportDelay_IWORK.CircularBufSize = 8192;
    pBuffer[0] = my_pv_system_1_sm_master_P.TransportDelay_InitOutput;
    pBuffer[8192] = my_pv_system_1_sm_master_M->Timing.t[0];
    my_pv_system_1_sm_master_DW.TransportDelay_PWORK.TUbufferPtrs[0] = (void *)
      &pBuffer[0];
    my_pv_system_1_sm_master_DW.TransportDelay_PWORK.TUbufferPtrs[1] = (void *)
      &pBuffer[8192];
  }

  /* Start for TransportDelay: '<S30>/Transport Delay' */
  {
    real_T *pBuffer =
      &my_pv_system_1_sm_master_DW.TransportDelay_RWORK_o.TUbufferArea[0];
    my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.Tail = 0;
    my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.Head = 0;
    my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.Last = 0;
    my_pv_system_1_sm_master_DW.TransportDelay_IWORK_i.CircularBufSize = 8192;
    pBuffer[0] = my_pv_system_1_sm_master_P.TransportDelay_InitOutput_i;
    pBuffer[8192] = my_pv_system_1_sm_master_M->Timing.t[0];
    my_pv_system_1_sm_master_DW.TransportDelay_PWORK_a.TUbufferPtrs[0] = (void *)
      &pBuffer[0];
    my_pv_system_1_sm_master_DW.TransportDelay_PWORK_a.TUbufferPtrs[1] = (void *)
      &pBuffer[8192];
  }

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

  /* InitializeConditions for Memory: '<S2>/S-Function' */
  my_pv_system_1_sm_master_DW.SFunction_PreviousInput =
    my_pv_system_1_sm_master_P.SFunction_X0;

  /* InitializeConditions for Memory: '<S3>/Memory' */
  my_pv_system_1_sm_master_DW.Memory_1_PreviousInput =
    my_pv_system_1_sm_master_P.Memory_1_X0;

  /* InitializeConditions for Memory: '<S3>/Memory' */
  my_pv_system_1_sm_master_DW.Memory_2_PreviousInput =
    my_pv_system_1_sm_master_P.Memory_2_X0;

  /* InitializeConditions for Memory: '<S3>/Memory' */
  my_pv_system_1_sm_master_DW.Memory_3_PreviousInput =
    my_pv_system_1_sm_master_P.Memory_3_X0;

  /* InitializeConditions for Memory: '<S3>/Memory' */
  my_pv_system_1_sm_master_DW.Memory_4_PreviousInput =
    my_pv_system_1_sm_master_P.Memory_4_X0;

  /* InitializeConditions for Memory: '<S3>/Memory' */
  my_pv_system_1_sm_master_DW.Memory_5_PreviousInput =
    my_pv_system_1_sm_master_P.Memory_5_X0;

  /* InitializeConditions for Memory: '<S3>/Memory' */
  my_pv_system_1_sm_master_DW.Memory_6_PreviousInput =
    my_pv_system_1_sm_master_P.Memory_6_X0;

  /* InitializeConditions for Memory: '<S3>/Memory1' */
  my_pv_system_1_sm_master_DW.Memory1_1_PreviousInput =
    my_pv_system_1_sm_master_P.Memory1_1_X0;

  /* InitializeConditions for Memory: '<S3>/Memory1' */
  my_pv_system_1_sm_master_DW.Memory1_2_PreviousInput =
    my_pv_system_1_sm_master_P.Memory1_2_X0;

  /* InitializeConditions for Memory: '<S3>/Memory1' */
  my_pv_system_1_sm_master_DW.Memory1_3_PreviousInput =
    my_pv_system_1_sm_master_P.Memory1_3_X0;

  /* InitializeConditions for Memory: '<S3>/Memory1' */
  my_pv_system_1_sm_master_DW.Memory1_4_PreviousInput =
    my_pv_system_1_sm_master_P.Memory1_4_X0;

  /* InitializeConditions for Memory: '<S3>/Memory1' */
  my_pv_system_1_sm_master_DW.Memory1_5_PreviousInput =
    my_pv_system_1_sm_master_P.Memory1_5_X0;

  /* InitializeConditions for Memory: '<S3>/Memory1' */
  my_pv_system_1_sm_master_DW.Memory1_6_PreviousInput =
    my_pv_system_1_sm_master_P.Memory1_6_X0;

  /* InitializeConditions for Memory: '<S3>/Memory1' */
  my_pv_system_1_sm_master_DW.Memory1_7_PreviousInput =
    my_pv_system_1_sm_master_P.Memory1_7_X0;

  /* InitializeConditions for Memory: '<S3>/Memory1' */
  my_pv_system_1_sm_master_DW.Memory1_8_PreviousInput[0] =
    my_pv_system_1_sm_master_P.Memory1_8_X0;
  my_pv_system_1_sm_master_DW.Memory1_8_PreviousInput[1] =
    my_pv_system_1_sm_master_P.Memory1_8_X0;
  my_pv_system_1_sm_master_DW.Memory1_8_PreviousInput[2] =
    my_pv_system_1_sm_master_P.Memory1_8_X0;
  my_pv_system_1_sm_master_DW.Memory1_8_PreviousInput[3] =
    my_pv_system_1_sm_master_P.Memory1_8_X0;

  /* InitializeConditions for Atomic SubSystem: '<S16>/Subsystem5' */

  /* Level2 S-Function Block: '<S50>/S-Function' (send_rt) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[0];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* End of InitializeConditions for SubSystem: '<S16>/Subsystem5' */

  /* Level2 S-Function Block: '<S48>/S-Function' (OP_SEND) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[1];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Level2 S-Function Block: '<S3>/OpMonitor' (opmonitor) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[2];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Level2 S-Function Block: '<S49>/S-Function' (OP_SEND) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[3];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* InitializeConditions for UnitDelay: '<S24>/Unit Delay' */
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[0] =
    my_pv_system_1_sm_master_P.UnitDelay_InitialCondition_p;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[1] =
    my_pv_system_1_sm_master_P.UnitDelay_InitialCondition_p;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[2] =
    my_pv_system_1_sm_master_P.UnitDelay_InitialCondition_p;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[3] =
    my_pv_system_1_sm_master_P.UnitDelay_InitialCondition_p;

  /* Level2 S-Function Block: '<S51>/S-Function' (recv_rt) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[4];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Level2 S-Function Block: '<S54>/State-Space' (fts5abcd_noncomp) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[5];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* InitializeConditions for Integrator: '<S29>/integrator' */
  my_pv_system_1_sm_master_X.integrator_CSTATE =
    my_pv_system_1_sm_master_P.integrator_IC;

  /* InitializeConditions for Memory: '<S29>/Memory' */
  my_pv_system_1_sm_master_DW.Memory_PreviousInput =
    my_pv_system_1_sm_master_P.Memory_X0;

  /* InitializeConditions for Integrator: '<S30>/integrator' */
  my_pv_system_1_sm_master_X.integrator_CSTATE_a =
    my_pv_system_1_sm_master_P.integrator_IC_i;

  /* InitializeConditions for Memory: '<S30>/Memory' */
  my_pv_system_1_sm_master_DW.Memory_PreviousInput_h =
    my_pv_system_1_sm_master_P.Memory_X0_k;

  /* Level2 S-Function Block: '<S31>/S-Function' (RECV_Param) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[6];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* InitializeConditions for RateLimiter: '<S12>/Rate Limiter' */
  my_pv_system_1_sm_master_DW.PrevY = my_pv_system_1_sm_master_P.RateLimiter_IC;
}

/* Model terminate function */
void my_pv_system_1_sm_master_terminate(void)
{
  /* Terminate for Atomic SubSystem: '<S16>/Subsystem5' */

  /* Level2 S-Function Block: '<S50>/S-Function' (send_rt) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[0];
    sfcnTerminate(rts);
  }

  /* End of Terminate for SubSystem: '<S16>/Subsystem5' */

  /* Level2 S-Function Block: '<S48>/S-Function' (OP_SEND) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[1];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S3>/OpMonitor' (opmonitor) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[2];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S49>/S-Function' (OP_SEND) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[3];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S51>/S-Function' (recv_rt) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[4];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S54>/State-Space' (fts5abcd_noncomp) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[5];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S31>/S-Function' (RECV_Param) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[6];
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
  my_pv_system_1_sm_master_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  my_pv_system_1_sm_master_update();
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
  my_pv_system_1_sm_master_initialize();
}

void MdlTerminate(void)
{
  my_pv_system_1_sm_master_terminate();
}

/* Registration function */
RT_MODEL_my_pv_system_1_sm_master_T *my_pv_system_1_sm_master(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  my_pv_system_1_sm_master_P.Saturation_UpperSat = rtInf;

  /* initialize real-time model */
  (void) memset((void *)my_pv_system_1_sm_master_M, 0,
                sizeof(RT_MODEL_my_pv_system_1_sm_master_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&my_pv_system_1_sm_master_M->solverInfo,
                          &my_pv_system_1_sm_master_M->Timing.simTimeStep);
    rtsiSetTPtr(&my_pv_system_1_sm_master_M->solverInfo, &rtmGetTPtr
                (my_pv_system_1_sm_master_M));
    rtsiSetStepSizePtr(&my_pv_system_1_sm_master_M->solverInfo,
                       &my_pv_system_1_sm_master_M->Timing.stepSize0);
    rtsiSetdXPtr(&my_pv_system_1_sm_master_M->solverInfo,
                 &my_pv_system_1_sm_master_M->ModelData.derivs);
    rtsiSetContStatesPtr(&my_pv_system_1_sm_master_M->solverInfo, (real_T **)
                         &my_pv_system_1_sm_master_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&my_pv_system_1_sm_master_M->solverInfo,
      &my_pv_system_1_sm_master_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&my_pv_system_1_sm_master_M->solverInfo,
                          (&rtmGetErrorStatus(my_pv_system_1_sm_master_M)));
    rtsiSetRTModelPtr(&my_pv_system_1_sm_master_M->solverInfo,
                      my_pv_system_1_sm_master_M);
  }

  rtsiSetSimTimeStep(&my_pv_system_1_sm_master_M->solverInfo, MAJOR_TIME_STEP);
  my_pv_system_1_sm_master_M->ModelData.intgData.y =
    my_pv_system_1_sm_master_M->ModelData.odeY;
  my_pv_system_1_sm_master_M->ModelData.intgData.f[0] =
    my_pv_system_1_sm_master_M->ModelData.odeF[0];
  my_pv_system_1_sm_master_M->ModelData.intgData.f[1] =
    my_pv_system_1_sm_master_M->ModelData.odeF[1];
  my_pv_system_1_sm_master_M->ModelData.intgData.f[2] =
    my_pv_system_1_sm_master_M->ModelData.odeF[2];
  my_pv_system_1_sm_master_M->ModelData.contStates = ((real_T *)
    &my_pv_system_1_sm_master_X);
  rtsiSetSolverData(&my_pv_system_1_sm_master_M->solverInfo, (void *)
                    &my_pv_system_1_sm_master_M->ModelData.intgData);
  rtsiSetSolverName(&my_pv_system_1_sm_master_M->solverInfo,"ode3");
  my_pv_system_1_sm_master_M->solverInfoPtr =
    (&my_pv_system_1_sm_master_M->solverInfo);

  /* Initialize timing info */
  {
    int_T *mdlTsMap = my_pv_system_1_sm_master_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    my_pv_system_1_sm_master_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    my_pv_system_1_sm_master_M->Timing.sampleTimes =
      (&my_pv_system_1_sm_master_M->Timing.sampleTimesArray[0]);
    my_pv_system_1_sm_master_M->Timing.offsetTimes =
      (&my_pv_system_1_sm_master_M->Timing.offsetTimesArray[0]);

    /* task periods */
    my_pv_system_1_sm_master_M->Timing.sampleTimes[0] = (0.0);
    my_pv_system_1_sm_master_M->Timing.sampleTimes[1] = (5.0E-5);

    /* task offsets */
    my_pv_system_1_sm_master_M->Timing.offsetTimes[0] = (0.0);
    my_pv_system_1_sm_master_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(my_pv_system_1_sm_master_M,
             &my_pv_system_1_sm_master_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = my_pv_system_1_sm_master_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    my_pv_system_1_sm_master_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(my_pv_system_1_sm_master_M, 10.0);
  my_pv_system_1_sm_master_M->Timing.stepSize0 = 5.0E-5;
  my_pv_system_1_sm_master_M->Timing.stepSize1 = 5.0E-5;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    my_pv_system_1_sm_master_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(my_pv_system_1_sm_master_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(my_pv_system_1_sm_master_M->rtwLogInfo, (NULL));
    rtliSetLogT(my_pv_system_1_sm_master_M->rtwLogInfo, "");
    rtliSetLogX(my_pv_system_1_sm_master_M->rtwLogInfo, "");
    rtliSetLogXFinal(my_pv_system_1_sm_master_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(my_pv_system_1_sm_master_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(my_pv_system_1_sm_master_M->rtwLogInfo, 0);
    rtliSetLogMaxRows(my_pv_system_1_sm_master_M->rtwLogInfo, 1000);
    rtliSetLogDecimation(my_pv_system_1_sm_master_M->rtwLogInfo, 1);
    rtliSetLogY(my_pv_system_1_sm_master_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(my_pv_system_1_sm_master_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(my_pv_system_1_sm_master_M->rtwLogInfo, (NULL));
  }

  my_pv_system_1_sm_master_M->solverInfoPtr =
    (&my_pv_system_1_sm_master_M->solverInfo);
  my_pv_system_1_sm_master_M->Timing.stepSize = (5.0E-5);
  rtsiSetFixedStepSize(&my_pv_system_1_sm_master_M->solverInfo, 5.0E-5);
  rtsiSetSolverMode(&my_pv_system_1_sm_master_M->solverInfo,
                    SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  my_pv_system_1_sm_master_M->ModelData.blockIO = ((void *)
    &my_pv_system_1_sm_master_B);
  (void) memset(((void *) &my_pv_system_1_sm_master_B), 0,
                sizeof(B_my_pv_system_1_sm_master_T));

  {
    int_T i;
    for (i = 0; i < 6; i++) {
      my_pv_system_1_sm_master_B.SFunction_c[i] = 0.0;
    }

    for (i = 0; i < 12; i++) {
      my_pv_system_1_sm_master_B.StateSpace_o1[i] = 0.0;
    }

    my_pv_system_1_sm_master_B.SFunction = 0.0;
    my_pv_system_1_sm_master_B.Sum = 0.0;
    my_pv_system_1_sm_master_B.IrradianceWm2 = 0.0;
    my_pv_system_1_sm_master_B.I_PV = 0.0;
    my_pv_system_1_sm_master_B.V_PV = 0.0;
    my_pv_system_1_sm_master_B.I_Diode = 0.0;
    my_pv_system_1_sm_master_B.Temperature = 0.0;
    my_pv_system_1_sm_master_B.Vph = 0.0;
    my_pv_system_1_sm_master_B.Vinv = 0.0;
    my_pv_system_1_sm_master_B.Vdc = 0.0;
    my_pv_system_1_sm_master_B.Igrid = 0.0;
    my_pv_system_1_sm_master_B.Vgrid = 0.0;
    my_pv_system_1_sm_master_B.Ig = 0.0;
    my_pv_system_1_sm_master_B.P_PV = 0.0;
    my_pv_system_1_sm_master_B.P_PV_MEAN = 0.0;
    my_pv_system_1_sm_master_B.Gate[0] = 0.0;
    my_pv_system_1_sm_master_B.Gate[1] = 0.0;
    my_pv_system_1_sm_master_B.Gate[2] = 0.0;
    my_pv_system_1_sm_master_B.Gate[3] = 0.0;
    my_pv_system_1_sm_master_B.Computationtime = 0.0;
    my_pv_system_1_sm_master_B.Realstepsize = 0.0;
    my_pv_system_1_sm_master_B.Idletime = 0.0;
    my_pv_system_1_sm_master_B.Overruntimes = 0.0;
    my_pv_system_1_sm_master_B.UnitDelay[0] = 0.0;
    my_pv_system_1_sm_master_B.UnitDelay[1] = 0.0;
    my_pv_system_1_sm_master_B.UnitDelay[2] = 0.0;
    my_pv_system_1_sm_master_B.UnitDelay[3] = 0.0;
    my_pv_system_1_sm_master_B.DataTypeConversion[0] = 0.0;
    my_pv_system_1_sm_master_B.DataTypeConversion[1] = 0.0;
    my_pv_system_1_sm_master_B.DataTypeConversion[2] = 0.0;
    my_pv_system_1_sm_master_B.DataTypeConversion[3] = 0.0;
    my_pv_system_1_sm_master_B.AC = 0.0;
    my_pv_system_1_sm_master_B.StateSpace_o2[0] = 0.0;
    my_pv_system_1_sm_master_B.StateSpace_o2[1] = 0.0;
    my_pv_system_1_sm_master_B.StateSpace_o2[2] = 0.0;
    my_pv_system_1_sm_master_B.StateSpace_o2[3] = 0.0;
    my_pv_system_1_sm_master_B.Ron[0] = 0.0;
    my_pv_system_1_sm_master_B.Ron[1] = 0.0;
    my_pv_system_1_sm_master_B.Ron[2] = 0.0;
    my_pv_system_1_sm_master_B.Ron[3] = 0.0;
    my_pv_system_1_sm_master_B.Switch[0] = 0.0;
    my_pv_system_1_sm_master_B.Switch[1] = 0.0;
    my_pv_system_1_sm_master_B.Switch[2] = 0.0;
    my_pv_system_1_sm_master_B.Switch[3] = 0.0;
    my_pv_system_1_sm_master_B.Saturation[0] = 0.0;
    my_pv_system_1_sm_master_B.Saturation[1] = 0.0;
    my_pv_system_1_sm_master_B.Saturation[2] = 0.0;
    my_pv_system_1_sm_master_B.Saturation[3] = 0.0;
    my_pv_system_1_sm_master_B.Switch_e[0] = 0.0;
    my_pv_system_1_sm_master_B.Switch_e[1] = 0.0;
    my_pv_system_1_sm_master_B.Switch_e[2] = 0.0;
    my_pv_system_1_sm_master_B.Switch_e[3] = 0.0;
    my_pv_system_1_sm_master_B.donotdeletethisgain = 0.0;
    my_pv_system_1_sm_master_B.Clock = 0.0;
    my_pv_system_1_sm_master_B.integrator = 0.0;
    my_pv_system_1_sm_master_B.TransportDelay = 0.0;
    my_pv_system_1_sm_master_B.Memory = 0.0;
    my_pv_system_1_sm_master_B.Switch_d = 0.0;
    my_pv_system_1_sm_master_B.Clock_c = 0.0;
    my_pv_system_1_sm_master_B.integrator_i = 0.0;
    my_pv_system_1_sm_master_B.TransportDelay_j = 0.0;
    my_pv_system_1_sm_master_B.Memory_j = 0.0;
    my_pv_system_1_sm_master_B.Switch_p = 0.0;
    my_pv_system_1_sm_master_B.SFunction_f[0] = 0.0;
    my_pv_system_1_sm_master_B.SFunction_f[1] = 0.0;
    my_pv_system_1_sm_master_B.Product = 0.0;
    my_pv_system_1_sm_master_B.donotdeletethisgain_d = 0.0;
    my_pv_system_1_sm_master_B.donotdeletethisgain_a = 0.0;
    my_pv_system_1_sm_master_B.Product1 = 0.0;
    my_pv_system_1_sm_master_B.IrradianceWm2_j = 0.0;
    my_pv_system_1_sm_master_B.Temperature_l = 0.0;
    my_pv_system_1_sm_master_B.donotdeletethisgain_db = 0.0;
    my_pv_system_1_sm_master_B.donotdeletethisgain_e = 0.0;
    my_pv_system_1_sm_master_B.donotdeletethisgain_g = 0.0;
    my_pv_system_1_sm_master_B.donotdeletethisgain_i = 0.0;
    my_pv_system_1_sm_master_B.donotdeletethisgain_de = 0.0;
    my_pv_system_1_sm_master_B.Sum_p = 0.0;
    my_pv_system_1_sm_master_B.Gain = 0.0;
    my_pv_system_1_sm_master_B.Sum_g = 0.0;
    my_pv_system_1_sm_master_B.Gain_f = 0.0;
    my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[0] = 0.0;
    my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[1] = 0.0;
    my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[2] = 0.0;
    my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[3] = 0.0;
    my_pv_system_1_sm_master_B.u9Tf[0] = 0.0;
    my_pv_system_1_sm_master_B.u9Tf[1] = 0.0;
    my_pv_system_1_sm_master_B.u9Tf[2] = 0.0;
    my_pv_system_1_sm_master_B.u9Tf[3] = 0.0;
    my_pv_system_1_sm_master_B.Add[0] = 0.0;
    my_pv_system_1_sm_master_B.Add[1] = 0.0;
    my_pv_system_1_sm_master_B.Add[2] = 0.0;
    my_pv_system_1_sm_master_B.Add[3] = 0.0;
    my_pv_system_1_sm_master_B.Saturation1[0] = 0.0;
    my_pv_system_1_sm_master_B.Saturation1[1] = 0.0;
    my_pv_system_1_sm_master_B.Saturation1[2] = 0.0;
    my_pv_system_1_sm_master_B.Saturation1[3] = 0.0;
    my_pv_system_1_sm_master_B.Add1[0] = 0.0;
    my_pv_system_1_sm_master_B.Add1[1] = 0.0;
    my_pv_system_1_sm_master_B.Add1[2] = 0.0;
    my_pv_system_1_sm_master_B.Add1[3] = 0.0;
    my_pv_system_1_sm_master_B.uTt[0] = 0.0;
    my_pv_system_1_sm_master_B.uTt[1] = 0.0;
    my_pv_system_1_sm_master_B.uTt[2] = 0.0;
    my_pv_system_1_sm_master_B.uTt[3] = 0.0;
    my_pv_system_1_sm_master_B.Saturation2[0] = 0.0;
    my_pv_system_1_sm_master_B.Saturation2[1] = 0.0;
    my_pv_system_1_sm_master_B.Saturation2[2] = 0.0;
    my_pv_system_1_sm_master_B.Saturation2[3] = 0.0;
    my_pv_system_1_sm_master_B.Add2[0] = 0.0;
    my_pv_system_1_sm_master_B.Add2[1] = 0.0;
    my_pv_system_1_sm_master_B.Add2[2] = 0.0;
    my_pv_system_1_sm_master_B.Add2[3] = 0.0;
    my_pv_system_1_sm_master_B.UnitDelay_l[0] = 0.0;
    my_pv_system_1_sm_master_B.UnitDelay_l[1] = 0.0;
    my_pv_system_1_sm_master_B.UnitDelay_l[2] = 0.0;
    my_pv_system_1_sm_master_B.UnitDelay_l[3] = 0.0;
    my_pv_system_1_sm_master_B.Switch_j[0] = 0.0;
    my_pv_system_1_sm_master_B.Switch_j[1] = 0.0;
    my_pv_system_1_sm_master_B.Switch_j[2] = 0.0;
    my_pv_system_1_sm_master_B.Switch_j[3] = 0.0;
    my_pv_system_1_sm_master_B.Product_d[0] = 0.0;
    my_pv_system_1_sm_master_B.Product_d[1] = 0.0;
    my_pv_system_1_sm_master_B.Product_d[2] = 0.0;
    my_pv_system_1_sm_master_B.Product_d[3] = 0.0;
  }

  /* parameters */
  my_pv_system_1_sm_master_M->ModelData.defaultParam = ((real_T *)
    &my_pv_system_1_sm_master_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &my_pv_system_1_sm_master_X;
    my_pv_system_1_sm_master_M->ModelData.contStates = (x);
    (void) memset((void *)&my_pv_system_1_sm_master_X, 0,
                  sizeof(X_my_pv_system_1_sm_master_T));
  }

  /* states (dwork) */
  my_pv_system_1_sm_master_M->ModelData.dwork = ((void *)
    &my_pv_system_1_sm_master_DW);
  (void) memset((void *)&my_pv_system_1_sm_master_DW, 0,
                sizeof(DW_my_pv_system_1_sm_master_T));
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[0] = 0.0;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[1] = 0.0;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[2] = 0.0;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[3] = 0.0;

  {
    int_T i;
    for (i = 0; i < 20; i++) {
      my_pv_system_1_sm_master_DW.StateSpace_DSTATE[i] = 0.0;
    }
  }

  my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[0] = 0.0;
  my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[1] = 0.0;
  my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[2] = 0.0;
  my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[3] = 0.0;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[0] = 0.0;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[1] = 0.0;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[2] = 0.0;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[3] = 0.0;
  my_pv_system_1_sm_master_DW.SFunction_PreviousInput = 0.0;
  my_pv_system_1_sm_master_DW.Memory_1_PreviousInput = 0.0;
  my_pv_system_1_sm_master_DW.Memory_2_PreviousInput = 0.0;
  my_pv_system_1_sm_master_DW.Memory_3_PreviousInput = 0.0;
  my_pv_system_1_sm_master_DW.Memory_4_PreviousInput = 0.0;
  my_pv_system_1_sm_master_DW.Memory_5_PreviousInput = 0.0;
  my_pv_system_1_sm_master_DW.Memory_6_PreviousInput = 0.0;
  my_pv_system_1_sm_master_DW.Memory1_1_PreviousInput = 0.0;
  my_pv_system_1_sm_master_DW.Memory1_2_PreviousInput = 0.0;
  my_pv_system_1_sm_master_DW.Memory1_3_PreviousInput = 0.0;
  my_pv_system_1_sm_master_DW.Memory1_4_PreviousInput = 0.0;
  my_pv_system_1_sm_master_DW.Memory1_5_PreviousInput = 0.0;
  my_pv_system_1_sm_master_DW.Memory1_6_PreviousInput = 0.0;
  my_pv_system_1_sm_master_DW.Memory1_7_PreviousInput = 0.0;
  my_pv_system_1_sm_master_DW.Memory1_8_PreviousInput[0] = 0.0;
  my_pv_system_1_sm_master_DW.Memory1_8_PreviousInput[1] = 0.0;
  my_pv_system_1_sm_master_DW.Memory1_8_PreviousInput[2] = 0.0;
  my_pv_system_1_sm_master_DW.Memory1_8_PreviousInput[3] = 0.0;
  my_pv_system_1_sm_master_DW.Memory_PreviousInput = 0.0;
  my_pv_system_1_sm_master_DW.Memory_PreviousInput_h = 0.0;
  my_pv_system_1_sm_master_DW.PrevY = 0.0;

  {
    int_T i;
    for (i = 0; i < 74; i++) {
      my_pv_system_1_sm_master_DW.StateSpace_RWORK[i] = 0.0;
    }
  }

  my_pv_system_1_sm_master_DW.TransportDelay_RWORK.modelTStart = 0.0;

  {
    int_T i;
    for (i = 0; i < 16384; i++) {
      my_pv_system_1_sm_master_DW.TransportDelay_RWORK.TUbufferArea[i] = 0.0;
    }
  }

  my_pv_system_1_sm_master_DW.TransportDelay_RWORK_o.modelTStart = 0.0;

  {
    int_T i;
    for (i = 0; i < 16384; i++) {
      my_pv_system_1_sm_master_DW.TransportDelay_RWORK_o.TUbufferArea[i] = 0.0;
    }
  }

  /* child S-Function registration */
  {
    RTWSfcnInfo *sfcnInfo =
      &my_pv_system_1_sm_master_M->NonInlinedSFcns.sfcnInfo;
    my_pv_system_1_sm_master_M->sfcnInfo = (sfcnInfo);
    rtssSetErrorStatusPtr(sfcnInfo, (&rtmGetErrorStatus
      (my_pv_system_1_sm_master_M)));
    rtssSetNumRootSampTimesPtr(sfcnInfo,
      &my_pv_system_1_sm_master_M->Sizes.numSampTimes);
    my_pv_system_1_sm_master_M->NonInlinedSFcns.taskTimePtrs[0] = &(rtmGetTPtr
      (my_pv_system_1_sm_master_M)[0]);
    my_pv_system_1_sm_master_M->NonInlinedSFcns.taskTimePtrs[1] = &(rtmGetTPtr
      (my_pv_system_1_sm_master_M)[1]);
    rtssSetTPtrPtr(sfcnInfo,
                   my_pv_system_1_sm_master_M->NonInlinedSFcns.taskTimePtrs);
    rtssSetTStartPtr(sfcnInfo, &rtmGetTStart(my_pv_system_1_sm_master_M));
    rtssSetTFinalPtr(sfcnInfo, &rtmGetTFinal(my_pv_system_1_sm_master_M));
    rtssSetTimeOfLastOutputPtr(sfcnInfo, &rtmGetTimeOfLastOutput
      (my_pv_system_1_sm_master_M));
    rtssSetStepSizePtr(sfcnInfo, &my_pv_system_1_sm_master_M->Timing.stepSize);
    rtssSetStopRequestedPtr(sfcnInfo, &rtmGetStopRequested
      (my_pv_system_1_sm_master_M));
    rtssSetDerivCacheNeedsResetPtr(sfcnInfo,
      &my_pv_system_1_sm_master_M->ModelData.derivCacheNeedsReset);
    rtssSetZCCacheNeedsResetPtr(sfcnInfo,
      &my_pv_system_1_sm_master_M->ModelData.zCCacheNeedsReset);
    rtssSetBlkStateChangePtr(sfcnInfo,
      &my_pv_system_1_sm_master_M->ModelData.blkStateChange);
    rtssSetSampleHitsPtr(sfcnInfo,
                         &my_pv_system_1_sm_master_M->Timing.sampleHits);
    rtssSetPerTaskSampleHitsPtr(sfcnInfo,
      &my_pv_system_1_sm_master_M->Timing.perTaskSampleHits);
    rtssSetSimModePtr(sfcnInfo, &my_pv_system_1_sm_master_M->simMode);
    rtssSetSolverInfoPtr(sfcnInfo, &my_pv_system_1_sm_master_M->solverInfoPtr);
  }

  my_pv_system_1_sm_master_M->Sizes.numSFcns = (7);

  /* register each child */
  {
    (void) memset((void *)
                  &my_pv_system_1_sm_master_M->NonInlinedSFcns.childSFunctions[0],
                  0,
                  7*sizeof(SimStruct));
    my_pv_system_1_sm_master_M->childSfunctions =
      (&my_pv_system_1_sm_master_M->NonInlinedSFcns.childSFunctionPtrs[0]);

    {
      int_T i;
      for (i = 0; i < 7; i++) {
        my_pv_system_1_sm_master_M->childSfunctions[i] =
          (&my_pv_system_1_sm_master_M->NonInlinedSFcns.childSFunctions[i]);
      }
    }

    /* Level2 S-Function Block: my_pv_system_1_sm_master/<S50>/S-Function (send_rt) */
    {
      SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[0];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn0.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn0.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn0.sfcnTsMap;
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
                         &my_pv_system_1_sm_master_M->NonInlinedSFcns.blkInfo2[0]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_1_sm_master_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_1_sm_master_M->
                           NonInlinedSFcns.methods2[0]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_1_sm_master_M->
                           NonInlinedSFcns.methods3[0]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_1_sm_master_M->NonInlinedSFcns.statesInfo2
                         [0]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn0.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn0.UPtrs0;
          sfcnUPtrs[0] = &my_pv_system_1_sm_master_B.IrradianceWm2;
          sfcnUPtrs[1] = &my_pv_system_1_sm_master_B.I_PV;
          sfcnUPtrs[2] = &my_pv_system_1_sm_master_B.V_PV;
          sfcnUPtrs[3] = &my_pv_system_1_sm_master_B.I_Diode;
          sfcnUPtrs[4] = &my_pv_system_1_sm_master_B.Temperature;
          sfcnUPtrs[5] = &my_pv_system_1_sm_master_B.Vph;
          sfcnUPtrs[6] = &my_pv_system_1_sm_master_B.Vinv;
          sfcnUPtrs[7] = &my_pv_system_1_sm_master_B.Vdc;
          sfcnUPtrs[8] = &my_pv_system_1_sm_master_B.Igrid;
          sfcnUPtrs[9] = &my_pv_system_1_sm_master_B.Vgrid;
          sfcnUPtrs[10] = &my_pv_system_1_sm_master_B.Ig;
          sfcnUPtrs[11] = &my_pv_system_1_sm_master_B.P_PV;
          sfcnUPtrs[12] = &my_pv_system_1_sm_master_B.P_PV_MEAN;
          sfcnUPtrs[13] = &my_pv_system_1_sm_master_B.Gate[0];
          sfcnUPtrs[14] = &my_pv_system_1_sm_master_B.Gate[1];
          sfcnUPtrs[15] = &my_pv_system_1_sm_master_B.Gate[2];
          sfcnUPtrs[16] = &my_pv_system_1_sm_master_B.Gate[3];
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 17);
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "my_pv_system_1_sm_master/SM_master/rtlab_send_subsystem/Subsystem5/Send5/S-Function");
      ssSetRTModel(rts,my_pv_system_1_sm_master_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn0.params;
        ssSetSFcnParamsCount(rts, 3);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_1_sm_master_P.SFunction_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_1_sm_master_P.SFunction_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       my_pv_system_1_sm_master_P.SFunction_P3_Size);
      }

      /* work vectors */
      ssSetIWork(rts, (int_T *) &my_pv_system_1_sm_master_DW.SFunction_IWORK_c[0]);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn0.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn0.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* IWORK */
        ssSetDWorkWidth(rts, 0, 5);
        ssSetDWorkDataType(rts, 0,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_1_sm_master_DW.SFunction_IWORK_c[0]);
      }

      /* registration */
      send_rt(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 17);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: my_pv_system_1_sm_master/<S48>/S-Function (OP_SEND) */
    {
      SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[1];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn1.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn1.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn1.sfcnTsMap;
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
                         &my_pv_system_1_sm_master_M->NonInlinedSFcns.blkInfo2[1]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_1_sm_master_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_1_sm_master_M->
                           NonInlinedSFcns.methods2[1]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_1_sm_master_M->
                           NonInlinedSFcns.methods3[1]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_1_sm_master_M->NonInlinedSFcns.statesInfo2
                         [1]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn1.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn1.UPtrs0;
          sfcnUPtrs[0] = &my_pv_system_1_sm_master_B.IrradianceWm2;
          sfcnUPtrs[1] = &my_pv_system_1_sm_master_B.I_PV;
          sfcnUPtrs[2] = &my_pv_system_1_sm_master_B.V_PV;
          sfcnUPtrs[3] = &my_pv_system_1_sm_master_B.I_Diode;
          sfcnUPtrs[4] = &my_pv_system_1_sm_master_B.Temperature;
          sfcnUPtrs[5] = &my_pv_system_1_sm_master_B.Vph;
          sfcnUPtrs[6] = &my_pv_system_1_sm_master_B.Vinv;
          sfcnUPtrs[7] = &my_pv_system_1_sm_master_B.Vdc;
          sfcnUPtrs[8] = &my_pv_system_1_sm_master_B.Igrid;
          sfcnUPtrs[9] = &my_pv_system_1_sm_master_B.Vgrid;
          sfcnUPtrs[10] = &my_pv_system_1_sm_master_B.Ig;
          sfcnUPtrs[11] = &my_pv_system_1_sm_master_B.P_PV;
          sfcnUPtrs[12] = &my_pv_system_1_sm_master_B.P_PV_MEAN;
          sfcnUPtrs[13] = &my_pv_system_1_sm_master_B.Gate[0];
          sfcnUPtrs[14] = &my_pv_system_1_sm_master_B.Gate[1];
          sfcnUPtrs[15] = &my_pv_system_1_sm_master_B.Gate[2];
          sfcnUPtrs[16] = &my_pv_system_1_sm_master_B.Gate[3];
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 17);
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "my_pv_system_1_sm_master/SM_master/rtlab_send_subsystem/Subsystem1/Send1/S-Function");
      ssSetRTModel(rts,my_pv_system_1_sm_master_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn1.params;
        ssSetSFcnParamsCount(rts, 1);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_1_sm_master_P.SFunction_P1_Size_c);
      }

      /* work vectors */
      ssSetIWork(rts, (int_T *) &my_pv_system_1_sm_master_DW.SFunction_IWORK[0]);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn1.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn1.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* IWORK */
        ssSetDWorkWidth(rts, 0, 5);
        ssSetDWorkDataType(rts, 0,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_1_sm_master_DW.SFunction_IWORK[0]);
      }

      /* registration */
      OP_SEND(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 17);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: my_pv_system_1_sm_master/<S3>/OpMonitor (opmonitor) */
    {
      SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[2];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn2.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn2.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn2.sfcnTsMap;
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
                         &my_pv_system_1_sm_master_M->NonInlinedSFcns.blkInfo2[2]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_1_sm_master_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_1_sm_master_M->
                           NonInlinedSFcns.methods2[2]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_1_sm_master_M->
                           NonInlinedSFcns.methods3[2]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_1_sm_master_M->NonInlinedSFcns.statesInfo2
                         [2]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn2.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn2.UPtrs0;
          sfcnUPtrs[0] = (real_T*)&my_pv_system_1_sm_master_RGND;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn2.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 4);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &my_pv_system_1_sm_master_B.Computationtime));
        }

        /* port 1 */
        {
          _ssSetOutputPortNumDimensions(rts, 1, 1);
          ssSetOutputPortWidth(rts, 1, 1);
          ssSetOutputPortSignal(rts, 1, ((real_T *)
            &my_pv_system_1_sm_master_B.Realstepsize));
        }

        /* port 2 */
        {
          _ssSetOutputPortNumDimensions(rts, 2, 1);
          ssSetOutputPortWidth(rts, 2, 1);
          ssSetOutputPortSignal(rts, 2, ((real_T *)
            &my_pv_system_1_sm_master_B.Idletime));
        }

        /* port 3 */
        {
          _ssSetOutputPortNumDimensions(rts, 3, 1);
          ssSetOutputPortWidth(rts, 3, 1);
          ssSetOutputPortSignal(rts, 3, ((real_T *)
            &my_pv_system_1_sm_master_B.Overruntimes));
        }
      }

      /* path info */
      ssSetModelName(rts, "OpMonitor");
      ssSetPath(rts, "my_pv_system_1_sm_master/SM_master/OpMonitor");
      ssSetRTModel(rts,my_pv_system_1_sm_master_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn2.params;
        ssSetSFcnParamsCount(rts, 6);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_1_sm_master_P.OpMonitor_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_1_sm_master_P.OpMonitor_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       my_pv_system_1_sm_master_P.OpMonitor_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       my_pv_system_1_sm_master_P.OpMonitor_P4_Size);
        ssSetSFcnParam(rts, 4, (mxArray*)
                       my_pv_system_1_sm_master_P.OpMonitor_P5_Size);
        ssSetSFcnParam(rts, 5, (mxArray*)
                       my_pv_system_1_sm_master_P.OpMonitor_P6_Size);
      }

      /* work vectors */
      ssSetPWork(rts, (void **) &my_pv_system_1_sm_master_DW.OpMonitor_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn2.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn2.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_1_sm_master_DW.OpMonitor_PWORK);
      }

      /* registration */
      opmonitor(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
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

    /* Level2 S-Function Block: my_pv_system_1_sm_master/<S49>/S-Function (OP_SEND) */
    {
      SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[3];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn3.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn3.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn3.sfcnTsMap;
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
                         &my_pv_system_1_sm_master_M->NonInlinedSFcns.blkInfo2[3]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_1_sm_master_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_1_sm_master_M->
                           NonInlinedSFcns.methods2[3]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_1_sm_master_M->
                           NonInlinedSFcns.methods3[3]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_1_sm_master_M->NonInlinedSFcns.statesInfo2
                         [3]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn3.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn3.UPtrs0;
          sfcnUPtrs[0] = &my_pv_system_1_sm_master_B.Computationtime;
          sfcnUPtrs[1] = &my_pv_system_1_sm_master_B.Realstepsize;
          sfcnUPtrs[2] = &my_pv_system_1_sm_master_B.Idletime;
          sfcnUPtrs[3] = &my_pv_system_1_sm_master_B.Overruntimes;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 4);
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "my_pv_system_1_sm_master/SM_master/rtlab_send_subsystem/Subsystem2/Send2/S-Function");
      ssSetRTModel(rts,my_pv_system_1_sm_master_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn3.params;
        ssSetSFcnParamsCount(rts, 1);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_1_sm_master_P.SFunction_P1_Size_b);
      }

      /* work vectors */
      ssSetIWork(rts, (int_T *) &my_pv_system_1_sm_master_DW.SFunction_IWORK_n[0]);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn3.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn3.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* IWORK */
        ssSetDWorkWidth(rts, 0, 5);
        ssSetDWorkDataType(rts, 0,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_1_sm_master_DW.SFunction_IWORK_n[0]);
      }

      /* registration */
      OP_SEND(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
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

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: my_pv_system_1_sm_master/<S51>/S-Function (recv_rt) */
    {
      SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[4];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn4.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn4.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn4.sfcnTsMap;
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
                         &my_pv_system_1_sm_master_M->NonInlinedSFcns.blkInfo2[4]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_1_sm_master_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_1_sm_master_M->
                           NonInlinedSFcns.methods2[4]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_1_sm_master_M->
                           NonInlinedSFcns.methods3[4]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_1_sm_master_M->NonInlinedSFcns.statesInfo2
                         [4]);
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn4.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 6);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            my_pv_system_1_sm_master_B.SFunction_c));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "my_pv_system_1_sm_master/SM_master/zzzOpComm1/Receive_1/S-Function");
      ssSetRTModel(rts,my_pv_system_1_sm_master_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn4.params;
        ssSetSFcnParamsCount(rts, 3);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_1_sm_master_P.SFunction_P1_Size_b3);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_1_sm_master_P.SFunction_P2_Size_o);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       my_pv_system_1_sm_master_P.SFunction_P3_Size_c);
      }

      /* work vectors */
      ssSetPWork(rts, (void **) &my_pv_system_1_sm_master_DW.SFunction_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn4.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn4.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_1_sm_master_DW.SFunction_PWORK);
      }

      /* registration */
      recv_rt(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
    }

    /* Level2 S-Function Block: my_pv_system_1_sm_master/<S54>/State-Space (fts5abcd_noncomp) */
    {
      SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[5];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn5.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn5.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn5.sfcnTsMap;
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
                         &my_pv_system_1_sm_master_M->NonInlinedSFcns.blkInfo2[5]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_1_sm_master_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_1_sm_master_M->
                           NonInlinedSFcns.methods2[5]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_1_sm_master_M->
                           NonInlinedSFcns.methods3[5]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_1_sm_master_M->NonInlinedSFcns.statesInfo2
                         [5]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn5.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn5.UPtrs0;
          sfcnUPtrs[0] = &my_pv_system_1_sm_master_B.Product_d[0];
          sfcnUPtrs[1] = &my_pv_system_1_sm_master_B.Product_d[1];
          sfcnUPtrs[2] = &my_pv_system_1_sm_master_B.Product_d[2];
          sfcnUPtrs[3] = &my_pv_system_1_sm_master_B.Product_d[3];
          sfcnUPtrs[4] = &my_pv_system_1_sm_master_B.AC;
          sfcnUPtrs[5] = &my_pv_system_1_sm_master_B.SFunction_c[5];
          sfcnUPtrs[6] = &my_pv_system_1_sm_master_B.SFunction_c[4];
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 7);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn5.UPtrs1;
          sfcnUPtrs[0] = my_pv_system_1_sm_master_B.DataTypeConversion;
          sfcnUPtrs[1] = &my_pv_system_1_sm_master_B.DataTypeConversion[1];
          sfcnUPtrs[2] = &my_pv_system_1_sm_master_B.DataTypeConversion[2];
          sfcnUPtrs[3] = &my_pv_system_1_sm_master_B.DataTypeConversion[3];
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 4);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn5.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 2);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 12);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            my_pv_system_1_sm_master_B.StateSpace_o1));
        }

        /* port 1 */
        {
          _ssSetOutputPortNumDimensions(rts, 1, 1);
          ssSetOutputPortWidth(rts, 1, 4);
          ssSetOutputPortSignal(rts, 1, ((real_T *)
            my_pv_system_1_sm_master_B.StateSpace_o2));
        }
      }

      /* states */
      ssSetDiscStates(rts, (real_T *)
                      &my_pv_system_1_sm_master_DW.StateSpace_DSTATE[0]);

      /* path info */
      ssSetModelName(rts, "State-Space");
      ssSetPath(rts,
                "my_pv_system_1_sm_master/powergui/EquivalentModel1/State-Space");
      ssSetRTModel(rts,my_pv_system_1_sm_master_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn5.params;
        ssSetSFcnParamsCount(rts, 23);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P4_Size);
        ssSetSFcnParam(rts, 4, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P5_Size);
        ssSetSFcnParam(rts, 5, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P6_Size);
        ssSetSFcnParam(rts, 6, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P7_Size);
        ssSetSFcnParam(rts, 7, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P8_Size);
        ssSetSFcnParam(rts, 8, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P9_Size);
        ssSetSFcnParam(rts, 9, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P10_Size);
        ssSetSFcnParam(rts, 10, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P11_Size);
        ssSetSFcnParam(rts, 11, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P12_Size);
        ssSetSFcnParam(rts, 12, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P13_Size);
        ssSetSFcnParam(rts, 13, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P14_Size);
        ssSetSFcnParam(rts, 14, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P15_Size);
        ssSetSFcnParam(rts, 15, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P16_Size);
        ssSetSFcnParam(rts, 16, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P17_Size);
        ssSetSFcnParam(rts, 17, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P18_Size);
        ssSetSFcnParam(rts, 18, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P19_Size);
        ssSetSFcnParam(rts, 19, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P20_Size);
        ssSetSFcnParam(rts, 20, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P21_Size);
        ssSetSFcnParam(rts, 21, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P22_Size);
        ssSetSFcnParam(rts, 22, (mxArray*)
                       my_pv_system_1_sm_master_P.StateSpace_P23_Size);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *) &my_pv_system_1_sm_master_DW.StateSpace_RWORK[0]);
      ssSetIWork(rts, (int_T *) &my_pv_system_1_sm_master_DW.StateSpace_IWORK[0]);
      ssSetPWork(rts, (void **) &my_pv_system_1_sm_master_DW.StateSpace_PWORK[0]);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn5.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn5.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 4);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 74);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_1_sm_master_DW.StateSpace_RWORK[0]);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 317);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &my_pv_system_1_sm_master_DW.StateSpace_IWORK[0]);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 27);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &my_pv_system_1_sm_master_DW.StateSpace_PWORK[0]);

        /* DSTATE */
        ssSetDWorkWidth(rts, 3, 20);
        ssSetDWorkDataType(rts, 3,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 3, 0);
        ssSetDWorkUsedAsDState(rts, 3, 1);
        ssSetDWork(rts, 3, &my_pv_system_1_sm_master_DW.StateSpace_DSTATE[0]);
      }

      /* registration */
      fts5abcd_noncomp(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 1, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);
      _ssSetOutputPortBeingMerged(rts, 1, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
    }

    /* Level2 S-Function Block: my_pv_system_1_sm_master/<S31>/S-Function (RECV_Param) */
    {
      SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[6];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn6.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn6.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn6.sfcnTsMap;
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
                         &my_pv_system_1_sm_master_M->NonInlinedSFcns.blkInfo2[6]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_1_sm_master_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_1_sm_master_M->
                           NonInlinedSFcns.methods2[6]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_1_sm_master_M->
                           NonInlinedSFcns.methods3[6]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_1_sm_master_M->NonInlinedSFcns.statesInfo2
                         [6]);
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn6.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 2);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            my_pv_system_1_sm_master_B.SFunction_f));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "my_pv_system_1_sm_master/SM_master/OpComm/Receive/S-Function");
      ssSetRTModel(rts,my_pv_system_1_sm_master_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn6.params;
        ssSetSFcnParamsCount(rts, 2);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_1_sm_master_P.SFunction_P1_Size_k);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_1_sm_master_P.SFunction_P2_Size_k);
      }

      /* registration */
      RECV_Param(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
    }
  }

  /* Initialize Sizes */
  my_pv_system_1_sm_master_M->Sizes.numContStates = (2);/* Number of continuous states */
  my_pv_system_1_sm_master_M->Sizes.numY = (0);/* Number of model outputs */
  my_pv_system_1_sm_master_M->Sizes.numU = (0);/* Number of model inputs */
  my_pv_system_1_sm_master_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  my_pv_system_1_sm_master_M->Sizes.numSampTimes = (2);/* Number of sample times */
  my_pv_system_1_sm_master_M->Sizes.numBlocks = (153);/* Number of blocks */
  my_pv_system_1_sm_master_M->Sizes.numBlockIO = (70);/* Number of block outputs */
  my_pv_system_1_sm_master_M->Sizes.numBlockPrms = (703);/* Sum of parameter "widths" */
  return my_pv_system_1_sm_master_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
