/*
 * my_pv_system_3_ss_controller.c
 *
 * Code generation for model "my_pv_system_3_ss_controller".
 *
 * Model version              : 1.207
 * Simulink Coder version : 8.7 (R2014b) 08-Sep-2014
 * C source code generated on : Wed Sep 07 16:46:24 2016
 *
 * Target selection: rtlab_rtmodel.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "my_pv_system_3_ss_controller.h"
#include "my_pv_system_3_ss_controller_private.h"

const real_T my_pv_system_3_ss_controller_RGND = 0.0;/* real_T ground */

/* Block signals (auto storage) */
B_my_pv_system_3_ss_controller_T my_pv_system_3_ss_controller_B;

/* Block states (auto storage) */
DW_my_pv_system_3_ss_controller_T my_pv_system_3_ss_controller_DW;

/* Real-time model */
RT_MODEL_my_pv_system_3_ss_controller_T my_pv_system_3_ss_controller_M_;
RT_MODEL_my_pv_system_3_ss_controller_T *const my_pv_system_3_ss_controller_M =
  &my_pv_system_3_ss_controller_M_;

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
static void my_pv_system_3_ss_controller_output(void)
{
  real_T P;
  real_T dV;
  real_T dP;
  real_T u0;
  real_T u2;

  /* Memory: '<S1>/S-Function' */
  my_pv_system_3_ss_controller_B.SFunction =
    my_pv_system_3_ss_controller_DW.SFunction_PreviousInput;

  /* Sum: '<S1>/Sum' incorporates:
   *  Constant: '<S1>/S-Function1'
   */
  my_pv_system_3_ss_controller_B.Sum =
    my_pv_system_3_ss_controller_P.SFunction1_Value +
    my_pv_system_3_ss_controller_B.SFunction;

  /* Stop: '<S1>/Stop Simulation' */
  if (my_pv_system_3_ss_controller_B.Sum != 0.0) {
    rtmSetStopRequested(my_pv_system_3_ss_controller_M, 1);
  }

  /* End of Stop: '<S1>/Stop Simulation' */

  /* Memory: '<S2>/Memory' */
  my_pv_system_3_ss_controller_B.Gate[0] =
    my_pv_system_3_ss_controller_DW.Memory_1_PreviousInput[0];
  my_pv_system_3_ss_controller_B.Gate[1] =
    my_pv_system_3_ss_controller_DW.Memory_1_PreviousInput[1];
  my_pv_system_3_ss_controller_B.Gate[2] =
    my_pv_system_3_ss_controller_DW.Memory_1_PreviousInput[2];
  my_pv_system_3_ss_controller_B.Gate[3] =
    my_pv_system_3_ss_controller_DW.Memory_1_PreviousInput[3];

  /* Memory: '<S2>/Memory' */
  my_pv_system_3_ss_controller_B.Iph =
    my_pv_system_3_ss_controller_DW.Memory_2_PreviousInput;

  /* Memory: '<S2>/Memory' */
  my_pv_system_3_ss_controller_B.Io =
    my_pv_system_3_ss_controller_DW.Memory_3_PreviousInput;

  /* Outputs for Atomic SubSystem: '<S6>/Subsystem4' */

  /* Level2 S-Function Block: '<S90>/S-Function' (send_rt) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[2];
    sfcnOutputs(rts, 0);
  }

  /* End of Outputs for SubSystem: '<S6>/Subsystem4' */

  /* Level2 S-Function Block: '<S2>/OpMonitor' (opmonitor) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[3];
    sfcnOutputs(rts, 0);
  }

  /* Level2 S-Function Block: '<S89>/S-Function' (OP_SEND) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[4];
    sfcnOutputs(rts, 0);
  }

  /* Level2 S-Function Block: '<S91>/S-Function' (recv_rt) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[5];
    sfcnOutputs(rts, 0);
  }

  /* RateTransition: '<S2>/Rate Transition1' */
  my_pv_system_3_ss_controller_B.RateTransition1 =
    my_pv_system_3_ss_controller_B.SFunction_g[8];

  /* Gain: '<S11>/A->pu' */
  P = my_pv_system_3_ss_controller_P.InverterControl_Vnom_prim;
  dP = my_pv_system_3_ss_controller_P.InverterControl_Pnom;
  dV = P / dP;
  dV /= 1.4142135623730951;
  my_pv_system_3_ss_controller_B.Apu = dV *
    my_pv_system_3_ss_controller_B.RateTransition1;

  /* UnitDelay: '<S31>/Unit Delay' */
  my_pv_system_3_ss_controller_B.UnitDelay =
    my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE;

  /* Saturate: '<S20>/avoid division by zero' */
  u0 = my_pv_system_3_ss_controller_B.UnitDelay;
  dP = my_pv_system_3_ss_controller_P.avoiddivisionbyzero_LowerSat;
  u2 = my_pv_system_3_ss_controller_P.avoiddivisionbyzero_UpperSat;
  if (u0 > u2) {
    my_pv_system_3_ss_controller_B.avoiddivisionbyzero = u2;
  } else if (u0 < dP) {
    my_pv_system_3_ss_controller_B.avoiddivisionbyzero = dP;
  } else {
    my_pv_system_3_ss_controller_B.avoiddivisionbyzero = u0;
  }

  /* End of Saturate: '<S20>/avoid division by zero' */

  /* Math: '<S20>/Math Function'
   *
   * About '<S20>/Math Function':
   *  Operator: reciprocal
   */
  P = my_pv_system_3_ss_controller_B.avoiddivisionbyzero;
  my_pv_system_3_ss_controller_B.MathFunction = 1.0 / P;

  /* Gain: '<S20>/Gain' */
  my_pv_system_3_ss_controller_B.Gain =
    my_pv_system_3_ss_controller_P.Gain_Gain_c *
    my_pv_system_3_ss_controller_B.MathFunction;

  /* Level2 S-Function Block: '<S57>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[6];
    sfcnOutputs(rts, 0);
  }

  /* DiscreteIntegrator: '<S31>/Discrete-Time Integrator' */
  my_pv_system_3_ss_controller_B.DiscreteTimeIntegrator =
    my_pv_system_3_ss_controller_DW.DiscreteTimeIntegrator_DSTATE;

  /* Math: '<S31>/Math Function' incorporates:
   *  Constant: '<S31>/Constant4'
   */
  my_pv_system_3_ss_controller_B.MathFunction_i = rt_modd_snf
    (my_pv_system_3_ss_controller_B.DiscreteTimeIntegrator,
     my_pv_system_3_ss_controller_P.Constant4_Value_i);

  /* RelationalOperator: '<S58>/Compare' incorporates:
   *  Constant: '<S56>/Constant'
   *  Constant: '<S58>/Constant'
   */
  my_pv_system_3_ss_controller_B.Compare = (uint8_T)
    (my_pv_system_3_ss_controller_P.AlphaBetaZerotodq0_Alignment ==
     my_pv_system_3_ss_controller_P.CompareToConstant_const);

  /* Outputs for Enabled SubSystem: '<S56>/Subsystem1' incorporates:
   *  EnablePort: '<S61>/Enable'
   */
  if (my_pv_system_3_ss_controller_B.Compare > 0) {
    /* Fcn: '<S61>/Fcn' */
    my_pv_system_3_ss_controller_B.Fcn = my_pv_system_3_ss_controller_B.Apu *
      cos(my_pv_system_3_ss_controller_B.MathFunction_i) +
      my_pv_system_3_ss_controller_B.SFunction_o * sin
      (my_pv_system_3_ss_controller_B.MathFunction_i);

    /* Fcn: '<S61>/Fcn1' */
    my_pv_system_3_ss_controller_B.Fcn1 = -my_pv_system_3_ss_controller_B.Apu *
      sin(my_pv_system_3_ss_controller_B.MathFunction_i) +
      my_pv_system_3_ss_controller_B.SFunction_o * cos
      (my_pv_system_3_ss_controller_B.MathFunction_i);
  }

  /* End of Outputs for SubSystem: '<S56>/Subsystem1' */

  /* RelationalOperator: '<S59>/Compare' incorporates:
   *  Constant: '<S56>/Constant'
   *  Constant: '<S59>/Constant'
   */
  my_pv_system_3_ss_controller_B.Compare_e = (uint8_T)
    (my_pv_system_3_ss_controller_P.AlphaBetaZerotodq0_Alignment ==
     my_pv_system_3_ss_controller_P.CompareToConstant1_const);

  /* Outputs for Enabled SubSystem: '<S56>/Subsystem - pi//2 delay' incorporates:
   *  EnablePort: '<S60>/Enable'
   */
  if (my_pv_system_3_ss_controller_B.Compare_e > 0) {
    /* Fcn: '<S60>/Fcn' */
    my_pv_system_3_ss_controller_B.Fcn_l = my_pv_system_3_ss_controller_B.Apu *
      sin(my_pv_system_3_ss_controller_B.MathFunction_i) -
      my_pv_system_3_ss_controller_B.SFunction_o * cos
      (my_pv_system_3_ss_controller_B.MathFunction_i);

    /* Fcn: '<S60>/Fcn1' */
    my_pv_system_3_ss_controller_B.Fcn1_f = my_pv_system_3_ss_controller_B.Apu *
      cos(my_pv_system_3_ss_controller_B.MathFunction_i) +
      my_pv_system_3_ss_controller_B.SFunction_o * sin
      (my_pv_system_3_ss_controller_B.MathFunction_i);
  }

  /* End of Outputs for SubSystem: '<S56>/Subsystem - pi//2 delay' */

  /* Step: '<S20>/First cycle of simulation Id=0.92, Iq=0' */
  u0 = my_pv_system_3_ss_controller_M->Timing.t[0];
  dP = my_pv_system_3_ss_controller_P.InverterControl_Fnom;
  dV = 1.0 / dP;
  if (u0 < dV) {
    my_pv_system_3_ss_controller_B.FirstcycleofsimulationId092Iq0 =
      my_pv_system_3_ss_controller_P.FirstcycleofsimulationId092Iq0_;
  } else {
    my_pv_system_3_ss_controller_B.FirstcycleofsimulationId092Iq0 =
      my_pv_system_3_ss_controller_P.FirstcycleofsimulationId092Iq_o;
  }

  /* End of Step: '<S20>/First cycle of simulation Id=0.92, Iq=0' */

  /* Switch: '<S20>/Switch' incorporates:
   *  Constant: '<S20>/Constant'
   */
  if (my_pv_system_3_ss_controller_B.FirstcycleofsimulationId092Iq0 != 0.0) {
    /* Switch: '<S56>/Switch' */
    if (my_pv_system_3_ss_controller_B.Compare != 0) {
      my_pv_system_3_ss_controller_B.Switch_p[0] =
        my_pv_system_3_ss_controller_B.Fcn;
    } else {
      my_pv_system_3_ss_controller_B.Switch_p[0] =
        my_pv_system_3_ss_controller_B.Fcn_l;
    }

    if (my_pv_system_3_ss_controller_B.Compare != 0) {
      my_pv_system_3_ss_controller_B.Switch_p[1] =
        my_pv_system_3_ss_controller_B.Fcn1;
    } else {
      my_pv_system_3_ss_controller_B.Switch_p[1] =
        my_pv_system_3_ss_controller_B.Fcn1_f;
    }

    /* End of Switch: '<S56>/Switch' */
    my_pv_system_3_ss_controller_B.Switch[0] =
      my_pv_system_3_ss_controller_B.Switch_p[0];
    my_pv_system_3_ss_controller_B.Switch[1] =
      my_pv_system_3_ss_controller_B.Switch_p[1];
  } else {
    my_pv_system_3_ss_controller_B.Switch[0] =
      my_pv_system_3_ss_controller_P.Constant_Value_a[0];
    my_pv_system_3_ss_controller_B.Switch[1] =
      my_pv_system_3_ss_controller_P.Constant_Value_a[1];
  }

  /* End of Switch: '<S20>/Switch' */

  /* Gain: '<S52>/D*u(k)' */
  my_pv_system_3_ss_controller_B.Duk[0] =
    my_pv_system_3_ss_controller_P.Duk_Gain *
    my_pv_system_3_ss_controller_B.Switch[0];
  my_pv_system_3_ss_controller_B.Duk[1] =
    my_pv_system_3_ss_controller_P.Duk_Gain *
    my_pv_system_3_ss_controller_B.Switch[1];

  /* UnitDelay: '<S52>/Delay_x1' */
  my_pv_system_3_ss_controller_B.x1k[0] =
    my_pv_system_3_ss_controller_DW.Delay_x1_DSTATE[0];
  my_pv_system_3_ss_controller_B.x1k[1] =
    my_pv_system_3_ss_controller_DW.Delay_x1_DSTATE[1];

  /* Gain: '<S55>/C11' */
  my_pv_system_3_ss_controller_B.C11[0] =
    my_pv_system_3_ss_controller_P.C11_Gain *
    my_pv_system_3_ss_controller_B.x1k[0];
  my_pv_system_3_ss_controller_B.C11[1] =
    my_pv_system_3_ss_controller_P.C11_Gain *
    my_pv_system_3_ss_controller_B.x1k[1];

  /* UnitDelay: '<S52>/Delay_x2' */
  my_pv_system_3_ss_controller_B.x2k[0] =
    my_pv_system_3_ss_controller_DW.Delay_x2_DSTATE[0];
  my_pv_system_3_ss_controller_B.x2k[1] =
    my_pv_system_3_ss_controller_DW.Delay_x2_DSTATE[1];

  /* Gain: '<S55>/C12' */
  my_pv_system_3_ss_controller_B.C12[0] =
    my_pv_system_3_ss_controller_P.C12_Gain *
    my_pv_system_3_ss_controller_B.x2k[0];
  my_pv_system_3_ss_controller_B.C12[1] =
    my_pv_system_3_ss_controller_P.C12_Gain *
    my_pv_system_3_ss_controller_B.x2k[1];

  /* Sum: '<S55>/sum2' */
  my_pv_system_3_ss_controller_B.sum2[0] = my_pv_system_3_ss_controller_B.C11[0]
    + my_pv_system_3_ss_controller_B.C12[0];
  my_pv_system_3_ss_controller_B.sum2[1] = my_pv_system_3_ss_controller_B.C11[1]
    + my_pv_system_3_ss_controller_B.C12[1];

  /* Sum: '<S52>/C*X(k)+D*u(k)' */
  my_pv_system_3_ss_controller_B.yk[0] = my_pv_system_3_ss_controller_B.Duk[0] +
    my_pv_system_3_ss_controller_B.sum2[0];
  my_pv_system_3_ss_controller_B.yk[1] = my_pv_system_3_ss_controller_B.Duk[1] +
    my_pv_system_3_ss_controller_B.sum2[1];

  /* UnitDelay: '<S4>/Unit Delay2' */
  my_pv_system_3_ss_controller_B.UnitDelay2 =
    my_pv_system_3_ss_controller_DW.UnitDelay2_DSTATE;

  /* Sum: '<S8>/Sum' incorporates:
   *  Constant: '<S4>/Iq_ref'
   */
  my_pv_system_3_ss_controller_B.Sum_e[0] =
    my_pv_system_3_ss_controller_B.UnitDelay2 -
    my_pv_system_3_ss_controller_B.yk[0];
  my_pv_system_3_ss_controller_B.Sum_e[1] =
    my_pv_system_3_ss_controller_P.Iq_ref_Value -
    my_pv_system_3_ss_controller_B.yk[1];

  /* Gain: '<S15>/Proportional Gain' */
  my_pv_system_3_ss_controller_B.ProportionalGain[0] =
    my_pv_system_3_ss_controller_P.InverterControl_Kp_Ireg *
    my_pv_system_3_ss_controller_B.Sum_e[0];
  my_pv_system_3_ss_controller_B.ProportionalGain[1] =
    my_pv_system_3_ss_controller_P.InverterControl_Kp_Ireg *
    my_pv_system_3_ss_controller_B.Sum_e[1];

  /* DiscreteIntegrator: '<S15>/Integrator' */
  my_pv_system_3_ss_controller_B.Integrator[0] =
    my_pv_system_3_ss_controller_DW.Integrator_DSTATE[0];
  my_pv_system_3_ss_controller_B.Integrator[1] =
    my_pv_system_3_ss_controller_DW.Integrator_DSTATE[1];

  /* Sum: '<S15>/Sum' */
  my_pv_system_3_ss_controller_B.Sum_g[0] =
    my_pv_system_3_ss_controller_B.ProportionalGain[0] +
    my_pv_system_3_ss_controller_B.Integrator[0];
  my_pv_system_3_ss_controller_B.Sum_g[1] =
    my_pv_system_3_ss_controller_B.ProportionalGain[1] +
    my_pv_system_3_ss_controller_B.Integrator[1];

  /* Saturate: '<S15>/Saturate' */
  u0 = my_pv_system_3_ss_controller_B.Sum_g[0];
  dP = my_pv_system_3_ss_controller_P.PI_LowerSaturationLimit;
  u2 = my_pv_system_3_ss_controller_P.PI_UpperSaturationLimit;
  if (u0 > u2) {
    u0 = u2;
  } else {
    if (u0 < dP) {
      u0 = dP;
    }
  }

  my_pv_system_3_ss_controller_B.Saturate[0] = u0;
  u0 = my_pv_system_3_ss_controller_B.Sum_g[1];
  dP = my_pv_system_3_ss_controller_P.PI_LowerSaturationLimit;
  u2 = my_pv_system_3_ss_controller_P.PI_UpperSaturationLimit;
  if (u0 > u2) {
    u0 = u2;
  } else {
    if (u0 < dP) {
      u0 = dP;
    }
  }

  my_pv_system_3_ss_controller_B.Saturate[1] = u0;

  /* End of Saturate: '<S15>/Saturate' */

  /* RateTransition: '<S2>/Rate Transition' */
  my_pv_system_3_ss_controller_B.RateTransition =
    my_pv_system_3_ss_controller_B.SFunction_g[9];

  /* Gain: '<S11>/V->pu' */
  dP = my_pv_system_3_ss_controller_P.InverterControl_Vnom_prim *
    1.4142135623730951;
  dV = 1.0 / dP;
  my_pv_system_3_ss_controller_B.Vpu = dV *
    my_pv_system_3_ss_controller_B.RateTransition;

  /* Trigonometry: '<S16>/Trigonometric Function' */
  my_pv_system_3_ss_controller_B.TrigonometricFunction = sin
    (my_pv_system_3_ss_controller_B.MathFunction_i);

  /* Gain: '<S16>/Gain1' */
  my_pv_system_3_ss_controller_B.Gain1 =
    my_pv_system_3_ss_controller_P.Gain1_Gain_l *
    my_pv_system_3_ss_controller_B.TrigonometricFunction;

  /* Product: '<S16>/Product1' */
  my_pv_system_3_ss_controller_B.Product1 = my_pv_system_3_ss_controller_B.Vpu *
    my_pv_system_3_ss_controller_B.Gain1;

  /* DiscreteIntegrator: '<S23>/Integ4' */
  if (my_pv_system_3_ss_controller_DW.Integ4_SYSTEM_ENABLE != 0) {
    my_pv_system_3_ss_controller_B.Integ4 =
      my_pv_system_3_ss_controller_DW.Integ4_DSTATE;
  } else {
    my_pv_system_3_ss_controller_B.Integ4 =
      my_pv_system_3_ss_controller_P.Integ4_gainval_o *
      my_pv_system_3_ss_controller_B.Product1 +
      my_pv_system_3_ss_controller_DW.Integ4_DSTATE;
  }

  /* End of DiscreteIntegrator: '<S23>/Integ4' */

  /* Saturate: '<S23>/To avoid division  by zero' */
  u0 = my_pv_system_3_ss_controller_B.UnitDelay;
  dP = my_pv_system_3_ss_controller_P.Toavoiddivisionbyzero_LowerSa_l;
  u2 = my_pv_system_3_ss_controller_P.Toavoiddivisionbyzero_UpperSa_j;
  if (u0 > u2) {
    my_pv_system_3_ss_controller_B.Freq = u2;
  } else if (u0 < dP) {
    my_pv_system_3_ss_controller_B.Freq = dP;
  } else {
    my_pv_system_3_ss_controller_B.Freq = u0;
  }

  /* End of Saturate: '<S23>/To avoid division  by zero' */

  /* Fcn: '<S23>/Number of samples per cycle' */
  my_pv_system_3_ss_controller_B.Numberofsamplespercycle = 1.0 /
    my_pv_system_3_ss_controller_B.Freq / 5.0e-5;

  /* Rounding: '<S23>/Rounding Function' */
  my_pv_system_3_ss_controller_B.RoundingFunction = ceil
    (my_pv_system_3_ss_controller_B.Numberofsamplespercycle);

  /* Gain: '<S23>/Gain' */
  my_pv_system_3_ss_controller_B.Delay =
    my_pv_system_3_ss_controller_P.InverterControl_Ts_Control *
    my_pv_system_3_ss_controller_B.RoundingFunction;

  /* Level2 S-Function Block: '<S25>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[7];
    sfcnOutputs(rts, 0);
  }

  /* UnitDelay: '<S24>/Unit Delay' */
  my_pv_system_3_ss_controller_B.UnitDelay_f =
    my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_g;

  /* DigitalClock: '<S23>/Digital  Clock' */
  my_pv_system_3_ss_controller_B.DigitalClock =
    my_pv_system_3_ss_controller_M->Timing.t[0];

  /* RelationalOperator: '<S23>/Relational Operator' incorporates:
   *  Constant: '<S23>/Constant'
   */
  my_pv_system_3_ss_controller_B.RelationalOperator =
    (my_pv_system_3_ss_controller_B.DigitalClock >=
     my_pv_system_3_ss_controller_P.Constant_Value_p);

  /* UnitDelay: '<S23>/Unit Delay1' */
  my_pv_system_3_ss_controller_B.UnitDelay1 =
    my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE;

  /* Switch: '<S23>/Switch' */
  if (my_pv_system_3_ss_controller_B.RelationalOperator) {
    /* Sum: '<S24>/Sum1' */
    my_pv_system_3_ss_controller_B.Sum1_c3 =
      my_pv_system_3_ss_controller_B.Product1 -
      my_pv_system_3_ss_controller_B.UnitDelay_f;

    /* Sum: '<S24>/Sum5' */
    my_pv_system_3_ss_controller_B.Sum5_cr =
      my_pv_system_3_ss_controller_B.Numberofsamplespercycle -
      my_pv_system_3_ss_controller_B.RoundingFunction;

    /* Product: '<S24>/Product5' */
    my_pv_system_3_ss_controller_B.Product5_j =
      my_pv_system_3_ss_controller_B.Sum5_cr *
      my_pv_system_3_ss_controller_B.Sum1_c3;

    /* Gain: '<S24>/Gain1' */
    my_pv_system_3_ss_controller_B.Gain1_b =
      my_pv_system_3_ss_controller_P.Gain1_Gain *
      my_pv_system_3_ss_controller_B.Product5_j;

    /* Sum: '<S24>/Sum4' */
    my_pv_system_3_ss_controller_B.Sum4_e =
      my_pv_system_3_ss_controller_B.Gain1_b +
      my_pv_system_3_ss_controller_B.Product1;

    /* Product: '<S24>/Product2' */
    my_pv_system_3_ss_controller_B.Product2_f =
      my_pv_system_3_ss_controller_B.Sum5_cr /
      my_pv_system_3_ss_controller_B.Numberofsamplespercycle;

    /* Product: '<S24>/Product4' */
    my_pv_system_3_ss_controller_B.Product4_c =
      my_pv_system_3_ss_controller_B.Product2_f *
      my_pv_system_3_ss_controller_B.Sum4_e;

    /* Sum: '<S23>/Sum7' */
    my_pv_system_3_ss_controller_B.Sum7_eg =
      my_pv_system_3_ss_controller_B.Integ4 -
      my_pv_system_3_ss_controller_B.SFunction_p;

    /* Product: '<S23>/Product' */
    my_pv_system_3_ss_controller_B.Meanvalue_c =
      my_pv_system_3_ss_controller_B.Sum7_eg *
      my_pv_system_3_ss_controller_B.UnitDelay;

    /* Sum: '<S23>/Sum5' */
    my_pv_system_3_ss_controller_B.Sum5_m =
      my_pv_system_3_ss_controller_B.Meanvalue_c +
      my_pv_system_3_ss_controller_B.Product4_c;
    my_pv_system_3_ss_controller_B.Switch_g =
      my_pv_system_3_ss_controller_B.Sum5_m;
  } else {
    my_pv_system_3_ss_controller_B.Switch_g =
      my_pv_system_3_ss_controller_B.UnitDelay1;
  }

  /* End of Switch: '<S23>/Switch' */

  /* Trigonometry: '<S16>/Trigonometric Function3' */
  my_pv_system_3_ss_controller_B.TrigonometricFunction3 = cos
    (my_pv_system_3_ss_controller_B.MathFunction_i);

  /* Gain: '<S16>/Gain3' */
  my_pv_system_3_ss_controller_B.Gain3 =
    my_pv_system_3_ss_controller_P.Gain3_Gain_i *
    my_pv_system_3_ss_controller_B.TrigonometricFunction3;

  /* Product: '<S16>/Product2' */
  my_pv_system_3_ss_controller_B.Product2 = my_pv_system_3_ss_controller_B.Vpu *
    my_pv_system_3_ss_controller_B.Gain3;

  /* DiscreteIntegrator: '<S26>/Integ4' */
  if (my_pv_system_3_ss_controller_DW.Integ4_SYSTEM_ENABLE_e != 0) {
    my_pv_system_3_ss_controller_B.Integ4_h =
      my_pv_system_3_ss_controller_DW.Integ4_DSTATE_d;
  } else {
    my_pv_system_3_ss_controller_B.Integ4_h =
      my_pv_system_3_ss_controller_P.Integ4_gainval_f *
      my_pv_system_3_ss_controller_B.Product2 +
      my_pv_system_3_ss_controller_DW.Integ4_DSTATE_d;
  }

  /* End of DiscreteIntegrator: '<S26>/Integ4' */

  /* Saturate: '<S26>/To avoid division  by zero' */
  u0 = my_pv_system_3_ss_controller_B.UnitDelay;
  dP = my_pv_system_3_ss_controller_P.Toavoiddivisionbyzero_LowerSa_b;
  u2 = my_pv_system_3_ss_controller_P.Toavoiddivisionbyzero_UpperSa_c;
  if (u0 > u2) {
    my_pv_system_3_ss_controller_B.Freq_m = u2;
  } else if (u0 < dP) {
    my_pv_system_3_ss_controller_B.Freq_m = dP;
  } else {
    my_pv_system_3_ss_controller_B.Freq_m = u0;
  }

  /* End of Saturate: '<S26>/To avoid division  by zero' */

  /* Fcn: '<S26>/Number of samples per cycle' */
  my_pv_system_3_ss_controller_B.Numberofsamplespercycle_m = 1.0 /
    my_pv_system_3_ss_controller_B.Freq_m / 5.0e-5;

  /* Rounding: '<S26>/Rounding Function' */
  my_pv_system_3_ss_controller_B.RoundingFunction_p = ceil
    (my_pv_system_3_ss_controller_B.Numberofsamplespercycle_m);

  /* Gain: '<S26>/Gain' */
  my_pv_system_3_ss_controller_B.Delay_d =
    my_pv_system_3_ss_controller_P.InverterControl_Ts_Control *
    my_pv_system_3_ss_controller_B.RoundingFunction_p;

  /* Level2 S-Function Block: '<S28>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[8];
    sfcnOutputs(rts, 0);
  }

  /* UnitDelay: '<S27>/Unit Delay' */
  my_pv_system_3_ss_controller_B.UnitDelay_l =
    my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_i;

  /* DigitalClock: '<S26>/Digital  Clock' */
  my_pv_system_3_ss_controller_B.DigitalClock_m =
    my_pv_system_3_ss_controller_M->Timing.t[0];

  /* RelationalOperator: '<S26>/Relational Operator' incorporates:
   *  Constant: '<S26>/Constant'
   */
  my_pv_system_3_ss_controller_B.RelationalOperator_m =
    (my_pv_system_3_ss_controller_B.DigitalClock_m >=
     my_pv_system_3_ss_controller_P.Constant_Value_m);

  /* UnitDelay: '<S26>/Unit Delay1' */
  my_pv_system_3_ss_controller_B.UnitDelay1_e =
    my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_b;

  /* Switch: '<S26>/Switch' */
  if (my_pv_system_3_ss_controller_B.RelationalOperator_m) {
    /* Sum: '<S27>/Sum1' */
    my_pv_system_3_ss_controller_B.Sum1_cr =
      my_pv_system_3_ss_controller_B.Product2 -
      my_pv_system_3_ss_controller_B.UnitDelay_l;

    /* Sum: '<S27>/Sum5' */
    my_pv_system_3_ss_controller_B.Sum5_g =
      my_pv_system_3_ss_controller_B.Numberofsamplespercycle_m -
      my_pv_system_3_ss_controller_B.RoundingFunction_p;

    /* Product: '<S27>/Product5' */
    my_pv_system_3_ss_controller_B.Product5_a =
      my_pv_system_3_ss_controller_B.Sum5_g *
      my_pv_system_3_ss_controller_B.Sum1_cr;

    /* Gain: '<S27>/Gain1' */
    my_pv_system_3_ss_controller_B.Gain1_if =
      my_pv_system_3_ss_controller_P.Gain1_Gain_h *
      my_pv_system_3_ss_controller_B.Product5_a;

    /* Sum: '<S27>/Sum4' */
    my_pv_system_3_ss_controller_B.Sum4_b =
      my_pv_system_3_ss_controller_B.Gain1_if +
      my_pv_system_3_ss_controller_B.Product2;

    /* Product: '<S27>/Product2' */
    my_pv_system_3_ss_controller_B.Product2_dp =
      my_pv_system_3_ss_controller_B.Sum5_g /
      my_pv_system_3_ss_controller_B.Numberofsamplespercycle_m;

    /* Product: '<S27>/Product4' */
    my_pv_system_3_ss_controller_B.Product4_g =
      my_pv_system_3_ss_controller_B.Product2_dp *
      my_pv_system_3_ss_controller_B.Sum4_b;

    /* Sum: '<S26>/Sum7' */
    my_pv_system_3_ss_controller_B.Sum7_i1 =
      my_pv_system_3_ss_controller_B.Integ4_h -
      my_pv_system_3_ss_controller_B.SFunction_i;

    /* Product: '<S26>/Product' */
    my_pv_system_3_ss_controller_B.Meanvalue_e =
      my_pv_system_3_ss_controller_B.Sum7_i1 *
      my_pv_system_3_ss_controller_B.UnitDelay;

    /* Sum: '<S26>/Sum5' */
    my_pv_system_3_ss_controller_B.Sum5_g2 =
      my_pv_system_3_ss_controller_B.Meanvalue_e +
      my_pv_system_3_ss_controller_B.Product4_g;
    my_pv_system_3_ss_controller_B.Switch_k =
      my_pv_system_3_ss_controller_B.Sum5_g2;
  } else {
    my_pv_system_3_ss_controller_B.Switch_k =
      my_pv_system_3_ss_controller_B.UnitDelay1_e;
  }

  /* End of Switch: '<S26>/Switch' */

  /* RealImagToComplex: '<S16>/Real-Imag to Complex' */
  my_pv_system_3_ss_controller_B.RealImagtoComplex.re =
    my_pv_system_3_ss_controller_B.Switch_g;
  my_pv_system_3_ss_controller_B.RealImagtoComplex.im =
    my_pv_system_3_ss_controller_B.Switch_k;

  /* ComplexToMagnitudeAngle: '<S16>/Complex to Magnitude-Angle' */
  my_pv_system_3_ss_controller_B.ComplextoMagnitudeAngle_o1 = rt_hypotd_snf
    (my_pv_system_3_ss_controller_B.RealImagtoComplex.re,
     my_pv_system_3_ss_controller_B.RealImagtoComplex.im);
  my_pv_system_3_ss_controller_B.ComplextoMagnitudeAngle_o2 = rt_atan2d_snf
    (my_pv_system_3_ss_controller_B.RealImagtoComplex.im,
     my_pv_system_3_ss_controller_B.RealImagtoComplex.re);

  /* Gain: '<S16>/Rad->Deg.' */
  my_pv_system_3_ss_controller_B.RadDeg =
    my_pv_system_3_ss_controller_P.RadDeg_Gain_f *
    my_pv_system_3_ss_controller_B.ComplextoMagnitudeAngle_o2;

  /* Gain: '<S11>/to-rad' */
  my_pv_system_3_ss_controller_B.torad =
    my_pv_system_3_ss_controller_P.torad_Gain *
    my_pv_system_3_ss_controller_B.RadDeg;

  /* MagnitudeAngleToComplex: '<S11>/Magnitude-Angle to Complex' */
  P = my_pv_system_3_ss_controller_B.torad;
  dV = my_pv_system_3_ss_controller_B.ComplextoMagnitudeAngle_o1;
  u2 = cos(P);
  u0 = sin(P);
  my_pv_system_3_ss_controller_B.MagnitudeAngletoComplex.re = dV * u2;
  my_pv_system_3_ss_controller_B.MagnitudeAngletoComplex.im = dV * u0;

  /* ComplexToRealImag: '<S11>/Complex to Real-Imag' */
  my_pv_system_3_ss_controller_B.ComplextoRealImag_o1 =
    my_pv_system_3_ss_controller_B.MagnitudeAngletoComplex.re;
  my_pv_system_3_ss_controller_B.ComplextoRealImag_o2 =
    my_pv_system_3_ss_controller_B.MagnitudeAngletoComplex.im;

  /* Gain: '<S8>/Rff ' */
  my_pv_system_3_ss_controller_B.Rff = my_pv_system_3_ss_controller_P.Rff_Gain *
    my_pv_system_3_ss_controller_B.UnitDelay2;

  /* Gain: '<S8>/Lff  ' incorporates:
   *  Constant: '<S4>/Iq_ref'
   */
  my_pv_system_3_ss_controller_B.Lff = my_pv_system_3_ss_controller_P.Lff_Gain *
    my_pv_system_3_ss_controller_P.Iq_ref_Value;

  /* Sum: '<S8>/Add1' */
  my_pv_system_3_ss_controller_B.Feedforward =
    (my_pv_system_3_ss_controller_B.ComplextoRealImag_o1 +
     my_pv_system_3_ss_controller_B.Rff) - my_pv_system_3_ss_controller_B.Lff;

  /* Gain: '<S8>/Rff' incorporates:
   *  Constant: '<S4>/Iq_ref'
   */
  my_pv_system_3_ss_controller_B.Rff_e =
    my_pv_system_3_ss_controller_P.Rff_Gain_a *
    my_pv_system_3_ss_controller_P.Iq_ref_Value;

  /* Gain: '<S8>/Lff' */
  my_pv_system_3_ss_controller_B.Lff_i =
    my_pv_system_3_ss_controller_P.Lff_Gain_c *
    my_pv_system_3_ss_controller_B.UnitDelay2;

  /* Sum: '<S8>/Add3' */
  my_pv_system_3_ss_controller_B.Add3 =
    (my_pv_system_3_ss_controller_B.ComplextoRealImag_o2 +
     my_pv_system_3_ss_controller_B.Rff_e) +
    my_pv_system_3_ss_controller_B.Lff_i;

  /* Sum: '<S8>/Add2' */
  my_pv_system_3_ss_controller_B.Add2[0] =
    my_pv_system_3_ss_controller_B.Saturate[0] +
    my_pv_system_3_ss_controller_B.Feedforward;
  my_pv_system_3_ss_controller_B.Add2[1] =
    my_pv_system_3_ss_controller_B.Saturate[1] +
    my_pv_system_3_ss_controller_B.Add3;

  /* Gain: '<S15>/Integral Gain' */
  my_pv_system_3_ss_controller_B.IntegralGain[0] =
    my_pv_system_3_ss_controller_P.InverterControl_Ki_Ireg *
    my_pv_system_3_ss_controller_B.Sum_e[0];
  my_pv_system_3_ss_controller_B.IntegralGain[1] =
    my_pv_system_3_ss_controller_P.InverterControl_Ki_Ireg *
    my_pv_system_3_ss_controller_B.Sum_e[1];

  /* Saturate: '<S8>/Saturation' */
  u0 = my_pv_system_3_ss_controller_B.Add2[0];
  dP = my_pv_system_3_ss_controller_P.Saturation_LowerSat_k;
  u2 = my_pv_system_3_ss_controller_P.Saturation_UpperSat_c;
  if (u0 > u2) {
    u0 = u2;
  } else {
    if (u0 < dP) {
      u0 = dP;
    }
  }

  my_pv_system_3_ss_controller_B.Saturation[0] = u0;
  u0 = my_pv_system_3_ss_controller_B.Add2[1];
  dP = my_pv_system_3_ss_controller_P.Saturation_LowerSat_k;
  u2 = my_pv_system_3_ss_controller_P.Saturation_UpperSat_c;
  if (u0 > u2) {
    u0 = u2;
  } else {
    if (u0 < dP) {
      u0 = dP;
    }
  }

  my_pv_system_3_ss_controller_B.Saturation[1] = u0;

  /* End of Saturate: '<S8>/Saturation' */

  /* RateTransition: '<S2>/Rate Transition3' */
  my_pv_system_3_ss_controller_B.RateTransition3 =
    my_pv_system_3_ss_controller_B.SFunction_g[2];

  /* RateTransition: '<S2>/Rate Transition4' */
  my_pv_system_3_ss_controller_B.RateTransition4 =
    my_pv_system_3_ss_controller_B.SFunction_g[1];

  /* SignalConversion: '<S10>/TmpSignal ConversionAt SFunction Inport1' incorporates:
   *  Constant: '<S9>/Iph_'
   *  Constant: '<S9>/Iph_1'
   *  Constant: '<S9>/Iph_2'
   *  Constant: '<S9>/Iph_3'
   *  MATLAB Function: '<S4>/MPPT Controller using Perturbe  & Observe technique  '
   */
  my_pv_system_3_ss_controller_B.TmpSignalConversionAtSFunctionI[0] =
    my_pv_system_3_ss_controller_P.InverterControl_Vdc_ref_Init;
  my_pv_system_3_ss_controller_B.TmpSignalConversionAtSFunctionI[1] =
    my_pv_system_3_ss_controller_P.Iph_1_Value;
  my_pv_system_3_ss_controller_B.TmpSignalConversionAtSFunctionI[2] =
    my_pv_system_3_ss_controller_P.Iph_2_Value;
  my_pv_system_3_ss_controller_B.TmpSignalConversionAtSFunctionI[3] =
    my_pv_system_3_ss_controller_P.InverterControl_Increment_MPPT;

  /* MATLAB Function: '<S4>/MPPT Controller using Perturbe  & Observe technique  ' incorporates:
   *  Constant: '<S2>/MPPT_On'
   */
  /* MATLAB Function 'SS_controller/Inverter Control/MPPT Controller using Perturbe  & Observe technique  ': '<S10>:1' */
  /*  MPPT controller based on the Perturb & Observe algorithm. */
  /*  D output = Reference for DC link voltage (Vdc_ref)  */
  /*  */
  /*  Enabled input = 1 to enable the MPPT controller */
  /*  V input = PV array terminal voltage (V) */
  /*  I input = PV array current (A) */
  /*  */
  /*  Param input: */
  /* Initial value for Vdc_ref */
  /* '<S10>:1:13' */
  /* Maximum value for Vdc_ref */
  /* '<S10>:1:14' */
  /* Minimum value for Vdc_ref */
  /* '<S10>:1:15' */
  /* Increment value used to increase/decrease Vdc_ref */
  /*    */
  if (!my_pv_system_3_ss_controller_DW.Vold_not_empty) {
    /* '<S10>:1:22' */
    my_pv_system_3_ss_controller_DW.Vold_not_empty = true;

    /* '<S10>:1:25' */
    my_pv_system_3_ss_controller_DW.Dold =
      my_pv_system_3_ss_controller_B.TmpSignalConversionAtSFunctionI[0];
  }

  /* '<S10>:1:27' */
  P = my_pv_system_3_ss_controller_B.RateTransition3 *
    my_pv_system_3_ss_controller_B.RateTransition4;

  /* '<S10>:1:28' */
  dV = my_pv_system_3_ss_controller_B.RateTransition3 -
    my_pv_system_3_ss_controller_DW.Vold;

  /* '<S10>:1:29' */
  dP = P - my_pv_system_3_ss_controller_DW.Pold;
  if ((dP != 0.0) && (my_pv_system_3_ss_controller_P.MPPT_On_Value != 0.0)) {
    /* '<S10>:1:31' */
    if (dP < 0.0) {
      /* '<S10>:1:32' */
      if (dV < 0.0) {
        /* '<S10>:1:33' */
        /* '<S10>:1:34' */
        my_pv_system_3_ss_controller_B.D = my_pv_system_3_ss_controller_DW.Dold
          + my_pv_system_3_ss_controller_B.TmpSignalConversionAtSFunctionI[3];
      } else {
        /* '<S10>:1:36' */
        my_pv_system_3_ss_controller_B.D = my_pv_system_3_ss_controller_DW.Dold
          - my_pv_system_3_ss_controller_B.TmpSignalConversionAtSFunctionI[3];
      }
    } else if (dV < 0.0) {
      /* '<S10>:1:39' */
      /* '<S10>:1:40' */
      my_pv_system_3_ss_controller_B.D = my_pv_system_3_ss_controller_DW.Dold -
        my_pv_system_3_ss_controller_B.TmpSignalConversionAtSFunctionI[3];
    } else {
      /* '<S10>:1:42' */
      my_pv_system_3_ss_controller_B.D = my_pv_system_3_ss_controller_DW.Dold +
        my_pv_system_3_ss_controller_B.TmpSignalConversionAtSFunctionI[3];
    }
  } else {
    /* '<S10>:1:45' */
    my_pv_system_3_ss_controller_B.D = my_pv_system_3_ss_controller_DW.Dold;
  }

  if ((my_pv_system_3_ss_controller_B.D >=
       my_pv_system_3_ss_controller_B.TmpSignalConversionAtSFunctionI[1]) ||
      (my_pv_system_3_ss_controller_B.D <=
       my_pv_system_3_ss_controller_B.TmpSignalConversionAtSFunctionI[2])) {
    /* '<S10>:1:48' */
    /* '<S10>:1:49' */
    my_pv_system_3_ss_controller_B.D = my_pv_system_3_ss_controller_DW.Dold;
  }

  /* '<S10>:1:52' */
  my_pv_system_3_ss_controller_DW.Dold = my_pv_system_3_ss_controller_B.D;

  /* '<S10>:1:53' */
  my_pv_system_3_ss_controller_DW.Vold =
    my_pv_system_3_ss_controller_B.RateTransition3;

  /* '<S10>:1:54' */
  my_pv_system_3_ss_controller_DW.Pold = P;

  /* DigitalClock: '<S29>/Digital  Clock' */
  my_pv_system_3_ss_controller_B.DigitalClock_a =
    my_pv_system_3_ss_controller_M->Timing.t[0];

  /* RateTransition: '<S2>/Rate Transition2' */
  my_pv_system_3_ss_controller_B.RateTransition2 =
    my_pv_system_3_ss_controller_B.SFunction_g[7];

  /* DiscreteIntegrator: '<S29>/Integ4' */
  if (my_pv_system_3_ss_controller_DW.Integ4_SYSTEM_ENABLE_d != 0) {
    my_pv_system_3_ss_controller_B.Integ4_c =
      my_pv_system_3_ss_controller_DW.Integ4_DSTATE_d5;
  } else {
    my_pv_system_3_ss_controller_B.Integ4_c =
      my_pv_system_3_ss_controller_P.Integ4_gainval_j *
      my_pv_system_3_ss_controller_B.RateTransition2 +
      my_pv_system_3_ss_controller_DW.Integ4_DSTATE_d5;
  }

  /* End of DiscreteIntegrator: '<S29>/Integ4' */

  /* Constant: '<S29>/K1' */
  my_pv_system_3_ss_controller_B.K1 = my_pv_system_3_ss_controller_P.K1_Value;

  /* Level2 S-Function Block: '<S30>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[9];
    sfcnOutputs(rts, 0);
  }

  /* UnitDelay: '<S29>/Unit Delay' */
  my_pv_system_3_ss_controller_B.UnitDelay_h =
    my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_iq;

  /* RelationalOperator: '<S29>/Relational Operator' */
  my_pv_system_3_ss_controller_B.RelationalOperator_j =
    (my_pv_system_3_ss_controller_B.DigitalClock_a >=
     my_pv_system_3_ss_controller_B.K1);

  /* UnitDelay: '<S29>/Unit Delay1' */
  my_pv_system_3_ss_controller_B.UnitDelay1_p =
    my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_a;

  /* Switch: '<S29>/Switch' */
  if (my_pv_system_3_ss_controller_B.RelationalOperator_j) {
    /* Gain: '<S29>/Gain1' */
    my_pv_system_3_ss_controller_B.Gain1_c =
      my_pv_system_3_ss_controller_P.Gain1_Gain_o *
      my_pv_system_3_ss_controller_B.RateTransition2;

    /* Gain: '<S29>/Gain' */
    my_pv_system_3_ss_controller_B.Gain_f =
      my_pv_system_3_ss_controller_P.Gain_Gain *
      my_pv_system_3_ss_controller_B.UnitDelay_h;

    /* Sum: '<S29>/Sum1' */
    my_pv_system_3_ss_controller_B.Correction =
      my_pv_system_3_ss_controller_B.Gain1_c -
      my_pv_system_3_ss_controller_B.Gain_f;

    /* Sum: '<S29>/Sum7' */
    my_pv_system_3_ss_controller_B.Sum7_i =
      my_pv_system_3_ss_controller_B.Integ4_c -
      my_pv_system_3_ss_controller_B.SFunction_pd;

    /* Product: '<S29>/Product' incorporates:
     *  Constant: '<S29>/K2'
     */
    my_pv_system_3_ss_controller_B.Mean = my_pv_system_3_ss_controller_B.Sum7_i *
      my_pv_system_3_ss_controller_P.K2_Value;

    /* Sum: '<S29>/Sum5' */
    my_pv_system_3_ss_controller_B.Sum5_c = my_pv_system_3_ss_controller_B.Mean
      + my_pv_system_3_ss_controller_B.Correction;
    my_pv_system_3_ss_controller_B.Switch_c =
      my_pv_system_3_ss_controller_B.Sum5_c;
  } else {
    my_pv_system_3_ss_controller_B.Switch_c =
      my_pv_system_3_ss_controller_B.UnitDelay1_p;
  }

  /* End of Switch: '<S29>/Switch' */

  /* Outputs for Enabled SubSystem: '<S31>/Automatic Gain Control' incorporates:
   *  EnablePort: '<S32>/Enable'
   */
  /* Constant: '<S31>/Constant1' */
  if (my_pv_system_3_ss_controller_P.PLL_AGC > 0.0) {
    if (!my_pv_system_3_ss_controller_DW.AutomaticGainControl_MODE) {
      /* Enable for DiscreteIntegrator: '<S39>/Integ4' */
      my_pv_system_3_ss_controller_DW.Integ4_SYSTEM_ENABLE_g = 1U;

      /* Enable for DiscreteIntegrator: '<S42>/Integ4' */
      my_pv_system_3_ss_controller_DW.Integ4_SYSTEM_ENABLE_gx = 1U;
      my_pv_system_3_ss_controller_DW.AutomaticGainControl_MODE = true;
    }

    /* Trigonometry: '<S36>/Trigonometric Function' */
    my_pv_system_3_ss_controller_B.TrigonometricFunction_p = sin
      (my_pv_system_3_ss_controller_B.MathFunction_i);

    /* Gain: '<S36>/Gain1' */
    my_pv_system_3_ss_controller_B.Gain1_i =
      my_pv_system_3_ss_controller_P.Gain1_Gain_e *
      my_pv_system_3_ss_controller_B.TrigonometricFunction_p;

    /* Product: '<S36>/Product1' */
    my_pv_system_3_ss_controller_B.Product1_o =
      my_pv_system_3_ss_controller_B.Vpu *
      my_pv_system_3_ss_controller_B.Gain1_i;

    /* DiscreteIntegrator: '<S39>/Integ4' */
    if (my_pv_system_3_ss_controller_DW.Integ4_SYSTEM_ENABLE_g != 0) {
      my_pv_system_3_ss_controller_B.Integ4_d =
        my_pv_system_3_ss_controller_DW.Integ4_DSTATE_ky;
    } else {
      my_pv_system_3_ss_controller_B.Integ4_d =
        my_pv_system_3_ss_controller_P.Integ4_gainval *
        my_pv_system_3_ss_controller_B.Product1_o +
        my_pv_system_3_ss_controller_DW.Integ4_DSTATE_ky;
    }

    /* End of DiscreteIntegrator: '<S39>/Integ4' */

    /* Saturate: '<S39>/To avoid division  by zero' */
    u0 = my_pv_system_3_ss_controller_B.UnitDelay;
    dP = my_pv_system_3_ss_controller_P.Toavoiddivisionbyzero_LowerSat;
    u2 = my_pv_system_3_ss_controller_P.Toavoiddivisionbyzero_UpperSat;
    if (u0 > u2) {
      my_pv_system_3_ss_controller_B.Freq_h = u2;
    } else if (u0 < dP) {
      my_pv_system_3_ss_controller_B.Freq_h = dP;
    } else {
      my_pv_system_3_ss_controller_B.Freq_h = u0;
    }

    /* End of Saturate: '<S39>/To avoid division  by zero' */

    /* Fcn: '<S39>/Number of samples per cycle' */
    my_pv_system_3_ss_controller_B.Numberofsamplespercycle_h = 1.0 /
      my_pv_system_3_ss_controller_B.Freq_h / 5.0e-5;

    /* Rounding: '<S39>/Rounding Function' */
    my_pv_system_3_ss_controller_B.RoundingFunction_e = ceil
      (my_pv_system_3_ss_controller_B.Numberofsamplespercycle_h);

    /* Gain: '<S39>/Gain' */
    my_pv_system_3_ss_controller_B.Delay_db =
      my_pv_system_3_ss_controller_P.InverterControl_Ts_Control *
      my_pv_system_3_ss_controller_B.RoundingFunction_e;

    /* Level2 S-Function Block: '<S41>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[0];
      sfcnOutputs(rts, 0);
    }

    /* UnitDelay: '<S40>/Unit Delay' */
    my_pv_system_3_ss_controller_B.UnitDelay_b =
      my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_l;

    /* DigitalClock: '<S39>/Digital  Clock' */
    my_pv_system_3_ss_controller_B.DigitalClock_o =
      my_pv_system_3_ss_controller_M->Timing.t[0];

    /* RelationalOperator: '<S39>/Relational Operator' incorporates:
     *  Constant: '<S39>/Constant'
     */
    my_pv_system_3_ss_controller_B.RelationalOperator_d =
      (my_pv_system_3_ss_controller_B.DigitalClock_o >=
       my_pv_system_3_ss_controller_P.Constant_Value);

    /* UnitDelay: '<S39>/Unit Delay1' */
    my_pv_system_3_ss_controller_B.UnitDelay1_c2 =
      my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_o;

    /* Switch: '<S39>/Switch' */
    if (my_pv_system_3_ss_controller_B.RelationalOperator_d) {
      /* Sum: '<S40>/Sum1' */
      my_pv_system_3_ss_controller_B.Sum1_p =
        my_pv_system_3_ss_controller_B.Product1_o -
        my_pv_system_3_ss_controller_B.UnitDelay_b;

      /* Sum: '<S40>/Sum5' */
      my_pv_system_3_ss_controller_B.Sum5_k =
        my_pv_system_3_ss_controller_B.Numberofsamplespercycle_h -
        my_pv_system_3_ss_controller_B.RoundingFunction_e;

      /* Product: '<S40>/Product5' */
      my_pv_system_3_ss_controller_B.Product5_n =
        my_pv_system_3_ss_controller_B.Sum5_k *
        my_pv_system_3_ss_controller_B.Sum1_p;

      /* Gain: '<S40>/Gain1' */
      my_pv_system_3_ss_controller_B.Gain1_o =
        my_pv_system_3_ss_controller_P.Gain1_Gain_f *
        my_pv_system_3_ss_controller_B.Product5_n;

      /* Sum: '<S40>/Sum4' */
      my_pv_system_3_ss_controller_B.Sum4_n =
        my_pv_system_3_ss_controller_B.Gain1_o +
        my_pv_system_3_ss_controller_B.Product1_o;

      /* Product: '<S40>/Product2' */
      my_pv_system_3_ss_controller_B.Product2_d =
        my_pv_system_3_ss_controller_B.Sum5_k /
        my_pv_system_3_ss_controller_B.Numberofsamplespercycle_h;

      /* Product: '<S40>/Product4' */
      my_pv_system_3_ss_controller_B.Product4_o =
        my_pv_system_3_ss_controller_B.Product2_d *
        my_pv_system_3_ss_controller_B.Sum4_n;

      /* Sum: '<S39>/Sum7' */
      my_pv_system_3_ss_controller_B.Sum7_e =
        my_pv_system_3_ss_controller_B.Integ4_d -
        my_pv_system_3_ss_controller_B.SFunction_ij;

      /* Product: '<S39>/Product' */
      my_pv_system_3_ss_controller_B.Meanvalue_h =
        my_pv_system_3_ss_controller_B.Sum7_e *
        my_pv_system_3_ss_controller_B.UnitDelay;

      /* Sum: '<S39>/Sum5' */
      my_pv_system_3_ss_controller_B.Sum5_f =
        my_pv_system_3_ss_controller_B.Meanvalue_h +
        my_pv_system_3_ss_controller_B.Product4_o;
      my_pv_system_3_ss_controller_B.Switch_j =
        my_pv_system_3_ss_controller_B.Sum5_f;
    } else {
      my_pv_system_3_ss_controller_B.Switch_j =
        my_pv_system_3_ss_controller_B.UnitDelay1_c2;
    }

    /* End of Switch: '<S39>/Switch' */

    /* Trigonometry: '<S36>/Trigonometric Function3' */
    my_pv_system_3_ss_controller_B.TrigonometricFunction3_p = cos
      (my_pv_system_3_ss_controller_B.MathFunction_i);

    /* Gain: '<S36>/Gain3' */
    my_pv_system_3_ss_controller_B.Gain3_c =
      my_pv_system_3_ss_controller_P.Gain3_Gain *
      my_pv_system_3_ss_controller_B.TrigonometricFunction3_p;

    /* Product: '<S36>/Product2' */
    my_pv_system_3_ss_controller_B.Product2_i =
      my_pv_system_3_ss_controller_B.Vpu *
      my_pv_system_3_ss_controller_B.Gain3_c;

    /* DiscreteIntegrator: '<S42>/Integ4' */
    if (my_pv_system_3_ss_controller_DW.Integ4_SYSTEM_ENABLE_gx != 0) {
      my_pv_system_3_ss_controller_B.Integ4_dy =
        my_pv_system_3_ss_controller_DW.Integ4_DSTATE_j;
    } else {
      my_pv_system_3_ss_controller_B.Integ4_dy =
        my_pv_system_3_ss_controller_P.Integ4_gainval_h *
        my_pv_system_3_ss_controller_B.Product2_i +
        my_pv_system_3_ss_controller_DW.Integ4_DSTATE_j;
    }

    /* End of DiscreteIntegrator: '<S42>/Integ4' */

    /* Saturate: '<S42>/To avoid division  by zero' */
    u0 = my_pv_system_3_ss_controller_B.UnitDelay;
    dP = my_pv_system_3_ss_controller_P.Toavoiddivisionbyzero_LowerSa_m;
    u2 = my_pv_system_3_ss_controller_P.Toavoiddivisionbyzero_UpperSa_n;
    if (u0 > u2) {
      my_pv_system_3_ss_controller_B.Freq_mz = u2;
    } else if (u0 < dP) {
      my_pv_system_3_ss_controller_B.Freq_mz = dP;
    } else {
      my_pv_system_3_ss_controller_B.Freq_mz = u0;
    }

    /* End of Saturate: '<S42>/To avoid division  by zero' */

    /* Fcn: '<S42>/Number of samples per cycle' */
    my_pv_system_3_ss_controller_B.Numberofsamplespercycle_d = 1.0 /
      my_pv_system_3_ss_controller_B.Freq_mz / 5.0e-5;

    /* Rounding: '<S42>/Rounding Function' */
    my_pv_system_3_ss_controller_B.RoundingFunction_h = ceil
      (my_pv_system_3_ss_controller_B.Numberofsamplespercycle_d);

    /* Gain: '<S42>/Gain' */
    my_pv_system_3_ss_controller_B.Delay_i =
      my_pv_system_3_ss_controller_P.InverterControl_Ts_Control *
      my_pv_system_3_ss_controller_B.RoundingFunction_h;

    /* Level2 S-Function Block: '<S44>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[1];
      sfcnOutputs(rts, 0);
    }

    /* UnitDelay: '<S43>/Unit Delay' */
    my_pv_system_3_ss_controller_B.UnitDelay_fb =
      my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_d;

    /* DigitalClock: '<S42>/Digital  Clock' */
    my_pv_system_3_ss_controller_B.DigitalClock_l =
      my_pv_system_3_ss_controller_M->Timing.t[0];

    /* RelationalOperator: '<S42>/Relational Operator' incorporates:
     *  Constant: '<S42>/Constant'
     */
    my_pv_system_3_ss_controller_B.RelationalOperator_k =
      (my_pv_system_3_ss_controller_B.DigitalClock_l >=
       my_pv_system_3_ss_controller_P.Constant_Value_c);

    /* UnitDelay: '<S42>/Unit Delay1' */
    my_pv_system_3_ss_controller_B.UnitDelay1_j =
      my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_k;

    /* Switch: '<S42>/Switch' */
    if (my_pv_system_3_ss_controller_B.RelationalOperator_k) {
      /* Sum: '<S43>/Sum1' */
      my_pv_system_3_ss_controller_B.Sum1_c =
        my_pv_system_3_ss_controller_B.Product2_i -
        my_pv_system_3_ss_controller_B.UnitDelay_fb;

      /* Sum: '<S43>/Sum5' */
      my_pv_system_3_ss_controller_B.Sum5_o =
        my_pv_system_3_ss_controller_B.Numberofsamplespercycle_d -
        my_pv_system_3_ss_controller_B.RoundingFunction_h;

      /* Product: '<S43>/Product5' */
      my_pv_system_3_ss_controller_B.Product5_e =
        my_pv_system_3_ss_controller_B.Sum5_o *
        my_pv_system_3_ss_controller_B.Sum1_c;

      /* Gain: '<S43>/Gain1' */
      my_pv_system_3_ss_controller_B.Gain1_n =
        my_pv_system_3_ss_controller_P.Gain1_Gain_d *
        my_pv_system_3_ss_controller_B.Product5_e;

      /* Sum: '<S43>/Sum4' */
      my_pv_system_3_ss_controller_B.Sum4_g =
        my_pv_system_3_ss_controller_B.Gain1_n +
        my_pv_system_3_ss_controller_B.Product2_i;

      /* Product: '<S43>/Product2' */
      my_pv_system_3_ss_controller_B.Product2_m =
        my_pv_system_3_ss_controller_B.Sum5_o /
        my_pv_system_3_ss_controller_B.Numberofsamplespercycle_d;

      /* Product: '<S43>/Product4' */
      my_pv_system_3_ss_controller_B.Product4_d =
        my_pv_system_3_ss_controller_B.Product2_m *
        my_pv_system_3_ss_controller_B.Sum4_g;

      /* Sum: '<S42>/Sum7' */
      my_pv_system_3_ss_controller_B.Sum7_m =
        my_pv_system_3_ss_controller_B.Integ4_dy -
        my_pv_system_3_ss_controller_B.SFunction_ed;

      /* Product: '<S42>/Product' */
      my_pv_system_3_ss_controller_B.Meanvalue_g =
        my_pv_system_3_ss_controller_B.Sum7_m *
        my_pv_system_3_ss_controller_B.UnitDelay;

      /* Sum: '<S42>/Sum5' */
      my_pv_system_3_ss_controller_B.Sum5_i =
        my_pv_system_3_ss_controller_B.Meanvalue_g +
        my_pv_system_3_ss_controller_B.Product4_d;
      my_pv_system_3_ss_controller_B.Switch_m =
        my_pv_system_3_ss_controller_B.Sum5_i;
    } else {
      my_pv_system_3_ss_controller_B.Switch_m =
        my_pv_system_3_ss_controller_B.UnitDelay1_j;
    }

    /* End of Switch: '<S42>/Switch' */

    /* RealImagToComplex: '<S36>/Real-Imag to Complex' */
    my_pv_system_3_ss_controller_B.RealImagtoComplex_l.re =
      my_pv_system_3_ss_controller_B.Switch_j;
    my_pv_system_3_ss_controller_B.RealImagtoComplex_l.im =
      my_pv_system_3_ss_controller_B.Switch_m;

    /* ComplexToMagnitudeAngle: '<S36>/Complex to Magnitude-Angle' */
    my_pv_system_3_ss_controller_B.ComplextoMagnitudeAngle_o1_d = rt_hypotd_snf
      (my_pv_system_3_ss_controller_B.RealImagtoComplex_l.re,
       my_pv_system_3_ss_controller_B.RealImagtoComplex_l.im);
    my_pv_system_3_ss_controller_B.ComplextoMagnitudeAngle_o2_i = rt_atan2d_snf
      (my_pv_system_3_ss_controller_B.RealImagtoComplex_l.im,
       my_pv_system_3_ss_controller_B.RealImagtoComplex_l.re);

    /* Gain: '<S36>/Rad->Deg.' */
    my_pv_system_3_ss_controller_B.RadDeg_h =
      my_pv_system_3_ss_controller_P.RadDeg_Gain *
      my_pv_system_3_ss_controller_B.ComplextoMagnitudeAngle_o2_i;

    /* Saturate: '<S32>/Saturation' */
    u0 = my_pv_system_3_ss_controller_B.ComplextoMagnitudeAngle_o1_d;
    dP = my_pv_system_3_ss_controller_P.Saturation_LowerSat;
    u2 = my_pv_system_3_ss_controller_P.Saturation_UpperSat;
    if (u0 > u2) {
      my_pv_system_3_ss_controller_B.Saturation_n = u2;
    } else if (u0 < dP) {
      my_pv_system_3_ss_controller_B.Saturation_n = dP;
    } else {
      my_pv_system_3_ss_controller_B.Saturation_n = u0;
    }

    /* End of Saturate: '<S32>/Saturation' */

    /* Math: '<S32>/Math Function'
     *
     * About '<S32>/Math Function':
     *  Operator: reciprocal
     */
    P = my_pv_system_3_ss_controller_B.Saturation_n;
    my_pv_system_3_ss_controller_B.MathFunction_k = 1.0 / P;
  } else {
    if (my_pv_system_3_ss_controller_DW.AutomaticGainControl_MODE) {
      my_pv_system_3_ss_controller_DW.AutomaticGainControl_MODE = false;
    }
  }

  /* End of Constant: '<S31>/Constant1' */
  /* End of Outputs for SubSystem: '<S31>/Automatic Gain Control' */

  /* Trigonometry: '<S31>/Trigonometric Function2' */
  my_pv_system_3_ss_controller_B.TrigonometricFunction2 = cos
    (my_pv_system_3_ss_controller_B.MathFunction_i);

  /* Product: '<S31>/Product1' */
  my_pv_system_3_ss_controller_B.Product1_c = my_pv_system_3_ss_controller_B.Vpu
    * my_pv_system_3_ss_controller_B.TrigonometricFunction2;

  /* DiscreteIntegrator: '<S45>/Integ4' */
  if (my_pv_system_3_ss_controller_DW.Integ4_SYSTEM_ENABLE_k != 0) {
    my_pv_system_3_ss_controller_B.Integ4_e =
      my_pv_system_3_ss_controller_DW.Integ4_DSTATE_k;
  } else {
    my_pv_system_3_ss_controller_B.Integ4_e =
      my_pv_system_3_ss_controller_P.Integ4_gainval_jy *
      my_pv_system_3_ss_controller_B.Product1_c +
      my_pv_system_3_ss_controller_DW.Integ4_DSTATE_k;
  }

  /* End of DiscreteIntegrator: '<S45>/Integ4' */

  /* Saturate: '<S45>/To avoid division  by zero' */
  u0 = my_pv_system_3_ss_controller_B.UnitDelay;
  dP = my_pv_system_3_ss_controller_P.Toavoiddivisionbyzero_LowerSa_a;
  u2 = my_pv_system_3_ss_controller_P.Toavoiddivisionbyzero_UpperSa_e;
  if (u0 > u2) {
    my_pv_system_3_ss_controller_B.Freq_c = u2;
  } else if (u0 < dP) {
    my_pv_system_3_ss_controller_B.Freq_c = dP;
  } else {
    my_pv_system_3_ss_controller_B.Freq_c = u0;
  }

  /* End of Saturate: '<S45>/To avoid division  by zero' */

  /* Fcn: '<S45>/Number of samples per cycle' */
  my_pv_system_3_ss_controller_B.Numberofsamplespercycle_f = 1.0 /
    my_pv_system_3_ss_controller_B.Freq_c / 5.0e-5;

  /* Rounding: '<S45>/Rounding Function' */
  my_pv_system_3_ss_controller_B.RoundingFunction_b = ceil
    (my_pv_system_3_ss_controller_B.Numberofsamplespercycle_f);

  /* Gain: '<S45>/Gain' */
  my_pv_system_3_ss_controller_B.Delay_m =
    my_pv_system_3_ss_controller_P.InverterControl_Ts_Control *
    my_pv_system_3_ss_controller_B.RoundingFunction_b;

  /* Level2 S-Function Block: '<S47>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[10];
    sfcnOutputs(rts, 0);
  }

  /* UnitDelay: '<S46>/Unit Delay' */
  my_pv_system_3_ss_controller_B.UnitDelay_o =
    my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_a;

  /* DigitalClock: '<S45>/Digital  Clock' */
  my_pv_system_3_ss_controller_B.DigitalClock_e =
    my_pv_system_3_ss_controller_M->Timing.t[0];

  /* RelationalOperator: '<S45>/Relational Operator' incorporates:
   *  Constant: '<S45>/Constant'
   */
  my_pv_system_3_ss_controller_B.RelationalOperator_l =
    (my_pv_system_3_ss_controller_B.DigitalClock_e >=
     my_pv_system_3_ss_controller_P.Constant_Value_f);

  /* UnitDelay: '<S45>/Unit Delay1' */
  my_pv_system_3_ss_controller_B.UnitDelay1_a =
    my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_i;

  /* Switch: '<S45>/Switch' */
  if (my_pv_system_3_ss_controller_B.RelationalOperator_l) {
    /* Sum: '<S46>/Sum1' */
    my_pv_system_3_ss_controller_B.Sum1 =
      my_pv_system_3_ss_controller_B.Product1_c -
      my_pv_system_3_ss_controller_B.UnitDelay_o;

    /* Sum: '<S46>/Sum5' */
    my_pv_system_3_ss_controller_B.Sum5 =
      my_pv_system_3_ss_controller_B.Numberofsamplespercycle_f -
      my_pv_system_3_ss_controller_B.RoundingFunction_b;

    /* Product: '<S46>/Product5' */
    my_pv_system_3_ss_controller_B.Product5 =
      my_pv_system_3_ss_controller_B.Sum5 * my_pv_system_3_ss_controller_B.Sum1;

    /* Gain: '<S46>/Gain1' */
    my_pv_system_3_ss_controller_B.Gain1_pw =
      my_pv_system_3_ss_controller_P.Gain1_Gain_hn *
      my_pv_system_3_ss_controller_B.Product5;

    /* Sum: '<S46>/Sum4' */
    my_pv_system_3_ss_controller_B.Sum4 =
      my_pv_system_3_ss_controller_B.Gain1_pw +
      my_pv_system_3_ss_controller_B.Product1_c;

    /* Product: '<S46>/Product2' */
    my_pv_system_3_ss_controller_B.Product2_b =
      my_pv_system_3_ss_controller_B.Sum5 /
      my_pv_system_3_ss_controller_B.Numberofsamplespercycle_f;

    /* Product: '<S46>/Product4' */
    my_pv_system_3_ss_controller_B.Product4 =
      my_pv_system_3_ss_controller_B.Product2_b *
      my_pv_system_3_ss_controller_B.Sum4;

    /* Sum: '<S45>/Sum7' */
    my_pv_system_3_ss_controller_B.Sum7 =
      my_pv_system_3_ss_controller_B.Integ4_e -
      my_pv_system_3_ss_controller_B.SFunction_e;

    /* Product: '<S45>/Product' */
    my_pv_system_3_ss_controller_B.Meanvalue =
      my_pv_system_3_ss_controller_B.Sum7 *
      my_pv_system_3_ss_controller_B.UnitDelay;

    /* Sum: '<S45>/Sum5' */
    my_pv_system_3_ss_controller_B.Sum5_b =
      my_pv_system_3_ss_controller_B.Meanvalue +
      my_pv_system_3_ss_controller_B.Product4;
    my_pv_system_3_ss_controller_B.Switch_d =
      my_pv_system_3_ss_controller_B.Sum5_b;
  } else {
    my_pv_system_3_ss_controller_B.Switch_d =
      my_pv_system_3_ss_controller_B.UnitDelay1_a;
  }

  /* End of Switch: '<S45>/Switch' */

  /* Product: '<S31>/Divide' */
  my_pv_system_3_ss_controller_B.Divide =
    my_pv_system_3_ss_controller_B.Switch_d *
    my_pv_system_3_ss_controller_B.MathFunction_k;

  /* DiscreteTransferFcn: '<S33>/Discrete Derivative ' */
  P = my_pv_system_3_ss_controller_B.Divide;
  P -= my_pv_system_3_ss_controller_P.DiscreteDerivative_DenCoef[1] *
    my_pv_system_3_ss_controller_DW.DiscreteDerivative_states;
  P /= my_pv_system_3_ss_controller_P.DiscreteDerivative_DenCoef[0];
  my_pv_system_3_ss_controller_DW.DiscreteDerivative_tmp = P;
  u0 = 1.0;
  P = my_pv_system_3_ss_controller_P.Discrete_Kd;
  u0 *= P;
  dV = u0 * my_pv_system_3_ss_controller_DW.DiscreteDerivative_tmp;
  u0 = (-1.0);
  P = my_pv_system_3_ss_controller_P.Discrete_Kd;
  u0 *= P;
  dV += u0 * my_pv_system_3_ss_controller_DW.DiscreteDerivative_states;
  my_pv_system_3_ss_controller_B.DiscreteDerivative = dV;

  /* DiscreteIntegrator: '<S33>/Discrete-Time Integrator' */
  my_pv_system_3_ss_controller_B.DiscreteTimeIntegrator_g =
    my_pv_system_3_ss_controller_DW.DiscreteTimeIntegrator_DSTATE_d;

  /* Gain: '<S33>/Kp4' */
  my_pv_system_3_ss_controller_B.Kp4 =
    my_pv_system_3_ss_controller_P.Discrete_Kp *
    my_pv_system_3_ss_controller_B.Divide;

  /* Sum: '<S33>/Sum6' */
  my_pv_system_3_ss_controller_B.Sum6 = (my_pv_system_3_ss_controller_B.Kp4 +
    my_pv_system_3_ss_controller_B.DiscreteTimeIntegrator_g) +
    my_pv_system_3_ss_controller_B.DiscreteDerivative;

  /* Saturate: '<S33>/Saturation1' */
  u0 = my_pv_system_3_ss_controller_B.Sum6;
  dP = my_pv_system_3_ss_controller_P.Saturation1_LowerSat;
  u2 = my_pv_system_3_ss_controller_P.Saturation1_UpperSat;
  if (u0 > u2) {
    my_pv_system_3_ss_controller_B.Saturation1 = u2;
  } else if (u0 < dP) {
    my_pv_system_3_ss_controller_B.Saturation1 = dP;
  } else {
    my_pv_system_3_ss_controller_B.Saturation1 = u0;
  }

  /* End of Saturate: '<S33>/Saturation1' */

  /* Gain: '<S31>/Gain10' */
  my_pv_system_3_ss_controller_B.Gain10 =
    my_pv_system_3_ss_controller_P.Gain10_Gain *
    my_pv_system_3_ss_controller_B.Saturation1;

  /* RateLimiter: '<S31>/Rate Limiter' */
  P = my_pv_system_3_ss_controller_B.Gain10 -
    my_pv_system_3_ss_controller_DW.PrevY;
  if (P > my_pv_system_3_ss_controller_P.RateLimiter_RisingLim) {
    my_pv_system_3_ss_controller_B.RateLimiter =
      my_pv_system_3_ss_controller_DW.PrevY +
      my_pv_system_3_ss_controller_P.RateLimiter_RisingLim;
  } else if (P < my_pv_system_3_ss_controller_P.RateLimiter_FallingLim) {
    my_pv_system_3_ss_controller_B.RateLimiter =
      my_pv_system_3_ss_controller_DW.PrevY +
      my_pv_system_3_ss_controller_P.RateLimiter_FallingLim;
  } else {
    my_pv_system_3_ss_controller_B.RateLimiter =
      my_pv_system_3_ss_controller_B.Gain10;
  }

  my_pv_system_3_ss_controller_DW.PrevY =
    my_pv_system_3_ss_controller_B.RateLimiter;

  /* End of RateLimiter: '<S31>/Rate Limiter' */

  /* UnitDelay: '<S48>/Delay_x1' */
  my_pv_system_3_ss_controller_B.x1k_d =
    my_pv_system_3_ss_controller_DW.Delay_x1_DSTATE_d;

  /* Gain: '<S49>/A11' */
  my_pv_system_3_ss_controller_B.A11 = my_pv_system_3_ss_controller_P.A11_Gain *
    my_pv_system_3_ss_controller_B.x1k_d;

  /* UnitDelay: '<S48>/Delay_x2' */
  my_pv_system_3_ss_controller_B.x2k_k =
    my_pv_system_3_ss_controller_DW.Delay_x2_DSTATE_m;

  /* Gain: '<S49>/A12' */
  my_pv_system_3_ss_controller_B.A12 = my_pv_system_3_ss_controller_P.A12_Gain *
    my_pv_system_3_ss_controller_B.x2k_k;

  /* Gain: '<S49>/A21' */
  my_pv_system_3_ss_controller_B.A21 = my_pv_system_3_ss_controller_P.A21_Gain *
    my_pv_system_3_ss_controller_B.x1k_d;

  /* Gain: '<S49>/A22' */
  my_pv_system_3_ss_controller_B.A22 = my_pv_system_3_ss_controller_P.A22_Gain *
    my_pv_system_3_ss_controller_B.x2k_k;

  /* Sum: '<S49>/sum2' */
  my_pv_system_3_ss_controller_B.sum2_c = my_pv_system_3_ss_controller_B.A11 +
    my_pv_system_3_ss_controller_B.A12;

  /* Sum: '<S49>/sum3' */
  my_pv_system_3_ss_controller_B.sum3 = my_pv_system_3_ss_controller_B.A21 +
    my_pv_system_3_ss_controller_B.A22;

  /* Gain: '<S50>/B11' */
  my_pv_system_3_ss_controller_B.B11 = my_pv_system_3_ss_controller_P.B11_Gain *
    my_pv_system_3_ss_controller_B.RateLimiter;

  /* Sum: '<S48>/A*x1(k) + B*u1(k) ' */
  my_pv_system_3_ss_controller_B.x1k1 = my_pv_system_3_ss_controller_B.sum2_c +
    my_pv_system_3_ss_controller_B.B11;

  /* Gain: '<S50>/B21' */
  my_pv_system_3_ss_controller_B.B21 = my_pv_system_3_ss_controller_P.B21_Gain *
    my_pv_system_3_ss_controller_B.RateLimiter;

  /* Sum: '<S48>/A*x2(k) + B*u2(k)' */
  my_pv_system_3_ss_controller_B.x2k1 = my_pv_system_3_ss_controller_B.sum3 +
    my_pv_system_3_ss_controller_B.B21;

  /* Gain: '<S48>/D*u(k)' */
  my_pv_system_3_ss_controller_B.Duk_b =
    my_pv_system_3_ss_controller_P.Duk_Gain_p *
    my_pv_system_3_ss_controller_B.RateLimiter;

  /* Gain: '<S51>/C11' */
  my_pv_system_3_ss_controller_B.C11_d =
    my_pv_system_3_ss_controller_P.C11_Gain_l *
    my_pv_system_3_ss_controller_B.x1k_d;

  /* Gain: '<S51>/C12' */
  my_pv_system_3_ss_controller_B.C12_g =
    my_pv_system_3_ss_controller_P.C12_Gain_k *
    my_pv_system_3_ss_controller_B.x2k_k;

  /* Sum: '<S51>/sum2' */
  my_pv_system_3_ss_controller_B.sum2_o = my_pv_system_3_ss_controller_B.C11_d +
    my_pv_system_3_ss_controller_B.C12_g;

  /* Sum: '<S48>/C*X(k)+D*u(k)' */
  my_pv_system_3_ss_controller_B.yk_e = my_pv_system_3_ss_controller_B.Duk_b +
    my_pv_system_3_ss_controller_B.sum2_o;

  /* Gain: '<S53>/A11' */
  my_pv_system_3_ss_controller_B.A11_e[0] =
    my_pv_system_3_ss_controller_P.A11_Gain_g *
    my_pv_system_3_ss_controller_B.x1k[0];
  my_pv_system_3_ss_controller_B.A11_e[1] =
    my_pv_system_3_ss_controller_P.A11_Gain_g *
    my_pv_system_3_ss_controller_B.x1k[1];

  /* Gain: '<S53>/A12' */
  my_pv_system_3_ss_controller_B.A12_i[0] =
    my_pv_system_3_ss_controller_P.A12_Gain_l *
    my_pv_system_3_ss_controller_B.x2k[0];
  my_pv_system_3_ss_controller_B.A12_i[1] =
    my_pv_system_3_ss_controller_P.A12_Gain_l *
    my_pv_system_3_ss_controller_B.x2k[1];

  /* Gain: '<S53>/A21' */
  my_pv_system_3_ss_controller_B.A21_o[0] =
    my_pv_system_3_ss_controller_P.A21_Gain_k *
    my_pv_system_3_ss_controller_B.x1k[0];
  my_pv_system_3_ss_controller_B.A21_o[1] =
    my_pv_system_3_ss_controller_P.A21_Gain_k *
    my_pv_system_3_ss_controller_B.x1k[1];

  /* Gain: '<S53>/A22' */
  my_pv_system_3_ss_controller_B.A22_p[0] =
    my_pv_system_3_ss_controller_P.A22_Gain_p *
    my_pv_system_3_ss_controller_B.x2k[0];
  my_pv_system_3_ss_controller_B.A22_p[1] =
    my_pv_system_3_ss_controller_P.A22_Gain_p *
    my_pv_system_3_ss_controller_B.x2k[1];

  /* Sum: '<S53>/sum2' */
  my_pv_system_3_ss_controller_B.sum2_e[0] =
    my_pv_system_3_ss_controller_B.A11_e[0] +
    my_pv_system_3_ss_controller_B.A12_i[0];
  my_pv_system_3_ss_controller_B.sum2_e[1] =
    my_pv_system_3_ss_controller_B.A11_e[1] +
    my_pv_system_3_ss_controller_B.A12_i[1];

  /* Sum: '<S53>/sum3' */
  my_pv_system_3_ss_controller_B.sum3_g[0] =
    my_pv_system_3_ss_controller_B.A21_o[0] +
    my_pv_system_3_ss_controller_B.A22_p[0];
  my_pv_system_3_ss_controller_B.sum3_g[1] =
    my_pv_system_3_ss_controller_B.A21_o[1] +
    my_pv_system_3_ss_controller_B.A22_p[1];

  /* Gain: '<S54>/B11' */
  my_pv_system_3_ss_controller_B.B11_o[0] =
    my_pv_system_3_ss_controller_P.B11_Gain_n *
    my_pv_system_3_ss_controller_B.Switch[0];
  my_pv_system_3_ss_controller_B.B11_o[1] =
    my_pv_system_3_ss_controller_P.B11_Gain_n *
    my_pv_system_3_ss_controller_B.Switch[1];

  /* Sum: '<S52>/A*x1(k) + B*u1(k) ' */
  my_pv_system_3_ss_controller_B.x1k1_l[0] =
    my_pv_system_3_ss_controller_B.sum2_e[0] +
    my_pv_system_3_ss_controller_B.B11_o[0];
  my_pv_system_3_ss_controller_B.x1k1_l[1] =
    my_pv_system_3_ss_controller_B.sum2_e[1] +
    my_pv_system_3_ss_controller_B.B11_o[1];

  /* Gain: '<S54>/B21' */
  my_pv_system_3_ss_controller_B.B21_p[0] =
    my_pv_system_3_ss_controller_P.B21_Gain_m *
    my_pv_system_3_ss_controller_B.Switch[0];
  my_pv_system_3_ss_controller_B.B21_p[1] =
    my_pv_system_3_ss_controller_P.B21_Gain_m *
    my_pv_system_3_ss_controller_B.Switch[1];

  /* Sum: '<S52>/A*x2(k) + B*u2(k)' */
  my_pv_system_3_ss_controller_B.x2k1_k[0] =
    my_pv_system_3_ss_controller_B.sum3_g[0] +
    my_pv_system_3_ss_controller_B.B21_p[0];
  my_pv_system_3_ss_controller_B.x2k1_k[1] =
    my_pv_system_3_ss_controller_B.sum3_g[1] +
    my_pv_system_3_ss_controller_B.B21_p[1];

  /* Constant: '<S20>/Constant1' */
  my_pv_system_3_ss_controller_B.Constant1 =
    my_pv_system_3_ss_controller_P.Constant1_Value_p;

  /* Sum: '<S62>/Add3' incorporates:
   *  Constant: '<S12>/Constant10'
   */
  my_pv_system_3_ss_controller_B.Add3_n =
    my_pv_system_3_ss_controller_P.PWM_Generator_MinMax[1] -
    my_pv_system_3_ss_controller_P.PWM_Generator_MinMax[0];

  /* DigitalClock: '<S84>/Digital Clock' */
  my_pv_system_3_ss_controller_B.DigitalClock_mz =
    my_pv_system_3_ss_controller_M->Timing.t[0];

  /* Sum: '<S84>/Add1' incorporates:
   *  Constant: '<S84>/Constant3'
   */
  my_pv_system_3_ss_controller_B.Add1 =
    my_pv_system_3_ss_controller_B.DigitalClock_mz +
    my_pv_system_3_ss_controller_P.Constant3_Value_c;

  /* Math: '<S84>/Math Function' incorporates:
   *  Constant: '<S84>/Constant1'
   */
  my_pv_system_3_ss_controller_B.MathFunction_b = rt_remd_snf
    (my_pv_system_3_ss_controller_B.Add1,
     my_pv_system_3_ss_controller_P.Constant1_Value_h);

  /* Gain: '<S84>/1\ib1' */
  my_pv_system_3_ss_controller_B.ib1 = my_pv_system_3_ss_controller_P.ib1_Gain *
    my_pv_system_3_ss_controller_B.MathFunction_b;

  /* Lookup: '<S84>/Lookup Table'
   * About '<S84>/Lookup Table':
   * Input0  Data Type:  Floating Point real_T
   * Output0 Data Type:  Floating Point real_T
   * Lookup Method: Linear_Endpoint
   *
   * XData parameter uses the same data type and scaling as Input0
   * YData parameter uses the same data type and scaling as Output0
   */
  LookUp_real_T_real_T( &(my_pv_system_3_ss_controller_B.LookupTable),
                       my_pv_system_3_ss_controller_P.LookupTable_YData,
                       my_pv_system_3_ss_controller_B.ib1,
                       my_pv_system_3_ss_controller_P.LookupTable_XData, 2U);

  /* Sum: '<S84>/Add3' incorporates:
   *  Constant: '<S84>/Constant2'
   */
  my_pv_system_3_ss_controller_B.Add3_o =
    my_pv_system_3_ss_controller_B.LookupTable -
    my_pv_system_3_ss_controller_P.Constant2_Value_a;

  /* Gain: '<S62>/Gain1' */
  my_pv_system_3_ss_controller_B.Gain1_p =
    my_pv_system_3_ss_controller_P.Gain1_Gain_du *
    my_pv_system_3_ss_controller_B.Add3_n;

  /* Product: '<S62>/MUL1' */
  my_pv_system_3_ss_controller_B.MUL1 = my_pv_system_3_ss_controller_B.Add3_o *
    my_pv_system_3_ss_controller_B.Gain1_p;

  /* Sum: '<S62>/Add4' incorporates:
   *  Constant: '<S12>/Constant10'
   */
  my_pv_system_3_ss_controller_B.Add4 =
    (my_pv_system_3_ss_controller_P.PWM_Generator_MinMax[0] +
     my_pv_system_3_ss_controller_B.MUL1) +
    my_pv_system_3_ss_controller_B.Gain1_p;

  /* UnitDelay: '<S4>/Unit Delay' */
  my_pv_system_3_ss_controller_B.UnitDelay_m =
    my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_il;

  /* RelationalOperator: '<S67>/Relational Operator1' */
  my_pv_system_3_ss_controller_B.RelationalOperator1 =
    (my_pv_system_3_ss_controller_B.UnitDelay_m >=
     my_pv_system_3_ss_controller_B.Add4);

  /* Gain: '<S67>/Gain' */
  my_pv_system_3_ss_controller_B.Gain_a =
    my_pv_system_3_ss_controller_P.Gain_Gain_cv *
    my_pv_system_3_ss_controller_B.UnitDelay_m;

  /* RelationalOperator: '<S67>/Relational Operator3' */
  my_pv_system_3_ss_controller_B.RelationalOperator3 =
    (my_pv_system_3_ss_controller_B.Gain_a >=
     my_pv_system_3_ss_controller_B.Add4);

  /* Logic: '<S12>/Logical Operator4' */
  my_pv_system_3_ss_controller_B.LogicalOperator4[0] =
    !my_pv_system_3_ss_controller_B.RelationalOperator1;
  my_pv_system_3_ss_controller_B.LogicalOperator4[1] =
    !my_pv_system_3_ss_controller_B.RelationalOperator3;

  /* DataTypeConversion: '<S12>/Data Type Conversion' */
  my_pv_system_3_ss_controller_B.DataTypeConversion[0] =
    my_pv_system_3_ss_controller_B.RelationalOperator1;
  my_pv_system_3_ss_controller_B.DataTypeConversion[1] =
    my_pv_system_3_ss_controller_B.LogicalOperator4[0];
  my_pv_system_3_ss_controller_B.DataTypeConversion[2] =
    my_pv_system_3_ss_controller_B.RelationalOperator3;
  my_pv_system_3_ss_controller_B.DataTypeConversion[3] =
    my_pv_system_3_ss_controller_B.LogicalOperator4[1];

  /* Switch: '<S4>/Switch' incorporates:
   *  Constant: '<S2>/MPPT_On'
   *  Constant: '<S4>/Vnom_dc1'
   */
  if (my_pv_system_3_ss_controller_P.MPPT_On_Value != 0.0) {
    my_pv_system_3_ss_controller_B.Switch_d2 = my_pv_system_3_ss_controller_B.D;
  } else {
    my_pv_system_3_ss_controller_B.Switch_d2 =
      my_pv_system_3_ss_controller_P.InverterControl_Vdc_ref_Init;
  }

  /* End of Switch: '<S4>/Switch' */

  /* Sum: '<S13>/Add1' incorporates:
   *  Constant: '<S13>/Constant2'
   *  Constant: '<S13>/Constant4'
   */
  dV = my_pv_system_3_ss_controller_P.InverterControl_Ts_Control *
    my_pv_system_3_ss_controller_P.InverterControl_Fnom * 6.2831853071795862;
  my_pv_system_3_ss_controller_B.Add1_a =
    (my_pv_system_3_ss_controller_B.MathFunction_i +
     my_pv_system_3_ss_controller_P.Constant2_Value_b) + dV;

  /* UnitDelay: '<S4>/Unit Delay3' */
  my_pv_system_3_ss_controller_B.UnitDelay3[0] =
    my_pv_system_3_ss_controller_DW.UnitDelay3_DSTATE[0];
  my_pv_system_3_ss_controller_B.UnitDelay3[1] =
    my_pv_system_3_ss_controller_DW.UnitDelay3_DSTATE[1];

  /* Gain: '<S13>/Gain1' */
  my_pv_system_3_ss_controller_B.Gain1_g =
    my_pv_system_3_ss_controller_P.Gain1_Gain_p *
    my_pv_system_3_ss_controller_B.Switch_c;

  /* Product: '<S13>/Product' incorporates:
   *  Constant: '<S13>/Constant3'
   */
  dV = my_pv_system_3_ss_controller_P.InverterControl_Vnom_prim *
    1.4142135623730951;
  my_pv_system_3_ss_controller_B.Product =
    my_pv_system_3_ss_controller_B.Gain1_g / dV;

  /* Product: '<S13>/Product1' */
  my_pv_system_3_ss_controller_B.Product1_d[0] =
    my_pv_system_3_ss_controller_B.UnitDelay3[0] /
    my_pv_system_3_ss_controller_B.Product;
  my_pv_system_3_ss_controller_B.Product1_d[1] =
    my_pv_system_3_ss_controller_B.UnitDelay3[1] /
    my_pv_system_3_ss_controller_B.Product;

  /* RealImagToComplex: '<S13>/Real-Imag to Complex' */
  my_pv_system_3_ss_controller_B.RealImagtoComplex_n.re =
    my_pv_system_3_ss_controller_B.Product1_d[0];
  my_pv_system_3_ss_controller_B.RealImagtoComplex_n.im =
    my_pv_system_3_ss_controller_B.Product1_d[1];

  /* ComplexToMagnitudeAngle: '<S13>/Complex to Magnitude-Angle' */
  my_pv_system_3_ss_controller_B.ComplextoMagnitudeAngle_o1_n = rt_hypotd_snf
    (my_pv_system_3_ss_controller_B.RealImagtoComplex_n.re,
     my_pv_system_3_ss_controller_B.RealImagtoComplex_n.im);
  my_pv_system_3_ss_controller_B.ComplextoMagnitudeAngle_o2_p = rt_atan2d_snf
    (my_pv_system_3_ss_controller_B.RealImagtoComplex_n.im,
     my_pv_system_3_ss_controller_B.RealImagtoComplex_n.re);

  /* Sum: '<S13>/Add2' */
  my_pv_system_3_ss_controller_B.Add2_f =
    my_pv_system_3_ss_controller_B.ComplextoMagnitudeAngle_o2_p +
    my_pv_system_3_ss_controller_B.Add1_a;

  /* Trigonometry: '<S13>/Trigonometric Function' */
  my_pv_system_3_ss_controller_B.TrigonometricFunction_h = sin
    (my_pv_system_3_ss_controller_B.Add2_f);

  /* Product: '<S13>/Product2' */
  my_pv_system_3_ss_controller_B.Product2_n =
    my_pv_system_3_ss_controller_B.ComplextoMagnitudeAngle_o1_n *
    my_pv_system_3_ss_controller_B.TrigonometricFunction_h;

  /* UnitDelay: '<S4>/Unit Delay1' */
  my_pv_system_3_ss_controller_B.UnitDelay1_c =
    my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_e;

  /* Sum: '<S14>/Sum' */
  my_pv_system_3_ss_controller_B.Sum_k = my_pv_system_3_ss_controller_B.Switch_c
    - my_pv_system_3_ss_controller_B.UnitDelay1_c;

  /* Gain: '<S14>/Rtot_pu2' */
  dP = my_pv_system_3_ss_controller_P.InverterControl_Vnom_dc;
  dV = 1.0 / dP;
  my_pv_system_3_ss_controller_B.Rtot_pu2 = dV *
    my_pv_system_3_ss_controller_B.Sum_k;

  /* Gain: '<S86>/Integral Gain' */
  my_pv_system_3_ss_controller_B.IntegralGain_a =
    my_pv_system_3_ss_controller_P.InverterControl_Ki_VDCreg *
    my_pv_system_3_ss_controller_B.Rtot_pu2;

  /* DiscreteIntegrator: '<S86>/Integrator' */
  my_pv_system_3_ss_controller_B.Integrator_e =
    my_pv_system_3_ss_controller_P.Integrator_gainval_i *
    my_pv_system_3_ss_controller_B.IntegralGain_a +
    my_pv_system_3_ss_controller_DW.Integrator_DSTATE_d;

  /* Gain: '<S86>/Proportional Gain' */
  my_pv_system_3_ss_controller_B.ProportionalGain_b =
    my_pv_system_3_ss_controller_P.InverterControl_Kp_VDCreg *
    my_pv_system_3_ss_controller_B.Rtot_pu2;

  /* Sum: '<S86>/Sum' */
  my_pv_system_3_ss_controller_B.Sum_i =
    my_pv_system_3_ss_controller_B.ProportionalGain_b +
    my_pv_system_3_ss_controller_B.Integrator_e;

  /* Saturate: '<S86>/Saturate' */
  u0 = my_pv_system_3_ss_controller_B.Sum_i;
  dP = my_pv_system_3_ss_controller_P.PI_LowerSaturationLimit_p;
  u2 = my_pv_system_3_ss_controller_P.PI_UpperSaturationLimit_f;
  if (u0 > u2) {
    my_pv_system_3_ss_controller_B.Saturate_p = u2;
  } else if (u0 < dP) {
    my_pv_system_3_ss_controller_B.Saturate_p = dP;
  } else {
    my_pv_system_3_ss_controller_B.Saturate_p = u0;
  }

  /* End of Saturate: '<S86>/Saturate' */

  /* MATLAB Function: '<S2>/MATLAB Function' incorporates:
   *  Constant: '<S2>/Constant1'
   *  Constant: '<S2>/Constant2'
   *  Constant: '<S2>/Constant3'
   *  Constant: '<S2>/Constant4'
   */
  /* MATLAB Function 'SS_controller/MATLAB Function': '<S5>:1' */
  /* '<S5>:1:3' */
  /* '<S5>:1:4' */
  /* '<S5>:1:6' */
  /* '<S5>:1:7' */
  /* '<S5>:1:8' */
  /*  bandgap energy(1.12 eV) */
  /*  ideal factor */
  /* '<S5>:1:13' */
  /* '<S5>:1:14' */
  /* '<S5>:1:17' */
  my_pv_system_3_ss_controller_B.Iph_p =
    (((my_pv_system_3_ss_controller_B.SFunction_g[4] + 273.0) - 298.0) *
     (my_pv_system_3_ss_controller_P.Constant3_Value / 100.0) + 1.0) *
    (my_pv_system_3_ss_controller_B.SFunction_g[0] / 1000.0 *
     my_pv_system_3_ss_controller_P.Constant1_Value);

  /* '<S5>:1:18' */
  /* '<S5>:1:19' */
  /* '<S5>:1:21' */
  my_pv_system_3_ss_controller_B.Io_c =
    my_pv_system_3_ss_controller_P.Constant1_Value / (exp
    (my_pv_system_3_ss_controller_P.Constant2_Value /
     (my_pv_system_3_ss_controller_P.Constant4_Value * 1.381E-23 * 1.3 * 298.0 /
      1.6E-19)) - 1.0) * rt_powd_snf
    ((my_pv_system_3_ss_controller_B.SFunction_g[4] + 273.0) / 298.0, 3.0) * exp
    ((0.0033557046979865771 - 1.0 / (my_pv_system_3_ss_controller_B.SFunction_g
       [4] + 273.0)) * 9981.61867097421) * (exp
    (my_pv_system_3_ss_controller_B.SFunction_g[5] /
     (my_pv_system_3_ss_controller_P.Constant4_Value * 1.381E-23 * 1.3 *
      (my_pv_system_3_ss_controller_B.SFunction_g[4] + 273.0) / 1.6E-19)) - 1.0);
}

/* Model update function */
static void my_pv_system_3_ss_controller_update(void)
{
  /* Update for Memory: '<S1>/S-Function' */
  my_pv_system_3_ss_controller_DW.SFunction_PreviousInput =
    my_pv_system_3_ss_controller_B.Sum;

  /* Update for Memory: '<S2>/Memory' */
  my_pv_system_3_ss_controller_DW.Memory_1_PreviousInput[0] =
    my_pv_system_3_ss_controller_B.DataTypeConversion[0];
  my_pv_system_3_ss_controller_DW.Memory_1_PreviousInput[1] =
    my_pv_system_3_ss_controller_B.DataTypeConversion[1];
  my_pv_system_3_ss_controller_DW.Memory_1_PreviousInput[2] =
    my_pv_system_3_ss_controller_B.DataTypeConversion[2];
  my_pv_system_3_ss_controller_DW.Memory_1_PreviousInput[3] =
    my_pv_system_3_ss_controller_B.DataTypeConversion[3];

  /* Update for Memory: '<S2>/Memory' */
  my_pv_system_3_ss_controller_DW.Memory_2_PreviousInput =
    my_pv_system_3_ss_controller_B.Iph_p;

  /* Update for Memory: '<S2>/Memory' */
  my_pv_system_3_ss_controller_DW.Memory_3_PreviousInput =
    my_pv_system_3_ss_controller_B.Io_c;

  /* Update for UnitDelay: '<S31>/Unit Delay' */
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE =
    my_pv_system_3_ss_controller_B.yk_e;

  /* Level2 S-Function Block: '<S57>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[6];
    sfcnUpdate(rts, 0);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Update for DiscreteIntegrator: '<S31>/Discrete-Time Integrator' */
  my_pv_system_3_ss_controller_DW.DiscreteTimeIntegrator_DSTATE +=
    my_pv_system_3_ss_controller_P.DiscreteTimeIntegrator_gainval *
    my_pv_system_3_ss_controller_B.Saturation1;

  /* Update for UnitDelay: '<S52>/Delay_x1' */
  my_pv_system_3_ss_controller_DW.Delay_x1_DSTATE[0] =
    my_pv_system_3_ss_controller_B.x1k1_l[0];
  my_pv_system_3_ss_controller_DW.Delay_x1_DSTATE[1] =
    my_pv_system_3_ss_controller_B.x1k1_l[1];

  /* Update for UnitDelay: '<S52>/Delay_x2' */
  my_pv_system_3_ss_controller_DW.Delay_x2_DSTATE[0] =
    my_pv_system_3_ss_controller_B.x2k1_k[0];
  my_pv_system_3_ss_controller_DW.Delay_x2_DSTATE[1] =
    my_pv_system_3_ss_controller_B.x2k1_k[1];

  /* Update for UnitDelay: '<S4>/Unit Delay2' */
  my_pv_system_3_ss_controller_DW.UnitDelay2_DSTATE =
    my_pv_system_3_ss_controller_B.Saturate_p;

  /* Update for DiscreteIntegrator: '<S15>/Integrator' */
  my_pv_system_3_ss_controller_DW.Integrator_DSTATE[0] +=
    my_pv_system_3_ss_controller_P.Integrator_gainval *
    my_pv_system_3_ss_controller_B.IntegralGain[0];
  my_pv_system_3_ss_controller_DW.Integrator_DSTATE[1] +=
    my_pv_system_3_ss_controller_P.Integrator_gainval *
    my_pv_system_3_ss_controller_B.IntegralGain[1];

  /* Update for DiscreteIntegrator: '<S23>/Integ4' */
  my_pv_system_3_ss_controller_DW.Integ4_SYSTEM_ENABLE = 0U;
  my_pv_system_3_ss_controller_DW.Integ4_DSTATE =
    my_pv_system_3_ss_controller_P.Integ4_gainval_o *
    my_pv_system_3_ss_controller_B.Product1 +
    my_pv_system_3_ss_controller_B.Integ4;

  /* Level2 S-Function Block: '<S25>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[7];
    sfcnUpdate(rts, 0);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Update for UnitDelay: '<S24>/Unit Delay' */
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_g =
    my_pv_system_3_ss_controller_B.Product1;

  /* Update for UnitDelay: '<S23>/Unit Delay1' */
  my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE =
    my_pv_system_3_ss_controller_B.Switch_g;

  /* Update for DiscreteIntegrator: '<S26>/Integ4' */
  my_pv_system_3_ss_controller_DW.Integ4_SYSTEM_ENABLE_e = 0U;
  my_pv_system_3_ss_controller_DW.Integ4_DSTATE_d =
    my_pv_system_3_ss_controller_P.Integ4_gainval_f *
    my_pv_system_3_ss_controller_B.Product2 +
    my_pv_system_3_ss_controller_B.Integ4_h;

  /* Level2 S-Function Block: '<S28>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[8];
    sfcnUpdate(rts, 0);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Update for UnitDelay: '<S27>/Unit Delay' */
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_i =
    my_pv_system_3_ss_controller_B.Product2;

  /* Update for UnitDelay: '<S26>/Unit Delay1' */
  my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_b =
    my_pv_system_3_ss_controller_B.Switch_k;

  /* Update for DiscreteIntegrator: '<S29>/Integ4' */
  my_pv_system_3_ss_controller_DW.Integ4_SYSTEM_ENABLE_d = 0U;
  my_pv_system_3_ss_controller_DW.Integ4_DSTATE_d5 =
    my_pv_system_3_ss_controller_P.Integ4_gainval_j *
    my_pv_system_3_ss_controller_B.RateTransition2 +
    my_pv_system_3_ss_controller_B.Integ4_c;

  /* Level2 S-Function Block: '<S30>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[9];
    sfcnUpdate(rts, 0);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Update for UnitDelay: '<S29>/Unit Delay' */
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_iq =
    my_pv_system_3_ss_controller_B.RateTransition2;

  /* Update for UnitDelay: '<S29>/Unit Delay1' */
  my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_a =
    my_pv_system_3_ss_controller_B.Switch_c;

  /* Update for Enabled SubSystem: '<S31>/Automatic Gain Control' incorporates:
   *  Update for EnablePort: '<S32>/Enable'
   */
  if (my_pv_system_3_ss_controller_DW.AutomaticGainControl_MODE) {
    /* Update for DiscreteIntegrator: '<S39>/Integ4' */
    my_pv_system_3_ss_controller_DW.Integ4_SYSTEM_ENABLE_g = 0U;
    my_pv_system_3_ss_controller_DW.Integ4_DSTATE_ky =
      my_pv_system_3_ss_controller_P.Integ4_gainval *
      my_pv_system_3_ss_controller_B.Product1_o +
      my_pv_system_3_ss_controller_B.Integ4_d;

    /* Level2 S-Function Block: '<S41>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[0];
      sfcnUpdate(rts, 0);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Update for UnitDelay: '<S40>/Unit Delay' */
    my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_l =
      my_pv_system_3_ss_controller_B.Product1_o;

    /* Update for UnitDelay: '<S39>/Unit Delay1' */
    my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_o =
      my_pv_system_3_ss_controller_B.Switch_j;

    /* Update for DiscreteIntegrator: '<S42>/Integ4' */
    my_pv_system_3_ss_controller_DW.Integ4_SYSTEM_ENABLE_gx = 0U;
    my_pv_system_3_ss_controller_DW.Integ4_DSTATE_j =
      my_pv_system_3_ss_controller_P.Integ4_gainval_h *
      my_pv_system_3_ss_controller_B.Product2_i +
      my_pv_system_3_ss_controller_B.Integ4_dy;

    /* Level2 S-Function Block: '<S44>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[1];
      sfcnUpdate(rts, 0);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Update for UnitDelay: '<S43>/Unit Delay' */
    my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_d =
      my_pv_system_3_ss_controller_B.Product2_i;

    /* Update for UnitDelay: '<S42>/Unit Delay1' */
    my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_k =
      my_pv_system_3_ss_controller_B.Switch_m;
  }

  /* End of Update for SubSystem: '<S31>/Automatic Gain Control' */

  /* Update for DiscreteIntegrator: '<S45>/Integ4' */
  my_pv_system_3_ss_controller_DW.Integ4_SYSTEM_ENABLE_k = 0U;
  my_pv_system_3_ss_controller_DW.Integ4_DSTATE_k =
    my_pv_system_3_ss_controller_P.Integ4_gainval_jy *
    my_pv_system_3_ss_controller_B.Product1_c +
    my_pv_system_3_ss_controller_B.Integ4_e;

  /* Level2 S-Function Block: '<S47>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[10];
    sfcnUpdate(rts, 0);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Update for UnitDelay: '<S46>/Unit Delay' */
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_a =
    my_pv_system_3_ss_controller_B.Product1_c;

  /* Update for UnitDelay: '<S45>/Unit Delay1' */
  my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_i =
    my_pv_system_3_ss_controller_B.Switch_d;

  /* Update for DiscreteTransferFcn: '<S33>/Discrete Derivative ' */
  my_pv_system_3_ss_controller_DW.DiscreteDerivative_states =
    my_pv_system_3_ss_controller_DW.DiscreteDerivative_tmp;

  /* Update for DiscreteIntegrator: '<S33>/Discrete-Time Integrator' */
  my_pv_system_3_ss_controller_DW.DiscreteTimeIntegrator_DSTATE_d +=
    my_pv_system_3_ss_controller_P.DiscreteTimeIntegrator_gainva_j *
    my_pv_system_3_ss_controller_B.Divide;
  if (my_pv_system_3_ss_controller_DW.DiscreteTimeIntegrator_DSTATE_d >=
      my_pv_system_3_ss_controller_P.DiscreteTimeIntegrator_UpperSat) {
    my_pv_system_3_ss_controller_DW.DiscreteTimeIntegrator_DSTATE_d =
      my_pv_system_3_ss_controller_P.DiscreteTimeIntegrator_UpperSat;
  } else {
    if (my_pv_system_3_ss_controller_DW.DiscreteTimeIntegrator_DSTATE_d <=
        my_pv_system_3_ss_controller_P.DiscreteTimeIntegrator_LowerSat) {
      my_pv_system_3_ss_controller_DW.DiscreteTimeIntegrator_DSTATE_d =
        my_pv_system_3_ss_controller_P.DiscreteTimeIntegrator_LowerSat;
    }
  }

  /* End of Update for DiscreteIntegrator: '<S33>/Discrete-Time Integrator' */

  /* Update for UnitDelay: '<S48>/Delay_x1' */
  my_pv_system_3_ss_controller_DW.Delay_x1_DSTATE_d =
    my_pv_system_3_ss_controller_B.x1k1;

  /* Update for UnitDelay: '<S48>/Delay_x2' */
  my_pv_system_3_ss_controller_DW.Delay_x2_DSTATE_m =
    my_pv_system_3_ss_controller_B.x2k1;

  /* Update for UnitDelay: '<S4>/Unit Delay' */
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_il =
    my_pv_system_3_ss_controller_B.Product2_n;

  /* Update for UnitDelay: '<S4>/Unit Delay3' */
  my_pv_system_3_ss_controller_DW.UnitDelay3_DSTATE[0] =
    my_pv_system_3_ss_controller_B.Saturation[0];
  my_pv_system_3_ss_controller_DW.UnitDelay3_DSTATE[1] =
    my_pv_system_3_ss_controller_B.Saturation[1];

  /* Update for UnitDelay: '<S4>/Unit Delay1' */
  my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_e =
    my_pv_system_3_ss_controller_B.Switch_d2;

  /* Update for DiscreteIntegrator: '<S86>/Integrator' */
  my_pv_system_3_ss_controller_DW.Integrator_DSTATE_d =
    my_pv_system_3_ss_controller_P.Integrator_gainval_i *
    my_pv_system_3_ss_controller_B.IntegralGain_a +
    my_pv_system_3_ss_controller_B.Integrator_e;

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++my_pv_system_3_ss_controller_M->Timing.clockTick0)) {
    ++my_pv_system_3_ss_controller_M->Timing.clockTickH0;
  }

  my_pv_system_3_ss_controller_M->Timing.t[0] =
    my_pv_system_3_ss_controller_M->Timing.clockTick0 *
    my_pv_system_3_ss_controller_M->Timing.stepSize0 +
    my_pv_system_3_ss_controller_M->Timing.clockTickH0 *
    my_pv_system_3_ss_controller_M->Timing.stepSize0 * 4294967296.0;
}

/* Model initialize function */
static void my_pv_system_3_ss_controller_initialize(void)
{
  /* Level2 S-Function Block: '<S57>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[6];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Start for Enabled SubSystem: '<S56>/Subsystem1' */
  /* VirtualOutportStart for Outport: '<S61>/dq' */
  my_pv_system_3_ss_controller_B.Fcn = my_pv_system_3_ss_controller_P.dq_Y0_p[0];
  my_pv_system_3_ss_controller_B.Fcn1 = my_pv_system_3_ss_controller_P.dq_Y0_p[1];

  /* End of Start for SubSystem: '<S56>/Subsystem1' */

  /* Start for Enabled SubSystem: '<S56>/Subsystem - pi//2 delay' */
  /* VirtualOutportStart for Outport: '<S60>/dq' */
  my_pv_system_3_ss_controller_B.Fcn_l = my_pv_system_3_ss_controller_P.dq_Y0[0];
  my_pv_system_3_ss_controller_B.Fcn1_f = my_pv_system_3_ss_controller_P.dq_Y0[1];

  /* End of Start for SubSystem: '<S56>/Subsystem - pi//2 delay' */
  /* Level2 S-Function Block: '<S25>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[7];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Level2 S-Function Block: '<S28>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[8];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Start for Constant: '<S29>/K1' */
  my_pv_system_3_ss_controller_B.K1 = my_pv_system_3_ss_controller_P.K1_Value;

  /* Level2 S-Function Block: '<S30>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[9];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Start for Enabled SubSystem: '<S31>/Automatic Gain Control' */

  /* Level2 S-Function Block: '<S41>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[0];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Level2 S-Function Block: '<S44>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[1];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* End of Start for SubSystem: '<S31>/Automatic Gain Control' */

  /* InitializeConditions for Enabled SubSystem: '<S31>/Automatic Gain Control' */
  /* InitializeConditions for DiscreteIntegrator: '<S39>/Integ4' */
  my_pv_system_3_ss_controller_DW.Integ4_DSTATE_ky =
    my_pv_system_3_ss_controller_P.Integ4_IC;

  /* Level2 S-Function Block: '<S41>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[0];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* InitializeConditions for UnitDelay: '<S40>/Unit Delay' */
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_l =
    my_pv_system_3_ss_controller_P.UnitDelay_InitialCondition;

  /* InitializeConditions for UnitDelay: '<S39>/Unit Delay1' */
  my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_o =
    my_pv_system_3_ss_controller_P.UnitDelay1_InitialCondition;

  /* InitializeConditions for DiscreteIntegrator: '<S42>/Integ4' */
  my_pv_system_3_ss_controller_DW.Integ4_DSTATE_j =
    my_pv_system_3_ss_controller_P.Integ4_IC_k;

  /* Level2 S-Function Block: '<S44>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[1];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* InitializeConditions for UnitDelay: '<S43>/Unit Delay' */
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_d =
    my_pv_system_3_ss_controller_P.UnitDelay_InitialCondition_c;

  /* InitializeConditions for UnitDelay: '<S42>/Unit Delay1' */
  my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_k =
    my_pv_system_3_ss_controller_P.UnitDelay1_InitialCondition_c;

  /* End of InitializeConditions for SubSystem: '<S31>/Automatic Gain Control' */

  /* Start for Enabled SubSystem: '<S31>/Automatic Gain Control' */
  /* VirtualOutportStart for Outport: '<S32>/Gain' */
  my_pv_system_3_ss_controller_B.MathFunction_k =
    my_pv_system_3_ss_controller_P.Gain_Y0;

  /* End of Start for SubSystem: '<S31>/Automatic Gain Control' */
  /* Level2 S-Function Block: '<S47>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[10];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
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

  /* InitializeConditions for Memory: '<S1>/S-Function' */
  my_pv_system_3_ss_controller_DW.SFunction_PreviousInput =
    my_pv_system_3_ss_controller_P.SFunction_X0;

  /* InitializeConditions for Memory: '<S2>/Memory' */
  my_pv_system_3_ss_controller_DW.Memory_1_PreviousInput[0] =
    my_pv_system_3_ss_controller_P.Memory_1_X0;
  my_pv_system_3_ss_controller_DW.Memory_1_PreviousInput[1] =
    my_pv_system_3_ss_controller_P.Memory_1_X0;
  my_pv_system_3_ss_controller_DW.Memory_1_PreviousInput[2] =
    my_pv_system_3_ss_controller_P.Memory_1_X0;
  my_pv_system_3_ss_controller_DW.Memory_1_PreviousInput[3] =
    my_pv_system_3_ss_controller_P.Memory_1_X0;

  /* InitializeConditions for Memory: '<S2>/Memory' */
  my_pv_system_3_ss_controller_DW.Memory_2_PreviousInput =
    my_pv_system_3_ss_controller_P.Memory_2_X0;

  /* InitializeConditions for Memory: '<S2>/Memory' */
  my_pv_system_3_ss_controller_DW.Memory_3_PreviousInput =
    my_pv_system_3_ss_controller_P.Memory_3_X0;

  /* InitializeConditions for Atomic SubSystem: '<S6>/Subsystem4' */

  /* Level2 S-Function Block: '<S90>/S-Function' (send_rt) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[2];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* End of InitializeConditions for SubSystem: '<S6>/Subsystem4' */

  /* Level2 S-Function Block: '<S2>/OpMonitor' (opmonitor) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[3];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Level2 S-Function Block: '<S89>/S-Function' (OP_SEND) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[4];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Level2 S-Function Block: '<S91>/S-Function' (recv_rt) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[5];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* InitializeConditions for UnitDelay: '<S31>/Unit Delay' */
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE =
    my_pv_system_3_ss_controller_P.UnitDelay_InitialCondition_a;

  /* Level2 S-Function Block: '<S57>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[6];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* InitializeConditions for DiscreteIntegrator: '<S31>/Discrete-Time Integrator' */
  my_pv_system_3_ss_controller_DW.DiscreteTimeIntegrator_DSTATE =
    my_pv_system_3_ss_controller_P.DiscreteTimeIntegrator_IC;

  /* InitializeConditions for UnitDelay: '<S52>/Delay_x1' */
  my_pv_system_3_ss_controller_DW.Delay_x1_DSTATE[0] =
    my_pv_system_3_ss_controller_P.Delay_x1_InitialCondition[0];
  my_pv_system_3_ss_controller_DW.Delay_x1_DSTATE[1] =
    my_pv_system_3_ss_controller_P.Delay_x1_InitialCondition[1];

  /* InitializeConditions for UnitDelay: '<S52>/Delay_x2' */
  my_pv_system_3_ss_controller_DW.Delay_x2_DSTATE[0] =
    my_pv_system_3_ss_controller_P.Delay_x2_InitialCondition[0];
  my_pv_system_3_ss_controller_DW.Delay_x2_DSTATE[1] =
    my_pv_system_3_ss_controller_P.Delay_x2_InitialCondition[1];

  /* InitializeConditions for UnitDelay: '<S4>/Unit Delay2' */
  my_pv_system_3_ss_controller_DW.UnitDelay2_DSTATE =
    my_pv_system_3_ss_controller_P.UnitDelay2_InitialCondition;

  /* InitializeConditions for DiscreteIntegrator: '<S15>/Integrator' */
  my_pv_system_3_ss_controller_DW.Integrator_DSTATE[0] =
    my_pv_system_3_ss_controller_P.Integrator_IC;
  my_pv_system_3_ss_controller_DW.Integrator_DSTATE[1] =
    my_pv_system_3_ss_controller_P.Integrator_IC;

  /* InitializeConditions for DiscreteIntegrator: '<S23>/Integ4' */
  my_pv_system_3_ss_controller_DW.Integ4_DSTATE =
    my_pv_system_3_ss_controller_P.Integ4_IC_a;

  /* Level2 S-Function Block: '<S25>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[7];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* InitializeConditions for UnitDelay: '<S24>/Unit Delay' */
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_g =
    my_pv_system_3_ss_controller_P.UnitDelay_InitialCondition_cj;

  /* InitializeConditions for UnitDelay: '<S23>/Unit Delay1' */
  my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE =
    my_pv_system_3_ss_controller_P.UnitDelay1_InitialCondition_p;

  /* InitializeConditions for DiscreteIntegrator: '<S26>/Integ4' */
  my_pv_system_3_ss_controller_DW.Integ4_DSTATE_d =
    my_pv_system_3_ss_controller_P.Integ4_IC_l;

  /* Level2 S-Function Block: '<S28>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[8];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* InitializeConditions for UnitDelay: '<S27>/Unit Delay' */
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_i =
    my_pv_system_3_ss_controller_P.UnitDelay_InitialCondition_l;

  /* InitializeConditions for UnitDelay: '<S26>/Unit Delay1' */
  my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_b =
    my_pv_system_3_ss_controller_P.UnitDelay1_InitialCondition_k;

  /* InitializeConditions for MATLAB Function: '<S4>/MPPT Controller using Perturbe  & Observe technique  ' */
  my_pv_system_3_ss_controller_DW.Vold_not_empty = false;
  my_pv_system_3_ss_controller_DW.Vold = 0.0;
  my_pv_system_3_ss_controller_DW.Pold = 0.0;

  /* InitializeConditions for DiscreteIntegrator: '<S29>/Integ4' */
  my_pv_system_3_ss_controller_DW.Integ4_DSTATE_d5 =
    my_pv_system_3_ss_controller_P.Integ4_IC_i;

  /* Level2 S-Function Block: '<S30>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[9];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* InitializeConditions for UnitDelay: '<S29>/Unit Delay' */
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_iq =
    my_pv_system_3_ss_controller_P.UnitDelay_InitialCondition_b;

  /* InitializeConditions for UnitDelay: '<S29>/Unit Delay1' */
  my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_a =
    my_pv_system_3_ss_controller_P.UnitDelay1_InitialCondition_kn;

  /* InitializeConditions for DiscreteIntegrator: '<S45>/Integ4' */
  my_pv_system_3_ss_controller_DW.Integ4_DSTATE_k =
    my_pv_system_3_ss_controller_P.Integ4_IC_h;

  /* Level2 S-Function Block: '<S47>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[10];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* InitializeConditions for UnitDelay: '<S46>/Unit Delay' */
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_a =
    my_pv_system_3_ss_controller_P.UnitDelay_InitialCondition_j;

  /* InitializeConditions for UnitDelay: '<S45>/Unit Delay1' */
  my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_i =
    my_pv_system_3_ss_controller_P.UnitDelay1_InitialCondition_j;

  /* InitializeConditions for DiscreteTransferFcn: '<S33>/Discrete Derivative ' */
  my_pv_system_3_ss_controller_DW.DiscreteDerivative_states =
    my_pv_system_3_ss_controller_P.DiscreteDerivative_InitialState;

  /* InitializeConditions for DiscreteIntegrator: '<S33>/Discrete-Time Integrator' */
  my_pv_system_3_ss_controller_DW.DiscreteTimeIntegrator_DSTATE_d =
    my_pv_system_3_ss_controller_P.Discrete_Init;

  /* InitializeConditions for RateLimiter: '<S31>/Rate Limiter' */
  my_pv_system_3_ss_controller_DW.PrevY =
    my_pv_system_3_ss_controller_P.RateLimiter_IC;

  /* InitializeConditions for UnitDelay: '<S48>/Delay_x1' */
  my_pv_system_3_ss_controller_DW.Delay_x1_DSTATE_d =
    my_pv_system_3_ss_controller_P.Delay_x1_InitialCondition_i;

  /* InitializeConditions for UnitDelay: '<S48>/Delay_x2' */
  my_pv_system_3_ss_controller_DW.Delay_x2_DSTATE_m =
    my_pv_system_3_ss_controller_P.Delay_x2_InitialCondition_g;

  /* InitializeConditions for UnitDelay: '<S4>/Unit Delay' */
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_il =
    my_pv_system_3_ss_controller_P.UnitDelay_InitialCondition_i;

  /* InitializeConditions for UnitDelay: '<S4>/Unit Delay3' */
  my_pv_system_3_ss_controller_DW.UnitDelay3_DSTATE[0] =
    my_pv_system_3_ss_controller_P.UnitDelay3_InitialCondition;
  my_pv_system_3_ss_controller_DW.UnitDelay3_DSTATE[1] =
    my_pv_system_3_ss_controller_P.UnitDelay3_InitialCondition;

  /* InitializeConditions for UnitDelay: '<S4>/Unit Delay1' */
  my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_e =
    my_pv_system_3_ss_controller_P.UnitDelay1_InitialCondition_n;

  /* InitializeConditions for DiscreteIntegrator: '<S86>/Integrator' */
  my_pv_system_3_ss_controller_DW.Integrator_DSTATE_d =
    my_pv_system_3_ss_controller_P.Integrator_IC_i;

  /* Enable for DiscreteIntegrator: '<S23>/Integ4' */
  my_pv_system_3_ss_controller_DW.Integ4_SYSTEM_ENABLE = 1U;

  /* Enable for DiscreteIntegrator: '<S26>/Integ4' */
  my_pv_system_3_ss_controller_DW.Integ4_SYSTEM_ENABLE_e = 1U;

  /* Enable for DiscreteIntegrator: '<S29>/Integ4' */
  my_pv_system_3_ss_controller_DW.Integ4_SYSTEM_ENABLE_d = 1U;

  /* Enable for DiscreteIntegrator: '<S45>/Integ4' */
  my_pv_system_3_ss_controller_DW.Integ4_SYSTEM_ENABLE_k = 1U;
}

/* Model terminate function */
static void my_pv_system_3_ss_controller_terminate(void)
{
  /* Terminate for Atomic SubSystem: '<S6>/Subsystem4' */

  /* Level2 S-Function Block: '<S90>/S-Function' (send_rt) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[2];
    sfcnTerminate(rts);
  }

  /* End of Terminate for SubSystem: '<S6>/Subsystem4' */

  /* Level2 S-Function Block: '<S2>/OpMonitor' (opmonitor) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[3];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S89>/S-Function' (OP_SEND) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[4];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S91>/S-Function' (recv_rt) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[5];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S57>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[6];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S25>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[7];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S28>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[8];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S30>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[9];
    sfcnTerminate(rts);
  }

  /* Terminate for Enabled SubSystem: '<S31>/Automatic Gain Control' */

  /* Level2 S-Function Block: '<S41>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[0];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S44>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[1];
    sfcnTerminate(rts);
  }

  /* End of Terminate for SubSystem: '<S31>/Automatic Gain Control' */

  /* Level2 S-Function Block: '<S47>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[10];
    sfcnTerminate(rts);
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/
void MdlOutputs(int_T tid)
{
  my_pv_system_3_ss_controller_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  my_pv_system_3_ss_controller_update();
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
  my_pv_system_3_ss_controller_initialize();
}

void MdlTerminate(void)
{
  my_pv_system_3_ss_controller_terminate();
}

/* Registration function */
RT_MODEL_my_pv_system_3_ss_controller_T *my_pv_system_3_ss_controller(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  my_pv_system_3_ss_controller_P.Saturation_UpperSat = rtInf;
  my_pv_system_3_ss_controller_P.DiscreteTimeIntegrator_UpperSat = rtInf;
  my_pv_system_3_ss_controller_P.Saturation1_UpperSat = rtInf;

  /* initialize real-time model */
  (void) memset((void *)my_pv_system_3_ss_controller_M, 0,
                sizeof(RT_MODEL_my_pv_system_3_ss_controller_T));
  rtsiSetSolverName(&my_pv_system_3_ss_controller_M->solverInfo,
                    "FixedStepDiscrete");
  my_pv_system_3_ss_controller_M->solverInfoPtr =
    (&my_pv_system_3_ss_controller_M->solverInfo);

  /* Initialize timing info */
  {
    int_T *mdlTsMap =
      my_pv_system_3_ss_controller_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    my_pv_system_3_ss_controller_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    my_pv_system_3_ss_controller_M->Timing.sampleTimes =
      (&my_pv_system_3_ss_controller_M->Timing.sampleTimesArray[0]);
    my_pv_system_3_ss_controller_M->Timing.offsetTimes =
      (&my_pv_system_3_ss_controller_M->Timing.offsetTimesArray[0]);

    /* task periods */
    my_pv_system_3_ss_controller_M->Timing.sampleTimes[0] = (5.0E-5);

    /* task offsets */
    my_pv_system_3_ss_controller_M->Timing.offsetTimes[0] = (0.0);
  }

  rtmSetTPtr(my_pv_system_3_ss_controller_M,
             &my_pv_system_3_ss_controller_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = my_pv_system_3_ss_controller_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    my_pv_system_3_ss_controller_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(my_pv_system_3_ss_controller_M, -1);
  my_pv_system_3_ss_controller_M->Timing.stepSize0 = 5.0E-5;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    my_pv_system_3_ss_controller_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(my_pv_system_3_ss_controller_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(my_pv_system_3_ss_controller_M->rtwLogInfo, (NULL));
    rtliSetLogT(my_pv_system_3_ss_controller_M->rtwLogInfo, "");
    rtliSetLogX(my_pv_system_3_ss_controller_M->rtwLogInfo, "");
    rtliSetLogXFinal(my_pv_system_3_ss_controller_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(my_pv_system_3_ss_controller_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(my_pv_system_3_ss_controller_M->rtwLogInfo, 0);
    rtliSetLogMaxRows(my_pv_system_3_ss_controller_M->rtwLogInfo, 1000);
    rtliSetLogDecimation(my_pv_system_3_ss_controller_M->rtwLogInfo, 1);
    rtliSetLogY(my_pv_system_3_ss_controller_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(my_pv_system_3_ss_controller_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(my_pv_system_3_ss_controller_M->rtwLogInfo, (NULL));
  }

  my_pv_system_3_ss_controller_M->solverInfoPtr =
    (&my_pv_system_3_ss_controller_M->solverInfo);
  my_pv_system_3_ss_controller_M->Timing.stepSize = (5.0E-5);
  rtsiSetFixedStepSize(&my_pv_system_3_ss_controller_M->solverInfo, 5.0E-5);
  rtsiSetSolverMode(&my_pv_system_3_ss_controller_M->solverInfo,
                    SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  my_pv_system_3_ss_controller_M->ModelData.blockIO = ((void *)
    &my_pv_system_3_ss_controller_B);
  (void) memset(((void *) &my_pv_system_3_ss_controller_B), 0,
                sizeof(B_my_pv_system_3_ss_controller_T));

  {
    int_T i;
    for (i = 0; i < 12; i++) {
      my_pv_system_3_ss_controller_B.SFunction_g[i] = 0.0;
    }

    my_pv_system_3_ss_controller_B.RealImagtoComplex.re = 0.0;
    my_pv_system_3_ss_controller_B.RealImagtoComplex.im = 0.0;
    my_pv_system_3_ss_controller_B.MagnitudeAngletoComplex.re = 0.0;
    my_pv_system_3_ss_controller_B.MagnitudeAngletoComplex.im = 0.0;
    my_pv_system_3_ss_controller_B.RealImagtoComplex_n.re = 0.0;
    my_pv_system_3_ss_controller_B.RealImagtoComplex_n.im = 0.0;
    my_pv_system_3_ss_controller_B.RealImagtoComplex_l.re = 0.0;
    my_pv_system_3_ss_controller_B.RealImagtoComplex_l.im = 0.0;
    my_pv_system_3_ss_controller_B.SFunction = 0.0;
    my_pv_system_3_ss_controller_B.Sum = 0.0;
    my_pv_system_3_ss_controller_B.Gate[0] = 0.0;
    my_pv_system_3_ss_controller_B.Gate[1] = 0.0;
    my_pv_system_3_ss_controller_B.Gate[2] = 0.0;
    my_pv_system_3_ss_controller_B.Gate[3] = 0.0;
    my_pv_system_3_ss_controller_B.Iph = 0.0;
    my_pv_system_3_ss_controller_B.Io = 0.0;
    my_pv_system_3_ss_controller_B.Computationtime = 0.0;
    my_pv_system_3_ss_controller_B.Realstepsize = 0.0;
    my_pv_system_3_ss_controller_B.Idletime = 0.0;
    my_pv_system_3_ss_controller_B.Overruntimes = 0.0;
    my_pv_system_3_ss_controller_B.RateTransition1 = 0.0;
    my_pv_system_3_ss_controller_B.Apu = 0.0;
    my_pv_system_3_ss_controller_B.UnitDelay = 0.0;
    my_pv_system_3_ss_controller_B.avoiddivisionbyzero = 0.0;
    my_pv_system_3_ss_controller_B.MathFunction = 0.0;
    my_pv_system_3_ss_controller_B.Gain = 0.0;
    my_pv_system_3_ss_controller_B.SFunction_o = 0.0;
    my_pv_system_3_ss_controller_B.DiscreteTimeIntegrator = 0.0;
    my_pv_system_3_ss_controller_B.MathFunction_i = 0.0;
    my_pv_system_3_ss_controller_B.FirstcycleofsimulationId092Iq0 = 0.0;
    my_pv_system_3_ss_controller_B.Switch[0] = 0.0;
    my_pv_system_3_ss_controller_B.Switch[1] = 0.0;
    my_pv_system_3_ss_controller_B.Duk[0] = 0.0;
    my_pv_system_3_ss_controller_B.Duk[1] = 0.0;
    my_pv_system_3_ss_controller_B.x1k[0] = 0.0;
    my_pv_system_3_ss_controller_B.x1k[1] = 0.0;
    my_pv_system_3_ss_controller_B.C11[0] = 0.0;
    my_pv_system_3_ss_controller_B.C11[1] = 0.0;
    my_pv_system_3_ss_controller_B.x2k[0] = 0.0;
    my_pv_system_3_ss_controller_B.x2k[1] = 0.0;
    my_pv_system_3_ss_controller_B.C12[0] = 0.0;
    my_pv_system_3_ss_controller_B.C12[1] = 0.0;
    my_pv_system_3_ss_controller_B.sum2[0] = 0.0;
    my_pv_system_3_ss_controller_B.sum2[1] = 0.0;
    my_pv_system_3_ss_controller_B.yk[0] = 0.0;
    my_pv_system_3_ss_controller_B.yk[1] = 0.0;
    my_pv_system_3_ss_controller_B.UnitDelay2 = 0.0;
    my_pv_system_3_ss_controller_B.Sum_e[0] = 0.0;
    my_pv_system_3_ss_controller_B.Sum_e[1] = 0.0;
    my_pv_system_3_ss_controller_B.ProportionalGain[0] = 0.0;
    my_pv_system_3_ss_controller_B.ProportionalGain[1] = 0.0;
    my_pv_system_3_ss_controller_B.Integrator[0] = 0.0;
    my_pv_system_3_ss_controller_B.Integrator[1] = 0.0;
    my_pv_system_3_ss_controller_B.Sum_g[0] = 0.0;
    my_pv_system_3_ss_controller_B.Sum_g[1] = 0.0;
    my_pv_system_3_ss_controller_B.Saturate[0] = 0.0;
    my_pv_system_3_ss_controller_B.Saturate[1] = 0.0;
    my_pv_system_3_ss_controller_B.RateTransition = 0.0;
    my_pv_system_3_ss_controller_B.Vpu = 0.0;
    my_pv_system_3_ss_controller_B.TrigonometricFunction = 0.0;
    my_pv_system_3_ss_controller_B.Gain1 = 0.0;
    my_pv_system_3_ss_controller_B.Product1 = 0.0;
    my_pv_system_3_ss_controller_B.Integ4 = 0.0;
    my_pv_system_3_ss_controller_B.Freq = 0.0;
    my_pv_system_3_ss_controller_B.Numberofsamplespercycle = 0.0;
    my_pv_system_3_ss_controller_B.RoundingFunction = 0.0;
    my_pv_system_3_ss_controller_B.Delay = 0.0;
    my_pv_system_3_ss_controller_B.SFunction_p = 0.0;
    my_pv_system_3_ss_controller_B.UnitDelay_f = 0.0;
    my_pv_system_3_ss_controller_B.DigitalClock = 0.0;
    my_pv_system_3_ss_controller_B.UnitDelay1 = 0.0;
    my_pv_system_3_ss_controller_B.Switch_g = 0.0;
    my_pv_system_3_ss_controller_B.TrigonometricFunction3 = 0.0;
    my_pv_system_3_ss_controller_B.Gain3 = 0.0;
    my_pv_system_3_ss_controller_B.Product2 = 0.0;
    my_pv_system_3_ss_controller_B.Integ4_h = 0.0;
    my_pv_system_3_ss_controller_B.Freq_m = 0.0;
    my_pv_system_3_ss_controller_B.Numberofsamplespercycle_m = 0.0;
    my_pv_system_3_ss_controller_B.RoundingFunction_p = 0.0;
    my_pv_system_3_ss_controller_B.Delay_d = 0.0;
    my_pv_system_3_ss_controller_B.SFunction_i = 0.0;
    my_pv_system_3_ss_controller_B.UnitDelay_l = 0.0;
    my_pv_system_3_ss_controller_B.DigitalClock_m = 0.0;
    my_pv_system_3_ss_controller_B.UnitDelay1_e = 0.0;
    my_pv_system_3_ss_controller_B.Switch_k = 0.0;
    my_pv_system_3_ss_controller_B.ComplextoMagnitudeAngle_o1 = 0.0;
    my_pv_system_3_ss_controller_B.ComplextoMagnitudeAngle_o2 = 0.0;
    my_pv_system_3_ss_controller_B.RadDeg = 0.0;
    my_pv_system_3_ss_controller_B.torad = 0.0;
    my_pv_system_3_ss_controller_B.ComplextoRealImag_o1 = 0.0;
    my_pv_system_3_ss_controller_B.ComplextoRealImag_o2 = 0.0;
    my_pv_system_3_ss_controller_B.Rff = 0.0;
    my_pv_system_3_ss_controller_B.Lff = 0.0;
    my_pv_system_3_ss_controller_B.Feedforward = 0.0;
    my_pv_system_3_ss_controller_B.Rff_e = 0.0;
    my_pv_system_3_ss_controller_B.Lff_i = 0.0;
    my_pv_system_3_ss_controller_B.Add3 = 0.0;
    my_pv_system_3_ss_controller_B.Add2[0] = 0.0;
    my_pv_system_3_ss_controller_B.Add2[1] = 0.0;
    my_pv_system_3_ss_controller_B.IntegralGain[0] = 0.0;
    my_pv_system_3_ss_controller_B.IntegralGain[1] = 0.0;
    my_pv_system_3_ss_controller_B.Saturation[0] = 0.0;
    my_pv_system_3_ss_controller_B.Saturation[1] = 0.0;
    my_pv_system_3_ss_controller_B.RateTransition3 = 0.0;
    my_pv_system_3_ss_controller_B.RateTransition4 = 0.0;
    my_pv_system_3_ss_controller_B.DigitalClock_a = 0.0;
    my_pv_system_3_ss_controller_B.RateTransition2 = 0.0;
    my_pv_system_3_ss_controller_B.Integ4_c = 0.0;
    my_pv_system_3_ss_controller_B.K1 = 0.0;
    my_pv_system_3_ss_controller_B.SFunction_pd = 0.0;
    my_pv_system_3_ss_controller_B.UnitDelay_h = 0.0;
    my_pv_system_3_ss_controller_B.UnitDelay1_p = 0.0;
    my_pv_system_3_ss_controller_B.Switch_c = 0.0;
    my_pv_system_3_ss_controller_B.TrigonometricFunction2 = 0.0;
    my_pv_system_3_ss_controller_B.Product1_c = 0.0;
    my_pv_system_3_ss_controller_B.Integ4_e = 0.0;
    my_pv_system_3_ss_controller_B.Freq_c = 0.0;
    my_pv_system_3_ss_controller_B.Numberofsamplespercycle_f = 0.0;
    my_pv_system_3_ss_controller_B.RoundingFunction_b = 0.0;
    my_pv_system_3_ss_controller_B.Delay_m = 0.0;
    my_pv_system_3_ss_controller_B.SFunction_e = 0.0;
    my_pv_system_3_ss_controller_B.UnitDelay_o = 0.0;
    my_pv_system_3_ss_controller_B.DigitalClock_e = 0.0;
    my_pv_system_3_ss_controller_B.UnitDelay1_a = 0.0;
    my_pv_system_3_ss_controller_B.Switch_d = 0.0;
    my_pv_system_3_ss_controller_B.Divide = 0.0;
    my_pv_system_3_ss_controller_B.DiscreteDerivative = 0.0;
    my_pv_system_3_ss_controller_B.DiscreteTimeIntegrator_g = 0.0;
    my_pv_system_3_ss_controller_B.Kp4 = 0.0;
    my_pv_system_3_ss_controller_B.Sum6 = 0.0;
    my_pv_system_3_ss_controller_B.Saturation1 = 0.0;
    my_pv_system_3_ss_controller_B.Gain10 = 0.0;
    my_pv_system_3_ss_controller_B.RateLimiter = 0.0;
    my_pv_system_3_ss_controller_B.x1k_d = 0.0;
    my_pv_system_3_ss_controller_B.A11 = 0.0;
    my_pv_system_3_ss_controller_B.x2k_k = 0.0;
    my_pv_system_3_ss_controller_B.A12 = 0.0;
    my_pv_system_3_ss_controller_B.A21 = 0.0;
    my_pv_system_3_ss_controller_B.A22 = 0.0;
    my_pv_system_3_ss_controller_B.sum2_c = 0.0;
    my_pv_system_3_ss_controller_B.sum3 = 0.0;
    my_pv_system_3_ss_controller_B.B11 = 0.0;
    my_pv_system_3_ss_controller_B.x1k1 = 0.0;
    my_pv_system_3_ss_controller_B.B21 = 0.0;
    my_pv_system_3_ss_controller_B.x2k1 = 0.0;
    my_pv_system_3_ss_controller_B.Duk_b = 0.0;
    my_pv_system_3_ss_controller_B.C11_d = 0.0;
    my_pv_system_3_ss_controller_B.C12_g = 0.0;
    my_pv_system_3_ss_controller_B.sum2_o = 0.0;
    my_pv_system_3_ss_controller_B.yk_e = 0.0;
    my_pv_system_3_ss_controller_B.A11_e[0] = 0.0;
    my_pv_system_3_ss_controller_B.A11_e[1] = 0.0;
    my_pv_system_3_ss_controller_B.A12_i[0] = 0.0;
    my_pv_system_3_ss_controller_B.A12_i[1] = 0.0;
    my_pv_system_3_ss_controller_B.A21_o[0] = 0.0;
    my_pv_system_3_ss_controller_B.A21_o[1] = 0.0;
    my_pv_system_3_ss_controller_B.A22_p[0] = 0.0;
    my_pv_system_3_ss_controller_B.A22_p[1] = 0.0;
    my_pv_system_3_ss_controller_B.sum2_e[0] = 0.0;
    my_pv_system_3_ss_controller_B.sum2_e[1] = 0.0;
    my_pv_system_3_ss_controller_B.sum3_g[0] = 0.0;
    my_pv_system_3_ss_controller_B.sum3_g[1] = 0.0;
    my_pv_system_3_ss_controller_B.B11_o[0] = 0.0;
    my_pv_system_3_ss_controller_B.B11_o[1] = 0.0;
    my_pv_system_3_ss_controller_B.x1k1_l[0] = 0.0;
    my_pv_system_3_ss_controller_B.x1k1_l[1] = 0.0;
    my_pv_system_3_ss_controller_B.B21_p[0] = 0.0;
    my_pv_system_3_ss_controller_B.B21_p[1] = 0.0;
    my_pv_system_3_ss_controller_B.x2k1_k[0] = 0.0;
    my_pv_system_3_ss_controller_B.x2k1_k[1] = 0.0;
    my_pv_system_3_ss_controller_B.Constant1 = 0.0;
    my_pv_system_3_ss_controller_B.Add3_n = 0.0;
    my_pv_system_3_ss_controller_B.DigitalClock_mz = 0.0;
    my_pv_system_3_ss_controller_B.Add1 = 0.0;
    my_pv_system_3_ss_controller_B.MathFunction_b = 0.0;
    my_pv_system_3_ss_controller_B.ib1 = 0.0;
    my_pv_system_3_ss_controller_B.LookupTable = 0.0;
    my_pv_system_3_ss_controller_B.Add3_o = 0.0;
    my_pv_system_3_ss_controller_B.Gain1_p = 0.0;
    my_pv_system_3_ss_controller_B.MUL1 = 0.0;
    my_pv_system_3_ss_controller_B.Add4 = 0.0;
    my_pv_system_3_ss_controller_B.UnitDelay_m = 0.0;
    my_pv_system_3_ss_controller_B.Gain_a = 0.0;
    my_pv_system_3_ss_controller_B.DataTypeConversion[0] = 0.0;
    my_pv_system_3_ss_controller_B.DataTypeConversion[1] = 0.0;
    my_pv_system_3_ss_controller_B.DataTypeConversion[2] = 0.0;
    my_pv_system_3_ss_controller_B.DataTypeConversion[3] = 0.0;
    my_pv_system_3_ss_controller_B.Switch_d2 = 0.0;
    my_pv_system_3_ss_controller_B.Add1_a = 0.0;
    my_pv_system_3_ss_controller_B.UnitDelay3[0] = 0.0;
    my_pv_system_3_ss_controller_B.UnitDelay3[1] = 0.0;
    my_pv_system_3_ss_controller_B.Gain1_g = 0.0;
    my_pv_system_3_ss_controller_B.Product = 0.0;
    my_pv_system_3_ss_controller_B.Product1_d[0] = 0.0;
    my_pv_system_3_ss_controller_B.Product1_d[1] = 0.0;
    my_pv_system_3_ss_controller_B.ComplextoMagnitudeAngle_o1_n = 0.0;
    my_pv_system_3_ss_controller_B.ComplextoMagnitudeAngle_o2_p = 0.0;
    my_pv_system_3_ss_controller_B.Add2_f = 0.0;
    my_pv_system_3_ss_controller_B.TrigonometricFunction_h = 0.0;
    my_pv_system_3_ss_controller_B.Product2_n = 0.0;
    my_pv_system_3_ss_controller_B.UnitDelay1_c = 0.0;
    my_pv_system_3_ss_controller_B.Sum_k = 0.0;
    my_pv_system_3_ss_controller_B.Rtot_pu2 = 0.0;
    my_pv_system_3_ss_controller_B.IntegralGain_a = 0.0;
    my_pv_system_3_ss_controller_B.Integrator_e = 0.0;
    my_pv_system_3_ss_controller_B.ProportionalGain_b = 0.0;
    my_pv_system_3_ss_controller_B.Sum_i = 0.0;
    my_pv_system_3_ss_controller_B.Saturate_p = 0.0;
    my_pv_system_3_ss_controller_B.Iph_p = 0.0;
    my_pv_system_3_ss_controller_B.Io_c = 0.0;
    my_pv_system_3_ss_controller_B.Fcn = 0.0;
    my_pv_system_3_ss_controller_B.Fcn1 = 0.0;
    my_pv_system_3_ss_controller_B.Fcn_l = 0.0;
    my_pv_system_3_ss_controller_B.Fcn1_f = 0.0;
    my_pv_system_3_ss_controller_B.Switch_p[0] = 0.0;
    my_pv_system_3_ss_controller_B.Switch_p[1] = 0.0;
    my_pv_system_3_ss_controller_B.Sum1 = 0.0;
    my_pv_system_3_ss_controller_B.Sum5 = 0.0;
    my_pv_system_3_ss_controller_B.Product5 = 0.0;
    my_pv_system_3_ss_controller_B.Gain1_pw = 0.0;
    my_pv_system_3_ss_controller_B.Sum4 = 0.0;
    my_pv_system_3_ss_controller_B.Product2_b = 0.0;
    my_pv_system_3_ss_controller_B.Product4 = 0.0;
    my_pv_system_3_ss_controller_B.Sum7 = 0.0;
    my_pv_system_3_ss_controller_B.Meanvalue = 0.0;
    my_pv_system_3_ss_controller_B.Sum5_b = 0.0;
    my_pv_system_3_ss_controller_B.TrigonometricFunction_p = 0.0;
    my_pv_system_3_ss_controller_B.Gain1_i = 0.0;
    my_pv_system_3_ss_controller_B.Product1_o = 0.0;
    my_pv_system_3_ss_controller_B.Integ4_d = 0.0;
    my_pv_system_3_ss_controller_B.Freq_h = 0.0;
    my_pv_system_3_ss_controller_B.Numberofsamplespercycle_h = 0.0;
    my_pv_system_3_ss_controller_B.RoundingFunction_e = 0.0;
    my_pv_system_3_ss_controller_B.Delay_db = 0.0;
    my_pv_system_3_ss_controller_B.SFunction_ij = 0.0;
    my_pv_system_3_ss_controller_B.UnitDelay_b = 0.0;
    my_pv_system_3_ss_controller_B.DigitalClock_o = 0.0;
    my_pv_system_3_ss_controller_B.UnitDelay1_c2 = 0.0;
    my_pv_system_3_ss_controller_B.Switch_j = 0.0;
    my_pv_system_3_ss_controller_B.TrigonometricFunction3_p = 0.0;
    my_pv_system_3_ss_controller_B.Gain3_c = 0.0;
    my_pv_system_3_ss_controller_B.Product2_i = 0.0;
    my_pv_system_3_ss_controller_B.Integ4_dy = 0.0;
    my_pv_system_3_ss_controller_B.Freq_mz = 0.0;
    my_pv_system_3_ss_controller_B.Numberofsamplespercycle_d = 0.0;
    my_pv_system_3_ss_controller_B.RoundingFunction_h = 0.0;
    my_pv_system_3_ss_controller_B.Delay_i = 0.0;
    my_pv_system_3_ss_controller_B.SFunction_ed = 0.0;
    my_pv_system_3_ss_controller_B.UnitDelay_fb = 0.0;
    my_pv_system_3_ss_controller_B.DigitalClock_l = 0.0;
    my_pv_system_3_ss_controller_B.UnitDelay1_j = 0.0;
    my_pv_system_3_ss_controller_B.Switch_m = 0.0;
    my_pv_system_3_ss_controller_B.ComplextoMagnitudeAngle_o1_d = 0.0;
    my_pv_system_3_ss_controller_B.ComplextoMagnitudeAngle_o2_i = 0.0;
    my_pv_system_3_ss_controller_B.RadDeg_h = 0.0;
    my_pv_system_3_ss_controller_B.Saturation_n = 0.0;
    my_pv_system_3_ss_controller_B.MathFunction_k = 0.0;
    my_pv_system_3_ss_controller_B.Sum1_c = 0.0;
    my_pv_system_3_ss_controller_B.Sum5_o = 0.0;
    my_pv_system_3_ss_controller_B.Product5_e = 0.0;
    my_pv_system_3_ss_controller_B.Gain1_n = 0.0;
    my_pv_system_3_ss_controller_B.Sum4_g = 0.0;
    my_pv_system_3_ss_controller_B.Product2_m = 0.0;
    my_pv_system_3_ss_controller_B.Product4_d = 0.0;
    my_pv_system_3_ss_controller_B.Sum7_m = 0.0;
    my_pv_system_3_ss_controller_B.Meanvalue_g = 0.0;
    my_pv_system_3_ss_controller_B.Sum5_i = 0.0;
    my_pv_system_3_ss_controller_B.Sum1_p = 0.0;
    my_pv_system_3_ss_controller_B.Sum5_k = 0.0;
    my_pv_system_3_ss_controller_B.Product5_n = 0.0;
    my_pv_system_3_ss_controller_B.Gain1_o = 0.0;
    my_pv_system_3_ss_controller_B.Sum4_n = 0.0;
    my_pv_system_3_ss_controller_B.Product2_d = 0.0;
    my_pv_system_3_ss_controller_B.Product4_o = 0.0;
    my_pv_system_3_ss_controller_B.Sum7_e = 0.0;
    my_pv_system_3_ss_controller_B.Meanvalue_h = 0.0;
    my_pv_system_3_ss_controller_B.Sum5_f = 0.0;
    my_pv_system_3_ss_controller_B.Gain1_c = 0.0;
    my_pv_system_3_ss_controller_B.Gain_f = 0.0;
    my_pv_system_3_ss_controller_B.Correction = 0.0;
    my_pv_system_3_ss_controller_B.Sum7_i = 0.0;
    my_pv_system_3_ss_controller_B.Mean = 0.0;
    my_pv_system_3_ss_controller_B.Sum5_c = 0.0;
    my_pv_system_3_ss_controller_B.Sum1_cr = 0.0;
    my_pv_system_3_ss_controller_B.Sum5_g = 0.0;
    my_pv_system_3_ss_controller_B.Product5_a = 0.0;
    my_pv_system_3_ss_controller_B.Gain1_if = 0.0;
    my_pv_system_3_ss_controller_B.Sum4_b = 0.0;
    my_pv_system_3_ss_controller_B.Product2_dp = 0.0;
    my_pv_system_3_ss_controller_B.Product4_g = 0.0;
    my_pv_system_3_ss_controller_B.Sum7_i1 = 0.0;
    my_pv_system_3_ss_controller_B.Meanvalue_e = 0.0;
    my_pv_system_3_ss_controller_B.Sum5_g2 = 0.0;
    my_pv_system_3_ss_controller_B.Sum1_c3 = 0.0;
    my_pv_system_3_ss_controller_B.Sum5_cr = 0.0;
    my_pv_system_3_ss_controller_B.Product5_j = 0.0;
    my_pv_system_3_ss_controller_B.Gain1_b = 0.0;
    my_pv_system_3_ss_controller_B.Sum4_e = 0.0;
    my_pv_system_3_ss_controller_B.Product2_f = 0.0;
    my_pv_system_3_ss_controller_B.Product4_c = 0.0;
    my_pv_system_3_ss_controller_B.Sum7_eg = 0.0;
    my_pv_system_3_ss_controller_B.Meanvalue_c = 0.0;
    my_pv_system_3_ss_controller_B.Sum5_m = 0.0;
    my_pv_system_3_ss_controller_B.TmpSignalConversionAtSFunctionI[0] = 0.0;
    my_pv_system_3_ss_controller_B.TmpSignalConversionAtSFunctionI[1] = 0.0;
    my_pv_system_3_ss_controller_B.TmpSignalConversionAtSFunctionI[2] = 0.0;
    my_pv_system_3_ss_controller_B.TmpSignalConversionAtSFunctionI[3] = 0.0;
    my_pv_system_3_ss_controller_B.D = 0.0;
  }

  /* parameters */
  my_pv_system_3_ss_controller_M->ModelData.defaultParam = ((real_T *)
    &my_pv_system_3_ss_controller_P);

  /* states (dwork) */
  my_pv_system_3_ss_controller_M->ModelData.dwork = ((void *)
    &my_pv_system_3_ss_controller_DW);
  (void) memset((void *)&my_pv_system_3_ss_controller_DW, 0,
                sizeof(DW_my_pv_system_3_ss_controller_T));
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE = 0.0;
  my_pv_system_3_ss_controller_DW.DiscreteTimeIntegrator_DSTATE = 0.0;
  my_pv_system_3_ss_controller_DW.Delay_x1_DSTATE[0] = 0.0;
  my_pv_system_3_ss_controller_DW.Delay_x1_DSTATE[1] = 0.0;
  my_pv_system_3_ss_controller_DW.Delay_x2_DSTATE[0] = 0.0;
  my_pv_system_3_ss_controller_DW.Delay_x2_DSTATE[1] = 0.0;
  my_pv_system_3_ss_controller_DW.UnitDelay2_DSTATE = 0.0;
  my_pv_system_3_ss_controller_DW.Integrator_DSTATE[0] = 0.0;
  my_pv_system_3_ss_controller_DW.Integrator_DSTATE[1] = 0.0;
  my_pv_system_3_ss_controller_DW.Integ4_DSTATE = 0.0;
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_g = 0.0;
  my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE = 0.0;
  my_pv_system_3_ss_controller_DW.Integ4_DSTATE_d = 0.0;
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_i = 0.0;
  my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_b = 0.0;
  my_pv_system_3_ss_controller_DW.Integ4_DSTATE_d5 = 0.0;
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_iq = 0.0;
  my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_a = 0.0;
  my_pv_system_3_ss_controller_DW.Integ4_DSTATE_k = 0.0;
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_a = 0.0;
  my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_i = 0.0;
  my_pv_system_3_ss_controller_DW.DiscreteDerivative_states = 0.0;
  my_pv_system_3_ss_controller_DW.DiscreteTimeIntegrator_DSTATE_d = 0.0;
  my_pv_system_3_ss_controller_DW.Delay_x1_DSTATE_d = 0.0;
  my_pv_system_3_ss_controller_DW.Delay_x2_DSTATE_m = 0.0;
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_il = 0.0;
  my_pv_system_3_ss_controller_DW.UnitDelay3_DSTATE[0] = 0.0;
  my_pv_system_3_ss_controller_DW.UnitDelay3_DSTATE[1] = 0.0;
  my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_e = 0.0;
  my_pv_system_3_ss_controller_DW.Integrator_DSTATE_d = 0.0;
  my_pv_system_3_ss_controller_DW.Integ4_DSTATE_ky = 0.0;
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_l = 0.0;
  my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_o = 0.0;
  my_pv_system_3_ss_controller_DW.Integ4_DSTATE_j = 0.0;
  my_pv_system_3_ss_controller_DW.UnitDelay_DSTATE_d = 0.0;
  my_pv_system_3_ss_controller_DW.UnitDelay1_DSTATE_k = 0.0;
  my_pv_system_3_ss_controller_DW.SFunction_PreviousInput = 0.0;
  my_pv_system_3_ss_controller_DW.Memory_1_PreviousInput[0] = 0.0;
  my_pv_system_3_ss_controller_DW.Memory_1_PreviousInput[1] = 0.0;
  my_pv_system_3_ss_controller_DW.Memory_1_PreviousInput[2] = 0.0;
  my_pv_system_3_ss_controller_DW.Memory_1_PreviousInput[3] = 0.0;
  my_pv_system_3_ss_controller_DW.Memory_2_PreviousInput = 0.0;
  my_pv_system_3_ss_controller_DW.Memory_3_PreviousInput = 0.0;
  my_pv_system_3_ss_controller_DW.DiscreteDerivative_tmp = 0.0;
  my_pv_system_3_ss_controller_DW.PrevY = 0.0;
  my_pv_system_3_ss_controller_DW.Vold = 0.0;
  my_pv_system_3_ss_controller_DW.Pold = 0.0;
  my_pv_system_3_ss_controller_DW.Dold = 0.0;
  my_pv_system_3_ss_controller_DW.SFunction_RWORK = 0.0;
  my_pv_system_3_ss_controller_DW.SFunction_RWORK_a = 0.0;
  my_pv_system_3_ss_controller_DW.SFunction_RWORK_o = 0.0;
  my_pv_system_3_ss_controller_DW.SFunction_RWORK_k = 0.0;
  my_pv_system_3_ss_controller_DW.SFunction_RWORK_j = 0.0;
  my_pv_system_3_ss_controller_DW.SFunction_RWORK_kx = 0.0;
  my_pv_system_3_ss_controller_DW.SFunction_RWORK_m = 0.0;

  /* child S-Function registration */
  {
    RTWSfcnInfo *sfcnInfo =
      &my_pv_system_3_ss_controller_M->NonInlinedSFcns.sfcnInfo;
    my_pv_system_3_ss_controller_M->sfcnInfo = (sfcnInfo);
    rtssSetErrorStatusPtr(sfcnInfo, (&rtmGetErrorStatus
      (my_pv_system_3_ss_controller_M)));
    rtssSetNumRootSampTimesPtr(sfcnInfo,
      &my_pv_system_3_ss_controller_M->Sizes.numSampTimes);
    my_pv_system_3_ss_controller_M->NonInlinedSFcns.taskTimePtrs[0] =
      &(rtmGetTPtr(my_pv_system_3_ss_controller_M)[0]);
    rtssSetTPtrPtr(sfcnInfo,
                   my_pv_system_3_ss_controller_M->NonInlinedSFcns.taskTimePtrs);
    rtssSetTStartPtr(sfcnInfo, &rtmGetTStart(my_pv_system_3_ss_controller_M));
    rtssSetTFinalPtr(sfcnInfo, &rtmGetTFinal(my_pv_system_3_ss_controller_M));
    rtssSetTimeOfLastOutputPtr(sfcnInfo, &rtmGetTimeOfLastOutput
      (my_pv_system_3_ss_controller_M));
    rtssSetStepSizePtr(sfcnInfo,
                       &my_pv_system_3_ss_controller_M->Timing.stepSize);
    rtssSetStopRequestedPtr(sfcnInfo, &rtmGetStopRequested
      (my_pv_system_3_ss_controller_M));
    rtssSetDerivCacheNeedsResetPtr(sfcnInfo,
      &my_pv_system_3_ss_controller_M->ModelData.derivCacheNeedsReset);
    rtssSetZCCacheNeedsResetPtr(sfcnInfo,
      &my_pv_system_3_ss_controller_M->ModelData.zCCacheNeedsReset);
    rtssSetBlkStateChangePtr(sfcnInfo,
      &my_pv_system_3_ss_controller_M->ModelData.blkStateChange);
    rtssSetSampleHitsPtr(sfcnInfo,
                         &my_pv_system_3_ss_controller_M->Timing.sampleHits);
    rtssSetPerTaskSampleHitsPtr(sfcnInfo,
      &my_pv_system_3_ss_controller_M->Timing.perTaskSampleHits);
    rtssSetSimModePtr(sfcnInfo, &my_pv_system_3_ss_controller_M->simMode);
    rtssSetSolverInfoPtr(sfcnInfo,
                         &my_pv_system_3_ss_controller_M->solverInfoPtr);
  }

  my_pv_system_3_ss_controller_M->Sizes.numSFcns = (11);

  /* register each child */
  {
    (void) memset((void *)
                  &my_pv_system_3_ss_controller_M->NonInlinedSFcns.childSFunctions
                  [0], 0,
                  11*sizeof(SimStruct));
    my_pv_system_3_ss_controller_M->childSfunctions =
      (&my_pv_system_3_ss_controller_M->NonInlinedSFcns.childSFunctionPtrs[0]);

    {
      int_T i;
      for (i = 0; i < 11; i++) {
        my_pv_system_3_ss_controller_M->childSfunctions[i] =
          (&my_pv_system_3_ss_controller_M->NonInlinedSFcns.childSFunctions[i]);
      }
    }

    /* Level2 S-Function Block: my_pv_system_3_ss_controller/<S41>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[0];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn0.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn0.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn0.sfcnTsMap;
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
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.blkInfo2
                         [0]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_3_ss_controller_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods2
                           [0]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods3
                           [0]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.statesInfo2
                         [0]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn0.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn0.UPtrs0;
          sfcnUPtrs[0] = &my_pv_system_3_ss_controller_B.Integ4_d;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn0.UPtrs1;
          sfcnUPtrs[0] = &my_pv_system_3_ss_controller_B.Delay_db;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn0.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &my_pv_system_3_ss_controller_B.SFunction_ij));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,my_pv_system_3_ss_controller_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn0.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P4_Size);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *)
                 &my_pv_system_3_ss_controller_DW.SFunction_RWORK_kx);
      ssSetIWork(rts, (int_T *)
                 &my_pv_system_3_ss_controller_DW.SFunction_IWORK_pd);
      ssSetPWork(rts, (void **)
                 &my_pv_system_3_ss_controller_DW.SFunction_PWORK_n);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn0.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn0.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_3_ss_controller_DW.SFunction_RWORK_kx);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &my_pv_system_3_ss_controller_DW.SFunction_IWORK_pd);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &my_pv_system_3_ss_controller_DW.SFunction_PWORK_n);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

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

    /* Level2 S-Function Block: my_pv_system_3_ss_controller/<S44>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[1];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn1.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn1.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn1.sfcnTsMap;
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
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.blkInfo2
                         [1]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_3_ss_controller_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods2
                           [1]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods3
                           [1]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.statesInfo2
                         [1]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn1.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn1.UPtrs0;
          sfcnUPtrs[0] = &my_pv_system_3_ss_controller_B.Integ4_dy;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn1.UPtrs1;
          sfcnUPtrs[0] = &my_pv_system_3_ss_controller_B.Delay_i;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn1.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &my_pv_system_3_ss_controller_B.SFunction_ed));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,my_pv_system_3_ss_controller_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn1.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P1_Size_a);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P2_Size_k);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P3_Size_h);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P4_Size_l);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *)
                 &my_pv_system_3_ss_controller_DW.SFunction_RWORK_m);
      ssSetIWork(rts, (int_T *)
                 &my_pv_system_3_ss_controller_DW.SFunction_IWORK_c4);
      ssSetPWork(rts, (void **)
                 &my_pv_system_3_ss_controller_DW.SFunction_PWORK_p);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn1.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn1.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_3_ss_controller_DW.SFunction_RWORK_m);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &my_pv_system_3_ss_controller_DW.SFunction_IWORK_c4);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &my_pv_system_3_ss_controller_DW.SFunction_PWORK_p);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

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

    /* Level2 S-Function Block: my_pv_system_3_ss_controller/<S90>/S-Function (send_rt) */
    {
      SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[2];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn2.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn2.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn2.sfcnTsMap;
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
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.blkInfo2
                         [2]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_3_ss_controller_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods2
                           [2]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods3
                           [2]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.statesInfo2
                         [2]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn2.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn2.UPtrs0;
          sfcnUPtrs[0] = &my_pv_system_3_ss_controller_B.Gate[0];
          sfcnUPtrs[1] = &my_pv_system_3_ss_controller_B.Gate[1];
          sfcnUPtrs[2] = &my_pv_system_3_ss_controller_B.Gate[2];
          sfcnUPtrs[3] = &my_pv_system_3_ss_controller_B.Gate[3];
          sfcnUPtrs[4] = &my_pv_system_3_ss_controller_B.Iph;
          sfcnUPtrs[5] = &my_pv_system_3_ss_controller_B.Io;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 6);
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "my_pv_system_3_ss_controller/SS_controller/rtlab_send_subsystem/Subsystem4/Send4/S-Function");
      ssSetRTModel(rts,my_pv_system_3_ss_controller_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn2.params;
        ssSetSFcnParamsCount(rts, 3);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P1_Size_e);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P2_Size_p);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P3_Size_e);
      }

      /* work vectors */
      ssSetIWork(rts, (int_T *)
                 &my_pv_system_3_ss_controller_DW.SFunction_IWORK_i[0]);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn2.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn2.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* IWORK */
        ssSetDWorkWidth(rts, 0, 5);
        ssSetDWorkDataType(rts, 0,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_3_ss_controller_DW.SFunction_IWORK_i[0]);
      }

      /* registration */
      send_rt(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 6);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: my_pv_system_3_ss_controller/<S2>/OpMonitor (opmonitor) */
    {
      SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[3];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn3.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn3.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn3.sfcnTsMap;
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
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.blkInfo2
                         [3]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_3_ss_controller_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods2
                           [3]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods3
                           [3]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.statesInfo2
                         [3]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn3.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn3.UPtrs0;
          sfcnUPtrs[0] = (real_T*)&my_pv_system_3_ss_controller_RGND;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn3.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 4);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &my_pv_system_3_ss_controller_B.Computationtime));
        }

        /* port 1 */
        {
          _ssSetOutputPortNumDimensions(rts, 1, 1);
          ssSetOutputPortWidth(rts, 1, 1);
          ssSetOutputPortSignal(rts, 1, ((real_T *)
            &my_pv_system_3_ss_controller_B.Realstepsize));
        }

        /* port 2 */
        {
          _ssSetOutputPortNumDimensions(rts, 2, 1);
          ssSetOutputPortWidth(rts, 2, 1);
          ssSetOutputPortSignal(rts, 2, ((real_T *)
            &my_pv_system_3_ss_controller_B.Idletime));
        }

        /* port 3 */
        {
          _ssSetOutputPortNumDimensions(rts, 3, 1);
          ssSetOutputPortWidth(rts, 3, 1);
          ssSetOutputPortSignal(rts, 3, ((real_T *)
            &my_pv_system_3_ss_controller_B.Overruntimes));
        }
      }

      /* path info */
      ssSetModelName(rts, "OpMonitor");
      ssSetPath(rts, "my_pv_system_3_ss_controller/SS_controller/OpMonitor");
      ssSetRTModel(rts,my_pv_system_3_ss_controller_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn3.params;
        ssSetSFcnParamsCount(rts, 6);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_3_ss_controller_P.OpMonitor_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_3_ss_controller_P.OpMonitor_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       my_pv_system_3_ss_controller_P.OpMonitor_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       my_pv_system_3_ss_controller_P.OpMonitor_P4_Size);
        ssSetSFcnParam(rts, 4, (mxArray*)
                       my_pv_system_3_ss_controller_P.OpMonitor_P5_Size);
        ssSetSFcnParam(rts, 5, (mxArray*)
                       my_pv_system_3_ss_controller_P.OpMonitor_P6_Size);
      }

      /* work vectors */
      ssSetPWork(rts, (void **) &my_pv_system_3_ss_controller_DW.OpMonitor_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn3.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn3.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_3_ss_controller_DW.OpMonitor_PWORK);
      }

      /* registration */
      opmonitor(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

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

    /* Level2 S-Function Block: my_pv_system_3_ss_controller/<S89>/S-Function (OP_SEND) */
    {
      SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[4];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn4.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn4.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn4.sfcnTsMap;
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
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.blkInfo2
                         [4]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_3_ss_controller_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods2
                           [4]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods3
                           [4]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.statesInfo2
                         [4]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn4.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn4.UPtrs0;
          sfcnUPtrs[0] = &my_pv_system_3_ss_controller_B.Computationtime;
          sfcnUPtrs[1] = &my_pv_system_3_ss_controller_B.Realstepsize;
          sfcnUPtrs[2] = &my_pv_system_3_ss_controller_B.Idletime;
          sfcnUPtrs[3] = &my_pv_system_3_ss_controller_B.Overruntimes;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 4);
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "my_pv_system_3_ss_controller/SS_controller/rtlab_send_subsystem/Subsystem2/Send2/S-Function");
      ssSetRTModel(rts,my_pv_system_3_ss_controller_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn4.params;
        ssSetSFcnParamsCount(rts, 1);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P1_Size_h);
      }

      /* work vectors */
      ssSetIWork(rts, (int_T *)
                 &my_pv_system_3_ss_controller_DW.SFunction_IWORK[0]);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn4.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn4.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* IWORK */
        ssSetDWorkWidth(rts, 0, 5);
        ssSetDWorkDataType(rts, 0,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_3_ss_controller_DW.SFunction_IWORK[0]);
      }

      /* registration */
      OP_SEND(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

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

    /* Level2 S-Function Block: my_pv_system_3_ss_controller/<S91>/S-Function (recv_rt) */
    {
      SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[5];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn5.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn5.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn5.sfcnTsMap;
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
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.blkInfo2
                         [5]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_3_ss_controller_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods2
                           [5]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods3
                           [5]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.statesInfo2
                         [5]);
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn5.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 12);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            my_pv_system_3_ss_controller_B.SFunction_g));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "my_pv_system_3_ss_controller/SS_controller/zzzOpComm/Receive_1/S-Function");
      ssSetRTModel(rts,my_pv_system_3_ss_controller_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn5.params;
        ssSetSFcnParamsCount(rts, 3);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P1_Size_o);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P2_Size_m);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P3_Size_b);
      }

      /* work vectors */
      ssSetPWork(rts, (void **) &my_pv_system_3_ss_controller_DW.SFunction_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn5.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn5.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_3_ss_controller_DW.SFunction_PWORK);
      }

      /* registration */
      recv_rt(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
    }

    /* Level2 S-Function Block: my_pv_system_3_ss_controller/<S57>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[6];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn6.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn6.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn6.sfcnTsMap;
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
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.blkInfo2
                         [6]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_3_ss_controller_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods2
                           [6]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods3
                           [6]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.statesInfo2
                         [6]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn6.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn6.UPtrs0;
          sfcnUPtrs[0] = &my_pv_system_3_ss_controller_B.Apu;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn6.UPtrs1;
          sfcnUPtrs[0] = &my_pv_system_3_ss_controller_B.Gain;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn6.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &my_pv_system_3_ss_controller_B.SFunction_o));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,my_pv_system_3_ss_controller_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn6.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P1_Size_f);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P2_Size_j);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P3_Size_h5);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P4_Size_o);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *)
                 &my_pv_system_3_ss_controller_DW.SFunction_RWORK);
      ssSetIWork(rts, (int_T *)
                 &my_pv_system_3_ss_controller_DW.SFunction_IWORK_p);
      ssSetPWork(rts, (void **)
                 &my_pv_system_3_ss_controller_DW.SFunction_PWORK_a);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn6.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn6.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_3_ss_controller_DW.SFunction_RWORK);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &my_pv_system_3_ss_controller_DW.SFunction_IWORK_p);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &my_pv_system_3_ss_controller_DW.SFunction_PWORK_a);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

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

    /* Level2 S-Function Block: my_pv_system_3_ss_controller/<S25>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[7];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn7.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn7.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn7.sfcnTsMap;
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
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.blkInfo2
                         [7]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_3_ss_controller_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods2
                           [7]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods3
                           [7]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.statesInfo2
                         [7]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn7.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn7.UPtrs0;
          sfcnUPtrs[0] = &my_pv_system_3_ss_controller_B.Integ4;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn7.UPtrs1;
          sfcnUPtrs[0] = &my_pv_system_3_ss_controller_B.Delay;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn7.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &my_pv_system_3_ss_controller_B.SFunction_p));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,my_pv_system_3_ss_controller_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn7.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P1_Size_fw);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P2_Size_a);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P3_Size_l);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P4_Size_n);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *)
                 &my_pv_system_3_ss_controller_DW.SFunction_RWORK_a);
      ssSetIWork(rts, (int_T *)
                 &my_pv_system_3_ss_controller_DW.SFunction_IWORK_f);
      ssSetPWork(rts, (void **)
                 &my_pv_system_3_ss_controller_DW.SFunction_PWORK_c);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn7.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn7.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_3_ss_controller_DW.SFunction_RWORK_a);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &my_pv_system_3_ss_controller_DW.SFunction_IWORK_f);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &my_pv_system_3_ss_controller_DW.SFunction_PWORK_c);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

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

    /* Level2 S-Function Block: my_pv_system_3_ss_controller/<S28>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[8];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn8.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn8.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn8.sfcnTsMap;
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
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.blkInfo2
                         [8]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_3_ss_controller_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods2
                           [8]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods3
                           [8]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.statesInfo2
                         [8]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn8.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn8.UPtrs0;
          sfcnUPtrs[0] = &my_pv_system_3_ss_controller_B.Integ4_h;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn8.UPtrs1;
          sfcnUPtrs[0] = &my_pv_system_3_ss_controller_B.Delay_d;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn8.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &my_pv_system_3_ss_controller_B.SFunction_i));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,my_pv_system_3_ss_controller_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn8.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P1_Size_l);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P2_Size_f);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P3_Size_f);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P4_Size_c);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *)
                 &my_pv_system_3_ss_controller_DW.SFunction_RWORK_o);
      ssSetIWork(rts, (int_T *)
                 &my_pv_system_3_ss_controller_DW.SFunction_IWORK_pq);
      ssSetPWork(rts, (void **)
                 &my_pv_system_3_ss_controller_DW.SFunction_PWORK_f);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn8.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn8.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_3_ss_controller_DW.SFunction_RWORK_o);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &my_pv_system_3_ss_controller_DW.SFunction_IWORK_pq);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &my_pv_system_3_ss_controller_DW.SFunction_PWORK_f);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

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

    /* Level2 S-Function Block: my_pv_system_3_ss_controller/<S30>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[9];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn9.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn9.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn9.sfcnTsMap;
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
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.blkInfo2
                         [9]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_3_ss_controller_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods2
                           [9]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods3
                           [9]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.statesInfo2
                         [9]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn9.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn9.UPtrs0;
          sfcnUPtrs[0] = &my_pv_system_3_ss_controller_B.Integ4_c;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn9.UPtrs1;
          sfcnUPtrs[0] = &my_pv_system_3_ss_controller_B.K1;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn9.outputPortInfo
          [0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &my_pv_system_3_ss_controller_B.SFunction_pd));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Mean/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,my_pv_system_3_ss_controller_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn9.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P1_Size_p);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P2_Size_ar);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P3_Size_ep);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P4_Size_c0);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *)
                 &my_pv_system_3_ss_controller_DW.SFunction_RWORK_k);
      ssSetIWork(rts, (int_T *)
                 &my_pv_system_3_ss_controller_DW.SFunction_IWORK_b);
      ssSetPWork(rts, (void **)
                 &my_pv_system_3_ss_controller_DW.SFunction_PWORK_i);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn9.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn9.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_3_ss_controller_DW.SFunction_RWORK_k);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &my_pv_system_3_ss_controller_DW.SFunction_IWORK_b);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &my_pv_system_3_ss_controller_DW.SFunction_PWORK_i);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

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

    /* Level2 S-Function Block: my_pv_system_3_ss_controller/<S47>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = my_pv_system_3_ss_controller_M->childSfunctions[10];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn10.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn10.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn10.sfcnTsMap;
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
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.blkInfo2
                         [10]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_3_ss_controller_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods2
                           [10]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_3_ss_controller_M->NonInlinedSFcns.methods3
                           [10]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_3_ss_controller_M->NonInlinedSFcns.statesInfo2
                         [10]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn10.inputPortInfo
          [0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn10.UPtrs0;
          sfcnUPtrs[0] = &my_pv_system_3_ss_controller_B.Integ4_e;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn10.UPtrs1;
          sfcnUPtrs[0] = &my_pv_system_3_ss_controller_B.Delay_m;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &my_pv_system_3_ss_controller_M->
          NonInlinedSFcns.Sfcn10.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &my_pv_system_3_ss_controller_B.SFunction_e));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Mean (Variable Frequency)/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,my_pv_system_3_ss_controller_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn10.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P1_Size_ol);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P2_Size_h);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P3_Size_j);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       my_pv_system_3_ss_controller_P.SFunction_P4_Size_d);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *)
                 &my_pv_system_3_ss_controller_DW.SFunction_RWORK_j);
      ssSetIWork(rts, (int_T *)
                 &my_pv_system_3_ss_controller_DW.SFunction_IWORK_c);
      ssSetPWork(rts, (void **)
                 &my_pv_system_3_ss_controller_DW.SFunction_PWORK_a5);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn10.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_3_ss_controller_M->NonInlinedSFcns.Sfcn10.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_3_ss_controller_DW.SFunction_RWORK_j);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &my_pv_system_3_ss_controller_DW.SFunction_IWORK_c);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &my_pv_system_3_ss_controller_DW.SFunction_PWORK_a5);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

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
  }

  /* Initialize Sizes */
  my_pv_system_3_ss_controller_M->Sizes.numContStates = (0);/* Number of continuous states */
  my_pv_system_3_ss_controller_M->Sizes.numY = (0);/* Number of model outputs */
  my_pv_system_3_ss_controller_M->Sizes.numU = (0);/* Number of model inputs */
  my_pv_system_3_ss_controller_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  my_pv_system_3_ss_controller_M->Sizes.numSampTimes = (1);/* Number of sample times */
  my_pv_system_3_ss_controller_M->Sizes.numBlocks = (320);/* Number of blocks */
  my_pv_system_3_ss_controller_M->Sizes.numBlockIO = (278);/* Number of block outputs */
  my_pv_system_3_ss_controller_M->Sizes.numBlockPrms = (329);/* Sum of parameter "widths" */
  return my_pv_system_3_ss_controller_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
