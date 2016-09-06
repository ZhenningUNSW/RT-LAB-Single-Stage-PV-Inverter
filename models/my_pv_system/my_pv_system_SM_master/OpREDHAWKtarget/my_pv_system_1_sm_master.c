/*
 * my_pv_system_1_sm_master.c
 *
 * Code generation for model "my_pv_system_1_sm_master".
 *
 * Model version              : 1.177
 * Simulink Coder version : 8.7 (R2014b) 08-Sep-2014
 * C source code generated on : Tue Sep 06 10:33:15 2016
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

/* Block states (auto storage) */
DW_my_pv_system_1_sm_master_T my_pv_system_1_sm_master_DW;

/* Real-time model */
RT_MODEL_my_pv_system_1_sm_master_T my_pv_system_1_sm_master_M_;
RT_MODEL_my_pv_system_1_sm_master_T *const my_pv_system_1_sm_master_M =
  &my_pv_system_1_sm_master_M_;
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
static void my_pv_system_1_sm_master_output(void)
{
  boolean_T y;
  real_T y_0;
  real_T u1;
  real_T u2;

  /* Memory: '<S1>/S-Function' */
  my_pv_system_1_sm_master_B.SFunction =
    my_pv_system_1_sm_master_DW.SFunction_PreviousInput;

  /* Sum: '<S1>/Sum' incorporates:
   *  Constant: '<S1>/S-Function1'
   */
  my_pv_system_1_sm_master_B.Sum = my_pv_system_1_sm_master_P.SFunction1_Value +
    my_pv_system_1_sm_master_B.SFunction;

  /* Stop: '<S1>/Stop Simulation' */
  if (my_pv_system_1_sm_master_B.Sum != 0.0) {
    rtmSetStopRequested(my_pv_system_1_sm_master_M, 1);
  }

  /* End of Stop: '<S1>/Stop Simulation' */

  /* Memory: '<S2>/Memory' */
  my_pv_system_1_sm_master_B.IrradianceWm2 =
    my_pv_system_1_sm_master_DW.Memory_1_PreviousInput;

  /* Memory: '<S2>/Memory' */
  my_pv_system_1_sm_master_B.I_PV =
    my_pv_system_1_sm_master_DW.Memory_2_PreviousInput;

  /* Memory: '<S2>/Memory' */
  my_pv_system_1_sm_master_B.V_PV =
    my_pv_system_1_sm_master_DW.Memory_3_PreviousInput;

  /* Memory: '<S2>/Memory' */
  my_pv_system_1_sm_master_B.I_Diode =
    my_pv_system_1_sm_master_DW.Memory_4_PreviousInput;

  /* Memory: '<S2>/Memory1' */
  my_pv_system_1_sm_master_B.Vinv =
    my_pv_system_1_sm_master_DW.Memory1_1_PreviousInput;

  /* Memory: '<S2>/Memory1' */
  my_pv_system_1_sm_master_B.Vdc =
    my_pv_system_1_sm_master_DW.Memory1_2_PreviousInput;

  /* Memory: '<S2>/Memory1' */
  my_pv_system_1_sm_master_B.Igrid =
    my_pv_system_1_sm_master_DW.Memory1_3_PreviousInput;

  /* Memory: '<S2>/Memory1' */
  my_pv_system_1_sm_master_B.Vgrid =
    my_pv_system_1_sm_master_DW.Memory1_4_PreviousInput;

  /* Memory: '<S2>/Memory1' */
  my_pv_system_1_sm_master_B.Ig =
    my_pv_system_1_sm_master_DW.Memory1_5_PreviousInput;

  /* Outputs for Atomic SubSystem: '<S13>/Subsystem5' */

  /* Level2 S-Function Block: '<S46>/S-Function' (send_rt) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[0];
    sfcnOutputs(rts, 1);
  }

  /* End of Outputs for SubSystem: '<S13>/Subsystem5' */

  /* Level2 S-Function Block: '<S44>/S-Function' (OP_SEND) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[1];
    sfcnOutputs(rts, 1);
  }

  /* Level2 S-Function Block: '<S2>/OpMonitor' (opmonitor) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[2];
    sfcnOutputs(rts, 1);
  }

  /* Level2 S-Function Block: '<S45>/S-Function' (OP_SEND) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[3];
    sfcnOutputs(rts, 1);
  }

  /* UnitDelay: '<S21>/Unit Delay' */
  my_pv_system_1_sm_master_B.UnitDelay[0] =
    my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[0];
  my_pv_system_1_sm_master_B.UnitDelay[1] =
    my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[1];
  my_pv_system_1_sm_master_B.UnitDelay[2] =
    my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[2];
  my_pv_system_1_sm_master_B.UnitDelay[3] =
    my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[3];

  /* Level2 S-Function Block: '<S47>/S-Function' (recv_rt) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[4];
    sfcnOutputs(rts, 1);
  }

  /* Outputs for Enabled SubSystem: '<S21>/Tail' incorporates:
   *  EnablePort: '<S22>/Enable'
   */
  /* Constant: '<S21>/2' */
  if (my_pv_system_1_sm_master_P._Value_g) {
    if (!my_pv_system_1_sm_master_DW.Tail_MODE) {
      /* InitializeConditions for DiscreteIntegrator: '<S22>/Discrete-Time Integrator' */
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

      /* InitializeConditions for UnitDelay: '<S22>/Unit Delay' */
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

    /* Constant: '<S22>/2' */
    y = (my_pv_system_1_sm_master_P.Tail_Tf + my_pv_system_1_sm_master_P.Tail_Tt
         > 0.0);
    my_pv_system_1_sm_master_B.u = y;

    /* DiscreteIntegrator: '<S22>/Discrete-Time Integrator' */
    if ((my_pv_system_1_sm_master_B.SFunction_p[0] <= 0.0) &&
        (my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[0] == 1)) {
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[0] =
        my_pv_system_1_sm_master_P.DiscreteTimeIntegrator_IC;
    }

    if ((my_pv_system_1_sm_master_B.SFunction_p[1] <= 0.0) &&
        (my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[1] == 1)) {
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[1] =
        my_pv_system_1_sm_master_P.DiscreteTimeIntegrator_IC;
    }

    if ((my_pv_system_1_sm_master_B.SFunction_p[2] <= 0.0) &&
        (my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[2] == 1)) {
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[2] =
        my_pv_system_1_sm_master_P.DiscreteTimeIntegrator_IC;
    }

    if ((my_pv_system_1_sm_master_B.SFunction_p[3] <= 0.0) &&
        (my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[3] == 1)) {
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

    /* End of DiscreteIntegrator: '<S22>/Discrete-Time Integrator' */

    /* Gain: '<S24>/-0.9//Tf' */
    y_0 = my_pv_system_1_sm_master_P.Tail_Tf + 2.2204460492503131E-16;
    y_0 = -0.9 / y_0;
    my_pv_system_1_sm_master_B.u9Tf[0] = y_0 *
      my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[0];
    my_pv_system_1_sm_master_B.u9Tf[1] = y_0 *
      my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[1];
    my_pv_system_1_sm_master_B.u9Tf[2] = y_0 *
      my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[2];
    my_pv_system_1_sm_master_B.u9Tf[3] = y_0 *
      my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[3];

    /* Sum: '<S24>/Add' incorporates:
     *  Constant: '<S24>/Constant'
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

    /* Saturate: '<S24>/Saturation1' */
    y_0 = my_pv_system_1_sm_master_B.Add[0];
    u1 = my_pv_system_1_sm_master_P.Saturation1_LowerSat;
    u2 = my_pv_system_1_sm_master_P.Saturation1_UpperSat;
    if (y_0 > u2) {
      y_0 = u2;
    } else {
      if (y_0 < u1) {
        y_0 = u1;
      }
    }

    my_pv_system_1_sm_master_B.Saturation1[0] = y_0;
    y_0 = my_pv_system_1_sm_master_B.Add[1];
    u1 = my_pv_system_1_sm_master_P.Saturation1_LowerSat;
    u2 = my_pv_system_1_sm_master_P.Saturation1_UpperSat;
    if (y_0 > u2) {
      y_0 = u2;
    } else {
      if (y_0 < u1) {
        y_0 = u1;
      }
    }

    my_pv_system_1_sm_master_B.Saturation1[1] = y_0;
    y_0 = my_pv_system_1_sm_master_B.Add[2];
    u1 = my_pv_system_1_sm_master_P.Saturation1_LowerSat;
    u2 = my_pv_system_1_sm_master_P.Saturation1_UpperSat;
    if (y_0 > u2) {
      y_0 = u2;
    } else {
      if (y_0 < u1) {
        y_0 = u1;
      }
    }

    my_pv_system_1_sm_master_B.Saturation1[2] = y_0;
    y_0 = my_pv_system_1_sm_master_B.Add[3];
    u1 = my_pv_system_1_sm_master_P.Saturation1_LowerSat;
    u2 = my_pv_system_1_sm_master_P.Saturation1_UpperSat;
    if (y_0 > u2) {
      y_0 = u2;
    } else {
      if (y_0 < u1) {
        y_0 = u1;
      }
    }

    my_pv_system_1_sm_master_B.Saturation1[3] = y_0;

    /* End of Saturate: '<S24>/Saturation1' */

    /* Sum: '<S24>/Add1' incorporates:
     *  Constant: '<S24>/Constant2'
     */
    y_0 = my_pv_system_1_sm_master_P.Tail_Tf +
      my_pv_system_1_sm_master_P.Tail_Tt;
    my_pv_system_1_sm_master_B.Add1[0] = y_0 -
      my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[0];
    my_pv_system_1_sm_master_B.Add1[1] = y_0 -
      my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[1];
    my_pv_system_1_sm_master_B.Add1[2] = y_0 -
      my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[2];
    my_pv_system_1_sm_master_B.Add1[3] = y_0 -
      my_pv_system_1_sm_master_B.DiscreteTimeIntegrator[3];

    /* Gain: '<S24>/0.1//Tt' */
    y_0 = my_pv_system_1_sm_master_P.Tail_Tt + 2.2204460492503131E-16;
    y_0 = 0.1 / y_0;
    my_pv_system_1_sm_master_B.uTt[0] = y_0 * my_pv_system_1_sm_master_B.Add1[0];
    my_pv_system_1_sm_master_B.uTt[1] = y_0 * my_pv_system_1_sm_master_B.Add1[1];
    my_pv_system_1_sm_master_B.uTt[2] = y_0 * my_pv_system_1_sm_master_B.Add1[2];
    my_pv_system_1_sm_master_B.uTt[3] = y_0 * my_pv_system_1_sm_master_B.Add1[3];

    /* Saturate: '<S24>/Saturation2' */
    y_0 = my_pv_system_1_sm_master_B.uTt[0];
    u1 = my_pv_system_1_sm_master_P.Saturation2_LowerSat;
    u2 = my_pv_system_1_sm_master_P.Saturation2_UpperSat;
    if (y_0 > u2) {
      y_0 = u2;
    } else {
      if (y_0 < u1) {
        y_0 = u1;
      }
    }

    my_pv_system_1_sm_master_B.Saturation2[0] = y_0;
    y_0 = my_pv_system_1_sm_master_B.uTt[1];
    u1 = my_pv_system_1_sm_master_P.Saturation2_LowerSat;
    u2 = my_pv_system_1_sm_master_P.Saturation2_UpperSat;
    if (y_0 > u2) {
      y_0 = u2;
    } else {
      if (y_0 < u1) {
        y_0 = u1;
      }
    }

    my_pv_system_1_sm_master_B.Saturation2[1] = y_0;
    y_0 = my_pv_system_1_sm_master_B.uTt[2];
    u1 = my_pv_system_1_sm_master_P.Saturation2_LowerSat;
    u2 = my_pv_system_1_sm_master_P.Saturation2_UpperSat;
    if (y_0 > u2) {
      y_0 = u2;
    } else {
      if (y_0 < u1) {
        y_0 = u1;
      }
    }

    my_pv_system_1_sm_master_B.Saturation2[2] = y_0;
    y_0 = my_pv_system_1_sm_master_B.uTt[3];
    u1 = my_pv_system_1_sm_master_P.Saturation2_LowerSat;
    u2 = my_pv_system_1_sm_master_P.Saturation2_UpperSat;
    if (y_0 > u2) {
      y_0 = u2;
    } else {
      if (y_0 < u1) {
        y_0 = u1;
      }
    }

    my_pv_system_1_sm_master_B.Saturation2[3] = y_0;

    /* End of Saturate: '<S24>/Saturation2' */

    /* Sum: '<S24>/Add2' */
    my_pv_system_1_sm_master_B.Add2[0] = my_pv_system_1_sm_master_B.Saturation1
      [0] + my_pv_system_1_sm_master_B.Saturation2[0];
    my_pv_system_1_sm_master_B.Add2[1] = my_pv_system_1_sm_master_B.Saturation1
      [1] + my_pv_system_1_sm_master_B.Saturation2[1];
    my_pv_system_1_sm_master_B.Add2[2] = my_pv_system_1_sm_master_B.Saturation1
      [2] + my_pv_system_1_sm_master_B.Saturation2[2];
    my_pv_system_1_sm_master_B.Add2[3] = my_pv_system_1_sm_master_B.Saturation1
      [3] + my_pv_system_1_sm_master_B.Saturation2[3];

    /* UnitDelay: '<S22>/Unit Delay' */
    my_pv_system_1_sm_master_B.UnitDelay_l[0] =
      my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[0];
    my_pv_system_1_sm_master_B.UnitDelay_l[1] =
      my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[1];
    my_pv_system_1_sm_master_B.UnitDelay_l[2] =
      my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[2];
    my_pv_system_1_sm_master_B.UnitDelay_l[3] =
      my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[3];

    /* Switch: '<S22>/Switch' */
    if (my_pv_system_1_sm_master_B.SFunction_p[0] >=
        my_pv_system_1_sm_master_P.Switch_Threshold) {
      my_pv_system_1_sm_master_B.Switch_j[0] =
        my_pv_system_1_sm_master_B.UnitDelay[0];
    } else {
      my_pv_system_1_sm_master_B.Switch_j[0] =
        my_pv_system_1_sm_master_B.UnitDelay_l[0];
    }

    if (my_pv_system_1_sm_master_B.SFunction_p[1] >=
        my_pv_system_1_sm_master_P.Switch_Threshold) {
      my_pv_system_1_sm_master_B.Switch_j[1] =
        my_pv_system_1_sm_master_B.UnitDelay[1];
    } else {
      my_pv_system_1_sm_master_B.Switch_j[1] =
        my_pv_system_1_sm_master_B.UnitDelay_l[1];
    }

    if (my_pv_system_1_sm_master_B.SFunction_p[2] >=
        my_pv_system_1_sm_master_P.Switch_Threshold) {
      my_pv_system_1_sm_master_B.Switch_j[2] =
        my_pv_system_1_sm_master_B.UnitDelay[2];
    } else {
      my_pv_system_1_sm_master_B.Switch_j[2] =
        my_pv_system_1_sm_master_B.UnitDelay_l[2];
    }

    if (my_pv_system_1_sm_master_B.SFunction_p[3] >=
        my_pv_system_1_sm_master_P.Switch_Threshold) {
      my_pv_system_1_sm_master_B.Switch_j[3] =
        my_pv_system_1_sm_master_B.UnitDelay[3];
    } else {
      my_pv_system_1_sm_master_B.Switch_j[3] =
        my_pv_system_1_sm_master_B.UnitDelay_l[3];
    }

    /* End of Switch: '<S22>/Switch' */

    /* Product: '<S22>/Product' */
    my_pv_system_1_sm_master_B.Product[0] = my_pv_system_1_sm_master_B.Add2[0] *
      my_pv_system_1_sm_master_B.Switch_j[0] * my_pv_system_1_sm_master_B.u;
    my_pv_system_1_sm_master_B.Product[1] = my_pv_system_1_sm_master_B.Add2[1] *
      my_pv_system_1_sm_master_B.Switch_j[1] * my_pv_system_1_sm_master_B.u;
    my_pv_system_1_sm_master_B.Product[2] = my_pv_system_1_sm_master_B.Add2[2] *
      my_pv_system_1_sm_master_B.Switch_j[2] * my_pv_system_1_sm_master_B.u;
    my_pv_system_1_sm_master_B.Product[3] = my_pv_system_1_sm_master_B.Add2[3] *
      my_pv_system_1_sm_master_B.Switch_j[3] * my_pv_system_1_sm_master_B.u;
  } else {
    if (my_pv_system_1_sm_master_DW.Tail_MODE) {
      my_pv_system_1_sm_master_DW.Tail_MODE = false;
    }
  }

  /* End of Constant: '<S21>/2' */
  /* End of Outputs for SubSystem: '<S21>/Tail' */

  /* Sin: '<S40>/AC' */
  u1 = my_pv_system_1_sm_master_P.Vs14400V_Phase * 3.1415926535897931;
  y_0 = u1 / 180.0;
  u1 = 6.2831853071795862 * my_pv_system_1_sm_master_P.Vs14400V_Frequency;
  my_pv_system_1_sm_master_B.AC = sin(u1 * my_pv_system_1_sm_master_M->Timing.t
    [0] + y_0) * my_pv_system_1_sm_master_P.Vs14400V_Amplitude +
    my_pv_system_1_sm_master_P.AC_Bias;

  /* Delay: '<S9>/Delay2' */
  my_pv_system_1_sm_master_B.I_Diode_k =
    my_pv_system_1_sm_master_DW.Delay2_DSTATE;

  /* Delay: '<S9>/Delay1' */
  my_pv_system_1_sm_master_B.Delay1 = my_pv_system_1_sm_master_DW.Delay1_DSTATE;

  /* S-Function block: <S50>/State-Space */
  {
    real_T accum;

    /* Circuit has switches */
    int_T *switch_status = (int_T*)
      my_pv_system_1_sm_master_DW.StateSpace_PWORK.SWITCH_STATUS;
    int_T *switch_status_init = (int_T*)
      my_pv_system_1_sm_master_DW.StateSpace_PWORK.SWITCH_STATUS_INIT;
    int_T *SwitchChange = (int_T*)
      my_pv_system_1_sm_master_DW.StateSpace_PWORK.SW_CHG;
    int_T *Chopper = (int_T*)
      my_pv_system_1_sm_master_DW.StateSpace_PWORK.CHOPPER;
    int_T *gState = (int_T*)
      my_pv_system_1_sm_master_DW.StateSpace_PWORK.G_STATE;
    real_T *yswitch = (real_T*)
      my_pv_system_1_sm_master_DW.StateSpace_PWORK.Y_SWITCH;
    int_T *switchTypes = (int_T*)
      my_pv_system_1_sm_master_DW.StateSpace_PWORK.SWITCH_TYPES;
    int_T *idxOutSw = (int_T*)
      my_pv_system_1_sm_master_DW.StateSpace_PWORK.IDX_OUT_SW;
    real_T *DxCol = (real_T*)my_pv_system_1_sm_master_DW.StateSpace_PWORK.DX_COL;
    real_T *tmp2 = (real_T*)my_pv_system_1_sm_master_DW.StateSpace_PWORK.TMP2;
    real_T *BDcol = (real_T*)my_pv_system_1_sm_master_DW.StateSpace_PWORK.BD_COL;
    real_T *tmp1 = (real_T*)my_pv_system_1_sm_master_DW.StateSpace_PWORK.TMP1;
    int_T newState;
    int_T swChanged = 0;
    int loopsToDo = 20;
    real_T temp;

    /* keep an initial copy of switch_status*/
    memcpy(switch_status_init, switch_status, 4 * sizeof(int_T));
    do {
      if (loopsToDo == 1) {            /* Need to reset some variables: */
        swChanged = 0;

        /* return to the original switch status*/
        {
          int_T i1;
          for (i1=0; i1 < 4; i1++) {
            swChanged = ((SwitchChange[i1] = switch_status_init[i1] -
                          switch_status[i1]) != 0) ? 1 : swChanged;
            switch_status[i1] = switch_status_init[i1];
          }
        }
      } else {
        /*
         * Compute outputs:
         * ---------------
         */

        /*
         * Chopper parameter will force zero current (y[i])
         * for an open switch.
         */
        real_T *Cs = (real_T*)my_pv_system_1_sm_master_DW.StateSpace_PWORK.CS;
        real_T *Ds = (real_T*)my_pv_system_1_sm_master_DW.StateSpace_PWORK.DS;

        {
          int_T i1;
          real_T *y0 = &my_pv_system_1_sm_master_B.StateSpace_o1[0];
          for (i1=0; i1 < 12; i1++) {
            accum = 0.0;

            {
              int_T i2;
              real_T *xd = &my_pv_system_1_sm_master_DW.StateSpace_DSTATE[0];
              for (i2=0; i2 < 10; i2++) {
                accum += *(Cs++) * xd[i2];
              }
            }

            accum += *(Ds++) * my_pv_system_1_sm_master_B.Product[0];
            accum += *(Ds++) * my_pv_system_1_sm_master_B.Product[1];
            accum += *(Ds++) * my_pv_system_1_sm_master_B.Product[2];
            accum += *(Ds++) * my_pv_system_1_sm_master_B.Product[3];
            accum += *(Ds++) * my_pv_system_1_sm_master_B.AC;
            accum += *(Ds++) * my_pv_system_1_sm_master_B.I_Diode_k;
            accum += *(Ds++) * my_pv_system_1_sm_master_B.Delay1;
            y0[i1] = accum * Chopper[i1];
          }
        }

        swChanged = 0;

        /* MOSFETs or IGBT/Diode pairs */
        newState = ((my_pv_system_1_sm_master_B.StateSpace_o1[0] > 0.0) &&
                    (gState[0] > 0)) ||
          (my_pv_system_1_sm_master_B.StateSpace_o1[0] < 0.0) ? 1 :
          (((my_pv_system_1_sm_master_B.StateSpace_o1[0] > 0.0) && gState[0] ==
            0) ? 0 : switch_status[0]);
        swChanged = ((SwitchChange[0] = newState - switch_status[0]) != 0) ? 1 :
          swChanged;
        switch_status[0] = newState;   /* Keep new state */

        /* MOSFETs or IGBT/Diode pairs */
        newState = ((my_pv_system_1_sm_master_B.StateSpace_o1[1] > 0.0) &&
                    (gState[1] > 0)) ||
          (my_pv_system_1_sm_master_B.StateSpace_o1[1] < 0.0) ? 1 :
          (((my_pv_system_1_sm_master_B.StateSpace_o1[1] > 0.0) && gState[1] ==
            0) ? 0 : switch_status[1]);
        swChanged = ((SwitchChange[1] = newState - switch_status[1]) != 0) ? 1 :
          swChanged;
        switch_status[1] = newState;   /* Keep new state */

        /* MOSFETs or IGBT/Diode pairs */
        newState = ((my_pv_system_1_sm_master_B.StateSpace_o1[2] > 0.0) &&
                    (gState[2] > 0)) ||
          (my_pv_system_1_sm_master_B.StateSpace_o1[2] < 0.0) ? 1 :
          (((my_pv_system_1_sm_master_B.StateSpace_o1[2] > 0.0) && gState[2] ==
            0) ? 0 : switch_status[2]);
        swChanged = ((SwitchChange[2] = newState - switch_status[2]) != 0) ? 1 :
          swChanged;
        switch_status[2] = newState;   /* Keep new state */

        /* MOSFETs or IGBT/Diode pairs */
        newState = ((my_pv_system_1_sm_master_B.StateSpace_o1[3] > 0.0) &&
                    (gState[3] > 0)) ||
          (my_pv_system_1_sm_master_B.StateSpace_o1[3] < 0.0) ? 1 :
          (((my_pv_system_1_sm_master_B.StateSpace_o1[3] > 0.0) && gState[3] ==
            0) ? 0 : switch_status[3]);
        swChanged = ((SwitchChange[3] = newState - switch_status[3]) != 0) ? 1 :
          swChanged;
        switch_status[3] = newState;   /* Keep new state */
      }

      /*
       * Compute new As, Bs, Cs and Ds matrixes:
       * --------------------------------------
       */
      if (swChanged) {
        real_T *As = (real_T*)my_pv_system_1_sm_master_DW.StateSpace_PWORK.AS;
        real_T *Cs = (real_T*)my_pv_system_1_sm_master_DW.StateSpace_PWORK.CS;
        real_T *Bs = (real_T*)my_pv_system_1_sm_master_DW.StateSpace_PWORK.BS;
        real_T *Ds = (real_T*)my_pv_system_1_sm_master_DW.StateSpace_PWORK.DS;
        real_T a1;

        {
          int_T i1;
          for (i1=0; i1 < 4; i1++) {
            if (SwitchChange[i1] != 0) {
              if (idxOutSw[i1] > -1 ) {/* A positive index points to a switch actual measure output */
                Chopper[idxOutSw[i1]] = switch_status[i1];
              }

              a1 = 1000.0*SwitchChange[i1];
              temp = 1/(1-Ds[i1*8]*a1);

              {
                int_T i2;
                for (i2=0; i2 < 12; i2++) {
                  DxCol[i2]= Ds[i2 * 7 + i1]*temp*a1;
                }
              }

              DxCol[i1] = temp;

              {
                int_T i2;
                for (i2=0; i2 < 10; i2++) {
                  BDcol[i2]= Bs[i2 * 7 + i1]*a1;
                }
              }

              /* Copy row nSw of Cs into tmp1 and zero it out in Cs */
              memcpy(tmp1, &Cs[i1 * 10], 10 * sizeof(real_T));
              memset(&Cs[i1 * 10], '\0', 10 * sizeof(real_T));

              /* Copy row nSw of Ds into tmp2 and zero it out in Ds */
              memcpy(tmp2, &Ds[i1 * 7], 7 * sizeof(real_T));
              memset(&Ds[i1 * 7], '\0', 7 * sizeof(real_T));

              /* Cs = Cs + DxCol * tmp1, Ds = Ds + DxCol * tmp2 *******************/
              {
                int_T i2;
                for (i2=0; i2 < 12; i2++) {
                  a1 = DxCol[i2];

                  {
                    int_T i3;
                    for (i3=0; i3 < 10; i3++) {
                      Cs[i2 * 10 + i3] += a1 * tmp1[i3];
                    }
                  }

                  {
                    int_T i3;
                    for (i3=0; i3 < 7; i3++) {
                      Ds[i2 * 7 + i3] += a1 * tmp2[i3];
                    }
                  }
                }
              }

              /* As = As + BdCol*Cs(nSw,:), Bs = Bs + BdCol*Ds(nSw,:) *************/
              {
                int_T i2;
                for (i2=0; i2 < 10; i2++) {
                  a1 = BDcol[i2];

                  {
                    int_T i3;
                    for (i3=0; i3 < 10; i3++) {
                      As[i2 * 10 + i3] += a1 * Cs[i1 * 10 + i3];
                    }
                  }

                  {
                    int_T i3;
                    for (i3=0; i3 < 7; i3++) {
                      Bs[i2 * 7 + i3] += a1 * Ds[i1 * 7 + i3];
                    }
                  }
                }
              }
            }
          }
        }
      }                                /* if (swChanged) */
    } while (swChanged > 0 && --loopsToDo > 0);

    if (loopsToDo == 0) {
      real_T *Cs = (real_T*)my_pv_system_1_sm_master_DW.StateSpace_PWORK.CS;
      real_T *Ds = (real_T*)my_pv_system_1_sm_master_DW.StateSpace_PWORK.DS;

      {
        int_T i1;
        real_T *y0 = &my_pv_system_1_sm_master_B.StateSpace_o1[0];
        for (i1=0; i1 < 12; i1++) {
          accum = 0.0;

          {
            int_T i2;
            real_T *xd = &my_pv_system_1_sm_master_DW.StateSpace_DSTATE[0];
            for (i2=0; i2 < 10; i2++) {
              accum += *(Cs++) * xd[i2];
            }
          }

          accum += *(Ds++) * my_pv_system_1_sm_master_B.Product[0];
          accum += *(Ds++) * my_pv_system_1_sm_master_B.Product[1];
          accum += *(Ds++) * my_pv_system_1_sm_master_B.Product[2];
          accum += *(Ds++) * my_pv_system_1_sm_master_B.Product[3];
          accum += *(Ds++) * my_pv_system_1_sm_master_B.AC;
          accum += *(Ds++) * my_pv_system_1_sm_master_B.I_Diode_k;
          accum += *(Ds++) * my_pv_system_1_sm_master_B.Delay1;
          y0[i1] = accum * Chopper[i1];
        }
      }
    }

    /* Output new switches states */
    my_pv_system_1_sm_master_B.StateSpace_o2[0] = (real_T)switch_status[0];
    my_pv_system_1_sm_master_B.StateSpace_o2[1] = (real_T)switch_status[1];
    my_pv_system_1_sm_master_B.StateSpace_o2[2] = (real_T)switch_status[2];
    my_pv_system_1_sm_master_B.StateSpace_o2[3] = (real_T)switch_status[3];
  }

  /* Gain: '<S5>/do not delete this gain' */
  my_pv_system_1_sm_master_B.donotdeletethisgain =
    my_pv_system_1_sm_master_P.donotdeletethisgain_Gain *
    my_pv_system_1_sm_master_B.StateSpace_o1[9];

  /* Gain: '<S10>/do not delete this gain' */
  my_pv_system_1_sm_master_B.donotdeletethisgain_e =
    my_pv_system_1_sm_master_P.donotdeletethisgain_Gain_g *
    my_pv_system_1_sm_master_B.StateSpace_o1[7];

  /* Gain: '<S11>/do not delete this gain' */
  my_pv_system_1_sm_master_B.donotdeletethisgain_g =
    my_pv_system_1_sm_master_P.donotdeletethisgain_Gain_h *
    my_pv_system_1_sm_master_B.StateSpace_o1[8];

  /* Level2 S-Function Block: '<S25>/S-Function' (RECV_Param) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[5];
    sfcnOutputs(rts, 1);
  }

  /* Gain: '<S29>/do not delete this gain' */
  my_pv_system_1_sm_master_B.donotdeletethisgain_a =
    my_pv_system_1_sm_master_P.donotdeletethisgain_Gain_j *
    my_pv_system_1_sm_master_B.StateSpace_o1[11];

  /* Delay: '<S9>/Delay' */
  my_pv_system_1_sm_master_B.Delay = my_pv_system_1_sm_master_DW.Delay_DSTATE;

  /* MATLAB Function: '<S9>/MATLAB Function' incorporates:
   *  Constant: '<S9>/Constant1'
   *  Constant: '<S9>/Constant2'
   *  Constant: '<S9>/Constant3'
   *  Constant: '<S9>/Constant4'
   */
  /* MATLAB Function 'SM_master/Subsystem/MATLAB Function': '<S30>:1' */
  /* '<S30>:1:3' */
  /* '<S30>:1:4' */
  /* '<S30>:1:6' */
  /* '<S30>:1:7' */
  /* '<S30>:1:8' */
  /*  bandgap energy(1.12 eV) */
  /*  ideal factor */
  /* '<S30>:1:13' */
  /* '<S30>:1:14' */
  /* '<S30>:1:17' */
  my_pv_system_1_sm_master_B.Iph = (((my_pv_system_1_sm_master_B.SFunction_d[1]
    + 273.0) - 298.0) * (my_pv_system_1_sm_master_P.Constant3_Value / 100.0) +
    1.0) * (my_pv_system_1_sm_master_B.SFunction_d[0] / 1000.0 *
            my_pv_system_1_sm_master_P.Constant1_Value);

  /* '<S30>:1:18' */
  /* '<S30>:1:19' */
  /* '<S30>:1:21' */
  my_pv_system_1_sm_master_B.Io = my_pv_system_1_sm_master_P.Constant1_Value /
    (exp(my_pv_system_1_sm_master_P.Constant2_Value /
         (my_pv_system_1_sm_master_P.Constant4_Value * 1.381E-23 * 1.3 * 298.0 /
          1.6E-19)) - 1.0) * rt_powd_snf
    ((my_pv_system_1_sm_master_B.SFunction_d[1] + 273.0) / 298.0, 3.0) * exp
    ((0.0033557046979865771 - 1.0 / (my_pv_system_1_sm_master_B.SFunction_d[1] +
       273.0)) * 9981.61867097421) * (exp(my_pv_system_1_sm_master_B.Delay /
    (my_pv_system_1_sm_master_P.Constant4_Value * 1.381E-23 * 1.3 *
     (my_pv_system_1_sm_master_B.SFunction_d[1] + 273.0) / 1.6E-19)) - 1.0);

  /* Gain: '<S31>/do not delete this gain' */
  my_pv_system_1_sm_master_B.donotdeletethisgain_d =
    my_pv_system_1_sm_master_P.donotdeletethisgain_Gain_n *
    my_pv_system_1_sm_master_B.StateSpace_o1[5];

  /* Gain: '<S32>/do not delete this gain' */
  my_pv_system_1_sm_master_B.donotdeletethisgain_db =
    my_pv_system_1_sm_master_P.donotdeletethisgain_Gain_jy *
    my_pv_system_1_sm_master_B.StateSpace_o1[6];

  /* Gain: '<S21>/1//Ron' */
  my_pv_system_1_sm_master_B.Ron[0] = my_pv_system_1_sm_master_P.Ron_Gain *
    my_pv_system_1_sm_master_B.StateSpace_o1[0];
  my_pv_system_1_sm_master_B.Ron[1] = my_pv_system_1_sm_master_P.Ron_Gain *
    my_pv_system_1_sm_master_B.StateSpace_o1[1];
  my_pv_system_1_sm_master_B.Ron[2] = my_pv_system_1_sm_master_P.Ron_Gain *
    my_pv_system_1_sm_master_B.StateSpace_o1[2];
  my_pv_system_1_sm_master_B.Ron[3] = my_pv_system_1_sm_master_P.Ron_Gain *
    my_pv_system_1_sm_master_B.StateSpace_o1[3];

  /* DataTypeConversion: '<S21>/Data Type Conversion' */
  my_pv_system_1_sm_master_B.DataTypeConversion[0] =
    my_pv_system_1_sm_master_B.SFunction_p[0];
  my_pv_system_1_sm_master_B.DataTypeConversion[1] =
    my_pv_system_1_sm_master_B.SFunction_p[1];
  my_pv_system_1_sm_master_B.DataTypeConversion[2] =
    my_pv_system_1_sm_master_B.SFunction_p[2];
  my_pv_system_1_sm_master_B.DataTypeConversion[3] =
    my_pv_system_1_sm_master_B.SFunction_p[3];

  /* Switch: '<S21>/Switch' incorporates:
   *  Constant: '<S21>/0 4'
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

  /* End of Switch: '<S21>/Switch' */

  /* Saturate: '<S21>/Saturation' */
  y_0 = my_pv_system_1_sm_master_B.Switch[0];
  u1 = my_pv_system_1_sm_master_P.Saturation_LowerSat;
  u2 = my_pv_system_1_sm_master_P.Saturation_UpperSat;
  if (y_0 > u2) {
    y_0 = u2;
  } else {
    if (y_0 < u1) {
      y_0 = u1;
    }
  }

  my_pv_system_1_sm_master_B.Saturation[0] = y_0;
  y_0 = my_pv_system_1_sm_master_B.Switch[1];
  u1 = my_pv_system_1_sm_master_P.Saturation_LowerSat;
  u2 = my_pv_system_1_sm_master_P.Saturation_UpperSat;
  if (y_0 > u2) {
    y_0 = u2;
  } else {
    if (y_0 < u1) {
      y_0 = u1;
    }
  }

  my_pv_system_1_sm_master_B.Saturation[1] = y_0;
  y_0 = my_pv_system_1_sm_master_B.Switch[2];
  u1 = my_pv_system_1_sm_master_P.Saturation_LowerSat;
  u2 = my_pv_system_1_sm_master_P.Saturation_UpperSat;
  if (y_0 > u2) {
    y_0 = u2;
  } else {
    if (y_0 < u1) {
      y_0 = u1;
    }
  }

  my_pv_system_1_sm_master_B.Saturation[2] = y_0;
  y_0 = my_pv_system_1_sm_master_B.Switch[3];
  u1 = my_pv_system_1_sm_master_P.Saturation_LowerSat;
  u2 = my_pv_system_1_sm_master_P.Saturation_UpperSat;
  if (y_0 > u2) {
    y_0 = u2;
  } else {
    if (y_0 < u1) {
      y_0 = u1;
    }
  }

  my_pv_system_1_sm_master_B.Saturation[3] = y_0;

  /* End of Saturate: '<S21>/Saturation' */

  /* Switch: '<S23>/Switch' incorporates:
   *  Constant: '<S23>/Vf Devices & Clamping Diodes'
   *  Constant: '<S23>/Vf Diodes'
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

  /* End of Switch: '<S23>/Switch' */

  /* Gain: '<S17>/do not delete this gain' */
  my_pv_system_1_sm_master_B.donotdeletethisgain_i =
    my_pv_system_1_sm_master_P.donotdeletethisgain_Gain_c *
    my_pv_system_1_sm_master_B.StateSpace_o1[10];

  /* Gain: '<S18>/do not delete this gain' */
  my_pv_system_1_sm_master_B.donotdeletethisgain_de =
    my_pv_system_1_sm_master_P.donotdeletethisgain_Gain_k *
    my_pv_system_1_sm_master_B.StateSpace_o1[4];
}

/* Model update function */
static void my_pv_system_1_sm_master_update(void)
{
  real_T tmp;

  /* Update for Memory: '<S1>/S-Function' */
  my_pv_system_1_sm_master_DW.SFunction_PreviousInput =
    my_pv_system_1_sm_master_B.Sum;

  /* Update for Memory: '<S2>/Memory' */
  my_pv_system_1_sm_master_DW.Memory_1_PreviousInput =
    my_pv_system_1_sm_master_B.SFunction_d[0];

  /* Update for Memory: '<S2>/Memory' */
  my_pv_system_1_sm_master_DW.Memory_2_PreviousInput =
    my_pv_system_1_sm_master_B.donotdeletethisgain_a;

  /* Update for Memory: '<S2>/Memory' */
  my_pv_system_1_sm_master_DW.Memory_3_PreviousInput =
    my_pv_system_1_sm_master_B.donotdeletethisgain_d;

  /* Update for Memory: '<S2>/Memory' */
  my_pv_system_1_sm_master_DW.Memory_4_PreviousInput =
    my_pv_system_1_sm_master_B.I_Diode_k;

  /* Update for Memory: '<S2>/Memory1' */
  my_pv_system_1_sm_master_DW.Memory1_1_PreviousInput =
    my_pv_system_1_sm_master_B.donotdeletethisgain_g;

  /* Update for Memory: '<S2>/Memory1' */
  my_pv_system_1_sm_master_DW.Memory1_2_PreviousInput =
    my_pv_system_1_sm_master_B.donotdeletethisgain_e;

  /* Update for Memory: '<S2>/Memory1' */
  my_pv_system_1_sm_master_DW.Memory1_3_PreviousInput =
    my_pv_system_1_sm_master_B.donotdeletethisgain_i;

  /* Update for Memory: '<S2>/Memory1' */
  my_pv_system_1_sm_master_DW.Memory1_4_PreviousInput =
    my_pv_system_1_sm_master_B.donotdeletethisgain_de;

  /* Update for Memory: '<S2>/Memory1' */
  my_pv_system_1_sm_master_DW.Memory1_5_PreviousInput =
    my_pv_system_1_sm_master_B.donotdeletethisgain;

  /* Update for UnitDelay: '<S21>/Unit Delay' */
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[0] =
    my_pv_system_1_sm_master_B.Saturation[0];
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[1] =
    my_pv_system_1_sm_master_B.Saturation[1];
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[2] =
    my_pv_system_1_sm_master_B.Saturation[2];
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[3] =
    my_pv_system_1_sm_master_B.Saturation[3];

  /* Update for Enabled SubSystem: '<S21>/Tail' incorporates:
   *  Update for EnablePort: '<S22>/Enable'
   */
  if (my_pv_system_1_sm_master_DW.Tail_MODE) {
    /* Update for DiscreteIntegrator: '<S22>/Discrete-Time Integrator' incorporates:
     *  Constant: '<S22>/1'
     */
    tmp = my_pv_system_1_sm_master_P.DiscreteTimeIntegrator_gainval *
      my_pv_system_1_sm_master_P._Value;
    my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[0] += tmp;
    my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[1] += tmp;
    my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[2] += tmp;
    my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_DSTATE[3] += tmp;
    if (my_pv_system_1_sm_master_B.SFunction_p[0] > 0.0) {
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[0] = 1;
    } else if (my_pv_system_1_sm_master_B.SFunction_p[0] < 0.0) {
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[0] = -1;
    } else if (my_pv_system_1_sm_master_B.SFunction_p[0] == 0.0) {
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[0] = 0;
    } else {
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[0] = 2;
    }

    if (my_pv_system_1_sm_master_B.SFunction_p[1] > 0.0) {
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[1] = 1;
    } else if (my_pv_system_1_sm_master_B.SFunction_p[1] < 0.0) {
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[1] = -1;
    } else if (my_pv_system_1_sm_master_B.SFunction_p[1] == 0.0) {
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[1] = 0;
    } else {
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[1] = 2;
    }

    if (my_pv_system_1_sm_master_B.SFunction_p[2] > 0.0) {
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[2] = 1;
    } else if (my_pv_system_1_sm_master_B.SFunction_p[2] < 0.0) {
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[2] = -1;
    } else if (my_pv_system_1_sm_master_B.SFunction_p[2] == 0.0) {
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[2] = 0;
    } else {
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[2] = 2;
    }

    if (my_pv_system_1_sm_master_B.SFunction_p[3] > 0.0) {
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[3] = 1;
    } else if (my_pv_system_1_sm_master_B.SFunction_p[3] < 0.0) {
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[3] = -1;
    } else if (my_pv_system_1_sm_master_B.SFunction_p[3] == 0.0) {
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[3] = 0;
    } else {
      my_pv_system_1_sm_master_DW.DiscreteTimeIntegrator_PrevRese[3] = 2;
    }

    /* End of Update for DiscreteIntegrator: '<S22>/Discrete-Time Integrator' */

    /* Update for UnitDelay: '<S22>/Unit Delay' */
    my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[0] =
      my_pv_system_1_sm_master_B.Switch_j[0];
    my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[1] =
      my_pv_system_1_sm_master_B.Switch_j[1];
    my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[2] =
      my_pv_system_1_sm_master_B.Switch_j[2];
    my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[3] =
      my_pv_system_1_sm_master_B.Switch_j[3];
  }

  /* End of Update for SubSystem: '<S21>/Tail' */

  /* Update for Delay: '<S9>/Delay2' */
  my_pv_system_1_sm_master_DW.Delay2_DSTATE = my_pv_system_1_sm_master_B.Io;

  /* Update for Delay: '<S9>/Delay1' */
  my_pv_system_1_sm_master_DW.Delay1_DSTATE = my_pv_system_1_sm_master_B.Iph;

  /* S-Function block: <S50>/State-Space */
  {
    const real_T *As = (real_T*)my_pv_system_1_sm_master_DW.StateSpace_PWORK.AS;
    const real_T *Bs = (real_T*)my_pv_system_1_sm_master_DW.StateSpace_PWORK.BS;
    real_T *xtmp = (real_T*)my_pv_system_1_sm_master_DW.StateSpace_PWORK.XTMP;
    real_T accum;

    /* Calculate new states... */
    {
      int_T i1;
      real_T *xd = &my_pv_system_1_sm_master_DW.StateSpace_DSTATE[0];
      for (i1=0; i1 < 10; i1++) {
        accum = 0.0;

        {
          int_T i2;
          real_T *xd = &my_pv_system_1_sm_master_DW.StateSpace_DSTATE[0];
          for (i2=0; i2 < 10; i2++) {
            accum += *(As++) * xd[i2];
          }
        }

        accum += *(Bs++) * my_pv_system_1_sm_master_B.Product[0];
        accum += *(Bs++) * my_pv_system_1_sm_master_B.Product[1];
        accum += *(Bs++) * my_pv_system_1_sm_master_B.Product[2];
        accum += *(Bs++) * my_pv_system_1_sm_master_B.Product[3];
        accum += *(Bs++) * my_pv_system_1_sm_master_B.AC;
        accum += *(Bs++) * my_pv_system_1_sm_master_B.I_Diode_k;
        accum += *(Bs++) * my_pv_system_1_sm_master_B.Delay1;
        xtmp[i1] = accum;
      }
    }

    {
      int_T i1;
      real_T *xd = &my_pv_system_1_sm_master_DW.StateSpace_DSTATE[0];
      for (i1=0; i1 < 10; i1++) {
        xd[i1] = xtmp[i1];
      }
    }

    {
      int_T *gState = (int_T*)
        my_pv_system_1_sm_master_DW.StateSpace_PWORK.G_STATE;

      /* Store switch gates values for next step */
      *(gState++) = (int_T) my_pv_system_1_sm_master_B.DataTypeConversion[0];
      *(gState++) = (int_T) my_pv_system_1_sm_master_B.DataTypeConversion[1];
      *(gState++) = (int_T) my_pv_system_1_sm_master_B.DataTypeConversion[2];
      *(gState++) = (int_T) my_pv_system_1_sm_master_B.DataTypeConversion[3];
    }
  }

  /* Update for Delay: '<S9>/Delay' */
  my_pv_system_1_sm_master_DW.Delay_DSTATE =
    my_pv_system_1_sm_master_B.donotdeletethisgain_db;

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

  my_pv_system_1_sm_master_M->Timing.t[0] =
    my_pv_system_1_sm_master_M->Timing.clockTick0 *
    my_pv_system_1_sm_master_M->Timing.stepSize0 +
    my_pv_system_1_sm_master_M->Timing.clockTickH0 *
    my_pv_system_1_sm_master_M->Timing.stepSize0 * 4294967296.0;

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

/* Model initialize function */
static void my_pv_system_1_sm_master_initialize(void)
{
  /* InitializeConditions for Enabled SubSystem: '<S21>/Tail' */
  /* InitializeConditions for DiscreteIntegrator: '<S22>/Discrete-Time Integrator' */
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

  /* InitializeConditions for UnitDelay: '<S22>/Unit Delay' */
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[0] =
    my_pv_system_1_sm_master_P.UnitDelay_InitialCondition;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[1] =
    my_pv_system_1_sm_master_P.UnitDelay_InitialCondition;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[2] =
    my_pv_system_1_sm_master_P.UnitDelay_InitialCondition;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE_e[3] =
    my_pv_system_1_sm_master_P.UnitDelay_InitialCondition;

  /* End of InitializeConditions for SubSystem: '<S21>/Tail' */

  /* Start for Enabled SubSystem: '<S21>/Tail' */
  /* VirtualOutportStart for Outport: '<S22>/itail' */
  my_pv_system_1_sm_master_B.Product[0] = my_pv_system_1_sm_master_P.itail_Y0;
  my_pv_system_1_sm_master_B.Product[1] = my_pv_system_1_sm_master_P.itail_Y0;
  my_pv_system_1_sm_master_B.Product[2] = my_pv_system_1_sm_master_P.itail_Y0;
  my_pv_system_1_sm_master_B.Product[3] = my_pv_system_1_sm_master_P.itail_Y0;

  /* End of Start for SubSystem: '<S21>/Tail' */

  /* S-Function block: <S50>/State-Space */
  {
    my_pv_system_1_sm_master_DW.StateSpace_PWORK.AS = (real_T*)calloc(10 * 10,
      sizeof(real_T));
    my_pv_system_1_sm_master_DW.StateSpace_PWORK.BS = (real_T*)calloc(10 * 7,
      sizeof(real_T));
    my_pv_system_1_sm_master_DW.StateSpace_PWORK.CS = (real_T*)calloc(12 * 10,
      sizeof(real_T));
    my_pv_system_1_sm_master_DW.StateSpace_PWORK.DS = (real_T*)calloc(12 * 7,
      sizeof(real_T));
    my_pv_system_1_sm_master_DW.StateSpace_PWORK.DX_COL = (real_T*)calloc(12,
      sizeof(real_T));
    my_pv_system_1_sm_master_DW.StateSpace_PWORK.TMP2 = (real_T*)calloc(7,
      sizeof(real_T));
    my_pv_system_1_sm_master_DW.StateSpace_PWORK.BD_COL = (real_T*)calloc(10,
      sizeof(real_T));
    my_pv_system_1_sm_master_DW.StateSpace_PWORK.TMP1 = (real_T*)calloc(10,
      sizeof(real_T));
    my_pv_system_1_sm_master_DW.StateSpace_PWORK.XTMP = (real_T*)calloc(10,
      sizeof(real_T));
    my_pv_system_1_sm_master_DW.StateSpace_PWORK.CHOPPER = (int_T*)calloc(12,
      sizeof(int_T));
    my_pv_system_1_sm_master_DW.StateSpace_PWORK.SWITCH_STATUS = (int_T*)calloc
      (4, sizeof(int_T));
    my_pv_system_1_sm_master_DW.StateSpace_PWORK.SW_CHG = (int_T*)calloc(4,
      sizeof(int_T));
    my_pv_system_1_sm_master_DW.StateSpace_PWORK.G_STATE = (int_T*)calloc(4,
      sizeof(int_T));
    my_pv_system_1_sm_master_DW.StateSpace_PWORK.Y_SWITCH = (real_T*)calloc(4,
      sizeof(real_T));
    my_pv_system_1_sm_master_DW.StateSpace_PWORK.SWITCH_TYPES = (int_T*)calloc(4,
      sizeof(int_T));
    my_pv_system_1_sm_master_DW.StateSpace_PWORK.IDX_OUT_SW = (int_T*)calloc(4,
      sizeof(int_T));
    my_pv_system_1_sm_master_DW.StateSpace_PWORK.SWITCH_STATUS_INIT = (int_T*)
      calloc(4, sizeof(int_T));
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
  my_pv_system_1_sm_master_DW.SFunction_PreviousInput =
    my_pv_system_1_sm_master_P.SFunction_X0;

  /* InitializeConditions for Memory: '<S2>/Memory' */
  my_pv_system_1_sm_master_DW.Memory_1_PreviousInput =
    my_pv_system_1_sm_master_P.Memory_1_X0;

  /* InitializeConditions for Memory: '<S2>/Memory' */
  my_pv_system_1_sm_master_DW.Memory_2_PreviousInput =
    my_pv_system_1_sm_master_P.Memory_2_X0;

  /* InitializeConditions for Memory: '<S2>/Memory' */
  my_pv_system_1_sm_master_DW.Memory_3_PreviousInput =
    my_pv_system_1_sm_master_P.Memory_3_X0;

  /* InitializeConditions for Memory: '<S2>/Memory' */
  my_pv_system_1_sm_master_DW.Memory_4_PreviousInput =
    my_pv_system_1_sm_master_P.Memory_4_X0;

  /* InitializeConditions for Memory: '<S2>/Memory1' */
  my_pv_system_1_sm_master_DW.Memory1_1_PreviousInput =
    my_pv_system_1_sm_master_P.Memory1_1_X0;

  /* InitializeConditions for Memory: '<S2>/Memory1' */
  my_pv_system_1_sm_master_DW.Memory1_2_PreviousInput =
    my_pv_system_1_sm_master_P.Memory1_2_X0;

  /* InitializeConditions for Memory: '<S2>/Memory1' */
  my_pv_system_1_sm_master_DW.Memory1_3_PreviousInput =
    my_pv_system_1_sm_master_P.Memory1_3_X0;

  /* InitializeConditions for Memory: '<S2>/Memory1' */
  my_pv_system_1_sm_master_DW.Memory1_4_PreviousInput =
    my_pv_system_1_sm_master_P.Memory1_4_X0;

  /* InitializeConditions for Memory: '<S2>/Memory1' */
  my_pv_system_1_sm_master_DW.Memory1_5_PreviousInput =
    my_pv_system_1_sm_master_P.Memory1_5_X0;

  /* InitializeConditions for Atomic SubSystem: '<S13>/Subsystem5' */

  /* Level2 S-Function Block: '<S46>/S-Function' (send_rt) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[0];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* End of InitializeConditions for SubSystem: '<S13>/Subsystem5' */

  /* Level2 S-Function Block: '<S44>/S-Function' (OP_SEND) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[1];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Level2 S-Function Block: '<S2>/OpMonitor' (opmonitor) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[2];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Level2 S-Function Block: '<S45>/S-Function' (OP_SEND) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[3];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* InitializeConditions for UnitDelay: '<S21>/Unit Delay' */
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[0] =
    my_pv_system_1_sm_master_P.UnitDelay_InitialCondition_p;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[1] =
    my_pv_system_1_sm_master_P.UnitDelay_InitialCondition_p;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[2] =
    my_pv_system_1_sm_master_P.UnitDelay_InitialCondition_p;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[3] =
    my_pv_system_1_sm_master_P.UnitDelay_InitialCondition_p;

  /* Level2 S-Function Block: '<S47>/S-Function' (recv_rt) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[4];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* InitializeConditions for Delay: '<S9>/Delay2' */
  my_pv_system_1_sm_master_DW.Delay2_DSTATE =
    my_pv_system_1_sm_master_P.Delay2_InitialCondition;

  /* InitializeConditions for Delay: '<S9>/Delay1' */
  my_pv_system_1_sm_master_DW.Delay1_DSTATE =
    my_pv_system_1_sm_master_P.Delay1_InitialCondition;

  /* InitializeConditions for S-Function (sfun_spssw_discc): '<S50>/State-Space' */
  {
    real_T *As = (real_T*)my_pv_system_1_sm_master_DW.StateSpace_PWORK.AS;
    real_T *Bs = (real_T*)my_pv_system_1_sm_master_DW.StateSpace_PWORK.BS;
    real_T *Cs = (real_T*)my_pv_system_1_sm_master_DW.StateSpace_PWORK.CS;
    real_T *Ds = (real_T*)my_pv_system_1_sm_master_DW.StateSpace_PWORK.DS;
    real_T *X0 = (real_T*)&my_pv_system_1_sm_master_DW.StateSpace_DSTATE[0];
    int_T *Chopper = (int_T*)
      my_pv_system_1_sm_master_DW.StateSpace_PWORK.CHOPPER;
    X0[0] = 8.500016821605308E+6;
    X0[1] = -76192.288557594118;
    X0[2] = -24335.850149922797;
    X0[3] = 27520.735648827354;
    X0[4] = -8.41080432615598;
    X0[5] = 187428.76021043005;
    X0[6] = 469206.94944928697;
    X0[7] = -469206.94944928709;
    X0[8] = -10716.825349448978;
    X0[9] = -3523.7302283658632;

    /* Copy and transpose A and B */
    As[0] = 0.99999604198305492;
    As[1] = 0.0;
    As[2] = 0.0;
    As[3] = 0.0;
    As[4] = 0.0;
    As[5] = 0.0;
    As[6] = 0.0;
    As[7] = 0.0;
    As[8] = 0.0;
    As[9] = 0.0;
    As[10] = 1.9746258200928638E-6;
    As[11] = -0.026898014275293838;
    As[12] = -6.2770915771320074E-6;
    As[13] = -0.972406011326852;
    As[14] = 3.94925164016541E-6;
    As[15] = -0.011140069399679358;
    As[16] = -0.0091796697954719675;
    As[17] = 0.0091796546506022689;
    As[18] = -1.0944611153800048;
    As[19] = 1.0978994663797219;
    As[20] = 1.9746258200921968E-6;
    As[21] = 0.00016710115735295506;
    As[22] = -0.99965078209984171;
    As[23] = 0.0001733782489305182;
    As[24] = 3.9492516401844172E-6;
    As[25] = -1.9637013775810864E-6;
    As[26] = 2.2045320775900863E-8;
    As[27] = -3.7190190475927092E-8;
    As[28] = 3.5312281726489978E-6;
    As[29] = -3.5423218531343406E-6;
    As[30] = 1.9746258200911295E-6;
    As[31] = -0.97241228841845151;
    As[32] = 6.27709157713114E-6;
    As[33] = -0.026904291366893391;
    As[34] = 3.9492516402535995E-6;
    As[35] = 0.011140069399679361;
    As[36] = 0.0091796546506022758;
    As[37] = -0.0091796697954719779;
    As[38] = 1.0944611153800055;
    As[39] = -1.0978994663797232;
    As[40] = -0.01233943236710369;
    As[41] = -2.1553040826300554;
    As[42] = 2.7328017233509883E-16;
    As[43] = -2.1553040826318886;
    As[44] = 0.97531717724884748;
    As[45] = -1.1981500649522748E-14;
    As[46] = 4.7327717802710895E-5;
    As[47] = 4.7327717822455572E-5;
    As[48] = -1.1774144643136774E-12;
    As[49] = 1.1811138885264033E-12;
    As[50] = 1.5177839283960908E-18;
    As[51] = 1.0052768611564225;
    As[52] = 0.00035447025025181489;
    As[53] = -1.0056313314066079;
    As[54] = -9.8863144504680173E-17;
    As[55] = 0.988483379163918;
    As[56] = -0.0094882457421109744;
    As[57] = 0.0094882457421109779;
    As[58] = -1.1312525194935341;
    As[59] = 1.13480645409813;
    As[60] = -8.260245034785625E-12;
    As[61] = -0.00087437942724596473;
    As[62] = 5.6422620546020625E-9;
    As[63] = 0.0008743708998453204;
    As[64] = -1.652048987904197E-11;
    As[65] = 1.0013425817849203E-5;
    As[66] = 0.99910520242946776;
    As[67] = 0.00088643118041807677;
    As[68] = -0.10618504026323157;
    As[69] = 0.10651863040564366;
    As[70] = -8.2602450341504376E-12;
    As[71] = 0.00087437654210751887;
    As[72] = -5.6422620546054506E-9;
    As[73] = -0.00087437378498376626;
    As[74] = -1.6520489966154015E-11;
    As[75] = -1.0013425817849198E-5;
    As[76] = 0.00088643118041807655;
    As[77] = 0.99910520242946754;
    As[78] = 0.10618504026323154;
    As[79] = -0.10651863040564363;
    As[80] = -1.806704629483595E-20;
    As[81] = -0.010859294632359674;
    As[82] = 7.00737977369823E-8;
    As[83] = 0.010859224558563311;
    As[84] = 1.0891469547897699E-18;
    As[85] = 0.00012436125238849216;
    As[86] = -0.011060941694086621;
    As[87] = 0.011060941694086619;
    As[88] = -0.33124657050830419;
    As[89] = 1.3291456196471505;
    As[90] = 1.38368853348407E-23;
    As[91] = 6.5360460675609146E-6;
    As[92] = -4.2176364639105434E-11;
    As[93] = -6.5360038911961563E-6;
    As[94] = -1.2197121191924856E-21;
    As[95] = -7.4851166871232261E-8;
    As[96] = 6.6574144003527313E-6;
    As[97] = -6.6574144003527322E-6;
    As[98] = 0.00079748737178829;
    As[99] = 0.99920000724774327;
    Bs[0] = -166.6662257211041;
    Bs[1] = -166.6662257211041;
    Bs[2] = -166.6662257211041;
    Bs[3] = -166.6662257211041;
    Bs[4] = 0.0;
    Bs[5] = -333.05986765390293;
    Bs[6] = 333.05986765390293;
    Bs[7] = 20105.688547766218;
    Bs[8] = -20105.688547769179;
    Bs[9] = 19386.827854080508;
    Bs[10] = -19386.827854083469;
    Bs[11] = -2.4872330727123226;
    Bs[12] = -3.385312863368739E-15;
    Bs[13] = 3.385312863368739E-15;
    Bs[14] = -250.18592174786178;
    Bs[15] = 250.18592174786093;
    Bs[16] = 39742.7023235918;
    Bs[17] = -39742.70232359181;
    Bs[18] = 8.0249424807171557E-6;
    Bs[19] = 3.8489177123235407E-18;
    Bs[20] = -3.8489177123235407E-18;
    Bs[21] = 19386.827854073057;
    Bs[22] = -19386.827854076018;
    Bs[23] = 20105.688547758768;
    Bs[24] = -20105.688547761729;
    Bs[25] = 2.4872330727123226;
    Bs[26] = 1.1878627421968435E-14;
    Bs[27] = -1.1878627421968435E-14;
    Bs[28] = -1.234140304226495E+8;
    Bs[29] = 1.234141970888752E+8;
    Bs[30] = -1.2341403042264949E+8;
    Bs[31] = 1.234141970888752E+8;
    Bs[32] = -2.6760885609852398E-12;
    Bs[33] = 166.52993382695146;
    Bs[34] = -166.52993382695146;
    Bs[35] = -20297.197105579078;
    Bs[36] = 20297.197105579078;
    Bs[37] = 20297.197105597705;
    Bs[38] = -20297.197105597705;
    Bs[39] = -2.5708438980004567;
    Bs[40] = -1.6969252419492379E-14;
    Bs[41] = 1.6969252419492379E-14;
    Bs[42] = -0.40568204578162514;
    Bs[43] = 0.40568204578007461;
    Bs[44] = 0.24047714509210669;
    Bs[45] = -0.24047714509365728;
    Bs[46] = -0.24131231366616357;
    Bs[47] = 3.1759779530756728E-17;
    Bs[48] = -3.1759779530756728E-17;
    Bs[49] = 0.24047714508846871;
    Bs[50] = -0.2404771450900193;
    Bs[51] = -0.40568204578526312;
    Bs[52] = 0.40568204578371259;
    Bs[53] = 0.24131231366616357;
    Bs[54] = 1.7041402926250065E-17;
    Bs[55] = -1.7041402926250065E-17;
    Bs[56] = -4.0124712402876721;
    Bs[57] = 4.0124712402876721;
    Bs[58] = 4.0124712400204343;
    Bs[59] = -4.0124712400204343;
    Bs[60] = 1.5197850558116974;
    Bs[61] = 1.8739371941459857E-16;
    Bs[62] = -1.8739371941459857E-16;
    Bs[63] = 0.0024150460743174208;
    Bs[64] = -0.0024150460743174208;
    Bs[65] = -0.0024150460740591895;
    Bs[66] = 0.0024150460740591895;
    Bs[67] = 0.0018123411954142964;
    Bs[68] = -2.0772791276375983E-19;
    Bs[69] = 2.0772791276375983E-19;
    Cs[0] = 4.9682773275935779E-5;
    Cs[1] = 0.006532451863604365;
    Cs[2] = -0.0043652237519788244;
    Cs[3] = 0.0064539882188939376;
    Cs[4] = 4.9365645503569615E-5;
    Cs[5] = 2.4546267219807305E-5;
    Cs[6] = -4.6487738081850693E-7;
    Cs[7] = 2.7556650959468242E-7;
    Cs[8] = -4.4140352157410235E-5;
    Cs[9] = 4.4279023152193989E-5;
    Cs[10] = 3.1712777364060197E-7;
    Cs[11] = -0.006532451863604365;
    Cs[12] = 0.0043652237519788244;
    Cs[13] = -0.0064539882188939376;
    Cs[14] = -4.9365645503569615E-5;
    Cs[15] = -2.4546267219807305E-5;
    Cs[16] = 4.6487738081850693E-7;
    Cs[17] = -2.7556650959468242E-7;
    Cs[18] = 4.4140352157410235E-5;
    Cs[19] = -4.4279023152193989E-5;
    Cs[20] = 4.9682773275940652E-5;
    Cs[21] = 0.0020887644669121976;
    Cs[22] = 0.0043652237519788235;
    Cs[23] = 0.0021672281116315111;
    Cs[24] = 4.936564550230551E-5;
    Cs[25] = -2.4546267219763581E-5;
    Cs[26] = 2.7556650969875963E-7;
    Cs[27] = -4.6487738094908696E-7;
    Cs[28] = 4.4140352158112459E-5;
    Cs[29] = -4.4279023164179248E-5;
    Cs[30] = 3.1712777363572253E-7;
    Cs[31] = -0.0020887644669121976;
    Cs[32] = -0.0043652237519788235;
    Cs[33] = -0.0021672281116315111;
    Cs[34] = -4.936564550230551E-5;
    Cs[35] = 2.4546267219763581E-5;
    Cs[36] = -2.7556650969875963E-7;
    Cs[37] = 4.6487738094908696E-7;
    Cs[38] = -4.4140352158112459E-5;
    Cs[39] = 4.4279023164179248E-5;
    Cs[40] = -4.4219665156954679E-23;
    Cs[41] = 4.0078534576901142E-5;
    Cs[42] = -2.5862224210854778E-10;
    Cs[43] = -4.0078275954659044E-5;
    Cs[44] = 5.369481714531313E-21;
    Cs[45] = -4.589816302640751E-7;
    Cs[46] = 4.0822755910727191E-5;
    Cs[47] = -4.0822755910727225E-5;
    Cs[48] = 0.0048671637236079061;
    Cs[49] = -0.0048824543694058094;
    Cs[50] = 4.9999901049576376E-5;
    Cs[51] = 0.0;
    Cs[52] = 0.0;
    Cs[53] = 0.0;
    Cs[54] = 0.0;
    Cs[55] = 0.0;
    Cs[56] = 0.0;
    Cs[57] = 0.0;
    Cs[58] = 0.0;
    Cs[59] = 0.0;
    Cs[60] = 4.9959013454072214E-5;
    Cs[61] = 0.0;
    Cs[62] = 0.0;
    Cs[63] = 0.0;
    Cs[64] = 0.0;
    Cs[65] = 0.0;
    Cs[66] = 0.0;
    Cs[67] = 0.0;
    Cs[68] = 0.0;
    Cs[69] = 0.0;
    Cs[70] = 4.9999901049576376E-5;
    Cs[71] = 0.0;
    Cs[72] = 0.0;
    Cs[73] = 0.0;
    Cs[74] = 0.0;
    Cs[75] = 0.0;
    Cs[76] = 0.0;
    Cs[77] = 0.0;
    Cs[78] = 0.0;
    Cs[79] = 0.0;
    Cs[80] = 4.8794301932275633E-18;
    Cs[81] = -0.0044436873966929857;
    Cs[82] = 0.0087304475039576487;
    Cs[83] = -0.0042867601072633993;
    Cs[84] = -1.2641054869533265E-15;
    Cs[85] = -4.9092534439569135E-5;
    Cs[86] = 7.4044389050413888E-7;
    Cs[87] = -7.4044389052687621E-7;
    Cs[88] = 8.828070431482047E-5;
    Cs[89] = -8.85580463160295E-5;
    Cs[90] = 9.87312910046001E-11;
    Cs[91] = 1.7242432661035423E-8;
    Cs[92] = -2.1811098391798234E-24;
    Cs[93] = 1.7242432661052076E-8;
    Cs[94] = 1.9746258201175027E-10;
    Cs[95] = 8.6736173798840359E-23;
    Cs[96] = -3.7862174220512636E-13;
    Cs[97] = -3.7862174272554343E-13;
    Cs[98] = 0.0;
    Cs[99] = -2.2204460492503132E-20;
    Cs[100] = 4.9365645502321738E-11;
    Cs[101] = 2.4327549643119992E-5;
    Cs[102] = -1.5692728942835346E-10;
    Cs[103] = -2.4310150283169517E-5;
    Cs[104] = 9.8731291003627025E-11;
    Cs[105] = -2.78501734991984E-7;
    Cs[106] = -2.2949174488679916E-7;
    Cs[107] = 2.2949136626505678E-7;
    Cs[108] = -2.7361527884500131E-5;
    Cs[109] = 2.7447486659493056E-5;
    Cs[110] = -1.1824058850249033E-8;
    Cs[111] = 0.0;
    Cs[112] = 0.0;
    Cs[113] = 0.0;
    Cs[114] = 0.0;
    Cs[115] = 0.0;
    Cs[116] = 0.0;
    Cs[117] = 0.0;
    Cs[118] = 0.0;
    Cs[119] = 0.0;
    Ds[0] = -3216.2230385176372;
    Ds[1] = 3216.2188718621037;
    Ds[2] = -3127.3261052635435;
    Ds[3] = 3127.3219386079072;
    Ds[4] = -0.0001003117810073686;
    Ds[5] = -0.0041632483455676713;
    Ds[6] = 0.0041632483455676713;
    Ds[7] = 3216.2188718621037;
    Ds[8] = -3216.2230385176372;
    Ds[9] = 3127.3219386079077;
    Ds[10] = -3127.326105263543;
    Ds[11] = 0.0001003117810073686;
    Ds[12] = -0.0041632483457799043;
    Ds[13] = 0.0041632483457799043;
    Ds[14] = -3127.3261051760919;
    Ds[15] = 3127.3219385204552;
    Ds[16] = -3216.2230384302675;
    Ds[17] = 3216.2188717745594;
    Ds[18] = 0.00010031178100896446;
    Ds[19] = -0.0041632483456737384;
    Ds[20] = 0.0041632483456737384;
    Ds[21] = 3127.3219385204557;
    Ds[22] = -3127.326105176091;
    Ds[23] = 3216.2188717746758;
    Ds[24] = -3216.2230384302093;
    Ds[25] = -0.00010031178100896446;
    Ds[26] = -0.0041632483456738355;
    Ds[27] = 0.0041632483456738355;
    Ds[28] = 0.014808877811953426;
    Ds[29] = -0.014808877811953426;
    Ds[30] = -0.014808877812325955;
    Ds[31] = 0.014808877812325955;
    Ds[32] = 0.011060941694086616;
    Ds[33] = 4.5445451458473518E-19;
    Ds[34] = -4.5445451458473518E-19;
    Ds[35] = -0.0041666556430276025;
    Ds[36] = -0.0041666556430276025;
    Ds[37] = -0.0041666556430276025;
    Ds[38] = -0.0041666556430276025;
    Ds[39] = 0.0;
    Ds[40] = -0.0083264966913475739;
    Ds[41] = 0.0083264966913475739;
    Ds[42] = -0.0041632483456737878;
    Ds[43] = -0.0041632483456737878;
    Ds[44] = -0.0041632483456737878;
    Ds[45] = -0.0041632483456737878;
    Ds[46] = 0.0;
    Ds[47] = -3.4634918959680077;
    Ds[48] = 3.4634918959680077;
    Ds[49] = -0.0041666556430276025;
    Ds[50] = -0.0041666556430276025;
    Ds[51] = -0.0041666556430276025;
    Ds[52] = -0.0041666556430276025;
    Ds[53] = 0.0;
    Ds[54] = -0.0083264966913475739;
    Ds[55] = 0.0083264966913475739;
    Ds[56] = 88.896933290292509;
    Ds[57] = -88.896933290292509;
    Ds[58] = -88.896933217882179;
    Ds[59] = 88.896933217882179;
    Ds[60] = 0.00020062356201473719;
    Ds[61] = -1.060684826439291E-13;
    Ds[62] = 1.060684826439291E-13;
    Ds[63] = 0.9873129100460234;
    Ds[64] = -0.9873129100460234;
    Ds[65] = 0.9873129100460234;
    Ds[66] = -0.9873129100460234;
    Ds[67] = 0.0;
    Ds[68] = 2.1233286396499241E-19;
    Ds[69] = -2.1233286396499241E-19;
    Ds[70] = 0.50264221369415552;
    Ds[71] = -0.50264221369422946;
    Ds[72] = 0.48467069635201271;
    Ds[73] = -0.48467069635208676;
    Ds[74] = -6.2180826817808072E-5;
    Ds[75] = -8.4632821584218484E-20;
    Ds[76] = 8.4632821584218484E-20;
    Ds[77] = 9.8533758062903236E-7;
    Ds[78] = 9.8533758062903236E-7;
    Ds[79] = 9.8533758062903236E-7;
    Ds[80] = 9.8533758062903236E-7;
    Ds[81] = 0.0;
    Ds[82] = -0.99918027740794091;
    Ds[83] = 0.99918027740794091;

    {
      int_T i1;
      for (i1=0; i1 < 12; i1++) {
        Chopper[i1] = 1;
      }
    }

    {
      /* Switches work vectors */
      int_T *switch_status = (int_T*)
        my_pv_system_1_sm_master_DW.StateSpace_PWORK.SWITCH_STATUS;
      int_T *gState = (int_T*)
        my_pv_system_1_sm_master_DW.StateSpace_PWORK.G_STATE;
      real_T *yswitch = (real_T*)
        my_pv_system_1_sm_master_DW.StateSpace_PWORK.Y_SWITCH;
      int_T *switchTypes = (int_T*)
        my_pv_system_1_sm_master_DW.StateSpace_PWORK.SWITCH_TYPES;
      int_T *idxOutSw = (int_T*)
        my_pv_system_1_sm_master_DW.StateSpace_PWORK.IDX_OUT_SW;
      int_T *switch_status_init = (int_T*)
        my_pv_system_1_sm_master_DW.StateSpace_PWORK.SWITCH_STATUS_INIT;

      /* Initialize work vectors */
      switch_status[0] = 0;
      switch_status_init[0] = 0;
      gState[0] = (int_T) 0.0;
      yswitch[0] = 1/0.001;
      switchTypes[0] = (int_T)7.0;
      idxOutSw[0] = ((int_T)0.0) - 1;
      switch_status[1] = 0;
      switch_status_init[1] = 0;
      gState[1] = (int_T) 0.0;
      yswitch[1] = 1/0.001;
      switchTypes[1] = (int_T)7.0;
      idxOutSw[1] = ((int_T)0.0) - 1;
      switch_status[2] = 0;
      switch_status_init[2] = 0;
      gState[2] = (int_T) 0.0;
      yswitch[2] = 1/0.001;
      switchTypes[2] = (int_T)7.0;
      idxOutSw[2] = ((int_T)0.0) - 1;
      switch_status[3] = 0;
      switch_status_init[3] = 0;
      gState[3] = (int_T) 0.0;
      yswitch[3] = 1/0.001;
      switchTypes[3] = (int_T)7.0;
      idxOutSw[3] = ((int_T)0.0) - 1;
    }
  }

  /* Level2 S-Function Block: '<S25>/S-Function' (RECV_Param) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[5];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* InitializeConditions for Delay: '<S9>/Delay' */
  my_pv_system_1_sm_master_DW.Delay_DSTATE =
    my_pv_system_1_sm_master_P.Delay_InitialCondition;
}

/* Model terminate function */
static void my_pv_system_1_sm_master_terminate(void)
{
  /* Terminate for Atomic SubSystem: '<S13>/Subsystem5' */

  /* Level2 S-Function Block: '<S46>/S-Function' (send_rt) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[0];
    sfcnTerminate(rts);
  }

  /* End of Terminate for SubSystem: '<S13>/Subsystem5' */

  /* Level2 S-Function Block: '<S44>/S-Function' (OP_SEND) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[1];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S2>/OpMonitor' (opmonitor) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[2];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S45>/S-Function' (OP_SEND) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[3];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S47>/S-Function' (recv_rt) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[4];
    sfcnTerminate(rts);
  }

  /* S-Function block: <S50>/State-Space */
  {
    /* Free memory */
    free(my_pv_system_1_sm_master_DW.StateSpace_PWORK.AS);
    free(my_pv_system_1_sm_master_DW.StateSpace_PWORK.BS);
    free(my_pv_system_1_sm_master_DW.StateSpace_PWORK.CS);
    free(my_pv_system_1_sm_master_DW.StateSpace_PWORK.DS);
    free(my_pv_system_1_sm_master_DW.StateSpace_PWORK.DX_COL);
    free(my_pv_system_1_sm_master_DW.StateSpace_PWORK.TMP2);
    free(my_pv_system_1_sm_master_DW.StateSpace_PWORK.BD_COL);
    free(my_pv_system_1_sm_master_DW.StateSpace_PWORK.TMP1);
    free(my_pv_system_1_sm_master_DW.StateSpace_PWORK.XTMP);

    /*
     * Circuit has switches*/
    free(my_pv_system_1_sm_master_DW.StateSpace_PWORK.CHOPPER);
    free(my_pv_system_1_sm_master_DW.StateSpace_PWORK.G_STATE);
    free(my_pv_system_1_sm_master_DW.StateSpace_PWORK.SWITCH_STATUS);
    free(my_pv_system_1_sm_master_DW.StateSpace_PWORK.SW_CHG);
    free(my_pv_system_1_sm_master_DW.StateSpace_PWORK.SWITCH_STATUS_INIT);
  }

  /* Level2 S-Function Block: '<S25>/S-Function' (RECV_Param) */
  {
    SimStruct *rts = my_pv_system_1_sm_master_M->childSfunctions[5];
    sfcnTerminate(rts);
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/
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
    rtsiSetErrorStatusPtr(&my_pv_system_1_sm_master_M->solverInfo,
                          (&rtmGetErrorStatus(my_pv_system_1_sm_master_M)));
    rtsiSetRTModelPtr(&my_pv_system_1_sm_master_M->solverInfo,
                      my_pv_system_1_sm_master_M);
  }

  rtsiSetSimTimeStep(&my_pv_system_1_sm_master_M->solverInfo, MAJOR_TIME_STEP);
  rtsiSetSolverName(&my_pv_system_1_sm_master_M->solverInfo,"FixedStepDiscrete");
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

  rtmSetTFinal(my_pv_system_1_sm_master_M, -1);
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

  {
    int_T i;
    for (i = 0; i < 12; i++) {
      my_pv_system_1_sm_master_B.StateSpace_o1[i] = 0.0;
    }

    my_pv_system_1_sm_master_B.SFunction = 0.0;
    my_pv_system_1_sm_master_B.Sum = 0.0;
    my_pv_system_1_sm_master_B.IrradianceWm2 = 0.0;
    my_pv_system_1_sm_master_B.I_PV = 0.0;
    my_pv_system_1_sm_master_B.V_PV = 0.0;
    my_pv_system_1_sm_master_B.I_Diode = 0.0;
    my_pv_system_1_sm_master_B.Vinv = 0.0;
    my_pv_system_1_sm_master_B.Vdc = 0.0;
    my_pv_system_1_sm_master_B.Igrid = 0.0;
    my_pv_system_1_sm_master_B.Vgrid = 0.0;
    my_pv_system_1_sm_master_B.Ig = 0.0;
    my_pv_system_1_sm_master_B.OpMonitor_o1 = 0.0;
    my_pv_system_1_sm_master_B.OpMonitor_o2 = 0.0;
    my_pv_system_1_sm_master_B.OpMonitor_o3 = 0.0;
    my_pv_system_1_sm_master_B.OpMonitor_o4 = 0.0;
    my_pv_system_1_sm_master_B.UnitDelay[0] = 0.0;
    my_pv_system_1_sm_master_B.UnitDelay[1] = 0.0;
    my_pv_system_1_sm_master_B.UnitDelay[2] = 0.0;
    my_pv_system_1_sm_master_B.UnitDelay[3] = 0.0;
    my_pv_system_1_sm_master_B.SFunction_p[0] = 0.0;
    my_pv_system_1_sm_master_B.SFunction_p[1] = 0.0;
    my_pv_system_1_sm_master_B.SFunction_p[2] = 0.0;
    my_pv_system_1_sm_master_B.SFunction_p[3] = 0.0;
    my_pv_system_1_sm_master_B.AC = 0.0;
    my_pv_system_1_sm_master_B.I_Diode_k = 0.0;
    my_pv_system_1_sm_master_B.Delay1 = 0.0;
    my_pv_system_1_sm_master_B.StateSpace_o2[0] = 0.0;
    my_pv_system_1_sm_master_B.StateSpace_o2[1] = 0.0;
    my_pv_system_1_sm_master_B.StateSpace_o2[2] = 0.0;
    my_pv_system_1_sm_master_B.StateSpace_o2[3] = 0.0;
    my_pv_system_1_sm_master_B.donotdeletethisgain = 0.0;
    my_pv_system_1_sm_master_B.donotdeletethisgain_e = 0.0;
    my_pv_system_1_sm_master_B.donotdeletethisgain_g = 0.0;
    my_pv_system_1_sm_master_B.SFunction_d[0] = 0.0;
    my_pv_system_1_sm_master_B.SFunction_d[1] = 0.0;
    my_pv_system_1_sm_master_B.donotdeletethisgain_a = 0.0;
    my_pv_system_1_sm_master_B.Delay = 0.0;
    my_pv_system_1_sm_master_B.donotdeletethisgain_d = 0.0;
    my_pv_system_1_sm_master_B.donotdeletethisgain_db = 0.0;
    my_pv_system_1_sm_master_B.Ron[0] = 0.0;
    my_pv_system_1_sm_master_B.Ron[1] = 0.0;
    my_pv_system_1_sm_master_B.Ron[2] = 0.0;
    my_pv_system_1_sm_master_B.Ron[3] = 0.0;
    my_pv_system_1_sm_master_B.DataTypeConversion[0] = 0.0;
    my_pv_system_1_sm_master_B.DataTypeConversion[1] = 0.0;
    my_pv_system_1_sm_master_B.DataTypeConversion[2] = 0.0;
    my_pv_system_1_sm_master_B.DataTypeConversion[3] = 0.0;
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
    my_pv_system_1_sm_master_B.donotdeletethisgain_i = 0.0;
    my_pv_system_1_sm_master_B.donotdeletethisgain_de = 0.0;
    my_pv_system_1_sm_master_B.Iph = 0.0;
    my_pv_system_1_sm_master_B.Io = 0.0;
    my_pv_system_1_sm_master_B.u = 0.0;
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
    my_pv_system_1_sm_master_B.Product[0] = 0.0;
    my_pv_system_1_sm_master_B.Product[1] = 0.0;
    my_pv_system_1_sm_master_B.Product[2] = 0.0;
    my_pv_system_1_sm_master_B.Product[3] = 0.0;
  }

  /* parameters */
  my_pv_system_1_sm_master_M->ModelData.defaultParam = ((real_T *)
    &my_pv_system_1_sm_master_P);

  /* states (dwork) */
  my_pv_system_1_sm_master_M->ModelData.dwork = ((void *)
    &my_pv_system_1_sm_master_DW);
  (void) memset((void *)&my_pv_system_1_sm_master_DW, 0,
                sizeof(DW_my_pv_system_1_sm_master_T));
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[0] = 0.0;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[1] = 0.0;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[2] = 0.0;
  my_pv_system_1_sm_master_DW.UnitDelay_DSTATE[3] = 0.0;
  my_pv_system_1_sm_master_DW.Delay2_DSTATE = 0.0;
  my_pv_system_1_sm_master_DW.Delay1_DSTATE = 0.0;

  {
    int_T i;
    for (i = 0; i < 10; i++) {
      my_pv_system_1_sm_master_DW.StateSpace_DSTATE[i] = 0.0;
    }
  }

  my_pv_system_1_sm_master_DW.Delay_DSTATE = 0.0;
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
  my_pv_system_1_sm_master_DW.Memory1_1_PreviousInput = 0.0;
  my_pv_system_1_sm_master_DW.Memory1_2_PreviousInput = 0.0;
  my_pv_system_1_sm_master_DW.Memory1_3_PreviousInput = 0.0;
  my_pv_system_1_sm_master_DW.Memory1_4_PreviousInput = 0.0;
  my_pv_system_1_sm_master_DW.Memory1_5_PreviousInput = 0.0;

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

  my_pv_system_1_sm_master_M->Sizes.numSFcns = (6);

  /* register each child */
  {
    (void) memset((void *)
                  &my_pv_system_1_sm_master_M->NonInlinedSFcns.childSFunctions[0],
                  0,
                  6*sizeof(SimStruct));
    my_pv_system_1_sm_master_M->childSfunctions =
      (&my_pv_system_1_sm_master_M->NonInlinedSFcns.childSFunctionPtrs[0]);

    {
      int_T i;
      for (i = 0; i < 6; i++) {
        my_pv_system_1_sm_master_M->childSfunctions[i] =
          (&my_pv_system_1_sm_master_M->NonInlinedSFcns.childSFunctions[i]);
      }
    }

    /* Level2 S-Function Block: my_pv_system_1_sm_master/<S46>/S-Function (send_rt) */
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
          sfcnUPtrs[4] = &my_pv_system_1_sm_master_B.Vinv;
          sfcnUPtrs[5] = &my_pv_system_1_sm_master_B.Vdc;
          sfcnUPtrs[6] = &my_pv_system_1_sm_master_B.Igrid;
          sfcnUPtrs[7] = &my_pv_system_1_sm_master_B.Vgrid;
          sfcnUPtrs[8] = &my_pv_system_1_sm_master_B.Ig;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 9);
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
      ssSetIWork(rts, (int_T *) &my_pv_system_1_sm_master_DW.SFunction_IWORK_g[0]);

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
        ssSetDWork(rts, 0, &my_pv_system_1_sm_master_DW.SFunction_IWORK_g[0]);
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
      ssSetInputPortWidth(rts, 0, 9);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: my_pv_system_1_sm_master/<S44>/S-Function (OP_SEND) */
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
          sfcnUPtrs[4] = &my_pv_system_1_sm_master_B.Vinv;
          sfcnUPtrs[5] = &my_pv_system_1_sm_master_B.Vdc;
          sfcnUPtrs[6] = &my_pv_system_1_sm_master_B.Igrid;
          sfcnUPtrs[7] = &my_pv_system_1_sm_master_B.Vgrid;
          sfcnUPtrs[8] = &my_pv_system_1_sm_master_B.Ig;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 9);
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
                       my_pv_system_1_sm_master_P.SFunction_P1_Size_b);
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
      ssSetInputPortWidth(rts, 0, 9);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: my_pv_system_1_sm_master/<S2>/OpMonitor (opmonitor) */
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
            &my_pv_system_1_sm_master_B.OpMonitor_o1));
        }

        /* port 1 */
        {
          _ssSetOutputPortNumDimensions(rts, 1, 1);
          ssSetOutputPortWidth(rts, 1, 1);
          ssSetOutputPortSignal(rts, 1, ((real_T *)
            &my_pv_system_1_sm_master_B.OpMonitor_o2));
        }

        /* port 2 */
        {
          _ssSetOutputPortNumDimensions(rts, 2, 1);
          ssSetOutputPortWidth(rts, 2, 1);
          ssSetOutputPortSignal(rts, 2, ((real_T *)
            &my_pv_system_1_sm_master_B.OpMonitor_o3));
        }

        /* port 3 */
        {
          _ssSetOutputPortNumDimensions(rts, 3, 1);
          ssSetOutputPortWidth(rts, 3, 1);
          ssSetOutputPortSignal(rts, 3, ((real_T *)
            &my_pv_system_1_sm_master_B.OpMonitor_o4));
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

    /* Level2 S-Function Block: my_pv_system_1_sm_master/<S45>/S-Function (OP_SEND) */
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
          sfcnUPtrs[0] = &my_pv_system_1_sm_master_B.OpMonitor_o1;
          sfcnUPtrs[1] = &my_pv_system_1_sm_master_B.OpMonitor_o2;
          sfcnUPtrs[2] = &my_pv_system_1_sm_master_B.OpMonitor_o3;
          sfcnUPtrs[3] = &my_pv_system_1_sm_master_B.OpMonitor_o4;
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
                       my_pv_system_1_sm_master_P.SFunction_P1_Size_l);
      }

      /* work vectors */
      ssSetIWork(rts, (int_T *) &my_pv_system_1_sm_master_DW.SFunction_IWORK_j[0]);

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
        ssSetDWork(rts, 0, &my_pv_system_1_sm_master_DW.SFunction_IWORK_j[0]);
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

    /* Level2 S-Function Block: my_pv_system_1_sm_master/<S47>/S-Function (recv_rt) */
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
          ssSetOutputPortWidth(rts, 0, 4);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            my_pv_system_1_sm_master_B.SFunction_p));
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
                       my_pv_system_1_sm_master_P.SFunction_P1_Size_o);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_1_sm_master_P.SFunction_P2_Size_h);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       my_pv_system_1_sm_master_P.SFunction_P3_Size_j);
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

    /* Level2 S-Function Block: my_pv_system_1_sm_master/<S25>/S-Function (RECV_Param) */
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

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn5.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 2);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            my_pv_system_1_sm_master_B.SFunction_d));
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
          &my_pv_system_1_sm_master_M->NonInlinedSFcns.Sfcn5.params;
        ssSetSFcnParamsCount(rts, 2);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_1_sm_master_P.SFunction_P1_Size_bw);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_1_sm_master_P.SFunction_P2_Size_d);
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
  my_pv_system_1_sm_master_M->Sizes.numContStates = (0);/* Number of continuous states */
  my_pv_system_1_sm_master_M->Sizes.numY = (0);/* Number of model outputs */
  my_pv_system_1_sm_master_M->Sizes.numU = (0);/* Number of model inputs */
  my_pv_system_1_sm_master_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  my_pv_system_1_sm_master_M->Sizes.numSampTimes = (2);/* Number of sample times */
  my_pv_system_1_sm_master_M->Sizes.numBlocks = (133);/* Number of blocks */
  my_pv_system_1_sm_master_M->Sizes.numBlockIO = (51);/* Number of block outputs */
  my_pv_system_1_sm_master_M->Sizes.numBlockPrms = (141);/* Sum of parameter "widths" */
  return my_pv_system_1_sm_master_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
