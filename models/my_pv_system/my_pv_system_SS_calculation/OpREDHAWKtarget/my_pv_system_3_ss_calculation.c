/*
 * my_pv_system_3_ss_calculation.c
 *
 * Code generation for model "my_pv_system_3_ss_calculation".
 *
 * Model version              : 1.157
 * Simulink Coder version : 8.7 (R2014b) 08-Sep-2014
 * C source code generated on : Mon Sep 05 22:15:04 2016
 *
 * Target selection: rtlab_rtmodel.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "my_pv_system_3_ss_calculation.h"
#include "my_pv_system_3_ss_calculation_private.h"

/* Block signals (auto storage) */
B_my_pv_system_3_ss_calculation_T my_pv_system_3_ss_calculation_B;

/* Block states (auto storage) */
DW_my_pv_system_3_ss_calculation_T my_pv_system_3_ss_calculation_DW;

/* Real-time model */
RT_MODEL_my_pv_system_3_ss_calculation_T my_pv_system_3_ss_calculation_M_;
RT_MODEL_my_pv_system_3_ss_calculation_T *const my_pv_system_3_ss_calculation_M =
  &my_pv_system_3_ss_calculation_M_;

/* Model output function */
static void my_pv_system_3_ss_calculation_output(void)
{
  /* Memory: '<S1>/S-Function' */
  my_pv_system_3_ss_calculation_B.SFunction =
    my_pv_system_3_ss_calculation_DW.SFunction_PreviousInput;

  /* Sum: '<S1>/Sum' incorporates:
   *  Constant: '<S1>/S-Function1'
   */
  my_pv_system_3_ss_calculation_B.Sum =
    my_pv_system_3_ss_calculation_P.SFunction1_Value +
    my_pv_system_3_ss_calculation_B.SFunction;

  /* Stop: '<S1>/Stop Simulation' */
  if (my_pv_system_3_ss_calculation_B.Sum != 0.0) {
    rtmSetStopRequested(my_pv_system_3_ss_calculation_M, 1);
  }

  /* End of Stop: '<S1>/Stop Simulation' */

  /* Memory: '<S2>/Memory' */
  my_pv_system_3_ss_calculation_B.Memory[0] =
    my_pv_system_3_ss_calculation_DW.Memory_PreviousInput[0];
  my_pv_system_3_ss_calculation_B.Memory[1] =
    my_pv_system_3_ss_calculation_DW.Memory_PreviousInput[1];
  my_pv_system_3_ss_calculation_B.Memory[2] =
    my_pv_system_3_ss_calculation_DW.Memory_PreviousInput[2];

  /* Level2 S-Function Block: '<S13>/S-Function' (OP_SEND) */
  {
    SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[0];
    sfcnOutputs(rts, 0);
  }

  /* DigitalClock: '<S8>/Digital  Clock' */
  my_pv_system_3_ss_calculation_B.DigitalClock =
    my_pv_system_3_ss_calculation_M->Timing.t[0];

  /* Level2 S-Function Block: '<S14>/S-Function' (recv_rt) */
  {
    SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[1];
    sfcnOutputs(rts, 0);
  }

  /* DiscreteIntegrator: '<S8>/Integ4' */
  if (my_pv_system_3_ss_calculation_DW.Integ4_SYSTEM_ENABLE != 0) {
    my_pv_system_3_ss_calculation_B.Integ4 =
      my_pv_system_3_ss_calculation_DW.Integ4_DSTATE;
  } else {
    my_pv_system_3_ss_calculation_B.Integ4 =
      my_pv_system_3_ss_calculation_P.Integ4_gainval *
      my_pv_system_3_ss_calculation_B.SFunction_d[2] +
      my_pv_system_3_ss_calculation_DW.Integ4_DSTATE;
  }

  /* End of DiscreteIntegrator: '<S8>/Integ4' */

  /* Constant: '<S8>/K1' */
  my_pv_system_3_ss_calculation_B.K1 = my_pv_system_3_ss_calculation_P.K1_Value;

  /* Level2 S-Function Block: '<S9>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[2];
    sfcnOutputs(rts, 0);
  }

  /* UnitDelay: '<S8>/Unit Delay' */
  my_pv_system_3_ss_calculation_B.UnitDelay =
    my_pv_system_3_ss_calculation_DW.UnitDelay_DSTATE;

  /* RelationalOperator: '<S8>/Relational Operator' */
  my_pv_system_3_ss_calculation_B.RelationalOperator =
    (my_pv_system_3_ss_calculation_B.DigitalClock >=
     my_pv_system_3_ss_calculation_B.K1);

  /* UnitDelay: '<S8>/Unit Delay1' */
  my_pv_system_3_ss_calculation_B.UnitDelay1 =
    my_pv_system_3_ss_calculation_DW.UnitDelay1_DSTATE;

  /* Switch: '<S8>/Switch' */
  if (my_pv_system_3_ss_calculation_B.RelationalOperator) {
    /* Gain: '<S8>/Gain1' */
    my_pv_system_3_ss_calculation_B.Gain1_m =
      my_pv_system_3_ss_calculation_P.Gain1_Gain *
      my_pv_system_3_ss_calculation_B.SFunction_d[2];

    /* Gain: '<S8>/Gain' */
    my_pv_system_3_ss_calculation_B.Gain_a =
      my_pv_system_3_ss_calculation_P.Gain_Gain *
      my_pv_system_3_ss_calculation_B.UnitDelay;

    /* Sum: '<S8>/Sum1' */
    my_pv_system_3_ss_calculation_B.Correction_e =
      my_pv_system_3_ss_calculation_B.Gain1_m -
      my_pv_system_3_ss_calculation_B.Gain_a;

    /* Sum: '<S8>/Sum7' */
    my_pv_system_3_ss_calculation_B.Sum7_g =
      my_pv_system_3_ss_calculation_B.Integ4 -
      my_pv_system_3_ss_calculation_B.SFunction_a;

    /* Product: '<S8>/Product' incorporates:
     *  Constant: '<S8>/K2'
     */
    my_pv_system_3_ss_calculation_B.Mean_n =
      my_pv_system_3_ss_calculation_B.Sum7_g *
      my_pv_system_3_ss_calculation_P.K2_Value;

    /* Sum: '<S8>/Sum5' */
    my_pv_system_3_ss_calculation_B.Sum5_p =
      my_pv_system_3_ss_calculation_B.Mean_n +
      my_pv_system_3_ss_calculation_B.Correction_e;
    my_pv_system_3_ss_calculation_B.Switch =
      my_pv_system_3_ss_calculation_B.Sum5_p;
  } else {
    my_pv_system_3_ss_calculation_B.Switch =
      my_pv_system_3_ss_calculation_B.UnitDelay1;
  }

  /* End of Switch: '<S8>/Switch' */

  /* DigitalClock: '<S10>/Digital  Clock' */
  my_pv_system_3_ss_calculation_B.DigitalClock_h =
    my_pv_system_3_ss_calculation_M->Timing.t[0];

  /* DiscreteIntegrator: '<S10>/Integ4' */
  if (my_pv_system_3_ss_calculation_DW.Integ4_SYSTEM_ENABLE_f != 0) {
    my_pv_system_3_ss_calculation_B.Integ4_b =
      my_pv_system_3_ss_calculation_DW.Integ4_DSTATE_e;
  } else {
    my_pv_system_3_ss_calculation_B.Integ4_b =
      my_pv_system_3_ss_calculation_P.Integ4_gainval_n *
      my_pv_system_3_ss_calculation_B.SFunction_d[1] +
      my_pv_system_3_ss_calculation_DW.Integ4_DSTATE_e;
  }

  /* End of DiscreteIntegrator: '<S10>/Integ4' */

  /* Constant: '<S10>/K1' */
  my_pv_system_3_ss_calculation_B.K1_i =
    my_pv_system_3_ss_calculation_P.K1_Value_g;

  /* Level2 S-Function Block: '<S11>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[3];
    sfcnOutputs(rts, 0);
  }

  /* UnitDelay: '<S10>/Unit Delay' */
  my_pv_system_3_ss_calculation_B.UnitDelay_m =
    my_pv_system_3_ss_calculation_DW.UnitDelay_DSTATE_l;

  /* RelationalOperator: '<S10>/Relational Operator' */
  my_pv_system_3_ss_calculation_B.RelationalOperator_n =
    (my_pv_system_3_ss_calculation_B.DigitalClock_h >=
     my_pv_system_3_ss_calculation_B.K1_i);

  /* UnitDelay: '<S10>/Unit Delay1' */
  my_pv_system_3_ss_calculation_B.UnitDelay1_k =
    my_pv_system_3_ss_calculation_DW.UnitDelay1_DSTATE_p;

  /* Switch: '<S10>/Switch' */
  if (my_pv_system_3_ss_calculation_B.RelationalOperator_n) {
    /* Gain: '<S10>/Gain1' */
    my_pv_system_3_ss_calculation_B.Gain1 =
      my_pv_system_3_ss_calculation_P.Gain1_Gain_o *
      my_pv_system_3_ss_calculation_B.SFunction_d[1];

    /* Gain: '<S10>/Gain' */
    my_pv_system_3_ss_calculation_B.Gain =
      my_pv_system_3_ss_calculation_P.Gain_Gain_o *
      my_pv_system_3_ss_calculation_B.UnitDelay_m;

    /* Sum: '<S10>/Sum1' */
    my_pv_system_3_ss_calculation_B.Correction =
      my_pv_system_3_ss_calculation_B.Gain1 -
      my_pv_system_3_ss_calculation_B.Gain;

    /* Sum: '<S10>/Sum7' */
    my_pv_system_3_ss_calculation_B.Sum7 =
      my_pv_system_3_ss_calculation_B.Integ4_b -
      my_pv_system_3_ss_calculation_B.SFunction_f;

    /* Product: '<S10>/Product' incorporates:
     *  Constant: '<S10>/K2'
     */
    my_pv_system_3_ss_calculation_B.Mean = my_pv_system_3_ss_calculation_B.Sum7 *
      my_pv_system_3_ss_calculation_P.K2_Value_o;

    /* Sum: '<S10>/Sum5' */
    my_pv_system_3_ss_calculation_B.Sum5 = my_pv_system_3_ss_calculation_B.Mean
      + my_pv_system_3_ss_calculation_B.Correction;
    my_pv_system_3_ss_calculation_B.Switch_j =
      my_pv_system_3_ss_calculation_B.Sum5;
  } else {
    my_pv_system_3_ss_calculation_B.Switch_j =
      my_pv_system_3_ss_calculation_B.UnitDelay1_k;
  }

  /* End of Switch: '<S10>/Switch' */

  /* Product: '<S2>/Product' */
  my_pv_system_3_ss_calculation_B.Pmean = my_pv_system_3_ss_calculation_B.Switch
    * my_pv_system_3_ss_calculation_B.Switch_j;
}

/* Model update function */
static void my_pv_system_3_ss_calculation_update(void)
{
  /* Update for Memory: '<S1>/S-Function' */
  my_pv_system_3_ss_calculation_DW.SFunction_PreviousInput =
    my_pv_system_3_ss_calculation_B.Sum;

  /* Update for Memory: '<S2>/Memory' */
  my_pv_system_3_ss_calculation_DW.Memory_PreviousInput[0] =
    my_pv_system_3_ss_calculation_B.SFunction_d[0];
  my_pv_system_3_ss_calculation_DW.Memory_PreviousInput[1] =
    my_pv_system_3_ss_calculation_B.Switch;
  my_pv_system_3_ss_calculation_DW.Memory_PreviousInput[2] =
    my_pv_system_3_ss_calculation_B.Pmean;

  /* Update for DiscreteIntegrator: '<S8>/Integ4' */
  my_pv_system_3_ss_calculation_DW.Integ4_SYSTEM_ENABLE = 0U;
  my_pv_system_3_ss_calculation_DW.Integ4_DSTATE =
    my_pv_system_3_ss_calculation_P.Integ4_gainval *
    my_pv_system_3_ss_calculation_B.SFunction_d[2] +
    my_pv_system_3_ss_calculation_B.Integ4;

  /* Level2 S-Function Block: '<S9>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[2];
    sfcnUpdate(rts, 0);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Update for UnitDelay: '<S8>/Unit Delay' */
  my_pv_system_3_ss_calculation_DW.UnitDelay_DSTATE =
    my_pv_system_3_ss_calculation_B.SFunction_d[2];

  /* Update for UnitDelay: '<S8>/Unit Delay1' */
  my_pv_system_3_ss_calculation_DW.UnitDelay1_DSTATE =
    my_pv_system_3_ss_calculation_B.Switch;

  /* Update for DiscreteIntegrator: '<S10>/Integ4' */
  my_pv_system_3_ss_calculation_DW.Integ4_SYSTEM_ENABLE_f = 0U;
  my_pv_system_3_ss_calculation_DW.Integ4_DSTATE_e =
    my_pv_system_3_ss_calculation_P.Integ4_gainval_n *
    my_pv_system_3_ss_calculation_B.SFunction_d[1] +
    my_pv_system_3_ss_calculation_B.Integ4_b;

  /* Level2 S-Function Block: '<S11>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[3];
    sfcnUpdate(rts, 0);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Update for UnitDelay: '<S10>/Unit Delay' */
  my_pv_system_3_ss_calculation_DW.UnitDelay_DSTATE_l =
    my_pv_system_3_ss_calculation_B.SFunction_d[1];

  /* Update for UnitDelay: '<S10>/Unit Delay1' */
  my_pv_system_3_ss_calculation_DW.UnitDelay1_DSTATE_p =
    my_pv_system_3_ss_calculation_B.Switch_j;

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++my_pv_system_3_ss_calculation_M->Timing.clockTick0)) {
    ++my_pv_system_3_ss_calculation_M->Timing.clockTickH0;
  }

  my_pv_system_3_ss_calculation_M->Timing.t[0] =
    my_pv_system_3_ss_calculation_M->Timing.clockTick0 *
    my_pv_system_3_ss_calculation_M->Timing.stepSize0 +
    my_pv_system_3_ss_calculation_M->Timing.clockTickH0 *
    my_pv_system_3_ss_calculation_M->Timing.stepSize0 * 4294967296.0;
}

/* Model initialize function */
static void my_pv_system_3_ss_calculation_initialize(void)
{
  /* Start for Constant: '<S8>/K1' */
  my_pv_system_3_ss_calculation_B.K1 = my_pv_system_3_ss_calculation_P.K1_Value;

  /* Level2 S-Function Block: '<S9>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[2];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Start for Constant: '<S10>/K1' */
  my_pv_system_3_ss_calculation_B.K1_i =
    my_pv_system_3_ss_calculation_P.K1_Value_g;

  /* Level2 S-Function Block: '<S11>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[3];
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
  my_pv_system_3_ss_calculation_DW.SFunction_PreviousInput =
    my_pv_system_3_ss_calculation_P.SFunction_X0;

  /* InitializeConditions for Memory: '<S2>/Memory' */
  my_pv_system_3_ss_calculation_DW.Memory_PreviousInput[0] =
    my_pv_system_3_ss_calculation_P.Memory_X0;
  my_pv_system_3_ss_calculation_DW.Memory_PreviousInput[1] =
    my_pv_system_3_ss_calculation_P.Memory_X0;
  my_pv_system_3_ss_calculation_DW.Memory_PreviousInput[2] =
    my_pv_system_3_ss_calculation_P.Memory_X0;

  /* Level2 S-Function Block: '<S13>/S-Function' (OP_SEND) */
  {
    SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[0];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Level2 S-Function Block: '<S14>/S-Function' (recv_rt) */
  {
    SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[1];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* InitializeConditions for DiscreteIntegrator: '<S8>/Integ4' */
  my_pv_system_3_ss_calculation_DW.Integ4_DSTATE =
    my_pv_system_3_ss_calculation_P.Integ4_IC;

  /* Level2 S-Function Block: '<S9>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[2];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* InitializeConditions for UnitDelay: '<S8>/Unit Delay' */
  my_pv_system_3_ss_calculation_DW.UnitDelay_DSTATE =
    my_pv_system_3_ss_calculation_P.UnitDelay_InitialCondition;

  /* InitializeConditions for UnitDelay: '<S8>/Unit Delay1' */
  my_pv_system_3_ss_calculation_DW.UnitDelay1_DSTATE =
    my_pv_system_3_ss_calculation_P.UnitDelay1_InitialCondition;

  /* InitializeConditions for DiscreteIntegrator: '<S10>/Integ4' */
  my_pv_system_3_ss_calculation_DW.Integ4_DSTATE_e =
    my_pv_system_3_ss_calculation_P.Integ4_IC_p;

  /* Level2 S-Function Block: '<S11>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[3];
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* InitializeConditions for UnitDelay: '<S10>/Unit Delay' */
  my_pv_system_3_ss_calculation_DW.UnitDelay_DSTATE_l =
    my_pv_system_3_ss_calculation_P.UnitDelay_InitialCondition_g;

  /* InitializeConditions for UnitDelay: '<S10>/Unit Delay1' */
  my_pv_system_3_ss_calculation_DW.UnitDelay1_DSTATE_p =
    my_pv_system_3_ss_calculation_P.UnitDelay1_InitialCondition_h;

  /* Enable for DiscreteIntegrator: '<S8>/Integ4' */
  my_pv_system_3_ss_calculation_DW.Integ4_SYSTEM_ENABLE = 1U;

  /* Enable for DiscreteIntegrator: '<S10>/Integ4' */
  my_pv_system_3_ss_calculation_DW.Integ4_SYSTEM_ENABLE_f = 1U;
}

/* Model terminate function */
static void my_pv_system_3_ss_calculation_terminate(void)
{
  /* Level2 S-Function Block: '<S13>/S-Function' (OP_SEND) */
  {
    SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[0];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S14>/S-Function' (recv_rt) */
  {
    SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[1];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S9>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[2];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S11>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[3];
    sfcnTerminate(rts);
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/
void MdlOutputs(int_T tid)
{
  my_pv_system_3_ss_calculation_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  my_pv_system_3_ss_calculation_update();
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
  my_pv_system_3_ss_calculation_initialize();
}

void MdlTerminate(void)
{
  my_pv_system_3_ss_calculation_terminate();
}

/* Registration function */
RT_MODEL_my_pv_system_3_ss_calculation_T *my_pv_system_3_ss_calculation(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)my_pv_system_3_ss_calculation_M, 0,
                sizeof(RT_MODEL_my_pv_system_3_ss_calculation_T));
  rtsiSetSolverName(&my_pv_system_3_ss_calculation_M->solverInfo,
                    "FixedStepDiscrete");
  my_pv_system_3_ss_calculation_M->solverInfoPtr =
    (&my_pv_system_3_ss_calculation_M->solverInfo);

  /* Initialize timing info */
  {
    int_T *mdlTsMap =
      my_pv_system_3_ss_calculation_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    my_pv_system_3_ss_calculation_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    my_pv_system_3_ss_calculation_M->Timing.sampleTimes =
      (&my_pv_system_3_ss_calculation_M->Timing.sampleTimesArray[0]);
    my_pv_system_3_ss_calculation_M->Timing.offsetTimes =
      (&my_pv_system_3_ss_calculation_M->Timing.offsetTimesArray[0]);

    /* task periods */
    my_pv_system_3_ss_calculation_M->Timing.sampleTimes[0] = (5.0E-5);

    /* task offsets */
    my_pv_system_3_ss_calculation_M->Timing.offsetTimes[0] = (0.0);
  }

  rtmSetTPtr(my_pv_system_3_ss_calculation_M,
             &my_pv_system_3_ss_calculation_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits =
      my_pv_system_3_ss_calculation_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    my_pv_system_3_ss_calculation_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(my_pv_system_3_ss_calculation_M, -1);
  my_pv_system_3_ss_calculation_M->Timing.stepSize0 = 5.0E-5;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    my_pv_system_3_ss_calculation_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(my_pv_system_3_ss_calculation_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(my_pv_system_3_ss_calculation_M->rtwLogInfo, (NULL));
    rtliSetLogT(my_pv_system_3_ss_calculation_M->rtwLogInfo, "");
    rtliSetLogX(my_pv_system_3_ss_calculation_M->rtwLogInfo, "");
    rtliSetLogXFinal(my_pv_system_3_ss_calculation_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(my_pv_system_3_ss_calculation_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(my_pv_system_3_ss_calculation_M->rtwLogInfo, 0);
    rtliSetLogMaxRows(my_pv_system_3_ss_calculation_M->rtwLogInfo, 1000);
    rtliSetLogDecimation(my_pv_system_3_ss_calculation_M->rtwLogInfo, 1);
    rtliSetLogY(my_pv_system_3_ss_calculation_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(my_pv_system_3_ss_calculation_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(my_pv_system_3_ss_calculation_M->rtwLogInfo, (NULL));
  }

  my_pv_system_3_ss_calculation_M->solverInfoPtr =
    (&my_pv_system_3_ss_calculation_M->solverInfo);
  my_pv_system_3_ss_calculation_M->Timing.stepSize = (5.0E-5);
  rtsiSetFixedStepSize(&my_pv_system_3_ss_calculation_M->solverInfo, 5.0E-5);
  rtsiSetSolverMode(&my_pv_system_3_ss_calculation_M->solverInfo,
                    SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  my_pv_system_3_ss_calculation_M->ModelData.blockIO = ((void *)
    &my_pv_system_3_ss_calculation_B);
  (void) memset(((void *) &my_pv_system_3_ss_calculation_B), 0,
                sizeof(B_my_pv_system_3_ss_calculation_T));

  {
    my_pv_system_3_ss_calculation_B.SFunction = 0.0;
    my_pv_system_3_ss_calculation_B.Sum = 0.0;
    my_pv_system_3_ss_calculation_B.Memory[0] = 0.0;
    my_pv_system_3_ss_calculation_B.Memory[1] = 0.0;
    my_pv_system_3_ss_calculation_B.Memory[2] = 0.0;
    my_pv_system_3_ss_calculation_B.DigitalClock = 0.0;
    my_pv_system_3_ss_calculation_B.SFunction_d[0] = 0.0;
    my_pv_system_3_ss_calculation_B.SFunction_d[1] = 0.0;
    my_pv_system_3_ss_calculation_B.SFunction_d[2] = 0.0;
    my_pv_system_3_ss_calculation_B.SFunction_d[3] = 0.0;
    my_pv_system_3_ss_calculation_B.Integ4 = 0.0;
    my_pv_system_3_ss_calculation_B.K1 = 0.0;
    my_pv_system_3_ss_calculation_B.SFunction_a = 0.0;
    my_pv_system_3_ss_calculation_B.UnitDelay = 0.0;
    my_pv_system_3_ss_calculation_B.UnitDelay1 = 0.0;
    my_pv_system_3_ss_calculation_B.Switch = 0.0;
    my_pv_system_3_ss_calculation_B.DigitalClock_h = 0.0;
    my_pv_system_3_ss_calculation_B.Integ4_b = 0.0;
    my_pv_system_3_ss_calculation_B.K1_i = 0.0;
    my_pv_system_3_ss_calculation_B.SFunction_f = 0.0;
    my_pv_system_3_ss_calculation_B.UnitDelay_m = 0.0;
    my_pv_system_3_ss_calculation_B.UnitDelay1_k = 0.0;
    my_pv_system_3_ss_calculation_B.Switch_j = 0.0;
    my_pv_system_3_ss_calculation_B.Pmean = 0.0;
    my_pv_system_3_ss_calculation_B.Gain1 = 0.0;
    my_pv_system_3_ss_calculation_B.Gain = 0.0;
    my_pv_system_3_ss_calculation_B.Correction = 0.0;
    my_pv_system_3_ss_calculation_B.Sum7 = 0.0;
    my_pv_system_3_ss_calculation_B.Mean = 0.0;
    my_pv_system_3_ss_calculation_B.Sum5 = 0.0;
    my_pv_system_3_ss_calculation_B.Gain1_m = 0.0;
    my_pv_system_3_ss_calculation_B.Gain_a = 0.0;
    my_pv_system_3_ss_calculation_B.Correction_e = 0.0;
    my_pv_system_3_ss_calculation_B.Sum7_g = 0.0;
    my_pv_system_3_ss_calculation_B.Mean_n = 0.0;
    my_pv_system_3_ss_calculation_B.Sum5_p = 0.0;
  }

  /* parameters */
  my_pv_system_3_ss_calculation_M->ModelData.defaultParam = ((real_T *)
    &my_pv_system_3_ss_calculation_P);

  /* states (dwork) */
  my_pv_system_3_ss_calculation_M->ModelData.dwork = ((void *)
    &my_pv_system_3_ss_calculation_DW);
  (void) memset((void *)&my_pv_system_3_ss_calculation_DW, 0,
                sizeof(DW_my_pv_system_3_ss_calculation_T));
  my_pv_system_3_ss_calculation_DW.Integ4_DSTATE = 0.0;
  my_pv_system_3_ss_calculation_DW.UnitDelay_DSTATE = 0.0;
  my_pv_system_3_ss_calculation_DW.UnitDelay1_DSTATE = 0.0;
  my_pv_system_3_ss_calculation_DW.Integ4_DSTATE_e = 0.0;
  my_pv_system_3_ss_calculation_DW.UnitDelay_DSTATE_l = 0.0;
  my_pv_system_3_ss_calculation_DW.UnitDelay1_DSTATE_p = 0.0;
  my_pv_system_3_ss_calculation_DW.SFunction_PreviousInput = 0.0;
  my_pv_system_3_ss_calculation_DW.Memory_PreviousInput[0] = 0.0;
  my_pv_system_3_ss_calculation_DW.Memory_PreviousInput[1] = 0.0;
  my_pv_system_3_ss_calculation_DW.Memory_PreviousInput[2] = 0.0;
  my_pv_system_3_ss_calculation_DW.SFunction_RWORK = 0.0;
  my_pv_system_3_ss_calculation_DW.SFunction_RWORK_f = 0.0;

  /* child S-Function registration */
  {
    RTWSfcnInfo *sfcnInfo =
      &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.sfcnInfo;
    my_pv_system_3_ss_calculation_M->sfcnInfo = (sfcnInfo);
    rtssSetErrorStatusPtr(sfcnInfo, (&rtmGetErrorStatus
      (my_pv_system_3_ss_calculation_M)));
    rtssSetNumRootSampTimesPtr(sfcnInfo,
      &my_pv_system_3_ss_calculation_M->Sizes.numSampTimes);
    my_pv_system_3_ss_calculation_M->NonInlinedSFcns.taskTimePtrs[0] =
      &(rtmGetTPtr(my_pv_system_3_ss_calculation_M)[0]);
    rtssSetTPtrPtr(sfcnInfo,
                   my_pv_system_3_ss_calculation_M->NonInlinedSFcns.taskTimePtrs);
    rtssSetTStartPtr(sfcnInfo, &rtmGetTStart(my_pv_system_3_ss_calculation_M));
    rtssSetTFinalPtr(sfcnInfo, &rtmGetTFinal(my_pv_system_3_ss_calculation_M));
    rtssSetTimeOfLastOutputPtr(sfcnInfo, &rtmGetTimeOfLastOutput
      (my_pv_system_3_ss_calculation_M));
    rtssSetStepSizePtr(sfcnInfo,
                       &my_pv_system_3_ss_calculation_M->Timing.stepSize);
    rtssSetStopRequestedPtr(sfcnInfo, &rtmGetStopRequested
      (my_pv_system_3_ss_calculation_M));
    rtssSetDerivCacheNeedsResetPtr(sfcnInfo,
      &my_pv_system_3_ss_calculation_M->ModelData.derivCacheNeedsReset);
    rtssSetZCCacheNeedsResetPtr(sfcnInfo,
      &my_pv_system_3_ss_calculation_M->ModelData.zCCacheNeedsReset);
    rtssSetBlkStateChangePtr(sfcnInfo,
      &my_pv_system_3_ss_calculation_M->ModelData.blkStateChange);
    rtssSetSampleHitsPtr(sfcnInfo,
                         &my_pv_system_3_ss_calculation_M->Timing.sampleHits);
    rtssSetPerTaskSampleHitsPtr(sfcnInfo,
      &my_pv_system_3_ss_calculation_M->Timing.perTaskSampleHits);
    rtssSetSimModePtr(sfcnInfo, &my_pv_system_3_ss_calculation_M->simMode);
    rtssSetSolverInfoPtr(sfcnInfo,
                         &my_pv_system_3_ss_calculation_M->solverInfoPtr);
  }

  my_pv_system_3_ss_calculation_M->Sizes.numSFcns = (4);

  /* register each child */
  {
    (void) memset((void *)
                  &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.childSFunctions
                  [0], 0,
                  4*sizeof(SimStruct));
    my_pv_system_3_ss_calculation_M->childSfunctions =
      (&my_pv_system_3_ss_calculation_M->NonInlinedSFcns.childSFunctionPtrs[0]);
    my_pv_system_3_ss_calculation_M->childSfunctions[0] =
      (&my_pv_system_3_ss_calculation_M->NonInlinedSFcns.childSFunctions[0]);
    my_pv_system_3_ss_calculation_M->childSfunctions[1] =
      (&my_pv_system_3_ss_calculation_M->NonInlinedSFcns.childSFunctions[1]);
    my_pv_system_3_ss_calculation_M->childSfunctions[2] =
      (&my_pv_system_3_ss_calculation_M->NonInlinedSFcns.childSFunctions[2]);
    my_pv_system_3_ss_calculation_M->childSfunctions[3] =
      (&my_pv_system_3_ss_calculation_M->NonInlinedSFcns.childSFunctions[3]);

    /* Level2 S-Function Block: my_pv_system_3_ss_calculation/<S13>/S-Function (OP_SEND) */
    {
      SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[0];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn0.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn0.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn0.sfcnTsMap;
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
                         &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.blkInfo2
                         [0]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_3_ss_calculation_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.methods2
                           [0]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.methods3
                           [0]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.statesInfo2
                         [0]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn0.inputPortInfo
          [0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn0.UPtrs0;
          sfcnUPtrs[0] = my_pv_system_3_ss_calculation_B.Memory;
          sfcnUPtrs[1] = &my_pv_system_3_ss_calculation_B.Memory[1];
          sfcnUPtrs[2] = &my_pv_system_3_ss_calculation_B.Memory[2];
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 3);
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "my_pv_system_3_ss_calculation/SS_calculation/rtlab_send_subsystem/Subsystem1/Send1/S-Function");
      ssSetRTModel(rts,my_pv_system_3_ss_calculation_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn0.params;
        ssSetSFcnParamsCount(rts, 1);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_3_ss_calculation_P.SFunction_P1_Size);
      }

      /* work vectors */
      ssSetIWork(rts, (int_T *)
                 &my_pv_system_3_ss_calculation_DW.SFunction_IWORK[0]);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn0.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn0.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* IWORK */
        ssSetDWorkWidth(rts, 0, 5);
        ssSetDWorkDataType(rts, 0,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_3_ss_calculation_DW.SFunction_IWORK[0]);
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
      ssSetInputPortWidth(rts, 0, 3);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: my_pv_system_3_ss_calculation/<S14>/S-Function (recv_rt) */
    {
      SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[1];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn1.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn1.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn1.sfcnTsMap;
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
                         &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.blkInfo2
                         [1]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_3_ss_calculation_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.methods2
                           [1]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.methods3
                           [1]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.statesInfo2
                         [1]);
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &my_pv_system_3_ss_calculation_M->
          NonInlinedSFcns.Sfcn1.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 4);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            my_pv_system_3_ss_calculation_B.SFunction_d));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "my_pv_system_3_ss_calculation/SS_calculation/zzzOpComm/Receive_1/S-Function");
      ssSetRTModel(rts,my_pv_system_3_ss_calculation_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn1.params;
        ssSetSFcnParamsCount(rts, 3);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_3_ss_calculation_P.SFunction_P1_Size_g);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_3_ss_calculation_P.SFunction_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       my_pv_system_3_ss_calculation_P.SFunction_P3_Size);
      }

      /* work vectors */
      ssSetPWork(rts, (void **)
                 &my_pv_system_3_ss_calculation_DW.SFunction_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn1.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn1.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_3_ss_calculation_DW.SFunction_PWORK);
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

    /* Level2 S-Function Block: my_pv_system_3_ss_calculation/<S9>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[2];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn2.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn2.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn2.sfcnTsMap;
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
                         &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.blkInfo2
                         [2]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_3_ss_calculation_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.methods2
                           [2]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.methods3
                           [2]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.statesInfo2
                         [2]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn2.inputPortInfo
          [0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn2.UPtrs0;
          sfcnUPtrs[0] = &my_pv_system_3_ss_calculation_B.Integ4;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn2.UPtrs1;
          sfcnUPtrs[0] = &my_pv_system_3_ss_calculation_B.K1;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &my_pv_system_3_ss_calculation_M->
          NonInlinedSFcns.Sfcn2.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &my_pv_system_3_ss_calculation_B.SFunction_a));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "my_pv_system_3_ss_calculation/SS_calculation/Mean/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,my_pv_system_3_ss_calculation_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn2.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_3_ss_calculation_P.SFunction_P1_Size_e);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_3_ss_calculation_P.SFunction_P2_Size_e);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       my_pv_system_3_ss_calculation_P.SFunction_P3_Size_a);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       my_pv_system_3_ss_calculation_P.SFunction_P4_Size);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *)
                 &my_pv_system_3_ss_calculation_DW.SFunction_RWORK);
      ssSetIWork(rts, (int_T *)
                 &my_pv_system_3_ss_calculation_DW.SFunction_IWORK_f);
      ssSetPWork(rts, (void **)
                 &my_pv_system_3_ss_calculation_DW.SFunction_PWORK_o);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn2.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn2.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_3_ss_calculation_DW.SFunction_RWORK);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &my_pv_system_3_ss_calculation_DW.SFunction_IWORK_f);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &my_pv_system_3_ss_calculation_DW.SFunction_PWORK_o);
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

    /* Level2 S-Function Block: my_pv_system_3_ss_calculation/<S11>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = my_pv_system_3_ss_calculation_M->childSfunctions[3];

      /* timing info */
      time_T *sfcnPeriod =
        my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn3.sfcnPeriod;
      time_T *sfcnOffset =
        my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn3.sfcnOffset;
      int_T *sfcnTsMap =
        my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn3.sfcnTsMap;
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
                         &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.blkInfo2
                         [3]);
      }

      ssSetRTWSfcnInfo(rts, my_pv_system_3_ss_calculation_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.methods2
                           [3]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.methods3
                           [3]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.statesInfo2
                         [3]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn3.inputPortInfo
          [0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn3.UPtrs0;
          sfcnUPtrs[0] = &my_pv_system_3_ss_calculation_B.Integ4_b;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn3.UPtrs1;
          sfcnUPtrs[0] = &my_pv_system_3_ss_calculation_B.K1_i;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &my_pv_system_3_ss_calculation_M->
          NonInlinedSFcns.Sfcn3.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &my_pv_system_3_ss_calculation_B.SFunction_f));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "my_pv_system_3_ss_calculation/SS_calculation/Mean1/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,my_pv_system_3_ss_calculation_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn3.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       my_pv_system_3_ss_calculation_P.SFunction_P1_Size_ea);
        ssSetSFcnParam(rts, 1, (mxArray*)
                       my_pv_system_3_ss_calculation_P.SFunction_P2_Size_l);
        ssSetSFcnParam(rts, 2, (mxArray*)
                       my_pv_system_3_ss_calculation_P.SFunction_P3_Size_m);
        ssSetSFcnParam(rts, 3, (mxArray*)
                       my_pv_system_3_ss_calculation_P.SFunction_P4_Size_k);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *)
                 &my_pv_system_3_ss_calculation_DW.SFunction_RWORK_f);
      ssSetIWork(rts, (int_T *)
                 &my_pv_system_3_ss_calculation_DW.SFunction_IWORK_n);
      ssSetPWork(rts, (void **)
                 &my_pv_system_3_ss_calculation_DW.SFunction_PWORK_k);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn3.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &my_pv_system_3_ss_calculation_M->NonInlinedSFcns.Sfcn3.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &my_pv_system_3_ss_calculation_DW.SFunction_RWORK_f);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &my_pv_system_3_ss_calculation_DW.SFunction_IWORK_n);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &my_pv_system_3_ss_calculation_DW.SFunction_PWORK_k);
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
  my_pv_system_3_ss_calculation_M->Sizes.numContStates = (0);/* Number of continuous states */
  my_pv_system_3_ss_calculation_M->Sizes.numY = (0);/* Number of model outputs */
  my_pv_system_3_ss_calculation_M->Sizes.numU = (0);/* Number of model inputs */
  my_pv_system_3_ss_calculation_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  my_pv_system_3_ss_calculation_M->Sizes.numSampTimes = (1);/* Number of sample times */
  my_pv_system_3_ss_calculation_M->Sizes.numBlocks = (40);/* Number of blocks */
  my_pv_system_3_ss_calculation_M->Sizes.numBlockIO = (33);/* Number of block outputs */
  my_pv_system_3_ss_calculation_M->Sizes.numBlockPrms = (55);/* Sum of parameter "widths" */
  return my_pv_system_3_ss_calculation_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
