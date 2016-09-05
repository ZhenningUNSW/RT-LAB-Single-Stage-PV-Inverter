/*
 * my_pv_system_1_sm_master_data.c
 *
 * Code generation for model "my_pv_system_1_sm_master".
 *
 * Model version              : 1.157
 * Simulink Coder version : 8.7 (R2014b) 08-Sep-2014
 * C source code generated on : Mon Sep 05 22:14:56 2016
 *
 * Target selection: rtlab_rtmodel.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "my_pv_system_1_sm_master.h"
#include "my_pv_system_1_sm_master_private.h"

/* Block parameters (auto storage) */
P_my_pv_system_1_sm_master_T my_pv_system_1_sm_master_P = {
  20364.675298172569,                  /* Mask Parameter: Vs14400V_Amplitude
                                        * Referenced by: '<S40>/AC'
                                        */
  50.0,                                /* Mask Parameter: Vs14400V_Frequency
                                        * Referenced by: '<S40>/AC'
                                        */
  0.0,                                 /* Mask Parameter: Vs14400V_Phase
                                        * Referenced by: '<S40>/AC'
                                        */
  0.0,                                 /* Mask Parameter: Tail_Tf
                                        * Referenced by:
                                        *   '<S22>/2'
                                        *   '<S24>/Constant2'
                                        *   '<S24>/-0.9//Tf'
                                        */
  0.0,                                 /* Mask Parameter: Tail_Tt
                                        * Referenced by:
                                        *   '<S22>/2'
                                        *   '<S24>/Constant2'
                                        *   '<S24>/0.1//Tt'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S22>/itail'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S22>/1'
                                        */
  5.0E-5,                              /* Computed Parameter: DiscreteTimeIntegrator_gainval
                                        * Referenced by: '<S22>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S22>/Discrete-Time Integrator'
                                        */
  0.9,                                 /* Expression: 0.9
                                        * Referenced by: '<S24>/Constant'
                                        */
  0.9,                                 /* Expression: 0.9
                                        * Referenced by: '<S24>/Saturation1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S24>/Saturation1'
                                        */
  0.1,                                 /* Expression: 0.1
                                        * Referenced by: '<S24>/Saturation2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S24>/Saturation2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S22>/Unit Delay'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S22>/Switch'
                                        */

  /*  Computed Parameter: SFunction_P1_Size
   * Referenced by: '<S45>/S-Function'
   */
  { 1.0, 1.0 },
  3.0,                                 /* Expression: dest
                                        * Referenced by: '<S45>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size
   * Referenced by: '<S45>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: priority2
                                        * Referenced by: '<S45>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size
   * Referenced by: '<S45>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: st
                                        * Referenced by: '<S45>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_k
   * Referenced by: '<S46>/S-Function'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: dest
                                        * Referenced by: '<S46>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_c
   * Referenced by: '<S46>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: priority2
                                        * Referenced by: '<S46>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_i
   * Referenced by: '<S46>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: st
                                        * Referenced by: '<S46>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S1>/S-Function1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S1>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory1'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_c
   * Referenced by: '<S44>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: Acqu_group
                                        * Referenced by: '<S44>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S21>/Unit Delay'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_cu
   * Referenced by: '<S47>/S-Function'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: src
                                        * Referenced by: '<S47>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_f
   * Referenced by: '<S47>/S-Function'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: Data_width
                                        * Referenced by: '<S47>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_b
   * Referenced by: '<S47>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: st
                                        * Referenced by: '<S47>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S40>/AC'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S9>/Delay2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S9>/Delay1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S5>/do not delete this gain'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_l
   * Referenced by: '<S25>/S-Function'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: Data_width
                                        * Referenced by: '<S25>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_m
   * Referenced by: '<S25>/S-Function'
   */
  { 2.0, 1.0 },

  /*  Expression: InitialConditions
   * Referenced by: '<S25>/S-Function'
   */
  { 1000.0, 25.0 },
  8.55,                                /* Expression: 8.55
                                        * Referenced by: '<S9>/Constant1'
                                        */
  523.6,                               /* Expression: 37.4*14
                                        * Referenced by: '<S9>/Constant2'
                                        */
  0.06,                                /* Expression: 0.06
                                        * Referenced by: '<S9>/Constant3'
                                        */
  840.0,                               /* Expression: 60*14
                                        * Referenced by: '<S9>/Constant4'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S29>/do not delete this gain'
                                        */
  523.6,                               /* Expression: 37.4*14
                                        * Referenced by: '<S9>/Delay'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S31>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S32>/do not delete this gain'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S21>/0 4'
                                        */
  1000.0,                              /* Expression: 1./Ron
                                        * Referenced by: '<S21>/1//Ron'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S21>/Switch'
                                        */
  0.0,                                 /* Expression: inf
                                        * Referenced by: '<S21>/Saturation'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S21>/Saturation'
                                        */

  /*  Expression: Vf_SwitchOn
   * Referenced by: '<S23>/Vf Devices & Clamping Diodes'
   */
  { 0.0, 0.0, 0.0, 0.0 },

  /*  Expression: Vf_DiodeOn
   * Referenced by: '<S23>/Vf Diodes'
   */
  { -0.0, -0.0, -0.0, -0.0 },
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S10>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S11>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S17>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S18>/do not delete this gain'
                                        */
  1U,                                  /* Computed Parameter: Delay2_DelayLength
                                        * Referenced by: '<S9>/Delay2'
                                        */
  1U,                                  /* Computed Parameter: Delay1_DelayLength
                                        * Referenced by: '<S9>/Delay1'
                                        */
  1U,                                  /* Computed Parameter: Delay_DelayLength
                                        * Referenced by: '<S9>/Delay'
                                        */
  0                                    /* Expression: Tf_sps>0 | Tt_sps>0
                                        * Referenced by: '<S21>/2'
                                        */
};
