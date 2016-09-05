/*
 * my_pv_system_3_ss_calculation_data.c
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

/* Block parameters (auto storage) */
P_my_pv_system_3_ss_calculation_T my_pv_system_3_ss_calculation_P = {
  0.0,                                 /* Expression: sps.K2
                                        * Referenced by: '<S8>/Gain1'
                                        */
  0.0,                                 /* Expression: sps.K1
                                        * Referenced by: '<S8>/Gain'
                                        */
  0.0,                                 /* Expression: sps.K2
                                        * Referenced by: '<S10>/Gain1'
                                        */
  0.0,                                 /* Expression: sps.K1
                                        * Referenced by: '<S10>/Gain'
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

  /*  Computed Parameter: SFunction_P1_Size
   * Referenced by: '<S13>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: Acqu_group
                                        * Referenced by: '<S13>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_g
   * Referenced by: '<S14>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: src
                                        * Referenced by: '<S14>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size
   * Referenced by: '<S14>/S-Function'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: Data_width
                                        * Referenced by: '<S14>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size
   * Referenced by: '<S14>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: st
                                        * Referenced by: '<S14>/S-Function'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval
                                        * Referenced by: '<S8>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S8>/Integ4'
                                        */
  0.02,                                /* Expression: sps.Delay
                                        * Referenced by: '<S8>/K1'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_e
   * Referenced by: '<S9>/S-Function'
   */
  { 1.0, 1.0 },
  0.0201,                              /* Expression: MaxDelay
                                        * Referenced by: '<S9>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_e
   * Referenced by: '<S9>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S9>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_a
   * Referenced by: '<S9>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S9>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size
   * Referenced by: '<S9>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S9>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S8>/Unit Delay'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S8>/K2'
                                        */
  406.0,                               /* Expression: sps.Vinit
                                        * Referenced by: '<S8>/Unit Delay1'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval_n
                                        * Referenced by: '<S10>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S10>/Integ4'
                                        */
  0.02,                                /* Expression: sps.Delay
                                        * Referenced by: '<S10>/K1'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_ea
   * Referenced by: '<S11>/S-Function'
   */
  { 1.0, 1.0 },
  0.0201,                              /* Expression: MaxDelay
                                        * Referenced by: '<S11>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_l
   * Referenced by: '<S11>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S11>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_m
   * Referenced by: '<S11>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S11>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_k
   * Referenced by: '<S11>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S11>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S10>/Unit Delay'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S10>/K2'
                                        */
  8.05                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S10>/Unit Delay1'
                                        */
};
