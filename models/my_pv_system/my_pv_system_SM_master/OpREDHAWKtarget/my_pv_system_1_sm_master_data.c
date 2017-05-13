/*
 * my_pv_system_1_sm_master_data.c
 *
 * Code generation for model "my_pv_system_1_sm_master".
 *
 * Model version              : 1.225
 * Simulink Coder version : 8.7 (R2014b) 08-Sep-2014
 * C source code generated on : Mon Feb 27 11:22:27 2017
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
                                        * Referenced by: '<S44>/AC'
                                        */
  50.0,                                /* Mask Parameter: Vs14400V_Frequency
                                        * Referenced by: '<S44>/AC'
                                        */
  0.0,                                 /* Mask Parameter: Vs14400V_Phase
                                        * Referenced by: '<S44>/AC'
                                        */
  0.0,                                 /* Mask Parameter: Tail_Tf
                                        * Referenced by:
                                        *   '<S28>/Constant2'
                                        *   '<S28>/-0.9//Tf'
                                        */
  0.0,                                 /* Mask Parameter: Tail_Tt
                                        * Referenced by:
                                        *   '<S28>/Constant2'
                                        *   '<S28>/0.1//Tt'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S25>/itail'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S25>/1'
                                        */
  0.0,                                 /* Computed Parameter: _Value_a
                                        * Referenced by: '<S25>/2'
                                        */
  5.0E-5,                              /* Computed Parameter: DiscreteTimeIntegrator_gainval
                                        * Referenced by: '<S25>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S25>/Discrete-Time Integrator'
                                        */
  0.9,                                 /* Expression: 0.9
                                        * Referenced by: '<S28>/Constant'
                                        */
  0.9,                                 /* Expression: 0.9
                                        * Referenced by: '<S28>/Saturation1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S28>/Saturation1'
                                        */
  0.1,                                 /* Expression: 0.1
                                        * Referenced by: '<S28>/Saturation2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S28>/Saturation2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S25>/Unit Delay'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S25>/Switch'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S29>/Gain'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S30>/Gain'
                                        */

  /*  Computed Parameter: SFunction_P1_Size
   * Referenced by: '<S50>/S-Function'
   */
  { 1.0, 1.0 },
  3.0,                                 /* Expression: dest
                                        * Referenced by: '<S50>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size
   * Referenced by: '<S50>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: priority2
                                        * Referenced by: '<S50>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size
   * Referenced by: '<S50>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: st
                                        * Referenced by: '<S50>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/S-Function1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S3>/Memory'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S3>/Memory'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S3>/Memory'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S3>/Memory'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S3>/Memory'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S3>/Memory'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S3>/Memory1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S3>/Memory1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S3>/Memory1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S3>/Memory1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S3>/Memory1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S3>/Memory1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S3>/Memory1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S3>/Memory1'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_i
   * Referenced by: '<S48>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: Acqu_group
                                        * Referenced by: '<S48>/S-Function'
                                        */

  /*  Computed Parameter: OpMonitor_P1_Size
   * Referenced by: '<S3>/OpMonitor'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: compute_time
                                        * Referenced by: '<S3>/OpMonitor'
                                        */

  /*  Computed Parameter: OpMonitor_P2_Size
   * Referenced by: '<S3>/OpMonitor'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: real_step
                                        * Referenced by: '<S3>/OpMonitor'
                                        */

  /*  Computed Parameter: OpMonitor_P3_Size
   * Referenced by: '<S3>/OpMonitor'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: idle_time
                                        * Referenced by: '<S3>/OpMonitor'
                                        */

  /*  Computed Parameter: OpMonitor_P4_Size
   * Referenced by: '<S3>/OpMonitor'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: nb_overruns
                                        * Referenced by: '<S3>/OpMonitor'
                                        */

  /*  Computed Parameter: OpMonitor_P5_Size
   * Referenced by: '<S3>/OpMonitor'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: user_time
                                        * Referenced by: '<S3>/OpMonitor'
                                        */

  /*  Computed Parameter: OpMonitor_P6_Size
   * Referenced by: '<S3>/OpMonitor'
   */
  { 1.0, 32.0 },

  /*  Computed Parameter: OpMonitor_P6
   * Referenced by: '<S3>/OpMonitor'
   */
  { 109.0, 121.0, 95.0, 101.0, 118.0, 101.0, 110.0, 116.0, 95.0, 110.0, 97.0,
    109.0, 101.0, 44.0, 97.0, 110.0, 111.0, 116.0, 104.0, 101.0, 114.0, 95.0,
    101.0, 118.0, 101.0, 110.0, 116.0, 95.0, 110.0, 97.0, 109.0, 101.0 },

  /*  Computed Parameter: SFunction_P1_Size_p
   * Referenced by: '<S49>/S-Function'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: Acqu_group
                                        * Referenced by: '<S49>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S24>/0 4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S24>/Unit Delay'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_d
   * Referenced by: '<S51>/S-Function'
   */
  { 1.0, 1.0 },
  3.0,                                 /* Expression: src
                                        * Referenced by: '<S51>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_b
   * Referenced by: '<S51>/S-Function'
   */
  { 1.0, 1.0 },
  6.0,                                 /* Expression: Data_width
                                        * Referenced by: '<S51>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_d
   * Referenced by: '<S51>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: st
                                        * Referenced by: '<S51>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S24>/Not in ARTEMIS'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S44>/AC'
                                        */

  /*  Computed Parameter: StateSpace_P1_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 1.0, 19.0 },

  /*  Expression: sizeABCD
   * Referenced by: '<S54>/State-Space'
   */
  { 10.0, 10.0, 16.0, 10.0, 7.0, 16.0, 12.0, 10.0, 16.0, 12.0, 7.0, 16.0, 100.0,
    70.0, 29.0, 4.0, 5.0, 6.0, 12.0 },

  /*  Computed Parameter: StateSpace_P2_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 10.0, 1.0 },

  /*  Expression: x0
   * Referenced by: '<S54>/State-Space'
   */
  { 425.0, -2.5932275980290846, -1.0034857626241627E-6, 2.5932275980290864,
    -8.8987152882012675E-11, 6.6910267898712252, 23.460832962498461,
    -23.460832962498465, -0.52701625488899384, -0.17619081302360645 },

  /*  Computed Parameter: StateSpace_P3_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: sw_initcond_decimal
                                        * Referenced by: '<S54>/State-Space'
                                        */

  /*  Computed Parameter: StateSpace_P4_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 7.0, 1.0 },

  /*  Expression: u0
   * Referenced by: '<S54>/State-Space'
   */
  { 0.0, 0.0, 0.0, 0.0, -319.87441686059287, 8.05, 8.05 },

  /*  Computed Parameter: StateSpace_P5_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S54>/State-Space'
                                        */

  /*  Computed Parameter: StateSpace_P6_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 1.0, 11.0 },

  /*  Expression: holdtype
   * Referenced by: '<S54>/State-Space'
   */
  { 1.0, 1.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.25, 1000.0, 0.0, 1.0 },

  /*  Computed Parameter: StateSpace_P7_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 1.0, 42.0 },

  /*  Expression: Cspa_param
   * Referenced by: '<S54>/State-Space'
   */
  { 29.0, 4.0, 4.0, 2.0, 2.0, 6.0, 1.0, 1.0, 1.0, 4.0, 2.0, 1.0, 1.0, 1.0, 2.0,
    3.0, 4.0, 1.0, 2.0, 3.0, 4.0, 1.0, 3.0, 1.0, 3.0, 2.0, 4.0, 7.0, 8.0, 9.0,
    10.0, 1.0, 1.0, 1.0, 1.0, 2.0, 3.0, 4.0, 2.0, 4.0, 2.0, 1.0 },

  /*  Computed Parameter: StateSpace_P8_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 1.0, 17.0 },

  /*  Expression: Dspa_param
   * Referenced by: '<S54>/State-Space'
   */
  { 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 6.0, 7.0,
    6.0, 7.0 },

  /*  Computed Parameter: StateSpace_P9_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: Bspa_param
                                        * Referenced by: '<S54>/State-Space'
                                        */

  /*  Computed Parameter: StateSpace_P10_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: nb_sw
                                        * Referenced by: '<S54>/State-Space'
                                        */

  /*  Computed Parameter: StateSpace_P11_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: matfile
                                        * Referenced by: '<S54>/State-Space'
                                        */

  /*  Computed Parameter: StateSpace_P12_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 1.0, 1.0 },
  5.0,                                 /* Expression: methode
                                        * Referenced by: '<S54>/State-Space'
                                        */

  /*  Computed Parameter: StateSpace_P13_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: subnum
                                        * Referenced by: '<S54>/State-Space'
                                        */

  /*  Computed Parameter: StateSpace_P14_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 4.0, 1.0 },

  /*  Expression: Rswitch
   * Referenced by: '<S54>/State-Space'
   */
  { 0.001, 0.001, 0.001, 0.001 },

  /*  Computed Parameter: StateSpace_P15_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: dyna
                                        * Referenced by: '<S54>/State-Space'
                                        */

  /*  Computed Parameter: StateSpace_P16_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 1.0, 3.0 },

  /*  Expression: FOHinput
   * Referenced by: '<S54>/State-Space'
   */
  { 5.0, 6.0, 7.0 },

  /*  Computed Parameter: StateSpace_P17_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 1.0, 4.0 },

  /*  Expression: TypeSwitch
   * Referenced by: '<S54>/State-Space'
   */
  { 7.0, 7.0, 7.0, 7.0 },

  /*  Computed Parameter: StateSpace_P18_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 10.0, 10.0 },

  /*  Expression: Aart
   * Referenced by: '<S54>/State-Space'
   */
  { -0.079160495561056213, 114.52130096197894, 114.52130096197894,
    114.52130096197892, 0.039580247780528106, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -5.7261274116030477E+7, 5.7260604233292848E+7, -5.7260123111382917E+7,
    -1.2499999999999999E+8, 0.0, -54.721483276631467, 54.554154773111371,
    -682.83617057677736, 0.4097017023460664, 0.0, -502.38395577669144,
    -2.2904310807795495E+8, 502.38395574688911, 0.0, 41344.773369914241, 0.0,
    0.0, 0.0, 0.0, 0.0, -5.7259620727427147E+7, 5.7261106617248647E+7,
    -5.7261776499986216E+7, -1.2499999999999999E+8, -41344.773369914241,
    54.554154773111371, -54.721483276631467, 682.83617057677736,
    -0.4097017023460664, 0.0, 229.04260192395787, 229.04260192395787,
    229.04260192395785, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -229.04260192395785,
    -229.04260192395787, 229.04260192395782, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -574.05647937438619, -0.43917622739078865, 573.17812691960467, 0.0, 0.0,
    -54.721483276631467, 54.554154773111371, -682.83617057677736,
    0.4097017023460664, 0.0, 573.17812691960455, -0.43917622739081708,
    -574.0564793743863, 0.0, 0.0, 54.554154773111371, -54.721483276631467,
    682.83617057677736, -0.4097017023460664, 0.0, -68820.3123151764, 0.0,
    68820.312315176416, 0.0, 0.0, -6555.2272375370621, 6555.2272375370621,
    -82175.25796264915, 49.229756553903329, 0.0, 68820.3123151764, 0.0,
    -68820.312315176416, 0.0, 0.0, 6555.2272375370621, -6555.2272375370621,
    82049.594256505559, -49.229756553903329 },

  /*  Computed Parameter: StateSpace_P19_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 10.0, 7.0 },

  /*  Expression: Bart
   * Referenced by: '<S54>/State-Space'
   */
  { -166.66655555562963, 5.7260650480989471E+7, -5.7260650480989479E+7,
    5.7260650480989456E+7, 83.333277777814814, 0.0, 0.0, 0.0, 0.0, 0.0,
    -166.66655555562963, -5.7260650480989471E+7, 5.7260650480989479E+7,
    -5.7260650480989456E+7, 83.333277777814814, 0.0, 0.0, 0.0, 0.0, 0.0,
    -166.66655555562963, 5.7260650480989471E+7, 1.7178195144296843E+8,
    5.7260650480989479E+7, 83.333277777814814, 0.0, 0.0, 0.0, 0.0, 0.0,
    -166.66655555562963, -5.7260650480989471E+7, -1.7178195144296843E+8,
    -5.7260650480989479E+7, 83.333277777814814, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.5451282604018992, 0.0,
    -333.06052678350733, 0.0, 0.0, 0.0, 166.53026339175366, 0.0, 0.0, 0.0, 0.0,
    0.0, 333.06052678350733, 0.0, 0.0, 0.0, -166.53026339175366, 0.0, 0.0, 0.0,
    0.0, 0.0 },

  /*  Computed Parameter: StateSpace_P20_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 12.0, 10.0 },

  /*  Expression: Cart
   * Referenced by: '<S54>/State-Space'
   */
  { 0.5, 0.5, 0.5, 0.5, 0.0, 1.0, 0.99918224647157561, 1.0, 0.0, 0.0, 0.0,
    -0.00023648164500415973, 500000.0, -500000.0, 0.0, 0.0, 2.504413145539782,
    0.0, 0.0, 0.0, -500000.0, 1.0, 1.0, 0.0, -500000.0, 500000.0, 500000.0,
    -500000.0, 0.0, 0.0, 0.0, 0.0, 1.0E+6, 0.0, 0.0, 0.0, 500000.0, -500000.0,
    0.0, 0.0, -2.504413145539782, 0.0, 0.0, 0.0, -500000.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    2.504413145539782, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -2.504413145539782, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    300.46948356806018, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -300.46948356806018, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /*  Computed Parameter: StateSpace_P21_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 12.0, 7.0 },

  /*  Expression: Dart
   * Referenced by: '<S54>/State-Space'
   */
  { -500000.0, 500000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 500000.0, 0.0, 0.0, 0.0,
    500000.0, -500000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -500000.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -500000.0, 500000.0, 0.0, 0.0, 0.0, 0.0, -500000.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 500000.0, -500000.0, 0.0, 0.0, 0.0, 0.0, 500000.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -3.4551722082987086, 0.0, 0.0, 0.0, 0.0, -0.99918224647157561,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.4551722082987086, 0.0, 0.0, 0.0, 0.0,
    0.99918224647157561 },

  /*  Computed Parameter: StateSpace_P22_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: Opaldtc
                                        * Referenced by: '<S54>/State-Space'
                                        */

  /*  Computed Parameter: StateSpace_P23_Size
   * Referenced by: '<S54>/State-Space'
   */
  { 1.0, 3.0 },

  /*  Expression: Thy_ind_src
   * Referenced by: '<S54>/State-Space'
   */
  { 3.0, 0.0, 0.0 },
  1000.0,                              /* Expression: 1./Ron
                                        * Referenced by: '<S24>/1//Ron'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S24>/Switch'
                                        */
  0.0,                                 /* Expression: inf
                                        * Referenced by: '<S24>/Saturation'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S24>/Saturation'
                                        */

  /*  Expression: Vf_SwitchOn
   * Referenced by: '<S26>/Vf Devices & Clamping Diodes'
   */
  { 0.0, 0.0, 0.0, 0.0 },

  /*  Expression: Vf_DiodeOn
   * Referenced by: '<S26>/Vf Diodes'
   */
  { -0.0, -0.0, -0.0, -0.0 },
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S6>/do not delete this gain'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S29>/integrator'
                                        */
  0.02,                                /* Expression: 1./sps.Freq
                                        * Referenced by: '<S29>/Transport Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S29>/Transport Delay'
                                        */
  0.02,                                /* Expression: 1./sps.Freq
                                        * Referenced by: '<S29>/K1'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S29>/Memory'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S30>/integrator'
                                        */
  0.02,                                /* Expression: 1./sps.Freq
                                        * Referenced by: '<S30>/Transport Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S30>/Transport Delay'
                                        */
  0.02,                                /* Expression: 1./sps.Freq
                                        * Referenced by: '<S30>/K1'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S30>/Memory'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_f
   * Referenced by: '<S31>/S-Function'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: Data_width
                                        * Referenced by: '<S31>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_m
   * Referenced by: '<S31>/S-Function'
   */
  { 2.0, 1.0 },

  /*  Expression: InitialConditions
   * Referenced by: '<S31>/S-Function'
   */
  { 1000.0, 25.0 },
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S35>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S34>/do not delete this gain'
                                        */
  0.4,                                 /* Computed Parameter: RateLimiter_RisingLim
                                        * Referenced by: '<S12>/Rate Limiter'
                                        */
  -0.4,                                /* Computed Parameter: RateLimiter_FallingLim
                                        * Referenced by: '<S12>/Rate Limiter'
                                        */
  1000.0,                              /* Expression: 1000
                                        * Referenced by: '<S12>/Rate Limiter'
                                        */
  55.0,                                /* Expression: 55
                                        * Referenced by: '<S12>/Saturation'
                                        */
  5.0,                                 /* Expression: 5
                                        * Referenced by: '<S12>/Saturation'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S36>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S13>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S14>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S20>/do not delete this gain'
                                        */
  1.0                                  /* Expression: 1
                                        * Referenced by: '<S21>/do not delete this gain'
                                        */
};
