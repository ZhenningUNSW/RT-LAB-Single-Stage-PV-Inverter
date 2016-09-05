/*
 * my_pv_system_4_ss_controller_data.c
 *
 * Code generation for model "my_pv_system_4_ss_controller".
 *
 * Model version              : 1.157
 * Simulink Coder version : 8.7 (R2014b) 08-Sep-2014
 * C source code generated on : Mon Sep 05 22:15:14 2016
 *
 * Target selection: rtlab_rtmodel.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "my_pv_system_4_ss_controller.h"
#include "my_pv_system_4_ss_controller_private.h"

/* Block parameters (auto storage) */
P_my_pv_system_4_ss_controller_T my_pv_system_4_ss_controller_P = {
  1.0,                                 /* Mask Parameter: PLL_AGC
                                        * Referenced by: '<S30>/Constant1'
                                        */
  2.0,                                 /* Mask Parameter: AlphaBetaZerotodq0_Alignment
                                        * Referenced by: '<S55>/Constant'
                                        */
  50.0,                                /* Mask Parameter: InverterControl_Fnom
                                        * Referenced by:
                                        *   '<S12>/Constant4'
                                        *   '<S19>/First cycle of simulation Id=0.92, Iq=0'
                                        */
  0.003,                               /* Mask Parameter: InverterControl_Increment_MPPT
                                        * Referenced by: '<S8>/Iph_3'
                                        */
  314.15926535897933,                  /* Mask Parameter: Discrete_Init
                                        * Referenced by: '<S32>/Discrete-Time Integrator'
                                        */
  1.0,                                 /* Mask Parameter: Discrete_Kd
                                        * Referenced by: '<S32>/Discrete Derivative '
                                        */
  6.6,                                 /* Mask Parameter: InverterControl_Ki_Ireg
                                        * Referenced by: '<S14>/Integral Gain'
                                        */
  200.0,                               /* Mask Parameter: InverterControl_Ki_VDCreg
                                        * Referenced by: '<S85>/Integral Gain'
                                        */
  180.0,                               /* Mask Parameter: Discrete_Kp
                                        * Referenced by: '<S32>/Kp4'
                                        */
  0.15,                                /* Mask Parameter: InverterControl_Kp_Ireg
                                        * Referenced by: '<S14>/Proportional Gain'
                                        */
  12.0,                                /* Mask Parameter: InverterControl_Kp_VDCreg
                                        * Referenced by: '<S85>/Proportional Gain'
                                        */
  -1.5,                                /* Mask Parameter: PI_LowerSaturationLimit
                                        * Referenced by: '<S14>/Saturate'
                                        */
  -1.5,                                /* Mask Parameter: PI_LowerSaturationLimit_p
                                        * Referenced by: '<S85>/Saturate'
                                        */

  /*  Mask Parameter: PWM_Generator_MinMax
   * Referenced by: '<S11>/Constant10'
   */
  { -1.0, 1.0 },
  3500.0,                              /* Mask Parameter: InverterControl_Pnom
                                        * Referenced by: '<S10>/A->pu'
                                        */
  5.0E-5,                              /* Mask Parameter: InverterControl_Ts_Control
                                        * Referenced by:
                                        *   '<S12>/Constant4'
                                        *   '<S22>/Gain'
                                        *   '<S25>/Gain'
                                        *   '<S44>/Gain'
                                        *   '<S38>/Gain'
                                        *   '<S41>/Gain'
                                        */
  1.5,                                 /* Mask Parameter: PI_UpperSaturationLimit
                                        * Referenced by: '<S14>/Saturate'
                                        */
  1.5,                                 /* Mask Parameter: PI_UpperSaturationLimit_f
                                        * Referenced by: '<S85>/Saturate'
                                        */
  425.0,                               /* Mask Parameter: InverterControl_Vdc_ref_Init
                                        * Referenced by:
                                        *   '<S4>/Vnom_dc1'
                                        *   '<S8>/Iph_'
                                        */
  400.0,                               /* Mask Parameter: InverterControl_Vnom_dc
                                        * Referenced by: '<S13>/Rtot_pu2'
                                        */
  240.0,                               /* Mask Parameter: InverterControl_Vnom_prim
                                        * Referenced by:
                                        *   '<S10>/A->pu'
                                        *   '<S10>/V->pu'
                                        *   '<S12>/Constant3'
                                        */
  1.0,                                 /* Mask Parameter: CompareToConstant_const
                                        * Referenced by: '<S57>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: CompareToConstant1_const
                                        * Referenced by: '<S58>/Constant'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S23>/Gain1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S26>/Gain1'
                                        */
  0.0,                                 /* Expression: sps.K2
                                        * Referenced by: '<S28>/Gain1'
                                        */
  0.0,                                 /* Expression: sps.K1
                                        * Referenced by: '<S28>/Gain'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S39>/Gain1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S42>/Gain1'
                                        */
  1.0,                                 /* Expression: [1]
                                        * Referenced by: '<S31>/Gain'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S35>/Gain1'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval
                                        * Referenced by: '<S38>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S38>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S38>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S38>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size
   * Referenced by: '<S40>/S-Function'
   */
  { 1.0, 1.0 },
  0.025050000000000003,                /* Expression: MaxDelay
                                        * Referenced by: '<S40>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size
   * Referenced by: '<S40>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S40>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size
   * Referenced by: '<S40>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S40>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size
   * Referenced by: '<S40>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S40>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S39>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S38>/Constant'
                                        */
  1.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S38>/Unit Delay1'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S35>/Gain3'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval_h
                                        * Referenced by: '<S41>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S41>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S41>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S41>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_a
   * Referenced by: '<S43>/S-Function'
   */
  { 1.0, 1.0 },
  0.025050000000000003,                /* Expression: MaxDelay
                                        * Referenced by: '<S43>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_k
   * Referenced by: '<S43>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S43>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_h
   * Referenced by: '<S43>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S43>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_l
   * Referenced by: '<S43>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S43>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S42>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S41>/Constant'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S41>/Unit Delay1'
                                        */
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<S35>/Rad->Deg.'
                                        */
  0.0,                                 /* Expression: inf
                                        * Referenced by: '<S31>/Saturation'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S31>/Saturation'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S45>/Gain1'
                                        */

  /*  Expression: [0.92 0]
   * Referenced by: '<S19>/Constant'
   */
  { 0.92, 0.0 },

  /*  Expression: [0,0]
   * Referenced by: '<S59>/dq'
   */
  { 0.0, 0.0 },

  /*  Expression: [0,0]
   * Referenced by: '<S60>/dq'
   */
  { 0.0, 0.0 },

  /*  Computed Parameter: SFunction_P1_Size_b
   * Referenced by: '<S87>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: dest
                                        * Referenced by: '<S87>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_d
   * Referenced by: '<S87>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: priority2
                                        * Referenced by: '<S87>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_d
   * Referenced by: '<S87>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: st
                                        * Referenced by: '<S87>/S-Function'
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

  /*  Computed Parameter: SFunction_P1_Size_c
   * Referenced by: '<S88>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: src
                                        * Referenced by: '<S88>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_kw
   * Referenced by: '<S88>/S-Function'
   */
  { 1.0, 1.0 },
  7.0,                                 /* Expression: Data_width
                                        * Referenced by: '<S88>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_n
   * Referenced by: '<S88>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: st
                                        * Referenced by: '<S88>/S-Function'
                                        */
  50.0,                                /* Expression: sps.Finit
                                        * Referenced by: '<S30>/Unit Delay'
                                        */
  70.0,                                /* Expression: 70
                                        * Referenced by: '<S19>/avoid division by zero'
                                        */
  40.0,                                /* Expression: 40
                                        * Referenced by: '<S19>/avoid division by zero'
                                        */
  0.25,                                /* Expression: 1/4
                                        * Referenced by: '<S19>/Gain'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_f
   * Referenced by: '<S56>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: MaxDelay
                                        * Referenced by: '<S56>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_j
   * Referenced by: '<S56>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S56>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_h5
   * Referenced by: '<S56>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S56>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_o
   * Referenced by: '<S56>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S56>/S-Function'
                                        */
  5.0E-5,                              /* Computed Parameter: DiscreteTimeIntegrator_gainval
                                        * Referenced by: '<S30>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: sps.Phase_Init*pi/180
                                        * Referenced by: '<S30>/Discrete-Time Integrator'
                                        */
  6.2831853071795862,                  /* Expression: 2*pi
                                        * Referenced by: '<S30>/Constant4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S19>/First cycle of simulation Id=0.92, Iq=0'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S19>/First cycle of simulation Id=0.92, Iq=0'
                                        */
  0.060141321173410534,                /* Expression: sps.D
                                        * Referenced by: '<S51>/D*u(k)'
                                        */

  /*  Expression: sps.x0(1,:)
   * Referenced by: '<S51>/Delay_x1'
   */
  { 0.00010656152949374528, 0.0 },
  7937.8816629464736,                  /* Expression: sps.C11
                                        * Referenced by: '<S54>/C11'
                                        */

  /*  Expression: sps.x0(2,:)
   * Referenced by: '<S51>/Delay_x2'
   */
  { 0.0, 0.0 },
  0.12028264234682107,                 /* Expression: sps.C12
                                        * Referenced by: '<S54>/C12'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S4>/Unit Delay2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S4>/Iq_ref'
                                        */
  5.0E-5,                              /* Computed Parameter: Integrator_gainval
                                        * Referenced by: '<S14>/Integrator'
                                        */
  0.0,                                 /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S14>/Integrator'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S15>/Gain1'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval_o
                                        * Referenced by: '<S22>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S22>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S22>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S22>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_fw
   * Referenced by: '<S24>/S-Function'
   */
  { 1.0, 1.0 },
  0.022272222222222225,                /* Expression: MaxDelay
                                        * Referenced by: '<S24>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_a
   * Referenced by: '<S24>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S24>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_l
   * Referenced by: '<S24>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S24>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_n
   * Referenced by: '<S24>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S24>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S23>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S22>/Constant'
                                        */
  1.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S22>/Unit Delay1'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S15>/Gain3'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval_f
                                        * Referenced by: '<S25>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S25>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S25>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S25>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_l
   * Referenced by: '<S27>/S-Function'
   */
  { 1.0, 1.0 },
  0.022272222222222225,                /* Expression: MaxDelay
                                        * Referenced by: '<S27>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_f
   * Referenced by: '<S27>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S27>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_f
   * Referenced by: '<S27>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S27>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_c
   * Referenced by: '<S27>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S27>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S26>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S25>/Constant'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S25>/Unit Delay1'
                                        */
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<S15>/Rad->Deg.'
                                        */
  0.017453292519943295,                /* Expression: pi/180
                                        * Referenced by: '<S10>/to-rad'
                                        */
  0.002,                               /* Expression: RLff(1)
                                        * Referenced by: '<S7>/Rff '
                                        */
  0.2,                                 /* Expression: RLff(2)
                                        * Referenced by: '<S7>/Lff  '
                                        */
  0.002,                               /* Expression: RLff(1)
                                        * Referenced by: '<S7>/Rff'
                                        */
  0.2,                                 /* Expression: RLff(2)
                                        * Referenced by: '<S7>/Lff'
                                        */
  1.5,                                 /* Expression: 1.5
                                        * Referenced by: '<S7>/Saturation'
                                        */
  -1.5,                                /* Expression: -1.5
                                        * Referenced by: '<S7>/Saturation'
                                        */
  450.0,                               /* Expression: Limits_MPPT(1)
                                        * Referenced by: '<S8>/Iph_1'
                                        */
  375.0,                               /* Expression: Limits_MPPT(2)
                                        * Referenced by: '<S8>/Iph_2'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S2>/MPPT_On'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval_j
                                        * Referenced by: '<S28>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S28>/Integ4'
                                        */
  0.02,                                /* Expression: sps.Delay
                                        * Referenced by: '<S28>/K1'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_p
   * Referenced by: '<S29>/S-Function'
   */
  { 1.0, 1.0 },
  0.0201,                              /* Expression: MaxDelay
                                        * Referenced by: '<S29>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_ar
   * Referenced by: '<S29>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S29>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_e
   * Referenced by: '<S29>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S29>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_c0
   * Referenced by: '<S29>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S29>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S28>/Unit Delay'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S28>/K2'
                                        */
  425.0,                               /* Expression: sps.Vinit
                                        * Referenced by: '<S28>/Unit Delay1'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval_jy
                                        * Referenced by: '<S44>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S44>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S44>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S44>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_o
   * Referenced by: '<S46>/S-Function'
   */
  { 1.0, 1.0 },
  0.025050000000000003,                /* Expression: MaxDelay
                                        * Referenced by: '<S46>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_h
   * Referenced by: '<S46>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S46>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_j
   * Referenced by: '<S46>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S46>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_d
   * Referenced by: '<S46>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S46>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S45>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S44>/Constant'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S44>/Unit Delay1'
                                        */

  /*  Expression: [ TcD  Ts-TcD ]
   * Referenced by: '<S32>/Discrete Derivative '
   */
  { 0.0001, -5.0E-5 },
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S32>/Discrete Derivative '
                                        */
  0.16,                                /* Computed Parameter: DiscreteTimeIntegrator_gainva_j
                                        * Referenced by: '<S32>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: Par_Limits(1)
                                        * Referenced by: '<S32>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: Par_Limits(2)
                                        * Referenced by: '<S32>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: Par_Limits(1)
                                        * Referenced by: '<S32>/Saturation1'
                                        */
  0.0,                                 /* Expression: Par_Limits(2)
                                        * Referenced by: '<S32>/Saturation1'
                                        */
  0.15915494309189535,                 /* Expression: 1/2/pi
                                        * Referenced by: '<S30>/Gain10'
                                        */
  0.00060000000000000006,              /* Computed Parameter: RateLimiter_RisingLim
                                        * Referenced by: '<S30>/Rate Limiter'
                                        */
  -0.00060000000000000006,             /* Computed Parameter: RateLimiter_FallingLim
                                        * Referenced by: '<S30>/Rate Limiter'
                                        */
  50.0,                                /* Expression: sps.Finit
                                        * Referenced by: '<S30>/Rate Limiter'
                                        */
  40.528056790910988,                  /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S47>/Delay_x1'
                                        */
  0.99996932795768223,                 /* Expression: sps.A11
                                        * Referenced by: '<S48>/A11'
                                        */
  0.0,                                 /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S47>/Delay_x2'
                                        */
  4.9723130917851353E-5,               /* Expression: sps.A12
                                        * Referenced by: '<S48>/A12'
                                        */
  -1.226881692709296,                  /* Expression: sps.A21
                                        * Referenced by: '<S48>/A21'
                                        */
  0.988925236714054,                   /* Expression: sps.A22
                                        * Referenced by: '<S48>/A22'
                                        */
  2.4861565458925677E-5,               /* Expression: sps.B11
                                        * Referenced by: '<S49>/B11'
                                        */
  0.994462618357027,                   /* Expression: sps.B21
                                        * Referenced by: '<S49>/B21'
                                        */
  1.53360211588662E-5,                 /* Expression: sps.D
                                        * Referenced by: '<S47>/D*u(k)'
                                        */
  1.233694313470147,                   /* Expression: sps.C11
                                        * Referenced by: '<S50>/C11'
                                        */
  3.06720423177324E-5,                 /* Expression: sps.C12
                                        * Referenced by: '<S50>/C12'
                                        */
  0.87971735765317893,                 /* Expression: sps.A11
                                        * Referenced by: '<S52>/A11'
                                        */
  2.8483338533391971E-5,               /* Expression: sps.A12
                                        * Referenced by: '<S52>/A12'
                                        */
  -4811.3056938728414,                 /* Expression: sps.A21
                                        * Referenced by: '<S52>/A21'
                                        */
  0.13933354133567871,                 /* Expression: sps.A22
                                        * Referenced by: '<S52>/A22'
                                        */
  1.4241669266695987E-5,               /* Expression: sps.B11
                                        * Referenced by: '<S53>/B11'
                                        */
  0.56966677066783933,                 /* Expression: sps.B21
                                        * Referenced by: '<S53>/B21'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S19>/Constant1'
                                        */
  6.25E-5,                             /* Expression: sps.Delay
                                        * Referenced by: '<S83>/Constant3'
                                        */
  0.00025,                             /* Expression: sps.Period
                                        * Referenced by: '<S83>/Constant1'
                                        */
  4000.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S83>/1\ib1'
                                        */

  /*  Expression: [0 .5 1]
   * Referenced by: '<S83>/Lookup Table'
   */
  { 0.0, 0.5, 1.0 },

  /*  Expression: [0 2 0]
   * Referenced by: '<S83>/Lookup Table'
   */
  { 0.0, 2.0, 0.0 },
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S83>/Constant2'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S61>/Gain1'
                                        */
  0.1684,                              /* Expression: 0.1684
                                        * Referenced by: '<S4>/Unit Delay'
                                        */
  -1.0,                                /* Expression: -1
                                        * Referenced by: '<S66>/Gain'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S12>/Constant2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S4>/Unit Delay3'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S12>/Gain1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S4>/Unit Delay1'
                                        */
  2.5E-5,                              /* Computed Parameter: Integrator_gainval_i
                                        * Referenced by: '<S85>/Integrator'
                                        */
  0.9226                               /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S85>/Integrator'
                                        */
};
