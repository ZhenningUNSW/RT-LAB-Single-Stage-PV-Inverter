/*
 * my_pv_system_3_ss_controller_data.c
 *
 * Code generation for model "my_pv_system_3_ss_controller".
 *
 * Model version              : 1.225
 * Simulink Coder version : 8.7 (R2014b) 08-Sep-2014
 * C source code generated on : Mon Feb 27 11:22:43 2017
 *
 * Target selection: rtlab_rtmodel.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "my_pv_system_3_ss_controller.h"
#include "my_pv_system_3_ss_controller_private.h"

/* Block parameters (auto storage) */
P_my_pv_system_3_ss_controller_T my_pv_system_3_ss_controller_P = {
  1.0,                                 /* Mask Parameter: PLL_AGC
                                        * Referenced by: '<S32>/Constant1'
                                        */
  2.0,                                 /* Mask Parameter: AlphaBetaZerotodq0_Alignment
                                        * Referenced by: '<S57>/Constant'
                                        */
  50.0,                                /* Mask Parameter: InverterControl_Fnom
                                        * Referenced by:
                                        *   '<S14>/Constant4'
                                        *   '<S21>/First cycle of simulation Id=0.92, Iq=0'
                                        */
  0.003,                               /* Mask Parameter: InverterControl_Increment_MPPT
                                        * Referenced by: '<S10>/Iph_3'
                                        */
  314.15926535897933,                  /* Mask Parameter: Discrete_Init
                                        * Referenced by: '<S34>/Discrete-Time Integrator'
                                        */
  1.0,                                 /* Mask Parameter: Discrete_Kd
                                        * Referenced by: '<S34>/Discrete Derivative '
                                        */
  6.6,                                 /* Mask Parameter: InverterControl_Ki_Ireg
                                        * Referenced by: '<S16>/Integral Gain'
                                        */
  200.0,                               /* Mask Parameter: InverterControl_Ki_VDCreg
                                        * Referenced by: '<S87>/Integral Gain'
                                        */
  180.0,                               /* Mask Parameter: Discrete_Kp
                                        * Referenced by: '<S34>/Kp4'
                                        */
  0.15,                                /* Mask Parameter: InverterControl_Kp_Ireg
                                        * Referenced by: '<S16>/Proportional Gain'
                                        */
  12.0,                                /* Mask Parameter: InverterControl_Kp_VDCreg
                                        * Referenced by: '<S87>/Proportional Gain'
                                        */
  -1.5,                                /* Mask Parameter: PI_LowerSaturationLimit
                                        * Referenced by: '<S16>/Saturate'
                                        */
  -1.5,                                /* Mask Parameter: PI_LowerSaturationLimit_p
                                        * Referenced by: '<S87>/Saturate'
                                        */

  /*  Mask Parameter: PWM_Generator_MinMax
   * Referenced by: '<S13>/Constant10'
   */
  { -1.0, 1.0 },
  3500.0,                              /* Mask Parameter: InverterControl_Pnom
                                        * Referenced by: '<S12>/A->pu'
                                        */
  5.0E-5,                              /* Mask Parameter: InverterControl_Ts_Control
                                        * Referenced by:
                                        *   '<S14>/Constant4'
                                        *   '<S24>/Gain'
                                        *   '<S27>/Gain'
                                        *   '<S46>/Gain'
                                        *   '<S40>/Gain'
                                        *   '<S43>/Gain'
                                        */
  1.5,                                 /* Mask Parameter: PI_UpperSaturationLimit
                                        * Referenced by: '<S16>/Saturate'
                                        */
  1.5,                                 /* Mask Parameter: PI_UpperSaturationLimit_f
                                        * Referenced by: '<S87>/Saturate'
                                        */
  425.0,                               /* Mask Parameter: InverterControl_Vdc_ref_Init
                                        * Referenced by:
                                        *   '<S5>/Vnom_dc1'
                                        *   '<S10>/Iph_'
                                        */
  400.0,                               /* Mask Parameter: InverterControl_Vnom_dc
                                        * Referenced by: '<S15>/Rtot_pu2'
                                        */
  240.0,                               /* Mask Parameter: InverterControl_Vnom_prim
                                        * Referenced by:
                                        *   '<S12>/A->pu'
                                        *   '<S12>/V->pu'
                                        *   '<S14>/Constant3'
                                        */
  1.0,                                 /* Mask Parameter: CompareToConstant_const
                                        * Referenced by: '<S59>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: CompareToConstant1_const
                                        * Referenced by: '<S60>/Constant'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S25>/Gain1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S28>/Gain1'
                                        */
  0.0,                                 /* Expression: sps.K2
                                        * Referenced by: '<S30>/Gain1'
                                        */
  0.0,                                 /* Expression: sps.K1
                                        * Referenced by: '<S30>/Gain'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S41>/Gain1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S44>/Gain1'
                                        */
  1.0,                                 /* Expression: [1]
                                        * Referenced by: '<S33>/Gain'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S37>/Gain1'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval
                                        * Referenced by: '<S40>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S40>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S40>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S40>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size
   * Referenced by: '<S42>/S-Function'
   */
  { 1.0, 1.0 },
  0.025050000000000003,                /* Expression: MaxDelay
                                        * Referenced by: '<S42>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size
   * Referenced by: '<S42>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S42>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size
   * Referenced by: '<S42>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S42>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size
   * Referenced by: '<S42>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S42>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S41>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S40>/Constant'
                                        */
  1.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S40>/Unit Delay1'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S37>/Gain3'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval_h
                                        * Referenced by: '<S43>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S43>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S43>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S43>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_a
   * Referenced by: '<S45>/S-Function'
   */
  { 1.0, 1.0 },
  0.025050000000000003,                /* Expression: MaxDelay
                                        * Referenced by: '<S45>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_k
   * Referenced by: '<S45>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S45>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_h
   * Referenced by: '<S45>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S45>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_l
   * Referenced by: '<S45>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S45>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S44>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S43>/Constant'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S43>/Unit Delay1'
                                        */
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<S37>/Rad->Deg.'
                                        */
  0.0,                                 /* Expression: inf
                                        * Referenced by: '<S33>/Saturation'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S33>/Saturation'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S47>/Gain1'
                                        */

  /*  Expression: [0.92 0]
   * Referenced by: '<S21>/Constant'
   */
  { 0.92, 0.0 },

  /*  Expression: [0,0]
   * Referenced by: '<S61>/dq'
   */
  { 0.0, 0.0 },

  /*  Expression: [0,0]
   * Referenced by: '<S62>/dq'
   */
  { 0.0, 0.0 },

  /*  Computed Parameter: SFunction_P1_Size_g
   * Referenced by: '<S91>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: dest
                                        * Referenced by: '<S91>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_h
   * Referenced by: '<S91>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: priority2
                                        * Referenced by: '<S91>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_f
   * Referenced by: '<S91>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: st
                                        * Referenced by: '<S91>/S-Function'
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
   * Referenced by: '<S90>/S-Function'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: Acqu_group
                                        * Referenced by: '<S90>/S-Function'
                                        */
  8.55,                                /* Expression: 8.55
                                        * Referenced by: '<S3>/Constant1'
                                        */
  523.6,                               /* Expression: 37.4*14
                                        * Referenced by: '<S3>/Constant2'
                                        */
  0.06,                                /* Expression: 0.06
                                        * Referenced by: '<S3>/Constant3'
                                        */
  840.0,                               /* Expression: 60*14
                                        * Referenced by: '<S3>/Constant4'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_pd
   * Referenced by: '<S92>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: src
                                        * Referenced by: '<S92>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_p
   * Referenced by: '<S92>/S-Function'
   */
  { 1.0, 1.0 },
  17.0,                                /* Expression: Data_width
                                        * Referenced by: '<S92>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_m
   * Referenced by: '<S92>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: st
                                        * Referenced by: '<S92>/S-Function'
                                        */
  50.0,                                /* Expression: sps.Finit
                                        * Referenced by: '<S32>/Unit Delay'
                                        */
  70.0,                                /* Expression: 70
                                        * Referenced by: '<S21>/avoid division by zero'
                                        */
  40.0,                                /* Expression: 40
                                        * Referenced by: '<S21>/avoid division by zero'
                                        */
  0.25,                                /* Expression: 1/4
                                        * Referenced by: '<S21>/Gain'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_f
   * Referenced by: '<S58>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: MaxDelay
                                        * Referenced by: '<S58>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_j
   * Referenced by: '<S58>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S58>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_h5
   * Referenced by: '<S58>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S58>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_o
   * Referenced by: '<S58>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S58>/S-Function'
                                        */
  5.0E-5,                              /* Computed Parameter: DiscreteTimeIntegrator_gainval
                                        * Referenced by: '<S32>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: sps.Phase_Init*pi/180
                                        * Referenced by: '<S32>/Discrete-Time Integrator'
                                        */
  6.2831853071795862,                  /* Expression: 2*pi
                                        * Referenced by: '<S32>/Constant4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S21>/First cycle of simulation Id=0.92, Iq=0'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S21>/First cycle of simulation Id=0.92, Iq=0'
                                        */
  0.060141321173410534,                /* Expression: sps.D
                                        * Referenced by: '<S53>/D*u(k)'
                                        */

  /*  Expression: sps.x0(1,:)
   * Referenced by: '<S53>/Delay_x1'
   */
  { 0.00010656152949374528, 0.0 },
  7937.8816629464736,                  /* Expression: sps.C11
                                        * Referenced by: '<S56>/C11'
                                        */

  /*  Expression: sps.x0(2,:)
   * Referenced by: '<S53>/Delay_x2'
   */
  { 0.0, 0.0 },
  0.12028264234682107,                 /* Expression: sps.C12
                                        * Referenced by: '<S56>/C12'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S5>/Unit Delay2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S5>/Iq_ref'
                                        */
  5.0E-5,                              /* Computed Parameter: Integrator_gainval
                                        * Referenced by: '<S16>/Integrator'
                                        */
  0.0,                                 /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S16>/Integrator'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S17>/Gain1'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval_o
                                        * Referenced by: '<S24>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S24>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S24>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S24>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_fw
   * Referenced by: '<S26>/S-Function'
   */
  { 1.0, 1.0 },
  0.022272222222222225,                /* Expression: MaxDelay
                                        * Referenced by: '<S26>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_a
   * Referenced by: '<S26>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S26>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_l
   * Referenced by: '<S26>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S26>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_n
   * Referenced by: '<S26>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S26>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S25>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S24>/Constant'
                                        */
  1.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S24>/Unit Delay1'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S17>/Gain3'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval_f
                                        * Referenced by: '<S27>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S27>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S27>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S27>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_l
   * Referenced by: '<S29>/S-Function'
   */
  { 1.0, 1.0 },
  0.022272222222222225,                /* Expression: MaxDelay
                                        * Referenced by: '<S29>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_f
   * Referenced by: '<S29>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S29>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_fo
   * Referenced by: '<S29>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S29>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_c
   * Referenced by: '<S29>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S29>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S28>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S27>/Constant'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S27>/Unit Delay1'
                                        */
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<S17>/Rad->Deg.'
                                        */
  0.017453292519943295,                /* Expression: pi/180
                                        * Referenced by: '<S12>/to-rad'
                                        */
  0.002,                               /* Expression: RLff(1)
                                        * Referenced by: '<S9>/Rff '
                                        */
  0.2,                                 /* Expression: RLff(2)
                                        * Referenced by: '<S9>/Lff  '
                                        */
  0.002,                               /* Expression: RLff(1)
                                        * Referenced by: '<S9>/Rff'
                                        */
  0.2,                                 /* Expression: RLff(2)
                                        * Referenced by: '<S9>/Lff'
                                        */
  1.5,                                 /* Expression: 1.5
                                        * Referenced by: '<S9>/Saturation'
                                        */
  -1.5,                                /* Expression: -1.5
                                        * Referenced by: '<S9>/Saturation'
                                        */
  450.0,                               /* Expression: Limits_MPPT(1)
                                        * Referenced by: '<S10>/Iph_1'
                                        */
  375.0,                               /* Expression: Limits_MPPT(2)
                                        * Referenced by: '<S10>/Iph_2'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S3>/MPPT_On'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval_j
                                        * Referenced by: '<S30>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S30>/Integ4'
                                        */
  0.02,                                /* Expression: sps.Delay
                                        * Referenced by: '<S30>/K1'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_pj
   * Referenced by: '<S31>/S-Function'
   */
  { 1.0, 1.0 },
  0.0201,                              /* Expression: MaxDelay
                                        * Referenced by: '<S31>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_ar
   * Referenced by: '<S31>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S31>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_e
   * Referenced by: '<S31>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S31>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_c0
   * Referenced by: '<S31>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S31>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S30>/Unit Delay'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S30>/K2'
                                        */
  425.0,                               /* Expression: sps.Vinit
                                        * Referenced by: '<S30>/Unit Delay1'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval_jy
                                        * Referenced by: '<S46>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S46>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S46>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S46>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_o
   * Referenced by: '<S48>/S-Function'
   */
  { 1.0, 1.0 },
  0.025050000000000003,                /* Expression: MaxDelay
                                        * Referenced by: '<S48>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_hd
   * Referenced by: '<S48>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S48>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_j
   * Referenced by: '<S48>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S48>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_d
   * Referenced by: '<S48>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S48>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S47>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S46>/Constant'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S46>/Unit Delay1'
                                        */

  /*  Expression: [ TcD  Ts-TcD ]
   * Referenced by: '<S34>/Discrete Derivative '
   */
  { 0.0001, -5.0E-5 },
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S34>/Discrete Derivative '
                                        */
  0.16,                                /* Computed Parameter: DiscreteTimeIntegrator_gainva_j
                                        * Referenced by: '<S34>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: Par_Limits(1)
                                        * Referenced by: '<S34>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: Par_Limits(2)
                                        * Referenced by: '<S34>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: Par_Limits(1)
                                        * Referenced by: '<S34>/Saturation1'
                                        */
  0.0,                                 /* Expression: Par_Limits(2)
                                        * Referenced by: '<S34>/Saturation1'
                                        */
  0.15915494309189535,                 /* Expression: 1/2/pi
                                        * Referenced by: '<S32>/Gain10'
                                        */
  0.00060000000000000006,              /* Computed Parameter: RateLimiter_RisingLim
                                        * Referenced by: '<S32>/Rate Limiter'
                                        */
  -0.00060000000000000006,             /* Computed Parameter: RateLimiter_FallingLim
                                        * Referenced by: '<S32>/Rate Limiter'
                                        */
  50.0,                                /* Expression: sps.Finit
                                        * Referenced by: '<S32>/Rate Limiter'
                                        */
  40.528056790910988,                  /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S49>/Delay_x1'
                                        */
  0.99996932795768223,                 /* Expression: sps.A11
                                        * Referenced by: '<S50>/A11'
                                        */
  0.0,                                 /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S49>/Delay_x2'
                                        */
  4.9723130917851353E-5,               /* Expression: sps.A12
                                        * Referenced by: '<S50>/A12'
                                        */
  -1.226881692709296,                  /* Expression: sps.A21
                                        * Referenced by: '<S50>/A21'
                                        */
  0.988925236714054,                   /* Expression: sps.A22
                                        * Referenced by: '<S50>/A22'
                                        */
  2.4861565458925677E-5,               /* Expression: sps.B11
                                        * Referenced by: '<S51>/B11'
                                        */
  0.994462618357027,                   /* Expression: sps.B21
                                        * Referenced by: '<S51>/B21'
                                        */
  1.53360211588662E-5,                 /* Expression: sps.D
                                        * Referenced by: '<S49>/D*u(k)'
                                        */
  1.233694313470147,                   /* Expression: sps.C11
                                        * Referenced by: '<S52>/C11'
                                        */
  3.06720423177324E-5,                 /* Expression: sps.C12
                                        * Referenced by: '<S52>/C12'
                                        */
  0.87971735765317893,                 /* Expression: sps.A11
                                        * Referenced by: '<S54>/A11'
                                        */
  2.8483338533391971E-5,               /* Expression: sps.A12
                                        * Referenced by: '<S54>/A12'
                                        */
  -4811.3056938728414,                 /* Expression: sps.A21
                                        * Referenced by: '<S54>/A21'
                                        */
  0.13933354133567871,                 /* Expression: sps.A22
                                        * Referenced by: '<S54>/A22'
                                        */
  1.4241669266695987E-5,               /* Expression: sps.B11
                                        * Referenced by: '<S55>/B11'
                                        */
  0.56966677066783933,                 /* Expression: sps.B21
                                        * Referenced by: '<S55>/B21'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S21>/Constant1'
                                        */
  3.125E-5,                            /* Expression: sps.Delay
                                        * Referenced by: '<S85>/Constant3'
                                        */
  0.000125,                            /* Expression: sps.Period
                                        * Referenced by: '<S85>/Constant1'
                                        */
  8000.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S85>/1\ib1'
                                        */

  /*  Expression: [0 .5 1]
   * Referenced by: '<S85>/Lookup Table'
   */
  { 0.0, 0.5, 1.0 },

  /*  Expression: [0 2 0]
   * Referenced by: '<S85>/Lookup Table'
   */
  { 0.0, 2.0, 0.0 },
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S85>/Constant2'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S63>/Gain1'
                                        */
  0.1684,                              /* Expression: 0.1684
                                        * Referenced by: '<S5>/Unit Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S14>/Constant2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S5>/Unit Delay3'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S14>/Gain1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S5>/Unit Delay1'
                                        */
  2.5E-5,                              /* Computed Parameter: Integrator_gainval_i
                                        * Referenced by: '<S87>/Integrator'
                                        */
  0.9226                               /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S87>/Integrator'
                                        */
};
