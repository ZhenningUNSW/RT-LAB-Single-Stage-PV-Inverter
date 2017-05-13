/*
 * boost_and_two_level__1_sm_ehs_data.c
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

/* Block parameters (auto storage) */
P_boost_and_two_level__1_sm_ehs_T boost_and_two_level__1_sm_ehs_P = {
  2.0E-5,                              /* Variable: Ts
                                        * Referenced by:
                                        *   '<S38>/Constant4'
                                        *   '<S48>/Gain'
                                        *   '<S51>/Gain'
                                        *   '<S70>/Gain'
                                        *   '<S64>/Gain'
                                        *   '<S67>/Gain'
                                        */
  1.0,                                 /* Mask Parameter: PLL_AGC
                                        * Referenced by: '<S56>/Constant1'
                                        */
  2.0,                                 /* Mask Parameter: AlphaBetaZerotodq0_Alignment
                                        * Referenced by: '<S81>/Constant'
                                        */
  50.0,                                /* Mask Parameter: InverterControl_Fnom
                                        * Referenced by:
                                        *   '<S38>/Constant4'
                                        *   '<S45>/First cycle of simulation Id=0.92, Iq=0'
                                        */
  0.003,                               /* Mask Parameter: InverterControl_Increment_MPPT
                                        * Referenced by: '<S34>/Iph_3'
                                        */
  314.15926535897933,                  /* Mask Parameter: Discrete_Init
                                        * Referenced by: '<S58>/Discrete-Time Integrator'
                                        */
  1.0,                                 /* Mask Parameter: Discrete_Kd
                                        * Referenced by: '<S58>/Discrete Derivative '
                                        */
  6.6,                                 /* Mask Parameter: InverterControl_Ki_Ireg
                                        * Referenced by: '<S40>/Integral Gain'
                                        */
  200.0,                               /* Mask Parameter: InverterControl_Ki_VDCreg
                                        * Referenced by: '<S111>/Integral Gain'
                                        */
  180.0,                               /* Mask Parameter: Discrete_Kp
                                        * Referenced by: '<S58>/Kp4'
                                        */
  0.15,                                /* Mask Parameter: InverterControl_Kp_Ireg
                                        * Referenced by: '<S40>/Proportional Gain'
                                        */
  12.0,                                /* Mask Parameter: InverterControl_Kp_VDCreg
                                        * Referenced by: '<S111>/Proportional Gain'
                                        */
  -1.5,                                /* Mask Parameter: PI_LowerSaturationLimit
                                        * Referenced by: '<S40>/Saturate'
                                        */
  -1.5,                                /* Mask Parameter: PI_LowerSaturationLimit_j
                                        * Referenced by: '<S111>/Saturate'
                                        */

  /*  Mask Parameter: PWM_Generator_MinMax
   * Referenced by: '<S37>/Constant10'
   */
  { -1.0, 1.0 },
  3500.0,                              /* Mask Parameter: InverterControl_Pnom
                                        * Referenced by: '<S36>/A->pu'
                                        */
  1.5,                                 /* Mask Parameter: PI_UpperSaturationLimit
                                        * Referenced by: '<S40>/Saturate'
                                        */
  1.5,                                 /* Mask Parameter: PI_UpperSaturationLimit_a
                                        * Referenced by: '<S111>/Saturate'
                                        */
  425.0,                               /* Mask Parameter: InverterControl_Vdc_ref_Init
                                        * Referenced by:
                                        *   '<S31>/Vnom_dc1'
                                        *   '<S34>/Iph_'
                                        */
  400.0,                               /* Mask Parameter: InverterControl_Vnom_dc
                                        * Referenced by: '<S39>/Rtot_pu2'
                                        */
  240.0,                               /* Mask Parameter: InverterControl_Vnom_prim
                                        * Referenced by:
                                        *   '<S36>/A->pu'
                                        *   '<S36>/V->pu'
                                        *   '<S38>/Constant3'
                                        */
  1.0,                                 /* Mask Parameter: CompareToConstant_const
                                        * Referenced by: '<S83>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: CompareToConstant1_const
                                        * Referenced by: '<S84>/Constant'
                                        */
  0.2,                                 /* Mask Parameter: CompareToConstant_const_d
                                        * Referenced by: '<S3>/Constant'
                                        */
  5.0,                                 /* Mask Parameter: Counter_max_count
                                        * Referenced by: '<S19>/Switch'
                                        */
  127U,                                /* Mask Parameter: BitwiseOperator_BitMask
                                        * Referenced by: '<S22>/Bitwise Operator'
                                        */
  65535U,                              /* Mask Parameter: BitwiseOperator_BitMask_j
                                        * Referenced by: '<S20>/Bitwise Operator'
                                        */
  65535U,                              /* Mask Parameter: BitwiseOperator_BitMask_h
                                        * Referenced by: '<S21>/Bitwise Operator'
                                        */
  3U,                                  /* Mask Parameter: DetectChange_vinit
                                        * Referenced by: '<S25>/Delay Input1'
                                        */
  0U,                                  /* Mask Parameter: DetectChange_vinit_i
                                        * Referenced by: '<S24>/Delay Input1'
                                        */
  0U,                                  /* Mask Parameter: DetectChange_vinit_g
                                        * Referenced by: '<S28>/Delay Input1'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S15>/Gain'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S16>/Gain'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S19>/Constant1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S19>/Constant2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S19>/Constant'
                                        */
  31.99951171875,                      /* Expression: 2^(N-Q)-2^(-Q)
                                        * Referenced by: '<S20>/Saturation1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S20>/Saturation1'
                                        */
  2048.0,                              /* Expression: 2^Q
                                        * Referenced by: '<S20>/Gain2'
                                        */
  15.99951171875,                      /* Expression: 2^(N-Q-1)-2^(-Q)
                                        * Referenced by: '<S20>/Saturation'
                                        */
  -16.0,                               /* Expression: -2^(N-Q-1)
                                        * Referenced by: '<S20>/Saturation'
                                        */
  2048.0,                              /* Expression: 2^Q
                                        * Referenced by: '<S20>/Gain1'
                                        */
  31.99951171875,                      /* Expression: 2^(N-Q)-2^(-Q)
                                        * Referenced by: '<S21>/Saturation1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S21>/Saturation1'
                                        */
  2048.0,                              /* Expression: 2^Q
                                        * Referenced by: '<S21>/Gain2'
                                        */
  15.99951171875,                      /* Expression: 2^(N-Q-1)-2^(-Q)
                                        * Referenced by: '<S21>/Saturation'
                                        */
  -16.0,                               /* Expression: -2^(N-Q-1)
                                        * Referenced by: '<S21>/Saturation'
                                        */
  2048.0,                              /* Expression: 2^Q
                                        * Referenced by: '<S21>/Gain1'
                                        */
  127.0,                               /* Expression: 2^(N-Q)-2^(-Q)
                                        * Referenced by: '<S22>/Saturation1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S22>/Saturation1'
                                        */
  1.0,                                 /* Expression: 2^Q
                                        * Referenced by: '<S22>/Gain2'
                                        */
  63.0,                                /* Expression: 2^(N-Q-1)-2^(-Q)
                                        * Referenced by: '<S22>/Saturation'
                                        */
  -64.0,                               /* Expression: -2^(N-Q-1)
                                        * Referenced by: '<S22>/Saturation'
                                        */
  1.0,                                 /* Expression: 2^Q
                                        * Referenced by: '<S22>/Gain1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S49>/Gain1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S52>/Gain1'
                                        */
  -1.1368683772160957E-16,             /* Expression: sps.K2
                                        * Referenced by: '<S54>/Gain1'
                                        */
  6.4623485355705288E-30,              /* Expression: sps.K1
                                        * Referenced by: '<S54>/Gain'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S65>/Gain1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S68>/Gain1'
                                        */
  1.0,                                 /* Expression: [1]
                                        * Referenced by: '<S57>/Gain'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S61>/Gain1'
                                        */
  1.0E-5,                              /* Computed Parameter: Integ4_gainval
                                        * Referenced by: '<S64>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S64>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S64>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S64>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size
   * Referenced by: '<S66>/S-Function'
   */
  { 1.0, 1.0 },
  0.02502,                             /* Expression: MaxDelay
                                        * Referenced by: '<S66>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size
   * Referenced by: '<S66>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S66>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size
   * Referenced by: '<S66>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S66>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size
   * Referenced by: '<S66>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S66>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S65>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S64>/Constant'
                                        */
  1.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S64>/Unit Delay1'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S61>/Gain3'
                                        */
  1.0E-5,                              /* Computed Parameter: Integ4_gainval_c
                                        * Referenced by: '<S67>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S67>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S67>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S67>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_h
   * Referenced by: '<S69>/S-Function'
   */
  { 1.0, 1.0 },
  0.02502,                             /* Expression: MaxDelay
                                        * Referenced by: '<S69>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_p
   * Referenced by: '<S69>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S69>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_p
   * Referenced by: '<S69>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S69>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_p
   * Referenced by: '<S69>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S69>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S68>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S67>/Constant'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S67>/Unit Delay1'
                                        */
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<S61>/Rad->Deg.'
                                        */
  0.0,                                 /* Expression: inf
                                        * Referenced by: '<S57>/Saturation'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S57>/Saturation'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S71>/Gain1'
                                        */

  /*  Expression: [0.92 0]
   * Referenced by: '<S45>/Constant'
   */
  { 0.92, 0.0 },

  /*  Expression: [0,0]
   * Referenced by: '<S85>/dq'
   */
  { 0.0, 0.0 },

  /*  Expression: [0,0]
   * Referenced by: '<S86>/dq'
   */
  { 0.0, 0.0 },
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S1>/S-Function1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S1>/S-Function'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P1_Size
   * Referenced by: '<S7>/Outputs_eHS1_Recv'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: Outputs_eHS1_Recv_P1
   * Referenced by: '<S7>/Outputs_eHS1_Recv'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: Outputs_eHS1_Recv_P2_Size
   * Referenced by: '<S7>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: FcnNos
                                        * Referenced by: '<S7>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P3_Size
   * Referenced by: '<S7>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  8.0,                                 /* Expression: width
                                        * Referenced by: '<S7>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P4_Size
   * Referenced by: '<S7>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S7>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P5_Size
   * Referenced by: '<S7>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: sampleTime
                                        * Referenced by: '<S7>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P6_Size
   * Referenced by: '<S7>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: checkVersion
                                        * Referenced by: '<S7>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P7_Size
   * Referenced by: '<S7>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedId
                                        * Referenced by: '<S7>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P8_Size
   * Referenced by: '<S7>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedVersion
                                        * Referenced by: '<S7>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P9_Size
   * Referenced by: '<S7>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: opComp
                                        * Referenced by: '<S7>/Outputs_eHS1_Recv'
                                        */
  0.1684,                              /* Expression: 0.1684
                                        * Referenced by: '<S31>/Unit Delay'
                                        */
  3.125E-5,                            /* Expression: sps.Delay
                                        * Referenced by: '<S109>/Constant3'
                                        */
  0.000125,                            /* Expression: sps.Period
                                        * Referenced by: '<S109>/Constant1'
                                        */
  8000.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S109>/1\ib1'
                                        */

  /*  Expression: [0 .5 1]
   * Referenced by: '<S109>/Lookup Table'
   */
  { 0.0, 0.5, 1.0 },

  /*  Expression: [0 2 0]
   * Referenced by: '<S109>/Lookup Table'
   */
  { 0.0, 2.0, 0.0 },
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S109>/Constant2'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S87>/Gain1'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_f
   * Referenced by: '<S10>/S-Function'
   */
  { 1.0, 1.0 },
  164.0,                               /* Expression: Data_width
                                        * Referenced by: '<S10>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_k
   * Referenced by: '<S10>/S-Function'
   */
  { 164.0, 1.0 },

  /*  Expression: InitialConditions
   * Referenced by: '<S10>/S-Function'
   */
  { 1000.0, 25.0, 1.0, 0.0, 48.0, 49.0, 52.0, 53.0, 53.0, 9.0, 1.0, 3.0, 12.0,
    13.0, 0.0, 2.0, 48.0, 49.0, 50.0, 51.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0,
    54.0, 54.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 0.0077647, 0.33,
    0.08, 0.0048529, 0.0048529, 1.0, 0.005, 0.005, 0.1, 0.1, 0.1, 0.1, 0.005,
    0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.65, 1.65, 1.65, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.3, 3.3, 3.3, 3.3, 3.3, 5.0, 5.0, 5.0, 5.0,
    5.0, 5.0, 5.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0,
    16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -16.0, -16.0, -16.0,
    -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0,
    -16.0, -16.0, -16.0, -16.0, -16.0, -16.0 },
  0.5,                                 /* Expression: .5
                                        * Referenced by: '<S7>/Relay1'
                                        */
  0.5,                                 /* Expression: .5
                                        * Referenced by: '<S7>/Relay1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S7>/Relay1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S7>/Relay1'
                                        */

  /*  Computed Parameter: eHS_rst_loadin_P1_Size
   * Referenced by: '<S7>/eHS_rst_loadin'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: eHS_rst_loadin_P1
   * Referenced by: '<S7>/eHS_rst_loadin'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: eHS_rst_loadin_P2_Size
   * Referenced by: '<S7>/eHS_rst_loadin'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: FcnNos
                                        * Referenced by: '<S7>/eHS_rst_loadin'
                                        */

  /*  Computed Parameter: eHS_rst_loadin_P3_Size
   * Referenced by: '<S7>/eHS_rst_loadin'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: width
                                        * Referenced by: '<S7>/eHS_rst_loadin'
                                        */

  /*  Computed Parameter: eHS_rst_loadin_P4_Size
   * Referenced by: '<S7>/eHS_rst_loadin'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S7>/eHS_rst_loadin'
                                        */

  /*  Computed Parameter: Automated_Solver_Mat_Initialisa
   * Referenced by: '<S7>/Automated_Solver_Mat_Initialisation_1'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: Automated_Solver_Mat_Initiali_c
   * Referenced by: '<S7>/Automated_Solver_Mat_Initialisation_1'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: Automated_Solver_Mat_Initiali_o
   * Referenced by: '<S7>/Automated_Solver_Mat_Initialisation_1'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: fpga_port_in
                                        * Referenced by: '<S7>/Automated_Solver_Mat_Initialisation_1'
                                        */

  /*  Computed Parameter: Automated_Solver_Mat_Initiali_e
   * Referenced by: '<S7>/Automated_Solver_Mat_Initialisation_1'
   */
  { 1.0, 31.0 },

  /*  Computed Parameter: Automated_Solver_Mat_Initial_er
   * Referenced by: '<S7>/Automated_Solver_Mat_Initialisation_1'
   */
  { 101.0, 104.0, 115.0, 95.0, 99.0, 111.0, 110.0, 102.0, 105.0, 103.0, 95.0,
    101.0, 72.0, 83.0, 120.0, 54.0, 52.0, 118.0, 51.0, 95.0, 83.0, 73.0, 68.0,
    52.0, 56.0, 57.0, 52.0, 46.0, 109.0, 97.0, 116.0 },
  20364.675298172569,                  /* Expression: 14400*sqrt(2)
                                        * Referenced by: '<S9>/Sine Wave Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S9>/Sine Wave Function'
                                        */
  314.15926535897933,                  /* Expression: 50*2*pi
                                        * Referenced by: '<S9>/Sine Wave Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S9>/Sine Wave Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S9>/Memory1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S9>/Memory'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P1_Size
   * Referenced by: '<S7>/Inputs_eHS1_Send'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: Inputs_eHS1_Send_P1
   * Referenced by: '<S7>/Inputs_eHS1_Send'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: Inputs_eHS1_Send_P2_Size
   * Referenced by: '<S7>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: FcnNos
                                        * Referenced by: '<S7>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P3_Size
   * Referenced by: '<S7>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  3.0,                                 /* Expression: width
                                        * Referenced by: '<S7>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P4_Size
   * Referenced by: '<S7>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S7>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P5_Size
   * Referenced by: '<S7>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: checkVersion
                                        * Referenced by: '<S7>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P6_Size
   * Referenced by: '<S7>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedId
                                        * Referenced by: '<S7>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P7_Size
   * Referenced by: '<S7>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedVersion
                                        * Referenced by: '<S7>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P8_Size
   * Referenced by: '<S7>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: opComp
                                        * Referenced by: '<S7>/Inputs_eHS1_Send'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S15>/integrator'
                                        */
  0.02,                                /* Expression: 1./sps.Freq
                                        * Referenced by: '<S15>/Transport Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S15>/Transport Delay'
                                        */
  0.02,                                /* Expression: 1./sps.Freq
                                        * Referenced by: '<S15>/K1'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S15>/Memory'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S16>/integrator'
                                        */
  0.02,                                /* Expression: 1./sps.Freq
                                        * Referenced by: '<S16>/Transport Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S16>/Transport Delay'
                                        */
  0.02,                                /* Expression: 1./sps.Freq
                                        * Referenced by: '<S16>/K1'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S16>/Memory'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P1_Size_e
   * Referenced by: '<S2>/Outputs_eHS1_Recv'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: Outputs_eHS1_Recv_P1_p
   * Referenced by: '<S2>/Outputs_eHS1_Recv'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: Outputs_eHS1_Recv_P2_Size_g
   * Referenced by: '<S2>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: FcnNos
                                        * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P3_Size_a
   * Referenced by: '<S2>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  7.0,                                 /* Expression: width
                                        * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P4_Size_b
   * Referenced by: '<S2>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P5_Size_m
   * Referenced by: '<S2>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: sampleTime
                                        * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P6_Size_b
   * Referenced by: '<S2>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: checkVersion
                                        * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P7_Size_j
   * Referenced by: '<S2>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedId
                                        * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P8_Size_p
   * Referenced by: '<S2>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedVersion
                                        * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P9_Size_m
   * Referenced by: '<S2>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: opComp
                                        * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: OpMonitor_P1_Size
   * Referenced by: '<S2>/OpMonitor'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: compute_time
                                        * Referenced by: '<S2>/OpMonitor'
                                        */

  /*  Computed Parameter: OpMonitor_P2_Size
   * Referenced by: '<S2>/OpMonitor'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: real_step
                                        * Referenced by: '<S2>/OpMonitor'
                                        */

  /*  Computed Parameter: OpMonitor_P3_Size
   * Referenced by: '<S2>/OpMonitor'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: idle_time
                                        * Referenced by: '<S2>/OpMonitor'
                                        */

  /*  Computed Parameter: OpMonitor_P4_Size
   * Referenced by: '<S2>/OpMonitor'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: nb_overruns
                                        * Referenced by: '<S2>/OpMonitor'
                                        */

  /*  Computed Parameter: OpMonitor_P5_Size
   * Referenced by: '<S2>/OpMonitor'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: user_time
                                        * Referenced by: '<S2>/OpMonitor'
                                        */

  /*  Computed Parameter: OpMonitor_P6_Size
   * Referenced by: '<S2>/OpMonitor'
   */
  { 1.0, 32.0 },

  /*  Computed Parameter: OpMonitor_P6
   * Referenced by: '<S2>/OpMonitor'
   */
  { 109.0, 121.0, 95.0, 101.0, 118.0, 101.0, 110.0, 116.0, 95.0, 110.0, 97.0,
    109.0, 101.0, 44.0, 97.0, 110.0, 111.0, 116.0, 104.0, 101.0, 114.0, 95.0,
    101.0, 118.0, 101.0, 110.0, 116.0, 95.0, 110.0, 97.0, 109.0, 101.0 },

  /*  Computed Parameter: SFunction_P1_Size_l
   * Referenced by: '<S30>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: Acqu_group
                                        * Referenced by: '<S30>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S19>/Memory1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S19>/Memory2'
                                        */

  /*  Computed Parameter: LoadIn_P1_Size
   * Referenced by: '<S17>/LoadIn'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: LoadIn_P1
   * Referenced by: '<S17>/LoadIn'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: LoadIn_P2_Size
   * Referenced by: '<S17>/LoadIn'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: FcnNos
                                        * Referenced by: '<S17>/LoadIn'
                                        */

  /*  Computed Parameter: LoadIn_P3_Size
   * Referenced by: '<S17>/LoadIn'
   */
  { 1.0, 1.0 },
  18.0,                                /* Expression: width
                                        * Referenced by: '<S17>/LoadIn'
                                        */

  /*  Computed Parameter: LoadIn_P4_Size
   * Referenced by: '<S17>/LoadIn'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S17>/LoadIn'
                                        */

  /*  Expression: ones(1,9)
   * Referenced by: '<S18>/Constant'
   */
  { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 },

  /*  Computed Parameter: DataInSend_P1_Size
   * Referenced by: '<S17>/DataIn Send'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: DataInSend_P1
   * Referenced by: '<S17>/DataIn Send'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: DataInSend_P2_Size
   * Referenced by: '<S17>/DataIn Send'
   */
  { 1.0, 1.0 },
  15.0,                                /* Expression: FcnNos
                                        * Referenced by: '<S17>/DataIn Send'
                                        */

  /*  Computed Parameter: DataInSend_P3_Size
   * Referenced by: '<S17>/DataIn Send'
   */
  { 1.0, 1.0 },
  16.0,                                /* Expression: width
                                        * Referenced by: '<S17>/DataIn Send'
                                        */

  /*  Computed Parameter: DataInSend_P4_Size
   * Referenced by: '<S17>/DataIn Send'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S17>/DataIn Send'
                                        */

  /*  Computed Parameter: DataInSend_P5_Size
   * Referenced by: '<S17>/DataIn Send'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: checkVersion
                                        * Referenced by: '<S17>/DataIn Send'
                                        */

  /*  Computed Parameter: DataInSend_P6_Size
   * Referenced by: '<S17>/DataIn Send'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedId
                                        * Referenced by: '<S17>/DataIn Send'
                                        */

  /*  Computed Parameter: DataInSend_P7_Size
   * Referenced by: '<S17>/DataIn Send'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedVersion
                                        * Referenced by: '<S17>/DataIn Send'
                                        */

  /*  Computed Parameter: DataInSend_P8_Size
   * Referenced by: '<S17>/DataIn Send'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: opComp
                                        * Referenced by: '<S17>/DataIn Send'
                                        */

  /*  Computed Parameter: OpWriteFile_P1_Size
   * Referenced by: '<S2>/OpWriteFile'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: Acq_Group
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */

  /*  Computed Parameter: OpWriteFile_P2_Size
   * Referenced by: '<S2>/OpWriteFile'
   */
  { 1.0, 17.0 },

  /*  Computed Parameter: OpWriteFile_P2
   * Referenced by: '<S2>/OpWriteFile'
   */
  { 109.0, 121.0, 80.0, 86.0, 95.0, 105.0, 110.0, 118.0, 101.0, 114.0, 116.0,
    101.0, 114.0, 46.0, 109.0, 97.0, 116.0 },

  /*  Computed Parameter: OpWriteFile_P3_Size
   * Referenced by: '<S2>/OpWriteFile'
   */
  { 1.0, 5.0 },

  /*  Computed Parameter: OpWriteFile_P3
   * Referenced by: '<S2>/OpWriteFile'
   */
  { 111.0, 112.0, 118.0, 97.0, 114.0 },

  /*  Computed Parameter: OpWriteFile_P4_Size
   * Referenced by: '<S2>/OpWriteFile'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: Decimation
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */

  /*  Computed Parameter: OpWriteFile_P5_Size
   * Referenced by: '<S2>/OpWriteFile'
   */
  { 1.0, 1.0 },
  400000.0,                            /* Expression: Nb_Samples
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */

  /*  Computed Parameter: OpWriteFile_P6_Size
   * Referenced by: '<S2>/OpWriteFile'
   */
  { 1.0, 1.0 },
  1.0E+6,                              /* Expression: Buffer_size
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */

  /*  Computed Parameter: OpWriteFile_P7_Size
   * Referenced by: '<S2>/OpWriteFile'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: Sim_Mode
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */

  /*  Computed Parameter: OpWriteFile_P8_Size
   * Referenced by: '<S2>/OpWriteFile'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: Static_File
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */

  /*  Computed Parameter: OpWriteFile_P9_Size
   * Referenced by: '<S2>/OpWriteFile'
   */
  { 1.0, 1.0 },
  1.4999999999999998E+9,               /* Expression: file_size
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */

  /*  Computed Parameter: OpWriteFile_P10_Size
   * Referenced by: '<S2>/OpWriteFile'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: write_offline
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */

  /*  Computed Parameter: OpWriteFile_P11_Size
   * Referenced by: '<S2>/OpWriteFile'
   */
  { 1.0, 2.0 },

  /*  Computed Parameter: OpWriteFile_P11
   * Referenced by: '<S2>/OpWriteFile'
   */
  { 46.0, 47.0 },

  /*  Computed Parameter: OpWriteFile_P12_Size
   * Referenced by: '<S2>/OpWriteFile'
   */
  { 1.0, 2.0 },

  /*  Computed Parameter: OpWriteFile_P12
   * Referenced by: '<S2>/OpWriteFile'
   */
  { 46.0, 92.0 },

  /*  Computed Parameter: RTEConversion1_P1_Size
   * Referenced by: '<S9>/RTE Conversion1'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: nbMaxEvents
                                        * Referenced by: '<S9>/RTE Conversion1'
                                        */

  /*  Computed Parameter: RTEConversion1_P2_Size
   * Referenced by: '<S9>/RTE Conversion1'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: inputdatatype
                                        * Referenced by: '<S9>/RTE Conversion1'
                                        */

  /*  Computed Parameter: RTEConversion1_P3_Size
   * Referenced by: '<S9>/RTE Conversion1'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: outputdatatype
                                        * Referenced by: '<S9>/RTE Conversion1'
                                        */

  /*  Computed Parameter: RTEConversion1_P4_Size
   * Referenced by: '<S9>/RTE Conversion1'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: compensation
                                        * Referenced by: '<S9>/RTE Conversion1'
                                        */

  /*  Computed Parameter: RTEConversion1_P5_Size
   * Referenced by: '<S9>/RTE Conversion1'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: sampleTime
                                        * Referenced by: '<S9>/RTE Conversion1'
                                        */
  0.01,                                /* Expression: 0.01
                                        * Referenced by: '<S9>/Step'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S9>/Step'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S9>/Step'
                                        */

  /*  Computed Parameter: RTEConversion_P1_Size
   * Referenced by: '<S9>/RTE Conversion'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: nbMaxEvents
                                        * Referenced by: '<S9>/RTE Conversion'
                                        */

  /*  Computed Parameter: RTEConversion_P2_Size
   * Referenced by: '<S9>/RTE Conversion'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: inputdatatype
                                        * Referenced by: '<S9>/RTE Conversion'
                                        */

  /*  Computed Parameter: RTEConversion_P3_Size
   * Referenced by: '<S9>/RTE Conversion'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: outputdatatype
                                        * Referenced by: '<S9>/RTE Conversion'
                                        */

  /*  Computed Parameter: RTEConversion_P4_Size
   * Referenced by: '<S9>/RTE Conversion'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: compensation
                                        * Referenced by: '<S9>/RTE Conversion'
                                        */

  /*  Computed Parameter: RTEConversion_P5_Size
   * Referenced by: '<S9>/RTE Conversion'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: sampleTime
                                        * Referenced by: '<S9>/RTE Conversion'
                                        */

  /*  Computed Parameter: RTELogicalOperator1_P1_Size
   * Referenced by: '<S9>/RTE Logical Operator1'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: LogicalOperator
                                        * Referenced by: '<S9>/RTE Logical Operator1'
                                        */

  /*  Computed Parameter: RTELogicalOperator1_P2_Size
   * Referenced by: '<S9>/RTE Logical Operator1'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: NbrInput
                                        * Referenced by: '<S9>/RTE Logical Operator1'
                                        */

  /*  Computed Parameter: RTELogicalOperator1_P3_Size
   * Referenced by: '<S9>/RTE Logical Operator1'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: NbrMaxEvents
                                        * Referenced by: '<S9>/RTE Logical Operator1'
                                        */

  /*  Computed Parameter: RTEGround_P1_Size
   * Referenced by: '<S9>/RTE Ground'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: SampleTime
                                        * Referenced by: '<S9>/RTE Ground'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P1_Size
   * Referenced by: '<S26>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: conversionType
                                        * Referenced by: '<S26>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P2_Size
   * Referenced by: '<S26>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  8.0,                                 /* Expression: nbLine
                                        * Referenced by: '<S26>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P3_Size
   * Referenced by: '<S26>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: nbEvents
                                        * Referenced by: '<S26>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P4_Size
   * Referenced by: '<S26>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: timeUnit
                                        * Referenced by: '<S26>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P5_Size
   * Referenced by: '<S26>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: enTimeFactor
                                        * Referenced by: '<S26>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P6_Size
   * Referenced by: '<S26>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: sampleTime
                                        * Referenced by: '<S26>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P7_Size
   * Referenced by: '<S26>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: initStates
                                        * Referenced by: '<S26>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: EventGen_eHS_1_P1_Size
   * Referenced by: '<S26>/EventGen_eHS_1'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: EventGen_eHS_1_P1
   * Referenced by: '<S26>/EventGen_eHS_1'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: EventGen_eHS_1_P2_Size
   * Referenced by: '<S26>/EventGen_eHS_1'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portNb
                                        * Referenced by: '<S26>/EventGen_eHS_1'
                                        */

  /*  Computed Parameter: EventGen_eHS_1_P3_Size
   * Referenced by: '<S26>/EventGen_eHS_1'
   */
  { 0.0, 0.0 },

  /*  Computed Parameter: EventGen_eHS_1_P4_Size
   * Referenced by: '<S26>/EventGen_eHS_1'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: size
                                        * Referenced by: '<S26>/EventGen_eHS_1'
                                        */

  /*  Computed Parameter: EventGen_eHS_1_P5_Size
   * Referenced by: '<S26>/EventGen_eHS_1'
   */
  { 1.0, 1.0 },
  8.0,                                 /* Expression: nbChannels
                                        * Referenced by: '<S26>/EventGen_eHS_1'
                                        */

  /*  Computed Parameter: EventGen_eHS_1_P6_Size
   * Referenced by: '<S26>/EventGen_eHS_1'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: numwidth
                                        * Referenced by: '<S26>/EventGen_eHS_1'
                                        */

  /*  Computed Parameter: EventGen_eHS_1_P7_Size
   * Referenced by: '<S26>/EventGen_eHS_1'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: timeUnit
                                        * Referenced by: '<S26>/EventGen_eHS_1'
                                        */
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
  50.0,                                /* Expression: sps.Finit
                                        * Referenced by: '<S56>/Unit Delay'
                                        */
  70.0,                                /* Expression: 70
                                        * Referenced by: '<S45>/avoid division by zero'
                                        */
  40.0,                                /* Expression: 40
                                        * Referenced by: '<S45>/avoid division by zero'
                                        */
  0.25,                                /* Expression: 1/4
                                        * Referenced by: '<S45>/Gain'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_i
   * Referenced by: '<S82>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: MaxDelay
                                        * Referenced by: '<S82>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_i
   * Referenced by: '<S82>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S82>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_i
   * Referenced by: '<S82>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S82>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_n
   * Referenced by: '<S82>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S82>/S-Function'
                                        */
  2.0E-5,                              /* Computed Parameter: DiscreteTimeIntegrator_gainval
                                        * Referenced by: '<S56>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: sps.Phase_Init*pi/180
                                        * Referenced by: '<S56>/Discrete-Time Integrator'
                                        */
  6.2831853071795862,                  /* Expression: 2*pi
                                        * Referenced by: '<S56>/Constant4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S45>/First cycle of simulation Id=0.92, Iq=0'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S45>/First cycle of simulation Id=0.92, Iq=0'
                                        */
  0.01257992044073508,                 /* Expression: sps.D
                                        * Referenced by: '<S77>/D*u(k)'
                                        */

  /*  Expression: sps.x0(1,:)
   * Referenced by: '<S77>/Delay_x1'
   */
  { 0.00028197057830512594, 0.0 },
  3151.6695002188571,                  /* Expression: sps.C11
                                        * Referenced by: '<S80>/C11'
                                        */

  /*  Expression: sps.x0(2,:)
   * Referenced by: '<S77>/Delay_x2'
   */
  { -5.5511151231257827E-17, 0.0 },
  0.02515984088147016,                 /* Expression: sps.C12
                                        * Referenced by: '<S80>/C12'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S31>/Unit Delay2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S31>/Iq_ref'
                                        */
  2.0E-5,                              /* Computed Parameter: Integrator_gainval
                                        * Referenced by: '<S40>/Integrator'
                                        */
  0.0,                                 /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S40>/Integrator'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S41>/Gain1'
                                        */
  1.0E-5,                              /* Computed Parameter: Integ4_gainval_l
                                        * Referenced by: '<S48>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S48>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S48>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S48>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_ll
   * Referenced by: '<S50>/S-Function'
   */
  { 1.0, 1.0 },
  0.022242222222222222,                /* Expression: MaxDelay
                                        * Referenced by: '<S50>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_c
   * Referenced by: '<S50>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S50>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_b
   * Referenced by: '<S50>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S50>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_f
   * Referenced by: '<S50>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S50>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S49>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S48>/Constant'
                                        */
  1.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S48>/Unit Delay1'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S41>/Gain3'
                                        */
  1.0E-5,                              /* Computed Parameter: Integ4_gainval_g
                                        * Referenced by: '<S51>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S51>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S51>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S51>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_hi
   * Referenced by: '<S53>/S-Function'
   */
  { 1.0, 1.0 },
  0.022242222222222222,                /* Expression: MaxDelay
                                        * Referenced by: '<S53>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_j
   * Referenced by: '<S53>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S53>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_c
   * Referenced by: '<S53>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S53>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_m
   * Referenced by: '<S53>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S53>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S52>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S51>/Constant'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S51>/Unit Delay1'
                                        */
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<S41>/Rad->Deg.'
                                        */
  0.017453292519943295,                /* Expression: pi/180
                                        * Referenced by: '<S36>/to-rad'
                                        */
  0.002,                               /* Expression: RLff(1)
                                        * Referenced by: '<S33>/Rff '
                                        */
  0.2,                                 /* Expression: RLff(2)
                                        * Referenced by: '<S33>/Lff  '
                                        */
  0.002,                               /* Expression: RLff(1)
                                        * Referenced by: '<S33>/Rff'
                                        */
  0.2,                                 /* Expression: RLff(2)
                                        * Referenced by: '<S33>/Lff'
                                        */
  1.5,                                 /* Expression: 1.5
                                        * Referenced by: '<S33>/Saturation'
                                        */
  -1.5,                                /* Expression: -1.5
                                        * Referenced by: '<S33>/Saturation'
                                        */
  450.0,                               /* Expression: Limits_MPPT(1)
                                        * Referenced by: '<S34>/Iph_1'
                                        */
  375.0,                               /* Expression: Limits_MPPT(2)
                                        * Referenced by: '<S34>/Iph_2'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S9>/MPPT_On'
                                        */
  1.0E-5,                              /* Computed Parameter: Integ4_gainval_o
                                        * Referenced by: '<S54>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S54>/Integ4'
                                        */
  0.02,                                /* Expression: sps.Delay
                                        * Referenced by: '<S54>/K1'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_o
   * Referenced by: '<S55>/S-Function'
   */
  { 1.0, 1.0 },
  0.02004,                             /* Expression: MaxDelay
                                        * Referenced by: '<S55>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_id
   * Referenced by: '<S55>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S55>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_l
   * Referenced by: '<S55>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S55>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_i
   * Referenced by: '<S55>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S55>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S54>/Unit Delay'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S54>/K2'
                                        */
  425.0,                               /* Expression: sps.Vinit
                                        * Referenced by: '<S54>/Unit Delay1'
                                        */
  1.0E-5,                              /* Computed Parameter: Integ4_gainval_b
                                        * Referenced by: '<S70>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S70>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S70>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S70>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_lc
   * Referenced by: '<S72>/S-Function'
   */
  { 1.0, 1.0 },
  0.02502,                             /* Expression: MaxDelay
                                        * Referenced by: '<S72>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_g
   * Referenced by: '<S72>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S72>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_f
   * Referenced by: '<S72>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S72>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_a
   * Referenced by: '<S72>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S72>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S71>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S70>/Constant'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S70>/Unit Delay1'
                                        */

  /*  Expression: [ TcD  Ts-TcD ]
   * Referenced by: '<S58>/Discrete Derivative '
   */
  { 0.0001, -8.0E-5 },
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S58>/Discrete Derivative '
                                        */
  0.064,                               /* Computed Parameter: DiscreteTimeIntegrator_gainva_a
                                        * Referenced by: '<S58>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: Par_Limits(1)
                                        * Referenced by: '<S58>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: Par_Limits(2)
                                        * Referenced by: '<S58>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: Par_Limits(1)
                                        * Referenced by: '<S58>/Saturation1'
                                        */
  0.0,                                 /* Expression: Par_Limits(2)
                                        * Referenced by: '<S58>/Saturation1'
                                        */
  0.15915494309189535,                 /* Expression: 1/2/pi
                                        * Referenced by: '<S56>/Gain10'
                                        */
  0.00024000000000000003,              /* Computed Parameter: RateLimiter_RisingLim
                                        * Referenced by: '<S56>/Rate Limiter'
                                        */
  -0.00024000000000000003,             /* Computed Parameter: RateLimiter_FallingLim
                                        * Referenced by: '<S56>/Rate Limiter'
                                        */
  50.0,                                /* Expression: sps.Finit
                                        * Referenced by: '<S56>/Rate Limiter'
                                        */
  101.32101697571224,                  /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S73>/Delay_x1'
                                        */
  0.999995076138259,                   /* Expression: sps.A11
                                        * Referenced by: '<S74>/A11'
                                        */
  0.0,                                 /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S73>/Delay_x2'
                                        */
  1.995562716158404E-5,                /* Expression: sps.A12
                                        * Referenced by: '<S74>/A12'
                                        */
  -0.49238617409376229,                /* Expression: sps.A21
                                        * Referenced by: '<S74>/A21'
                                        */
  0.99556271615840386,                 /* Expression: sps.A22
                                        * Referenced by: '<S74>/A22'
                                        */
  9.97781358079202E-6,                 /* Expression: sps.B11
                                        * Referenced by: '<S75>/B11'
                                        */
  0.99778135807920187,                 /* Expression: sps.B21
                                        * Referenced by: '<S75>/B21'
                                        */
  2.4619308704688116E-6,               /* Expression: sps.D
                                        * Referenced by: '<S73>/D*u(k)'
                                        */
  0.49347981688184189,                 /* Expression: sps.C11
                                        * Referenced by: '<S76>/C11'
                                        */
  4.9238617409376232E-6,               /* Expression: sps.C12
                                        * Referenced by: '<S76>/C12'
                                        */
  0.97484015911852984,                 /* Expression: sps.A11
                                        * Referenced by: '<S78>/A11'
                                        */
  1.576518862980693E-5,                /* Expression: sps.A12
                                        * Referenced by: '<S78>/A12'
                                        */
  -2515.9840881470318,                 /* Expression: sps.A21
                                        * Referenced by: '<S78>/A21'
                                        */
  0.57651886298069654,                 /* Expression: sps.A22
                                        * Referenced by: '<S78>/A22'
                                        */
  7.8825943149034344E-6,               /* Expression: sps.B11
                                        * Referenced by: '<S79>/B11'
                                        */
  0.78825943149034827,                 /* Expression: sps.B21
                                        * Referenced by: '<S79>/B21'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S45>/Constant1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S38>/Constant2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S31>/Unit Delay3'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S38>/Gain1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S31>/Unit Delay1'
                                        */
  1.0E-5,                              /* Computed Parameter: Integrator_gainval_e
                                        * Referenced by: '<S111>/Integrator'
                                        */
  0.9226,                              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S111>/Integrator'
                                        */
  10.0,                                /* Expression: 10
                                        * Referenced by: '<S9>/Saturation'
                                        */
  -10.0,                               /* Expression: -10
                                        * Referenced by: '<S9>/Saturation'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S2>/Constant1'
                                        */

  /*  Computed Parameter: OpTrigger_P1_Size
   * Referenced by: '<S2>/OpTrigger'
   */
  { 1.0, 1.0 },
  26.0,                                /* Expression: Acq_Group
                                        * Referenced by: '<S2>/OpTrigger'
                                        */

  /*  Computed Parameter: OpTrigger_P2_Size
   * Referenced by: '<S2>/OpTrigger'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: Trig_Type
                                        * Referenced by: '<S2>/OpTrigger'
                                        */

  /*  Computed Parameter: OpTrigger_P3_Size
   * Referenced by: '<S2>/OpTrigger'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: Trig_Offset
                                        * Referenced by: '<S2>/OpTrigger'
                                        */

  /*  Computed Parameter: OpCtrl_P1_Size
   * Referenced by: '<S2>/OpCtrl'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: OpCtrl_P1
   * Referenced by: '<S2>/OpCtrl'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: OpCtrl_P2_Size
   * Referenced by: '<S2>/OpCtrl'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: boardid
                                        * Referenced by: '<S2>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P3_Size
   * Referenced by: '<S2>/OpCtrl'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: mode
                                        * Referenced by: '<S2>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P4_Size
   * Referenced by: '<S2>/OpCtrl'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: externalClock
                                        * Referenced by: '<S2>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P5_Size
   * Referenced by: '<S2>/OpCtrl'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: decimRtsi
                                        * Referenced by: '<S2>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P6_Size
   * Referenced by: '<S2>/OpCtrl'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S2>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P7_Size
   * Referenced by: '<S2>/OpCtrl'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: SampleTime
                                        * Referenced by: '<S2>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P8_Size
   * Referenced by: '<S2>/OpCtrl'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: calibIO
                                        * Referenced by: '<S2>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P9_Size
   * Referenced by: '<S2>/OpCtrl'
   */
  { 1.0, 1.0 },
  -1.0,                                /* Expression: numconfig
                                        * Referenced by: '<S2>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P10_Size
   * Referenced by: '<S2>/OpCtrl'
   */
  { 1.0, 1.0 },
  -1.0,                                /* Expression: loadinport
                                        * Referenced by: '<S2>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P11_Size
   * Referenced by: '<S2>/OpCtrl'
   */
  { 1.0, 1.0 },
  6.0,                                 /* Expression: BoardType
                                        * Referenced by: '<S2>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P12_Size
   * Referenced by: '<S2>/OpCtrl'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: sync_type
                                        * Referenced by: '<S2>/OpCtrl'
                                        */

  /*  Computed Parameter: default_Value
   * Referenced by: '<S17>/default'
   */
  { 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U },
  336592896U,                          /* Computed Parameter: addr4_Value
                                        * Referenced by: '<S17>/addr4'
                                        */
  269484032U,                          /* Computed Parameter: addr3_Value
                                        * Referenced by: '<S17>/addr3'
                                        */
  202375168U,                          /* Computed Parameter: addr2_Value
                                        * Referenced by: '<S17>/addr2'
                                        */
  135266304U,                          /* Computed Parameter: addr1_Value
                                        * Referenced by: '<S17>/addr1'
                                        */
  68157440U,                           /* Computed Parameter: addr_Value
                                        * Referenced by: '<S17>/addr'
                                        */
  255U,                                /* Computed Parameter: Saturation_UpperSat_k
                                        * Referenced by: '<S7>/Saturation'
                                        */
  1U,                                  /* Computed Parameter: Saturation_LowerSat_f
                                        * Referenced by: '<S7>/Saturation'
                                        */
  101U,                                /* Computed Parameter: sat_scn_UpperSat
                                        * Referenced by: '<S7>/sat_scn'
                                        */
  0U,                                  /* Computed Parameter: sat_scn_LowerSat
                                        * Referenced by: '<S7>/sat_scn'
                                        */
  374865680U,                          /* Computed Parameter: blockID_Value
                                        * Referenced by: '<S17>/blockID'
                                        */

  /*  Computed Parameter: load_config1_Value
   * Referenced by: '<S7>/load_config1'
   */
  { 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U },
  2147483648U,                         /* Computed Parameter: shift_2bits_Gain
                                        * Referenced by: '<S7>/shift_2bits'
                                        */
  0,                                   /* Expression: ~strcmp(FormatType,'Unsigned')
                                        * Referenced by: '<S22>/Constant'
                                        */
  1,                                   /* Expression: ~strcmp(FormatType,'Unsigned')
                                        * Referenced by: '<S20>/Constant'
                                        */
  1                                    /* Expression: ~strcmp(FormatType,'Unsigned')
                                        * Referenced by: '<S21>/Constant'
                                        */
};
