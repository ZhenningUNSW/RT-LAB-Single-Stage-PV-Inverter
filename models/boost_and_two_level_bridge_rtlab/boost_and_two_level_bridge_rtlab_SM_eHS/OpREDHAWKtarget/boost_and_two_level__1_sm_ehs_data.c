/*
 * boost_and_two_level__1_sm_ehs_data.c
 *
 * Code generation for model "boost_and_two_level__1_sm_ehs".
 *
 * Model version              : 1.1107
 * Simulink Coder version : 8.7 (R2014b) 08-Sep-2014
 * C source code generated on : Sat May 20 23:01:43 2017
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
                                        *   '<S69>/Constant4'
                                        *   '<S79>/Gain'
                                        *   '<S82>/Gain'
                                        *   '<S101>/Gain'
                                        *   '<S95>/Gain'
                                        *   '<S98>/Gain'
                                        */
  1.0,                                 /* Mask Parameter: PLL_AGC
                                        * Referenced by: '<S87>/Constant1'
                                        */
  2.0,                                 /* Mask Parameter: AlphaBetaZerotodq0_Alignment
                                        * Referenced by: '<S112>/Constant'
                                        */
  50.0,                                /* Mask Parameter: InverterControl_Fnom
                                        * Referenced by:
                                        *   '<S69>/Constant4'
                                        *   '<S76>/First cycle of simulation Id=0.92, Iq=0'
                                        */
  0.003,                               /* Mask Parameter: InverterControl_Increment_MPPT
                                        * Referenced by: '<S65>/Iph_3'
                                        */
  314.15926535897933,                  /* Mask Parameter: Discrete_Init
                                        * Referenced by: '<S89>/Discrete-Time Integrator'
                                        */
  1.0,                                 /* Mask Parameter: Discrete_Kd
                                        * Referenced by: '<S89>/Discrete Derivative '
                                        */
  6.6,                                 /* Mask Parameter: InverterControl_Ki_Ireg
                                        * Referenced by: '<S71>/Integral Gain'
                                        */
  200.0,                               /* Mask Parameter: InverterControl_Ki_VDCreg
                                        * Referenced by: '<S142>/Integral Gain'
                                        */
  180.0,                               /* Mask Parameter: Discrete_Kp
                                        * Referenced by: '<S89>/Kp4'
                                        */
  0.15,                                /* Mask Parameter: InverterControl_Kp_Ireg
                                        * Referenced by: '<S71>/Proportional Gain'
                                        */
  12.0,                                /* Mask Parameter: InverterControl_Kp_VDCreg
                                        * Referenced by: '<S142>/Proportional Gain'
                                        */
  -1.5,                                /* Mask Parameter: PI_LowerSaturationLimit
                                        * Referenced by: '<S71>/Saturate'
                                        */
  -1.5,                                /* Mask Parameter: PI_LowerSaturationLimit_j
                                        * Referenced by: '<S142>/Saturate'
                                        */

  /*  Mask Parameter: PWM_Generator_MinMax
   * Referenced by: '<S68>/Constant10'
   */
  { -1.0, 1.0 },
  3500.0,                              /* Mask Parameter: InverterControl_Pnom
                                        * Referenced by: '<S67>/A->pu'
                                        */
  1.0,                                 /* Mask Parameter: RMS_TrueRMS
                                        * Referenced by: '<S6>/Constant'
                                        */
  1.0,                                 /* Mask Parameter: RMS1_TrueRMS
                                        * Referenced by: '<S7>/Constant'
                                        */
  1.5,                                 /* Mask Parameter: PI_UpperSaturationLimit
                                        * Referenced by: '<S71>/Saturate'
                                        */
  1.5,                                 /* Mask Parameter: PI_UpperSaturationLimit_a
                                        * Referenced by: '<S142>/Saturate'
                                        */
  425.0,                               /* Mask Parameter: InverterControl_Vdc_ref_Init
                                        * Referenced by:
                                        *   '<S62>/Vnom_dc1'
                                        *   '<S65>/Iph_'
                                        */
  400.0,                               /* Mask Parameter: InverterControl_Vnom_dc
                                        * Referenced by: '<S70>/Rtot_pu2'
                                        */
  240.0,                               /* Mask Parameter: InverterControl_Vnom_prim
                                        * Referenced by:
                                        *   '<S67>/A->pu'
                                        *   '<S67>/V->pu'
                                        *   '<S69>/Constant3'
                                        */
  1.0,                                 /* Mask Parameter: CompareToConstant_const
                                        * Referenced by: '<S114>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: CompareToConstant1_const
                                        * Referenced by: '<S115>/Constant'
                                        */
  1.0,                                 /* Mask Parameter: CompareToConstant_const_d
                                        * Referenced by: '<S3>/Constant'
                                        */
  5.0,                                 /* Mask Parameter: Counter_max_count
                                        * Referenced by: '<S42>/Switch'
                                        */
  127U,                                /* Mask Parameter: BitwiseOperator_BitMask
                                        * Referenced by: '<S45>/Bitwise Operator'
                                        */
  65535U,                              /* Mask Parameter: BitwiseOperator_BitMask_j
                                        * Referenced by: '<S43>/Bitwise Operator'
                                        */
  65535U,                              /* Mask Parameter: BitwiseOperator_BitMask_h
                                        * Referenced by: '<S44>/Bitwise Operator'
                                        */
  3U,                                  /* Mask Parameter: DetectChange_vinit
                                        * Referenced by: '<S48>/Delay Input1'
                                        */
  0U,                                  /* Mask Parameter: DetectChange_vinit_i
                                        * Referenced by: '<S47>/Delay Input1'
                                        */
  0U,                                  /* Mask Parameter: DetectChange_vinit_g
                                        * Referenced by: '<S51>/Delay Input1'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S38>/Gain'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S39>/Gain'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S42>/Constant1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S42>/Constant2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S42>/Constant'
                                        */
  31.99951171875,                      /* Expression: 2^(N-Q)-2^(-Q)
                                        * Referenced by: '<S43>/Saturation1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S43>/Saturation1'
                                        */
  2048.0,                              /* Expression: 2^Q
                                        * Referenced by: '<S43>/Gain2'
                                        */
  15.99951171875,                      /* Expression: 2^(N-Q-1)-2^(-Q)
                                        * Referenced by: '<S43>/Saturation'
                                        */
  -16.0,                               /* Expression: -2^(N-Q-1)
                                        * Referenced by: '<S43>/Saturation'
                                        */
  2048.0,                              /* Expression: 2^Q
                                        * Referenced by: '<S43>/Gain1'
                                        */
  31.99951171875,                      /* Expression: 2^(N-Q)-2^(-Q)
                                        * Referenced by: '<S44>/Saturation1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S44>/Saturation1'
                                        */
  2048.0,                              /* Expression: 2^Q
                                        * Referenced by: '<S44>/Gain2'
                                        */
  15.99951171875,                      /* Expression: 2^(N-Q-1)-2^(-Q)
                                        * Referenced by: '<S44>/Saturation'
                                        */
  -16.0,                               /* Expression: -2^(N-Q-1)
                                        * Referenced by: '<S44>/Saturation'
                                        */
  2048.0,                              /* Expression: 2^Q
                                        * Referenced by: '<S44>/Gain1'
                                        */
  127.0,                               /* Expression: 2^(N-Q)-2^(-Q)
                                        * Referenced by: '<S45>/Saturation1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S45>/Saturation1'
                                        */
  1.0,                                 /* Expression: 2^Q
                                        * Referenced by: '<S45>/Gain2'
                                        */
  63.0,                                /* Expression: 2^(N-Q-1)-2^(-Q)
                                        * Referenced by: '<S45>/Saturation'
                                        */
  -64.0,                               /* Expression: -2^(N-Q-1)
                                        * Referenced by: '<S45>/Saturation'
                                        */
  1.0,                                 /* Expression: 2^Q
                                        * Referenced by: '<S45>/Gain1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S80>/Gain1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S83>/Gain1'
                                        */
  -1.1368683772160957E-16,             /* Expression: sps.K2
                                        * Referenced by: '<S85>/Gain1'
                                        */
  6.4623485355705288E-30,              /* Expression: sps.K1
                                        * Referenced by: '<S85>/Gain'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S96>/Gain1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S99>/Gain1'
                                        */
  1.0,                                 /* Expression: [1]
                                        * Referenced by: '<S88>/Gain'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S92>/Gain1'
                                        */
  1.0E-5,                              /* Computed Parameter: Integ4_gainval
                                        * Referenced by: '<S95>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S95>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S95>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S95>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size
   * Referenced by: '<S97>/S-Function'
   */
  { 1.0, 1.0 },
  0.02502,                             /* Expression: MaxDelay
                                        * Referenced by: '<S97>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size
   * Referenced by: '<S97>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S97>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size
   * Referenced by: '<S97>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S97>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size
   * Referenced by: '<S97>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S97>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S96>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S95>/Constant'
                                        */
  1.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S95>/Unit Delay1'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S92>/Gain3'
                                        */
  1.0E-5,                              /* Computed Parameter: Integ4_gainval_c
                                        * Referenced by: '<S98>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S98>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S98>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S98>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_h
   * Referenced by: '<S100>/S-Function'
   */
  { 1.0, 1.0 },
  0.02502,                             /* Expression: MaxDelay
                                        * Referenced by: '<S100>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_p
   * Referenced by: '<S100>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S100>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_p
   * Referenced by: '<S100>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S100>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_p
   * Referenced by: '<S100>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S100>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S99>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S98>/Constant'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S98>/Unit Delay1'
                                        */
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<S92>/Rad->Deg.'
                                        */
  0.0,                                 /* Expression: inf
                                        * Referenced by: '<S88>/Saturation'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S88>/Saturation'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S102>/Gain1'
                                        */

  /*  Expression: [0.92 0]
   * Referenced by: '<S76>/Constant'
   */
  { 0.92, 0.0 },

  /*  Expression: [0,0]
   * Referenced by: '<S116>/dq'
   */
  { 0.0, 0.0 },

  /*  Expression: [0,0]
   * Referenced by: '<S117>/dq'
   */
  { 0.0, 0.0 },
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S1>/S-Function1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S1>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory1'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_m
   * Referenced by: '<S57>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: Acqu_group
                                        * Referenced by: '<S57>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_n
   * Referenced by: '<S15>/S-Function'
   */
  { 1.0, 1.0 },
  165.0,                               /* Expression: Data_width
                                        * Referenced by: '<S15>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_n
   * Referenced by: '<S15>/S-Function'
   */
  { 165.0, 1.0 },

  /*  Expression: InitialConditions
   * Referenced by: '<S15>/S-Function'
   */
  { 800.0, 25.0, 1.0, 0.0, 49.0, 48.0, 52.0, 53.0, 53.0, 9.0, 1.0, 3.0, 12.0,
    13.0, 0.0, 2.0, 48.0, 49.0, 50.0, 51.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0,
    54.0, 54.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 0.0066, 0.33,
    0.08, 0.004125, 0.0048529, 1.0, 0.005, 0.005, 0.1, 0.1, 0.1, 0.1, 0.005, 0.1,
    0.1, 0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.65, 1.65, 1.65, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 3.3, 3.3, 3.3, 3.3, 3.3, 5.0, 5.0, 5.0, 5.0, 5.0,
    5.0, 5.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0,
    16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 0.0, 0.0, 0.0, 0.0,
    0.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -16.0, -16.0, -16.0, -16.0,
    -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0,
    -16.0, -16.0, -16.0, -16.0, -16.0, 2.0 },
  0.5,                                 /* Expression: .5
                                        * Referenced by: '<S10>/Relay1'
                                        */
  0.5,                                 /* Expression: .5
                                        * Referenced by: '<S10>/Relay1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S10>/Relay1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S10>/Relay1'
                                        */

  /*  Computed Parameter: eHS_rst_loadin_P1_Size
   * Referenced by: '<S10>/eHS_rst_loadin'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: eHS_rst_loadin_P1
   * Referenced by: '<S10>/eHS_rst_loadin'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: eHS_rst_loadin_P2_Size
   * Referenced by: '<S10>/eHS_rst_loadin'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: FcnNos
                                        * Referenced by: '<S10>/eHS_rst_loadin'
                                        */

  /*  Computed Parameter: eHS_rst_loadin_P3_Size
   * Referenced by: '<S10>/eHS_rst_loadin'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: width
                                        * Referenced by: '<S10>/eHS_rst_loadin'
                                        */

  /*  Computed Parameter: eHS_rst_loadin_P4_Size
   * Referenced by: '<S10>/eHS_rst_loadin'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S10>/eHS_rst_loadin'
                                        */

  /*  Computed Parameter: Automated_Solver_Mat_Initialisa
   * Referenced by: '<S10>/Automated_Solver_Mat_Initialisation_1'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: Automated_Solver_Mat_Initiali_c
   * Referenced by: '<S10>/Automated_Solver_Mat_Initialisation_1'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: Automated_Solver_Mat_Initiali_o
   * Referenced by: '<S10>/Automated_Solver_Mat_Initialisation_1'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: fpga_port_in
                                        * Referenced by: '<S10>/Automated_Solver_Mat_Initialisation_1'
                                        */

  /*  Computed Parameter: Automated_Solver_Mat_Initiali_e
   * Referenced by: '<S10>/Automated_Solver_Mat_Initialisation_1'
   */
  { 1.0, 31.0 },

  /*  Computed Parameter: Automated_Solver_Mat_Initial_er
   * Referenced by: '<S10>/Automated_Solver_Mat_Initialisation_1'
   */
  { 101.0, 104.0, 115.0, 95.0, 99.0, 111.0, 110.0, 102.0, 105.0, 103.0, 95.0,
    101.0, 72.0, 83.0, 120.0, 54.0, 52.0, 118.0, 51.0, 95.0, 83.0, 73.0, 68.0,
    52.0, 56.0, 57.0, 52.0, 46.0, 109.0, 97.0, 116.0 },
  20364.675298172569,                  /* Expression: 14400*sqrt(2)
                                        * Referenced by: '<S12>/Sine Wave Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S12>/Sine Wave Function'
                                        */
  314.15926535897933,                  /* Expression: 50*2*pi
                                        * Referenced by: '<S12>/Sine Wave Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S12>/Sine Wave Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S12>/Memory1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S12>/Memory'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P1_Size
   * Referenced by: '<S10>/Inputs_eHS1_Send'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: Inputs_eHS1_Send_P1
   * Referenced by: '<S10>/Inputs_eHS1_Send'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: Inputs_eHS1_Send_P2_Size
   * Referenced by: '<S10>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: FcnNos
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P3_Size
   * Referenced by: '<S10>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  3.0,                                 /* Expression: width
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P4_Size
   * Referenced by: '<S10>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P5_Size
   * Referenced by: '<S10>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: checkVersion
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P6_Size
   * Referenced by: '<S10>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedId
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P7_Size
   * Referenced by: '<S10>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedVersion
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P8_Size
   * Referenced by: '<S10>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: opComp
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P1_Size
   * Referenced by: '<S10>/Outputs_eHS1_Recv'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: Outputs_eHS1_Recv_P1
   * Referenced by: '<S10>/Outputs_eHS1_Recv'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: Outputs_eHS1_Recv_P2_Size
   * Referenced by: '<S10>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: FcnNos
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P3_Size
   * Referenced by: '<S10>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  8.0,                                 /* Expression: width
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P4_Size
   * Referenced by: '<S10>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P5_Size
   * Referenced by: '<S10>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: sampleTime
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P6_Size
   * Referenced by: '<S10>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: checkVersion
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P7_Size
   * Referenced by: '<S10>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedId
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P8_Size
   * Referenced by: '<S10>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedVersion
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P9_Size
   * Referenced by: '<S10>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: opComp
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
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

  /*  Computed Parameter: SFunction_P1_Size_g
   * Referenced by: '<S58>/S-Function'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: Acqu_group
                                        * Referenced by: '<S58>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory2'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_b
   * Referenced by: '<S59>/S-Function'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: Acqu_group
                                        * Referenced by: '<S59>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_b0
   * Referenced by: '<S60>/S-Function'
   */
  { 1.0, 1.0 },
  3.0,                                 /* Expression: Acqu_group
                                        * Referenced by: '<S60>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory3'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory3'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_k
   * Referenced by: '<S61>/S-Function'
   */
  { 1.0, 1.0 },
  5.0,                                 /* Expression: Acqu_group
                                        * Referenced by: '<S61>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S42>/Memory1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S42>/Memory2'
                                        */

  /*  Computed Parameter: LoadIn_P1_Size
   * Referenced by: '<S40>/LoadIn'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: LoadIn_P1
   * Referenced by: '<S40>/LoadIn'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: LoadIn_P2_Size
   * Referenced by: '<S40>/LoadIn'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: FcnNos
                                        * Referenced by: '<S40>/LoadIn'
                                        */

  /*  Computed Parameter: LoadIn_P3_Size
   * Referenced by: '<S40>/LoadIn'
   */
  { 1.0, 1.0 },
  18.0,                                /* Expression: width
                                        * Referenced by: '<S40>/LoadIn'
                                        */

  /*  Computed Parameter: LoadIn_P4_Size
   * Referenced by: '<S40>/LoadIn'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S40>/LoadIn'
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

  /*  Expression: ones(1,9)
   * Referenced by: '<S41>/Constant'
   */
  { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 },

  /*  Computed Parameter: DataInSend_P1_Size
   * Referenced by: '<S40>/DataIn Send'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: DataInSend_P1
   * Referenced by: '<S40>/DataIn Send'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: DataInSend_P2_Size
   * Referenced by: '<S40>/DataIn Send'
   */
  { 1.0, 1.0 },
  15.0,                                /* Expression: FcnNos
                                        * Referenced by: '<S40>/DataIn Send'
                                        */

  /*  Computed Parameter: DataInSend_P3_Size
   * Referenced by: '<S40>/DataIn Send'
   */
  { 1.0, 1.0 },
  16.0,                                /* Expression: width
                                        * Referenced by: '<S40>/DataIn Send'
                                        */

  /*  Computed Parameter: DataInSend_P4_Size
   * Referenced by: '<S40>/DataIn Send'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S40>/DataIn Send'
                                        */

  /*  Computed Parameter: DataInSend_P5_Size
   * Referenced by: '<S40>/DataIn Send'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: checkVersion
                                        * Referenced by: '<S40>/DataIn Send'
                                        */

  /*  Computed Parameter: DataInSend_P6_Size
   * Referenced by: '<S40>/DataIn Send'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedId
                                        * Referenced by: '<S40>/DataIn Send'
                                        */

  /*  Computed Parameter: DataInSend_P7_Size
   * Referenced by: '<S40>/DataIn Send'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedVersion
                                        * Referenced by: '<S40>/DataIn Send'
                                        */

  /*  Computed Parameter: DataInSend_P8_Size
   * Referenced by: '<S40>/DataIn Send'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: opComp
                                        * Referenced by: '<S40>/DataIn Send'
                                        */

  /*  Computed Parameter: rtlab_io_block_P1_Size
   * Referenced by: '<S13>/rtlab_io_block'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: rtlab_io_block_P1
   * Referenced by: '<S13>/rtlab_io_block'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: rtlab_io_block_P2_Size
   * Referenced by: '<S13>/rtlab_io_block'
   */
  { 0.0, 0.0 },

  /*  Computed Parameter: rtlab_io_block_P3_Size
   * Referenced by: '<S13>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: sampleTime
                                        * Referenced by: '<S13>/rtlab_io_block'
                                        */

  /*  Computed Parameter: rtlab_io_block_P4_Size
   * Referenced by: '<S13>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  5.0,                                 /* Expression: portNb
                                        * Referenced by: '<S13>/rtlab_io_block'
                                        */

  /*  Computed Parameter: rtlab_io_block_P5_Size
   * Referenced by: '<S13>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  8.0,                                 /* Expression: numchan
                                        * Referenced by: '<S13>/rtlab_io_block'
                                        */

  /*  Computed Parameter: rtlab_io_block_P6_Size
   * Referenced by: '<S13>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  8.0,                                 /* Expression: maxCount
                                        * Referenced by: '<S13>/rtlab_io_block'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S13>/IOTypeSel1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S13>/Memory1'
                                        */

  /*  Computed Parameter: IOTypeSel_LoadIn_P1_Size
   * Referenced by: '<S13>/IOTypeSel_LoadIn'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: IOTypeSel_LoadIn_P1
   * Referenced by: '<S13>/IOTypeSel_LoadIn'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: IOTypeSel_LoadIn_P2_Size
   * Referenced by: '<S13>/IOTypeSel_LoadIn'
   */
  { 1.0, 1.0 },
  13.0,                                /* Expression: FcnNos
                                        * Referenced by: '<S13>/IOTypeSel_LoadIn'
                                        */

  /*  Computed Parameter: IOTypeSel_LoadIn_P3_Size
   * Referenced by: '<S13>/IOTypeSel_LoadIn'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: width
                                        * Referenced by: '<S13>/IOTypeSel_LoadIn'
                                        */

  /*  Computed Parameter: IOTypeSel_LoadIn_P4_Size
   * Referenced by: '<S13>/IOTypeSel_LoadIn'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S13>/IOTypeSel_LoadIn'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S13>/Memory'
                                        */

  /*  Computed Parameter: rtlab_io_block_P1_Size_b
   * Referenced by: '<S14>/rtlab_io_block'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: rtlab_io_block_P1_m
   * Referenced by: '<S14>/rtlab_io_block'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: rtlab_io_block_P2_Size_k
   * Referenced by: '<S14>/rtlab_io_block'
   */
  { 0.0, 0.0 },

  /*  Computed Parameter: rtlab_io_block_P3_Size_i
   * Referenced by: '<S14>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: sampleTime
                                        * Referenced by: '<S14>/rtlab_io_block'
                                        */

  /*  Computed Parameter: rtlab_io_block_P4_Size_g
   * Referenced by: '<S14>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  6.0,                                 /* Expression: portNb
                                        * Referenced by: '<S14>/rtlab_io_block'
                                        */

  /*  Computed Parameter: rtlab_io_block_P5_Size_j
   * Referenced by: '<S14>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  8.0,                                 /* Expression: numchan
                                        * Referenced by: '<S14>/rtlab_io_block'
                                        */

  /*  Computed Parameter: rtlab_io_block_P6_Size_k
   * Referenced by: '<S14>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  8.0,                                 /* Expression: maxCount
                                        * Referenced by: '<S14>/rtlab_io_block'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S14>/IOTypeSel1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S14>/Memory1'
                                        */

  /*  Computed Parameter: IOTypeSel_LoadIn_P1_Size_e
   * Referenced by: '<S14>/IOTypeSel_LoadIn'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: IOTypeSel_LoadIn_P1_l
   * Referenced by: '<S14>/IOTypeSel_LoadIn'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: IOTypeSel_LoadIn_P2_Size_o
   * Referenced by: '<S14>/IOTypeSel_LoadIn'
   */
  { 1.0, 1.0 },
  14.0,                                /* Expression: FcnNos
                                        * Referenced by: '<S14>/IOTypeSel_LoadIn'
                                        */

  /*  Computed Parameter: IOTypeSel_LoadIn_P3_Size_m
   * Referenced by: '<S14>/IOTypeSel_LoadIn'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: width
                                        * Referenced by: '<S14>/IOTypeSel_LoadIn'
                                        */

  /*  Computed Parameter: IOTypeSel_LoadIn_P4_Size_e
   * Referenced by: '<S14>/IOTypeSel_LoadIn'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S14>/IOTypeSel_LoadIn'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S14>/Memory'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S38>/integrator'
                                        */
  0.02,                                /* Expression: 1./sps.Freq
                                        * Referenced by: '<S38>/Transport Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S38>/Transport Delay'
                                        */
  0.02,                                /* Expression: 1./sps.Freq
                                        * Referenced by: '<S38>/K1'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S38>/Memory'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S39>/integrator'
                                        */
  0.02,                                /* Expression: 1./sps.Freq
                                        * Referenced by: '<S39>/Transport Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S39>/Transport Delay'
                                        */
  0.02,                                /* Expression: 1./sps.Freq
                                        * Referenced by: '<S39>/K1'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S39>/Memory'
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
  100.0,                               /* Expression: Decimation
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
  1.048576E+6,                         /* Expression: Buffer_size
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
  1.073741824E+9,                      /* Expression: file_size
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
  0.05,                                /* Expression: 0.05
                                        * Referenced by: '<S12>/Step'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S12>/Step'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S12>/Step'
                                        */

  /*  Computed Parameter: RTEConversion_P1_Size
   * Referenced by: '<S12>/RTE Conversion'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: nbMaxEvents
                                        * Referenced by: '<S12>/RTE Conversion'
                                        */

  /*  Computed Parameter: RTEConversion_P2_Size
   * Referenced by: '<S12>/RTE Conversion'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: inputdatatype
                                        * Referenced by: '<S12>/RTE Conversion'
                                        */

  /*  Computed Parameter: RTEConversion_P3_Size
   * Referenced by: '<S12>/RTE Conversion'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: outputdatatype
                                        * Referenced by: '<S12>/RTE Conversion'
                                        */

  /*  Computed Parameter: RTEConversion_P4_Size
   * Referenced by: '<S12>/RTE Conversion'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: compensation
                                        * Referenced by: '<S12>/RTE Conversion'
                                        */

  /*  Computed Parameter: RTEConversion_P5_Size
   * Referenced by: '<S12>/RTE Conversion'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: sampleTime
                                        * Referenced by: '<S12>/RTE Conversion'
                                        */
  0.1684,                              /* Expression: 0.1684
                                        * Referenced by: '<S62>/Unit Delay'
                                        */
  3.125E-5,                            /* Expression: sps.Delay
                                        * Referenced by: '<S140>/Constant3'
                                        */
  0.000125,                            /* Expression: sps.Period
                                        * Referenced by: '<S140>/Constant1'
                                        */
  8000.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S140>/1\ib1'
                                        */

  /*  Expression: [0 .5 1]
   * Referenced by: '<S140>/Lookup Table'
   */
  { 0.0, 0.5, 1.0 },

  /*  Expression: [0 2 0]
   * Referenced by: '<S140>/Lookup Table'
   */
  { 0.0, 2.0, 0.0 },
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S140>/Constant2'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S118>/Gain1'
                                        */

  /*  Computed Parameter: RTEConversion1_P1_Size
   * Referenced by: '<S12>/RTE Conversion1'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: nbMaxEvents
                                        * Referenced by: '<S12>/RTE Conversion1'
                                        */

  /*  Computed Parameter: RTEConversion1_P2_Size
   * Referenced by: '<S12>/RTE Conversion1'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: inputdatatype
                                        * Referenced by: '<S12>/RTE Conversion1'
                                        */

  /*  Computed Parameter: RTEConversion1_P3_Size
   * Referenced by: '<S12>/RTE Conversion1'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: outputdatatype
                                        * Referenced by: '<S12>/RTE Conversion1'
                                        */

  /*  Computed Parameter: RTEConversion1_P4_Size
   * Referenced by: '<S12>/RTE Conversion1'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: compensation
                                        * Referenced by: '<S12>/RTE Conversion1'
                                        */

  /*  Computed Parameter: RTEConversion1_P5_Size
   * Referenced by: '<S12>/RTE Conversion1'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: sampleTime
                                        * Referenced by: '<S12>/RTE Conversion1'
                                        */

  /*  Computed Parameter: RTELogicalOperator1_P1_Size
   * Referenced by: '<S12>/RTE Logical Operator1'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: LogicalOperator
                                        * Referenced by: '<S12>/RTE Logical Operator1'
                                        */

  /*  Computed Parameter: RTELogicalOperator1_P2_Size
   * Referenced by: '<S12>/RTE Logical Operator1'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: NbrInput
                                        * Referenced by: '<S12>/RTE Logical Operator1'
                                        */

  /*  Computed Parameter: RTELogicalOperator1_P3_Size
   * Referenced by: '<S12>/RTE Logical Operator1'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: NbrMaxEvents
                                        * Referenced by: '<S12>/RTE Logical Operator1'
                                        */

  /*  Computed Parameter: RTESPWM_P1_Size
   * Referenced by: '<S12>/RTE SPWM'
   */
  { 1.0, 1.0 },
  100.0,                               /* Expression: MaxEvents
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */

  /*  Computed Parameter: RTESPWM_P2_Size
   * Referenced by: '<S12>/RTE SPWM'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: SampleTime
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */

  /*  Computed Parameter: RTESPWM_P3_Size
   * Referenced by: '<S12>/RTE SPWM'
   */
  { 1.0, 1.0 },
  25000.0,                             /* Expression: MaxFreq
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */

  /*  Computed Parameter: RTESPWM_P4_Size
   * Referenced by: '<S12>/RTE SPWM'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: MinFreq
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */

  /*  Computed Parameter: RTESPWM_P5_Size
   * Referenced by: '<S12>/RTE SPWM'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: EnablingPort
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */

  /*  Computed Parameter: RTESPWM_P6_Size
   * Referenced by: '<S12>/RTE SPWM'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: NumberPhases
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */

  /*  Computed Parameter: RTESPWM_P7_Size
   * Referenced by: '<S12>/RTE SPWM'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: ComplementaryMode
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */

  /*  Computed Parameter: RTESPWM_P8_Size
   * Referenced by: '<S12>/RTE SPWM'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: RiseTimeDelay
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */

  /*  Computed Parameter: RTESPWM_P9_Size
   * Referenced by: '<S12>/RTE SPWM'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: CenterAlignmentMode
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */

  /*  Computed Parameter: RTESPWM_P10_Size
   * Referenced by: '<S12>/RTE SPWM'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: SpaceVector
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S12>/Switch'
                                        */

  /*  Computed Parameter: RTEGround_P1_Size
   * Referenced by: '<S12>/RTE Ground'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: SampleTime
                                        * Referenced by: '<S12>/RTE Ground'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P1_Size
   * Referenced by: '<S49>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: conversionType
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P2_Size
   * Referenced by: '<S49>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  8.0,                                 /* Expression: nbLine
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P3_Size
   * Referenced by: '<S49>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: nbEvents
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P4_Size
   * Referenced by: '<S49>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: timeUnit
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P5_Size
   * Referenced by: '<S49>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: enTimeFactor
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P6_Size
   * Referenced by: '<S49>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: sampleTime
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P7_Size
   * Referenced by: '<S49>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: initStates
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: EventGen_eHS_1_P1_Size
   * Referenced by: '<S49>/EventGen_eHS_1'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: EventGen_eHS_1_P1
   * Referenced by: '<S49>/EventGen_eHS_1'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: EventGen_eHS_1_P2_Size
   * Referenced by: '<S49>/EventGen_eHS_1'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portNb
                                        * Referenced by: '<S49>/EventGen_eHS_1'
                                        */

  /*  Computed Parameter: EventGen_eHS_1_P3_Size
   * Referenced by: '<S49>/EventGen_eHS_1'
   */
  { 0.0, 0.0 },

  /*  Computed Parameter: EventGen_eHS_1_P4_Size
   * Referenced by: '<S49>/EventGen_eHS_1'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: size
                                        * Referenced by: '<S49>/EventGen_eHS_1'
                                        */

  /*  Computed Parameter: EventGen_eHS_1_P5_Size
   * Referenced by: '<S49>/EventGen_eHS_1'
   */
  { 1.0, 1.0 },
  8.0,                                 /* Expression: nbChannels
                                        * Referenced by: '<S49>/EventGen_eHS_1'
                                        */

  /*  Computed Parameter: EventGen_eHS_1_P6_Size
   * Referenced by: '<S49>/EventGen_eHS_1'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: numwidth
                                        * Referenced by: '<S49>/EventGen_eHS_1'
                                        */

  /*  Computed Parameter: EventGen_eHS_1_P7_Size
   * Referenced by: '<S49>/EventGen_eHS_1'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: timeUnit
                                        * Referenced by: '<S49>/EventGen_eHS_1'
                                        */
  8.55,                                /* Expression: 8.55
                                        * Referenced by: '<S12>/Constant1'
                                        */
  523.6,                               /* Expression: 37.4*14
                                        * Referenced by: '<S12>/Constant2'
                                        */
  0.06,                                /* Expression: 0.06
                                        * Referenced by: '<S12>/Constant3'
                                        */
  840.0,                               /* Expression: 60*14
                                        * Referenced by: '<S12>/Constant4'
                                        */
  50.0,                                /* Expression: sps.Finit
                                        * Referenced by: '<S87>/Unit Delay'
                                        */
  70.0,                                /* Expression: 70
                                        * Referenced by: '<S76>/avoid division by zero'
                                        */
  40.0,                                /* Expression: 40
                                        * Referenced by: '<S76>/avoid division by zero'
                                        */
  0.25,                                /* Expression: 1/4
                                        * Referenced by: '<S76>/Gain'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_i
   * Referenced by: '<S113>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: MaxDelay
                                        * Referenced by: '<S113>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_i
   * Referenced by: '<S113>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S113>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_i
   * Referenced by: '<S113>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S113>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_n
   * Referenced by: '<S113>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S113>/S-Function'
                                        */
  2.0E-5,                              /* Computed Parameter: DiscreteTimeIntegrator_gainval
                                        * Referenced by: '<S87>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: sps.Phase_Init*pi/180
                                        * Referenced by: '<S87>/Discrete-Time Integrator'
                                        */
  6.2831853071795862,                  /* Expression: 2*pi
                                        * Referenced by: '<S87>/Constant4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S76>/First cycle of simulation Id=0.92, Iq=0'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S76>/First cycle of simulation Id=0.92, Iq=0'
                                        */
  0.01257992044073508,                 /* Expression: sps.D
                                        * Referenced by: '<S108>/D*u(k)'
                                        */

  /*  Expression: sps.x0(1,:)
   * Referenced by: '<S108>/Delay_x1'
   */
  { 0.00028197057830512594, 0.0 },
  3151.6695002188571,                  /* Expression: sps.C11
                                        * Referenced by: '<S111>/C11'
                                        */

  /*  Expression: sps.x0(2,:)
   * Referenced by: '<S108>/Delay_x2'
   */
  { -5.5511151231257827E-17, 0.0 },
  0.02515984088147016,                 /* Expression: sps.C12
                                        * Referenced by: '<S111>/C12'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S62>/Unit Delay2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S62>/Iq_ref'
                                        */
  2.0E-5,                              /* Computed Parameter: Integrator_gainval
                                        * Referenced by: '<S71>/Integrator'
                                        */
  0.0,                                 /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S71>/Integrator'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S72>/Gain1'
                                        */
  1.0E-5,                              /* Computed Parameter: Integ4_gainval_l
                                        * Referenced by: '<S79>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S79>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S79>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S79>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_l
   * Referenced by: '<S81>/S-Function'
   */
  { 1.0, 1.0 },
  0.022242222222222222,                /* Expression: MaxDelay
                                        * Referenced by: '<S81>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_c
   * Referenced by: '<S81>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S81>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_b
   * Referenced by: '<S81>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S81>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_f
   * Referenced by: '<S81>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S81>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S80>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S79>/Constant'
                                        */
  1.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S79>/Unit Delay1'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S72>/Gain3'
                                        */
  1.0E-5,                              /* Computed Parameter: Integ4_gainval_g
                                        * Referenced by: '<S82>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S82>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S82>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S82>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_hi
   * Referenced by: '<S84>/S-Function'
   */
  { 1.0, 1.0 },
  0.022242222222222222,                /* Expression: MaxDelay
                                        * Referenced by: '<S84>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_j
   * Referenced by: '<S84>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S84>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_c
   * Referenced by: '<S84>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S84>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_m
   * Referenced by: '<S84>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S84>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S83>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S82>/Constant'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S82>/Unit Delay1'
                                        */
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<S72>/Rad->Deg.'
                                        */
  0.017453292519943295,                /* Expression: pi/180
                                        * Referenced by: '<S67>/to-rad'
                                        */
  0.002,                               /* Expression: RLff(1)
                                        * Referenced by: '<S64>/Rff '
                                        */
  0.2,                                 /* Expression: RLff(2)
                                        * Referenced by: '<S64>/Lff  '
                                        */
  0.002,                               /* Expression: RLff(1)
                                        * Referenced by: '<S64>/Rff'
                                        */
  0.2,                                 /* Expression: RLff(2)
                                        * Referenced by: '<S64>/Lff'
                                        */
  1.5,                                 /* Expression: 1.5
                                        * Referenced by: '<S64>/Saturation'
                                        */
  -1.5,                                /* Expression: -1.5
                                        * Referenced by: '<S64>/Saturation'
                                        */
  450.0,                               /* Expression: Limits_MPPT(1)
                                        * Referenced by: '<S65>/Iph_1'
                                        */
  375.0,                               /* Expression: Limits_MPPT(2)
                                        * Referenced by: '<S65>/Iph_2'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S12>/MPPT_On'
                                        */
  1.0E-5,                              /* Computed Parameter: Integ4_gainval_o
                                        * Referenced by: '<S85>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S85>/Integ4'
                                        */
  0.02,                                /* Expression: sps.Delay
                                        * Referenced by: '<S85>/K1'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_o
   * Referenced by: '<S86>/S-Function'
   */
  { 1.0, 1.0 },
  0.02004,                             /* Expression: MaxDelay
                                        * Referenced by: '<S86>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_id
   * Referenced by: '<S86>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S86>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_l
   * Referenced by: '<S86>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S86>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_i
   * Referenced by: '<S86>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S86>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S85>/Unit Delay'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S85>/K2'
                                        */
  425.0,                               /* Expression: sps.Vinit
                                        * Referenced by: '<S85>/Unit Delay1'
                                        */
  1.0E-5,                              /* Computed Parameter: Integ4_gainval_b
                                        * Referenced by: '<S101>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S101>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S101>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S101>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_lc
   * Referenced by: '<S103>/S-Function'
   */
  { 1.0, 1.0 },
  0.02502,                             /* Expression: MaxDelay
                                        * Referenced by: '<S103>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_g
   * Referenced by: '<S103>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S103>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_f
   * Referenced by: '<S103>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S103>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_a
   * Referenced by: '<S103>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S103>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S102>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S101>/Constant'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S101>/Unit Delay1'
                                        */

  /*  Expression: [ TcD  Ts-TcD ]
   * Referenced by: '<S89>/Discrete Derivative '
   */
  { 0.0001, -8.0E-5 },
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S89>/Discrete Derivative '
                                        */
  0.064,                               /* Computed Parameter: DiscreteTimeIntegrator_gainva_a
                                        * Referenced by: '<S89>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: Par_Limits(1)
                                        * Referenced by: '<S89>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: Par_Limits(2)
                                        * Referenced by: '<S89>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: Par_Limits(1)
                                        * Referenced by: '<S89>/Saturation1'
                                        */
  0.0,                                 /* Expression: Par_Limits(2)
                                        * Referenced by: '<S89>/Saturation1'
                                        */
  0.15915494309189535,                 /* Expression: 1/2/pi
                                        * Referenced by: '<S87>/Gain10'
                                        */
  0.00024000000000000003,              /* Computed Parameter: RateLimiter_RisingLim
                                        * Referenced by: '<S87>/Rate Limiter'
                                        */
  -0.00024000000000000003,             /* Computed Parameter: RateLimiter_FallingLim
                                        * Referenced by: '<S87>/Rate Limiter'
                                        */
  50.0,                                /* Expression: sps.Finit
                                        * Referenced by: '<S87>/Rate Limiter'
                                        */
  101.32101697571224,                  /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S104>/Delay_x1'
                                        */
  0.999995076138259,                   /* Expression: sps.A11
                                        * Referenced by: '<S105>/A11'
                                        */
  0.0,                                 /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S104>/Delay_x2'
                                        */
  1.995562716158404E-5,                /* Expression: sps.A12
                                        * Referenced by: '<S105>/A12'
                                        */
  -0.49238617409376229,                /* Expression: sps.A21
                                        * Referenced by: '<S105>/A21'
                                        */
  0.99556271615840386,                 /* Expression: sps.A22
                                        * Referenced by: '<S105>/A22'
                                        */
  9.97781358079202E-6,                 /* Expression: sps.B11
                                        * Referenced by: '<S106>/B11'
                                        */
  0.99778135807920187,                 /* Expression: sps.B21
                                        * Referenced by: '<S106>/B21'
                                        */
  2.4619308704688116E-6,               /* Expression: sps.D
                                        * Referenced by: '<S104>/D*u(k)'
                                        */
  0.49347981688184189,                 /* Expression: sps.C11
                                        * Referenced by: '<S107>/C11'
                                        */
  4.9238617409376232E-6,               /* Expression: sps.C12
                                        * Referenced by: '<S107>/C12'
                                        */
  0.97484015911852984,                 /* Expression: sps.A11
                                        * Referenced by: '<S109>/A11'
                                        */
  1.576518862980693E-5,                /* Expression: sps.A12
                                        * Referenced by: '<S109>/A12'
                                        */
  -2515.9840881470318,                 /* Expression: sps.A21
                                        * Referenced by: '<S109>/A21'
                                        */
  0.57651886298069654,                 /* Expression: sps.A22
                                        * Referenced by: '<S109>/A22'
                                        */
  7.8825943149034344E-6,               /* Expression: sps.B11
                                        * Referenced by: '<S110>/B11'
                                        */
  0.78825943149034827,                 /* Expression: sps.B21
                                        * Referenced by: '<S110>/B21'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S76>/Constant1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S69>/Constant2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S62>/Unit Delay3'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S69>/Gain1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S62>/Unit Delay1'
                                        */
  1.0E-5,                              /* Computed Parameter: Integrator_gainval_e
                                        * Referenced by: '<S142>/Integrator'
                                        */
  0.9226,                              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S142>/Integrator'
                                        */

  /*  Computed Parameter: RTEPeriodMeter_P1_Size
   * Referenced by: '<S12>/RTE Period Meter'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: Period
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */

  /*  Computed Parameter: RTEPeriodMeter_P2_Size
   * Referenced by: '<S12>/RTE Period Meter'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: Frequency
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */

  /*  Computed Parameter: RTEPeriodMeter_P3_Size
   * Referenced by: '<S12>/RTE Period Meter'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: EdgeType
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */

  /*  Computed Parameter: RTEPeriodMeter_P4_Size
   * Referenced by: '<S12>/RTE Period Meter'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: TOn
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */

  /*  Computed Parameter: RTEPeriodMeter_P5_Size
   * Referenced by: '<S12>/RTE Period Meter'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: TOff
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */

  /*  Computed Parameter: RTEPeriodMeter_P6_Size
   * Referenced by: '<S12>/RTE Period Meter'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DutyCycle
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */

  /*  Computed Parameter: RTEPeriodMeter_P7_Size
   * Referenced by: '<S12>/RTE Period Meter'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: DuringType
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */

  /*  Computed Parameter: RTEPeriodMeter_P8_Size
   * Referenced by: '<S12>/RTE Period Meter'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: MinimumFrequency
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */
  10.0,                                /* Expression: 10
                                        * Referenced by: '<S12>/Saturation'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S12>/Saturation'
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
   * Referenced by: '<S40>/default'
   */
  { 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U },
  336592896U,                          /* Computed Parameter: addr4_Value
                                        * Referenced by: '<S40>/addr4'
                                        */
  269484032U,                          /* Computed Parameter: addr3_Value
                                        * Referenced by: '<S40>/addr3'
                                        */
  202375168U,                          /* Computed Parameter: addr2_Value
                                        * Referenced by: '<S40>/addr2'
                                        */
  135266304U,                          /* Computed Parameter: addr1_Value
                                        * Referenced by: '<S40>/addr1'
                                        */
  68157440U,                           /* Computed Parameter: addr_Value
                                        * Referenced by: '<S40>/addr'
                                        */
  101U,                                /* Computed Parameter: sat_scn_UpperSat
                                        * Referenced by: '<S10>/sat_scn'
                                        */
  0U,                                  /* Computed Parameter: sat_scn_LowerSat
                                        * Referenced by: '<S10>/sat_scn'
                                        */
  374865680U,                          /* Computed Parameter: blockID_Value
                                        * Referenced by: '<S40>/blockID'
                                        */
  2U,                                  /* Computed Parameter: IOTypeSel_Value
                                        * Referenced by: '<S13>/IOTypeSel'
                                        */
  2U,                                  /* Computed Parameter: IOTypeSel_Value_n
                                        * Referenced by: '<S14>/IOTypeSel'
                                        */
  255U,                                /* Computed Parameter: Saturation_UpperSat_k
                                        * Referenced by: '<S10>/Saturation'
                                        */
  1U,                                  /* Computed Parameter: Saturation_LowerSat_f
                                        * Referenced by: '<S10>/Saturation'
                                        */

  /*  Computed Parameter: load_config1_Value
   * Referenced by: '<S10>/load_config1'
   */
  { 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U },
  2147483648U,                         /* Computed Parameter: shift_2bits_Gain
                                        * Referenced by: '<S10>/shift_2bits'
                                        */
  0,                                   /* Expression: ~strcmp(FormatType,'Unsigned')
                                        * Referenced by: '<S45>/Constant'
                                        */
  1,                                   /* Expression: ~strcmp(FormatType,'Unsigned')
                                        * Referenced by: '<S43>/Constant'
                                        */
  1,                                   /* Expression: ~strcmp(FormatType,'Unsigned')
                                        * Referenced by: '<S44>/Constant'
                                        */

  /* Start of '<S7>/TrueRMS ' */
  {
    60.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S35>/Gain'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S35>/integrator'
                                        */
    0.016666666666666666,              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S35>/Transport Delay'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S35>/Transport Delay'
                                        */
    0.016666666666666666,              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S35>/K1'
                                        */
    14400.0,                           /* Expression: sps.Vinit
                                        * Referenced by: '<S35>/Memory'
                                        */
    0.0,                               /* Expression: Inf
                                        * Referenced by: '<S28>/Saturation to avoid negative sqrt'
                                        */
    0.0                                /* Expression: 0
                                        * Referenced by: '<S28>/Saturation to avoid negative sqrt'
                                        */
  }
  /* End of '<S7>/TrueRMS ' */
  ,

  /* Start of '<S7>/RMS ' */
  {
    60.0,                              /* Mask Parameter: Fourier1_Freq
                                        * Referenced by:
                                        *   '<S29>/cos(wt)'
                                        *   '<S29>/sin(wt)'
                                        */
    60.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S32>/Gain'
                                        */
    60.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S33>/Gain'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S33>/integrator'
                                        */
    0.016666666666666666,              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S33>/Transport Delay'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S33>/Transport Delay'
                                        */
    0.016666666666666666,              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S33>/K1'
                                        */
    169.70562748477141,                /* Expression: sps.Vinit
                                        * Referenced by: '<S33>/Memory'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S32>/integrator'
                                        */
    0.016666666666666666,              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S32>/Transport Delay'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S32>/Transport Delay'
                                        */
    0.016666666666666666,              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S32>/K1'
                                        */
    0.0,                               /* Expression: sps.Vinit
                                        * Referenced by: '<S32>/Memory'
                                        */
    2.0,                               /* Expression: sps.k
                                        * Referenced by: '<S29>/sin(wt)'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S29>/sin(wt)'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S29>/sin(wt)'
                                        */
    2.0,                               /* Expression: sps.k
                                        * Referenced by: '<S29>/cos(wt)'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S29>/cos(wt)'
                                        */
    1.5707963267948966,                /* Expression: pi/2
                                        * Referenced by: '<S29>/cos(wt)'
                                        */
    57.295779513082323,                /* Expression: 180/pi
                                        * Referenced by: '<S29>/Rad->Deg.'
                                        */
    0.70710678118654746                /* Expression: 1/sqrt(2)
                                        * Referenced by: '<S27>/Gain'
                                        */
  }
  /* End of '<S7>/RMS ' */
  ,

  /* Start of '<S6>/TrueRMS ' */
  {
    60.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S26>/Gain'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S26>/integrator'
                                        */
    0.016666666666666666,              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S26>/Transport Delay'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S26>/Transport Delay'
                                        */
    0.016666666666666666,              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S26>/K1'
                                        */
    14400.0,                           /* Expression: sps.Vinit
                                        * Referenced by: '<S26>/Memory'
                                        */
    0.0,                               /* Expression: Inf
                                        * Referenced by: '<S19>/Saturation to avoid negative sqrt'
                                        */
    0.0                                /* Expression: 0
                                        * Referenced by: '<S19>/Saturation to avoid negative sqrt'
                                        */
  }
  /* End of '<S6>/TrueRMS ' */
  ,

  /* Start of '<S6>/RMS ' */
  {
    60.0,                              /* Mask Parameter: Fourier1_Freq
                                        * Referenced by:
                                        *   '<S20>/cos(wt)'
                                        *   '<S20>/sin(wt)'
                                        */
    60.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S23>/Gain'
                                        */
    60.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S24>/Gain'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S24>/integrator'
                                        */
    0.016666666666666666,              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S24>/Transport Delay'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S24>/Transport Delay'
                                        */
    0.016666666666666666,              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S24>/K1'
                                        */
    169.70562748477141,                /* Expression: sps.Vinit
                                        * Referenced by: '<S24>/Memory'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S23>/integrator'
                                        */
    0.016666666666666666,              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S23>/Transport Delay'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S23>/Transport Delay'
                                        */
    0.016666666666666666,              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S23>/K1'
                                        */
    0.0,                               /* Expression: sps.Vinit
                                        * Referenced by: '<S23>/Memory'
                                        */
    2.0,                               /* Expression: sps.k
                                        * Referenced by: '<S20>/sin(wt)'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S20>/sin(wt)'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S20>/sin(wt)'
                                        */
    2.0,                               /* Expression: sps.k
                                        * Referenced by: '<S20>/cos(wt)'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S20>/cos(wt)'
                                        */
    1.5707963267948966,                /* Expression: pi/2
                                        * Referenced by: '<S20>/cos(wt)'
                                        */
    57.295779513082323,                /* Expression: 180/pi
                                        * Referenced by: '<S20>/Rad->Deg.'
                                        */
    0.70710678118654746                /* Expression: 1/sqrt(2)
                                        * Referenced by: '<S18>/Gain'
                                        */
  }
  /* End of '<S6>/RMS ' */
};
