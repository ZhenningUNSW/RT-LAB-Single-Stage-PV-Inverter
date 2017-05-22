/*
 * boost_and_two_level__1_sm_ehs_data.c
 *
 * Code generation for model "boost_and_two_level__1_sm_ehs".
 *
 * Model version              : 1.1155
 * Simulink Coder version : 8.7 (R2014b) 08-Sep-2014
 * C source code generated on : Mon May 22 20:50:22 2017
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
                                        *   '<S73>/Constant4'
                                        *   '<S83>/Gain'
                                        *   '<S86>/Gain'
                                        *   '<S105>/Gain'
                                        *   '<S99>/Gain'
                                        *   '<S102>/Gain'
                                        */
  1.0,                                 /* Mask Parameter: PLL_AGC
                                        * Referenced by: '<S91>/Constant1'
                                        */
  2.0,                                 /* Mask Parameter: AlphaBetaZerotodq0_Alignment
                                        * Referenced by: '<S116>/Constant'
                                        */
  50.0,                                /* Mask Parameter: InverterControl_Fnom
                                        * Referenced by:
                                        *   '<S73>/Constant4'
                                        *   '<S80>/First cycle of simulation Id=0.92, Iq=0'
                                        */
  0.003,                               /* Mask Parameter: InverterControl_Increment_MPPT
                                        * Referenced by: '<S69>/Iph_3'
                                        */
  314.15926535897933,                  /* Mask Parameter: Discrete_Init
                                        * Referenced by: '<S93>/Discrete-Time Integrator'
                                        */
  1.0,                                 /* Mask Parameter: Discrete_Kd
                                        * Referenced by: '<S93>/Discrete Derivative '
                                        */
  6.6,                                 /* Mask Parameter: InverterControl_Ki_Ireg
                                        * Referenced by: '<S75>/Integral Gain'
                                        */
  200.0,                               /* Mask Parameter: InverterControl_Ki_VDCreg
                                        * Referenced by: '<S146>/Integral Gain'
                                        */
  180.0,                               /* Mask Parameter: Discrete_Kp
                                        * Referenced by: '<S93>/Kp4'
                                        */
  0.15,                                /* Mask Parameter: InverterControl_Kp_Ireg
                                        * Referenced by: '<S75>/Proportional Gain'
                                        */
  12.0,                                /* Mask Parameter: InverterControl_Kp_VDCreg
                                        * Referenced by: '<S146>/Proportional Gain'
                                        */
  -1.5,                                /* Mask Parameter: PI_LowerSaturationLimit
                                        * Referenced by: '<S75>/Saturate'
                                        */
  -1.5,                                /* Mask Parameter: PI_LowerSaturationLimit_j
                                        * Referenced by: '<S146>/Saturate'
                                        */

  /*  Mask Parameter: PWM_Generator_MinMax
   * Referenced by: '<S72>/Constant10'
   */
  { -1.0, 1.0 },
  3500.0,                              /* Mask Parameter: InverterControl_Pnom
                                        * Referenced by: '<S71>/A->pu'
                                        */
  1.0,                                 /* Mask Parameter: RMS_TrueRMS
                                        * Referenced by: '<S7>/Constant'
                                        */
  1.0,                                 /* Mask Parameter: RMS1_TrueRMS
                                        * Referenced by: '<S8>/Constant'
                                        */
  1.0,                                 /* Mask Parameter: RMS1_TrueRMS_l
                                        * Referenced by: '<S67>/Constant'
                                        */
  1.5,                                 /* Mask Parameter: PI_UpperSaturationLimit
                                        * Referenced by: '<S75>/Saturate'
                                        */
  1.5,                                 /* Mask Parameter: PI_UpperSaturationLimit_a
                                        * Referenced by: '<S146>/Saturate'
                                        */
  425.0,                               /* Mask Parameter: InverterControl_Vdc_ref_Init
                                        * Referenced by:
                                        *   '<S65>/Vnom_dc1'
                                        *   '<S69>/Iph_'
                                        */
  400.0,                               /* Mask Parameter: InverterControl_Vnom_dc
                                        * Referenced by: '<S74>/Rtot_pu2'
                                        */
  240.0,                               /* Mask Parameter: InverterControl_Vnom_prim
                                        * Referenced by:
                                        *   '<S71>/A->pu'
                                        *   '<S71>/V->pu'
                                        *   '<S73>/Constant3'
                                        */
  1.0,                                 /* Mask Parameter: CompareToConstant_const
                                        * Referenced by: '<S118>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: CompareToConstant1_const
                                        * Referenced by: '<S119>/Constant'
                                        */
  1.0,                                 /* Mask Parameter: CompareToConstant_const_d
                                        * Referenced by: '<S4>/Constant'
                                        */
  5.0,                                 /* Mask Parameter: Counter_max_count
                                        * Referenced by: '<S45>/Switch'
                                        */
  127U,                                /* Mask Parameter: BitwiseOperator_BitMask
                                        * Referenced by: '<S48>/Bitwise Operator'
                                        */
  65535U,                              /* Mask Parameter: BitwiseOperator_BitMask_j
                                        * Referenced by: '<S46>/Bitwise Operator'
                                        */
  65535U,                              /* Mask Parameter: BitwiseOperator_BitMask_h
                                        * Referenced by: '<S47>/Bitwise Operator'
                                        */
  3U,                                  /* Mask Parameter: DetectChange_vinit
                                        * Referenced by: '<S51>/Delay Input1'
                                        */
  0U,                                  /* Mask Parameter: DetectChange_vinit_i
                                        * Referenced by: '<S50>/Delay Input1'
                                        */
  0U,                                  /* Mask Parameter: DetectChange_vinit_g
                                        * Referenced by: '<S54>/Delay Input1'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S41>/Gain'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S42>/Gain'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S45>/Constant1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S45>/Constant2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S45>/Constant'
                                        */
  31.99951171875,                      /* Expression: 2^(N-Q)-2^(-Q)
                                        * Referenced by: '<S46>/Saturation1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S46>/Saturation1'
                                        */
  2048.0,                              /* Expression: 2^Q
                                        * Referenced by: '<S46>/Gain2'
                                        */
  15.99951171875,                      /* Expression: 2^(N-Q-1)-2^(-Q)
                                        * Referenced by: '<S46>/Saturation'
                                        */
  -16.0,                               /* Expression: -2^(N-Q-1)
                                        * Referenced by: '<S46>/Saturation'
                                        */
  2048.0,                              /* Expression: 2^Q
                                        * Referenced by: '<S46>/Gain1'
                                        */
  31.99951171875,                      /* Expression: 2^(N-Q)-2^(-Q)
                                        * Referenced by: '<S47>/Saturation1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S47>/Saturation1'
                                        */
  2048.0,                              /* Expression: 2^Q
                                        * Referenced by: '<S47>/Gain2'
                                        */
  15.99951171875,                      /* Expression: 2^(N-Q-1)-2^(-Q)
                                        * Referenced by: '<S47>/Saturation'
                                        */
  -16.0,                               /* Expression: -2^(N-Q-1)
                                        * Referenced by: '<S47>/Saturation'
                                        */
  2048.0,                              /* Expression: 2^Q
                                        * Referenced by: '<S47>/Gain1'
                                        */
  127.0,                               /* Expression: 2^(N-Q)-2^(-Q)
                                        * Referenced by: '<S48>/Saturation1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S48>/Saturation1'
                                        */
  1.0,                                 /* Expression: 2^Q
                                        * Referenced by: '<S48>/Gain2'
                                        */
  63.0,                                /* Expression: 2^(N-Q-1)-2^(-Q)
                                        * Referenced by: '<S48>/Saturation'
                                        */
  -64.0,                               /* Expression: -2^(N-Q-1)
                                        * Referenced by: '<S48>/Saturation'
                                        */
  1.0,                                 /* Expression: 2^Q
                                        * Referenced by: '<S48>/Gain1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S13>/Constant5'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S13>/Constant'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S84>/Gain1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S87>/Gain1'
                                        */
  -1.1368683772160957E-16,             /* Expression: sps.K2
                                        * Referenced by: '<S89>/Gain1'
                                        */
  6.4623485355705288E-30,              /* Expression: sps.K1
                                        * Referenced by: '<S89>/Gain'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S100>/Gain1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S103>/Gain1'
                                        */
  1.0,                                 /* Expression: [1]
                                        * Referenced by: '<S92>/Gain'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S96>/Gain1'
                                        */
  1.0E-5,                              /* Computed Parameter: Integ4_gainval
                                        * Referenced by: '<S99>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S99>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S99>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S99>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size
   * Referenced by: '<S101>/S-Function'
   */
  { 1.0, 1.0 },
  0.02502,                             /* Expression: MaxDelay
                                        * Referenced by: '<S101>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size
   * Referenced by: '<S101>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S101>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size
   * Referenced by: '<S101>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S101>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size
   * Referenced by: '<S101>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S101>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S100>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S99>/Constant'
                                        */
  1.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S99>/Unit Delay1'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S96>/Gain3'
                                        */
  1.0E-5,                              /* Computed Parameter: Integ4_gainval_c
                                        * Referenced by: '<S102>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S102>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S102>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S102>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_h
   * Referenced by: '<S104>/S-Function'
   */
  { 1.0, 1.0 },
  0.02502,                             /* Expression: MaxDelay
                                        * Referenced by: '<S104>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_p
   * Referenced by: '<S104>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S104>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_p
   * Referenced by: '<S104>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S104>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_p
   * Referenced by: '<S104>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S104>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S103>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S102>/Constant'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S102>/Unit Delay1'
                                        */
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<S96>/Rad->Deg.'
                                        */
  0.0,                                 /* Expression: inf
                                        * Referenced by: '<S92>/Saturation'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S92>/Saturation'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S106>/Gain1'
                                        */

  /*  Expression: [0.92 0]
   * Referenced by: '<S80>/Constant'
   */
  { 0.92, 0.0 },

  /*  Expression: [0,0]
   * Referenced by: '<S120>/dq'
   */
  { 0.0, 0.0 },

  /*  Expression: [0,0]
   * Referenced by: '<S121>/dq'
   */
  { 0.0, 0.0 },
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/S-Function1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S2>/S-Function'
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

  /*  Computed Parameter: SFunction_P1_Size_a
   * Referenced by: '<S60>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: Acqu_group
                                        * Referenced by: '<S60>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_p
   * Referenced by: '<S18>/S-Function'
   */
  { 1.0, 1.0 },
  165.0,                               /* Expression: Data_width
                                        * Referenced by: '<S18>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_d
   * Referenced by: '<S18>/S-Function'
   */
  { 165.0, 1.0 },

  /*  Expression: InitialConditions
   * Referenced by: '<S18>/S-Function'
   */
  { 1000.0, 25.0, 1.0, 0.0, 49.0, 48.0, 52.0, 53.0, 53.0, 9.0, 1.0, 3.0, 12.0,
    13.0, 0.0, 2.0, 48.0, 49.0, 50.0, 51.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0,
    54.0, 54.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 0.0066, 0.33,
    0.066, 0.004125, 0.0048529, 1.0, 0.005, 0.005, 0.1, 0.1, 0.1, 0.1, 0.005,
    0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.65, 1.65, 1.65, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.3, 3.3, 3.3, 3.3, 3.3, 5.0, 5.0, 5.0, 5.0,
    5.0, 5.0, 5.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0,
    16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -16.0, -16.0, -16.0,
    -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0,
    -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, 2.0 },
  0.5,                                 /* Expression: .5
                                        * Referenced by: '<S11>/Relay1'
                                        */
  0.5,                                 /* Expression: .5
                                        * Referenced by: '<S11>/Relay1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S11>/Relay1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S11>/Relay1'
                                        */

  /*  Computed Parameter: eHS_rst_loadin_P1_Size
   * Referenced by: '<S11>/eHS_rst_loadin'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: eHS_rst_loadin_P1
   * Referenced by: '<S11>/eHS_rst_loadin'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: eHS_rst_loadin_P2_Size
   * Referenced by: '<S11>/eHS_rst_loadin'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: FcnNos
                                        * Referenced by: '<S11>/eHS_rst_loadin'
                                        */

  /*  Computed Parameter: eHS_rst_loadin_P3_Size
   * Referenced by: '<S11>/eHS_rst_loadin'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: width
                                        * Referenced by: '<S11>/eHS_rst_loadin'
                                        */

  /*  Computed Parameter: eHS_rst_loadin_P4_Size
   * Referenced by: '<S11>/eHS_rst_loadin'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S11>/eHS_rst_loadin'
                                        */

  /*  Computed Parameter: Automated_Solver_Mat_Initialisa
   * Referenced by: '<S11>/Automated_Solver_Mat_Initialisation_1'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: Automated_Solver_Mat_Initiali_c
   * Referenced by: '<S11>/Automated_Solver_Mat_Initialisation_1'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: Automated_Solver_Mat_Initiali_o
   * Referenced by: '<S11>/Automated_Solver_Mat_Initialisation_1'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: fpga_port_in
                                        * Referenced by: '<S11>/Automated_Solver_Mat_Initialisation_1'
                                        */

  /*  Computed Parameter: Automated_Solver_Mat_Initiali_e
   * Referenced by: '<S11>/Automated_Solver_Mat_Initialisation_1'
   */
  { 1.0, 31.0 },

  /*  Computed Parameter: Automated_Solver_Mat_Initial_er
   * Referenced by: '<S11>/Automated_Solver_Mat_Initialisation_1'
   */
  { 101.0, 104.0, 115.0, 95.0, 99.0, 111.0, 110.0, 102.0, 105.0, 103.0, 95.0,
    101.0, 72.0, 83.0, 120.0, 54.0, 52.0, 118.0, 51.0, 95.0, 83.0, 73.0, 68.0,
    52.0, 56.0, 57.0, 52.0, 46.0, 109.0, 97.0, 116.0 },
  20364.675298172569,                  /* Expression: 14400*sqrt(2)
                                        * Referenced by: '<S13>/Sine Wave Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S13>/Sine Wave Function'
                                        */
  314.15926535897933,                  /* Expression: 50*2*pi
                                        * Referenced by: '<S13>/Sine Wave Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S13>/Sine Wave Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S13>/Memory1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S13>/Memory'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P1_Size
   * Referenced by: '<S11>/Inputs_eHS1_Send'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: Inputs_eHS1_Send_P1
   * Referenced by: '<S11>/Inputs_eHS1_Send'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: Inputs_eHS1_Send_P2_Size
   * Referenced by: '<S11>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: FcnNos
                                        * Referenced by: '<S11>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P3_Size
   * Referenced by: '<S11>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  3.0,                                 /* Expression: width
                                        * Referenced by: '<S11>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P4_Size
   * Referenced by: '<S11>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S11>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P5_Size
   * Referenced by: '<S11>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: checkVersion
                                        * Referenced by: '<S11>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P6_Size
   * Referenced by: '<S11>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedId
                                        * Referenced by: '<S11>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P7_Size
   * Referenced by: '<S11>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedVersion
                                        * Referenced by: '<S11>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Inputs_eHS1_Send_P8_Size
   * Referenced by: '<S11>/Inputs_eHS1_Send'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: opComp
                                        * Referenced by: '<S11>/Inputs_eHS1_Send'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P1_Size
   * Referenced by: '<S11>/Outputs_eHS1_Recv'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: Outputs_eHS1_Recv_P1
   * Referenced by: '<S11>/Outputs_eHS1_Recv'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: Outputs_eHS1_Recv_P2_Size
   * Referenced by: '<S11>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: FcnNos
                                        * Referenced by: '<S11>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P3_Size
   * Referenced by: '<S11>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  8.0,                                 /* Expression: width
                                        * Referenced by: '<S11>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P4_Size
   * Referenced by: '<S11>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S11>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P5_Size
   * Referenced by: '<S11>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: sampleTime
                                        * Referenced by: '<S11>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P6_Size
   * Referenced by: '<S11>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: checkVersion
                                        * Referenced by: '<S11>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P7_Size
   * Referenced by: '<S11>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedId
                                        * Referenced by: '<S11>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P8_Size
   * Referenced by: '<S11>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedVersion
                                        * Referenced by: '<S11>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P9_Size
   * Referenced by: '<S11>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: opComp
                                        * Referenced by: '<S11>/Outputs_eHS1_Recv'
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

  /*  Computed Parameter: SFunction_P1_Size_pu
   * Referenced by: '<S61>/S-Function'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: Acqu_group
                                        * Referenced by: '<S61>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S3>/Memory2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S3>/Memory2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S3>/Memory2'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_b
   * Referenced by: '<S62>/S-Function'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: Acqu_group
                                        * Referenced by: '<S62>/S-Function'
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

  /*  Computed Parameter: SFunction_P1_Size_j
   * Referenced by: '<S63>/S-Function'
   */
  { 1.0, 1.0 },
  3.0,                                 /* Expression: Acqu_group
                                        * Referenced by: '<S63>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S3>/Memory3'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S3>/Memory3'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_ad
   * Referenced by: '<S64>/S-Function'
   */
  { 1.0, 1.0 },
  5.0,                                 /* Expression: Acqu_group
                                        * Referenced by: '<S64>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S45>/Memory1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S45>/Memory2'
                                        */

  /*  Computed Parameter: LoadIn_P1_Size
   * Referenced by: '<S43>/LoadIn'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: LoadIn_P1
   * Referenced by: '<S43>/LoadIn'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: LoadIn_P2_Size
   * Referenced by: '<S43>/LoadIn'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: FcnNos
                                        * Referenced by: '<S43>/LoadIn'
                                        */

  /*  Computed Parameter: LoadIn_P3_Size
   * Referenced by: '<S43>/LoadIn'
   */
  { 1.0, 1.0 },
  18.0,                                /* Expression: width
                                        * Referenced by: '<S43>/LoadIn'
                                        */

  /*  Computed Parameter: LoadIn_P4_Size
   * Referenced by: '<S43>/LoadIn'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S43>/LoadIn'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P1_Size_e
   * Referenced by: '<S3>/Outputs_eHS1_Recv'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: Outputs_eHS1_Recv_P1_p
   * Referenced by: '<S3>/Outputs_eHS1_Recv'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: Outputs_eHS1_Recv_P2_Size_g
   * Referenced by: '<S3>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: FcnNos
                                        * Referenced by: '<S3>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P3_Size_a
   * Referenced by: '<S3>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  7.0,                                 /* Expression: width
                                        * Referenced by: '<S3>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P4_Size_b
   * Referenced by: '<S3>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S3>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P5_Size_m
   * Referenced by: '<S3>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: sampleTime
                                        * Referenced by: '<S3>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P6_Size_b
   * Referenced by: '<S3>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: checkVersion
                                        * Referenced by: '<S3>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P7_Size_j
   * Referenced by: '<S3>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedId
                                        * Referenced by: '<S3>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P8_Size_p
   * Referenced by: '<S3>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedVersion
                                        * Referenced by: '<S3>/Outputs_eHS1_Recv'
                                        */

  /*  Computed Parameter: Outputs_eHS1_Recv_P9_Size_m
   * Referenced by: '<S3>/Outputs_eHS1_Recv'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: opComp
                                        * Referenced by: '<S3>/Outputs_eHS1_Recv'
                                        */

  /*  Expression: ones(1,9)
   * Referenced by: '<S44>/Constant'
   */
  { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 },

  /*  Computed Parameter: DataInSend_P1_Size
   * Referenced by: '<S43>/DataIn Send'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: DataInSend_P1
   * Referenced by: '<S43>/DataIn Send'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: DataInSend_P2_Size
   * Referenced by: '<S43>/DataIn Send'
   */
  { 1.0, 1.0 },
  15.0,                                /* Expression: FcnNos
                                        * Referenced by: '<S43>/DataIn Send'
                                        */

  /*  Computed Parameter: DataInSend_P3_Size
   * Referenced by: '<S43>/DataIn Send'
   */
  { 1.0, 1.0 },
  16.0,                                /* Expression: width
                                        * Referenced by: '<S43>/DataIn Send'
                                        */

  /*  Computed Parameter: DataInSend_P4_Size
   * Referenced by: '<S43>/DataIn Send'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S43>/DataIn Send'
                                        */

  /*  Computed Parameter: DataInSend_P5_Size
   * Referenced by: '<S43>/DataIn Send'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: checkVersion
                                        * Referenced by: '<S43>/DataIn Send'
                                        */

  /*  Computed Parameter: DataInSend_P6_Size
   * Referenced by: '<S43>/DataIn Send'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedId
                                        * Referenced by: '<S43>/DataIn Send'
                                        */

  /*  Computed Parameter: DataInSend_P7_Size
   * Referenced by: '<S43>/DataIn Send'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: expectedVersion
                                        * Referenced by: '<S43>/DataIn Send'
                                        */

  /*  Computed Parameter: DataInSend_P8_Size
   * Referenced by: '<S43>/DataIn Send'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: opComp
                                        * Referenced by: '<S43>/DataIn Send'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S5>/Constant'
                                        */

  /*  Computed Parameter: rtlab_io_block_P1_Size
   * Referenced by: '<S15>/rtlab_io_block'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: rtlab_io_block_P1
   * Referenced by: '<S15>/rtlab_io_block'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: rtlab_io_block_P2_Size
   * Referenced by: '<S15>/rtlab_io_block'
   */
  { 0.0, 0.0 },

  /*  Computed Parameter: rtlab_io_block_P3_Size
   * Referenced by: '<S15>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: sampleTime
                                        * Referenced by: '<S15>/rtlab_io_block'
                                        */

  /*  Computed Parameter: rtlab_io_block_P4_Size
   * Referenced by: '<S15>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  5.0,                                 /* Expression: portNb
                                        * Referenced by: '<S15>/rtlab_io_block'
                                        */

  /*  Computed Parameter: rtlab_io_block_P5_Size
   * Referenced by: '<S15>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  32.0,                                /* Expression: size
                                        * Referenced by: '<S15>/rtlab_io_block'
                                        */

  /*  Computed Parameter: rtlab_io_block_P6_Size
   * Referenced by: '<S15>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: nbChannels
                                        * Referenced by: '<S15>/rtlab_io_block'
                                        */

  /*  Computed Parameter: rtlab_io_block_P7_Size
   * Referenced by: '<S15>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: numwidth
                                        * Referenced by: '<S15>/rtlab_io_block'
                                        */

  /*  Computed Parameter: rtlab_io_block_P8_Size
   * Referenced by: '<S15>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: timeUnit
                                        * Referenced by: '<S15>/rtlab_io_block'
                                        */

  /*  Computed Parameter: rtlab_io_block_P9_Size
   * Referenced by: '<S15>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  8.0,                                 /* Expression: maxcount
                                        * Referenced by: '<S15>/rtlab_io_block'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S15>/IOTypeSel1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S15>/Memory1'
                                        */

  /*  Computed Parameter: IOTypeSel_LoadIn_P1_Size
   * Referenced by: '<S15>/IOTypeSel_LoadIn'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: IOTypeSel_LoadIn_P1
   * Referenced by: '<S15>/IOTypeSel_LoadIn'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: IOTypeSel_LoadIn_P2_Size
   * Referenced by: '<S15>/IOTypeSel_LoadIn'
   */
  { 1.0, 1.0 },
  13.0,                                /* Expression: FcnNos
                                        * Referenced by: '<S15>/IOTypeSel_LoadIn'
                                        */

  /*  Computed Parameter: IOTypeSel_LoadIn_P3_Size
   * Referenced by: '<S15>/IOTypeSel_LoadIn'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: width
                                        * Referenced by: '<S15>/IOTypeSel_LoadIn'
                                        */

  /*  Computed Parameter: IOTypeSel_LoadIn_P4_Size
   * Referenced by: '<S15>/IOTypeSel_LoadIn'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S15>/IOTypeSel_LoadIn'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S15>/Memory'
                                        */

  /*  Computed Parameter: rtlab_io_block_P1_Size_b
   * Referenced by: '<S16>/rtlab_io_block'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: rtlab_io_block_P1_m
   * Referenced by: '<S16>/rtlab_io_block'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: rtlab_io_block_P2_Size_k
   * Referenced by: '<S16>/rtlab_io_block'
   */
  { 0.0, 0.0 },

  /*  Computed Parameter: rtlab_io_block_P3_Size_i
   * Referenced by: '<S16>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: sampleTime
                                        * Referenced by: '<S16>/rtlab_io_block'
                                        */

  /*  Computed Parameter: rtlab_io_block_P4_Size_g
   * Referenced by: '<S16>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  6.0,                                 /* Expression: portNb
                                        * Referenced by: '<S16>/rtlab_io_block'
                                        */

  /*  Computed Parameter: rtlab_io_block_P5_Size_j
   * Referenced by: '<S16>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  32.0,                                /* Expression: size
                                        * Referenced by: '<S16>/rtlab_io_block'
                                        */

  /*  Computed Parameter: rtlab_io_block_P6_Size_k
   * Referenced by: '<S16>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  8.0,                                 /* Expression: nbChannels
                                        * Referenced by: '<S16>/rtlab_io_block'
                                        */

  /*  Computed Parameter: rtlab_io_block_P7_Size_a
   * Referenced by: '<S16>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: numwidth
                                        * Referenced by: '<S16>/rtlab_io_block'
                                        */

  /*  Computed Parameter: rtlab_io_block_P8_Size_b
   * Referenced by: '<S16>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: timeUnit
                                        * Referenced by: '<S16>/rtlab_io_block'
                                        */

  /*  Computed Parameter: rtlab_io_block_P9_Size_i
   * Referenced by: '<S16>/rtlab_io_block'
   */
  { 1.0, 1.0 },
  8.0,                                 /* Expression: maxcount
                                        * Referenced by: '<S16>/rtlab_io_block'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S16>/IOTypeSel1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S16>/Memory1'
                                        */

  /*  Computed Parameter: IOTypeSel_LoadIn_P1_Size_e
   * Referenced by: '<S16>/IOTypeSel_LoadIn'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: IOTypeSel_LoadIn_P1_l
   * Referenced by: '<S16>/IOTypeSel_LoadIn'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: IOTypeSel_LoadIn_P2_Size_o
   * Referenced by: '<S16>/IOTypeSel_LoadIn'
   */
  { 1.0, 1.0 },
  14.0,                                /* Expression: FcnNos
                                        * Referenced by: '<S16>/IOTypeSel_LoadIn'
                                        */

  /*  Computed Parameter: IOTypeSel_LoadIn_P3_Size_m
   * Referenced by: '<S16>/IOTypeSel_LoadIn'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: width
                                        * Referenced by: '<S16>/IOTypeSel_LoadIn'
                                        */

  /*  Computed Parameter: IOTypeSel_LoadIn_P4_Size_e
   * Referenced by: '<S16>/IOTypeSel_LoadIn'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portType
                                        * Referenced by: '<S16>/IOTypeSel_LoadIn'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S16>/Memory'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S41>/integrator'
                                        */
  0.02,                                /* Expression: 1./sps.Freq
                                        * Referenced by: '<S41>/Transport Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S41>/Transport Delay'
                                        */
  0.02,                                /* Expression: 1./sps.Freq
                                        * Referenced by: '<S41>/K1'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S41>/Memory'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S42>/integrator'
                                        */
  0.02,                                /* Expression: 1./sps.Freq
                                        * Referenced by: '<S42>/Transport Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S42>/Transport Delay'
                                        */
  0.02,                                /* Expression: 1./sps.Freq
                                        * Referenced by: '<S42>/K1'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S42>/Memory'
                                        */

  /*  Computed Parameter: OpWriteFile_P1_Size
   * Referenced by: '<S3>/OpWriteFile'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: Acq_Group
                                        * Referenced by: '<S3>/OpWriteFile'
                                        */

  /*  Computed Parameter: OpWriteFile_P2_Size
   * Referenced by: '<S3>/OpWriteFile'
   */
  { 1.0, 17.0 },

  /*  Computed Parameter: OpWriteFile_P2
   * Referenced by: '<S3>/OpWriteFile'
   */
  { 109.0, 121.0, 80.0, 86.0, 95.0, 105.0, 110.0, 118.0, 101.0, 114.0, 116.0,
    101.0, 114.0, 46.0, 109.0, 97.0, 116.0 },

  /*  Computed Parameter: OpWriteFile_P3_Size
   * Referenced by: '<S3>/OpWriteFile'
   */
  { 1.0, 5.0 },

  /*  Computed Parameter: OpWriteFile_P3
   * Referenced by: '<S3>/OpWriteFile'
   */
  { 111.0, 112.0, 118.0, 97.0, 114.0 },

  /*  Computed Parameter: OpWriteFile_P4_Size
   * Referenced by: '<S3>/OpWriteFile'
   */
  { 1.0, 1.0 },
  100.0,                               /* Expression: Decimation
                                        * Referenced by: '<S3>/OpWriteFile'
                                        */

  /*  Computed Parameter: OpWriteFile_P5_Size
   * Referenced by: '<S3>/OpWriteFile'
   */
  { 1.0, 1.0 },
  400000.0,                            /* Expression: Nb_Samples
                                        * Referenced by: '<S3>/OpWriteFile'
                                        */

  /*  Computed Parameter: OpWriteFile_P6_Size
   * Referenced by: '<S3>/OpWriteFile'
   */
  { 1.0, 1.0 },
  1.048576E+6,                         /* Expression: Buffer_size
                                        * Referenced by: '<S3>/OpWriteFile'
                                        */

  /*  Computed Parameter: OpWriteFile_P7_Size
   * Referenced by: '<S3>/OpWriteFile'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: Sim_Mode
                                        * Referenced by: '<S3>/OpWriteFile'
                                        */

  /*  Computed Parameter: OpWriteFile_P8_Size
   * Referenced by: '<S3>/OpWriteFile'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: Static_File
                                        * Referenced by: '<S3>/OpWriteFile'
                                        */

  /*  Computed Parameter: OpWriteFile_P9_Size
   * Referenced by: '<S3>/OpWriteFile'
   */
  { 1.0, 1.0 },
  1.073741824E+9,                      /* Expression: file_size
                                        * Referenced by: '<S3>/OpWriteFile'
                                        */

  /*  Computed Parameter: OpWriteFile_P10_Size
   * Referenced by: '<S3>/OpWriteFile'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: write_offline
                                        * Referenced by: '<S3>/OpWriteFile'
                                        */

  /*  Computed Parameter: OpWriteFile_P11_Size
   * Referenced by: '<S3>/OpWriteFile'
   */
  { 1.0, 2.0 },

  /*  Computed Parameter: OpWriteFile_P11
   * Referenced by: '<S3>/OpWriteFile'
   */
  { 46.0, 47.0 },

  /*  Computed Parameter: OpWriteFile_P12_Size
   * Referenced by: '<S3>/OpWriteFile'
   */
  { 1.0, 2.0 },

  /*  Computed Parameter: OpWriteFile_P12
   * Referenced by: '<S3>/OpWriteFile'
   */
  { 46.0, 92.0 },
  0.05,                                /* Expression: 0.05
                                        * Referenced by: '<S13>/Step'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S13>/Step'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S13>/Step'
                                        */
  -1.0,                                /* Expression: -1
                                        * Referenced by: '<S13>/Gain'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S13>/Switch1'
                                        */

  /*  Computed Parameter: RTEConversion_P1_Size
   * Referenced by: '<S13>/RTE Conversion'
   */
  { 1.0, 1.0 },
  100.0,                               /* Expression: nbMaxEvents
                                        * Referenced by: '<S13>/RTE Conversion'
                                        */

  /*  Computed Parameter: RTEConversion_P2_Size
   * Referenced by: '<S13>/RTE Conversion'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: inputdatatype
                                        * Referenced by: '<S13>/RTE Conversion'
                                        */

  /*  Computed Parameter: RTEConversion_P3_Size
   * Referenced by: '<S13>/RTE Conversion'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: outputdatatype
                                        * Referenced by: '<S13>/RTE Conversion'
                                        */

  /*  Computed Parameter: RTEConversion_P4_Size
   * Referenced by: '<S13>/RTE Conversion'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: compensation
                                        * Referenced by: '<S13>/RTE Conversion'
                                        */

  /*  Computed Parameter: RTEConversion_P5_Size
   * Referenced by: '<S13>/RTE Conversion'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: sampleTime
                                        * Referenced by: '<S13>/RTE Conversion'
                                        */
  0.1684,                              /* Expression: 0.1684
                                        * Referenced by: '<S65>/Unit Delay'
                                        */
  3.125E-5,                            /* Expression: sps.Delay
                                        * Referenced by: '<S144>/Constant3'
                                        */
  0.000125,                            /* Expression: sps.Period
                                        * Referenced by: '<S144>/Constant1'
                                        */
  8000.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S144>/1\ib1'
                                        */

  /*  Expression: [0 .5 1]
   * Referenced by: '<S144>/Lookup Table'
   */
  { 0.0, 0.5, 1.0 },

  /*  Expression: [0 2 0]
   * Referenced by: '<S144>/Lookup Table'
   */
  { 0.0, 2.0, 0.0 },
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S144>/Constant2'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S122>/Gain1'
                                        */

  /*  Computed Parameter: RTEConversion1_P1_Size
   * Referenced by: '<S13>/RTE Conversion1'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: nbMaxEvents
                                        * Referenced by: '<S13>/RTE Conversion1'
                                        */

  /*  Computed Parameter: RTEConversion1_P2_Size
   * Referenced by: '<S13>/RTE Conversion1'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: inputdatatype
                                        * Referenced by: '<S13>/RTE Conversion1'
                                        */

  /*  Computed Parameter: RTEConversion1_P3_Size
   * Referenced by: '<S13>/RTE Conversion1'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: outputdatatype
                                        * Referenced by: '<S13>/RTE Conversion1'
                                        */

  /*  Computed Parameter: RTEConversion1_P4_Size
   * Referenced by: '<S13>/RTE Conversion1'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: compensation
                                        * Referenced by: '<S13>/RTE Conversion1'
                                        */

  /*  Computed Parameter: RTEConversion1_P5_Size
   * Referenced by: '<S13>/RTE Conversion1'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: sampleTime
                                        * Referenced by: '<S13>/RTE Conversion1'
                                        */

  /*  Computed Parameter: RTEConversion2_P1_Size
   * Referenced by: '<S13>/RTE Conversion2'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: nbMaxEvents
                                        * Referenced by: '<S13>/RTE Conversion2'
                                        */

  /*  Computed Parameter: RTEConversion2_P2_Size
   * Referenced by: '<S13>/RTE Conversion2'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: inputdatatype
                                        * Referenced by: '<S13>/RTE Conversion2'
                                        */

  /*  Computed Parameter: RTEConversion2_P3_Size
   * Referenced by: '<S13>/RTE Conversion2'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: outputdatatype
                                        * Referenced by: '<S13>/RTE Conversion2'
                                        */

  /*  Computed Parameter: RTEConversion2_P4_Size
   * Referenced by: '<S13>/RTE Conversion2'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: compensation
                                        * Referenced by: '<S13>/RTE Conversion2'
                                        */

  /*  Computed Parameter: RTEConversion2_P5_Size
   * Referenced by: '<S13>/RTE Conversion2'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: sampleTime
                                        * Referenced by: '<S13>/RTE Conversion2'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S13>/Switch'
                                        */

  /*  Computed Parameter: RTELogicalOperator1_P1_Size
   * Referenced by: '<S13>/RTE Logical Operator1'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: LogicalOperator
                                        * Referenced by: '<S13>/RTE Logical Operator1'
                                        */

  /*  Computed Parameter: RTELogicalOperator1_P2_Size
   * Referenced by: '<S13>/RTE Logical Operator1'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: NbrInput
                                        * Referenced by: '<S13>/RTE Logical Operator1'
                                        */

  /*  Computed Parameter: RTELogicalOperator1_P3_Size
   * Referenced by: '<S13>/RTE Logical Operator1'
   */
  { 1.0, 1.0 },
  255.0,                               /* Expression: NbrMaxEvents
                                        * Referenced by: '<S13>/RTE Logical Operator1'
                                        */

  /*  Computed Parameter: RTEGround_P1_Size
   * Referenced by: '<S13>/RTE Ground'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: SampleTime
                                        * Referenced by: '<S13>/RTE Ground'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P1_Size
   * Referenced by: '<S52>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: conversionType
                                        * Referenced by: '<S52>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P2_Size
   * Referenced by: '<S52>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  8.0,                                 /* Expression: nbLine
                                        * Referenced by: '<S52>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P3_Size
   * Referenced by: '<S52>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: nbEvents
                                        * Referenced by: '<S52>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P4_Size
   * Referenced by: '<S52>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: timeUnit
                                        * Referenced by: '<S52>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P5_Size
   * Referenced by: '<S52>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: enTimeFactor
                                        * Referenced by: '<S52>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P6_Size
   * Referenced by: '<S52>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: sampleTime
                                        * Referenced by: '<S52>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: RTE_Conversion_1_P7_Size
   * Referenced by: '<S52>/RTE_Conversion_1'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: initStates
                                        * Referenced by: '<S52>/RTE_Conversion_1'
                                        */

  /*  Computed Parameter: EventGen_eHS_1_P1_Size
   * Referenced by: '<S52>/EventGen_eHS_1'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: EventGen_eHS_1_P1
   * Referenced by: '<S52>/EventGen_eHS_1'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: EventGen_eHS_1_P2_Size
   * Referenced by: '<S52>/EventGen_eHS_1'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: portNb
                                        * Referenced by: '<S52>/EventGen_eHS_1'
                                        */

  /*  Computed Parameter: EventGen_eHS_1_P3_Size
   * Referenced by: '<S52>/EventGen_eHS_1'
   */
  { 0.0, 0.0 },

  /*  Computed Parameter: EventGen_eHS_1_P4_Size
   * Referenced by: '<S52>/EventGen_eHS_1'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: size
                                        * Referenced by: '<S52>/EventGen_eHS_1'
                                        */

  /*  Computed Parameter: EventGen_eHS_1_P5_Size
   * Referenced by: '<S52>/EventGen_eHS_1'
   */
  { 1.0, 1.0 },
  8.0,                                 /* Expression: nbChannels
                                        * Referenced by: '<S52>/EventGen_eHS_1'
                                        */

  /*  Computed Parameter: EventGen_eHS_1_P6_Size
   * Referenced by: '<S52>/EventGen_eHS_1'
   */
  { 1.0, 1.0 },
  4.0,                                 /* Expression: numwidth
                                        * Referenced by: '<S52>/EventGen_eHS_1'
                                        */

  /*  Computed Parameter: EventGen_eHS_1_P7_Size
   * Referenced by: '<S52>/EventGen_eHS_1'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: timeUnit
                                        * Referenced by: '<S52>/EventGen_eHS_1'
                                        */
  8.55,                                /* Expression: 8.55
                                        * Referenced by: '<S13>/Constant1'
                                        */
  523.6,                               /* Expression: 37.4*14
                                        * Referenced by: '<S13>/Constant2'
                                        */
  0.06,                                /* Expression: 0.06
                                        * Referenced by: '<S13>/Constant3'
                                        */
  840.0,                               /* Expression: 60*14
                                        * Referenced by: '<S13>/Constant4'
                                        */
  50.0,                                /* Expression: sps.Finit
                                        * Referenced by: '<S91>/Unit Delay'
                                        */
  70.0,                                /* Expression: 70
                                        * Referenced by: '<S80>/avoid division by zero'
                                        */
  40.0,                                /* Expression: 40
                                        * Referenced by: '<S80>/avoid division by zero'
                                        */
  0.25,                                /* Expression: 1/4
                                        * Referenced by: '<S80>/Gain'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_i
   * Referenced by: '<S117>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: MaxDelay
                                        * Referenced by: '<S117>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_i
   * Referenced by: '<S117>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S117>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_i
   * Referenced by: '<S117>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S117>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_n
   * Referenced by: '<S117>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S117>/S-Function'
                                        */
  2.0E-5,                              /* Computed Parameter: DiscreteTimeIntegrator_gainval
                                        * Referenced by: '<S91>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: sps.Phase_Init*pi/180
                                        * Referenced by: '<S91>/Discrete-Time Integrator'
                                        */
  6.2831853071795862,                  /* Expression: 2*pi
                                        * Referenced by: '<S91>/Constant4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S80>/First cycle of simulation Id=0.92, Iq=0'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S80>/First cycle of simulation Id=0.92, Iq=0'
                                        */
  0.01257992044073508,                 /* Expression: sps.D
                                        * Referenced by: '<S112>/D*u(k)'
                                        */

  /*  Expression: sps.x0(1,:)
   * Referenced by: '<S112>/Delay_x1'
   */
  { 0.00028197057830512594, 0.0 },
  3151.6695002188571,                  /* Expression: sps.C11
                                        * Referenced by: '<S115>/C11'
                                        */

  /*  Expression: sps.x0(2,:)
   * Referenced by: '<S112>/Delay_x2'
   */
  { -5.5511151231257827E-17, 0.0 },
  0.02515984088147016,                 /* Expression: sps.C12
                                        * Referenced by: '<S115>/C12'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S65>/Unit Delay2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S65>/Iq_ref'
                                        */
  2.0E-5,                              /* Computed Parameter: Integrator_gainval
                                        * Referenced by: '<S75>/Integrator'
                                        */
  0.0,                                 /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S75>/Integrator'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S76>/Gain1'
                                        */
  1.0E-5,                              /* Computed Parameter: Integ4_gainval_l
                                        * Referenced by: '<S83>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S83>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S83>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S83>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_l
   * Referenced by: '<S85>/S-Function'
   */
  { 1.0, 1.0 },
  0.022242222222222222,                /* Expression: MaxDelay
                                        * Referenced by: '<S85>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_c
   * Referenced by: '<S85>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S85>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_b
   * Referenced by: '<S85>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S85>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_f
   * Referenced by: '<S85>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S85>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S84>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S83>/Constant'
                                        */
  1.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S83>/Unit Delay1'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S76>/Gain3'
                                        */
  1.0E-5,                              /* Computed Parameter: Integ4_gainval_g
                                        * Referenced by: '<S86>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S86>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S86>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S86>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_hi
   * Referenced by: '<S88>/S-Function'
   */
  { 1.0, 1.0 },
  0.022242222222222222,                /* Expression: MaxDelay
                                        * Referenced by: '<S88>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_j
   * Referenced by: '<S88>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S88>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_c
   * Referenced by: '<S88>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S88>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_m
   * Referenced by: '<S88>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S88>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S87>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S86>/Constant'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S86>/Unit Delay1'
                                        */
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<S76>/Rad->Deg.'
                                        */
  0.017453292519943295,                /* Expression: pi/180
                                        * Referenced by: '<S71>/to-rad'
                                        */
  0.002,                               /* Expression: RLff(1)
                                        * Referenced by: '<S68>/Rff '
                                        */
  0.2,                                 /* Expression: RLff(2)
                                        * Referenced by: '<S68>/Lff  '
                                        */
  0.002,                               /* Expression: RLff(1)
                                        * Referenced by: '<S68>/Rff'
                                        */
  0.2,                                 /* Expression: RLff(2)
                                        * Referenced by: '<S68>/Lff'
                                        */
  1.5,                                 /* Expression: 1.5
                                        * Referenced by: '<S68>/Saturation'
                                        */
  -1.5,                                /* Expression: -1.5
                                        * Referenced by: '<S68>/Saturation'
                                        */
  450.0,                               /* Expression: Limits_MPPT(1)
                                        * Referenced by: '<S69>/Iph_1'
                                        */
  375.0,                               /* Expression: Limits_MPPT(2)
                                        * Referenced by: '<S69>/Iph_2'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S13>/MPPT_On'
                                        */
  1.0E-5,                              /* Computed Parameter: Integ4_gainval_o
                                        * Referenced by: '<S89>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S89>/Integ4'
                                        */
  0.02,                                /* Expression: sps.Delay
                                        * Referenced by: '<S89>/K1'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_o
   * Referenced by: '<S90>/S-Function'
   */
  { 1.0, 1.0 },
  0.02004,                             /* Expression: MaxDelay
                                        * Referenced by: '<S90>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_id
   * Referenced by: '<S90>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S90>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_l
   * Referenced by: '<S90>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S90>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_i
   * Referenced by: '<S90>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S90>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S89>/Unit Delay'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S89>/K2'
                                        */
  425.0,                               /* Expression: sps.Vinit
                                        * Referenced by: '<S89>/Unit Delay1'
                                        */
  1.0E-5,                              /* Computed Parameter: Integ4_gainval_b
                                        * Referenced by: '<S105>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S105>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S105>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S105>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_lc
   * Referenced by: '<S107>/S-Function'
   */
  { 1.0, 1.0 },
  0.02502,                             /* Expression: MaxDelay
                                        * Referenced by: '<S107>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_g
   * Referenced by: '<S107>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S107>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_f
   * Referenced by: '<S107>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S107>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_a
   * Referenced by: '<S107>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S107>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S106>/Unit Delay'
                                        */
  0.02,                                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S105>/Constant'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S105>/Unit Delay1'
                                        */

  /*  Expression: [ TcD  Ts-TcD ]
   * Referenced by: '<S93>/Discrete Derivative '
   */
  { 0.0001, -8.0E-5 },
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S93>/Discrete Derivative '
                                        */
  0.064,                               /* Computed Parameter: DiscreteTimeIntegrator_gainva_a
                                        * Referenced by: '<S93>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: Par_Limits(1)
                                        * Referenced by: '<S93>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: Par_Limits(2)
                                        * Referenced by: '<S93>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: Par_Limits(1)
                                        * Referenced by: '<S93>/Saturation1'
                                        */
  0.0,                                 /* Expression: Par_Limits(2)
                                        * Referenced by: '<S93>/Saturation1'
                                        */
  0.15915494309189535,                 /* Expression: 1/2/pi
                                        * Referenced by: '<S91>/Gain10'
                                        */
  0.00024000000000000003,              /* Computed Parameter: RateLimiter_RisingLim
                                        * Referenced by: '<S91>/Rate Limiter'
                                        */
  -0.00024000000000000003,             /* Computed Parameter: RateLimiter_FallingLim
                                        * Referenced by: '<S91>/Rate Limiter'
                                        */
  50.0,                                /* Expression: sps.Finit
                                        * Referenced by: '<S91>/Rate Limiter'
                                        */
  101.32101697571224,                  /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S108>/Delay_x1'
                                        */
  0.999995076138259,                   /* Expression: sps.A11
                                        * Referenced by: '<S109>/A11'
                                        */
  0.0,                                 /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S108>/Delay_x2'
                                        */
  1.995562716158404E-5,                /* Expression: sps.A12
                                        * Referenced by: '<S109>/A12'
                                        */
  -0.49238617409376229,                /* Expression: sps.A21
                                        * Referenced by: '<S109>/A21'
                                        */
  0.99556271615840386,                 /* Expression: sps.A22
                                        * Referenced by: '<S109>/A22'
                                        */
  9.97781358079202E-6,                 /* Expression: sps.B11
                                        * Referenced by: '<S110>/B11'
                                        */
  0.99778135807920187,                 /* Expression: sps.B21
                                        * Referenced by: '<S110>/B21'
                                        */
  2.4619308704688116E-6,               /* Expression: sps.D
                                        * Referenced by: '<S108>/D*u(k)'
                                        */
  0.49347981688184189,                 /* Expression: sps.C11
                                        * Referenced by: '<S111>/C11'
                                        */
  4.9238617409376232E-6,               /* Expression: sps.C12
                                        * Referenced by: '<S111>/C12'
                                        */
  0.97484015911852984,                 /* Expression: sps.A11
                                        * Referenced by: '<S113>/A11'
                                        */
  1.576518862980693E-5,                /* Expression: sps.A12
                                        * Referenced by: '<S113>/A12'
                                        */
  -2515.9840881470318,                 /* Expression: sps.A21
                                        * Referenced by: '<S113>/A21'
                                        */
  0.57651886298069654,                 /* Expression: sps.A22
                                        * Referenced by: '<S113>/A22'
                                        */
  7.8825943149034344E-6,               /* Expression: sps.B11
                                        * Referenced by: '<S114>/B11'
                                        */
  0.78825943149034827,                 /* Expression: sps.B21
                                        * Referenced by: '<S114>/B21'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S80>/Constant1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S73>/Constant2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S65>/Unit Delay3'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S73>/Gain1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S65>/Unit Delay1'
                                        */
  1.0E-5,                              /* Computed Parameter: Integrator_gainval_e
                                        * Referenced by: '<S146>/Integrator'
                                        */
  0.9226,                              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S146>/Integrator'
                                        */

  /*  Computed Parameter: RTEPeriodMeter_P1_Size
   * Referenced by: '<S13>/RTE Period Meter'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: Period
                                        * Referenced by: '<S13>/RTE Period Meter'
                                        */

  /*  Computed Parameter: RTEPeriodMeter_P2_Size
   * Referenced by: '<S13>/RTE Period Meter'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: Frequency
                                        * Referenced by: '<S13>/RTE Period Meter'
                                        */

  /*  Computed Parameter: RTEPeriodMeter_P3_Size
   * Referenced by: '<S13>/RTE Period Meter'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: EdgeType
                                        * Referenced by: '<S13>/RTE Period Meter'
                                        */

  /*  Computed Parameter: RTEPeriodMeter_P4_Size
   * Referenced by: '<S13>/RTE Period Meter'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: TOn
                                        * Referenced by: '<S13>/RTE Period Meter'
                                        */

  /*  Computed Parameter: RTEPeriodMeter_P5_Size
   * Referenced by: '<S13>/RTE Period Meter'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: TOff
                                        * Referenced by: '<S13>/RTE Period Meter'
                                        */

  /*  Computed Parameter: RTEPeriodMeter_P6_Size
   * Referenced by: '<S13>/RTE Period Meter'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DutyCycle
                                        * Referenced by: '<S13>/RTE Period Meter'
                                        */

  /*  Computed Parameter: RTEPeriodMeter_P7_Size
   * Referenced by: '<S13>/RTE Period Meter'
   */
  { 1.0, 1.0 },
  2.0,                                 /* Expression: DuringType
                                        * Referenced by: '<S13>/RTE Period Meter'
                                        */

  /*  Computed Parameter: RTEPeriodMeter_P8_Size
   * Referenced by: '<S13>/RTE Period Meter'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: MinimumFrequency
                                        * Referenced by: '<S13>/RTE Period Meter'
                                        */
  10.0,                                /* Expression: 10
                                        * Referenced by: '<S13>/Saturation'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S13>/Saturation'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S3>/Constant1'
                                        */

  /*  Computed Parameter: OpTrigger_P1_Size
   * Referenced by: '<S3>/OpTrigger'
   */
  { 1.0, 1.0 },
  26.0,                                /* Expression: Acq_Group
                                        * Referenced by: '<S3>/OpTrigger'
                                        */

  /*  Computed Parameter: OpTrigger_P2_Size
   * Referenced by: '<S3>/OpTrigger'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: Trig_Type
                                        * Referenced by: '<S3>/OpTrigger'
                                        */

  /*  Computed Parameter: OpTrigger_P3_Size
   * Referenced by: '<S3>/OpTrigger'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: Trig_Offset
                                        * Referenced by: '<S3>/OpTrigger'
                                        */

  /*  Computed Parameter: OpCtrl_P1_Size
   * Referenced by: '<S3>/OpCtrl'
   */
  { 1.0, 6.0 },

  /*  Computed Parameter: OpCtrl_P1
   * Referenced by: '<S3>/OpCtrl'
   */
  { 79.0, 112.0, 67.0, 116.0, 114.0, 108.0 },

  /*  Computed Parameter: OpCtrl_P2_Size
   * Referenced by: '<S3>/OpCtrl'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: boardid
                                        * Referenced by: '<S3>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P3_Size
   * Referenced by: '<S3>/OpCtrl'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: mode
                                        * Referenced by: '<S3>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P4_Size
   * Referenced by: '<S3>/OpCtrl'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: externalClock
                                        * Referenced by: '<S3>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P5_Size
   * Referenced by: '<S3>/OpCtrl'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: decimRtsi
                                        * Referenced by: '<S3>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P6_Size
   * Referenced by: '<S3>/OpCtrl'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S3>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P7_Size
   * Referenced by: '<S3>/OpCtrl'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: SampleTime
                                        * Referenced by: '<S3>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P8_Size
   * Referenced by: '<S3>/OpCtrl'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: calibIO
                                        * Referenced by: '<S3>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P9_Size
   * Referenced by: '<S3>/OpCtrl'
   */
  { 1.0, 1.0 },
  -1.0,                                /* Expression: numconfig
                                        * Referenced by: '<S3>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P10_Size
   * Referenced by: '<S3>/OpCtrl'
   */
  { 1.0, 1.0 },
  -1.0,                                /* Expression: loadinport
                                        * Referenced by: '<S3>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P11_Size
   * Referenced by: '<S3>/OpCtrl'
   */
  { 1.0, 1.0 },
  6.0,                                 /* Expression: BoardType
                                        * Referenced by: '<S3>/OpCtrl'
                                        */

  /*  Computed Parameter: OpCtrl_P12_Size
   * Referenced by: '<S3>/OpCtrl'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: sync_type
                                        * Referenced by: '<S3>/OpCtrl'
                                        */

  /*  Computed Parameter: default_Value
   * Referenced by: '<S43>/default'
   */
  { 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U },
  336592896U,                          /* Computed Parameter: addr4_Value
                                        * Referenced by: '<S43>/addr4'
                                        */
  269484032U,                          /* Computed Parameter: addr3_Value
                                        * Referenced by: '<S43>/addr3'
                                        */
  202375168U,                          /* Computed Parameter: addr2_Value
                                        * Referenced by: '<S43>/addr2'
                                        */
  135266304U,                          /* Computed Parameter: addr1_Value
                                        * Referenced by: '<S43>/addr1'
                                        */
  68157440U,                           /* Computed Parameter: addr_Value
                                        * Referenced by: '<S43>/addr'
                                        */
  101U,                                /* Computed Parameter: sat_scn_UpperSat
                                        * Referenced by: '<S11>/sat_scn'
                                        */
  0U,                                  /* Computed Parameter: sat_scn_LowerSat
                                        * Referenced by: '<S11>/sat_scn'
                                        */
  374865680U,                          /* Computed Parameter: blockID_Value
                                        * Referenced by: '<S43>/blockID'
                                        */
  1U,                                  /* Computed Parameter: IOTypeSel_Value
                                        * Referenced by: '<S15>/IOTypeSel'
                                        */
  1U,                                  /* Computed Parameter: IOTypeSel_Value_n
                                        * Referenced by: '<S16>/IOTypeSel'
                                        */
  255U,                                /* Computed Parameter: Saturation_UpperSat_k
                                        * Referenced by: '<S11>/Saturation'
                                        */
  1U,                                  /* Computed Parameter: Saturation_LowerSat_f
                                        * Referenced by: '<S11>/Saturation'
                                        */

  /*  Computed Parameter: load_config1_Value
   * Referenced by: '<S11>/load_config1'
   */
  { 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U },
  2147483648U,                         /* Computed Parameter: shift_2bits_Gain
                                        * Referenced by: '<S11>/shift_2bits'
                                        */
  0,                                   /* Computed Parameter: Q_Y0
                                        * Referenced by: '<S17>/Q'
                                        */
  1,                                   /* Computed Parameter: Q_Y0_i
                                        * Referenced by: '<S17>/!Q'
                                        */
  0,                                   /* Expression: ~strcmp(FormatType,'Unsigned')
                                        * Referenced by: '<S48>/Constant'
                                        */
  1,                                   /* Expression: ~strcmp(FormatType,'Unsigned')
                                        * Referenced by: '<S46>/Constant'
                                        */
  1,                                   /* Expression: ~strcmp(FormatType,'Unsigned')
                                        * Referenced by: '<S47>/Constant'
                                        */

  /* Start of '<S67>/TrueRMS ' */
  {
    50.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S155>/Gain'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S155>/integrator'
                                        */
    0.02,                              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S155>/Transport Delay'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S155>/Transport Delay'
                                        */
    0.02,                              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S155>/K1'
                                        */
    0.0,                               /* Expression: sps.Vinit
                                        * Referenced by: '<S155>/Memory'
                                        */
    0.0,                               /* Expression: Inf
                                        * Referenced by: '<S148>/Saturation to avoid negative sqrt'
                                        */
    0.0                                /* Expression: 0
                                        * Referenced by: '<S148>/Saturation to avoid negative sqrt'
                                        */
  }
  /* End of '<S67>/TrueRMS ' */
  ,

  /* Start of '<S67>/RMS ' */
  {
    50.0,                              /* Mask Parameter: Fourier1_Freq
                                        * Referenced by:
                                        *   '<S149>/cos(wt)'
                                        *   '<S149>/sin(wt)'
                                        */
    50.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S152>/Gain'
                                        */
    50.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S153>/Gain'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S153>/integrator'
                                        */
    0.02,                              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S153>/Transport Delay'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S153>/Transport Delay'
                                        */
    0.02,                              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S153>/K1'
                                        */
    0.0,                               /* Expression: sps.Vinit
                                        * Referenced by: '<S153>/Memory'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S152>/integrator'
                                        */
    0.02,                              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S152>/Transport Delay'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S152>/Transport Delay'
                                        */
    0.02,                              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S152>/K1'
                                        */
    0.0,                               /* Expression: sps.Vinit
                                        * Referenced by: '<S152>/Memory'
                                        */
    2.0,                               /* Expression: sps.k
                                        * Referenced by: '<S149>/sin(wt)'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S149>/sin(wt)'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S149>/sin(wt)'
                                        */
    2.0,                               /* Expression: sps.k
                                        * Referenced by: '<S149>/cos(wt)'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S149>/cos(wt)'
                                        */
    1.5707963267948966,                /* Expression: pi/2
                                        * Referenced by: '<S149>/cos(wt)'
                                        */
    57.295779513082323,                /* Expression: 180/pi
                                        * Referenced by: '<S149>/Rad->Deg.'
                                        */
    0.70710678118654746                /* Expression: 1/sqrt(2)
                                        * Referenced by: '<S147>/Gain'
                                        */
  }
  /* End of '<S67>/RMS ' */
  ,

  /* Start of '<S8>/TrueRMS ' */
  {
    50.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S38>/Gain'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S38>/integrator'
                                        */
    0.02,                              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S38>/Transport Delay'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S38>/Transport Delay'
                                        */
    0.02,                              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S38>/K1'
                                        */
    0.0,                               /* Expression: sps.Vinit
                                        * Referenced by: '<S38>/Memory'
                                        */
    0.0,                               /* Expression: Inf
                                        * Referenced by: '<S31>/Saturation to avoid negative sqrt'
                                        */
    0.0                                /* Expression: 0
                                        * Referenced by: '<S31>/Saturation to avoid negative sqrt'
                                        */
  }
  /* End of '<S8>/TrueRMS ' */
  ,

  /* Start of '<S8>/RMS ' */
  {
    50.0,                              /* Mask Parameter: Fourier1_Freq
                                        * Referenced by:
                                        *   '<S32>/cos(wt)'
                                        *   '<S32>/sin(wt)'
                                        */
    50.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S35>/Gain'
                                        */
    50.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S36>/Gain'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S36>/integrator'
                                        */
    0.02,                              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S36>/Transport Delay'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S36>/Transport Delay'
                                        */
    0.02,                              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S36>/K1'
                                        */
    0.0,                               /* Expression: sps.Vinit
                                        * Referenced by: '<S36>/Memory'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S35>/integrator'
                                        */
    0.02,                              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S35>/Transport Delay'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S35>/Transport Delay'
                                        */
    0.02,                              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S35>/K1'
                                        */
    0.0,                               /* Expression: sps.Vinit
                                        * Referenced by: '<S35>/Memory'
                                        */
    2.0,                               /* Expression: sps.k
                                        * Referenced by: '<S32>/sin(wt)'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S32>/sin(wt)'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S32>/sin(wt)'
                                        */
    2.0,                               /* Expression: sps.k
                                        * Referenced by: '<S32>/cos(wt)'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S32>/cos(wt)'
                                        */
    1.5707963267948966,                /* Expression: pi/2
                                        * Referenced by: '<S32>/cos(wt)'
                                        */
    57.295779513082323,                /* Expression: 180/pi
                                        * Referenced by: '<S32>/Rad->Deg.'
                                        */
    0.70710678118654746                /* Expression: 1/sqrt(2)
                                        * Referenced by: '<S30>/Gain'
                                        */
  }
  /* End of '<S8>/RMS ' */
  ,

  /* Start of '<S7>/TrueRMS ' */
  {
    50.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S29>/Gain'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S29>/integrator'
                                        */
    0.02,                              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S29>/Transport Delay'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S29>/Transport Delay'
                                        */
    0.02,                              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S29>/K1'
                                        */
    0.0,                               /* Expression: sps.Vinit
                                        * Referenced by: '<S29>/Memory'
                                        */
    0.0,                               /* Expression: Inf
                                        * Referenced by: '<S22>/Saturation to avoid negative sqrt'
                                        */
    0.0                                /* Expression: 0
                                        * Referenced by: '<S22>/Saturation to avoid negative sqrt'
                                        */
  }
  /* End of '<S7>/TrueRMS ' */
  ,

  /* Start of '<S7>/RMS ' */
  {
    50.0,                              /* Mask Parameter: Fourier1_Freq
                                        * Referenced by:
                                        *   '<S23>/cos(wt)'
                                        *   '<S23>/sin(wt)'
                                        */
    50.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S26>/Gain'
                                        */
    50.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S27>/Gain'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S27>/integrator'
                                        */
    0.02,                              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S27>/Transport Delay'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S27>/Transport Delay'
                                        */
    0.02,                              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S27>/K1'
                                        */
    0.0,                               /* Expression: sps.Vinit
                                        * Referenced by: '<S27>/Memory'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S26>/integrator'
                                        */
    0.02,                              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S26>/Transport Delay'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S26>/Transport Delay'
                                        */
    0.02,                              /* Expression: 1./sps.Freq
                                        * Referenced by: '<S26>/K1'
                                        */
    0.0,                               /* Expression: sps.Vinit
                                        * Referenced by: '<S26>/Memory'
                                        */
    2.0,                               /* Expression: sps.k
                                        * Referenced by: '<S23>/sin(wt)'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S23>/sin(wt)'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S23>/sin(wt)'
                                        */
    2.0,                               /* Expression: sps.k
                                        * Referenced by: '<S23>/cos(wt)'
                                        */
    0.0,                               /* Expression: 0
                                        * Referenced by: '<S23>/cos(wt)'
                                        */
    1.5707963267948966,                /* Expression: pi/2
                                        * Referenced by: '<S23>/cos(wt)'
                                        */
    57.295779513082323,                /* Expression: 180/pi
                                        * Referenced by: '<S23>/Rad->Deg.'
                                        */
    0.70710678118654746                /* Expression: 1/sqrt(2)
                                        * Referenced by: '<S21>/Gain'
                                        */
  }
  /* End of '<S7>/RMS ' */
};
