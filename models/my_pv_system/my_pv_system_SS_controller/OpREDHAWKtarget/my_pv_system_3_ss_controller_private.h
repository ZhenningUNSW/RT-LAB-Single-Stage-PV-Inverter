/*
 * my_pv_system_3_ss_controller_private.h
 *
 * Code generation for model "my_pv_system_3_ss_controller".
 *
 * Model version              : 1.189
 * Simulink Coder version : 8.7 (R2014b) 08-Sep-2014
 * C source code generated on : Tue Sep 06 11:32:06 2016
 *
 * Target selection: rtlab_rtmodel.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#ifndef RTW_HEADER_my_pv_system_3_ss_controller_private_h_
#define RTW_HEADER_my_pv_system_3_ss_controller_private_h_
#include "rtwtypes.h"
#include "builtin_typeid_types.h"
#include "multiword_types.h"

extern real_T rt_roundd_snf(real_T u);
extern real_T rt_modd_snf(real_T u0, real_T u1);
extern real_T rt_hypotd_snf(real_T u0, real_T u1);
extern real_T rt_atan2d_snf(real_T u0, real_T u1);
extern real_T rt_remd_snf(real_T u0, real_T u1);
extern real_T rt_powd_snf(real_T u0, real_T u1);
void BINARYSEARCH_real_T(uint32_T *piLeft, uint32_T *piRght, real_T u, const
  real_T *pData, uint32_T iHi);
void LookUp_real_T_real_T(real_T *pY, const real_T *pYData, real_T u, const
  real_T *pUData, uint32_T iHi);
extern void sfun_discreteVariableDelay(SimStruct *rts);
extern void send_rt(SimStruct *rts);
extern void opmonitor(SimStruct *rts);
extern void OP_SEND(SimStruct *rts);
extern void recv_rt(SimStruct *rts);

#endif                                 /* RTW_HEADER_my_pv_system_3_ss_controller_private_h_ */
