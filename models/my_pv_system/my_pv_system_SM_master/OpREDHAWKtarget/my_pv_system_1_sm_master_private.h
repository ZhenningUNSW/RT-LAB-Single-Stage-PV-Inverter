/*
 * my_pv_system_1_sm_master_private.h
 *
 * Code generation for model "my_pv_system_1_sm_master".
 *
 * Model version              : 1.207
 * Simulink Coder version : 8.7 (R2014b) 08-Sep-2014
 * C source code generated on : Wed Sep 07 16:46:13 2016
 *
 * Target selection: rtlab_rtmodel.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#ifndef RTW_HEADER_my_pv_system_1_sm_master_private_h_
#define RTW_HEADER_my_pv_system_1_sm_master_private_h_
#include "rtwtypes.h"
#include "builtin_typeid_types.h"
#include "multiword_types.h"
#include <math.h>
#include <stdlib.h>
#ifndef CodeFormat
#define CodeFormat                     S-Function
#else
#undef CodeFormat
#define CodeFormat                     S-Function
#endif

#ifndef S_FUNCTION_NAME
#define S_FUNCTION_NAME                simulink_only_sfcn
#else
#undef S_FUNCTION_NAME
#define S_FUNCTION_NAME                simulink_only_sfcn
#endif

#ifndef S_FUNCTION_LEVEL
#define S_FUNCTION_LEVEL               2
#else
#undef S_FUNCTION_LEVEL
#define S_FUNCTION_LEVEL               2
#endif

#ifndef RTW_GENERATED_S_FUNCTION
#define RTW_GENERATED_S_FUNCTION
#endif

#ifndef rtmGetDataMapInfo
# define rtmGetDataMapInfo(rtm)        NULL
#endif

#ifndef rtmSetDataMapInfo
# define rtmSetDataMapInfo(rtm, val)
#endif

#if !defined(RTW_SFUNCTION_DEFINES)
#define RTW_SFUNCTION_DEFINES
#ifndef _RTW_COMMON_DEFINES_
#define _RTW_COMMON_DEFINES_
#endif
#endif

real_T rt_TDelayInterpolate(
  real_T tMinusDelay,                  /* tMinusDelay = currentSimTime - delay */
  real_T tStart,
  real_T *tBuf,
  real_T *uBuf,
  int_T bufSz,
  int_T *lastIdx,
  int_T oldestIdx,
  int_T newIdx,
  real_T initOutput,
  boolean_T discrete,
  boolean_T minorStepAndTAtLastMajorOutput)
  ;
extern void send_rt(SimStruct *rts);
extern void OP_SEND(SimStruct *rts);
extern void opmonitor(SimStruct *rts);
extern void recv_rt(SimStruct *rts);
extern void RECV_Param(SimStruct *rts);

/* private model entry point functions */
extern void my_pv_system_1_sm_master_derivatives(void);

#endif                                 /* RTW_HEADER_my_pv_system_1_sm_master_private_h_ */
