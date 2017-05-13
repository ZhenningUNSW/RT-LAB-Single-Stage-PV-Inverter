/*
 * boost_and_two_level__1_sm_ehs_private.h
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
#ifndef RTW_HEADER_boost_and_two_level__1_sm_ehs_private_h_
#define RTW_HEADER_boost_and_two_level__1_sm_ehs_private_h_
#include "rtwtypes.h"
#include "builtin_typeid_types.h"
#include "multiword_types.h"
#ifndef UCHAR_MAX
#include <limits.h>
#endif

#if ( UCHAR_MAX != (0xFFU) ) || ( SCHAR_MAX != (0x7F) )
#error Code was generated for compiler with different sized uchar/char. \
Consider adjusting Test hardware word size settings on the \
Hardware Implementation pane to match your compiler word sizes as \
defined in limits.h of the compiler. Alternatively, you can \
select the Test hardware is the same as production hardware option and \
select the Enable portable word sizes option on the Code Generation > \
Verification pane for ERT based targets, which will disable the \
preprocessor word size checks.
#endif

#if ( USHRT_MAX != (0xFFFFU) ) || ( SHRT_MAX != (0x7FFF) )
#error Code was generated for compiler with different sized ushort/short. \
Consider adjusting Test hardware word size settings on the \
Hardware Implementation pane to match your compiler word sizes as \
defined in limits.h of the compiler. Alternatively, you can \
select the Test hardware is the same as production hardware option and \
select the Enable portable word sizes option on the Code Generation > \
Verification pane for ERT based targets, which will disable the \
preprocessor word size checks.
#endif

#if ( UINT_MAX != (0xFFFFFFFFU) ) || ( INT_MAX != (0x7FFFFFFF) )
#error Code was generated for compiler with different sized uint/int. \
Consider adjusting Test hardware word size settings on the \
Hardware Implementation pane to match your compiler word sizes as \
defined in limits.h of the compiler. Alternatively, you can \
select the Test hardware is the same as production hardware option and \
select the Enable portable word sizes option on the Code Generation > \
Verification pane for ERT based targets, which will disable the \
preprocessor word size checks.
#endif

#if ( ULONG_MAX != (0xFFFFFFFFU) ) || ( LONG_MAX != (0x7FFFFFFF) )
#error Code was generated for compiler with different sized ulong/long. \
Consider adjusting Test hardware word size settings on the \
Hardware Implementation pane to match your compiler word sizes as \
defined in limits.h of the compiler. Alternatively, you can \
select the Test hardware is the same as production hardware option and \
select the Enable portable word sizes option on the Code Generation > \
Verification pane for ERT based targets, which will disable the \
preprocessor word size checks.
#endif

extern real_T rt_roundd_snf(real_T u);
extern real_T rt_remd_snf(real_T u0, real_T u1);
extern real_T rt_modd_snf(real_T u0, real_T u1);
extern real_T rt_hypotd_snf(real_T u0, real_T u1);
extern real_T rt_atan2d_snf(real_T u0, real_T u1);
extern real_T rt_powd_snf(real_T u0, real_T u1);
void BINARYSEARCH_real_T(uint32_T *piLeft, uint32_T *piRght, real_T u, const
  real_T *pData, uint32_T iHi);
void LookUp_real_T_real_T(real_T *pY, const real_T *pYData, real_T u, const
  real_T *pUData, uint32_T iHi);
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
extern void mul_wide_u32(uint32_T in0, uint32_T in1, uint32_T *ptrOutBitsHi,
  uint32_T *ptrOutBitsLo);
extern uint32_T mul_u32_u32_u32_sr29(uint32_T a, uint32_T b);
extern void sfun_discreteVariableDelay(SimStruct *rts);
extern void sfun_fct_op7160ex1_recv(SimStruct *rts);
extern void sfun_SFP2DBL(SimStruct *rts);
extern void RECV_Param(SimStruct *rts);
extern void sfun_fct_op7160ex1_load_in(SimStruct *rts);
extern void sfun_efs_solver_cfg(SimStruct *rts);
extern void sfun_DBL2SFP(SimStruct *rts);
extern void sfun_fct_op7160ex1_send(SimStruct *rts);
extern void opmonitor(SimStruct *rts);
extern void OP_SEND(SimStruct *rts);
extern void opwritefile(SimStruct *rts);
extern void rte_conversion(SimStruct *rts);
extern void rte_logical_operator(SimStruct *rts);
extern void rte_ground(SimStruct *rts);
extern void rte_conversion_ophsdio(SimStruct *rts);
extern void sfun_op7160ex1_event_generator(SimStruct *rts);
extern void optrigger(SimStruct *rts);
extern void sfun_ctrl_op7160ex1(SimStruct *rts);

/* private model entry point functions */
extern void boost_and_two_level__1_sm_ehs_derivatives(void);

#endif                                 /* RTW_HEADER_boost_and_two_level__1_sm_ehs_private_h_ */