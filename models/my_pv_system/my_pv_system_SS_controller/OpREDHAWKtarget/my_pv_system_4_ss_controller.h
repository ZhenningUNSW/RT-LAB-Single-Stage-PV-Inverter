/*
 * my_pv_system_4_ss_controller.h
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
#ifndef RTW_HEADER_my_pv_system_4_ss_controller_h_
#define RTW_HEADER_my_pv_system_4_ss_controller_h_
#include <stddef.h>
#include <string.h>
#include <float.h>
#include <math.h>
#ifndef my_pv_system_4_ss_controller_COMMON_INCLUDES_
# define my_pv_system_4_ss_controller_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "rt_logging.h"
#endif                                 /* my_pv_system_4_ss_controller_COMMON_INCLUDES_ */

#include "my_pv_system_4_ss_controller_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtGetInf.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetBlkStateChangeFlag
# define rtmGetBlkStateChangeFlag(rtm) ((rtm)->ModelData.blkStateChange)
#endif

#ifndef rtmSetBlkStateChangeFlag
# define rtmSetBlkStateChangeFlag(rtm, val) ((rtm)->ModelData.blkStateChange = (val))
#endif

#ifndef rtmGetBlockIO
# define rtmGetBlockIO(rtm)            ((rtm)->ModelData.blockIO)
#endif

#ifndef rtmSetBlockIO
# define rtmSetBlockIO(rtm, val)       ((rtm)->ModelData.blockIO = (val))
#endif

#ifndef rtmGetChecksums
# define rtmGetChecksums(rtm)          ((rtm)->Sizes.checksums)
#endif

#ifndef rtmSetChecksums
# define rtmSetChecksums(rtm, val)     ((rtm)->Sizes.checksums = (val))
#endif

#ifndef rtmGetConstBlockIO
# define rtmGetConstBlockIO(rtm)       ((rtm)->ModelData.constBlockIO)
#endif

#ifndef rtmSetConstBlockIO
# define rtmSetConstBlockIO(rtm, val)  ((rtm)->ModelData.constBlockIO = (val))
#endif

#ifndef rtmGetContStateDisabled
# define rtmGetContStateDisabled(rtm)  ((rtm)->ModelData.contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
# define rtmSetContStateDisabled(rtm, val) ((rtm)->ModelData.contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
# define rtmGetContStates(rtm)         ((rtm)->ModelData.contStates)
#endif

#ifndef rtmSetContStates
# define rtmSetContStates(rtm, val)    ((rtm)->ModelData.contStates = (val))
#endif

#ifndef rtmGetDataMapInfo
# define rtmGetDataMapInfo(rtm)        ()
#endif

#ifndef rtmSetDataMapInfo
# define rtmSetDataMapInfo(rtm, val)   ()
#endif

#ifndef rtmGetDefaultParam
# define rtmGetDefaultParam(rtm)       ((rtm)->ModelData.defaultParam)
#endif

#ifndef rtmSetDefaultParam
# define rtmSetDefaultParam(rtm, val)  ((rtm)->ModelData.defaultParam = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
# define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->ModelData.derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
# define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->ModelData.derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetDirectFeedThrough
# define rtmGetDirectFeedThrough(rtm)  ((rtm)->Sizes.sysDirFeedThru)
#endif

#ifndef rtmSetDirectFeedThrough
# define rtmSetDirectFeedThrough(rtm, val) ((rtm)->Sizes.sysDirFeedThru = (val))
#endif

#ifndef rtmGetErrorStatusFlag
# define rtmGetErrorStatusFlag(rtm)    ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatusFlag
# define rtmSetErrorStatusFlag(rtm, val) ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetFinalTime
# define rtmGetFinalTime(rtm)          ((rtm)->Timing.tFinal)
#endif

#ifndef rtmSetFinalTime
# define rtmSetFinalTime(rtm, val)     ((rtm)->Timing.tFinal = (val))
#endif

#ifndef rtmGetFirstInitCondFlag
# define rtmGetFirstInitCondFlag(rtm)  ()
#endif

#ifndef rtmSetFirstInitCondFlag
# define rtmSetFirstInitCondFlag(rtm, val) ()
#endif

#ifndef rtmGetIntgData
# define rtmGetIntgData(rtm)           ()
#endif

#ifndef rtmSetIntgData
# define rtmSetIntgData(rtm, val)      ()
#endif

#ifndef rtmGetMdlRefGlobalTID
# define rtmGetMdlRefGlobalTID(rtm)    ()
#endif

#ifndef rtmSetMdlRefGlobalTID
# define rtmSetMdlRefGlobalTID(rtm, val) ()
#endif

#ifndef rtmGetMdlRefTriggerTID
# define rtmGetMdlRefTriggerTID(rtm)   ()
#endif

#ifndef rtmSetMdlRefTriggerTID
# define rtmSetMdlRefTriggerTID(rtm, val) ()
#endif

#ifndef rtmGetModelMappingInfo
# define rtmGetModelMappingInfo(rtm)   ((rtm)->SpecialInfo.mappingInfo)
#endif

#ifndef rtmSetModelMappingInfo
# define rtmSetModelMappingInfo(rtm, val) ((rtm)->SpecialInfo.mappingInfo = (val))
#endif

#ifndef rtmGetModelName
# define rtmGetModelName(rtm)          ((rtm)->modelName)
#endif

#ifndef rtmSetModelName
# define rtmSetModelName(rtm, val)     ((rtm)->modelName = (val))
#endif

#ifndef rtmGetNonInlinedSFcns
# define rtmGetNonInlinedSFcns(rtm)    ((rtm)->NonInlinedSFcns)
#endif

#ifndef rtmSetNonInlinedSFcns
# define rtmSetNonInlinedSFcns(rtm, val) ((rtm)->NonInlinedSFcns = (val))
#endif

#ifndef rtmGetNumBlockIO
# define rtmGetNumBlockIO(rtm)         ((rtm)->Sizes.numBlockIO)
#endif

#ifndef rtmSetNumBlockIO
# define rtmSetNumBlockIO(rtm, val)    ((rtm)->Sizes.numBlockIO = (val))
#endif

#ifndef rtmGetNumBlockParams
# define rtmGetNumBlockParams(rtm)     ((rtm)->Sizes.numBlockPrms)
#endif

#ifndef rtmSetNumBlockParams
# define rtmSetNumBlockParams(rtm, val) ((rtm)->Sizes.numBlockPrms = (val))
#endif

#ifndef rtmGetNumBlocks
# define rtmGetNumBlocks(rtm)          ((rtm)->Sizes.numBlocks)
#endif

#ifndef rtmSetNumBlocks
# define rtmSetNumBlocks(rtm, val)     ((rtm)->Sizes.numBlocks = (val))
#endif

#ifndef rtmGetNumContStates
# define rtmGetNumContStates(rtm)      ((rtm)->Sizes.numContStates)
#endif

#ifndef rtmSetNumContStates
# define rtmSetNumContStates(rtm, val) ((rtm)->Sizes.numContStates = (val))
#endif

#ifndef rtmGetNumDWork
# define rtmGetNumDWork(rtm)           ((rtm)->Sizes.numDwork)
#endif

#ifndef rtmSetNumDWork
# define rtmSetNumDWork(rtm, val)      ((rtm)->Sizes.numDwork = (val))
#endif

#ifndef rtmGetNumInputPorts
# define rtmGetNumInputPorts(rtm)      ((rtm)->Sizes.numIports)
#endif

#ifndef rtmSetNumInputPorts
# define rtmSetNumInputPorts(rtm, val) ((rtm)->Sizes.numIports = (val))
#endif

#ifndef rtmGetNumNonSampledZCs
# define rtmGetNumNonSampledZCs(rtm)   ((rtm)->Sizes.numNonSampZCs)
#endif

#ifndef rtmSetNumNonSampledZCs
# define rtmSetNumNonSampledZCs(rtm, val) ((rtm)->Sizes.numNonSampZCs = (val))
#endif

#ifndef rtmGetNumOutputPorts
# define rtmGetNumOutputPorts(rtm)     ((rtm)->Sizes.numOports)
#endif

#ifndef rtmSetNumOutputPorts
# define rtmSetNumOutputPorts(rtm, val) ((rtm)->Sizes.numOports = (val))
#endif

#ifndef rtmGetNumSFcnParams
# define rtmGetNumSFcnParams(rtm)      ((rtm)->Sizes.numSFcnPrms)
#endif

#ifndef rtmSetNumSFcnParams
# define rtmSetNumSFcnParams(rtm, val) ((rtm)->Sizes.numSFcnPrms = (val))
#endif

#ifndef rtmGetNumSFunctions
# define rtmGetNumSFunctions(rtm)      ((rtm)->Sizes.numSFcns)
#endif

#ifndef rtmSetNumSFunctions
# define rtmSetNumSFunctions(rtm, val) ((rtm)->Sizes.numSFcns = (val))
#endif

#ifndef rtmGetNumSampleTimes
# define rtmGetNumSampleTimes(rtm)     ((rtm)->Sizes.numSampTimes)
#endif

#ifndef rtmSetNumSampleTimes
# define rtmSetNumSampleTimes(rtm, val) ((rtm)->Sizes.numSampTimes = (val))
#endif

#ifndef rtmGetNumU
# define rtmGetNumU(rtm)               ((rtm)->Sizes.numU)
#endif

#ifndef rtmSetNumU
# define rtmSetNumU(rtm, val)          ((rtm)->Sizes.numU = (val))
#endif

#ifndef rtmGetNumY
# define rtmGetNumY(rtm)               ((rtm)->Sizes.numY)
#endif

#ifndef rtmSetNumY
# define rtmSetNumY(rtm, val)          ((rtm)->Sizes.numY = (val))
#endif

#ifndef rtmGetOdeF
# define rtmGetOdeF(rtm)               ()
#endif

#ifndef rtmSetOdeF
# define rtmSetOdeF(rtm, val)          ()
#endif

#ifndef rtmGetOdeY
# define rtmGetOdeY(rtm)               ()
#endif

#ifndef rtmSetOdeY
# define rtmSetOdeY(rtm, val)          ()
#endif

#ifndef rtmGetOffsetTimeArray
# define rtmGetOffsetTimeArray(rtm)    ((rtm)->Timing.offsetTimesArray)
#endif

#ifndef rtmSetOffsetTimeArray
# define rtmSetOffsetTimeArray(rtm, val) ((rtm)->Timing.offsetTimesArray = (val))
#endif

#ifndef rtmGetOffsetTimePtr
# define rtmGetOffsetTimePtr(rtm)      ((rtm)->Timing.offsetTimes)
#endif

#ifndef rtmSetOffsetTimePtr
# define rtmSetOffsetTimePtr(rtm, val) ((rtm)->Timing.offsetTimes = (val))
#endif

#ifndef rtmGetOptions
# define rtmGetOptions(rtm)            ((rtm)->Sizes.options)
#endif

#ifndef rtmSetOptions
# define rtmSetOptions(rtm, val)       ((rtm)->Sizes.options = (val))
#endif

#ifndef rtmGetParamIsMalloced
# define rtmGetParamIsMalloced(rtm)    ()
#endif

#ifndef rtmSetParamIsMalloced
# define rtmSetParamIsMalloced(rtm, val) ()
#endif

#ifndef rtmGetPath
# define rtmGetPath(rtm)               ((rtm)->path)
#endif

#ifndef rtmSetPath
# define rtmSetPath(rtm, val)          ((rtm)->path = (val))
#endif

#ifndef rtmGetPerTaskSampleHits
# define rtmGetPerTaskSampleHits(rtm)  ()
#endif

#ifndef rtmSetPerTaskSampleHits
# define rtmSetPerTaskSampleHits(rtm, val) ()
#endif

#ifndef rtmGetPerTaskSampleHitsArray
# define rtmGetPerTaskSampleHitsArray(rtm) ((rtm)->Timing.perTaskSampleHitsArray)
#endif

#ifndef rtmSetPerTaskSampleHitsArray
# define rtmSetPerTaskSampleHitsArray(rtm, val) ((rtm)->Timing.perTaskSampleHitsArray = (val))
#endif

#ifndef rtmGetPerTaskSampleHitsPtr
# define rtmGetPerTaskSampleHitsPtr(rtm) ((rtm)->Timing.perTaskSampleHits)
#endif

#ifndef rtmSetPerTaskSampleHitsPtr
# define rtmSetPerTaskSampleHitsPtr(rtm, val) ((rtm)->Timing.perTaskSampleHits = (val))
#endif

#ifndef rtmGetPrevZCSigState
# define rtmGetPrevZCSigState(rtm)     ((rtm)->ModelData.prevZCSigState)
#endif

#ifndef rtmSetPrevZCSigState
# define rtmSetPrevZCSigState(rtm, val) ((rtm)->ModelData.prevZCSigState = (val))
#endif

#ifndef rtmGetRTWExtModeInfo
# define rtmGetRTWExtModeInfo(rtm)     ((rtm)->extModeInfo)
#endif

#ifndef rtmSetRTWExtModeInfo
# define rtmSetRTWExtModeInfo(rtm, val) ((rtm)->extModeInfo = (val))
#endif

#ifndef rtmGetRTWGeneratedSFcn
# define rtmGetRTWGeneratedSFcn(rtm)   ((rtm)->Sizes.rtwGenSfcn)
#endif

#ifndef rtmSetRTWGeneratedSFcn
# define rtmSetRTWGeneratedSFcn(rtm, val) ((rtm)->Sizes.rtwGenSfcn = (val))
#endif

#ifndef rtmGetRTWLogInfo
# define rtmGetRTWLogInfo(rtm)         ((rtm)->rtwLogInfo)
#endif

#ifndef rtmSetRTWLogInfo
# define rtmSetRTWLogInfo(rtm, val)    ((rtm)->rtwLogInfo = (val))
#endif

#ifndef rtmGetRTWRTModelMethodsInfo
# define rtmGetRTWRTModelMethodsInfo(rtm) ()
#endif

#ifndef rtmSetRTWRTModelMethodsInfo
# define rtmSetRTWRTModelMethodsInfo(rtm, val) ()
#endif

#ifndef rtmGetRTWSfcnInfo
# define rtmGetRTWSfcnInfo(rtm)        ((rtm)->sfcnInfo)
#endif

#ifndef rtmSetRTWSfcnInfo
# define rtmSetRTWSfcnInfo(rtm, val)   ((rtm)->sfcnInfo = (val))
#endif

#ifndef rtmGetRTWSolverInfo
# define rtmGetRTWSolverInfo(rtm)      ((rtm)->solverInfo)
#endif

#ifndef rtmSetRTWSolverInfo
# define rtmSetRTWSolverInfo(rtm, val) ((rtm)->solverInfo = (val))
#endif

#ifndef rtmGetRTWSolverInfoPtr
# define rtmGetRTWSolverInfoPtr(rtm)   ((rtm)->solverInfoPtr)
#endif

#ifndef rtmSetRTWSolverInfoPtr
# define rtmSetRTWSolverInfoPtr(rtm, val) ((rtm)->solverInfoPtr = (val))
#endif

#ifndef rtmGetReservedForXPC
# define rtmGetReservedForXPC(rtm)     ((rtm)->SpecialInfo.xpcData)
#endif

#ifndef rtmSetReservedForXPC
# define rtmSetReservedForXPC(rtm, val) ((rtm)->SpecialInfo.xpcData = (val))
#endif

#ifndef rtmGetRootDWork
# define rtmGetRootDWork(rtm)          ((rtm)->ModelData.dwork)
#endif

#ifndef rtmSetRootDWork
# define rtmSetRootDWork(rtm, val)     ((rtm)->ModelData.dwork = (val))
#endif

#ifndef rtmGetSFunctions
# define rtmGetSFunctions(rtm)         ((rtm)->childSfunctions)
#endif

#ifndef rtmSetSFunctions
# define rtmSetSFunctions(rtm, val)    ((rtm)->childSfunctions = (val))
#endif

#ifndef rtmGetSampleHitArray
# define rtmGetSampleHitArray(rtm)     ((rtm)->Timing.sampleHitArray)
#endif

#ifndef rtmSetSampleHitArray
# define rtmSetSampleHitArray(rtm, val) ((rtm)->Timing.sampleHitArray = (val))
#endif

#ifndef rtmGetSampleHitPtr
# define rtmGetSampleHitPtr(rtm)       ((rtm)->Timing.sampleHits)
#endif

#ifndef rtmSetSampleHitPtr
# define rtmSetSampleHitPtr(rtm, val)  ((rtm)->Timing.sampleHits = (val))
#endif

#ifndef rtmGetSampleTimeArray
# define rtmGetSampleTimeArray(rtm)    ((rtm)->Timing.sampleTimesArray)
#endif

#ifndef rtmSetSampleTimeArray
# define rtmSetSampleTimeArray(rtm, val) ((rtm)->Timing.sampleTimesArray = (val))
#endif

#ifndef rtmGetSampleTimePtr
# define rtmGetSampleTimePtr(rtm)      ((rtm)->Timing.sampleTimes)
#endif

#ifndef rtmSetSampleTimePtr
# define rtmSetSampleTimePtr(rtm, val) ((rtm)->Timing.sampleTimes = (val))
#endif

#ifndef rtmGetSampleTimeTaskIDArray
# define rtmGetSampleTimeTaskIDArray(rtm) ((rtm)->Timing.sampleTimeTaskIDArray)
#endif

#ifndef rtmSetSampleTimeTaskIDArray
# define rtmSetSampleTimeTaskIDArray(rtm, val) ((rtm)->Timing.sampleTimeTaskIDArray = (val))
#endif

#ifndef rtmGetSampleTimeTaskIDPtr
# define rtmGetSampleTimeTaskIDPtr(rtm) ((rtm)->Timing.sampleTimeTaskIDPtr)
#endif

#ifndef rtmSetSampleTimeTaskIDPtr
# define rtmSetSampleTimeTaskIDPtr(rtm, val) ((rtm)->Timing.sampleTimeTaskIDPtr = (val))
#endif

#ifndef rtmGetSimMode
# define rtmGetSimMode(rtm)            ((rtm)->simMode)
#endif

#ifndef rtmSetSimMode
# define rtmSetSimMode(rtm, val)       ((rtm)->simMode = (val))
#endif

#ifndef rtmGetSimTimeStep
# define rtmGetSimTimeStep(rtm)        ((rtm)->Timing.simTimeStep)
#endif

#ifndef rtmSetSimTimeStep
# define rtmSetSimTimeStep(rtm, val)   ((rtm)->Timing.simTimeStep = (val))
#endif

#ifndef rtmGetStartTime
# define rtmGetStartTime(rtm)          ((rtm)->Timing.tStart)
#endif

#ifndef rtmSetStartTime
# define rtmSetStartTime(rtm, val)     ((rtm)->Timing.tStart = (val))
#endif

#ifndef rtmGetStepSize
# define rtmGetStepSize(rtm)           ((rtm)->Timing.stepSize)
#endif

#ifndef rtmSetStepSize
# define rtmSetStepSize(rtm, val)      ((rtm)->Timing.stepSize = (val))
#endif

#ifndef rtmGetStopRequestedFlag
# define rtmGetStopRequestedFlag(rtm)  ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequestedFlag
# define rtmSetStopRequestedFlag(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetTaskCounters
# define rtmGetTaskCounters(rtm)       ()
#endif

#ifndef rtmSetTaskCounters
# define rtmSetTaskCounters(rtm, val)  ()
#endif

#ifndef rtmGetTaskTimeArray
# define rtmGetTaskTimeArray(rtm)      ((rtm)->Timing.tArray)
#endif

#ifndef rtmSetTaskTimeArray
# define rtmSetTaskTimeArray(rtm, val) ((rtm)->Timing.tArray = (val))
#endif

#ifndef rtmGetTimePtr
# define rtmGetTimePtr(rtm)            ((rtm)->Timing.t)
#endif

#ifndef rtmSetTimePtr
# define rtmSetTimePtr(rtm, val)       ((rtm)->Timing.t = (val))
#endif

#ifndef rtmGetTimingData
# define rtmGetTimingData(rtm)         ((rtm)->Timing.timingData)
#endif

#ifndef rtmSetTimingData
# define rtmSetTimingData(rtm, val)    ((rtm)->Timing.timingData = (val))
#endif

#ifndef rtmGetU
# define rtmGetU(rtm)                  ((rtm)->ModelData.inputs)
#endif

#ifndef rtmSetU
# define rtmSetU(rtm, val)             ((rtm)->ModelData.inputs = (val))
#endif

#ifndef rtmGetVarNextHitTimesListPtr
# define rtmGetVarNextHitTimesListPtr(rtm) ((rtm)->Timing.varNextHitTimesList)
#endif

#ifndef rtmSetVarNextHitTimesListPtr
# define rtmSetVarNextHitTimesListPtr(rtm, val) ((rtm)->Timing.varNextHitTimesList = (val))
#endif

#ifndef rtmGetY
# define rtmGetY(rtm)                  ((rtm)->ModelData.outputs)
#endif

#ifndef rtmSetY
# define rtmSetY(rtm, val)             ((rtm)->ModelData.outputs = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
# define rtmGetZCCacheNeedsReset(rtm)  ((rtm)->ModelData.zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
# define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->ModelData.zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetZCSignalValues
# define rtmGetZCSignalValues(rtm)     ((rtm)->ModelData.zcSignalValues)
#endif

#ifndef rtmSetZCSignalValues
# define rtmSetZCSignalValues(rtm, val) ((rtm)->ModelData.zcSignalValues = (val))
#endif

#ifndef rtmGet_TimeOfLastOutput
# define rtmGet_TimeOfLastOutput(rtm)  ((rtm)->Timing.timeOfLastOutput)
#endif

#ifndef rtmSet_TimeOfLastOutput
# define rtmSet_TimeOfLastOutput(rtm, val) ((rtm)->Timing.timeOfLastOutput = (val))
#endif

#ifndef rtmGetdX
# define rtmGetdX(rtm)                 ((rtm)->ModelData.derivs)
#endif

#ifndef rtmSetdX
# define rtmSetdX(rtm, val)            ((rtm)->ModelData.derivs = (val))
#endif

#ifndef rtmGetChecksumVal
# define rtmGetChecksumVal(rtm, idx)   ((rtm)->Sizes.checksums[idx])
#endif

#ifndef rtmSetChecksumVal
# define rtmSetChecksumVal(rtm, idx, val) ((rtm)->Sizes.checksums[idx] = (val))
#endif

#ifndef rtmGetDWork
# define rtmGetDWork(rtm, idx)         ((rtm)->ModelData.dwork[idx])
#endif

#ifndef rtmSetDWork
# define rtmSetDWork(rtm, idx, val)    ((rtm)->ModelData.dwork[idx] = (val))
#endif

#ifndef rtmGetOffsetTime
# define rtmGetOffsetTime(rtm, idx)    ((rtm)->Timing.offsetTimes[idx])
#endif

#ifndef rtmSetOffsetTime
# define rtmSetOffsetTime(rtm, idx, val) ((rtm)->Timing.offsetTimes[idx] = (val))
#endif

#ifndef rtmGetSFunction
# define rtmGetSFunction(rtm, idx)     ((rtm)->childSfunctions[idx])
#endif

#ifndef rtmSetSFunction
# define rtmSetSFunction(rtm, idx, val) ((rtm)->childSfunctions[idx] = (val))
#endif

#ifndef rtmGetSampleTime
# define rtmGetSampleTime(rtm, idx)    ((rtm)->Timing.sampleTimes[idx])
#endif

#ifndef rtmSetSampleTime
# define rtmSetSampleTime(rtm, idx, val) ((rtm)->Timing.sampleTimes[idx] = (val))
#endif

#ifndef rtmGetSampleTimeTaskID
# define rtmGetSampleTimeTaskID(rtm, idx) ((rtm)->Timing.sampleTimeTaskIDPtr[idx])
#endif

#ifndef rtmSetSampleTimeTaskID
# define rtmSetSampleTimeTaskID(rtm, idx, val) ((rtm)->Timing.sampleTimeTaskIDPtr[idx] = (val))
#endif

#ifndef rtmGetVarNextHitTimeList
# define rtmGetVarNextHitTimeList(rtm, idx) ((rtm)->Timing.varNextHitTimesList[idx])
#endif

#ifndef rtmSetVarNextHitTimeList
# define rtmSetVarNextHitTimeList(rtm, idx, val) ((rtm)->Timing.varNextHitTimesList[idx] = (val))
#endif

#ifndef rtmIsContinuousTask
# define rtmIsContinuousTask(rtm, tid) 0
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmIsSampleHit
# define rtmIsSampleHit(rtm, sti, tid) ((rtm)->Timing.sampleHits[(rtm)->Timing.sampleTimeTaskIDPtr[sti]])
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmSetT
# define rtmSetT(rtm, val)                                       /* Do Nothing */
#endif

#ifndef rtmGetTFinal
# define rtmGetTFinal(rtm)             ((rtm)->Timing.tFinal)
#endif

#ifndef rtmSetTFinal
# define rtmSetTFinal(rtm, val)        ((rtm)->Timing.tFinal = (val))
#endif

#ifndef rtmGetTPtr
# define rtmGetTPtr(rtm)               ((rtm)->Timing.t)
#endif

#ifndef rtmSetTPtr
# define rtmSetTPtr(rtm, val)          ((rtm)->Timing.t = (val))
#endif

#ifndef rtmGetTStart
# define rtmGetTStart(rtm)             ((rtm)->Timing.tStart)
#endif

#ifndef rtmSetTStart
# define rtmSetTStart(rtm, val)        ((rtm)->Timing.tStart = (val))
#endif

#ifndef rtmGetTaskTime
# define rtmGetTaskTime(rtm, sti)      (rtmGetTPtr((rtm))[(rtm)->Timing.sampleTimeTaskIDPtr[sti]])
#endif

#ifndef rtmSetTaskTime
# define rtmSetTaskTime(rtm, sti, val) (rtmGetTPtr((rtm))[sti] = (val))
#endif

#ifndef rtmGetTimeOfLastOutput
# define rtmGetTimeOfLastOutput(rtm)   ((rtm)->Timing.timeOfLastOutput)
#endif

#ifdef rtmGetRTWSolverInfo
#undef rtmGetRTWSolverInfo
#endif

#define rtmGetRTWSolverInfo(rtm)       &((rtm)->solverInfo)

/* Definition for use in the target main file */
#define my_pv_system_4_ss_controller_rtModel RT_MODEL_my_pv_system_4_ss_controller_T

/* user code (top of header file) */
/* System '<Root>' */
/* Opal-RT Technologies */
extern int opalSizeDwork;
extern int opalSizeBlockIO;
extern int opalSizeRTP;

#ifdef USE_RTMODEL

extern void * pRtModel;                //pointer on the RTmodel struc
int_T Opal_rtmGetNumBlockParams(void *ptr);
int_T Opal_rtmGetNumBlockIO(void *ptr);

#endif

/* Block signals (auto storage) */
typedef struct {
  creal_T RealImagtoComplex;           /* '<S15>/Real-Imag to Complex' */
  creal_T MagnitudeAngletoComplex;     /* '<S10>/Magnitude-Angle to Complex' */
  creal_T RealImagtoComplex_n;         /* '<S12>/Real-Imag to Complex' */
  creal_T RealImagtoComplex_l;         /* '<S35>/Real-Imag to Complex' */
  real_T SFunction;                    /* '<S1>/S-Function' */
  real_T Sum;                          /* '<S1>/Sum' */
  real_T Memory[4];                    /* '<S2>/Memory' */
  real_T SFunction_p[7];               /* '<S88>/S-Function' */
  real_T RateTransition1;              /* '<S2>/Rate Transition1' */
  real_T Apu;                          /* '<S10>/A->pu' */
  real_T UnitDelay;                    /* '<S30>/Unit Delay' */
  real_T avoiddivisionbyzero;          /* '<S19>/avoid division by zero' */
  real_T MathFunction;                 /* '<S19>/Math Function' */
  real_T Gain;                         /* '<S19>/Gain' */
  real_T SFunction_o;                  /* '<S56>/S-Function' */
  real_T DiscreteTimeIntegrator;       /* '<S30>/Discrete-Time Integrator' */
  real_T MathFunction_i;               /* '<S30>/Math Function' */
  real_T FirstcycleofsimulationId092Iq0;/* '<S19>/First cycle of simulation Id=0.92, Iq=0' */
  real_T Switch[2];                    /* '<S19>/Switch' */
  real_T Duk[2];                       /* '<S51>/D*u(k)' */
  real_T x1k[2];                       /* '<S51>/Delay_x1' */
  real_T C11[2];                       /* '<S54>/C11' */
  real_T x2k[2];                       /* '<S51>/Delay_x2' */
  real_T C12[2];                       /* '<S54>/C12' */
  real_T sum2[2];                      /* '<S54>/sum2' */
  real_T yk[2];                        /* '<S51>/C*X(k)+D*u(k)' */
  real_T UnitDelay2;                   /* '<S4>/Unit Delay2' */
  real_T Sum_e[2];                     /* '<S7>/Sum' */
  real_T ProportionalGain[2];          /* '<S14>/Proportional Gain' */
  real_T Integrator[2];                /* '<S14>/Integrator' */
  real_T Sum_g[2];                     /* '<S14>/Sum' */
  real_T Saturate[2];                  /* '<S14>/Saturate' */
  real_T RateTransition;               /* '<S2>/Rate Transition' */
  real_T Vpu;                          /* '<S10>/V->pu' */
  real_T TrigonometricFunction;        /* '<S15>/Trigonometric Function' */
  real_T Gain1;                        /* '<S15>/Gain1' */
  real_T Product1;                     /* '<S15>/Product1' */
  real_T Integ4;                       /* '<S22>/Integ4' */
  real_T Freq;                         /* '<S22>/To avoid division  by zero' */
  real_T Numberofsamplespercycle;      /* '<S22>/Number of samples per cycle' */
  real_T RoundingFunction;             /* '<S22>/Rounding Function' */
  real_T Delay;                        /* '<S22>/Gain' */
  real_T SFunction_pj;                 /* '<S24>/S-Function' */
  real_T UnitDelay_f;                  /* '<S23>/Unit Delay' */
  real_T DigitalClock;                 /* '<S22>/Digital  Clock' */
  real_T UnitDelay1;                   /* '<S22>/Unit Delay1' */
  real_T Switch_g;                     /* '<S22>/Switch' */
  real_T TrigonometricFunction3;       /* '<S15>/Trigonometric Function3' */
  real_T Gain3;                        /* '<S15>/Gain3' */
  real_T Product2;                     /* '<S15>/Product2' */
  real_T Integ4_h;                     /* '<S25>/Integ4' */
  real_T Freq_m;                       /* '<S25>/To avoid division  by zero' */
  real_T Numberofsamplespercycle_m;    /* '<S25>/Number of samples per cycle' */
  real_T RoundingFunction_p;           /* '<S25>/Rounding Function' */
  real_T Delay_d;                      /* '<S25>/Gain' */
  real_T SFunction_i;                  /* '<S27>/S-Function' */
  real_T UnitDelay_l;                  /* '<S26>/Unit Delay' */
  real_T DigitalClock_m;               /* '<S25>/Digital  Clock' */
  real_T UnitDelay1_e;                 /* '<S25>/Unit Delay1' */
  real_T Switch_k;                     /* '<S25>/Switch' */
  real_T ComplextoMagnitudeAngle_o1;   /* '<S15>/Complex to Magnitude-Angle' */
  real_T ComplextoMagnitudeAngle_o2;   /* '<S15>/Complex to Magnitude-Angle' */
  real_T RadDeg;                       /* '<S15>/Rad->Deg.' */
  real_T torad;                        /* '<S10>/to-rad' */
  real_T ComplextoRealImag_o1;         /* '<S10>/Complex to Real-Imag' */
  real_T ComplextoRealImag_o2;         /* '<S10>/Complex to Real-Imag' */
  real_T Rff;                          /* '<S7>/Rff ' */
  real_T Lff;                          /* '<S7>/Lff  ' */
  real_T Feedforward;                  /* '<S7>/Add1' */
  real_T Rff_e;                        /* '<S7>/Rff' */
  real_T Lff_i;                        /* '<S7>/Lff' */
  real_T Add3;                         /* '<S7>/Add3' */
  real_T Add2[2];                      /* '<S7>/Add2' */
  real_T IntegralGain[2];              /* '<S14>/Integral Gain' */
  real_T Saturation[2];                /* '<S7>/Saturation' */
  real_T RateTransition3;              /* '<S2>/Rate Transition3' */
  real_T RateTransition4;              /* '<S2>/Rate Transition4' */
  real_T DigitalClock_a;               /* '<S28>/Digital  Clock' */
  real_T RateTransition2;              /* '<S2>/Rate Transition2' */
  real_T Integ4_c;                     /* '<S28>/Integ4' */
  real_T K1;                           /* '<S28>/K1' */
  real_T SFunction_pd;                 /* '<S29>/S-Function' */
  real_T UnitDelay_h;                  /* '<S28>/Unit Delay' */
  real_T UnitDelay1_p;                 /* '<S28>/Unit Delay1' */
  real_T Switch_c;                     /* '<S28>/Switch' */
  real_T TrigonometricFunction2;       /* '<S30>/Trigonometric Function2' */
  real_T Product1_c;                   /* '<S30>/Product1' */
  real_T Integ4_e;                     /* '<S44>/Integ4' */
  real_T Freq_c;                       /* '<S44>/To avoid division  by zero' */
  real_T Numberofsamplespercycle_f;    /* '<S44>/Number of samples per cycle' */
  real_T RoundingFunction_b;           /* '<S44>/Rounding Function' */
  real_T Delay_m;                      /* '<S44>/Gain' */
  real_T SFunction_e;                  /* '<S46>/S-Function' */
  real_T UnitDelay_o;                  /* '<S45>/Unit Delay' */
  real_T DigitalClock_e;               /* '<S44>/Digital  Clock' */
  real_T UnitDelay1_a;                 /* '<S44>/Unit Delay1' */
  real_T Switch_d;                     /* '<S44>/Switch' */
  real_T Divide;                       /* '<S30>/Divide' */
  real_T DiscreteDerivative;           /* '<S32>/Discrete Derivative ' */
  real_T DiscreteTimeIntegrator_g;     /* '<S32>/Discrete-Time Integrator' */
  real_T Kp4;                          /* '<S32>/Kp4' */
  real_T Sum6;                         /* '<S32>/Sum6' */
  real_T Saturation1;                  /* '<S32>/Saturation1' */
  real_T Gain10;                       /* '<S30>/Gain10' */
  real_T RateLimiter;                  /* '<S30>/Rate Limiter' */
  real_T x1k_d;                        /* '<S47>/Delay_x1' */
  real_T A11;                          /* '<S48>/A11' */
  real_T x2k_k;                        /* '<S47>/Delay_x2' */
  real_T A12;                          /* '<S48>/A12' */
  real_T A21;                          /* '<S48>/A21' */
  real_T A22;                          /* '<S48>/A22' */
  real_T sum2_c;                       /* '<S48>/sum2' */
  real_T sum3;                         /* '<S48>/sum3' */
  real_T B11;                          /* '<S49>/B11' */
  real_T x1k1;                         /* '<S47>/A*x1(k) + B*u1(k) ' */
  real_T B21;                          /* '<S49>/B21' */
  real_T x2k1;                         /* '<S47>/A*x2(k) + B*u2(k)' */
  real_T Duk_b;                        /* '<S47>/D*u(k)' */
  real_T C11_d;                        /* '<S50>/C11' */
  real_T C12_g;                        /* '<S50>/C12' */
  real_T sum2_o;                       /* '<S50>/sum2' */
  real_T yk_e;                         /* '<S47>/C*X(k)+D*u(k)' */
  real_T A11_e[2];                     /* '<S52>/A11' */
  real_T A12_i[2];                     /* '<S52>/A12' */
  real_T A21_o[2];                     /* '<S52>/A21' */
  real_T A22_p[2];                     /* '<S52>/A22' */
  real_T sum2_e[2];                    /* '<S52>/sum2' */
  real_T sum3_g[2];                    /* '<S52>/sum3' */
  real_T B11_o[2];                     /* '<S53>/B11' */
  real_T x1k1_l[2];                    /* '<S51>/A*x1(k) + B*u1(k) ' */
  real_T B21_p[2];                     /* '<S53>/B21' */
  real_T x2k1_k[2];                    /* '<S51>/A*x2(k) + B*u2(k)' */
  real_T Constant1;                    /* '<S19>/Constant1' */
  real_T Add3_n;                       /* '<S61>/Add3' */
  real_T DigitalClock_mz;              /* '<S83>/Digital Clock' */
  real_T Add1;                         /* '<S83>/Add1' */
  real_T MathFunction_b;               /* '<S83>/Math Function' */
  real_T ib1;                          /* '<S83>/1\ib1' */
  real_T LookupTable;                  /* '<S83>/Lookup Table' */
  real_T Add3_o;                       /* '<S83>/Add3' */
  real_T Gain1_p;                      /* '<S61>/Gain1' */
  real_T MUL1;                         /* '<S61>/MUL1' */
  real_T Add4;                         /* '<S61>/Add4' */
  real_T UnitDelay_m;                  /* '<S4>/Unit Delay' */
  real_T Gain_a;                       /* '<S66>/Gain' */
  real_T DataTypeConversion[4];        /* '<S11>/Data Type Conversion' */
  real_T Switch_d2;                    /* '<S4>/Switch' */
  real_T Add1_a;                       /* '<S12>/Add1' */
  real_T UnitDelay3[2];                /* '<S4>/Unit Delay3' */
  real_T Gain1_g;                      /* '<S12>/Gain1' */
  real_T Product;                      /* '<S12>/Product' */
  real_T Product1_d[2];                /* '<S12>/Product1' */
  real_T ComplextoMagnitudeAngle_o1_n; /* '<S12>/Complex to Magnitude-Angle' */
  real_T ComplextoMagnitudeAngle_o2_p; /* '<S12>/Complex to Magnitude-Angle' */
  real_T Add2_f;                       /* '<S12>/Add2' */
  real_T TrigonometricFunction_h;      /* '<S12>/Trigonometric Function' */
  real_T Product2_n;                   /* '<S12>/Product2' */
  real_T UnitDelay1_c;                 /* '<S4>/Unit Delay1' */
  real_T Sum_k;                        /* '<S13>/Sum' */
  real_T Rtot_pu2;                     /* '<S13>/Rtot_pu2' */
  real_T IntegralGain_a;               /* '<S85>/Integral Gain' */
  real_T Integrator_e;                 /* '<S85>/Integrator' */
  real_T ProportionalGain_b;           /* '<S85>/Proportional Gain' */
  real_T Sum_i;                        /* '<S85>/Sum' */
  real_T Saturate_p;                   /* '<S85>/Saturate' */
  real_T Fcn;                          /* '<S60>/Fcn' */
  real_T Fcn1;                         /* '<S60>/Fcn1' */
  real_T Fcn_l;                        /* '<S59>/Fcn' */
  real_T Fcn1_f;                       /* '<S59>/Fcn1' */
  real_T Switch_p[2];                  /* '<S55>/Switch' */
  real_T Sum1;                         /* '<S45>/Sum1' */
  real_T Sum5;                         /* '<S45>/Sum5' */
  real_T Product5;                     /* '<S45>/Product5' */
  real_T Gain1_pw;                     /* '<S45>/Gain1' */
  real_T Sum4;                         /* '<S45>/Sum4' */
  real_T Product2_b;                   /* '<S45>/Product2' */
  real_T Product4;                     /* '<S45>/Product4' */
  real_T Sum7;                         /* '<S44>/Sum7' */
  real_T Meanvalue;                    /* '<S44>/Product' */
  real_T Sum5_b;                       /* '<S44>/Sum5' */
  real_T TrigonometricFunction_p;      /* '<S35>/Trigonometric Function' */
  real_T Gain1_i;                      /* '<S35>/Gain1' */
  real_T Product1_o;                   /* '<S35>/Product1' */
  real_T Integ4_d;                     /* '<S38>/Integ4' */
  real_T Freq_h;                       /* '<S38>/To avoid division  by zero' */
  real_T Numberofsamplespercycle_h;    /* '<S38>/Number of samples per cycle' */
  real_T RoundingFunction_e;           /* '<S38>/Rounding Function' */
  real_T Delay_db;                     /* '<S38>/Gain' */
  real_T SFunction_ij;                 /* '<S40>/S-Function' */
  real_T UnitDelay_b;                  /* '<S39>/Unit Delay' */
  real_T DigitalClock_o;               /* '<S38>/Digital  Clock' */
  real_T UnitDelay1_c2;                /* '<S38>/Unit Delay1' */
  real_T Switch_j;                     /* '<S38>/Switch' */
  real_T TrigonometricFunction3_p;     /* '<S35>/Trigonometric Function3' */
  real_T Gain3_c;                      /* '<S35>/Gain3' */
  real_T Product2_i;                   /* '<S35>/Product2' */
  real_T Integ4_dy;                    /* '<S41>/Integ4' */
  real_T Freq_mz;                      /* '<S41>/To avoid division  by zero' */
  real_T Numberofsamplespercycle_d;    /* '<S41>/Number of samples per cycle' */
  real_T RoundingFunction_h;           /* '<S41>/Rounding Function' */
  real_T Delay_i;                      /* '<S41>/Gain' */
  real_T SFunction_ed;                 /* '<S43>/S-Function' */
  real_T UnitDelay_fb;                 /* '<S42>/Unit Delay' */
  real_T DigitalClock_l;               /* '<S41>/Digital  Clock' */
  real_T UnitDelay1_j;                 /* '<S41>/Unit Delay1' */
  real_T Switch_m;                     /* '<S41>/Switch' */
  real_T ComplextoMagnitudeAngle_o1_d; /* '<S35>/Complex to Magnitude-Angle' */
  real_T ComplextoMagnitudeAngle_o2_i; /* '<S35>/Complex to Magnitude-Angle' */
  real_T RadDeg_h;                     /* '<S35>/Rad->Deg.' */
  real_T Saturation_n;                 /* '<S31>/Saturation' */
  real_T MathFunction_k;               /* '<S31>/Math Function' */
  real_T Sum1_c;                       /* '<S42>/Sum1' */
  real_T Sum5_o;                       /* '<S42>/Sum5' */
  real_T Product5_e;                   /* '<S42>/Product5' */
  real_T Gain1_n;                      /* '<S42>/Gain1' */
  real_T Sum4_g;                       /* '<S42>/Sum4' */
  real_T Product2_m;                   /* '<S42>/Product2' */
  real_T Product4_d;                   /* '<S42>/Product4' */
  real_T Sum7_m;                       /* '<S41>/Sum7' */
  real_T Meanvalue_g;                  /* '<S41>/Product' */
  real_T Sum5_i;                       /* '<S41>/Sum5' */
  real_T Sum1_p;                       /* '<S39>/Sum1' */
  real_T Sum5_k;                       /* '<S39>/Sum5' */
  real_T Product5_n;                   /* '<S39>/Product5' */
  real_T Gain1_o;                      /* '<S39>/Gain1' */
  real_T Sum4_n;                       /* '<S39>/Sum4' */
  real_T Product2_d;                   /* '<S39>/Product2' */
  real_T Product4_o;                   /* '<S39>/Product4' */
  real_T Sum7_e;                       /* '<S38>/Sum7' */
  real_T Meanvalue_h;                  /* '<S38>/Product' */
  real_T Sum5_f;                       /* '<S38>/Sum5' */
  real_T Gain1_c;                      /* '<S28>/Gain1' */
  real_T Gain_f;                       /* '<S28>/Gain' */
  real_T Correction;                   /* '<S28>/Sum1' */
  real_T Sum7_i;                       /* '<S28>/Sum7' */
  real_T Mean;                         /* '<S28>/Product' */
  real_T Sum5_c;                       /* '<S28>/Sum5' */
  real_T Sum1_cr;                      /* '<S26>/Sum1' */
  real_T Sum5_g;                       /* '<S26>/Sum5' */
  real_T Product5_a;                   /* '<S26>/Product5' */
  real_T Gain1_if;                     /* '<S26>/Gain1' */
  real_T Sum4_b;                       /* '<S26>/Sum4' */
  real_T Product2_dp;                  /* '<S26>/Product2' */
  real_T Product4_g;                   /* '<S26>/Product4' */
  real_T Sum7_i1;                      /* '<S25>/Sum7' */
  real_T Meanvalue_e;                  /* '<S25>/Product' */
  real_T Sum5_g2;                      /* '<S25>/Sum5' */
  real_T Sum1_c3;                      /* '<S23>/Sum1' */
  real_T Sum5_cr;                      /* '<S23>/Sum5' */
  real_T Product5_j;                   /* '<S23>/Product5' */
  real_T Gain1_b;                      /* '<S23>/Gain1' */
  real_T Sum4_e;                       /* '<S23>/Sum4' */
  real_T Product2_f;                   /* '<S23>/Product2' */
  real_T Product4_c;                   /* '<S23>/Product4' */
  real_T Sum7_eg;                      /* '<S22>/Sum7' */
  real_T Meanvalue_c;                  /* '<S22>/Product' */
  real_T Sum5_m;                       /* '<S22>/Sum5' */
  real_T TmpSignalConversionAtSFunctionI[4];/* '<S4>/MPPT Controller using Perturbe  & Observe technique  ' */
  real_T D;                            /* '<S4>/MPPT Controller using Perturbe  & Observe technique  ' */
  uint8_T Compare;                     /* '<S57>/Compare' */
  uint8_T Compare_e;                   /* '<S58>/Compare' */
  boolean_T RelationalOperator;        /* '<S22>/Relational Operator' */
  boolean_T RelationalOperator_m;      /* '<S25>/Relational Operator' */
  boolean_T RelationalOperator_j;      /* '<S28>/Relational Operator' */
  boolean_T RelationalOperator_l;      /* '<S44>/Relational Operator' */
  boolean_T RelationalOperator1;       /* '<S66>/Relational Operator1' */
  boolean_T RelationalOperator3;       /* '<S66>/Relational Operator3' */
  boolean_T LogicalOperator4[2];       /* '<S11>/Logical Operator4' */
  boolean_T RelationalOperator_d;      /* '<S38>/Relational Operator' */
  boolean_T RelationalOperator_k;      /* '<S41>/Relational Operator' */
} B_my_pv_system_4_ss_controller_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  real_T UnitDelay_DSTATE;             /* '<S30>/Unit Delay' */
  real_T DiscreteTimeIntegrator_DSTATE;/* '<S30>/Discrete-Time Integrator' */
  real_T Delay_x1_DSTATE[2];           /* '<S51>/Delay_x1' */
  real_T Delay_x2_DSTATE[2];           /* '<S51>/Delay_x2' */
  real_T UnitDelay2_DSTATE;            /* '<S4>/Unit Delay2' */
  real_T Integrator_DSTATE[2];         /* '<S14>/Integrator' */
  real_T Integ4_DSTATE;                /* '<S22>/Integ4' */
  real_T UnitDelay_DSTATE_g;           /* '<S23>/Unit Delay' */
  real_T UnitDelay1_DSTATE;            /* '<S22>/Unit Delay1' */
  real_T Integ4_DSTATE_d;              /* '<S25>/Integ4' */
  real_T UnitDelay_DSTATE_i;           /* '<S26>/Unit Delay' */
  real_T UnitDelay1_DSTATE_b;          /* '<S25>/Unit Delay1' */
  real_T Integ4_DSTATE_d5;             /* '<S28>/Integ4' */
  real_T UnitDelay_DSTATE_iq;          /* '<S28>/Unit Delay' */
  real_T UnitDelay1_DSTATE_a;          /* '<S28>/Unit Delay1' */
  real_T Integ4_DSTATE_k;              /* '<S44>/Integ4' */
  real_T UnitDelay_DSTATE_a;           /* '<S45>/Unit Delay' */
  real_T UnitDelay1_DSTATE_i;          /* '<S44>/Unit Delay1' */
  real_T DiscreteDerivative_states;    /* '<S32>/Discrete Derivative ' */
  real_T DiscreteTimeIntegrator_DSTATE_d;/* '<S32>/Discrete-Time Integrator' */
  real_T Delay_x1_DSTATE_d;            /* '<S47>/Delay_x1' */
  real_T Delay_x2_DSTATE_m;            /* '<S47>/Delay_x2' */
  real_T UnitDelay_DSTATE_il;          /* '<S4>/Unit Delay' */
  real_T UnitDelay3_DSTATE[2];         /* '<S4>/Unit Delay3' */
  real_T UnitDelay1_DSTATE_e;          /* '<S4>/Unit Delay1' */
  real_T Integrator_DSTATE_d;          /* '<S85>/Integrator' */
  real_T Integ4_DSTATE_ky;             /* '<S38>/Integ4' */
  real_T UnitDelay_DSTATE_l;           /* '<S39>/Unit Delay' */
  real_T UnitDelay1_DSTATE_o;          /* '<S38>/Unit Delay1' */
  real_T Integ4_DSTATE_j;              /* '<S41>/Integ4' */
  real_T UnitDelay_DSTATE_d;           /* '<S42>/Unit Delay' */
  real_T UnitDelay1_DSTATE_k;          /* '<S41>/Unit Delay1' */
  real_T SFunction_PreviousInput;      /* '<S1>/S-Function' */
  real_T Memory_PreviousInput[4];      /* '<S2>/Memory' */
  real_T DiscreteDerivative_tmp;       /* '<S32>/Discrete Derivative ' */
  real_T PrevY;                        /* '<S30>/Rate Limiter' */
  real_T Vold;                         /* '<S4>/MPPT Controller using Perturbe  & Observe technique  ' */
  real_T Pold;                         /* '<S4>/MPPT Controller using Perturbe  & Observe technique  ' */
  real_T Dold;                         /* '<S4>/MPPT Controller using Perturbe  & Observe technique  ' */
  real_T SFunction_RWORK;              /* '<S56>/S-Function' */
  real_T SFunction_RWORK_a;            /* '<S24>/S-Function' */
  real_T SFunction_RWORK_o;            /* '<S27>/S-Function' */
  real_T SFunction_RWORK_k;            /* '<S29>/S-Function' */
  real_T SFunction_RWORK_j;            /* '<S46>/S-Function' */
  real_T SFunction_RWORK_kx;           /* '<S40>/S-Function' */
  real_T SFunction_RWORK_m;            /* '<S43>/S-Function' */
  void *SFunction_PWORK;               /* '<S88>/S-Function' */
  void *SFunction_PWORK_a;             /* '<S56>/S-Function' */
  struct {
    void *LoggedData;
  } PI_Ireg1_PWORK;                    /* '<S7>/PI_Ireg1' */

  void *SFunction_PWORK_c;             /* '<S24>/S-Function' */
  void *SFunction_PWORK_f;             /* '<S27>/S-Function' */
  void *SFunction_PWORK_i;             /* '<S29>/S-Function' */
  void *SFunction_PWORK_a5;            /* '<S46>/S-Function' */
  void *SFunction_PWORK_n;             /* '<S40>/S-Function' */
  void *SFunction_PWORK_p;             /* '<S43>/S-Function' */
  int_T SFunction_IWORK;               /* '<S56>/S-Function' */
  int_T SFunction_IWORK_f;             /* '<S24>/S-Function' */
  int_T SFunction_IWORK_p;             /* '<S27>/S-Function' */
  int_T SFunction_IWORK_b;             /* '<S29>/S-Function' */
  int_T SFunction_IWORK_c;             /* '<S46>/S-Function' */
  int_T SFunction_IWORK_bf[5];         /* '<S87>/S-Function' */
  int_T SFunction_IWORK_pd;            /* '<S40>/S-Function' */
  int_T SFunction_IWORK_c4;            /* '<S43>/S-Function' */
  uint8_T Integ4_SYSTEM_ENABLE;        /* '<S22>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_e;      /* '<S25>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_d;      /* '<S28>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_k;      /* '<S44>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_g;      /* '<S38>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_gx;     /* '<S41>/Integ4' */
  boolean_T Vold_not_empty;            /* '<S4>/MPPT Controller using Perturbe  & Observe technique  ' */
  boolean_T AutomaticGainControl_MODE; /* '<S30>/Automatic Gain Control' */
} DW_my_pv_system_4_ss_controller_T;

/* Backward compatible GRT Identifiers */
#define rtB                            my_pv_system_4_ss_controller_B
#define BlockIO                        B_my_pv_system_4_ss_controller_T
#define rtP                            my_pv_system_4_ss_controller_P
#define Parameters                     P_my_pv_system_4_ss_controller_T
#define rtDWork                        my_pv_system_4_ss_controller_DW
#define D_Work                         DW_my_pv_system_4_ss_controller_T

/* Parameters (auto storage) */
struct P_my_pv_system_4_ss_controller_T_ {
  real_T PLL_AGC;                      /* Mask Parameter: PLL_AGC
                                        * Referenced by: '<S30>/Constant1'
                                        */
  real_T AlphaBetaZerotodq0_Alignment; /* Mask Parameter: AlphaBetaZerotodq0_Alignment
                                        * Referenced by: '<S55>/Constant'
                                        */
  real_T InverterControl_Fnom;         /* Mask Parameter: InverterControl_Fnom
                                        * Referenced by:
                                        *   '<S12>/Constant4'
                                        *   '<S19>/First cycle of simulation Id=0.92, Iq=0'
                                        */
  real_T InverterControl_Increment_MPPT;/* Mask Parameter: InverterControl_Increment_MPPT
                                         * Referenced by: '<S8>/Iph_3'
                                         */
  real_T Discrete_Init;                /* Mask Parameter: Discrete_Init
                                        * Referenced by: '<S32>/Discrete-Time Integrator'
                                        */
  real_T Discrete_Kd;                  /* Mask Parameter: Discrete_Kd
                                        * Referenced by: '<S32>/Discrete Derivative '
                                        */
  real_T InverterControl_Ki_Ireg;      /* Mask Parameter: InverterControl_Ki_Ireg
                                        * Referenced by: '<S14>/Integral Gain'
                                        */
  real_T InverterControl_Ki_VDCreg;    /* Mask Parameter: InverterControl_Ki_VDCreg
                                        * Referenced by: '<S85>/Integral Gain'
                                        */
  real_T Discrete_Kp;                  /* Mask Parameter: Discrete_Kp
                                        * Referenced by: '<S32>/Kp4'
                                        */
  real_T InverterControl_Kp_Ireg;      /* Mask Parameter: InverterControl_Kp_Ireg
                                        * Referenced by: '<S14>/Proportional Gain'
                                        */
  real_T InverterControl_Kp_VDCreg;    /* Mask Parameter: InverterControl_Kp_VDCreg
                                        * Referenced by: '<S85>/Proportional Gain'
                                        */
  real_T PI_LowerSaturationLimit;      /* Mask Parameter: PI_LowerSaturationLimit
                                        * Referenced by: '<S14>/Saturate'
                                        */
  real_T PI_LowerSaturationLimit_p;    /* Mask Parameter: PI_LowerSaturationLimit_p
                                        * Referenced by: '<S85>/Saturate'
                                        */
  real_T PWM_Generator_MinMax[2];      /* Mask Parameter: PWM_Generator_MinMax
                                        * Referenced by: '<S11>/Constant10'
                                        */
  real_T InverterControl_Pnom;         /* Mask Parameter: InverterControl_Pnom
                                        * Referenced by: '<S10>/A->pu'
                                        */
  real_T InverterControl_Ts_Control;   /* Mask Parameter: InverterControl_Ts_Control
                                        * Referenced by:
                                        *   '<S12>/Constant4'
                                        *   '<S22>/Gain'
                                        *   '<S25>/Gain'
                                        *   '<S44>/Gain'
                                        *   '<S38>/Gain'
                                        *   '<S41>/Gain'
                                        */
  real_T PI_UpperSaturationLimit;      /* Mask Parameter: PI_UpperSaturationLimit
                                        * Referenced by: '<S14>/Saturate'
                                        */
  real_T PI_UpperSaturationLimit_f;    /* Mask Parameter: PI_UpperSaturationLimit_f
                                        * Referenced by: '<S85>/Saturate'
                                        */
  real_T InverterControl_Vdc_ref_Init; /* Mask Parameter: InverterControl_Vdc_ref_Init
                                        * Referenced by:
                                        *   '<S4>/Vnom_dc1'
                                        *   '<S8>/Iph_'
                                        */
  real_T InverterControl_Vnom_dc;      /* Mask Parameter: InverterControl_Vnom_dc
                                        * Referenced by: '<S13>/Rtot_pu2'
                                        */
  real_T InverterControl_Vnom_prim;    /* Mask Parameter: InverterControl_Vnom_prim
                                        * Referenced by:
                                        *   '<S10>/A->pu'
                                        *   '<S10>/V->pu'
                                        *   '<S12>/Constant3'
                                        */
  real_T CompareToConstant_const;      /* Mask Parameter: CompareToConstant_const
                                        * Referenced by: '<S57>/Constant'
                                        */
  real_T CompareToConstant1_const;     /* Mask Parameter: CompareToConstant1_const
                                        * Referenced by: '<S58>/Constant'
                                        */
  real_T Gain1_Gain;                   /* Expression: 0.5
                                        * Referenced by: '<S23>/Gain1'
                                        */
  real_T Gain1_Gain_h;                 /* Expression: 0.5
                                        * Referenced by: '<S26>/Gain1'
                                        */
  real_T Gain1_Gain_o;                 /* Expression: sps.K2
                                        * Referenced by: '<S28>/Gain1'
                                        */
  real_T Gain_Gain;                    /* Expression: sps.K1
                                        * Referenced by: '<S28>/Gain'
                                        */
  real_T Gain1_Gain_f;                 /* Expression: 0.5
                                        * Referenced by: '<S39>/Gain1'
                                        */
  real_T Gain1_Gain_d;                 /* Expression: 0.5
                                        * Referenced by: '<S42>/Gain1'
                                        */
  real_T Gain_Y0;                      /* Expression: [1]
                                        * Referenced by: '<S31>/Gain'
                                        */
  real_T Gain1_Gain_e;                 /* Expression: 2
                                        * Referenced by: '<S35>/Gain1'
                                        */
  real_T Integ4_gainval;               /* Computed Parameter: Integ4_gainval
                                        * Referenced by: '<S38>/Integ4'
                                        */
  real_T Integ4_IC;                    /* Expression: 0
                                        * Referenced by: '<S38>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperSat;/* Expression: 1e6
                                         * Referenced by: '<S38>/To avoid division  by zero'
                                         */
  real_T Toavoiddivisionbyzero_LowerSat;/* Expression: eps
                                         * Referenced by: '<S38>/To avoid division  by zero'
                                         */
  real_T SFunction_P1_Size[2];         /* Computed Parameter: SFunction_P1_Size
                                        * Referenced by: '<S40>/S-Function'
                                        */
  real_T SFunction_P1;                 /* Expression: MaxDelay
                                        * Referenced by: '<S40>/S-Function'
                                        */
  real_T SFunction_P2_Size[2];         /* Computed Parameter: SFunction_P2_Size
                                        * Referenced by: '<S40>/S-Function'
                                        */
  real_T SFunction_P2;                 /* Expression: Ts
                                        * Referenced by: '<S40>/S-Function'
                                        */
  real_T SFunction_P3_Size[2];         /* Computed Parameter: SFunction_P3_Size
                                        * Referenced by: '<S40>/S-Function'
                                        */
  real_T SFunction_P3;                 /* Expression: InitialValue
                                        * Referenced by: '<S40>/S-Function'
                                        */
  real_T SFunction_P4_Size[2];         /* Computed Parameter: SFunction_P4_Size
                                        * Referenced by: '<S40>/S-Function'
                                        */
  real_T SFunction_P4;                 /* Expression: DFT
                                        * Referenced by: '<S40>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition;   /* Expression: 0
                                        * Referenced by: '<S39>/Unit Delay'
                                        */
  real_T Constant_Value;               /* Expression: 1/sps.Finit
                                        * Referenced by: '<S38>/Constant'
                                        */
  real_T UnitDelay1_InitialCondition;  /* Expression: sps.Vinit
                                        * Referenced by: '<S38>/Unit Delay1'
                                        */
  real_T Gain3_Gain;                   /* Expression: 2
                                        * Referenced by: '<S35>/Gain3'
                                        */
  real_T Integ4_gainval_h;             /* Computed Parameter: Integ4_gainval_h
                                        * Referenced by: '<S41>/Integ4'
                                        */
  real_T Integ4_IC_k;                  /* Expression: 0
                                        * Referenced by: '<S41>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperSa_n;/* Expression: 1e6
                                          * Referenced by: '<S41>/To avoid division  by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerSa_m;/* Expression: eps
                                          * Referenced by: '<S41>/To avoid division  by zero'
                                          */
  real_T SFunction_P1_Size_a[2];       /* Computed Parameter: SFunction_P1_Size_a
                                        * Referenced by: '<S43>/S-Function'
                                        */
  real_T SFunction_P1_a;               /* Expression: MaxDelay
                                        * Referenced by: '<S43>/S-Function'
                                        */
  real_T SFunction_P2_Size_k[2];       /* Computed Parameter: SFunction_P2_Size_k
                                        * Referenced by: '<S43>/S-Function'
                                        */
  real_T SFunction_P2_l;               /* Expression: Ts
                                        * Referenced by: '<S43>/S-Function'
                                        */
  real_T SFunction_P3_Size_h[2];       /* Computed Parameter: SFunction_P3_Size_h
                                        * Referenced by: '<S43>/S-Function'
                                        */
  real_T SFunction_P3_f;               /* Expression: InitialValue
                                        * Referenced by: '<S43>/S-Function'
                                        */
  real_T SFunction_P4_Size_l[2];       /* Computed Parameter: SFunction_P4_Size_l
                                        * Referenced by: '<S43>/S-Function'
                                        */
  real_T SFunction_P4_k;               /* Expression: DFT
                                        * Referenced by: '<S43>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_c; /* Expression: 0
                                        * Referenced by: '<S42>/Unit Delay'
                                        */
  real_T Constant_Value_c;             /* Expression: 1/sps.Finit
                                        * Referenced by: '<S41>/Constant'
                                        */
  real_T UnitDelay1_InitialCondition_c;/* Expression: sps.Vinit
                                        * Referenced by: '<S41>/Unit Delay1'
                                        */
  real_T RadDeg_Gain;                  /* Expression: 180/pi
                                        * Referenced by: '<S35>/Rad->Deg.'
                                        */
  real_T Saturation_UpperSat;          /* Expression: inf
                                        * Referenced by: '<S31>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: eps
                                        * Referenced by: '<S31>/Saturation'
                                        */
  real_T Gain1_Gain_hn;                /* Expression: 0.5
                                        * Referenced by: '<S45>/Gain1'
                                        */
  real_T Constant_Value_a[2];          /* Expression: [0.92 0]
                                        * Referenced by: '<S19>/Constant'
                                        */
  real_T dq_Y0[2];                     /* Expression: [0,0]
                                        * Referenced by: '<S59>/dq'
                                        */
  real_T dq_Y0_p[2];                   /* Expression: [0,0]
                                        * Referenced by: '<S60>/dq'
                                        */
  real_T SFunction_P1_Size_b[2];       /* Computed Parameter: SFunction_P1_Size_b
                                        * Referenced by: '<S87>/S-Function'
                                        */
  real_T SFunction_P1_b;               /* Expression: dest
                                        * Referenced by: '<S87>/S-Function'
                                        */
  real_T SFunction_P2_Size_d[2];       /* Computed Parameter: SFunction_P2_Size_d
                                        * Referenced by: '<S87>/S-Function'
                                        */
  real_T SFunction_P2_h;               /* Expression: priority2
                                        * Referenced by: '<S87>/S-Function'
                                        */
  real_T SFunction_P3_Size_d[2];       /* Computed Parameter: SFunction_P3_Size_d
                                        * Referenced by: '<S87>/S-Function'
                                        */
  real_T SFunction_P3_m;               /* Expression: st
                                        * Referenced by: '<S87>/S-Function'
                                        */
  real_T SFunction1_Value;             /* Expression: 0
                                        * Referenced by: '<S1>/S-Function1'
                                        */
  real_T SFunction_X0;                 /* Expression: 0
                                        * Referenced by: '<S1>/S-Function'
                                        */
  real_T Memory_X0;                    /* Expression: 0
                                        * Referenced by: '<S2>/Memory'
                                        */
  real_T SFunction_P1_Size_c[2];       /* Computed Parameter: SFunction_P1_Size_c
                                        * Referenced by: '<S88>/S-Function'
                                        */
  real_T SFunction_P1_i;               /* Expression: src
                                        * Referenced by: '<S88>/S-Function'
                                        */
  real_T SFunction_P2_Size_kw[2];      /* Computed Parameter: SFunction_P2_Size_kw
                                        * Referenced by: '<S88>/S-Function'
                                        */
  real_T SFunction_P2_o;               /* Expression: Data_width
                                        * Referenced by: '<S88>/S-Function'
                                        */
  real_T SFunction_P3_Size_n[2];       /* Computed Parameter: SFunction_P3_Size_n
                                        * Referenced by: '<S88>/S-Function'
                                        */
  real_T SFunction_P3_n;               /* Expression: st
                                        * Referenced by: '<S88>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_a; /* Expression: sps.Finit
                                        * Referenced by: '<S30>/Unit Delay'
                                        */
  real_T avoiddivisionbyzero_UpperSat; /* Expression: 70
                                        * Referenced by: '<S19>/avoid division by zero'
                                        */
  real_T avoiddivisionbyzero_LowerSat; /* Expression: 40
                                        * Referenced by: '<S19>/avoid division by zero'
                                        */
  real_T Gain_Gain_c;                  /* Expression: 1/4
                                        * Referenced by: '<S19>/Gain'
                                        */
  real_T SFunction_P1_Size_f[2];       /* Computed Parameter: SFunction_P1_Size_f
                                        * Referenced by: '<S56>/S-Function'
                                        */
  real_T SFunction_P1_k;               /* Expression: MaxDelay
                                        * Referenced by: '<S56>/S-Function'
                                        */
  real_T SFunction_P2_Size_j[2];       /* Computed Parameter: SFunction_P2_Size_j
                                        * Referenced by: '<S56>/S-Function'
                                        */
  real_T SFunction_P2_g;               /* Expression: Ts
                                        * Referenced by: '<S56>/S-Function'
                                        */
  real_T SFunction_P3_Size_h5[2];      /* Computed Parameter: SFunction_P3_Size_h5
                                        * Referenced by: '<S56>/S-Function'
                                        */
  real_T SFunction_P3_c;               /* Expression: InitialValue
                                        * Referenced by: '<S56>/S-Function'
                                        */
  real_T SFunction_P4_Size_o[2];       /* Computed Parameter: SFunction_P4_Size_o
                                        * Referenced by: '<S56>/S-Function'
                                        */
  real_T SFunction_P4_j;               /* Expression: DFT
                                        * Referenced by: '<S56>/S-Function'
                                        */
  real_T DiscreteTimeIntegrator_gainval;/* Computed Parameter: DiscreteTimeIntegrator_gainval
                                         * Referenced by: '<S30>/Discrete-Time Integrator'
                                         */
  real_T DiscreteTimeIntegrator_IC;    /* Expression: sps.Phase_Init*pi/180
                                        * Referenced by: '<S30>/Discrete-Time Integrator'
                                        */
  real_T Constant4_Value;              /* Expression: 2*pi
                                        * Referenced by: '<S30>/Constant4'
                                        */
  real_T FirstcycleofsimulationId092Iq0_;/* Expression: 0
                                          * Referenced by: '<S19>/First cycle of simulation Id=0.92, Iq=0'
                                          */
  real_T FirstcycleofsimulationId092Iq_o;/* Expression: 1
                                          * Referenced by: '<S19>/First cycle of simulation Id=0.92, Iq=0'
                                          */
  real_T Duk_Gain;                     /* Expression: sps.D
                                        * Referenced by: '<S51>/D*u(k)'
                                        */
  real_T Delay_x1_InitialCondition[2]; /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S51>/Delay_x1'
                                        */
  real_T C11_Gain;                     /* Expression: sps.C11
                                        * Referenced by: '<S54>/C11'
                                        */
  real_T Delay_x2_InitialCondition[2]; /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S51>/Delay_x2'
                                        */
  real_T C12_Gain;                     /* Expression: sps.C12
                                        * Referenced by: '<S54>/C12'
                                        */
  real_T UnitDelay2_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S4>/Unit Delay2'
                                        */
  real_T Iq_ref_Value;                 /* Expression: 0
                                        * Referenced by: '<S4>/Iq_ref'
                                        */
  real_T Integrator_gainval;           /* Computed Parameter: Integrator_gainval
                                        * Referenced by: '<S14>/Integrator'
                                        */
  real_T Integrator_IC;                /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S14>/Integrator'
                                        */
  real_T Gain1_Gain_l;                 /* Expression: 2
                                        * Referenced by: '<S15>/Gain1'
                                        */
  real_T Integ4_gainval_o;             /* Computed Parameter: Integ4_gainval_o
                                        * Referenced by: '<S22>/Integ4'
                                        */
  real_T Integ4_IC_a;                  /* Expression: 0
                                        * Referenced by: '<S22>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperSa_j;/* Expression: 1e6
                                          * Referenced by: '<S22>/To avoid division  by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerSa_l;/* Expression: eps
                                          * Referenced by: '<S22>/To avoid division  by zero'
                                          */
  real_T SFunction_P1_Size_fw[2];      /* Computed Parameter: SFunction_P1_Size_fw
                                        * Referenced by: '<S24>/S-Function'
                                        */
  real_T SFunction_P1_a2;              /* Expression: MaxDelay
                                        * Referenced by: '<S24>/S-Function'
                                        */
  real_T SFunction_P2_Size_a[2];       /* Computed Parameter: SFunction_P2_Size_a
                                        * Referenced by: '<S24>/S-Function'
                                        */
  real_T SFunction_P2_e;               /* Expression: Ts
                                        * Referenced by: '<S24>/S-Function'
                                        */
  real_T SFunction_P3_Size_l[2];       /* Computed Parameter: SFunction_P3_Size_l
                                        * Referenced by: '<S24>/S-Function'
                                        */
  real_T SFunction_P3_i;               /* Expression: InitialValue
                                        * Referenced by: '<S24>/S-Function'
                                        */
  real_T SFunction_P4_Size_n[2];       /* Computed Parameter: SFunction_P4_Size_n
                                        * Referenced by: '<S24>/S-Function'
                                        */
  real_T SFunction_P4_b;               /* Expression: DFT
                                        * Referenced by: '<S24>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_cj;/* Expression: 0
                                        * Referenced by: '<S23>/Unit Delay'
                                        */
  real_T Constant_Value_p;             /* Expression: 1/sps.Finit
                                        * Referenced by: '<S22>/Constant'
                                        */
  real_T UnitDelay1_InitialCondition_p;/* Expression: sps.Vinit
                                        * Referenced by: '<S22>/Unit Delay1'
                                        */
  real_T Gain3_Gain_i;                 /* Expression: 2
                                        * Referenced by: '<S15>/Gain3'
                                        */
  real_T Integ4_gainval_f;             /* Computed Parameter: Integ4_gainval_f
                                        * Referenced by: '<S25>/Integ4'
                                        */
  real_T Integ4_IC_l;                  /* Expression: 0
                                        * Referenced by: '<S25>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperSa_c;/* Expression: 1e6
                                          * Referenced by: '<S25>/To avoid division  by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerSa_b;/* Expression: eps
                                          * Referenced by: '<S25>/To avoid division  by zero'
                                          */
  real_T SFunction_P1_Size_l[2];       /* Computed Parameter: SFunction_P1_Size_l
                                        * Referenced by: '<S27>/S-Function'
                                        */
  real_T SFunction_P1_f;               /* Expression: MaxDelay
                                        * Referenced by: '<S27>/S-Function'
                                        */
  real_T SFunction_P2_Size_f[2];       /* Computed Parameter: SFunction_P2_Size_f
                                        * Referenced by: '<S27>/S-Function'
                                        */
  real_T SFunction_P2_d;               /* Expression: Ts
                                        * Referenced by: '<S27>/S-Function'
                                        */
  real_T SFunction_P3_Size_f[2];       /* Computed Parameter: SFunction_P3_Size_f
                                        * Referenced by: '<S27>/S-Function'
                                        */
  real_T SFunction_P3_h;               /* Expression: InitialValue
                                        * Referenced by: '<S27>/S-Function'
                                        */
  real_T SFunction_P4_Size_c[2];       /* Computed Parameter: SFunction_P4_Size_c
                                        * Referenced by: '<S27>/S-Function'
                                        */
  real_T SFunction_P4_g;               /* Expression: DFT
                                        * Referenced by: '<S27>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_l; /* Expression: 0
                                        * Referenced by: '<S26>/Unit Delay'
                                        */
  real_T Constant_Value_m;             /* Expression: 1/sps.Finit
                                        * Referenced by: '<S25>/Constant'
                                        */
  real_T UnitDelay1_InitialCondition_k;/* Expression: sps.Vinit
                                        * Referenced by: '<S25>/Unit Delay1'
                                        */
  real_T RadDeg_Gain_f;                /* Expression: 180/pi
                                        * Referenced by: '<S15>/Rad->Deg.'
                                        */
  real_T torad_Gain;                   /* Expression: pi/180
                                        * Referenced by: '<S10>/to-rad'
                                        */
  real_T Rff_Gain;                     /* Expression: RLff(1)
                                        * Referenced by: '<S7>/Rff '
                                        */
  real_T Lff_Gain;                     /* Expression: RLff(2)
                                        * Referenced by: '<S7>/Lff  '
                                        */
  real_T Rff_Gain_a;                   /* Expression: RLff(1)
                                        * Referenced by: '<S7>/Rff'
                                        */
  real_T Lff_Gain_c;                   /* Expression: RLff(2)
                                        * Referenced by: '<S7>/Lff'
                                        */
  real_T Saturation_UpperSat_c;        /* Expression: 1.5
                                        * Referenced by: '<S7>/Saturation'
                                        */
  real_T Saturation_LowerSat_k;        /* Expression: -1.5
                                        * Referenced by: '<S7>/Saturation'
                                        */
  real_T Iph_1_Value;                  /* Expression: Limits_MPPT(1)
                                        * Referenced by: '<S8>/Iph_1'
                                        */
  real_T Iph_2_Value;                  /* Expression: Limits_MPPT(2)
                                        * Referenced by: '<S8>/Iph_2'
                                        */
  real_T MPPT_On_Value;                /* Expression: 1
                                        * Referenced by: '<S2>/MPPT_On'
                                        */
  real_T Integ4_gainval_j;             /* Computed Parameter: Integ4_gainval_j
                                        * Referenced by: '<S28>/Integ4'
                                        */
  real_T Integ4_IC_i;                  /* Expression: 0
                                        * Referenced by: '<S28>/Integ4'
                                        */
  real_T K1_Value;                     /* Expression: sps.Delay
                                        * Referenced by: '<S28>/K1'
                                        */
  real_T SFunction_P1_Size_p[2];       /* Computed Parameter: SFunction_P1_Size_p
                                        * Referenced by: '<S29>/S-Function'
                                        */
  real_T SFunction_P1_km;              /* Expression: MaxDelay
                                        * Referenced by: '<S29>/S-Function'
                                        */
  real_T SFunction_P2_Size_ar[2];      /* Computed Parameter: SFunction_P2_Size_ar
                                        * Referenced by: '<S29>/S-Function'
                                        */
  real_T SFunction_P2_ge;              /* Expression: Ts
                                        * Referenced by: '<S29>/S-Function'
                                        */
  real_T SFunction_P3_Size_e[2];       /* Computed Parameter: SFunction_P3_Size_e
                                        * Referenced by: '<S29>/S-Function'
                                        */
  real_T SFunction_P3_p;               /* Expression: InitialValue
                                        * Referenced by: '<S29>/S-Function'
                                        */
  real_T SFunction_P4_Size_c0[2];      /* Computed Parameter: SFunction_P4_Size_c0
                                        * Referenced by: '<S29>/S-Function'
                                        */
  real_T SFunction_P4_f;               /* Expression: DFT
                                        * Referenced by: '<S29>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_b; /* Expression: 0
                                        * Referenced by: '<S28>/Unit Delay'
                                        */
  real_T K2_Value;                     /* Expression: sps.Freq
                                        * Referenced by: '<S28>/K2'
                                        */
  real_T UnitDelay1_InitialCondition_kn;/* Expression: sps.Vinit
                                         * Referenced by: '<S28>/Unit Delay1'
                                         */
  real_T Integ4_gainval_jy;            /* Computed Parameter: Integ4_gainval_jy
                                        * Referenced by: '<S44>/Integ4'
                                        */
  real_T Integ4_IC_h;                  /* Expression: 0
                                        * Referenced by: '<S44>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperSa_e;/* Expression: 1e6
                                          * Referenced by: '<S44>/To avoid division  by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerSa_a;/* Expression: eps
                                          * Referenced by: '<S44>/To avoid division  by zero'
                                          */
  real_T SFunction_P1_Size_o[2];       /* Computed Parameter: SFunction_P1_Size_o
                                        * Referenced by: '<S46>/S-Function'
                                        */
  real_T SFunction_P1_h;               /* Expression: MaxDelay
                                        * Referenced by: '<S46>/S-Function'
                                        */
  real_T SFunction_P2_Size_h[2];       /* Computed Parameter: SFunction_P2_Size_h
                                        * Referenced by: '<S46>/S-Function'
                                        */
  real_T SFunction_P2_hd;              /* Expression: Ts
                                        * Referenced by: '<S46>/S-Function'
                                        */
  real_T SFunction_P3_Size_j[2];       /* Computed Parameter: SFunction_P3_Size_j
                                        * Referenced by: '<S46>/S-Function'
                                        */
  real_T SFunction_P3_nb;              /* Expression: InitialValue
                                        * Referenced by: '<S46>/S-Function'
                                        */
  real_T SFunction_P4_Size_d[2];       /* Computed Parameter: SFunction_P4_Size_d
                                        * Referenced by: '<S46>/S-Function'
                                        */
  real_T SFunction_P4_a;               /* Expression: DFT
                                        * Referenced by: '<S46>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_j; /* Expression: 0
                                        * Referenced by: '<S45>/Unit Delay'
                                        */
  real_T Constant_Value_f;             /* Expression: 1/sps.Finit
                                        * Referenced by: '<S44>/Constant'
                                        */
  real_T UnitDelay1_InitialCondition_j;/* Expression: sps.Vinit
                                        * Referenced by: '<S44>/Unit Delay1'
                                        */
  real_T DiscreteDerivative_DenCoef[2];/* Expression: [ TcD  Ts-TcD ]
                                        * Referenced by: '<S32>/Discrete Derivative '
                                        */
  real_T DiscreteDerivative_InitialState;/* Expression: 0
                                          * Referenced by: '<S32>/Discrete Derivative '
                                          */
  real_T DiscreteTimeIntegrator_gainva_j;/* Computed Parameter: DiscreteTimeIntegrator_gainva_j
                                          * Referenced by: '<S32>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_UpperSat;/* Expression: Par_Limits(1)
                                          * Referenced by: '<S32>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_LowerSat;/* Expression: Par_Limits(2)
                                          * Referenced by: '<S32>/Discrete-Time Integrator'
                                          */
  real_T Saturation1_UpperSat;         /* Expression: Par_Limits(1)
                                        * Referenced by: '<S32>/Saturation1'
                                        */
  real_T Saturation1_LowerSat;         /* Expression: Par_Limits(2)
                                        * Referenced by: '<S32>/Saturation1'
                                        */
  real_T Gain10_Gain;                  /* Expression: 1/2/pi
                                        * Referenced by: '<S30>/Gain10'
                                        */
  real_T RateLimiter_RisingLim;        /* Computed Parameter: RateLimiter_RisingLim
                                        * Referenced by: '<S30>/Rate Limiter'
                                        */
  real_T RateLimiter_FallingLim;       /* Computed Parameter: RateLimiter_FallingLim
                                        * Referenced by: '<S30>/Rate Limiter'
                                        */
  real_T RateLimiter_IC;               /* Expression: sps.Finit
                                        * Referenced by: '<S30>/Rate Limiter'
                                        */
  real_T Delay_x1_InitialCondition_i;  /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S47>/Delay_x1'
                                        */
  real_T A11_Gain;                     /* Expression: sps.A11
                                        * Referenced by: '<S48>/A11'
                                        */
  real_T Delay_x2_InitialCondition_g;  /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S47>/Delay_x2'
                                        */
  real_T A12_Gain;                     /* Expression: sps.A12
                                        * Referenced by: '<S48>/A12'
                                        */
  real_T A21_Gain;                     /* Expression: sps.A21
                                        * Referenced by: '<S48>/A21'
                                        */
  real_T A22_Gain;                     /* Expression: sps.A22
                                        * Referenced by: '<S48>/A22'
                                        */
  real_T B11_Gain;                     /* Expression: sps.B11
                                        * Referenced by: '<S49>/B11'
                                        */
  real_T B21_Gain;                     /* Expression: sps.B21
                                        * Referenced by: '<S49>/B21'
                                        */
  real_T Duk_Gain_p;                   /* Expression: sps.D
                                        * Referenced by: '<S47>/D*u(k)'
                                        */
  real_T C11_Gain_l;                   /* Expression: sps.C11
                                        * Referenced by: '<S50>/C11'
                                        */
  real_T C12_Gain_k;                   /* Expression: sps.C12
                                        * Referenced by: '<S50>/C12'
                                        */
  real_T A11_Gain_g;                   /* Expression: sps.A11
                                        * Referenced by: '<S52>/A11'
                                        */
  real_T A12_Gain_l;                   /* Expression: sps.A12
                                        * Referenced by: '<S52>/A12'
                                        */
  real_T A21_Gain_k;                   /* Expression: sps.A21
                                        * Referenced by: '<S52>/A21'
                                        */
  real_T A22_Gain_p;                   /* Expression: sps.A22
                                        * Referenced by: '<S52>/A22'
                                        */
  real_T B11_Gain_n;                   /* Expression: sps.B11
                                        * Referenced by: '<S53>/B11'
                                        */
  real_T B21_Gain_m;                   /* Expression: sps.B21
                                        * Referenced by: '<S53>/B21'
                                        */
  real_T Constant1_Value;              /* Expression: 0
                                        * Referenced by: '<S19>/Constant1'
                                        */
  real_T Constant3_Value;              /* Expression: sps.Delay
                                        * Referenced by: '<S83>/Constant3'
                                        */
  real_T Constant1_Value_h;            /* Expression: sps.Period
                                        * Referenced by: '<S83>/Constant1'
                                        */
  real_T ib1_Gain;                     /* Expression: sps.Freq
                                        * Referenced by: '<S83>/1\ib1'
                                        */
  real_T LookupTable_XData[3];         /* Expression: [0 .5 1]
                                        * Referenced by: '<S83>/Lookup Table'
                                        */
  real_T LookupTable_YData[3];         /* Expression: [0 2 0]
                                        * Referenced by: '<S83>/Lookup Table'
                                        */
  real_T Constant2_Value;              /* Expression: 1
                                        * Referenced by: '<S83>/Constant2'
                                        */
  real_T Gain1_Gain_du;                /* Expression: 0.5
                                        * Referenced by: '<S61>/Gain1'
                                        */
  real_T UnitDelay_InitialCondition_i; /* Expression: 0.1684
                                        * Referenced by: '<S4>/Unit Delay'
                                        */
  real_T Gain_Gain_cv;                 /* Expression: -1
                                        * Referenced by: '<S66>/Gain'
                                        */
  real_T Constant2_Value_b;            /* Expression: 0
                                        * Referenced by: '<S12>/Constant2'
                                        */
  real_T UnitDelay3_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S4>/Unit Delay3'
                                        */
  real_T Gain1_Gain_p;                 /* Expression: 1
                                        * Referenced by: '<S12>/Gain1'
                                        */
  real_T UnitDelay1_InitialCondition_n;/* Expression: 0
                                        * Referenced by: '<S4>/Unit Delay1'
                                        */
  real_T Integrator_gainval_i;         /* Computed Parameter: Integrator_gainval_i
                                        * Referenced by: '<S85>/Integrator'
                                        */
  real_T Integrator_IC_i;              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S85>/Integrator'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_my_pv_system_4_ss_controller_T {
  const char_T *path;
  const char_T *modelName;
  struct SimStruct_tag * *childSfunctions;
  const char_T *errorStatus;
  SS_SimMode simMode;
  RTWLogInfo *rtwLogInfo;
  RTWExtModeInfo *extModeInfo;
  RTWSolverInfo solverInfo;
  RTWSolverInfo *solverInfoPtr;
  void *sfcnInfo;

  /*
   * NonInlinedSFcns:
   * The following substructure contains information regarding
   * non-inlined s-functions used in the model.
   */
  struct {
    RTWSfcnInfo sfcnInfo;
    time_T *taskTimePtrs[1];
    SimStruct childSFunctions[9];
    SimStruct *childSFunctionPtrs[9];
    struct _ssBlkInfo2 blkInfo2[9];
    struct _ssSFcnModelMethods2 methods2[9];
    struct _ssSFcnModelMethods3 methods3[9];
    struct _ssStatesInfo2 statesInfo2[9];
    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      real_T const *UPtrs0[1];
      real_T const *UPtrs1[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[4];
      mxArray *params[4];
      struct _ssDWorkRecord dWork[3];
      struct _ssDWorkAuxRecord dWorkAux[3];
    } Sfcn0;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      real_T const *UPtrs0[1];
      real_T const *UPtrs1[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[4];
      mxArray *params[4];
      struct _ssDWorkRecord dWork[3];
      struct _ssDWorkAuxRecord dWorkAux[3];
    } Sfcn1;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[4];
      uint_T attribs[3];
      mxArray *params[3];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn2;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[3];
      mxArray *params[3];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn3;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      real_T const *UPtrs0[1];
      real_T const *UPtrs1[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[4];
      mxArray *params[4];
      struct _ssDWorkRecord dWork[3];
      struct _ssDWorkAuxRecord dWorkAux[3];
    } Sfcn4;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      real_T const *UPtrs0[1];
      real_T const *UPtrs1[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[4];
      mxArray *params[4];
      struct _ssDWorkRecord dWork[3];
      struct _ssDWorkAuxRecord dWorkAux[3];
    } Sfcn5;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      real_T const *UPtrs0[1];
      real_T const *UPtrs1[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[4];
      mxArray *params[4];
      struct _ssDWorkRecord dWork[3];
      struct _ssDWorkAuxRecord dWorkAux[3];
    } Sfcn6;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      real_T const *UPtrs0[1];
      real_T const *UPtrs1[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[4];
      mxArray *params[4];
      struct _ssDWorkRecord dWork[3];
      struct _ssDWorkAuxRecord dWorkAux[3];
    } Sfcn7;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      real_T const *UPtrs0[1];
      real_T const *UPtrs1[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[4];
      mxArray *params[4];
      struct _ssDWorkRecord dWork[3];
      struct _ssDWorkAuxRecord dWorkAux[3];
    } Sfcn8;
  } NonInlinedSFcns;

  /*
   * ModelData:
   * The following substructure contains information regarding
   * the data used in the model.
   */
  struct {
    void *blockIO;
    const void *constBlockIO;
    void *defaultParam;
    ZCSigState *prevZCSigState;
    real_T *contStates;
    real_T *derivs;
    void *zcSignalValues;
    void *inputs;
    void *outputs;
    boolean_T *contStateDisabled;
    boolean_T zCCacheNeedsReset;
    boolean_T derivCacheNeedsReset;
    boolean_T blkStateChange;
    void *dwork;
  } ModelData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    uint32_T checksums[4];
    uint32_T options;
    int_T numContStates;
    int_T numU;
    int_T numY;
    int_T numSampTimes;
    int_T numBlocks;
    int_T numBlockIO;
    int_T numBlockPrms;
    int_T numDwork;
    int_T numSFcnPrms;
    int_T numSFcns;
    int_T numIports;
    int_T numOports;
    int_T numNonSampZCs;
    int_T sysDirFeedThru;
    int_T rtwGenSfcn;
  } Sizes;

  /*
   * SpecialInfo:
   * The following substructure contains special information
   * related to other components that are dependent on RTW.
   */
  struct {
    const void *mappingInfo;
    void *xpcData;
  } SpecialInfo;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    time_T stepSize;
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    time_T tStart;
    time_T tFinal;
    time_T timeOfLastOutput;
    void *timingData;
    real_T *varNextHitTimesList;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *sampleTimes;
    time_T *offsetTimes;
    int_T *sampleTimeTaskIDPtr;
    int_T *sampleHits;
    int_T *perTaskSampleHits;
    time_T *t;
    time_T sampleTimesArray[1];
    time_T offsetTimesArray[1];
    int_T sampleTimeTaskIDArray[1];
    int_T sampleHitArray[1];
    int_T perTaskSampleHitsArray[1];
    time_T tArray[1];
  } Timing;
};

/* Block parameters (auto storage) */
extern P_my_pv_system_4_ss_controller_T my_pv_system_4_ss_controller_P;

/* Block signals (auto storage) */
extern B_my_pv_system_4_ss_controller_T my_pv_system_4_ss_controller_B;

/* Block states (auto storage) */
extern DW_my_pv_system_4_ss_controller_T my_pv_system_4_ss_controller_DW;

/*====================*
 * External functions *
 *====================*/
extern my_pv_system_4_ss_controller_rtModel *my_pv_system_4_ss_controller(void);
extern void MdlInitializeSizes(void);
extern void MdlInitializeSampleTimes(void);
extern void MdlInitialize(void);
extern void MdlStart(void);
extern void MdlOutputs(int_T tid);
extern void MdlUpdate(int_T tid);
extern void MdlTerminate(void);

/* Real-time Model object */
extern RT_MODEL_my_pv_system_4_ss_controller_T *const
  my_pv_system_4_ss_controller_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'my_pv_system_4_ss_controller'
 * '<S1>'   : 'my_pv_system_4_ss_controller/OpCCode_do_not_touch'
 * '<S2>'   : 'my_pv_system_4_ss_controller/SS_controller'
 * '<S3>'   : 'my_pv_system_4_ss_controller/powergui'
 * '<S4>'   : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control'
 * '<S5>'   : 'my_pv_system_4_ss_controller/SS_controller/rtlab_send_subsystem'
 * '<S6>'   : 'my_pv_system_4_ss_controller/SS_controller/zzzOpComm'
 * '<S7>'   : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/Current Regulator'
 * '<S8>'   : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/MPPT  Parameters'
 * '<S9>'   : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/MPPT Controller using Perturbe  & Observe technique  '
 * '<S10>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements'
 * '<S11>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator'
 * '<S12>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/U_ref Generation '
 * '<S13>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/VDC Regulator'
 * '<S14>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/Current Regulator/PI'
 * '<S15>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)'
 * '<S16>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Mean'
 * '<S17>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL'
 * '<S18>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Second-Order Filter'
 * '<S19>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Single-Phase dq Transform'
 * '<S20>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)1'
 * '<S21>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)2'
 * '<S22>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model'
 * '<S23>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model/Correction subsystem'
 * '<S24>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model/Discrete Variable Time Delay'
 * '<S25>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model'
 * '<S26>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model/Correction subsystem'
 * '<S27>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model/Discrete Variable Time Delay'
 * '<S28>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Mean/Model'
 * '<S29>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Mean/Model/Discrete Variable Time Delay'
 * '<S30>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model'
 * '<S31>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control'
 * '<S32>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Discrete'
 * '<S33>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Mean (Variable Frequency)'
 * '<S34>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Second-Order Filter'
 * '<S35>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)'
 * '<S36>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)1'
 * '<S37>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)2'
 * '<S38>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model'
 * '<S39>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model/Correction subsystem'
 * '<S40>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model/Discrete Variable Time Delay'
 * '<S41>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model'
 * '<S42>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model/Correction subsystem'
 * '<S43>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model/Discrete Variable Time Delay'
 * '<S44>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Mean (Variable Frequency)/Model'
 * '<S45>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Mean (Variable Frequency)/Model/Correction subsystem'
 * '<S46>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Mean (Variable Frequency)/Model/Discrete Variable Time Delay'
 * '<S47>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Second-Order Filter/Model'
 * '<S48>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Second-Order Filter/Model/A*k(k-1)'
 * '<S49>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Second-Order Filter/Model/B*(u(k)+u(k-1))'
 * '<S50>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Second-Order Filter/Model/C*x(k)'
 * '<S51>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Second-Order Filter/Model'
 * '<S52>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Second-Order Filter/Model/A*k(k-1)'
 * '<S53>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Second-Order Filter/Model/B*(u(k)+u(k-1))'
 * '<S54>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Second-Order Filter/Model/C*x(k)'
 * '<S55>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Alpha-Beta-Zero to dq0'
 * '<S56>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Discrete Variable Time Delay'
 * '<S57>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S58>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S59>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S60>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S61>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Cr_MinMax'
 * '<S62>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Modulator type'
 * '<S63>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Reference signal'
 * '<S64>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling'
 * '<S65>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Modulator type/Full Bridge Bipolar'
 * '<S66>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Modulator type/Full Bridge Unipolar'
 * '<S67>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Modulator type/One Three Phase Bridge'
 * '<S68>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Reference signal/External'
 * '<S69>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Reference signal/Internal'
 * '<S70>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Sync Asymmetrical'
 * '<S71>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Sync Natural'
 * '<S72>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Sync Symmetrical'
 * '<S73>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Unsync Asymmetrical'
 * '<S74>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Unsync Natural'
 * '<S75>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Unsync Symmetrical'
 * '<S76>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Sync Asymmetrical/Sample & Hold'
 * '<S77>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Sync Natural/Sync_NaturalSampling'
 * '<S78>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Sync Symmetrical/Sync_SymmetricalSampling'
 * '<S79>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Sync Symmetrical/Sync_SymmetricalSampling/Sample & Hold'
 * '<S80>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Unsync Asymmetrical/Unsync_AsymmetricalSampling'
 * '<S81>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Unsync Natural/Unsync_NaturalSampling'
 * '<S82>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Unsync Natural/Unsync_NaturalSampling/Triangle Generator'
 * '<S83>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Unsync Natural/Unsync_NaturalSampling/Triangle Generator/Model'
 * '<S84>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Unsync Symmetrical/Unsync_SymmetricalSampling'
 * '<S85>'  : 'my_pv_system_4_ss_controller/SS_controller/Inverter Control/VDC Regulator/PI'
 * '<S86>'  : 'my_pv_system_4_ss_controller/SS_controller/rtlab_send_subsystem/Subsystem3'
 * '<S87>'  : 'my_pv_system_4_ss_controller/SS_controller/rtlab_send_subsystem/Subsystem3/Send3'
 * '<S88>'  : 'my_pv_system_4_ss_controller/SS_controller/zzzOpComm/Receive_1'
 * '<S89>'  : 'my_pv_system_4_ss_controller/SS_controller/zzzOpComm/busStruct'
 * '<S90>'  : 'my_pv_system_4_ss_controller/SS_controller/zzzOpComm/busStruct/Sub1'
 * '<S91>'  : 'my_pv_system_4_ss_controller/SS_controller/zzzOpComm/busStruct/Sub2'
 */
#endif                                 /* RTW_HEADER_my_pv_system_4_ss_controller_h_ */
